#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import itertools
import json
import copy
import argparse
import networkx as nx
from operator import methodcaller
#from networkx import Graph
#from networkx import readwrite
from Bio import bgzf
import DOMLight, json
from collections import deque
from collections import OrderedDict
from cStringIO import StringIO
# heap analysis from guppy import hpy
#requires 2.7 or greater
if sys.version_info < (2, 7):
    raise "must use python 2.7 or greater"

from lxml.etree import Element, ElementTree, tostring, fromstring, register_namespace, CDATA
#try:
#    from xml.etree.cElementTree import Element, ElementTree, tostring, fromstring, register_namespace
#except ImportError:
#    try:
#        from xml.etree.ElementTree import Element, ElementTree, tostring, fromstring, register_namespace
#    except ImportError:
#        pass

##This script expects files which contain a list of consecutive protein family ID's for the replicons in an organism
## and some kind of summary information about where the kmers come from

##NCBI_TAX_ID     RANK    TAX_PATH
##590     genus   220341,90370,59201,28901,590
#genome_name	genome_info_id	ncbi_tax_id	taxon_lineage_ids

##NAME    NCBI_TAX_ID     ACCESSION       START_MIN
##FIG01045527     946034  AERV01000001    507
#fam_id		gid	ncbi_tax_id	sequence_info_id	start_max	figfam_product

#SHOULD BE

ip={'org_id':0,'contig_id':1,'locus_id':2,'start':3, 'end':4, 'fam_id':5}
fi={'fam_id':0,'fam_description':1}

#Edge Classes by reverse status. Here indexed to zero. Class 1: Forward, Forward; Class2: Forward, Reverse; Class3:Reverse, Forward; Class4:Reverse, Reverse
edgeClass={(False,False):1,(False,True):2,(True,False):4,(True,True):8}
edgePossible=set([1,2,4,8])

def warning(*objs):
    for o in objs:
        print >> sys.stderr, o

##Class for storing information about the origin of a Kmer
class featureInfo():
        #expand to parse out this information from different sources
    def __init__(self, line=None):
        self.contig_id=None
        self.genome_id=None
        self.feature_id=None
        self.group_id=None
        self.group_num=None
        self.start=None
        self.end=None
        self.rf_forward=None # when this feature leaves out of its kmer window (left side) its in transition to this rf-node
        self.rf_reverse=None # when this feature leaves out of the kmer window (right side) its in transition to this rf-node
        self.pg_assignment=None
        self.repeat_num=1 #the local repeat number for this 'character'
        self.instance_key=None
        self.rf_ends=[]#Stores orieintation information as which rf-nodes this feature enters/exits

    def addRFPointer(self, direction, pointer):
        if direction=="increase":
            self.rf_forward=pointer
        else:
            self.rf_reverse=pointer

    def addEndInfo(self, rf_id, rhs, direction):
        self.rf_ends.append((rf_id, rhs, direction)) # storing the id of the feature on the rhs for convenience. effectively tells you what side *this* feature is on


    def compareInstance(self, other_feature):
        count1=self.instance_key.count(".")
        count2=other_feature.instance_key.count(".")
        short_str = long_str = None
        if count1 == count2 == self.ksize:
            return self.instance_key == other_feature.instance_key
        elif count1 >= count2:
            short_str = other_feature.instance_key
            long_str = self.instance_key
        else:
            short_str = self.instance_key
            long_str = other_feature.instance_key

        short_option2= ".".join(short_str.split(".")[::-1])  # reverse

        return long_str.startswith(short_str) or long_str.startswith(short_option2) or long_str.endswith(short_str) or long_str.endswith(short_option2)



    def getContextValue(self, context):
        if context=="genome":
            return self.genome_id
        elif context=="contig":
            return self.contig_id
        else:
            return None

    def getString(self):
        #result="|".join([self.contig_id,self.org_id,str(self.start),str(self.end),self.function,self.fam_id])
        #return result
        return line
    def getParts(self):
        return line.split("\t")

    def dict_to_line(self, target):
        line_list=[None for i in range(len(ip))]
        for k in ip:
            line_list[ip[k]]=str(target[k])
        return "\t".join(line_list)

    #calculate region between genes
    def getInterFeature(self,nxt_feature):
        #print "from "+self.fam_id+" to "+nxt_feature.fam_id
        cur_info=self.parse_line(self.line)
        nxt_info=self.parse_line(nxt_feature.line)
        edge_info={}
        edge_info['contig_id']=cur_info['contig_id']
        edge_info['org_id']=cur_info['org_id']
        edge_info['start']=int(cur_info['end'])+1
        edge_info['end']=int(nxt_info['start'])-1
        edge_info['fam_id']="EDGE"
        edge_info['locus_id']="None"
        result=geneInfo(self.dict_to_line(edge_info))
        return result

    def parse_line(self, line=None):
        result={}
        if line == None:
            line =self.line
        try:
            parts=line.strip().split("\t")
            result['fam_id']=parts[ip['fam_id']]
            result['contig_id']=parts[ip['contig_id']]
            result['org_id']=parts[ip['org_id']]
            result['start']=parts[ip['start']]
            result['end']=parts[ip['end']]
        except:
            warning("parsing problem. couldn't parse line: "+line)
        return result
    def getLocation(self):
        result=self.parse_line(self.line)
        return [result['contig_id'],result['start'],result['end']]
        #return [self.contig_id, self.start, self.end]
    def getLocationString(self):
        result=self.parse_line(self.line)
        return ":".join([result['contig_id'],str(result['start']),str(result['end'])])
        #return ":".join([self.contig_id, str(self.start), str(self.end)])
    def getFeatureString(self, delim=":"):
        result=self.parse_line(self.line)
        return delim.join([result['contig_id'],str(result['start']),str(result['end']), str(result['fam_id']), result['org_id']])
        #return delim.join([self.contig_id, str(self.start), str(self.end), str(fam_id), str(org_id)])
    def getReplicon(self):
        result=self.parse_line(self.line)
        return result['contig_id']
    def getOrganism(self):
        result=self.parse_line(self.line)
        return result['org_id']
    def getFam(self):
        result=self.parse_line(self.line)
        return result['fam_id']
        #return self.org_id

##Class for storing all the geneInfo in a particular node
##along with the kmer information. does not store information in direction
##specific way.
class rfNode():
    def __init__(self, nodeID, feature_list, ksize, reverse, palindrome):
        self.features=[set([]),set([])]# first position represents a k-lengthed series of features in the positive direction; the second, in reverse
        self.positive_features=[0]
        self.negative_features=[1]
        self.assigned_features=[set([]),set([])]# first position represents a k-lengthed series of features in the positive direction; the second, in reverse
        self.addFeatures(reverse, feature_list)
        self.duplicate=False #whether this node is duplicated in any context bin
        self.nodeID=nodeID
        self.palindrome=palindrome
        self.split=False
        self.has_forward=False
        self.has_reverse=False
        #dfs non-recursive variables
        #self.done=False
        #self.visited=False
        #self.descending=True
        #old variables
        self.weightLabel=None
        self.weight=None
        self.linkOut={}#four classes of edges
        self.visited=False
        self.queued=False
        self.self_edge=False
        #self.curRevStatus=rev_status
    
    def bidirectional(self):
        return (len(self.features[0]) > 0 or len(self.assigned_features[0]) > 0) and (len(self.features[1]) > 0 or len(self.assigned_features[1]) > 0)

    def anchorNode(self):
        return (not self.duplicate) and (not self.palindrome)
    def numFeatures(self):
        return len(self.features[0])+len(self.features[1])

    def addFeatures(self, reverse, feature_list):
        if(reverse):
            self.features[-1].add(feature_list[-1])#for space efficency only store right most feature in kmer
            self.has_reverse=True
        else:
            self.features[0].add(feature_list[-1])#for space efficency only store right most feature in kmer
            self.has_forward=True
    #each cell in list stores info[LetterOfKmer]=geneInfo()
    def addInfo(self, position, cur_fam, info):
        if self.infoList[position] != None:
            if self.infoList[position][0]!=cur_fam:
                sys.stderr.write("logical error: trying to insert information about wrong family\n")
                sys.exit()
            self.infoList[position][-1].append(info)
        else:
            self.infoList[position]=(cur_fam,[info])
    #add intergenic information to what will eventually become pan-genome edges
    def addPGEInfo(self, inter_info, position):
        if position > len(self.peInfo)-1:
            print "out of bounds "+str(position)+" for "+" ".join(self.peInfo)
        else:
            self.peInfo[position].add(inter_info)

    def addEdges(self,node_id,nxt_rev_status):
        #get class of edge type
        #if self.nodeID ==1 and node_id ==0:
        #	print "Debug: makes no sense for this to link backwards"
        edgeStatus=edgeClass[(self.curRevStatus,nxt_rev_status)]
        if node_id in self.linkOut:
            self.linkOut[node_id]=self.linkOut[node_id]|edgeStatus #bitwise OR to represent both multi status
        else:
            self.linkOut[node_id]=edgeStatus

    #if the node has been visited before update its references
    def updateNode(self, prev_node, in_edge_status, storage):
        update_pos=[] #ordered pg-node references to project onto current node	
        if (not in_edge_status in edgePossible):
            sys.stderr.write("unforseen case: transitioning from "+"|".join(prev_node.infoList.keys())+" to "+"|".join(self.infoList.keys()))
        #update references to pg-nodes from overlapping portion of previous k-mer
        if in_edge_status & 1:
            update_pos = range(1,len(prev_node.pgRefs),1)+[None]
        elif in_edge_status & 2:
            update_pos = [None]+range(len(prev_node.pgRefs)-1,0,-1)
        elif in_edge_status & 4:
            update_pos = range(len(prev_node.pgRefs)-2,-1,-1)+[None]
        elif in_edge_status & 8:
            update_pos = [None]+range(0,len(prev_node.pgRefs)-1,1)
        for cur_pos, prev_pos in enumerate(update_pos):
            if prev_pos != None:
                if self.pgRefs[cur_pos] == None: #happens if already queued. transfer the reference
                    self.pgRefs[cur_pos]=prev_node.pgRefs[prev_pos]
                elif prev_node.pgRefs[prev_pos] != self.pgRefs[cur_pos]:
                    storage.updatePGNode(prev_node.pgRefs[prev_pos], self.pgRefs[cur_pos])
                
            
    #if the node has not been visited before transfer previous references
    def transferRefs(self, prev_node, in_edge_status, storage):
        #add references to pg-nodes from overlapping portion of previous k-mer
        #if in_edge_status & 1:
        #	for n in prev_node.pgRefs[1:]:self.pgRefs.append(n)
        #elif in_edge_status & 2:
        #	for n in reversed(prev_node.pgRefs[1:]):self.pgRefs.append(n)
        #elif in_edge_status & 4:
        #	for n in reversed(prev_node.pgRefs[0:-1]):self.pgRefs.append(n)
        #elif in_edge_status & 8:
        #	for n in prev_node.pgRefs[0:-1]):self.pgRefs.append(n)
        #if self.nodeID == 1 or self.nodeID ==2 or (prev_node != None and (prev_node.nodeID ==1 or prev_node.nodeID ==2)):
        #	print "Debug: pgRefs and in_edge_status screwed up"
        update_pos=[] #ordered pg-node references to project onto current node	
        #update references to pg-nodes from overlapping portion of previous k-mer
        if in_edge_status & 1:
            update_pos = range(1,len(prev_node.pgRefs),1)+[None]
        elif in_edge_status & 2:
            update_pos = [None]+range(len(prev_node.pgRefs)-1,0,-1)
        elif in_edge_status & 4:
            update_pos = range(len(prev_node.pgRefs)-2,-1,-1)+[None]
        elif in_edge_status & 8:
            update_pos = [None]+range(0,len(prev_node.pgRefs)-1,1)
        for cur_pos, prev_pos in enumerate(update_pos):
            if prev_pos != None:
                if self.pgRefs[cur_pos] == None: #happens if already queued. transfer the reference
                    self.pgRefs[cur_pos]=prev_node.pgRefs[prev_pos]
                elif prev_node.pgRefs[prev_pos] != self.pgRefs[cur_pos]:
                    storage.updatePGNode(prev_node.pgRefs[prev_pos], self.pgRefs[cur_pos])

    #apply this kmers location info to current pg-node references
    def applyInfo(self,storage):
        #infoList is an OrderedDict
        #if self.nodeID == 1 or self.nodeID ==2:
        #	print "Debug: pgRefs and in_edge_status screwed up"
        for count,info in enumerate(self.infoList):
            nid=self.pgRefs[count]
            storage.addInfoPGNode(nid,info[-1])#adds the node. edges are implied within every k-mer 

    def addPGEdges(self,storage):
        #if self.nodeID == 1 or self.nodeID ==2 :
        #	print "Debug: pgRefs and in_edge_status screwed up"
        for i in range(0,len(self.pgRefs)-1,1):
            if self.pgRefs[i] == None or self.pgRefs[i+1] == None:
                print "missing pg-nodes in "+str(self.nodeID)
                sys.exit()
            if len(self.peInfo[i]):
                storage.getPGNode(self.pgRefs[i]).addEdge(self.pgRefs[i+1],self.peInfo[i])
        
    #1st process previous knode using incoming direction edge to put ref in this kmer. And add this kmers labels to previous references.
    #2nd Add edges to new family added in this kmer FOR ALL INCOMING EDGE TYPES
    #if there is no previous node just straight expand it
    #3rd Check outbound k nodes to see if identity process necessary (in BFS)
    #4th when checking outbound k nodes see if prev_node == next_node OR cur_node == next_node
    #NOTES direction does not matter at the pg-edge/node level
    #model letters in k-mer more explicitly than stupid | separated     	
    def visitNode(self, prev_node, in_edge_status, storage):
        #if self.nodeID == 1 or self.nodeID ==2 or (prev_node != None and (prev_node.nodeID ==1 or prev_node.nodeID ==2)):
        #	print "Debug: pgRefs and in_edge_status screwed up"
        #if self.nodeID ==3261:
        #	print "Debug: investigate here"
        if prev_node == None:
            for count,info in enumerate(self.infoList):
                g_id=storage.addPGNode(info[0],info[-1])#adds the node. edges are implied within every k-mer 
                self.pgRefs[count]=g_id
        else:
            
            if (not in_edge_status in edgePossible):
                sys.stderr.write("unforseen case: transitioning from "+"|".join([x[0] for x in prev_node.infoList])+" to "+"|".join([x[0] for x in self.infoList]))

            #if the beginnning of this kmer is new create a pg-node for it and a reference to it in this kmer
            #handle new portion exposed in this kmer
            #case 2|8 =10
            if (in_edge_status & 10):
                g_id=storage.addPGNode(self.infoList[0][0],self.infoList[0][-1])
                self.pgRefs[0]=g_id

            #transfer references from previous k-mer
            #if self.nodeID==3261:
            #	print "Debug: look at transfer of references to this node"
            self.transferRefs(prev_node, in_edge_status, storage)
    
            #if the end of *this* kmer is new create a PG-node for it and add the reference to this kmer
            #handle new portion exposed in this kmer
            #case 1|4 =5
            if (in_edge_status & 5):
                g_id=storage.addPGNode(self.infoList[-1][0],self.infoList[-1][-1])
                self.pgRefs[-1]=g_id
            #some information may be unique to this kmer. apply it to the pg-nodes
            self.applyInfo(storage)

        self.visited=True
                
    def getReplicons(self):
        result=set([])
        #all the replicons should be the same for each fam in this kmer
        for fam in self.infoList:
            for info in fam[-1]: 
                result.add(info.getReplicon())
        return result
    def testNode(self):
        #make sure that all the families in the kmer come from same replicons
        ref_set=set([info.getReplicon() for info in self.infoList[0][-1]])
        for tup in self.infoList:
            test_set=set([])
            for info in tup[-1]:
                test_set.add(info.getReplicon())
            if test_set != ref_set:
                warning("kmer "+self.nodeID+" has inconsistent replicons")
                
#pg-node "incubator" class
class pgShell():
    def __init__(self, nid,fid,gene_list):
        self.node_id=nid
        self.subsumed=False
        self.consumed_list=[]#ids of the things its consumed
        self.famSubset=famVersion(nid, fid,gene_list)
        self.edges={}#key is nodeRef, value is set of geneInfo intergenic
    def addEdge(self, nodeRef, e_info):
        if not nodeRef in self.edges:
            self.edges[nodeRef]=e_info.copy()
        else: self.edges[nodeRef].update(e_info)
    def addInfo(self, info_list):
        for i in info_list:
            self.famSubset.instances.add(i)
    def subsumeNode(self, target):
        for nid in target.edges:
            if nid in self.edges:
                self.edges[nid].update(target.edges[nid])
            else:
                self.edges[nid]=target.edges[nid] #Does this need to be copied?? It is a reference to a set after all...
        self.famSubset.instances.update(target.famSubset.instances)
        target.famSubset.subsumed=True
        self.consumed_list.append(target.node_id)
        self.consumed_list.extend(target.consumed_list)
        
#provides a summary of where this family occurs
#a family may be differentiated into multiple version depending on its ocurrence in kmers
class famVersion():
    def __init__(self, id, famID, id_list):
        self.summary_status=False
        self.id=id
        self.famID=famID
        self.instances=set(id_list) #set of locations that identify this version of family
        self.organisms=set()
        self.tax_summary=set()
        self.replicons=set()
        self.locations=set()
        self.functions=set()
    #returns of summary items
    def get_summary(self):
        if not self.summary_status:
            for i in self.instances:
                self.replicons.add(i.getReplicon())
                self.organisms.add(i.getOrganism())
                self.locations.add(i.getLocationString())
            self.summary_status=True
        result={"replicons":self.replicons, "organisms":self.organisms, "locations":self.locations, "functions":self.functions}
        return result
        
        
#storing information about each protein family
#function, name, locations
#organizes occurences into versions depending on kmer
class famInfo():
    def __init__(self, fID):
        self.fID=fID
        self.versions=[]#stores famSummary objects which detail locations
        self.label=""
        self.description=""
    #checks to see if the ID set that has changed now overlaps with any of the other sets
    #function returns the number to adjust original idx by to account for emptied sets
    def checkChainReaction(self, idx, fID, threshold, start=-1):
        debug=False
        num_adjust=0
        adjustment=True
        sets_merged=False #inefficient. should figure out which sets are merged and only updated those
        while adjustment:
            found=False
            for idx2, v in enumerate(self.versions):
                if idx2 != idx and idx2 > start: #id_set_list[0:idx]+id_set_list[idx+1:]:
                    intersect=self.versions[idx].instances.intersection(v.instances)
                    score=len(intersect)
                    if(score>=threshold):
                        if debug and fID == "FIG00638284":
                            warning("Merging groups for "+fID)
                            warning("Intersection", [':'.join(x.getLocation()) for x in intersect])
                            warning("Group1", [':'.join(x.getLocation()) for x in self.versions[idx].instances])
                            warning("Group2",[':'.join(x.getLocation()) for x in v.instances])
                        self.versions[idx].instances |= v.instances
                        v.instances=set([])
                        found=True
                        sets_merged=True
                        if idx2 < idx:
                            num_adjust=num_adjust+1
                        start=idx2
            adjustment= found
        #remove empty sets
        self.versions= [y for y in self.versions if len(y.instances)]				
        return (idx-num_adjust, sets_merged)

    #add id_set for an occurrence of the figfam in a kmer
    def add_instance(self, id_set, threshold, locationHash):
        matching_group=-1
        change_groups=False#keeps track of which groups need to be updated
        
        for idx, v in enumerate(self.versions):
            #bigset= id_set if len(id_set) > len(uid_set) else uid_set
            score=len(id_set.intersection(v.instances))
            if(score>=threshold):
                matching_group=idx
                #store the bigest id_set as the identifying one for this figfam
                v.instances |= id_set
                break
        if matching_group == -1:
            self.versions.append(famVersion(id_set))
            matching_group = len(self.versions)-1
        else:
            matching_group, change_groups=self.checkChainReaction(matching_group, self.fID, threshold)
        if not change_groups:
            for loc in self.versions[matching_group].instances:
                locationHash[loc]=(str(self.fID),str(matching_group))
        else:
            for idx_grp, v in enumerate(self.versions):
                for loc in v.instances:
                    locationHash[loc]=(str(self.fID),str(idx_grp))


        
class featureParser():
    def __init__(self, **kwargs):
        self.feature_file=kwargs['feature_file']
        self.file_type=kwargs['file_type']
        self.parse=None
        self.ip={'genome':0,'contig':1,'feature':2,'start':3, 'end':4, 'group':5}
        #self.ip={'taxid':2, 'genome':1, 'contig':3,'feature':2,'start':4, 'end':5, 'group':0}
        if self.file_type=="tab":
            self.parse=self.parseFeatureTab
    def parseFeatureTab(self):
        in_handle=open(self.feature_file)
        for line in in_handle:
            result=featureInfo()
            header=False
            try:
                header = line.startswith('#')
                if header:
                    #define column position based on header
                    parts=line.strip().replace("#","").split("\t")
                    x=0
                    while x < len(parts):
                        cur_part=parts[x].lower()
                        if cur_part in self.ip:
                            self.ip[cur_part]=x
                        x+=1
                    continue

                else:
                    parts=line.strip().split("\t")
                    result.group_id=parts[self.ip['group']]
                    result.contig_id=parts[self.ip['contig']]
                    result.genome_id=parts[self.ip['genome']]
                    result.start=int(parts[self.ip['start']])
                    #result.end=parts[ip['end']]
            except:
                warning("parsing problem. couldn't parse line: "+line)
                continue
            yield result


            
##CALCULATE DIVERSITY QUOTIENT!!! GENUS/TOTAL GENOMES
##CALCULATE NORMALIZED NUMBER WEIGHT of NUMBER OF genomes in edge/ total number of genomes

##This class is for storing dictionary structures that facilitate different pan-genome graph calculations
##Parameters are filepaths and the size of kmer to use
#GraphMaker(feature_tab=some_file, context="genome")
class GraphMaker():
    def __init__(self, **kwargs):#feature_file, family_file, summary_file, ksize, ignore_fams=set([])):
        #print feature_file
        #print summary_file
        #print str(ksize)
        self.feature_parser=None
        #convert option passed to file_type
        if "feature_tab" in kwargs:
            self.feature_parser=featureParser(feature_file=kwargs["feature_tab"], file_type="tab")
        self.context=kwargs["context"] #should be ["genome", "contig", "feature"]
        self.context_levels={"genome":0,"contig":1,"feature":2}
        self.ksize=kwargs["ksize"]
        self.break_conflict=kwargs["break_conflict"]
        self.num_pg_nodes=0
        self.rf_graph=nx.DiGraph()# the rf-graph (close to de bruijn) created from series of features with group designations
        self.pg_graph=nx.Graph()# pg-graph is an undirected grpah
        self.rf_node_index=[]
        self.replicon_map={}
        self.conflicts={}
        self.groups_seen={}
        self.group_index=[]
        self.context_bin=set([])
        self.feature_index=[]
        self.debug = False
        self.prebundle = False
        self.eat_repeats =False
        self.anchor_instance_keys={} #structured as {instance_key, [pg_assignment, [features with this instance key]]}
        self.non_anchor_guides={} # this is a lookup with the following structure [pg_id][rf_id]=feature_id. Allows looking of a guide_feature based on the pan-genome/transition to a particular rf_id. Should get limited use.
        #self.ignore_fams=ignore_fams
        self.kmerLevel=0 #the level of a kmer increases if it occurs in repeated series with itself
        self.kmerLookup={}#stores array for contig info and set for pointing to the next kmer 
        self.cur_rf_node=None
        self.prev_node=None
        self.pg_node_alt={}#for a given visit to an rf-node, for a given kmer position if there are multiple pg-nodes x=set([p1, p2, p3]) store it as pg_node_alt[pg1]=x, pg_node_alt[pg2]=x, etc.
        self.prev_indices=[]
        self.rf_starting_list=[]
        self.visit_number=0
        #based on the feature that is leaving the kmer: flip 0/1, orientation forward/reverse 0/1, leaving_position left/right 0/f1
        #gives position of the newest feature in the next kmer, the adjustment to the leaving feature to get the rhs of next kmer
        self.projection_table=[[[
            {"nxt_position":1,"rhs_adj":self.ksize,"feature_adj":self.ksize},# ++ lp=0
            {"nxt_position":0,"rhs_adj":-1,"feature_adj":-self.ksize} # ++ lp=k-1
            ],[
            {"nxt_position":1,"rhs_adj":-self.ksize,"feature_adj":-self.ksize},#projection_table[0][0][0] -- lp=0
            {"nxt_position":0,"rhs_adj":1,"feature_adj":self.ksize} #projection_table[0][0][1] -- lp=k-1
            ]],[[
            {"nxt_position":0,"rhs_adj":1,"feature_adj":self.ksize},# +- lp=0
            {"nxt_position":1,"rhs_adj":-self.ksize,"feature_adj":-self.ksize} # +- lp=k-1
            ],[
            {"nxt_position":0,"rhs_adj":-1,"feature_adj":-self.ksize},# -+ lp=0
            {"nxt_position":1,"rhs_adj":self.ksize,"feature_adj":self.ksize} # -+ lp=k-1
            ]]]

        #based on target info: orientation 0/1 forward/reverse, kmer_side 0/1 left/right 
        self.rhs_adj_table=[[
            {"new_feature_adj":-(self.ksize-1),"prev_feature_adj":-(self.ksize-2),"leaving_feature_adj":0},# increasing left side
            {"new_feature_adj":0,"prev_feature_adj":-1,"leaving_feature_adj":-(self.ksize-1)} # increasing right side
            ],[
            {"new_feature_adj":(self.ksize-1),"prev_feature_adj":(self.ksize-2),"leaving_feature_adj":0},#decreasing left side
            {"new_feature_adj":0,"prev_feature_adj":1,"leaving_feature_adj":self.ksize-1}#decreasing right side
            ]]

        #self.pg_initial=[] #initial node storage
        #self.pg_ptrs=[] #idx of nodes. for merging identity
        #self.figfamHash={}#stores sets of coordinates for each figfam used to distinguish between paralogs/orthologs/distant orthologs
        #self.summaryLookup={}
        #self.familyInfo={}
        #self.locationHash={}#stores the disambiguated 'version' of the protein family. hashed by (seq. accession, location)
        #self.geneHash={} #storing information about the individual genes
        #self.replicon_edges_dict={}#stores which replicons have which edges
        #self.summary_level=None#taxon level at which to summarize
        #self.ksize=ksize #size of the kmer to store
        #self.recentK=deque(maxlen=ksize-1)#used for elevating k-mers to the next level
        #self.replicon_map={}#stores relationships between org_ids and contig_ids (contig_ids)
        #self.parseFeatures(feature_file)
        #h=hpy()
        #print h.heap()	
        #self.parseSummary(summary_file)
        #self.parseFamilyInfo(family_file)



    def checkRFGraph(self):
        for r in self.rf_node_index:
            ambig=0
            if r.numFeatures() >0:
                ambig+=1
        print "rf-graph: "+str(ambig)+" nodes unexapanded"
                #assert LogicError("RFNode unexpanded")
    def checkPGGraph(self):
        for cnode in self.pg_graph.nodes_iter(data=True):
            if len(cnode[1]["features"]) == 0 :
                print "pg-graph node "+str(cnode[0])+"has no features"
            group_id=None
            for g in cnode[1]["features"]:
                for contig in cnode[1]["features"][g]:
                    for f in cnode[1]["features"][g][contig]:
                        if group_id == None:
                            group_id = self.feature_index[int(f)].group_id
                        elif group_id != self.feature_index[int(f)].group_id:
                            assert LogicError

    def calcStatistics(self):
        print "rf-graph:"
        print "nodes "+str(self.rf_graph.number_of_nodes())
        print "edges "+str(self.rf_graph.number_of_edges())
        print "pg-graph:"
        print "nodes "+str(self.pg_graph.number_of_nodes())
        print "edges "+str(self.pg_graph.number_of_edges())

    def finalizeGraphAttr(self):
        num_genomes=float(len(self.replicon_map.keys()))
        for e in self.pg_graph.edges_iter():
            attr=self.pg_graph.get_edge_data(*e)
            if "genomes" in attr:
                attr["weight"]=len(attr["genomes"])/num_genomes
            for a in attr:
                if type(attr[a])==set:
                    attr[a] = ','.join(attr[a])

    class taxInfo():
        def __init__(self, genome_name, summary_id):
            self.genome_name=genome_name
            self.summary_id=summary_id
        def get_summary_id(self):
            return self.summary_id
        
    ##adds a PGShell to pg_initial and a pointer in pg_ptrs
    def addPGNode(self,fid,gene_list):
        nid=len(self.pg_initial)
        self.pg_initial.append(pgShell(nid,fid,set(gene_list)))
        self.pg_ptrs.append(nid)
        return len(self.pg_ptrs)-1
    def getPGNode(self, node_idx):
        cur_ref=self.pg_ptrs[node_idx]
        result=self.pg_initial[cur_ref]
        #if result == None:
        #	print "Debug: None type"
        return result
    def addInfoPGNode(self, nid, gene_list):
        #if nid == None:
        #	print "Debug: whats going on?"
        cur_node=self.getPGNode(nid)
        cur_node.addInfo(gene_list)
    #using the idx provided make the main node subsume the target node
    #don't have to destroy target...
    def updatePGNode(self, main_idx, target_idx):
        main_idx2=self.pg_ptrs[main_idx]
        main_node = self.pg_initial[main_idx2]
        target_idx2 = self.pg_ptrs[target_idx]
        target_node = self.pg_initial[target_idx2]
        #if main_node == None or target_node == None or target_idx==18:
        #	print "Debug: None type here"
        if main_node.node_id != target_node.node_id:
            main_node.subsumeNode(target_node)
            #now point all future references to target_node at main_node
            for c in main_node.consumed_list:
                self.pg_ptrs[c]=main_idx2
            self.pg_initial[target_node.node_id]=None #destroy target
        

    #called after all features have been processed so that can flip them
    #and so that those that contain an anchor node can be added the anchor_instance_keys lookup
    def finalizeInstanceKeys(self):
        f=0
        while f < len(self.feature_index):
            if self.feature_index[f].instance_key !=None: #if its None then the contig wasn't big enough to create a node
                lv = [int(i) for i in self.feature_index[f].instance_key.split(".")]
                instance_reverse, instance_palindrome = self.flipKmer(lv)
                if instance_reverse:
                    lv.reverse()
                self.feature_index[f].instance_key_orig = ".".join(str(i) for i in lv)
                self.feature_index[f].instance_key = self.feature_index[f].instance_key_orig+"."+str(self.feature_index[f].repeat_num)+"r"
                # if the instance key has a single anchor node in it then it is an anchor instance key
                if self.feature_index[f].instance_key_orig in self.anchor_instance_keys:
                    self.anchor_instance_keys[self.feature_index[f].instance_key_orig][1].append(f)
                else:
                    #check to see if any of the instances are anchor nodes
                    for rf in lv:
                        if self.rf_node_index[rf].anchorNode():
                            self.anchor_instance_keys[self.feature_index[f].instance_key_orig]=[None,[f]]
                            break
            f+=1


    ##This function checks whether the kmer is in the graph
    #and links kmer graph data structure appropriately
    #store kmers according to the combined protein family ids, and a set of IDs for which kmer comes next
        #id of prev kmer, id of this kmer, information about this kmer, whether this kmer has been reversed
    def addRFNode(self, feature_list):
        reverse,palindrome,feature_indices,kmer_key=self.hashKmer(feature_list)#put IDs together to make kmer
        nodeID=None
        duplicate=False
        dup_number=0
        if not kmer_key in self.kmerLookup:
            nodeID = len(self.rf_node_index)
            self.rf_node_index.append(rfNode(nodeID, feature_indices, self.ksize, reverse, palindrome))
            self.cur_rf_node=self.rf_node_index[-1]
            self.kmerLookup[kmer_key]=self.rf_node_index[-1]
        else:
            duplicate=kmer_key in self.context_bin
            self.cur_rf_node=self.kmerLookup[kmer_key]
            self.cur_rf_node.addFeatures(reverse, feature_indices)
            if duplicate:
                self.cur_rf_node.duplicate=True
                dup_number=1
        self.context_bin.add(kmer_key)
        self.rf_graph.add_node(self.cur_rf_node.nodeID, label=kmer_key, duplicate=dup_number)

        #here we add a key that marks which kmers the feature occurs in for later disambiguation
        for f in feature_indices:
            if self.feature_index[f].instance_key == None:
                self.feature_index[f].instance_key= str(self.cur_rf_node.nodeID)
            else:
                self.feature_index[f].instance_key += "."+str(self.cur_rf_node.nodeID)

        #add information regarding which end of this rf-node the features are on
        self.feature_index[feature_indices[0]].addEndInfo(self.cur_rf_node.nodeID, feature_indices[self.ksize-1], reverse)
        self.feature_index[feature_indices[self.ksize-1]].addEndInfo(self.cur_rf_node.nodeID, feature_indices[self.ksize-1], reverse)

        #rf-edges. properties dictated by the relationship of the kmers (flipped or not)
        if self.prev_node!=None:
            if self.prev_reverse:
                if reverse: # -1 -1
                    leaving_position=self.ksize-1
                    reverse_lp=0
                else: # -1 +1
                    leaving_position=self.ksize-1
                    reverse_lp=self.ksize-1
            else:
                if reverse:# +1 -1
                    leaving_position=0
                    reverse_lp=0
                else:# +1 +1
                    leaving_position=0
                    reverse_lp=self.ksize-1

            self.feature_index[self.prev_indices[leaving_position]].addRFPointer(direction="increase", pointer=self.cur_rf_node.nodeID) #record which direction a feature is leaving the k-window and what rf-node it is traversing to
            self.feature_index[feature_indices[reverse_lp]].addRFPointer(direction="decrease", pointer=self.prev_node.nodeID) # to enable thread based navigation.
            if not self.rf_graph.has_edge(self.prev_node.nodeID, self.cur_rf_node.nodeID):
                rflip = fflip = self.prev_reverse ^ reverse #xor. if kmers are flipped relative to each other 
                if palindrome:
                    fflip = "no"
                if self.prev_node.palindrome:
                    rflip= "no"
                self.rf_graph.add_edge(self.prev_node.nodeID, self.cur_rf_node.nodeID, attr_dict={"flip":fflip,"leaving_position":leaving_position})
                self.rf_graph.add_edge(self.cur_rf_node.nodeID, self.prev_node.nodeID, attr_dict={"flip":rflip,"leaving_position":reverse_lp})
        self.prev_indices=feature_indices
        self.prev_node=self.cur_rf_node
        self.prev_reverse=reverse

    ##Create an ID for kmer
    ##In case directionality is flipped for entire genome I am flipping each kmer
    ##This shouldn't adversely affect inversions nor the overall result
    #takes a list of geneInfo objects
    def makeKey(self, k_info_list, prev_key):
        k_list=[]
        for k in k_info_list:
            k_list.append(k.getFam())
        id_sep="|"
        result=None
        rev_status=False
        if(k_list[0]> k_list[-1]):
            k_list.reverse()
            rev_status=True
        result=id_sep.join(k_list)
        if result in self.recentK:
            self.kmerLevel+=1
        else:
            self.kmerLevel=0
        self.recentK.append(result)
        return ((self.kmerLevel, result), rev_status)
                    
    ##Separate the kmer back into its parts
    def getParts(self, kmer):
        return(kmer.split('.'))
    
    #create a kmer based on the group ids in the features. Also append the feature to a feature_index
    def hashKmer(self, feature_list):
        kmer_hash=[]
        feature_indices=[]
        for feature in feature_list:
            if not feature.group_id in self.groups_seen:
                feature.group_num=len(self.group_index)
                self.groups_seen[feature.group_id]=feature.group_num
                self.group_index.append(feature.group_num)
            else:
                feature.group_num=self.groups_seen[feature.group_id]
            kmer_hash.append(feature.group_num)
            feature_indices.append(feature.feature_id)
        reverse, palindrome = self.flipKmer(kmer_hash)
        if reverse:
            feature_indices.reverse()
        return reverse, palindrome, feature_indices, ".".join([str(i) for i in kmer_hash])
            


    def flipKmer(self, feature_list):
        i=0
        k_size=len(feature_list)
        palindrome=0
        reverse=0
        while i<(k_size/2):
            if feature_list[i]< feature_list[k_size-(i+1)]:
                return (reverse, palindrome)
            elif feature_list[i] > feature_list[k_size-(i+1)]:
                reverse=1
                feature_list.reverse()
                return (reverse, palindrome)
            i+=1
        palindrome =1
        return (reverse, palindrome)


    def processFeatures(self):
        print "parsing features and constructing kmer graph"	
        kmer_q=deque()
        prev_feature=None
        #loop through figfams to create kmers
        repeat_num=1
        update_repeats=[]
        for feature in self.feature_parser.parse():
            if prev_feature and prev_feature.group_id == feature.group_id and prev_feature.contig_id == feature.contig_id:
                repeat_num+=1
                update_repeats.append(prev_feature)
                if self.eat_repeats:
                    continue
            else:
                if repeat_num >1:
                    update_repeats.append(prev_feature)
                    for to_up in update_repeats:
                        to_up.repeat_num =repeat_num
                repeat_num=1
                update_repeats=[]
            feature.feature_id= len(self.feature_index)
            self.feature_index.append(feature)
            if feature.genome_id not in self.replicon_map:
                self.replicon_map[feature.genome_id]=set()
            else:
                self.replicon_map[feature.genome_id].add(feature.contig_id)
            if(prev_feature and prev_feature.contig_id != feature.contig_id):
                kmer_q=deque()#clear kmer stack because switching replicons
                self.prev_node=None
                self.prev_indices=[]
            elif prev_feature and prev_feature.contig_id == feature.contig_id:
                if prev_feature.start > feature.start:
                    assert InputError
            #depending on the context populate the context bin with appropriate ids to detect duplicates
            if self.context:
                if(prev_feature and prev_feature.getContextValue(self.context) != feature.getContextValue(self.context)):
                    self.context_bin.clear()
            kmer_q.append(feature)#append the feature to the queue
            if(len(kmer_q)>self.ksize):
                kmer_q.popleft()
                self.addRFNode(kmer_q)
            elif(len(kmer_q)== self.ksize):
                self.addRFNode(kmer_q)#right now only passing in the last figfams information
            else:#kmer size is less than ksize
                kmer=None
            prev_feature=feature
        if self.debug: nx.set_node_attributes(self.rf_graph, "visit", "")


    def RF_to_PG(self):
        #determine list of starting nodes
        for node in self.rf_node_index:
            if node.anchorNode():
                self.rf_starting_list.append(node)
        self.finalizeInstanceKeys()
        #start with nodes that have the most features
        self.rf_starting_list.sort(key=lambda x: x.numFeatures(), reverse=True)
        for rf_node in self.rf_starting_list:
            self.tfs_expand_nr(None, rf_node, None, None)


        
    #for a given kmer return a set of the organisms involved
    def getOrgSummary(self, kmer):
        result=set()
        if kmer in self.kmerLookup:
            for i in self.kmerLookup[kmer][0].infoList.values()[0]:
                result.add(i.org_id)
        return result
        
    #for all the organisms of the kmer return a set of IDs at the taxonomy level
    #e.g. all the genus IDs stored in summaryLookup under various organism IDs
    def getTaxSummary(self,kmer):
        result=set()
        if kmer in self.kmerLookup:
            for i in self.kmerLookup[kmer][0].infoList.values()[0]:
                if(i.org_id in self.summaryLookup):
                    result.add(self.summaryLookup[i.org_id].get_summary_id())
        return result
    #for a given node return a set of the organisms involved
    def nodeOrgSummary(self,cnode):
        result=set()
        #print cnode.infoList
        for i in cnode.infoList.values()[0]:
            result.add(i.org_id)
        return result
    
    def nodeTaxSummary(self,cnode):
        result=set()
        for i in cnode.instances:
            org_id=i.getOrganism()
            if(org_id in self.summaryLookup):
                result.add(self.summaryLookup[org_id].get_summary_id())
        return result

    #Get the total number of unique taxonomy labels 
    def completeTaxSummary(self):
        result=set()
        for k in self.summaryLookup.keys():
            result.add(self.summaryLookup[k])
        return result

    #expects summary taxid, tax level, and the taxpath comma seperated values
    def parseSummary(self, summary_file):
        inHandle=open(summary_file, 'r')
        print "parsing taxonomy information and constructing taxon table"
        #header=inHandle.readline()
        for line in inHandle:
            if line.startswith('#'):
                continue
            summary_info=line.strip().split("\t")
            if(self.summary_level==None):
                self.summary_level=6
            ref_id=summary_info[1]
            genome_name=summary_info[0]
            summary_id=summary_info[-1].split(',')[self.summary_level]
            self.summaryLookup[ref_id]=self.taxInfo(genome_name,summary_id)
        inHandle.close()

    #expects two column family name information
    def parseFamilyInfo(self, family_file):
        in_handle=open(family_file, 'r')
        print "parsing family information table"
        #header=inHandle.readline()
        for line in in_handle:
            if line.startswith('#'):
                continue
            info_list=line.strip().split("\t")
            try:
                self.familyInfo[info_list[fi["fam_id"]]]=info_list[fi["fam_description"]]
            except:
                warning("problem parsing family info line: "+line)
                try:
                    self.familyInfo[info_list[fi["fam_id"]]]="No Label"
                except:
                    sys.exit()
                pass

        in_handle.close()
    def getFamilyInfo(self, fid):
        if fid in self.familyInfo:
            return self.familyInfo[fid]
        else: return None

    #i'm sorry universe this is going to be confusing
    #given the side of the kmer the feature is found on and its orientation relative to the kmer
    #return whether it leaves in the forward or reverse direction of the original kmer (before it was potentially flipped by hashing)
    #kmer_side=0 is left side, kmer_side=1 is right
    #orientation = 0 is increasing, orientation =1 is decreasing (feature series progression relative to kmer orientation) e.g 1,2,3 or 3,2,1
    def nextRFNode(self, kmer_side,orientation,feature):
        #left(0) and increasing(0) is the "back" of the kmer
        #right(1) and increasing(0) is the "back" of the kmer
        if (kmer_side==1 and orientation==0) or (kmer_side==0 and orientation==1):
            #progression="reverse"
            return self.feature_index[feature].rf_reverse

        #right(0) and decreasing(1) is the "front" of the kmer
        #left(1) and decreasing(1) is the "front" of the kmer
        else:
            #progression="forward"
            return self.feature_index[feature].rf_forward

    def projectFeature(self, cur_node, edge_data, kmer_side, orientation, leaving_feature, palindrome, nxt_node):
        to_queue=[] # used to process features with same instance key
        if (edge_data["leaving_position"]==self.ksize-1 and kmer_side!=1) or (edge_data["leaving_position"]==0 and kmer_side!=0):
            if not palindrome:
                sys.stderr.write("logic problem. calculated leaving side does not match")
                assert LogicError
        #project from leaving feature and orientation to what next feature should be next
        nxt_orientation=orientation
        flip=edge_data["flip"]
        if nxt_node.palindrome:#palindrome can't reliably use rf-edge info
            flip= orientation ^ 0 #xor. palindrome is always forward direction
        
        nxt_feature_info=self.projection_table[flip][orientation][kmer_side]
        nxt_feature=leaving_feature+nxt_feature_info["feature_adj"]
        
        if palindrome and nxt_node.bidirectional():
            #flip is uncertain
            k=self.ksize
            if nxt_feature > leaving_feature:
                if nxt_feature in nxt_node.features[0] or nxt_feature in nxt_node.assigned_features[0]:
                    nxt_orientation=0
                    nxt_position=1
                    nxt_target=nxt_feature
                elif (nxt_feature-(k-1)) in nxt_node.features[1] or nxt_feature-(k-1) in nxt_node.assigned_features[1]:
                    nxt_orientation=1
                    nxt_position=0
                    nxt_target=nxt_feature-(k-1)
            else:
                if (nxt_feature+(k-1)) in nxt_node.features[0] or nxt_feature+(k-1) in nxt_node.assigned_features[0]:
                    nxt_orientation = 0
                    nxt_position = 0
                    nxt_target=nxt_feature+(k-1)
                elif nxt_feature in nxt_node.features[1] or nxt_feature in nxt_node.assigned_features[1]:
                    nxt_orientation =1
                    nxt_position = 1
                    nxt_target=nxt_feature
                else:
                    assert LogicError
            to_queue.append((nxt_position, nxt_orientation, nxt_target))
        else:
            if flip==1:
                nxt_orientation = not orientation
            #projection_table: flip true/false, orientation forward/reverse true/false, leaving_position right/left true/false
            #"nxt_position" "feature_adj"
            nxt_target=leaving_feature+nxt_feature_info["rhs_adj"]
            to_queue.append((nxt_feature_info['nxt_position'], nxt_orientation, nxt_target))

        #expand projection/target selection based on instance key 
        instances_pkg = self.anchor_instance_keys.get(self.feature_index[nxt_feature].instance_key, [None,[]])
        #Turning off prebundling logic for now because it causes problems
        if self.prebundle and cur_node.nodeID != instances_pkg[0]:
            instances_pkg[0]=cur_node.nodeID
            for other_instance in instances_pkg[1]:
                instance_info=None
                #its gotta be one end or the other
                if self.feature_index[other_instance].rf_ends[0][0] == nxt_node.nodeID:
                    instance_info = self.feature_index[other_instance].rf_ends[0]
                else:
                    instance_info = self.feature_index[other_instance].rf_ends[1]
                rf_id, rhs, direction = instance_info
                position = 1 #set to rhs
                if other_instance!=rhs: # if its on the lhs
                    position = 0
                to_queue.append((position, direction, rhs))
        return to_queue


    #kmer_side=0 is left side, kmer_side=1 is right
    #orientation = 0 is increasing, orientation =1 is decreasing (feature series progression relative to kmer orientation) e.g 1,2,3 or 3,2,1
    #if a guide is passed, it is a feature from the leaving position. find its pg_node,
    #use that in combination with the nxt_rf_id to see if a guide/cat should be passed
    def queueFeature(self, cur_node, kmer_side, orientation, leaving_feature, node_queue, node_bundles, up_node=None):
        nxt_rf_id = self.nextRFNode(kmer_side, orientation, leaving_feature)
        prev_queued=True #has the rfid EVER been queued on THIS traversal
        up_queue= up_node!=None and nxt_rf_id == up_node.nodeID
        #If there are no unassigned features in the nxt node is there any point in visiting? conflict detection etc.?
        if nxt_rf_id!=None and self.rf_node_index[nxt_rf_id].numFeatures() > 0:
            bundle_id=nxt_rf_id
            cur_rf_id = cur_node.nodeID

            if up_queue:
                bundle_id=-1
            if(not bundle_id in node_bundles):
                node_bundles[bundle_id]=[[set([]),set([])],[set([]),set([])]]
                prev_queued=False
            #is the rfid currently queued on this traversal
            currently_queued=len(node_bundles[bundle_id][0][0])+len(node_bundles[bundle_id][0][1])+len(node_bundles[bundle_id][1][0])+len(node_bundles[bundle_id][1][1]) > 0
            #here look up edge information to project next
            edge_data=self.rf_graph[cur_node.nodeID][nxt_rf_id]
            to_queue = self.projectFeature(cur_node, edge_data,kmer_side,orientation,leaving_feature, cur_node.palindrome, self.rf_node_index[nxt_rf_id])
            for nxt_position, nxt_direction, nxt_target in to_queue:
                #structure for node_queue and node_bundles
                if not currently_queued and not up_queue: # if its being passed up (-1) then no need to queue
                    currently_queued=True
                    leaving_pg_id=self.feature_index[leaving_feature].pg_assignment
                    pg_node_id=self.feature_index[leaving_feature].pg_assignment
                    nxt_guide=nxt_guide_cat=nxt_guide_side=None
                    if prev_queued and leaving_pg_id in self.non_anchor_guides and nxt_rf_id in self.non_anchor_guides[leaving_pg_id]: #if rfnode previously been queued and not currently then its a re-descent and you need a guide to appropriately assign features to pg-nodes
                        nxt_guide, nxt_guide_cat, nxt_guide_side = self.non_anchor_guides[pg_node_id][nxt_rf_id]
                        #guide_obj = self.non_anchor_guides[leaving_pg_id][nxt_rf_id]
                    else:#first time queueing rfnode. store non_anchor_guides for later
                        if not pg_node_id in self.non_anchor_guides:
                            self.non_anchor_guides[pg_node_id]={nxt_rf_id:(nxt_target,nxt_direction,nxt_position)}
                        else:
                            self.non_anchor_guides[pg_node_id][nxt_rf_id]=(nxt_target, nxt_direction, nxt_position)
                    #no existing targets so needs to be queued. if its prev_queued then it will be queued with guide. else guide=None
                    if nxt_rf_id == cur_node.nodeID: #self loop goes first. 
                        node_queue.appendleft((nxt_rf_id,[nxt_guide, nxt_guide_cat, nxt_guide_side]))
                    else:
                        node_queue.append((nxt_rf_id,[nxt_guide, nxt_guide_cat, nxt_guide_side]))
                #node bundles exist separate from queue but are cleared out when rfid is taken from the queue
                node_bundles[bundle_id][nxt_position][nxt_direction].add(nxt_target)
            return nxt_rf_id
        return None

    #This should only be called on non-anchor nodes 
    def storeRevisitGuides(self, cur_rf_id, pre_assignments):
        #the leaving pg_node_id is used to distinguish the visit to an rf_node since an rf-node to rf-node transition can happen
        #multiple times. Instead of making this a pre-visit operation make it a visit operation. Where every
        #leaving pg-node points at the the same multi-pg-node array. That way no matter which is used
        #to look it up on re-descent it will give all pg-options.
        for i in [0,self.ksize-1]:
            for pg in pre_assignments[i]:
                for ik, features in pre_assignments[i][pg]["features"].iteritems():
                    for f in features:
                        cur_guide = f
                        cur_guide_cat = int(f < 0)
                        cur_guide_side = int(i == self.ksize-1)
                        #the same information that is used to queue a feature (that it is on the leaving side and what the next rf_node/feature will be 
                        #should be able to be used in REVERSE to project/find the feature that is used to come here. AND THUS the pg_node that 
                        #points here 
                        leaving_pg_id = None
                        multi_guide = {pg_id: {'nxt_guide':None,'nxt_guide_cat':None,'nxt_guide_side':None}}
                        self.non_anchor_guides.setdefault(leaving_pg_id,{}).setdefault(cur_rf_id,{}).setdefault(pg_id, dict(nxt_guide=cur_guide, nxt_guide_cat=cur_guide_cat, nxt_guide_side=cur_guide_side)) 


    def num_features_pg_node(self, node_id):
        num_features=0
        for g in self.pg_graph.node[node_id]['features']:
            for c in self.pg_graph.node[node_id]['features'][g]:
                num_features+=len(self.pg_graph.node[node_id]['features'][g][c])
        return num_features

    def merge_pg_node(self, node_id1, node_id2):
        print "merging "+str(node_id1)+" "+str(node_id2)
        conflict=False
        insert_level=None
        if node_id1 < node_id2:
            keep=node_id1
            remove=node_id2
        else:
            keep=node_id2
            remove=node_id1
        for g in self.pg_graph.node[remove]['features']:
            if not g in self.pg_graph.node[keep]['features']:
                self.pg_graph.node[keep]['features'][g]=self.pg_graph.node[remove]['features'][g]
                insert_level="genome"
            else:
                for c in self.pg_graph.node[remove]['features'][g]:
                    if not c in self.pg_graph.node[keep]['features'][g]:
                        self.pg_graph.node[keep]['features'][g][c] = self.pg_graph.node[remove]['features'][g][c]
                        insert_level="contig"
                    else:
                        insert_level="feature"
                        merge_set=set(self.pg_graph.node[keep]['features'][g][c])
                        for f in self.pg_graph.node[remove]['features'][g][c]:
                            if not f in merge_set:
                                merge_set.add(f)
                        self.pg_graph.node[keep]['features'][g][c]=list(merge_set)
            if self.context != "all" and insert_level != self.context:
                conflict=True
        for g in self.pg_graph.node[remove]['features']:
            for c in self.pg_graph.node[remove]['features'][g]:
                for f in self.pg_graph.node[remove]['features'][g][c]:
                    if keep == 74 and remove == 3161:
                        print "assigning "+str(f)+" to 74"
                    self.feature_index[f].pg_assignment=keep
        for e in self.pg_graph.edges(remove, data=True):
            if self.pg_graph.has_edge(keep, e[1]):
                cur_edge_data=self.pg_graph.get_edge_data(keep,e[1])
                for k in e[-1]:#edge dictionary
                    cur_edge_data[k].update(e[-1][k])
            else:
                self.pg_graph.add_edge(keep, e[1], attr_dict=e[-1])
        self.pg_graph.remove_node(remove)
        return (keep, conflict)


    def construct_pg_edge(self, prev_pg_id, cur_pg_id, genome_id, sequence_id):
        edge_data=self.pg_graph.get_edge_data(prev_pg_id, cur_pg_id, default=None)
        if edge_data != None:
            edge_data["genomes"].add(genome_id)
            edge_data["sequences"].add(sequence_id)
        else:
            self.pg_graph.add_edge(prev_pg_id, cur_pg_id, genomes=set([genome_id]), sequences=set([sequence_id]))

    def insert_feature(self, cur_pg_id, new_feature):
        #emit_extra=False
        genome_id=self.feature_index[new_feature].genome_id
        sequence_id=self.feature_index[new_feature].contig_id
        if not genome_id in self.pg_graph.node[cur_pg_id]['features']:
            self.pg_graph.node[cur_pg_id]['features'][genome_id]={sequence_id:[new_feature]}
            insert_level="genome"
        elif not sequence_id in self.pg_graph.node[cur_pg_id]['features'][genome_id]:
            self.pg_graph.node[cur_pg_id]['features'][genome_id][sequence_id]=[new_feature]
            insert_level="contig"
        else:
            #if self.context!="feature":
                #for cf in self.pg_graph.node[cur_pg_id]['features'][genome_id][sequence_id]:
                #    #if the distance is < k it is a special case of an 'extra character loop' which requires emitting an extra pg-node
                #    if abs(cf-new_feature)< self.ksize:
                #        emit_extra=True
                #        insert_level=self.context #so there won't be a problem
                #        break
            #if not emit_extra:
                self.pg_graph.node[cur_pg_id]['features'][genome_id][sequence_id].append(new_feature)
                insert_level="feature"
        return (insert_level)

    def detect_conflict(self, new_feature, cur_pg_id):
        genome_id=self.feature_index[new_feature].genome_id
        sequence_id=self.feature_index[new_feature].contig_id
        insert_level=None
        conflict=False
        end_fragments=False
        split=False
        if cur_pg_id == self.feature_index[new_feature].pg_assignment:
            return split, conflict, end_fragments
        if not genome_id in self.pg_graph.node[cur_pg_id]['features']:
            insert_level="genome"
        elif not sequence_id in self.pg_graph.node[cur_pg_id]['features'][genome_id]:
            insert_level="contig"
        else:
            insert_level="feature"
            split =True
        if self.context_levels[insert_level] > self.context_levels[self.context]:
            conflict=True
            cf = self.pg_graph.node[cur_pg_id]['features'][genome_id].values()[0][0]
            print "conflict between "+str(new_feature)+" and "+str(cf)+" in "+str(cur_pg_id)
            #determine if it is conflict class 1
            if self.context == "genome" and insert_level=="contig":
                nf_end=False
                cf_end=False
                i=0
                #check if both conflicting features are within k-distance from the end of a contig
                while i < self.ksize:
                    if self.feature_index[new_feature+i].rf_forward==None or self.feature_index[new_feature+i].rf_reverse==None:
                        nf_end=True
                    elif self.feature_index[new_feature-i].rf_forward==None or self.feature_index[new_feature-i].rf_reverse==None:
                        nf_end=True
                    if self.feature_index[cf+i].rf_forward==None or self.feature_index[cf+i].rf_reverse==None:
                        cf_end=True
                    elif self.feature_index[cf-i].rf_forward==None or self.feature_index[cf-i].rf_reverse==None:
                        cf_end=True
                    if nf_end and cf_end:
                        end_fragments=True
                        break
                    i+=1
            #if self.context!="feature" and \
            #genome_id in self.pg_graph.node[cur_pg_id]['features'] and \
            #sequence_id in self.pg_graph.node[cur_pg_id]['features'][genome_id]:
            #    for cf in self.pg_graph.node[cur_pg_id]['features'][genome_id][sequence_id]:
                    #if the distance is < k it is a special case of an 'extra character loop' which requires emitting an extra pg-node
            #        if abs(cf-new_feature)< self.ksize:
            #            split=True
            
        return split, conflict, end_fragments

    def detect_split(self, cur_pg_id, new_feature):
        genome_id=self.feature_index[new_feature].genome_id
        sequence_id=self.feature_index[new_feature].contig_id
        if self.context!="feature" and \
        genome_id in self.pg_graph.node[cur_pg_id]['features'] and \
        sequence_id in self.pg_graph.node[cur_pg_id]['features'][genome_id]:
            for cf in self.pg_graph.node[cur_pg_id]['features'][genome_id][sequence_id]:
                #if the distance is < k it is a special case of an 'extra character loop' which requires emitting an extra pg-node
                if abs(cf-new_feature)< self.ksize:
                    return True

    def move_features(self, target_pg_id, target_features):
        tn = self.pg_graph.node[target_pg_id]
        for t in target_features:
            rn_id = self.feature_index[t].pg_assignment
            if target_pg_id != rn_id:
                rn = self.pg_graph.node[rn_id]
                #print "moving feature"
                genome = self.feature_index[t].genome_id
                contig = self.feature_index[t].contig_id
                rn["features"][genome][contig].remove(t)
                if len(rn["features"][genome][contig]) == 0:
                    del rn["features"][genome][contig]
                if len(rn["features"][genome]) == 0:
                    del rn["features"][genome]
                tn["features"].setdefault(genome,{}).setdefault(contig,[]).append(t)
                self.feature_index[t].pg_assignment = target_pg_id


    def get_pg_id(self):
        return self.num_pg_nodes;

    def assign_pg_node(self, prev_feature, new_feature, pg_node):
        cur_pg_id=None
        #determine if there is a conflict based on mixed bundling
        conflict=False
        genome_id=self.feature_index[new_feature].genome_id
        sequence_id=self.feature_index[new_feature].contig_id
        #if there is already a node for this feature
        if self.feature_index[new_feature].pg_assignment != None:
            if (pg_node != None):
                if pg_node != self.feature_index[new_feature].pg_assignment:
                    cur_pg_id = self.feature_index[new_feature].pg_assignment
                    #cur_pg_id, conflict =self.merge_pg_node(self.feature_index[guide].pg_assignment, self.feature_index[new_feature].pg_assignment)
                #else they are EQUAL
                else:
                    cur_pg_id = pg_node
            #else there is no guide but this feature is already assigned
            else:
                cur_pg_id = self.feature_index[new_feature].pg_assignment

        #if this feature has yet to be assigned to a pg-node
        else:
            if pg_node!=None:
                #REPLACED insert_level=None
                cur_pg_id=pg_node
                #REPLACED insert_level = self.insert_feature(cur_pg_id, new_feature)
                #REPLACED if self.context_levels[insert_level] > self.context_levels[self.context]:
                #REPLACED    conflict=True

            else:
                cur_pg_id=self.num_pg_nodes
                self.num_pg_nodes+=1
                self.pg_graph.add_node(cur_pg_id, label=str(self.feature_index[new_feature].group_num), features={genome_id:{sequence_id:[new_feature]}})
            self.feature_index[new_feature].pg_assignment=cur_pg_id
        
        #REMOVE THIS
        #anomolous=set([2893])#, 1639, 1636, 3503, 2943, 3524, 3521, 1179, 1176])
        #if cur_pg_id in anomolous:
        #    print "hmmm"
        #END REMOVE

        if prev_feature != None:
            prev_pg_id = self.feature_index[prev_feature].pg_assignment
            self.construct_pg_edge(prev_pg_id, cur_pg_id, genome_id, sequence_id)
            if conflict:
                tier2=self.conflicts.setdefault(cur_pg_id, {})
                conflicted=tier2.setdefault(prev_pg_id, set([]))
                conflicted.add(new_feature)
        return cur_pg_id, conflict


    def maxInstanceOverlap(self, instance_key, comparison_list):
        instance_parts=set(instance_key.split('.'))
        max_intersect = 0
        max_key = []
        for c in comparison_list:
            cur_parts=set(c.split('.'))
            cur_intersect = len(instance_parts.intersection(cur_parts))
            if cur_intersect > max_intersect:
                max_interset = cur_intersect
                max_key = [c]
            elif cur_intersect > 0 and cur_intersect == max_intersect:
                max_key.append(c)
        if len(max_key) == 0:
            return [None]
        else:
            return max_key



    class featureContext():
        def __init__(self, kmer_side, direction, rhs_feature, prev_feature, new_feature, leaving_feature, conflict, split):
            self.kmer_side=kmer_side
            self.direction=direction
            self.rhs_feature=rhs_feature
            self.prev_feature=prev_feature
            self.new_feature=new_feature
            self.leaving_feature=leaving_feature
            self.conflict=conflict
            self.split=split


    def findKConflicts(self, pre_assignments, conflict_assignments):
        i=self.ksize-1
        while i > 0:
            for pg in pre_assignments[i]:
                for ik, features in pre_assignments[i][pg]["features"].iteritems():
                    for f in features:
                        #think about how the conflict should be represented overall
                        #look at the data structure required for determineAssignments
                        #start with that.
                        shift, conflict, end_fragments = self.detect_conflict(f, pg)
                        if shift:
                            conflict_assignments[i]["shift"].setdefault(pg, {'features':{}})
                            conflict_assignments[i]["shift"][pg]["features"].setdefault(ik,set([])).update(features)
                        elif end_fragments:
                            conflict_assignments[i]["c2_conflict"].setdefault(pg, {'features':{}})
                            conflict_assignments[i]["c2_conflict"][pg]["features"].setdefault(ik,set([])).update(features)
                        else:
                            conflict_assignments[i]["c1_conflict"].setdefault(pg, {'features':{}})
                            conflict_assignments[i]["c1_conflict"][pg]["features"].setdefault(ik,set([])).update(features)
            i-=1            


    #considering available guides/pg-nodes and conflicts to determine round of assignments
    #put new node assignments under None
    def determineAssignments(self, feature_pile, pre_assignments, conflicts=None):
        if conflicts == None:
            i=self.ksize-1
            while i >= 0:
                for ik, features in feature_pile[i]["by_instance"].iteritems():
                    max_keys = self.maxInstanceOverlap(ik, pre_assignments[i]["inst_key_map"].keys())
                    if max_keys[0] == None: #no guides exist
                        #create new node and put it in
                        tmp_id = "nn"+str(len(pre_assignments[i]["new_nodes"].keys()))
                        assign_block="new_nodes"
                        new_status=True
                    elif len(max_keys) == 1:
                        tmp_id = pre_assignments[i]["inst_key_map"][max_keys[0]]
                        new_status = (type(tmp_id) == str and tmp_id.startswith('nn'))
                        assign_block = "new_nodes" if new_status else "assignments"
                    elif len(max_keys) > 1:
                        #pick the one with the most features. this should be rare.
                        sys.stderr.write("Tie for pg-node selection based on instance keys "+" ".join(max_keys)+"\n")
                        max_features =0
                        tmp_id = None
                        mk = None
                        for k in max_keys:
                            pg_node = pre_assignments[i]["inst_key_map"][k]
                            num_features=self.num_features_pg_node(pg_node)
                            if num_features > max_features:
                                max_features=num_features
                                tmp_id = pg_node
                        new_status = (type(tmp_id) == str and tmp_id.startswith('nn'))
                        assign_block = "new_nodes" if new_status else "assignments"
                    if not ik in pre_assignments[i]["inst_key_map"]:
                        #track guide creates new_node position in pre_assignments
                        guide_feature=next(iter(features))
                        self.track_guide(pre_assignments, guide_feature.new_feature, i, negative=guide_feature.direction, target_pg=tmp_id, new_node=new_status)
                    pre_assignments[i][assign_block][tmp_id]["features"].setdefault(ik,set([])).update(features)
                i-=1
        else:
            i=self.ksize-1
            while i >= 0:
                #if there is a shift conflict create a new node. only need the pg_node and the offending pre_assignment instance key
                if "shift" in conflicts[i] and len(conflicts[i]["shift"]) > 0:
                    for pg in conflicts[i]["shift"]["features"]:
                        for ik in conflicts[i]["shift"]["features"][pg]:
                            tmp_id = len(pre_assignments[i]["new_nodes"].keys())
                            pre_assignments[i]["new_nodes"].setdefault(tmp_id, {'guides':{},'features':{}})
                            pre_assignments[i]["new_nodes"][tmp_id]["features"][ik]=conflicts[i]["shift"][pg]["features"][ik]
                            #take out features of pre_assignments
                            pre_assignments[i]["assignments"][pg]["features"].pop(ik)

    
    #store alternate pg-nodes at a given k-position
    def storeAlternates(self, pg_ids):
        assert type(pg_ids == set)
        for i in pg_ids:
            self.pg_node_alt.setdefault(i, set([])).update(pg_ids)

    #make assignments that are in pre_assignments
    def makeAssignments(self, pre_assignments):
        i=self.ksize-1
        while i > 0:
            new_nodes= set(pre_assignments[i]["assignments"].keys()+pre_assignments[i]["new_nodes"].keys())
            if len(new_nodes) > 1: storeAlternates(new_nodes)
            for pg in pre_assignments[i]["assignments"]:
                for ik, features in pre_assignments[i]["assignments"][pg]["features"].iteritems():
                    for f in features:
                        self.assign_pg_node(prev_feature=f.prev_feature, new_feature=f.new_feature, pg_node=pg)
            for pg in pre_assignments[i]["new_nodes"]:
                for ik, features in pre_assignments[i]["new_nodes"][pg]["features"].iteritems():
                    for f in features:
                        self.assign_pg_node(prev_feature=f.prev_feature, new_feature=f.new_feature, pg_node=None)
            i-=1
                        
                        

    def oldMethod(self):
        #PHASE 2:
        #Determine the best guide to use
        #Reorganize by target pg_node so that conflicts and splits can be detected per
        guide_to_assign={}
        for feature_tuple in to_assign:
            #unpack
            kmer_side, direction, rhs_feature, prev_feature, new_feature, leaving_feature, conflict, split= feature_tuple
            #determine the best guide to use
            if len(pg_set) > 1:
                instance_key = self.feature_index[new_feature].instance_key
                if instance_key in target_guides:
                    target_guide = target_guides[instance_key][0]
                else:
                    max_key = self.maxInstanceOverlap(instance_key, target_guides.keys())
                    if max_key != None:
                        target_guide = target_guides[max_key][0]
                    else:
                        target_guide = None
            elif len(target_guides.keys()) >= 1: #all instance keys refer to the same pg-node
                target_guide = target_guides.values()[0][0]
            else:
                target_guide = None
            if target_guide != None:
                pg = self.feature_index[target_guide].pg_assignment
                guide_to_assign.setdefault(pg,(target_guide, []))[1].append(feature_tuple)
            else:
                guide_to_assign.setdefault(None,(target_guide, []))[1].append(feature_tuple)

                #if (not cur_node.palindrome) and (kmer_side != rhs_guide_side):
                    #these should always be the same except for palindromes
                #    assert LogicError

            #PHASE 4:Make assignments based on the previous two phases
            for pg_node, assign_tuple in guide_to_assign.iteritems():
                #PHASE 3: Determine if there are ANY conflicts.
                target_guide, assign_list = assign_tuple
                num_conflict=num_split=0
                c1_conflict =False
                break_here=False # if break conflicts is true then want to use unassigned kmer as guide so that incoming features will be assigned to a new node
                if target_guide != None:
                    #conflict occurs when mixed bundling tries to violate synteny context
                    num_split, num_conflict, c1_conflict, assign_list = self.find_conflicts(assign_list, target_guide)
                if num_split > 0:
                    print "SPLIT CONFLICT: rf-node "+str(cur_node.nodeID)+" there are "+str(num_split)+" splits in a bundle of size "+str(num_targets)
                elif num_conflict > 0:
                    if c1_conflict:
                        print "C1 CONFLICT: rf-node "+str(cur_node.nodeID)+" there are "+str(num_conflict)+" conflicts in a bundle of size "+str(num_targets)
                    else:
                        if self.break_conflict:
                            break_here=True
                            target_guide=None
                        print "C2 CONFLICT: in rf-node "+str(cur_node.nodeID)+" there are "+str(num_conflict)+" conflicts in a bundle of size "+str(num_targets)
                default_guide = default_guide_cat = default_guide_side = None
                split_guide = split_guide_cat = split_guide_side = None
                for cur_tuple in assign_list:
                    #unpack
                    kmer_side, direction, rhs_feature, prev_feature, new_feature, leaving_feature, conflict, split= cur_tuple
                    cur_node.assigned_features[direction].add(rhs_feature)
                    
                    assigned =  self.feature_index[new_feature].pg_assignment != None
                    do_queue=True
                    #DON't DO THIS. Always queue.
                    #if break_here and assigned:
                    #    prev_feature = None
                    #    do_queue=False

                    #assign to pg-node
                    #conflict=False
                    if target_guide == None:
                        #if (not assigned) or (not break_here): #only care if it assigned or not if break_here aka conflict
                        if default_guide == None:
                            default_guide=new_feature
                            self.assign_pg_node(prev_feature=prev_feature, new_feature=new_feature, guide=None)
                        else:
                            new_guide=default_guide
                            assignment,conflict=self.assign_pg_node(prev_feature=prev_feature, new_feature=new_feature, guide=new_guide)

                    else:
                        if split:
                            if split_guide == None:
                                split_guide=new_feature
                                #IS THERE ANY GURANTEE THAT THE NODE YOU PICK WON'T BE A SPLIT CONFLICT AGAIN?
                                assignment, conflict = self.assign_pg_node(prev_feature=prev_feature, new_feature=new_feature, guide=None)
                            else:
                                assignment,conflict=self.assign_pg_node(prev_feature=prev_feature, new_feature=new_feature, guide=split_guide)
                        else:
                            assignment,conflict=self.assign_pg_node(prev_feature=prev_feature, new_feature=new_feature, guide=target_guide)
                            if up_targets and conflict:
                                print "conflict with up-targets! in pg-node "+str(assignment)+" from rf-node "+str(cur_node.nodeID)
                            elif conflict:
                                print "conflict in pg-node "+str(assignment)+" from rf-node "+str(cur_node.nodeID)
                    #make sure that where ever the feature is assigned, that all the other features in this visit, with the same instance key are there too.
                    instance_key = self.feature_index[new_feature].instance_key
                    #if instance_key in target_guides:
                    #    self.move_features(assignment, target_guides[instance_key])
                    #queue base on leaving feature
                    #if do_queue:

    #walk down the kmer and fill out to_assign
    #process features remaining in cur_node.features
    def kfill_feature_pile(self, feature_pile, pre_assignments, cur_node, prev_node, node_queue, node_bundles):
            i=0
            while i < self.ksize: #because any features remaining represent "new" threads need to assign the entire k-mer
                for direction in self.target_cat:
                    for rhs_feature in cur_node.features[direction]:
                        new_feature_adj= i
                        prev_feature_adj=i-1
                        #TODO Consider adjusting this based on signed features.
                        #TODO Consider what would be necessary to adapt this to consistency use for palindrome
                        if not direction:
                            new_feature_adj=new_feature_adj*-1
                            prev_feature_adj=prev_feature_adj*-1
                        new_feature=rhs_feature+new_feature_adj
                        instance_key = self.feature_index[new_feature].instance_key
                        prev_feature=rhs_feature+prev_feature_adj
                        if self.feature_index[new_feature].pg_assignment != None:
                            self.track_guide(pre_assignments, new_feature, self.ksize-(i+1), negative=direction)
                            #track feature to "add" in case of anchor instance key expansion AND/OR pg-edge addition
                            feature_info=self.featureContext(self.ksize-(i+1), direction, rhs_feature, prev_feature, new_feature, leaving_feature=None, conflict=False, split=False)
                            self.track_feature(self.ksize-(i+1), feature_pile, new_feature, feature_info, instance_key)
                        else:
                            #GET FEATURE_INFO HERE??????
                            feature_info=self.featureContext(self.ksize-(i+1), direction, rhs_feature, prev_feature, new_feature, leaving_feature=None, conflict=False, split=False)
                            self.track_feature(self.ksize-(i+1), feature_pile, new_feature, feature_info, instance_key)
                                    
                        cur_node.assigned_features[direction].add(rhs_feature) #only track rhs. IS THIS NECESSARY?
                        #there are two cases where the feature could be about to leave the kmer-frame. If they are on the left or right of the kmer
                        if i == 0:
                            #this feature is on the rhs of kmer
                            self.queueFeature(cur_node, 1, direction, new_feature, node_queue, node_bundles, up_node=prev_node)
                            prev_feature=None
                        if i == self.ksize-1:
                            #this feature is on the lhs of kmer
                            self.queueFeature(cur_node, 0, direction, new_feature, node_queue, node_bundles, up_node=prev_node)
                i+=1

    def getFeatureRecord(self, feature):
        return self.feature_index[abs(feature)]

    #for all features in all k positions. check if the instance key for the feature is an anchor instance key.
    #if so expand it through the whole kmer. this shouldn't be used for palindromes since threads here can lack exact orientation relative to the bundle.
    #the only time this can/needs to be used is when it is a duplicate node since anchor nodes will pick up all threads
    #ACTUALLY palindromes should be fine because the anchor instance key expansion will place it in the position it needs to be in.
    #if its a starting node then its an anchor node (no need to do middle, or lhs)
    #if its a palindrome can't do this
    #if its a duplicate node then the previous k-1 features will have been checked.
    #NOTE anchor instance keys are really just a move-ahead proceedure since eventually when that anchor rf-node is encountered
    #it will pick up all combine all threads and walk back the ones that haven't been added.
    #instance keys themselves will still be used to segregate assignments when multiple pg-nodes are available
    #also the concept of anchor keys is still valid since that guarantee is held up when the anchor rf-node is encountered.
    def anchorInstanceExpansion(self, feature_pile, cur_node):
        if cur_node.palindrome or (not cur_node.duplicate):
            return False #nothing expanded
        #expand projection/target selection based on instance key 
        instances_pkg = self.anchor_instance_keys.get(self.feature_index[nxt_feature].instance_key, [None,[]])
        

    
    def getColumnGuides(self, pre_assignments, col):
        result_guides =set([])
        for pg_node in pre_assignments[col]["assignments"]:
            for instance_key in pre_assignments[col]["assignments"][pg_node]:
                for instance_key in pre_assignments[col]["assignments"][pg_node]["guides"]:
                    result_guides.update(pre_assignments[col]["assignments"][pg_node]["guides"][instance_key])
        return result_guides


    #walk down the kmer and fill out guides based on those present on LHS or RHS
    #features are signed in pre_assignments
    def fill_guides(self, pre_assignments):
        lhs =self.getColumnGuides(pre_assignments, 0)
        rhs =self.getColumnGuides(pre_assignments, self.ksize-1)
        count = 1
        while count < self.ksize:
            for sf in lhs:
                cur_sf = sf + count
                self.track_guide(pre_assignments, abs(cur_sf), count, negative=(cur_sf < 0)) 
            for sf in rhs:
                cur_sf = sf - count
                self.track_guide(pre_assignments, abs(cur_sf), (self.ksize-1)-count, negative=(cur_sf < 0))
            count+=1


    #keep track of which nodes are available for assignment in this 'column' of the kmer
    def track_guide(self, pre_assignments, unsigned_guide, col, negative, target_pg=None, new_node=False):
        new_guide = self.sign_feature(unsigned_guide, negative)
        guide_ikey=self.getFeatureRecord(new_guide).instance_key
        if target_pg == None:
            target_pg =getFeatureRecord(new_guide).pg_assignment
            if target_pg == None:
                assert LogicError
        if not new_node:
            pre_assignments[col]["assignments"].setdefault(target_pg, {'guides':{},'features':{}})
            pre_assignments[col]["assignments"][target_pg]['guides'].setdefault(guide_ikey, set([])).add(new_guide)
            pre_assignments[col]["inst_key_map"].setdefault(guide_ikey, target_pg)
        else:
            pre_assignments[col]["new_nodes"].setdefault(target_pg, {'guides':{},'features':{}})
            pre_assignments[col]["new_nodes"][target_pg]['guides'].setdefault(guide_ikey, set([])).add(new_guide)
            pre_assignments[col]["inst_key_map"].setdefault(guide_ikey, target_pg)
            

    def track_feature(self, col, feature_pile, new_feature, feature_info, instance_key):
        if not new_feature in feature_pile[col]['all_features']:
            feature_pile[col]['by_instance'].setdefault(instance_key, set([])).add(feature_info)
            feature_pile[col]['all_features'].add(new_feature)

    #right now kmer_side is 0 to indicate left and 1 to indicate right
    def side_to_kcoord(self, kmer_side):
        return kmer_side*self.ksize

    #make a feature negative if it is present in a decreasing L to R series in this kmer
    def sign_feature(self, feature, negative):
        if negative:
            return feature * -1
        else:
            return feature


    def expand_features(self, prev_node, cur_node, targets, guide_array, node_queue, node_bundles,up_targets=False):
        pre_assignments = [dict(new_nodes=dict(), inst_key_map=dict(), assignments=dict()) for x in range(self.ksize)] # used to track which features are to be assigned to which pg-nodes. k dictionaries {pg_node_id: {"guides":{instance_key:set([feature_id])}, "features":{instance_key:set([feature_id])}}}
        conflict_assignments = [dict(shift=dict(), c1_conflict=dict(), c2_conflict=dict()) for x in range(self.ksize)] # used to track which features being assigned have conflicts
        feature_pile = [dict(by_instance=dict(), all_features=set([])) for x in range(self.ksize)] # which features to assign, organized by instance_key. k dictionaries {"by_instance": {instance_key:set([feature_id])}, "all_features":set([feature_id])}

        if self.debug:
            self.visit_number+=1
            existing_label = self.rf_graph.node[cur_node.nodeID]["visit"]
            self.rf_graph.node[cur_node.nodeID]["visit"] =  ",".join([str(self.visit_number),existing_label]) if len(existing_label) else str(self.visit_number)
        sys.stderr.write("number of pg-nodes is "+str(self.pg_graph.number_of_nodes())+"\n")
        num_targets=0
        if (targets!=None):
            num_targets=len(targets[0][0])+len(targets[0][1])+len(targets[1][0])+len(targets[1][1])
        rhs_guide=rhs_guide_cat=None
        self.target_cat=[0,1]
        target_guides={} #store available guides under their instance_key
        pg_set=set([]) # all the different pg_nodes these features could be assigned to
        
        
        #PHASE 1: Establish if any features have already been assigned, so that they can be used as guides for pg-assignment
        #TODO from the revist guide populate the to_assign structure
        #NOTE that the revist guide itself needs to be restructured to include the multiple choices that may have been present or created on the last visit
        #Assuming its possible for that to happen in a palindrome or duplicate node 
        #for guide_array in revisit_guides: #a guide is incoming when pass up new info to a non-anchor nonde via DFS ascending. and need to pass the information DOWN to a node that has already been visited
        if guide_array != None:
            rhs_guide=guide_array[0] #NOTE this could be combined with incoming guide parameter (maybe) since it will need a similar structure
            new_guide_cat=guide_array[1] #used if there are new things in this anchor node
            new_guide_side=guide_array[2]#only really needed if this is a palindrome THIS NEEDS TO BE RENAMED. Its actually the side of the kmer that the guide is good for.
            #this is a bit complicated. if the kmer is a palindrome you can't be sure what side of the kmer the feature will be on except from the previous kmer(s) and the features passed in 
            #from it (which is why a palindrome isn't an anchor node). In each kmer's data structure only the RHS feature is recorded. A revist guide is only used when this node has 
            #already been visited.
            #What happens next is the incoming guide is converted to the RHS of this kmer based on the information about what the edge type implies about the direction/side
            if rhs_guide != None:
                #TODO find the correct location to flag if there are multiple guides necessary in a duplicate or palindrome node
                #TODO bring the storage of guides in here from the parent function since need to have multiple
                #TODO NOTE can use the visit numbers as references to array of guides. though that supposes the current system isn't good enough
                rhs_adj_info=self.rhs_adj_table[rhs_guide_cat][new_guide_side]
                #the new guide is signed so can do simple increment/decrement depending on LHS or RHS. always take abs() to get actual ID
                new_guide= rhs_guide + rhs_adj_info['new_feature_adj'] #guides are aligned to the right side and walked left from there
                #REPLACED instance_key = self.feature_index[new_guide].instance_key
                self.track_guide(pre_assignments, new_guide, self.side_to_kcoord(new_guide_side), negative=new_guide_cat)
                # SEE IF YOU CAN GET AWAY WITHOUT DOING THIS TILL YOU NEED IT self.fill_guides(pre_assignments, rhs_guide, rhs_adj_info['new_feature_adj'])
                #REPLACED pg_set.add(self.feature_index[new_guide].pg_assignment)
                #REPLACED target_guides[instance_key]=[new_guide]

                        
        #whether this is an anchor or not there will be targets passed down if it is not the start of a traversal.
        if (num_targets>0):
            # if there are targets then this isn't the first node visited
            # this means only one new column aka 'character' in the kmer needs to be expanded 
            # no matter which feature is being 'targeted' only the RHS of this kmer exists in the features structure.
            #targets organized as targets["left" & "right" == 0 & 1][ "increasing" & "decreasing" == 0 & 1 ]
            #REPLACED currently_seen=set([])

            #REPLACED available_nodes=set([])
            #PROCESS TARGETS
            for kmer_side in target_cat:
                for direction in target_cat:
                    rhs_adj_info=self.rhs_adj_table[direction][kmer_side]
                    #get the targets that have been assigned
                    for rhs_feature in targets[kmer_side][direction]:
                        new_feature=rhs_feature + rhs_adj_info['new_feature_adj']
                        prev_feature= rhs_feature + rhs_adj_info['prev_feature_adj']
                        leaving_feature=rhs_feature + rhs_adj_info['leaving_feature_adj']
                        cur_instance_key = self.feature_index[new_feature].instance_key
                        #get the features with matching instance keys that can be co-assigned
                        instance_info = self.anchor_instance_keys.get(cur_instance_key, [None,[]])
                        if not rhs_feature in cur_node.features[direction]:
                            if not rhs_feature in cur_node.assigned_features[direction]:
                                print "missing projected "+str(rhs_feature)+" in "+str(cur_node.nodeID)+" from "+str(prev_node.nodeID)
                                assert LogicError
                            #NO PG-EDGE construction until node assignment finalized
                            #else:
                                #construct pg-edge to previously created node
                            #    self.construct_pg_edge(self.feature_index[prev_feature].pg_assignment, self.feature_index[new_feature].pg_assignment, self.feature_index[new_feature].genome_id, self.feature_index[new_feature].contig_id)
                        else:
                            #this initial loop through the targets is really just to see if any have already been assigned
                            if self.feature_index[new_feature].pg_assignment != None:
                                #REPLACED rhs_guide = rhs_feature
                                #REPLACED rhs_guide_cat=direction
                                #REPLACED rhs_guide_side=kmer_side #need this for palindromes
                                #REPLACED pg_set.add(self.feature_index[new_feature].pg_assignment)
                                #REPLACED target_guides.setdefault(instance_key,[]).append(new_feature)
                                self.track_guide(pre_assignments, new_feature, self.side_to_kcord(kmer_side), negative=direction)
                            else:
                                cur_node.features[direction].remove(rhs_feature)
                                cur_info=self.featureContext(kmer_side, direction, rhs_feature, prev_feature, new_feature, leaving_feature, False, False)
                                self.track_feature(self.side_to_kcord(kmer_side), feature_pile, new_feature, cur_info, cur_instance_key)
                        if up_targets:#if up_targets is true then these features were returned from a DFS exploration of an anchor node and passed here as target.
                            #when queueing based on return 'new' features make sure don't do a DFS "up"
                            q_rfid = self.queueFeature(cur_node, (not kmer_side), direction, leaving_feature, node_queue, node_bundles, up_node=prev_node) #no prevent_node
                        else:
                            q_rfid = self.queueFeature(cur_node, (not kmer_side), direction, leaving_feature, node_queue, node_bundles) #no prevent_node
                    #in the case of an anchor node also add non-target features as potential guides
                    #but to be used here they must be on the same side/direction as incoming targets
                    #REPLACED if cur_node.anchorNode() and len(targets[kmer_side][direction]):
                    #REPLACED    for rhs_feature in cur_node.features[direction]:
                    #REPLACED        potential_guide = rhs_feature + rhs_adj_info['new_feature_adj']
                    #REPLACED        if self.feature_index[potential_guide].pg_assignment != None:
                                #REPLACED pg_set.add(self.feature_index[potential_guide].pg_assignment)
                                #REPLACED instance_key= self.feature_index[potential_guide].instance_key
                                #REPLACED target_guides.setdefault(instance_key,[]).append(potential_guide)
                    #REPLACED            self.track_guide(pre_assignments, potential_guide, self.side_to_kcord(kmer_side))
                                
            #if there is nothing to assign leave
            #REPLACED num_features= sum(len(i["all_features"]) for i in to_assign)
            #REPLACED if num_features == 0:
            #REPLACED    return

            #PHASE 2:
            #Determine the best guide to use


        #TODO figure out consistentcy palindrome procedure. Flip LHS to RHS. Has to take place w/ targets.
        #TODO look at anchor instance keys. AND other instance keys. Do cardinality of instance keys

        #the only place where all k-features are added is at an anchor node if there are features that are not currently part of a bundle
        #if a "feature thread" is already in the bundle then the previous k-1 features have been processed
        #in that case the previous k-1 features should have already been checked for anchor instance expansion
        if cur_node.anchorNode():
            #if this is an anchor node and a starting node then everything needs to expanded/assigned to a pg-node
            #if this is an anchor node and had targets incoming then everything remaining is new and needs to be fully expanded
            #at this point anything remaining is regarded as 'new' and can be passed as targets up or down !!!!!
            #rhs_guide is used to track a feature "thread" that has already been assigned so that current features can be assigned to the correct pg_node
            self.kfill_feature_pile(feature_pile, pre_assignments, cur_node, prev_node, node_queue, node_bundles)
            self.fill_guides(pre_assignments)
            cur_node.features[0]=set([])#after assigning all features clear it out.
            cur_node.features[1]=set([])#after assigning all features clear it out.

        #anchor_expanded=False
        #anchor_expanded = self.anchorInstanceExpansion(feature_pile)

        self.determineAssignments(feature_pile, pre_assignments, conflicts=None)
        self.findKConflicts(pre_assignments, conflict_assignments)
        self.determineAssignments(feature_pile, pre_assignments, conflict_assignments)
        self.makeAssignments(pre_assignments)
        #self.storeGuides(pre_assignments)
        #end expand_features






    class VisitPack():
        def __init__(self, prev_node, cur_node, targets, guide):
            self.prev_node=prev_node
            self.cur_node=cur_node
            self.targets=targets
            self.guide=guide
            self.node_bundles={}#organized by rf-node id, values are next features "bundle" to look for (in case of palindrome or duplicate)
            self.node_queue=deque()#tuples of (next rf-node id to visit, the guide to send to it, and the next features to look for)
            self.new_targets=deque()
            self.visited=False
    
    #get a guide from a target bundle
    #this is used when new targets are returned from DFS/TFS.
    #Because of this the new target will be on the opposite side than the previous targets
    def getTargetGuide(self,targets):
        kmer_side=0
        guide=None
        while kmer_side < len(targets):
            if guide != None:
                break
            direction=0
            while direction < len(targets[kmer_side]):
                if len(targets[kmer_side][direction]) > 0:
                    guide= (iter(targets[kmer_side][direction]).next(),direction, not kmer_side) #can be any feature just assigned.
                    break
                direction+=1
            kmer_side+=1
        if guide == None:
            assert LogicError
        return guide



    #non-recursive version of tfs_expand
    def tfs_expand_nr(self, prev_node, cur_node,targets, guide):
        tfs_stack=[]
        tfs_stack.append(self.VisitPack(prev_node,cur_node,targets,guide))
        pv=None #pv previous visit
        while len(tfs_stack):
            cv=tfs_stack.pop() #cv current visit
            if not cv.visited:
                self.expand_features(cv.prev_node, cv.cur_node, cv.targets, cv.guide, cv.node_queue, cv.node_bundles)
                cv.visited=True
            next_visit=None # this is temp
            if len(cv.node_queue):
                next_node_id, next_guide= cv.node_queue.popleft()
                next_node=self.rf_node_index[next_node_id]
                next_targets=cv.node_bundles[next_node_id]
                cv.node_bundles[next_node_id]=[[set([]),set([])],[set([]),set([])]]
                next_visit=self.VisitPack(prev_node=cv.cur_node, cur_node=next_node, targets=next_targets, guide=next_guide)
                if next_node == cur_node.nodeID: #self loop special case. use same bundles
                    next_visit.node_bundles=cv.node_bundles
            if len(cv.node_queue) == 0 and next_visit==None:
                if pv!= None and -1 in cv.node_bundles: # -1 is used to track 'new' threads exposed by anchor node that are to be passed 'up'
                    pv.new_targets.append(cv.node_bundles[-1])
            else:
                tfs_stack.append(cv)
            if next_visit != None:
                tfs_stack.append(next_visit)
            if not cv.cur_node.anchorNode():
                while len(cv.new_targets):
                    cntargets=cv.new_targets.popleft()
                    new_guide= self.getTargetGuide(cv.targets)#can be any feature just assigned.
                    #after this or during this...need to think about the forking guide problem wrt restoring things into the queue
                    #if there are return targets and a guide for this node...it means a guide needs to be projected to go with all those nodes that have already been visited by TFS
                    #so if there is a guide: 
                    self.expand_features(cv.prev_node, cv.cur_node, targets=cntargets, guide=new_guide, node_queue=cv.node_queue, node_bundles=cv.node_bundles, up_targets=True)
            pv=cv

    #recursive function with a visit queue.
    def tfs_expand(self, prev_node, cur_node, targets, guide):
        node_bundles={}#organized by rf-node id, values are next features "bundle" to look for (in case of palindrome or duplicate)
        prev_bundles={}
        node_queue=deque()#tuples of (next rf-node id to visit, the guide to send to it, and the next features to look for)
        #seen_targets=set([])
        #put this conditional in expand_features
        #if cur_node.anchorNode():
            #expand_features assigns features to pg-nodes, queues rf-nodes for visiting, and organizes feature threads to pass to each
            #self.expand_features(cur_node, targets=None, guide=None, node_queue=node_queue)
            #figure out targets not seen before
        self.expand_features(prev_node, cur_node, targets=targets, guide=guide, node_queue=node_queue, node_bundles=node_bundles)
        #seen_targets.update(targets)
        #NOW need to think carefully about where targets will be set. I guess always right side? They need be divided into decreasing and increasing.
        #Also this means projection based on edge will need to be well thought out. 
        while len(node_queue):
            #next_guide needs to be populated IF AND WHEN new targets come back from the DFS AND 
            #they proceed to a Node that has already been visited WHICH surely has to be a 
            #duplicate/palindrome node else IT WOULD RETURN NEW TARGETS that would be passed down to the node from which NEW TARGETS comes
            next_node_id, next_guide= node_queue.popleft()
            #if prev_node and next_node_id != prev_node.nodeID: #prevent recursing "up"
            next_node=self.rf_node_index[next_node_id]
            next_targets=node_bundles[next_node_id]
            #careful here sets are passed by reference
            #EITHER the next rf-node is still in the queue in which case it is OK to accumulate bundle info from new_targets found in other visits
            #OR the bundle information is being transmitted right here in which case only the new stuff should be transmitted down next time
            #so clear out the bundle info here
            node_bundles[next_node_id]=[[set([]),set([])],[set([]),set([])]]
            #Recursion!
            new_targets = self.tfs_expand(prev_node=cur_node, cur_node=next_node, targets=next_targets, guide=next_guide)
            if not cur_node.anchorNode():
                #not needed because expand_features will put features in node bundles new_return_targets=return_targets.difference(seen_targets)
                #must process 
                #just got new targets returned from a DFS. expand them, and update queue based on them
                if(len(new_targets)):
                    new_guide= iter(targets).next() #can be any feature just assigned.
                    #after this or during this...need to think about the forking guide problem wrt restoring things into the queue
                    #if there are return targets and a guide for this node...it means a guide needs to be projected to go with all those nodes that have already been visited by TFS
                    #so if there is a guide: 
                    self.expand_features (prev_node, cur_node, targets=new_targets, guide=new_guide, node_queue=node_queue, node_bundles=node_bundles, up_targets=True)
                    #seen_targets.update(new_return_targets)
        if prev_node==None:
            return None
        elif -1 in node_bundles: # -1 is used to track 'new' threads exposed by anchor node that are to be passed 'up'
            return (node_bundles[-1])


        





        #transform the kmerNode graph (rf-graph) into a pg-graph
    #if the minOrg requirment is not met the node is added to the graph but is marked in active.
    #dfs still proceeds in case a node that does meet minOrg is encounterd (which will require considering prev. expanded nodes in identity resolution)
    def bfsExpand(self, minOrg):
        print "expanding kmer graph in to pg-graph total knodes: "+str(len(self.rf_node_list))
        for start_k_id, start_knode in enumerate(self.rf_node_list):
            if start_knode.visited:
                continue
            else:
                knode_q=deque()
                knode_q.append((start_k_id,None,None))
                prev_k_id=None
                in_edge_status=None # type of edge arrived by
                while len(knode_q) > 0:
                    visiting_k_id, prev_k_id, in_edge_status=knode_q.popleft()
                    cur_knode=self.rf_node_list[visiting_k_id]
                    #do work for expanding this kmer node into pg-graph nodes
                    #if prev_knode and incoming_status != None :
                    if prev_k_id != None:
                        cur_knode.visitNode(self.rf_node_list[prev_k_id], in_edge_status, self)#expand and store refs to pg-ndoes
                    else:
                        cur_knode.visitNode(None, None, self)
                    for k_id in cur_knode.linkOut:
                        if k_id == visiting_k_id:#self loop this should not happen because of kmer levels
                            continue
                            #something selfish
                            #cur_knode.self_edge=True
                        elif k_id == prev_k_id:#return loop
                            continue
                            #handle return loop. create single edge back and apply labels
                            #return_node=self.kmerList[k_id]
                            #return_node.applyInfo(self)
                        elif self.rf_node_list[k_id].visited or self.rf_node_list[k_id].queued:	
                            return_node=self.rf_node_list[k_id]
                            return_node.updateNode(cur_knode, cur_knode.linkOut[k_id], self)
                        else:
                            #if k_id ==208:
                            #	print "Debug: why is this being queued so much?"
                            knode_q.append((k_id, visiting_k_id, cur_knode.linkOut[k_id]))
                            self.rf_node_list[k_id].queued=True
                    cur_knode.addPGEdges(self)

# undirected weighted
class pFamGraph(nx.Graph):
    def __init__(self, storage, minOrg=2):
        #Graph.__init__(self, weighted=True)
        Graph.__init__(self)
        self.createGraph(storage, minOrg)
    def add_path_cumul_attr(self,nlist,**kwargs):
        edges=list(zip(nlist[:-1],nlist[1:]))#create list of edges
        edge_ids=[]
        for e in edges:
            if self.has_edge(*e):
                for k in kwargs:
                    if type(kwargs[k])==set:
                        try: self.adj[e[0]][e[1]][k] |= kwargs[k]#  union of attribute
                        except: 
                            try: self.adj[e[0]][e[1]][k]=kwargs[k].copy()
                            except: self.adj[e[0]][e[1]]=kwargs[k]
            else:
                kwargs['id']=str(self.number_of_edges())
                self.add_edge(e[0],e[1],kwargs)
                if kwargs['id']=="0":
                    warning("edge 0 is "+e[0]+" "+e[1])
            try: edge_ids.append(self.adj[e[0]][e[1]]['id'])
            except: warning("no ID for edge "+e[0]+" "+e[1])
        return edge_ids 

    #update the edge weight based on a designated attribute
    #also flatten to a string since writing list objects isn't supported
    #weight_attr has to be weight. label_attr = (what to get, and what to label it)
    #also setting ID so that it can be used in building map from sid to edge
    def update_edges(self, weight_attr='getOrganism', divisor=1, label_attr=('getReplicon','replicons'), remove_attrs=[]):
        edge_counter=itertools.count()
        for u,v,data in self.edges_iter(data=True):
            #try: self.adj[e[0]][e[1]][e_attr]=list(self.adj[e[0]][e[1]][e_attr])
            #except: pass
            data['label']=''
            weight_set=set()
            label_set=set()
            for i in data['instances']:
                weight_set.add(getattr(i,weight_attr))
                label_set.add(getattr(i,label_attr[0])())
            try: data['weight']=len(weight_set)/float(divisor)
            except:
                try:data['weight']=0
                except: pass
            if label_attr:
                try: data[label_attr[1]]=", ".join(list(label_set))
                except: pass
            for r in remove_attrs:
                try: data.pop(r,None)
                except: pass
            data['id']=next(edge_counter)
                
    def update_node_cumul_attr(self, n_id, **kwargs ):
        if n_id in self.node:
            for k in kwargs:
                try: self.node[n_id][k]=kwargs[k] | self.node[n_id][k]
                except:
                    try: self.node[n_id][k]=kwargs[k].copy()
                    except: print "cannot add attribute to node "+str(n_id)
    
    #calculate the node weight and change the set attributes to string
    #so that they can be written by graphml writer
    def update_node_attr_final(self, weight_func, family_func, divisor=1, remove_attrs=[], minOrg=2):
        remove_set=set()
        for n in self.nodes():
            weight_set=weight_func(n)
            node_summary=n.get_summary()
            if len(node_summary['organisms']) < minOrg:
                remove_set.add(n)
            try:
                self.node[n]['weight']=len(weight_set)/float(divisor)
                self.node[n]['id']=str(n.id)
                self.node[n]['familyID']=str(n.famID)
                self.node[n]['label']=family_func(n.famID)
                self.node[n]['locations']=','.join(list(node_summary['locations']))
                self.node[n]['organisms']=','.join(list(node_summary['organisms']))
                
            except: pass
            for r in remove_attrs:
                try: self.node[n].pop(r,None)
                except: pass
            for a in self.node[n]:
                if type(self.node[n][a])==set:
                    self.node[n][a] = ','.join(self.node[n][a])
        for n in remove_set:
            self.remove_node(n)

                                
                            
                    
                
                
    ##this function takes the storage class and constructs the graph from it
    def createGraph(self, storage, minOrg):
        num_orgs=len(storage.summaryLookup.keys())
        temp_size=len(storage.kmerLookup.keys())
        total_tax=len(storage.completeTaxSummary())
        for k in storage.replicon_map: storage.replicon_map[k]=list(storage.replicon_map[k])
        print " ".join(["starting",str(temp_size),str(total_tax),str(num_orgs)])
        for n in storage.pg_initial:
            if n != None:
                for e in n.edges:
                    n2=storage.getPGNode(e)
                    if n.subsumed or n2.subsumed:
                        sys.stderr.write("Logic Error: A node that should have been subsumed and removed is in the graph\n")
                        sys.exit()
                    self.add_edge(n.famSubset, n2.famSubset)
                    if not 'instances' in self[n.famSubset][n2.famSubset]:
                        self[n.famSubset][n2.famSubset]['instances']=set()
                    self[n.famSubset][n2.famSubset]['instances'].update(n.edges[e])
    
    def labelGraph(self, storage, minOrg):	
        num_orgs=len(storage.summaryLookup.keys())
        total_tax=len(storage.completeTaxSummary())
        self.update_edges(weight_attr='getOrganism',divisor=float(num_orgs), label_attr=('getReplicon','replicons'), remove_attrs=['instances'])
        self.update_node_attr_final(weight_func=storage.nodeTaxSummary, family_func=storage.getFamilyInfo, divisor=float(total_tax), remove_attrs=['instances'], minOrg=minOrg)
        
        #create attribute called paths which represents edges per replicon
        #self["paths"]=';'.join([k+':'+','.join(v) for k,v in storage.replicon_edges_dict.iteritems()])
            

        #get list of nodes and edges for testing
        #node_handle=open('new_loop_node_list.txt','w')
        #for n in self.nodes_iter():
        #	node_handle.write(n+"\n")
        #node_handle.close()
        #edge_handle=open('new_loop_edge_list.txt','w')
        #for e in self.edges_iter():
        #	edge_handle.write(str(e)+"\n")
        #edge_handle.close()



    def toXGMML(self, fhandle):
        xml = DOMLight.XMLMaker()
        fhandle.write("""<?xml version="1.0" encoding="UTF-8"?>
        <graph xmlns="http://www.cs.rpi.edu/XGMML" directed="0" label="PFam assembly">
        """)
        cid = 0
        cur_ids = {}
        for cn in self.nodes_iter():
            cur_ids[cn] = cid
            fhandle.write(str(xml.node({'id': cid, 'label': cn}, '<att type="real" name="weight" value="'+str(self.node[cn]['weight'])+'"/>')) + "\n")
            cid += 1
        count = 0
        for edge in self.edges_weight_iter():
            dom_edge=xml.edge()
            dom_edge.set({'weight': edge[2]['weight'], 'source': cur_ids[edge[0]], 'target': cur_ids[edge[1]], 'label': ""}, '<att type="real" name="weight" value="'+str(edge[2]['weight'])+'"/>')
            fhandle.write(str(dom_edge) + "\n")
            #if count == 1000:
            #	break
            count += 1
        fhandle.write("</graph>")
                            
    def toJSON(self, fhandle):
        cid = 0
        cur_ids = {}
        #fhandle.write("{\n\tnodes:[\n")
        results={"nodes" : [], "links" :[]}
        for cn in self.nodes_iter():
            cur_ids[cn] = cid
            results["nodes"].append({'id': cid, 'label': cn, 'weight': self.node[cn]['weight']})
            #fhandle.write(json.dumps({'id': cid, 'label': cn, 'weight': str(self.node[cn]['weight'])})+"\n")
            cid += 1
        #fhandle.write("\t],\n")
        #fhandle.write("\tlinks:[\n")
        count = 0
        for edge in self.edges_weight_iter():
            #fhandle.write(json.dumps({'source': cur_ids[edge[0]], 'target': cur_ids[edge[1]], 'weight': edge[2]['weight']})+"\n")
            results["links"].append({'source': cur_ids[edge[0]], 'target': cur_ids[edge[1]], 'weight': edge[2]['weight']})
            #if count == 1000:
            #	break
            count += 1
        #fhandle.write("\t]\n}")
        fhandle.write(json.dumps(results, indent=1))

    ## Get weighted edgesD from this graph.
    #def edges(self):
        # This is just the code from networkx.graph - except call our
    #	return list(self.edges_iter())
                

    ## Overwrite Graph edges_iter method so we get weight information too
    def edges_weight_iter(self, nbunch=None):
        for edge in Graph.edges_iter(self, nbunch, data=True):
            yield edge

def toGML(cur_graph, file_name):
        readwrite.graphml.write_graphml(cur_graph, file_name)

def node_to_gff(gff_handle, node, feature_counter, graphID):
    for i in node.instances:
        fid=next(feature_counter)
        contig_id, start, end = i.getLocation()
        #for now don't know direction, or feature id
        line="\t".join([contig_id, 'PanGraph', 'match', str(start), str(end), str(len(node.instances)), '+', '.', ";".join(["ID="+str(fid),"Name="+node.famID,"graphID="+str(graphID)])])
        gff_handle.write(line)

def edge_to_gff(gff_handle, edge):
    return 0

#create data maps for keeping track of relationships between graph entities
#finalaize node ids and edge ids
#create GFF records
#org_map genome id to name
#sid_to_edge sequence id to edge id
#takes number of nodes so that edges can have unique ids wrt to nodes
def create_indices(storage, pgraph, csize, gff_outfile):
    bgzf_handle = bgzf.BgzfWriter(gff_outfile,'wb')
    edge_counter=itertools.count()
    node_counter=itertools.count()
    feature_counter=itertools.count()
    storage.org_map={}#maps which genome ids have which names
    storage.sid_to_edge={}#maps which sequence ids have which edges
    storage.graph_to_offset={}#maps graph ID (node or edge) to offset location start,end
    for n in pgraph.nodes():
        ncount=str(next(node_counter))
        pgraph.node[n]['id']=ncount
        start_voff=bgzf_handle.tell()
        node_to_gff(gff_handle=bgzf_handle,node=n, feature_counter=feature_counter, graphID=ncount)
        end_voff=bgzf_handle.tell()
        storage.graph_to_offset[ncount]=[start_voff,end_voff]
    for e in pgraph.edges():
        ecount=next(edge_counter)
        pgraph.adj[e[0]][e[1]]['id']=str(ecount+csize)
        for r in (pgraph.adj[e[0]][e[1]]['replicons']).split(','):
            storage.sid_to_edge.setdefault(r,[]).append(pgraph.adj[e[0]][e[1]]['id'])
            #edge_to_gff(bgzf_handle)
    for k,v in storage.summaryLookup.iteritems():
        storage.org_map[k]=v.genome_name



#remove certain attributes that are no longer useful or prohibiitively large from graph output
def remove_attributes(pgraph, from_edges=[], from_nodes=[]):
    if from_edges:
        for e in pgraph.edges():
            for r in from_edges:
                try: pgraph.adj[e[0]][e[1]].pop(r,None)
                except: pass
    #from nodes need to implement
    if from_nodes:
        for n in pgraph.nodes():
            for r in from_nodes:
                try: delattr(n,r)
                except: pass

#possibly assembly mistakes
#or genome rearrangement points
def find_rearrangements(pgraph, storage, out_file, gminimum=None):
    out_handle=open(out_file,'w')
    if not gminimum:
        gminimum=len(storage.summaryLookup)#default to all genomes in 
    for u,v,data in pgraph.edges_iter(data=True):
        us=u.get_summary()
        vs=v.get_summary()
        if len(data['instances']) == 1 and len(us['organisms']) >= gminimum and len(vs['organisms']) >= gminimum:
            out_handle.write(next(iter(data['instances'])).getLocationString()+"\n")
    out_handle.close()
        
        
def modGexf(in_handle, out_file, k_size, minOrg, storage, pgraph):
    #register_namespace('',"http://www.gexf.net/1.1draft")
    encoding='utf-8'
    header='<?xml version="1.0" encoding="%s"?>'%encoding
    gexf_xml=ElementTree(fromstring(in_handle.getvalue()))
    metadata_element=Element("meta")
    metadata_element.append(Element("ksize",value=str(k_size)))
    metadata_element.append(Element("minorg",value=str(minOrg)))
    gn_element=Element("org_map")
    gn_element.text=CDATA(json.dumps(storage.org_map))
    metadata_element.append(gn_element)
    contig_element= Element("contig_map")
    contig_element.text = CDATA(json.dumps(storage.replicon_map))
    metadata_element.append(contig_element)
    sid_element = Element("edge_map")
    sid_element.text = CDATA(json.dumps(storage.sid_to_edge))
    metadata_element.append(sid_element)
    root=gexf_xml.getroot()
    root.insert(0, metadata_element)
    gexf_handle=open(out_file, 'w')
    gexf_handle.write(header.encode(encoding))
    gexf_xml.write(gexf_handle, encoding=encoding)
    gexf_handle.close()

#calculate graph statistics
def stats(graph):
    num_nodes=graph.order()
    num_edges=graph.size()
    avg_degree= float(num_edges)/num_nodes
    print "\t".join([str(num_nodes),str(num_edges),str(avg_degree)])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("feature_table", type=str, help="table specifying the group, genome, contig, feature, and start in sorted order")
    parser.add_argument("output_basename", type=str, help="the path and base name give to the output files")
    parser.add_argument("context", type=str, choices=["genome","contig","feature"], help="the synteny context")
    parser.add_argument("ksize", type=int, choices=range(3,10), help="the size of the kmer to use in constructing synteny")
    parser.add_argument('--break_conflict', help='Uses methods for dealing with latent updating to APIs', required=False, default=False, action='store_true')

    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    gmaker=GraphMaker(feature_tab=args.feature_table, context=args.context, ksize=args.ksize, break_conflict=args.break_conflict)
    gmaker.processFeatures()
    gmaker.RF_to_PG()
    nx.readwrite.write_gexf(gmaker.rf_graph, args.output_basename+".rf_graph.gexf")
    gmaker.checkPGGraph()
    gmaker.checkRFGraph()
    gmaker.calcStatistics()
    gmaker.finalizeGraphAttr()
    nx.readwrite.write_gexf(gmaker.pg_graph, args.output_basename+".gexf")


def old_main(init_args):
    if(len(init_args)<5):
        sys.stderr.write("Usage: figfam_to_graph.py feature_table family_table summary_table output_folder k-size minOrg\n")
        sys.exit()
    k_size=int(init_args[4])
    minOrg=int(init_args[5])
    if len(init_args)>=7:
        ignore_fams=init_args[6].replace(' ','').split(',')
    #ignore_fams=set(['FIG00638284','FIG01306568'])
    fstorage=Storage(init_args[0], init_args[1], init_args[2], k_size, ignore_fams=set(['FIG00638284','FIG01306568']))
    fstorage.bfsExpand(minOrg)
    out_basename=os.path.splitext(os.path.basename(init_args[0]))[0] #get basename of the file to name output
    out_folder=os.path.expanduser(init_args[3])
    out_file=os.path.join(out_folder,out_basename)
    pgraph=pFamGraph(fstorage,minOrg=minOrg)
    find_rearrangements(pgraph, fstorage, out_file+"_rearrangements.txt")
    pgraph.labelGraph(fstorage,minOrg=minOrg) #label/weight nodes and edges. also remove anything under minOrg
    csize=pgraph.order()
    create_indices(fstorage, pgraph, csize, out_file+".gff.gz")
    remove_attributes(pgraph, from_edges=["replicons"], from_nodes=["locations","organisms"])
    toGML(pgraph, out_file+".graphml")
    gexf_capture=StringIO()#lazy instead of patching NetworkX to include meta attribute. capture, mod xml.
    readwrite.write_gexf(pgraph, gexf_capture)
    modGexf(gexf_capture, out_file+".gexf", k_size, minOrg, fstorage, pgraph)
    result_handle=open(out_file+".xgmml", 'w')
    pgraph.toXGMML(result_handle)
    result_handle.close()
    stats(pgraph)
    #result_handle=open(out_file+".json", 'w')
    #pgraph.toJSON(result_handle)
    #result_handle.close()
    
if __name__ == "__main__":
    main()
