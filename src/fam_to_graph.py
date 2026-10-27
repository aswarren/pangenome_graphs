#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import itertools
import json
import re
import copy
import argparse
import networkx as nx
import fileinput
from operator import methodcaller
from operator import itemgetter
#from networkx import Graph
#from networkx import readwrite
#from Bio import bgzf
import json
from collections import deque
from collections import OrderedDict
from io import StringIO
import io
from subprocess import Popen, PIPE, STDOUT
import time
import requests
import logging
import urllib 
import collections

embed_info = False

##This script expects files which contain a list of consecutive protein family ID's for the replicons in an organism
## and some kind of summary information about where the kmers come from

##Class for storing information about the origin of a Kmer
class featureInfo():
        #expand to parse out this information from different sources
    def __init__(self, line=None):
        self.md5=None
        self.contig_id=None
        self.genome_id=None
        self.feature_id=None
        self.feature_ref=None
        self.group_id=None
        self.group_num=None
        self.start=None
        self.end=None
        self.strand='+'
        self.function=""
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

    def getContextValue(self, context):
        if context=="genome":
            return self.genome_id
        elif context=="contig":
            return self.contig_id
        else:
            return None


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
        self.numBins =0
    
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


        
class featureParser():
    def __init__(self, **kwargs):
        self.feature_files=kwargs['feature_files']
        self.feature_selector=kwargs['feature_selector']
        self.file_type=kwargs['file_type']
        self.parse_function=kwargs['parse_function']
        self.alpha=kwargs["alpha"]
        self.parse=None
        self.ip = None
        self.plaintab={'genome':0,'contig':1,'feature':2,'start':3, 'end':4, 'group':5}
        #self.ip={'taxid':2, 'genome':1, 'contig':3,'feature':2,'start':4, 'end':5, 'group':0}
        self.pc_figfam={'genome':0,'contig':2,'feature':5,'start':9, 'end':10, 'strand':11, 'group':15, 'function':14, 'organism':1}
        self.pc_plfam={'genome':0,'contig':2,'feature':5,'start':9, 'end':10, 'strand':11, 'group':16, 'function':14, 'organism':1}
        self.pc_pgfam={'genome':0,'contig':2,'feature':5,'start':9, 'end':10, 'strand':11, 'group':17, 'function':14, 'organism':1}
        if self.file_type=="patric_tab" or self.file_type=="patric_genomes":
            self.parse=self.parseFeatureTab
            self.ip = self.plaintab
            if self.alpha=="figfam_id":
                self.ip = self.pc_figfam
            if self.alpha=="plfam_id":
                self.ip = self.pc_plfam
            if self.alpha=="pgfam_id":
                self.ip = self.pc_pgfam
    

    def parseFeatureTab(self):
        #if parsing the feature_files to download patric genome ids, need to use the generator defined by genome_id_feature_gen
        if self.file_type=="patric_genomes":
            generator = self.genome_id_feature_gen()
        else:
            generator = fileinput.input(files=self.feature_files)
        #if files is empty it should read from stdin
#         for line in fileinput.input(files=self.feature_files):
        for line in generator:
            result=featureInfo()
            header=False
            try:
                header = line.startswith('#')
                if header:
                    #define column position based on header
                    parts=line.strip().replace("#","").split("\t")
                    valid_cols = {'genome', 'contig', 'feature', 'start', 'end', 'strand', 'group', 'function', 'organism'}
                    new_ip = {}
                    x=0
                    while x < len(parts):
                        cur_part=parts[x].lower()
                        if cur_part in valid_cols or cur_part in self.ip:
                            new_ip[cur_part]=x
                        x+=1
                    self.ip = new_ip
                    continue

                else:
                    parts =[p.replace('"', '') for p in line.strip().split("\t")]
                    result.group_id=parts[self.ip['group']]
                    result.contig_id=parts[self.ip['contig']]
                    result.genome_id=parts[self.ip['genome']]
                    result.start=int(parts[self.ip['start']])
                    result.end=int(parts[self.ip['end']])
                    if 'strand' in self.ip and self.ip['strand'] < len(parts):
                        result.strand = parts[self.ip['strand']]
                    result.feature_ref = parts[self.ip['feature']]
                    if 'organism' in self.ip and self.ip['organism'] < len(parts):
                        result.organism = parts[self.ip['organism']]
                    
                    if self.file_type=="patric_genomes" or self.alpha=="pgfam_id":
                        result.md5 = parts[-1]
                    
                    if self.parse_function and 'function' in self.ip and self.ip['function'] < len(parts):
                        result.function = parts[self.ip['function']]
            except:
                logging.debug("warning: couldn't parse line: "+line)
                continue
            yield result

    def chunker(self, seq, size):
        return (seq[pos:pos + size] for pos in range(0, len(seq), size))

    def genome_id_feature_gen(self, limit=100000):
        genome_ids =[]

        # Will open all files, or stdin if no arguments passed
        for line in fileinput.input(files=self.feature_files):
            delim = r'; |, |,|;| |\t'
            parsed = re.split(delim, line.strip())
            
            for l in parsed:
                # FIX 1: Strip single quotes AND double quotes
                clean_id = l.strip().replace('"', '').replace("'", "")
                if clean_id and clean_id != "genome_id":
                    genome_ids.append(clean_id)
                    
        logging.info("genome_ids: " + ",".join(genome_ids))
        
        # Process in chunks of 10 genomes
        for gids in self.chunker(genome_ids, 10):
            gids_str = ",".join(gids)
            selectors =[
                "ne(feature_type,source)",
                "eq(annotation,PATRIC)",
                f"in({self.feature_selector},({gids_str}))"
            ]
            genomes = f"and({','.join(selectors)})"   
            
            limit_str = f"limit({limit})"
            
            select = "select(genome_id,genome_name,accession,annotation,feature_type,patric_id,refseq_locus_tag,alt_locus_tag,uniprotkb_accession,start,end,strand,na_length,gene,product,figfam_id,plfam_id,pgfam_id,go,ec,pathway,aa_sequence_md5)&sort(+genome_id,+accession,+start)"
            base = "https://www.bv-brc.org/api/genome_feature/" 
            query = "&".join([genomes, limit_str, select])
            
            # --- FIX 2: URL Encode the query to protect the + signs! ---
            # We keep RQL syntax safe, but + will become %2B
            encoded_query = urllib.parse.quote(query, safe='&()=,')

            # ---------------------------------------------------
            # UNCOMMENT THIS BLOCK TO GET A BROWSER-TESTABLE GET URL
            #get_url = f"{base}?{urllib.parse.quote(query, safe='&()=+,')}"
            #logging.warning(f"\n[DEBUG] Browser GET URL:\n{get_url}\n")
            # ---------------------------------------------------
            
            headers = {"accept": "text/tsv", "content-type": "application/rqlquery+x-www-form-urlencoded"}

            # Send the ENCODED query in the POST body
            r = requests.post(url=base, data=encoded_query, headers=headers, stream=True) 
            if r.encoding is None:
                r.encoding = "utf-8"
                
            if not r.ok:
                logging.warning(f"Error in API request: {r.status_code} {r.reason}\n")
                
            for line in r.iter_lines(decode_unicode=True):
                yield line

##CALCULATE DIVERSITY QUOTIENT!!! GENUS/TOTAL GENOMES
##CALCULATE NORMALIZED NUMBER WEIGHT of NUMBER OF genomes in edge/ total number of genomes

##This class is for storing dictionary structures that facilitate different pan-genome graph calculations
##Parameters are filepaths and the size of kmer to use
#GraphMaker(feature_tab=some_file, context="genome")
class GraphMaker():
    def __init__(self, **kwargs):#feature_file, family_file, summary_file, ksize, ignore_fams=set([])):
        self.feature_parser=None
        label_func = kwargs.get("label_function", False)
        alpha = kwargs.get("alpha", "pgfam_id")
        #convert option passed to file_type
        self.feature_parser=featureParser(feature_files=kwargs["feature_files"], feature_selector=kwargs.get("feature_selector", "genome_id"), file_type=kwargs.get("file_type", "patric_tab"), parse_function=label_func, alpha=alpha)
        self.context=kwargs.get("context", "genome") #should be ["genome", "contig", "feature"]
        self.context_levels={"genome":0,"contig":1,"feature":2}
        self.ksize=kwargs.get("ksize", 3)
        self.break_conflict=kwargs.get("break_conflict", False)
        self.label_function=label_func
        self.diversity=kwargs.get("diversity", "genus")
        self.all_diversity={}#for tracking total taxa numbers in the graph
        self.num_pg_nodes=0
        self.rf_graph=nx.DiGraph()# the rf-graph (close to de bruijn) created from series of features with group designations
        self.pg_graph=pFamGraph()# pg-graph is an undirected grpah
        self.rf_node_index=[]
        self.replicon_map=OrderedDict()#stores {genome_id:OrderedDict(contig_id)}
        self.minSeq = kwargs.get("minSeq", 1)
        self.no_edge=set([])
        self.conflicts={}
        self.contig_order={}#stores {genome_id:OrderedDict{contig_id:[feature_number]}}
        self.contig_unorder={}#stores {genome_id:OrderedDict{contig_id:[feature_number]}}
        self.order_contigs=False
        self.contig_to_rf ={}# track which contigs have which rf-nodes
        self.contig_weight={} # for calculating which contigs should be prioritized in traversal to use in layout algorithm
        self.traverse_priority = None
        self.groups_seen={}
        self.group_index=[]
        self.context_bin=set([])
        self.feature_index=[]
        self.debug = False
        self.prebundle = False
        self.eat_repeats =False
        self.anchor_instance_keys={} #structured as {instance_key: [pg_assignment, [features with this instance key]]}
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
        self.alt_counter =0
        self.similarity_matrix = collections.defaultdict(lambda: collections.defaultdict(int))
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
        ambig=0
        for r in self.rf_node_index:
            if r.numFeatures() > 0:
                ambig+=1
        logging.info("rf-graph: "+str(ambig)+" nodes unexapanded")
                #assert LogicError("RFNode unexpanded")
    def checkPGGraph(self):
        for cnode in self.pg_graph.nodes(data=True):
            if len(cnode[1]["features"]) == 0 :
                logging.warning("pg-graph node "+str(cnode[0])+"has no features")
            group_id=None
            for g in cnode[1]["features"]:
                
                if g in ["md5", "start", "end", "info"]:
                    continue
                
                for contig in cnode[1]["features"][g]:
                    for f in cnode[1]["features"][g][contig]:
                        if group_id == None:
                            group_id = self.feature_index[int(f)].group_id
                        elif group_id != self.feature_index[int(f)].group_id:
                            raise Exception("LogicError")

    def calcStatistics(self):
        stats=[]
        stats.append("rf-graph:")
        stats.append("nodes "+str(self.rf_graph.number_of_nodes()))
        stats.append("edges "+str(self.rf_graph.number_of_edges()))
        stats.append("ps-graph:")
        stats.append("nodes "+str(self.pg_graph.number_of_nodes()))
        stats.append("edges "+str(self.pg_graph.number_of_edges()))
        stats.append("alt-nodes "+str(self.alt_counter))
        logging.warning("\n"+"\n".join(stats))

    def write_contigs(self, contig_file, unsorted_file=None):
            
        #first assure that all contigs are represented
        missing_contigs=0
        missing_genomes=0
        for k,v in self.replicon_map.items():
            if k not in self.contig_order:
                missing_genomes+=1
                logging.warning("WARNING: missing genome in contig order "+k+"\n")
                self.contig_order.setdefault(k, OrderedDict())
            for contig_id in v.keys():
                if contig_id not in self.contig_order[k]:
                    missing_contigs+=1
                    logging.warning("WARNING: missing contig in contig order "+" ".join([k,contig_id])+"\n")
                    if unsorted_file != None:
                        self.contig_unorder.setdefault(k,{}).setdefault(contig_id,[])
                    else:
                        self.contig_order[k].setdefault(contig_id,[])
        if missing_contigs or missing_genomes:
            logging.warning("WARNING: missing genomes count "+str(missing_genomes)+" missing contigs count "+str(missing_contigs)+"\n")
        with open(contig_file, 'w') as ch:
            for g in self.replicon_map.keys():
                ch.write("\t".join([g] + list(self.contig_order[g].keys())) + "\n")
        if unsorted_file != None:
            with open(unsorted_file, 'w') as ch:
                for g in self.contig_unorder.keys():
                    ch.write("\t".join([g] + list(self.contig_unorder[g].keys())) + "\n")

    def write_reports(self, base_path):
        import csv
        
        # Strip extension to use as a prefix
        prefix = base_path.replace(".gexf", "").replace(".gfa", "")
        
        logging.warning("Exporting 5 tabular TSV reports...")

        # 1. Synteny Similarity Matrix
        with open(f"{prefix}_synteny_similarity.tsv", 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            keys = sorted(self.similarity_matrix.keys())
            writer.writerow(["ID"] + keys)
            for k1 in keys:
                row = [k1] + [self.similarity_matrix[k1].get(k2, 0) for k2 in keys]
                writer.writerow(row)
                
        # 2. Scaffolding Recommendations
        with open(f"{prefix}_scaffold_recommendations.tsv", 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["Node_U_Family", "Node_V_Family", "Target_Broken_Genomes", "Consensus_Support_Weight"])
            for u, v, d in self.pg_graph.edges(data=True):
                if d.get("is_scaffold_path"):
                    u_fam = self.pg_graph.nodes[u].get('family', str(u))
                    v_fam = self.pg_graph.nodes[v].get('family', str(v))
                    targets = d.get("scaffold_entities") or d.get("sv_entities", "")
                    weight = d.get("weight", 0.0)
                    if targets: # Only write the edges that carry the target names
                        writer.writerow([u_fam, v_fam, targets, f"{weight:.4f}"])

        # 3. Mobilized Genes / Dispersed Paralogs
        with open(f"{prefix}_mobilized_gene_families.tsv", 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["CNV_Cluster_ID", "Protein_Family", "Function", "Total_Copies", "Replicons_Involved"])
            for n, d in self.pg_graph.nodes(data=True):
                if d.get("repeat_ambiguity_fragmentation"):
                    cluster_id = d.get("cnv_cluster_id", 0)
                    fam = d.get("family", "")
                    func = d.get("label", "")
                    
                    replicons = set()
                    copies = 0
                    for g_id, c_dict in d.get("features", {}).items():
                        if g_id in ["info", "md5", "start", "end"]: continue
                        for c_id, f_list in c_dict.items():
                            replicons.add(c_id)
                            copies += len(f_list)
                    writer.writerow([cluster_id, fam, func, copies, ",".join(replicons)])

        # 4. Structural Variants & Junctions
        sv_junction_flags = [
            "breakpoint_junction",
            "inter_contig_rearrangement",
            "intra_contig_rearrangement",
            "alt_path_junction",
            "repeat_fragmentation_junction",
            "repeat_ambiguity_detour"
        ]
        with open(f"{prefix}_structural_variants.tsv", 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["Junction_Type", "Node_U_Family", "Node_V_Family", "Divergent_Genomes", "Consensus_Support_Weight"])
            for u, v, d in self.pg_graph.edges(data=True):
                active_jts = [flag for flag in sv_junction_flags if d.get(flag)]
                if not active_jts and d.get("junction_type"):
                    active_jts = [d.get("junction_type")]
                if active_jts:
                    jt = ",".join(active_jts)
                    u_fam = self.pg_graph.nodes[u].get('family', str(u))
                    v_fam = self.pg_graph.nodes[v].get('family', str(v))
                    sv_ents = d.get("sv_entities", "")
                    weight = d.get("weight", 0.0)
                    writer.writerow([jt, u_fam, v_fam, sv_ents, f"{weight:.4f}"])

        # 5. Superbubbles Catalog
        if hasattr(self, 'completed_superbubbles'):
            with open(f"{prefix}_superbubbles.tsv", 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(["Superbubble_ID", "Entry_Node_Family", "Exit_Node_Family", "Total_Allelic_Paths", "Max_Path_Length", "Internal_Genes"])
                for bid, bubble in enumerate(self.completed_superbubbles, 1):
                    entry_n = bubble['origin']
                    paths = bubble['winning_paths']
                    if not paths: continue
                    exit_n = paths[0][-1]
                    
                    entry_fam = self.pg_graph.nodes[entry_n].get("family", str(entry_n))
                    exit_fam = self.pg_graph.nodes[exit_n].get("family", str(exit_n))
                    
                    max_len = max(len(p) for p in paths)
                    internal_genes = set()
                    for p in paths:
                        for internal_n in p[:-1]: # Don't list the exit node as internal
                            internal_genes.add(self.pg_graph.nodes[internal_n].get("family", str(internal_n)))
                            
                    writer.writerow([bid, entry_fam, exit_fam, len(paths), max_len, ",".join(internal_genes)])

    def getTaxaIndicator(self, feature_id, mode="name"):
            if mode=="name":
                string_part = 0
                if self.diversity == "genus":
                    string_part =0
                elif self.diversity == "species":
                    string_part = 1
                cur_taxa =self.feature_index[feature_id].organism.split()[string_part]
                return cur_taxa

    def trackDiversity(self, feature_id, tracking_dict):
        taxa = self.getTaxaIndicator(feature_id)
        tracking_dict.setdefault(taxa, 0)
        tracking_dict[taxa]+=1

    def calcDiversity(self, cur_profile):
        #taxa=self.getTaxaIndicator(feature_id)
        diversity = float(len(cur_profile.keys()))/float(len(self.all_diversity.keys()))
        return diversity

    def detect_and_annotate_superbubbles(self, max_depth=50):
        """
        TFS Superbubble Search.
        "Homebase" explores the graph using TFS (weight-prioritized DFS), 
        carrying "Hunts" (minority threads) as passengers.
        Prevents combinatorial explosion by stopping hunts upon reconvergence 
        and only spawning hunts in unvisited territory.
        """
        logging.warning("Detecting superbubbles via TFS Bounded Hunts...")
        UG = self.get_undirected_pg_graph()
        
        def get_node_genomes(n):
            return set(self.pg_graph.nodes[n].get('features', {}).keys()) - {'info', 'md5', 'start', 'end'}
            
        def get_edge_genomes(u, v):
            if UG.has_edge(u, v):
                return set(UG.edges[u, v].get('genomes', set()))
            return set()

        feat_to_idx = [0] * len(self.feature_index)
        for contigs in self.replicon_map.values():
            for f_arr in contigs.values():
                for idx, f_id in enumerate(f_arr):
                    feat_to_idx[f_id] = idx

        def get_projected_highway_nodes(cv_node, highway_genomes, budget=max_depth):
            import collections
            target_nodes = collections.defaultdict(set)
            feat_dict = self.pg_graph.nodes[cv_node].get('features', {})
            for g in highway_genomes:
                if g in feat_dict:
                    for c_id, f_list in feat_dict[g].items():
                        f_arr = self.replicon_map[g][c_id]
                        for f_id in f_list:
                            start_idx = feat_to_idx[f_id]
                            
                            # --- FIX: Project Forward AND Backward to handle Inversions! ---
                            # Project Forward
                            for i in range(start_idx + 1, min(len(f_arr), start_idx + budget + 1)):
                                nxt_pg = self.feature_index[f_arr[i]].pg_assignment
                                if nxt_pg is not None:
                                    target_nodes[nxt_pg].add(g)
                            # Project Backward
                            for i in range(start_idx - 1, max(-1, start_idx - budget - 1), -1):
                                nxt_pg = self.feature_index[f_arr[i]].pg_assignment
                                if nxt_pg is not None:
                                    target_nodes[nxt_pg].add(g)
                            # ---------------------------------------------------------------
            return target_nodes        

        class HuntState:
            def __init__(self, origin, highway_genomes, target_nodes, budget):
                self.origin = origin
                self.highway_genomes = highway_genomes
                self.target_nodes = target_nodes
                self.budget = budget
                self.best_score = 0
                self.current_path = []
                self.best_path = []
                
            def clone(self):
                new_hunt = HuntState(self.origin, self.highway_genomes, self.target_nodes, self.budget)
                new_hunt.best_score = self.best_score
                new_hunt.current_path = list(self.current_path)
                new_hunt.best_path = list(self.best_path)
                return new_hunt

        class VisitPack:
            def __init__(self, node, incoming_hunts, dfs_path):
                self.node = node
                self.incoming_hunts = incoming_hunts
                self.dfs_path = dfs_path
                self.visited = False
                self.children = []
                self.returned_hunts = []

        visited_global = set()
        hunt_memo = {} # (node, hunt_origin) -> max_budget_seen
        completed_bubbles = []
        
        # Sort nodes by weight to prioritize starting the Homebase on the Core Highways
        sorted_nodes = sorted(self.pg_graph.nodes(), key=lambda n: len(get_node_genomes(n)), reverse=True)
                              
        for start_node in sorted_nodes:
            if start_node in visited_global:
                continue
                
            stack = [VisitPack(start_node, [], set([start_node]))]
            
            while stack:
                cv = stack.pop()
                
                if not cv.visited:
                    cv.visited = True
                    
                    # USER FIX: Check if we are treading on previously mapped territory
                    is_new_territory = cv.node not in visited_global
                    if is_new_territory:
                        visited_global.add(cv.node)
                    
                    active_hunts_to_pass = []
                    node_gens = get_node_genomes(cv.node)
                    
                    # 1. Evaluate Incoming Hunts (The Passengers)
                    for hunt in cv.incoming_hunts:
                        hunt.current_path.append(cv.node)
                        hunt.budget -= 1
                        
                        # --- NEW: Check the Projected Target Mass ---
                        intersect_score = 0
                        if cv.node in hunt.target_nodes:
                            intersect_score = len(hunt.target_nodes[cv.node])
                            
                            # STRICTLY GREATER THAN: Freezes the path at the FIRST node 
                            # of maximum convergence, ignoring the rest of the highway tail.
                            if intersect_score > hunt.best_score:
                                hunt.best_score = intersect_score
                                hunt.best_path = list(hunt.current_path)
                                
                        # We only early-resolve if we hit exactly 100%. 
                        # Otherwise, keep exploring the budget to find the bulk merge!
                        is_100_percent = (intersect_score == len(hunt.highway_genomes))
                        
                        if hunt.budget > 0 and not is_100_percent:
                            memo_key = (cv.node, hunt.origin)
                            if hunt.budget > hunt_memo.get(memo_key, -1):
                                hunt_memo[memo_key] = hunt.budget
                                active_hunts_to_pass.append(hunt)
                            else:
                                cv.returned_hunts.append(hunt) # Pruned
                        else:
                            cv.returned_hunts.append(hunt) # Budget exhausted or 100% resolved
                            
                    # 2. Check for Divergence 
                    successors = list(UG.neighbors(cv.node))
                    succ_weights = [(succ, len(get_edge_genomes(cv.node, succ))) for succ in successors]
                    succ_weights.sort(key=lambda x: x[1], reverse=True)
                    
                    valid_succ_weights = [(s, w) for s, w in succ_weights if s not in cv.dfs_path]
                    stack.append(cv) # Re-push for bottom-up
                    
                    # 3. Prepare Children
                    for succ, weight in reversed(valid_succ_weights):
                        hunts_for_child = [h.clone() for h in active_hunts_to_pass]
                        
                        # USER FIX: ONLY spawn new hunts if we are in unvisited territory!
                        if is_new_territory and len(valid_succ_weights) > 1:
                            heaviest_succ_id = valid_succ_weights[0][0]
                            
                            # --- CRITICAL FIX: ONLY SPAWN SCOUTS DOWN THE DIRT ROAD! ---
                            # Do not let the Highway hunt for itself, or it creates 1-hop fake bubbles!
                            if succ != heaviest_succ_id:
                                highway_fingerprint = get_edge_genomes(cv.node, heaviest_succ_id)
                                
                                if len(highway_fingerprint) > 0:
                                    # Subset Checking 
                                    is_subset = False
                                    for active_hunt in active_hunts_to_pass:
                                        if highway_fingerprint.issubset(active_hunt.highway_genomes):
                                            is_subset = True
                                            break
                                    
                                    if not is_subset:
                                        target_nodes = get_projected_highway_nodes(cv.node, highway_fingerprint, max_depth)
                                        new_hunt = HuntState(origin=cv.node, highway_genomes=highway_fingerprint, target_nodes=target_nodes, budget=max_depth)
                                        hunts_for_child.append(new_hunt)
                            # -------------------------------------------------------------
                                        
                        if succ not in visited_global or hunts_for_child:
                            child_dfs_path = set(cv.dfs_path)
                            child_dfs_path.add(succ)
                            child_pack = VisitPack(succ, hunts_for_child, child_dfs_path)
                            cv.children.append(child_pack)
                            stack.append(child_pack)
                #vastly simplified            
                else:
                    # ---------------------------------------------------------
                    # BOTTOM-UP / TAIL-END RECURSION
                    # ---------------------------------------------------------
                    for child in cv.children:
                        cv.returned_hunts.extend(child.returned_hunts)
                        
                    hunts_to_pass_up = []
                    my_hunts = []
                    for h in cv.returned_hunts:
                        if h.origin == cv.node:
                            my_hunts.append(h)
                        else:
                            hunts_to_pass_up.append(h)
                            
                    if my_hunts:
                         # 1. Grab all hunts that successfully found at least 50% of the highway
                        winning_hunts = []
                        for h in my_hunts:
                            highway_size = len(h.highway_genomes)
                            if h.best_score > 0 and h.best_score >= (highway_size * 0.5):
                                winning_hunts.append(h)
                        
                        if winning_hunts:
                            # 2. What was the absolute best convergence any scout found?
                            max_score = max(h.best_score for h in winning_hunts)
                            best_hunts = [h for h in winning_hunts if h.best_score == max_score]
                            
                            # 3. Group them by their agreed Exit Node
                            exit_groups = {}
                            for h in best_hunts:
                                exit_node = h.best_path[-1]
                                exit_groups.setdefault(exit_node, []).append(h.best_path)
                            
                            # 4. Highway Backfill and Bubble Annotation
                            for exit_node, paths in exit_groups.items():
                                highway_path = []
                                curr = cv.node
                                highway_genomes = winning_hunts[0].highway_genomes
                                
                                # Trace the heaviest edges matching the highway genomes
                                for _ in range(max_depth):
                                    if curr == exit_node:
                                        break
                                    succs = [s for s in UG.neighbors(curr) if s != cv.node and s not in highway_path]
                                    if not succs: break
                                    
                                    best_succ = max(succs, key=lambda s: len(get_edge_genomes(curr, s).intersection(highway_genomes)))
                                    highway_path.append(best_succ)
                                    curr = best_succ
                                    
                                if highway_path and highway_path[-1] == exit_node:
                                    paths.append(highway_path)

                                completed_bubbles.append({
                                    'origin': cv.node,
                                    'winning_paths': paths
                                })
                            
                    cv.returned_hunts = hunts_to_pass_up

        # 5. Annotate Graph
        for n in self.pg_graph.nodes:
            self.pg_graph.nodes[n].setdefault('superbubble_id', set())
            self.pg_graph.nodes[n]['is_superbubble'] = False
        for u, v, d in self.pg_graph.edges(data=True):
            d.setdefault('superbubble_id', set())
            d['is_superbubble'] = False
            
        for bid, bubble in enumerate(completed_bubbles, 1):
                origin = bubble['origin']
                for path in bubble['winning_paths']:
                    
                    full_path = [origin] + path
                    
                    # Tag the Entry, Internal, and Exit nodes so the entire bubble is one unit
                    for n in full_path:
                        self.pg_graph.nodes[n]['superbubble_id'].add(bid)
                        self.pg_graph.nodes[n]['is_superbubble'] = True
                    for i in range(len(full_path)-1):
                        u, v = full_path[i], full_path[i+1]
                        # Tag both directions for undirected UI compatibility
                        if self.pg_graph.has_edge(u, v):
                            self.pg_graph.edges[u, v]['superbubble_id'].add(bid)
                            self.pg_graph.edges[u, v]['is_superbubble'] = True
                        if self.pg_graph.has_edge(v, u):
                            self.pg_graph.edges[v, u]['superbubble_id'].add(bid)
                            self.pg_graph.edges[v, u]['is_superbubble'] = True
        self.completed_superbubbles = completed_bubbles # Save for TSV export
        logging.warning(f"Superbubble Detection: Found and annotated {len(completed_bubbles)} superbubbles.")    


    def flag_bridges(self):
        """
        Scans all walks to find stealth bridges and scaffolding gaps.
        Pass 1: Logs unique topological jumps, including direct-edge translocations.
        Pass 2: Uses Relative Dominance Shielding to protect fractured highways.
        """
        def is_terminal(feat_id):
            feat = self.feature_index[feat_id]
            contig_array = self.replicon_map[feat.genome_id][feat.contig_id]
            return (contig_array[0] == feat_id) or (contig_array[-1] == feat_id)

        def target_is_terminal_in_node(target_genome, pg_node):
            feat_dict = self.pg_graph.nodes[pg_node].get('features', {})
            if target_genome in feat_dict:
                for c_id, f_list in feat_dict[target_genome].items():
                    if any(is_terminal(f) for f in f_list):
                        return True
            return False

        discovered_paths = {}
        

        # ---------------------------------------------------------
        # PASS 1: DISCOVERY
        # ---------------------------------------------------------
        for walker_genome, contigs in self.replicon_map.items():
            for contig_id, feature_array in contigs.items():
                
                last_seen = {} 

                for walk_idx, feature_id in enumerate(feature_array):
                    current_node = self.feature_index[feature_id].pg_assignment
                    if current_node is None: 
                        continue
                    
                    current_cnv = self.pg_graph.nodes[current_node].get('cnv_cluster_id', 0)
                    target_genomes = set(self.pg_graph.nodes[current_node].get('features', {}).keys()) - {'info', 'md5', 'start', 'end'}

                    for target_genome in target_genomes:
                        if target_genome == walker_genome: 
                            continue

                        if target_genome in last_seen:
                            prev = last_seen[target_genome]
                            gap_length = walk_idx - prev['walk_idx'] - 1
                            prev_node = prev['pg_node']
                            
                            # Check for Direct Edge Drop-outs
                            is_direct_dropout = False
                            if gap_length == 0:
                                edge_traversed = False
                                for dir_e in [(prev_node, current_node), (current_node, prev_node)]:
                                    if self.pg_graph.has_edge(*dir_e):
                                        if target_genome in self.pg_graph.edges[dir_e].get('genomes', set()):
                                            edge_traversed = True
                                            break
                                if not edge_traversed:
                                    is_direct_dropout = True
                            
                            # Target dropped out!
                            if gap_length > 0 or is_direct_dropout:
                                prev_cnv = self.pg_graph.nodes[prev_node].get('cnv_cluster_id', 0)

                                if current_cnv == prev_cnv and current_cnv != 0:
                                    continue

                                prev_is_term = target_is_terminal_in_node(target_genome, prev_node)
                                curr_is_term = target_is_terminal_in_node(target_genome, current_node)

                                bridge_features = feature_array[prev['walk_idx']+1 : walk_idx]
                                bridge_nodes = [self.feature_index[f].pg_assignment for f in bridge_features if self.feature_index[f].pg_assignment is not None]
                                
                                path_sequence = [prev_node] + bridge_nodes + [current_node]
                                path_edges = tuple((path_sequence[i], path_sequence[i+1]) for i in range(len(path_sequence)-1))

                                if path_edges not in discovered_paths:
                                    discovered_paths[path_edges] = {
                                        'u': prev_node,
                                        'w': current_node,
                                        'bridge_nodes': set(bridge_nodes),
                                        'path_edges': path_edges,
                                        'is_scaffold': True,
                                        'walkers': set()
                                    }
                                
                                discovered_paths[path_edges]['walkers'].add(walker_genome)
                                
                                if not (prev_is_term and curr_is_term):
                                    discovered_paths[path_edges]['is_scaffold'] = False

                        last_seen[target_genome] = {'pg_node': current_node, 'walk_idx': walk_idx}

        # ---------------------------------------------------------
        # PASS 2: RESOLUTION (DOMINANCE SHIELDING)
        # ---------------------------------------------------------
        alt_path_count = 0
        scaffold_count = 0
        alt_path_event_id = 1

        def normalize_path(path_tuple):
            # Normalize directional paths so Forward and Reverse count as the same Highway
            if path_tuple and path_tuple[0][0] > path_tuple[-1][-1]:
                return tuple((e[1], e[0]) for e in reversed(path_tuple))
            return path_tuple

        # 1. Pre-calculate the maximum combined (Fwd+Rev) support for every (U, W) bubble
        max_supports = {}
        path_supports = {} 
        
        for data in discovered_paths.values():
            u, w = data['u'], data['w']
            norm_uw = tuple(sorted([u, w]))
            norm_path = normalize_path(data['path_edges'])
            
            path_supports[norm_path] = path_supports.get(norm_path, 0) + len(data['walkers'])
            
            if path_supports[norm_path] > max_supports.get(norm_uw, 0):
                max_supports[norm_uw] = path_supports[norm_path]

        # Helper for direction-agnostic edge support
        def get_bidirectional_support(n1, n2):
            gens = set()
            if self.pg_graph.has_edge(n1, n2):
                gens.update(self.pg_graph.edges[n1, n2].get('genomes', set()))
            if self.pg_graph.has_edge(n2, n1):
                gens.update(self.pg_graph.edges[n2, n1].get('genomes', set()))
            return len(gens)

        # 2. Resolve the Bridges
        for path_edges, data in discovered_paths.items():
            u = data['u']
            w = data['w']
            bridge_nodes = data['bridge_nodes']
            is_scaffold = data['is_scaffold']
            
            norm_uw = tuple(sorted([u, w]))
            dominant_support = max_supports.get(norm_uw, 1)
            
            # The threshold is strictly less than the dominant path, unless the dominant path is just 1
            shield_threshold = dominant_support - 1 if dominant_support > 1 else 1

            # If there are no internal nodes between U and W, it is a direct jump (e.g., local deletion).
            # It is handled by tag_structural_variants. Do not flag it here as an MGE bridge!
            if not bridge_nodes:
                continue

            if is_scaffold:
                scaffold_count += 1
                # 1. Tag the nodes
                for n in bridge_nodes:
                    self.pg_graph.nodes[n]['is_scaffold_path'] = True
                
                # 2. Tag the edges of the path!
                # path_edges is a tuple of edges: ((u, x), (x, y), (y, w))
                for edge in path_edges:
                    for direction in [(edge[0], edge[1]), (edge[1], edge[0])]:
                        if self.pg_graph.has_edge(*direction):
                            self.pg_graph.edges[direction]['is_scaffold_path'] = True
                            self.pg_graph.edges[direction]['scaffold_entities'] = ",".join(data['walkers'])
            else:
                #superbubble check: if u and w share a superbubble, this is a local detour, not a translocation
                u_bubbles = self.pg_graph.nodes[u].get('superbubble_id', set())
                w_bubbles = self.pg_graph.nodes[w].get('superbubble_id', set())
                
                # If U and W share a superbubble ID, this path forms a closed
                # topological loop (a local detour/hotspot). It is NOT a translocation.
                if bool(u_bubbles & w_bubbles):
                    continue

                # Node-Level Check: Only flag nodes unique to the bridge
                true_bridge_nodes = []
                for n in bridge_nodes:
                    n_gens = set(self.pg_graph.nodes[n].get('features', {}).keys()) - {'info', 'md5', 'start', 'end'}
                    if len(n_gens) <= shield_threshold:
                        true_bridge_nodes.append(n)
                        self.pg_graph.nodes[n]['alternative_path'] = True
                        self.pg_graph.nodes[n]['alt_path_event_id'] = alt_path_event_id  # <-- Added Event ID
                # Edge-Level Check: Only flag edges unique to the dirt road
                true_bridge_edges = []
                for e in path_edges:
                    support = get_bidirectional_support(e[0], e[1])
                    if support <= shield_threshold:
                        true_bridge_edges.append(e)
                        
                        # Apply tags safely to whatever directional edges exist
                        for direction in [(e[0], e[1]), (e[1], e[0])]:
                            if self.pg_graph.has_edge(*direction):
                                # Internal edges get the boolean overlay and the event ID
                                self.pg_graph.edges[direction]['alternative_path'] = True
                                self.pg_graph.edges[direction]['alt_path_event_id'] = alt_path_event_id

                # Handle Entry and Exit classes (The actual junctions connecting to the backbone!)
                if true_bridge_edges:
                    alt_path_count += 1  # (Assuming you updated translocation_count to alt_path_count)
                    
                    # Tag the Entry Edge
                    edge_in = true_bridge_edges[0]
                    for direction in [(edge_in[0], edge_in[1]), (edge_in[1], edge_in[0])]:
                        if self.pg_graph.has_edge(*direction):
                            self.pg_graph.edges[direction]['alt_path_junction'] = True
                    
                    # Tag the Exit Edge
                    if len(true_bridge_edges) > 1:
                        edge_out = true_bridge_edges[-1]
                        for direction in [(edge_out[0], edge_out[1]), (edge_out[1], edge_out[0])]:
                            if self.pg_graph.has_edge(*direction):
                                self.pg_graph.edges[direction]['alt_path_junction'] = True
                
                # Increment the event ID only if we actually found unprotected elements
                if true_bridge_nodes or true_bridge_edges:
                    alt_path_event_id += 1

        logging.warning(f"Alternative Path Detection: Tagged {alt_path_count} alternative paths and {scaffold_count} scaffolding gaps.")
        return alt_path_count, scaffold_count   

    def compute_graph_metrics(self):
        """
        Calculates node metrics (length, diversity, labels, CNV clusters, and Macro-Conflicts).
        Leaves all features and edges as fast Python dicts and sets.
        """
        num_genomes = float(len(self.replicon_map.keys()))
        alt_group = {}
        processed_n = set([])
        grp_id = 1
        
        # 1. Build alt_groups mapping by protein family (plus any explicit pg_node_alt sets)
        family_to_nodes = collections.defaultdict(set)
        for n, d in self.pg_graph.nodes(data=True):
            first_group = next(i for key, i in d["features"].items() if type(i) == dict and key != "info")
            f_id = next(iter(first_group.values()))[0]
            fam_id = self.feature_index[f_id].group_id
            if fam_id:
                family_to_nodes[fam_id].add(n)

        for fam_id, nodes_set in family_to_nodes.items():
            if len(nodes_set) > 1:
                processed_n.update(nodes_set)
                for alt in nodes_set:
                    alt_group[alt] = grp_id
                grp_id += 1

        for n, alts in self.pg_node_alt.items():
            if n in processed_n:
                continue
            else:
                processed_n.update(alts)
                for alt in alts:
                    alt_group[alt] = grp_id
                grp_id += 1
                
        # 2. Macro-Node (CNV Cluster) Conflict Detection
        clusters = {}
        for node, c_id in alt_group.items():
            clusters.setdefault(c_id,[]).append(node)
            
        macro_conflicts = set()
        for c_id, nodes in clusters.items():
            cluster_genomes = {}
            for node in nodes:
                feat_dict = self.pg_graph.nodes[node].get("features", {})
                for gen_id, seq_dict in feat_dict.items():
                    if gen_id in ["info", "md5", "start", "end"]: continue
                    cluster_genomes.setdefault(gen_id, set()).update(seq_dict.keys())
            
            if any(len(contigs) > 1 for contigs in cluster_genomes.values()):
                macro_conflicts.update(nodes)

        # 3. Compute Node Metrics
        for n, d in self.pg_graph.nodes(data=True):
            label_set = set([])
            cur_diversity = {}

            first_group = next(i for key, i in d["features"].items() if type(i) == dict and key != "info")
            f_id = next(iter(first_group.values()))[0]
            d["label"] = self.feature_index[f_id].group_id
            d["family"] = d["label"]
            
            total_len = 0
            feat_count = 0
            for g_id, c_dict in d["features"].items():
                if g_id in ["info", "md5", "start", "end"]: continue
                for c_id, f_list in c_dict.items():
                    for f_idx in f_list:
                        f_obj = self.feature_index[f_idx]
                        total_len += abs(int(f_obj.end) - int(f_obj.start)) + 1
                        feat_count += 1
            
            d["length"] = int(total_len / feat_count) if feat_count > 0 else 1000

            for g in d["features"]:
                if g in ["md5", "start", "end", "info"]: continue
                
                f_id = next(iter(d["features"][g].values()))[0]
                self.trackDiversity(f_id, cur_diversity)
                for s in d["features"][g]:
                    feature_refs = []
                    for f in d["features"][g][s]:
                        feature_refs.append(self.feature_index[f].feature_ref)
                        if self.label_function:
                            label_set.add(self.feature_index[f].function)
                    
            d["diversity"] = self.calcDiversity(cur_diversity)            
            if self.label_function and label_set:
                d["label"] = list(label_set)[0]
                
            if n in macro_conflicts:
                d["repeat_ambiguity_fragmentation"] = True
                
            d["cnv_cluster_id"] = alt_group.get(n, 0)
            
        # 4. Compute Edge Weights (Leave them as sets!)
        total_contigs = float(sum(len(c) for c in self.replicon_map.values()))
        
        for u, v, attr in self.pg_graph.edges(data=True):
            if self.context == "contig" and "sequences" in attr:
                attr["weight"] = len(attr["sequences"]) / total_contigs if total_contigs > 0 else 0.0
            elif "genomes" in attr:
                attr["weight"] = len(attr["genomes"]) / num_genomes if num_genomes > 0 else 0.0

    def serialize_graph_for_gexf(self, target_graph=None):
            """
            Translates internal feature integers to external string references,
            and converts sets/dicts into JSON/CSV strings for GEXF XML export.
            MUST run after all graph analytics and GFA exports are complete.
            """
            if target_graph is None:
                target_graph = self.pg_graph
                
            for n, d in target_graph.nodes(data=True):
                if "features" in d and not isinstance(d["features"], str):
                    formatted_features = {}
                    for g_id, c_dict in d["features"].items():
                        if g_id in ["info", "md5", "start", "end"]:
                            formatted_features[g_id] = c_dict
                            continue
                        
                        formatted_features[g_id] = {}
                        for c_id, f_list in c_dict.items():
                            # Translate internal integer ID to external string reference
                            formatted_features[g_id][c_id] = [self.feature_index[f].feature_ref for f in f_list]
                            
                    d["features"] = json.dumps(formatted_features)
                
                for key, val in list(d.items()):
                    if isinstance(val, set):
                        d[key] = ','.join(map(str, val))
                    
            for u, v, attr in target_graph.edges(data=True):
                for a in list(attr.keys()):
                    if isinstance(attr[a], set):
                        attr[a] = ','.join(map(str, attr[a]))
            
    
    
    def tag_structural_variants(self):
        """
        Analyzes the generated PS-graph to explicitly tag SVs based on the chosen Context.
        Uses "Path Drop-Out" and "Dangling End" rules to precisely identify SVs while
        protecting against fragmented assemblies.
        """
        sv_edges = 0
        inv_graph = nx.Graph()
        UG = self.get_undirected_pg_graph()  # Needed for tracing precise CNV detours

        
        # --- NEW: Bring in the Terminal Check Helpers ---
        def is_terminal(feat_id):
            feat = self.feature_index[feat_id]
            contig_array = self.replicon_map[feat.genome_id][feat.contig_id]
            return (contig_array[0] == feat_id) or (contig_array[-1] == feat_id)

        def target_is_terminal_in_node(target_genome, pg_node):
            feat_dict = self.pg_graph.nodes[pg_node].get('features', {})
            if target_genome in feat_dict:
                for c_id, f_list in feat_dict[target_genome].items():
                    if any(is_terminal(f) for f in f_list):
                        return True
            return False
        # ------------------------------------------------
        
        def get_edge_genomes(n1, n2):
            if UG.has_edge(n1, n2):
                return set(UG.edges[n1, n2].get('genomes', set()))
            return set()
            
        for u, v, d in self.pg_graph.edges(data=True):
            d["inverted_block"] = False
            
            # Condition 1: Topological Inversion (Bidirectional edges)
            is_topological_inv = self.pg_graph.has_edge(v, u)
            
            # Condition 2: Edge-level Strand Flip (Multiple strands on the same edge)
            is_edge_strand_inv = len(d.get("strands", set())) > 1
            
            # Condition 3: Node-level Strand Flip (Both flanking genes are inverted!)
            is_node_inv = self.pg_graph.nodes[u].get("inverted_block") and self.pg_graph.nodes[v].get("inverted_block")
            
            if is_topological_inv or is_edge_strand_inv or is_node_inv:
                d["inverted_block"] = True
                inv_graph.add_edge(u, v)

            u_node = self.pg_graph.nodes[u]
            v_node = self.pg_graph.nodes[v]
            
            # --- REMOVED THE TARGET_CLASSES GATE ---
            # Now we evaluate all edges for drop-outs to catch clean SVs!
            drop_outs = set()
            
            u_feat = self.pg_graph.nodes[u].get('features', {})
            v_feat = self.pg_graph.nodes[v].get('features', {})
            
            u_cnv = self.pg_graph.nodes[u].get('cnv_cluster_id', 0)
            v_cnv = self.pg_graph.nodes[v].get('cnv_cluster_id', 0)
            
            if self.context == "genome":
                u_items = set(u_feat.keys()) - {'info', 'md5', 'start', 'end'}
                v_items = set(v_feat.keys()) - {'info', 'md5', 'start', 'end'}
                
                # Pool bidirectional support to prevent inversions from flagging as dropouts
                edge_items = set(d.get('genomes', set()))
                if self.pg_graph.has_edge(v, u):
                    edge_items.update(self.pg_graph.edges[v, u].get('genomes', set()))
                    
                shared_items = u_items.intersection(v_items)
                drop_outs = shared_items - edge_items
                
            elif self.context == "contig":
                u_items, v_items = set(), set()
                for gen_dict in u_feat.values():
                    if isinstance(gen_dict, dict): u_items.update(gen_dict.keys())
                for gen_dict in v_feat.values():
                    if isinstance(gen_dict, dict): v_items.update(gen_dict.keys())
                    
                edge_items = set(d.get('sequences', set()))
                if self.pg_graph.has_edge(v, u):
                    edge_items.update(self.pg_graph.edges[v, u].get('sequences', set()))
                    
                shared_items = u_items.intersection(v_items)
                drop_outs = shared_items - edge_items
                
            if drop_outs:
                    
                    # If this edge carries the majority of the shared genomes, it is the Highway.
                    # Do not impugn
                    edge_support = len(edge_items)
                    shared_support = len(shared_items)
                    
                    if edge_support > (shared_support / 2.0):
                        continue

                    true_rearrangements = set()
                    cnv_detours = set()
                    assembly_gaps = set()
                    
                    # Track the physical nature of the divergence
                    is_intra_contig = False
                    is_inter_contig = False
                    
                    for drop_gen in drop_outs:
                        u_term = target_is_terminal_in_node(drop_gen, u)
                        v_term = target_is_terminal_in_node(drop_gen, v)
                        
                        if u_term and v_term:
                            assembly_gaps.add(drop_gen)
                            continue
                            
                        # --- PRECISE CNV DETOUR LOGIC ---
                        is_detour = False
                        if u_cnv != 0 and u_cnv == v_cnv:
                            is_detour = True
                        if not is_detour and v_cnv != 0:
                            for neighbor in UG.neighbors(u):
                                if self.pg_graph.nodes[neighbor].get('cnv_cluster_id', 0) == v_cnv:
                                    if drop_gen in get_edge_genomes(u, neighbor):
                                        is_detour = True; break
                                        
                        if not is_detour and u_cnv != 0:
                            for neighbor in UG.neighbors(v):
                                if self.pg_graph.nodes[neighbor].get('cnv_cluster_id', 0) == u_cnv:
                                    if drop_gen in get_edge_genomes(v, neighbor):
                                        is_detour = True; break
                        # -------------------------------------
                        
                        if is_detour:
                            cnv_detours.add(drop_gen)
                        else:
                            true_rearrangements.add(drop_gen)
                            # DYNAMIC TOPOLOGY CHECK: Are U and V on the same contig for this genome?
                            u_contigs = set(self.pg_graph.nodes[u]['features'][drop_gen].keys())
                            v_contigs = set(self.pg_graph.nodes[v]['features'][drop_gen].keys())
                            if u_contigs & v_contigs:
                                is_intra_contig = True
                            else:
                                is_inter_contig = True
                            
                    # Apply final tags based on precisely parsed dropouts
                    if true_rearrangements:
                        assigned_jct = False
                        
                        if u_node.get("alternative_path") or v_node.get("alternative_path"):
                            d["alt_path_junction"] = True
                            assigned_jct = True
                        if u_node.get("synteny_breakpoint") or v_node.get("synteny_breakpoint"):
                            d["breakpoint_junction"] = True
                            assigned_jct = True
                        if u_node.get("repeat_ambiguity_fragmentation") or v_node.get("repeat_ambiguity_fragmentation"):
                            d["repeat_fragmentation_junction"] = True
                            assigned_jct = True
                            
                        # If none of the above, apply our dynamic structural divergence flags!
                        if not assigned_jct:
                            if is_inter_contig:
                                d["inter_contig_rearrangement"] = True
                            if is_intra_contig:
                                d["intra_contig_rearrangement"] = True
                            
                        d["sv_entities"] = ",".join(true_rearrangements)
                        sv_edges += 1
                    elif cnv_detours:
                        d["repeat_ambiguity_detour"] = True
                        d["sv_entities"] = ",".join(cnv_detours)
                    elif assembly_gaps:
                        d["is_scaffold_path"] = True
                        d["scaffold_entities"] = ",".join(assembly_gaps)
                        d["sv_entities"] = ",".join(assembly_gaps)

        inversion_events = nx.number_connected_components(inv_graph) if len(inv_graph) > 0 else 0
        logging.warning(f"SV Detection Context [{self.context}]: Tagged {inversion_events} inverted components and {sv_edges} path drop-outs.")
        return inversion_events, sv_edges
    
    def annotate_major_blocks_tfs(self, min_node_fraction=0.05):
        """
        Identifies major blocks using a Stack-Based Consensus Walk.
        Threads through Paralog Hubs by tracking physical sequence (contig) momentum.
        Assigns Block IDs to edges, and aggregates them on nodes.
        """
        # We walk the undirected graph to allow bidirectional replicon tracing
        undirected_pg = self.get_undirected_pg_graph()
        
        # Sort edges by weight (descending) so we seed the heaviest Highways first
        edges_sorted = sorted(undirected_pg.edges(data=True), key=lambda x: x[2].get('weight', 0), reverse=True)
        
        # Track unvisited edges (using sorted tuple keys to avoid directionality issues)
        unvisited_edges = {(min(u, v), max(u, v)): d for u, v, d in edges_sorted}
        
        block_id = 1
        block_manifest = []
        
        # Initialize default edge and node blocks on the Directed Graph
        for u, v in self.pg_graph.edges():
            self.pg_graph.edges[u, v]['block'] = "Fragment"
        for n in self.pg_graph.nodes():
            self.pg_graph.nodes[n]['block'] = set()
            
        total_nodes = self.pg_graph.number_of_nodes()
        min_nodes_req = total_nodes * min_node_fraction
        
        while unvisited_edges:
            # 1. Pop the heaviest unvisited edge to seed a new Block
            seed_edge, seed_data = next(iter(unvisited_edges.items()))
            del unvisited_edges[seed_edge]
            
            # RULE 2: Never seed on a known structural junction
            junction_flags = {'alt_path_junction', 'breakpoint_junction', 'repeat_fragmentation_junction', 'alternative_path',
                              'genomic_rearrangement', 'intra_contig_rearrangement', 'inter_contig_rearrangement'}
            if any(seed_data.get(flag) for flag in junction_flags):
                continue
                
            current_block = f"Block_{block_id}"
            current_nodes = set([seed_edge[0], seed_edge[1]])
            current_edges = [seed_edge]
            
            momentum_key = 'genomes' if self.context == 'feature' else 'sequences'
            momentum_set = set(seed_data.get(momentum_key, set()))
            
            stack = [(seed_edge[0], momentum_set), (seed_edge[1], momentum_set)]
            
            # 2. TFS Stack Traversal
            while stack:
                curr_node, momentum = stack.pop()
                
                # Grab all flags active on the current node
                curr_node_flags = self.pg_graph.nodes[curr_node]
                
                for nxt in undirected_pg.neighbors(curr_node):
                    edge_tuple = (min(curr_node, nxt), max(curr_node, nxt))
                    
                    if edge_tuple in unvisited_edges:
                        nxt_data = unvisited_edges[edge_tuple]
                        
                        # RULE 2: Never follow a structural junction edge
                        if any(nxt_data.get(flag) for flag in junction_flags):
                            continue
                            
                        # Rule A: Check Momentum Overlap
                        edge_threads = set(nxt_data.get(momentum_key, set()))
                        overlap = momentum.intersection(edge_threads)
                        
                        # If there is physical continuity (overlap > 0), thread through!
                        if len(overlap) > 0:
                            del unvisited_edges[edge_tuple]
                            current_edges.append(edge_tuple)
                            current_nodes.add(nxt)
                            
                            # RULE 1: MOMENTUM LOCK at Conflict Nodes
                            # Do not pick up new threads at a known evolutionary junction!
                            if any(curr_node_flags.get(flag) for flag in ["synteny_breakpoint", "repeat_ambiguity_fragmentation", "alternative_path"]):
                                new_momentum = momentum
                            else:
                                # Normal accumulation on safe backbone nodes
                                new_momentum = momentum.union(edge_threads)
                                
                            stack.append((nxt, new_momentum))
            
            # 3. Evaluate Block Size
            if len(current_nodes) >= min_nodes_req:
                final_block_label = current_block
                block_manifest.append(current_block)
                block_id += 1
            else:
                final_block_label = "Fragment"
                if "Fragment" not in block_manifest:
                    block_manifest.append("Fragment")
                    
            # 4. Apply Labels to the DIRECTED graph
            for u, v in current_edges:
                # Apply to both possible directions in the directed graph
                if self.pg_graph.has_edge(u, v):
                    self.pg_graph.edges[u, v]['block'] = final_block_label
                if self.pg_graph.has_edge(v, u):
                    self.pg_graph.edges[v, u]['block'] = final_block_label
                    
            # Nodes accumulate the blocks of their incident edges
            for n in current_nodes:
                if final_block_label != "Fragment":
                    self.pg_graph.nodes[n]['block'].add(final_block_label)

        # 5. Clean up node 'block' attributes for GEXF export
        for n, d in self.pg_graph.nodes(data=True):
            blocks = d.get('block', set())
            if not blocks:
                d['block'] = "Fragment"
            else:
                # e.g., "Block_1,Block_2" for Paralog Hubs!
                d['block'] = ",".join(sorted(list(blocks)))

        logging.warning(f"Block Detection (TFS Momentum): Found {block_id - 1} major blocks.")
        return block_manifest
    
    def export_gfa(self, gfa_file):
        """
        Exports the PS-graph to GFA v1.1 format.
        Nodes -> Segments (S), Edges -> Links (L), Sequences -> Walks (W).
        """
        logging.warning(f"Exporting GFA v1.1 to {gfa_file}")
        with open(gfa_file, 'w') as out:
            # Header
            out.write("H\tVN:Z:1.1\n")
            
            # 1. Write Segments (Nodes)
            for n, d in self.pg_graph.nodes(data=True):
                # In finalizeGraphAttr: 'family' holds the ID, 'label' holds the function string
                fam_id = d.get('family', 'unknown_fam')
                function_desc = d.get('label', 'hypothetical protein')
                dv = d.get('diversity', 0.0)
# Dynamically bundle all active boolean flags into a single comma-separated GFA tag
                node_bools = [k for k in ["repeat_ambiguity", "assembly_repeat_break", "synteny_breakpoint", "repeat_ambiguity_fragmentation", "alternative_path", "is_scaffold_path", "inverted_block"] if d.get(k)]
                nc_tag = f"\tnc:Z:{','.join(node_bools)}" if node_bools else ""
                al = d.get('alternate', 0)
                cnv = d.get('cnv_cluster_id', 0)

                node_len = d.get('length', 1)

                # Make a Bandage-friendly, space-free Segment Name (e.g., "0_PGF_12345")
                clean_fam_id = str(fam_id).replace(" ", "_")
                segment_name = f"{n}_{clean_fam_id}"
                
                # Tags: fm = family ID, fn = function (spaces are allowed in Z tags), LN = visual length
                tags = f"fm:Z:{fam_id}\tfn:Z:{function_desc}\tdv:f:{dv:.4f}\tcv:i:{cnv}\tLN:i:{node_len}{nc_tag}"
                al = d.get('alternate', 0)                
                out.write(f"S\t{segment_name}\t*\t{tags}\n")
                
            # 2. Write Links (Edges)
            for u, v, d in self.pg_graph.edges(data=True):
                seq_count = len(d.get('sequences', set()))
                gen_count = len(d.get('genomes', set()))
                
                # Format the u and v names to match the new Segment names
                u_node = self.pg_graph.nodes[u]
                v_node = self.pg_graph.nodes[v]
                u_fam = str(u_node.get('family', 'unknown')).replace(" ", "_")
                v_fam = str(v_node.get('family', 'unknown')).replace(" ", "_")
                
                u_name = f"{u}_{u_fam}"
                v_name = f"{v}_{v_fam}"
                
                # Assign GFA tags for overlapping SVs
                sv_tags = ""
                if d.get("inverted_block"):
                    sv_tags += "\tib:i:1"
                
                # Dynamically bundle all active boolean flags into a single comma-separated GFA tag
                edge_bools = [k for k in ['inter_contig_rearrangement', 'intra_contig_rearrangement', "alt_path_junction", "breakpoint_junction", "repeat_fragmentation_junction", "repeat_ambiguity_detour", "is_scaffold_path", "alternative_path"] if d.get(k)]
                if edge_bools:
                    sv_tags += f"\tjt:Z:{','.join(edge_bools)}"
                    
                # Output topology
                out.write(f"L\t{u_name}\t+\t{v_name}\t+\t0M\twc:i:{seq_count}\tgc:i:{gen_count}{sv_tags}\n")
                
            # 3. Write Walks (Genome Paths)
            for genome_id, contigs in self.replicon_map.items():
                for contig_id, features in contigs.items():
                    path =[]
                    prev_pg_id = None
                    for f_id in features:
                        pg_id = self.feature_index[f_id].pg_assignment
                        if pg_id is None:
                            continue
                            
                        # Format the path segment to match the new Segment names
                        pg_node = self.pg_graph.nodes[pg_id]
                        fam = str(pg_node.get('family', 'unknown')).replace(" ", "_")
                        segment_name = f"{pg_id}_{fam}"
                        
                        orient = ">"
                        if prev_pg_id is not None:
                            # Strand detection heuristic for inversion
                            if not self.pg_graph.has_edge(prev_pg_id, pg_id) and self.pg_graph.has_edge(pg_id, prev_pg_id):
                                orient = "<"
                                
                        path.append(f"{orient}{segment_name}")
                        prev_pg_id = pg_id
                        
                    if path:
                        path_str = "".join(path)
                        # W format: Sample Hap Index SeqID Start End Path
                        out.write(f"W\t{genome_id}\t1\t{contig_id}\t0\t{len(path)}\t{path_str}\n")

    #called after all features have been processed so that can flip them
    #and so that those that contain an anchor node can be added the anchor_instance_keys lookup
    def finalizeInstanceKeys(self):
        f=0
        while f < len(self.feature_index):
            if self.feature_index[f].instance_key !=None: #if its None then the contig wasn't big enough to create a node
                lv = [int(i) for i in self.feature_index[f].instance_key.split(".")]
                instance_reverse, instance_palindrome = self.flipKmer(lv)
                #the function actually flips the kmer
                #if instance_reverse:
                #    lv.reverse()
                instance_key_orig = ".".join(str(i) for i in lv)
                self.feature_index[f].instance_key = instance_key_orig+"."+str(self.feature_index[f].repeat_num)+"r"
                # if the instance key has a single anchor node in it then it is an anchor instance key
                if instance_key_orig in self.anchor_instance_keys:
                    self.anchor_instance_keys[instance_key_orig][1].append(f)
                else:
                    #check to see if any of the instances are anchor nodes
                    for rf in lv:
                        if self.rf_node_index[rf].anchorNode():
                            self.anchor_instance_keys[instance_key_orig]=[None,[f]]
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
        if self.context != "feature":
            self.context_bin.add(kmer_key)
        self.rf_graph.add_node(self.cur_rf_node.nodeID, label=kmer_key, duplicate=dup_number)
        if not duplicate:
            self.cur_rf_node.numBins+=1
        if self.traverse_priority == "area":
            self.contig_to_rf.setdefault(feature_list[0].contig_id,[]).append(self.cur_rf_node.nodeID)

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
                self.rf_graph.add_edge(self.prev_node.nodeID, self.cur_rf_node.nodeID, **{"flip":fflip,"leaving_position":leaving_position})
                self.rf_graph.add_edge(self.cur_rf_node.nodeID, self.prev_node.nodeID, **{"flip":rflip,"leaving_position":reverse_lp})
        self.prev_indices=feature_indices
        self.prev_node=self.cur_rf_node
        self.prev_reverse=reverse

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
        while i < (k_size // 2):
            if feature_list[i] < feature_list[k_size-(i+1)]:
                return (reverse, palindrome)
            elif feature_list[i] > feature_list[k_size-(i+1)]:
                reverse=1
                feature_list.reverse()
                return (reverse, palindrome)
            i+=1
        palindrome =1
        return (reverse, palindrome)


    def processFeatures(self):
        logging.warning("parsing features and constructing kmer graph")	
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
            if prev_feature == None or prev_feature.genome_id != feature.genome_id:
                self.trackDiversity(feature.feature_id, self.all_diversity)
            self.replicon_map.setdefault(feature.genome_id,OrderedDict()).setdefault(feature.contig_id,[]).append(feature.feature_id)
            if prev_feature and (prev_feature.contig_id != feature.contig_id or prev_feature.genome_id != feature.genome_id):
                kmer_q=deque()#clear kmer stack because switching replicons
                self.prev_node=None
                self.prev_indices=[]
            elif prev_feature and prev_feature.contig_id == feature.contig_id:
                if prev_feature.start > feature.start:
                    raise Exception("InputError")
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
        #sort by area that it will take up, then by whether it is earlier in the contig
        if self.traverse_priority == "area":
            #numBins will be the number of contexts the kmer will show up in (non-dup). this number will be inflated but correct for relative sorting.
            self.contig_weight={contig:sum([self.rf_node_index[r].numBins for r in rf_list]) for contig,rf_list in self.contig_to_rf.items()}
            temp_starting_list = [(self.rfNodeToMaxAlignArea(i.nodeID, start_tuple=True), i) for i in self.rf_starting_list] #make a list of ((weight, start), rfnode) tuples
            temp_starting_list.sort(key=lambda x: x[0][1])#first we sort by start so that lower coordinates come first (as a secondary sort criteria)
            temp_starting_list.sort(key=lambda x: x[0][0], reverse=True)#second we sort by area/weight so taht higher values come first (as a primary sort criteria)
            self.rf_starting_list=[i[1] for i in temp_starting_list]
        else:
            #start with nodes that have the most features
            self.rf_starting_list.sort(key=lambda x: x.numFeatures(), reverse=True)
        for rf_node in self.rf_starting_list:
            if rf_node.numFeatures() > 0:
                self.tfs_expand_nr(None, rf_node, None, None)

    def rfNodeToMaxAlignArea(self, rf_target_id, node_bundle=None, start_tuple=False):
        cur_weights=[]
        if node_bundle != None:
            features= sum([list(i) for i in sum(node_bundle[rf_target_id],[])], []) #a bunch of sets of features in the node bundle.
        else:
            features=sum([list(i) for i in self.rf_node_index[rf_target_id].features],[]) # get concatenated list of features
        for f in features:
            feature_obj = self.feature_index[f]
            cur_weights.append((self.contig_weight[feature_obj.contig_id], feature_obj.start))
        if start_tuple:
            cur_weights.sort(key=itemgetter(1))#first we sort by start so that lower coordinates come first (as a secondary sort criteria)
            cur_weights.sort(key=itemgetter(0), reverse=True)#second we sort by area/weight so that higher values come first (as a primary sort criteria)
            return cur_weights[0]
        else:
            return max(cur_weights, key=itemgetter(0))

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
                logging.warning("logic problem. calculated leaving side does not match")
                raise Exception("LogicError")
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
                    raise Exception("LogicError")
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
    def queueFeature(self, cur_node, kmer_side, orientation, leaving_feature, temp_self_queue, temp_reg_queue, node_bundles, up_node=None):
        nxt_rf_id = self.nextRFNode(kmer_side, orientation, leaving_feature)
        prev_queued=True #has the rfid EVER been queued on THIS traversal
        up_queue= up_node!=None and nxt_rf_id == up_node.nodeID
        #If there are no unassigned features in the nxt node is there any point in visiting? conflict detection etc.?
        if nxt_rf_id!=None and self.rf_node_index[nxt_rf_id].numFeatures() > 0 and self.rf_node_index[nxt_rf_id].numBins > self.minSeq:
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
                        temp_self_queue.append((nxt_rf_id,[nxt_guide, nxt_guide_cat, nxt_guide_side]))
                    else:
                        temp_reg_queue.append((nxt_rf_id,[nxt_guide, nxt_guide_cat, nxt_guide_side]))
                #node bundles exist separate from queue but are cleared out when rfid is taken from the queue
                node_bundles[bundle_id][nxt_position][nxt_direction].add(nxt_target)
            return nxt_rf_id
        return None


    def num_features_pg_node(self, node_id):
        num_features=0
        for g in self.pg_graph.nodes[node_id]['features']:
            for c in self.pg_graph.nodes[node_id]['features'][g]:
                num_features+=len(self.pg_graph.nodes[node_id]['features'][g][c])
        return num_features


    def construct_pg_edge(self, prev_pg_id, cur_pg_id, genome_id, sequence_id, strand='+'):
            # Determine the unit based on Context
            context_key = genome_id if self.context == "genome" else sequence_id
            
            edge_data = self.pg_graph.get_edge_data(prev_pg_id, cur_pg_id, default=None)
            
            if edge_data is not None:
                # Build the NxN Similarity Matrix for the Context Unit
                existing_entities = edge_data["genomes"] if self.context == "genome" else edge_data["sequences"]
                
                for existing_key in existing_entities:
                    if existing_key != context_key: # Skip self-loops
                        self.similarity_matrix[context_key][existing_key] += 1
                        self.similarity_matrix[existing_key][context_key] += 1
                
                edge_data["genomes"].add(genome_id)
                edge_data["sequences"].add(sequence_id)
                edge_data.setdefault("strands", set()).add(strand) # <-- Store the consensus strand!
            else:
                self.pg_graph.add_edge(prev_pg_id, cur_pg_id, genomes={genome_id}, sequences={sequence_id}, strands={strand})



    def get_undirected_pg_graph(self):
        """
        Safely converts the directed PS-graph to an undirected graph for GEXF export.
        Merges bidirectional edges (like inversions) and unifies their attributes.
        """
        undirected = nx.Graph()
        undirected.add_nodes_from(self.pg_graph.nodes(data=True))
        
        num_genomes = float(len(self.replicon_map.keys()))
        
        for u, v, d in self.pg_graph.edges(data=True):
            if undirected.has_edge(u, v):
                existing = undirected.edges[u, v]
                
                # Merge Genomes
                gen_set = set(existing.get('genomes', set())) | set(d.get('genomes', set()))
                gen_set.discard('')  # Remove empty strings
                existing['genomes'] = gen_set

                # Merge Sequences
                seq_set = set(existing.get('sequences', set())) | set(d.get('sequences', set()))
                seq_set.discard('')
                existing['sequences'] = seq_set
                                
                # Recalculate true combined weight
                existing['weight'] = len(gen_set) / num_genomes if num_genomes > 0 else 0.0
                
                # Merge SV Flags (if either direction was an inversion/translocation, the undirected edge is too)
                edge_bool_flags = ['inter_contig_rearrangement', 'intra_contig_rearrangement', 'inverted_block', 'alt_path_junction', 'breakpoint_junction', 'repeat_fragmentation_junction', 'repeat_ambiguity_detour', 'is_scaffold_path', 'alternative_path', 'is_superbubble']
                for flag in edge_bool_flags:
                    existing[flag] = existing.get(flag, False) or d.get(flag, False)
                if 'superbubble_id' in existing or 'superbubble_id' in d:
                    existing['superbubble_id'] = set(existing.get('superbubble_id', set())) | set(d.get('superbubble_id', set()))
                
            else:
                # Add new edge (make a copy of the dictionary to avoid mutating original)
                attr_dict = {k: set(v) if isinstance(v, set) else v for k, v in d.items()}
                undirected.add_edge(u, v, **attr_dict)
                
        return undirected

    def insert_feature(self, cur_pg_id, new_feature):
        #emit_extra=False
        md5=self.feature_index[new_feature].md5
        start=self.feature_index[new_feature].start
        end=self.feature_index[new_feature].end
#         group_id=self.feature_index[new_feature].group_id
        genome_id=self.feature_index[new_feature].genome_id
        sequence_id=self.feature_index[new_feature].contig_id
        if not genome_id in self.pg_graph.nodes[cur_pg_id]['features']:
            self.pg_graph.nodes[cur_pg_id]['features'][genome_id]={sequence_id:[new_feature]}
            if embed_info:
                self.pg_graph.nodes[cur_pg_id]['features']['info'][genome_id]={sequence_id:[{'md5':md5, 'start':start, 'end':end}]}
            insert_level="genome"
        elif not sequence_id in self.pg_graph.nodes[cur_pg_id]['features'][genome_id]:
            self.pg_graph.nodes[cur_pg_id]['features'][genome_id][sequence_id]=[new_feature]
            if embed_info:
                self.pg_graph.nodes[cur_pg_id]['features']['info'][genome_id][sequence_id]=[{'md5':md5, 'start':start, 'end':end}]
            insert_level="contig"
        else:
            #if self.context!="feature":
                #for cf in self.pg_graph.nodes[cur_pg_id]['features'][genome_id][sequence_id]:
                #    #if the distance is < k it is a special case of an 'extra character loop' which requires emitting an extra pg-node
                #    if abs(cf-new_feature)< self.ksize:
                #        emit_extra=True
                #        insert_level=self.context #so there won't be a problem
                #        break
            #if not emit_extra:
                self.pg_graph.nodes[cur_pg_id]['features'][genome_id][sequence_id].append(new_feature)
                if embed_info:
                    self.pg_graph.nodes[cur_pg_id]['features']['info'][genome_id][sequence_id].append({'md5':md5, 'start':start, 'end':end})
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
        if not genome_id in self.pg_graph.nodes[cur_pg_id]['features']:
            insert_level="genome"
        elif not sequence_id in self.pg_graph.nodes[cur_pg_id]['features'][genome_id]:
            insert_level="contig"
        else:
            insert_level="feature"
            split =True
        if self.context_levels[insert_level] > self.context_levels[self.context]:
            conflict=True
            cf = next(iter(self.pg_graph.nodes[cur_pg_id]['features'][genome_id].values()))[0]
            if not split:
                if self.debug:
                    logging.info("conflict between "+str(new_feature)+" and "+str(cf)+" in "+str(cur_pg_id))
                logging.info("conflict between "+self.feature_index[new_feature].feature_ref+" and "+self.feature_index[cf].feature_ref+" in "+str(cur_pg_id))
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
            #genome_id in self.pg_graph.nodes[cur_pg_id]['features'] and \
            #sequence_id in self.pg_graph.nodes[cur_pg_id]['features'][genome_id]:
            #    for cf in self.pg_graph.nodes[cur_pg_id]['features'][genome_id][sequence_id]:
                    #if the distance is < k it is a special case of an 'extra character loop' which requires emitting an extra pg-node
            #        if abs(cf-new_feature)< self.ksize:
            #            split=True
            
        return split, conflict, end_fragments

    def assign_pg_node(self, prev_feature, new_feature, pg_node, conflicts=None, relative_strand='+'):
        cur_pg_id=None
        #determine if there is a conflict based on mixed bundling
#         group_id=self.feature_index[new_feature].group_id
        genome_id=self.feature_index[new_feature].genome_id
        sequence_id=self.feature_index[new_feature].contig_id
        md5=self.feature_index[new_feature].md5
        start=self.feature_index[new_feature].start
        end=self.feature_index[new_feature].end
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
            if self.order_contigs:
                self.contig_order.setdefault(genome_id,OrderedDict()).setdefault(sequence_id, []).append(new_feature)
            if pg_node!=None:
                #REPLACED insert_level=None
                cur_pg_id=pg_node
                insert_level = self.insert_feature(cur_pg_id, new_feature)
                #REPLACED if self.context_levels[insert_level] > self.context_levels[self.context]:
                #REPLACED    conflict=True

            else:
                cur_pg_id=self.num_pg_nodes
                self.num_pg_nodes+=1
#                 self.pg_graph.add_node(cur_pg_id, label=str(self.feature_index[new_feature].group_num), features={genome_id:{sequence_id:[new_feature]}, 'md5':md5, 'start':start, 'end':end})
                if embed_info:
                    self.pg_graph.add_node(cur_pg_id, label=str(self.feature_index[new_feature].group_num), features={genome_id:{sequence_id:[new_feature]}, 'info':{genome_id:{sequence_id:[{'md5':md5, 'start':start, 'end':end}]}}})
                else:
                    self.pg_graph.add_node(cur_pg_id, label=str(self.feature_index[new_feature].group_num), features={genome_id:{sequence_id:[new_feature]}})
            self.feature_index[new_feature].pg_assignment=cur_pg_id
        
        
        if conflicts:
            for c in conflicts:
                self.pg_graph.nodes[cur_pg_id][c] = True
            # else where for more comprehensive marking
            #if conflict == "shift" and not "alternate" in self.pg_graph.nodes[cur_pg_id]:
            #    self.pg_graph.nodes[cur_pg_id]["alternate"]=1

        node_strands = self.pg_graph.nodes[cur_pg_id].setdefault("strands", set())
        node_strands.add(relative_strand)
        if len(node_strands) > 1:
            self.pg_graph.nodes[cur_pg_id]["inverted_block"] = True

        if prev_feature != None :
            prev_pg_id = self.feature_index[prev_feature].pg_assignment
            self.construct_pg_edge(prev_pg_id, cur_pg_id, genome_id, sequence_id, strand=relative_strand)
            #if conflict:
            #    tier2=self.conflicts.setdefault(cur_pg_id, {})
            #    conflicted=tier2.setdefault(prev_pg_id, set([]))
            #    conflicted.add(new_feature)
        return cur_pg_id


    def maxInstanceOverlap(self, instance_key, comparison_list):
        instance_parts=set(instance_key.split('.'))
        max_intersect = 0
        max_key = []
        for c in comparison_list:
            cur_parts=set(c.split('.'))
            cur_intersect = len(instance_parts.intersection(cur_parts))
            if cur_intersect > max_intersect:
                max_intersect = cur_intersect
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
        while i >= 0:
            #conflicts can only happen for pre-existing nodes.
            for pg in pre_assignments[i]["assignments"]:
                num_regular =0
                num_c1 =0
                for ik, features in pre_assignments[i]["assignments"][pg]["features"].items():
                    for f in features:
                        #think about how the conflict should be represented overall
                        #look at the data structure required for determineAssignments
                        #start with that.
                        shift, conflict, end_fragments = self.detect_conflict(f.new_feature, pg)
                        if shift:
                            conflict_assignments[i]["shift"].setdefault(pg, {'features':{}})
                            conflict_assignments[i]["shift"][pg]["features"].setdefault(ik,set([])).update(features)
                            #self.no_edge.add(pg)
                        elif end_fragments:
                            conflict_assignments[i]["c2_conflict"].setdefault(pg, {'features':{}})
                            conflict_assignments[i]["c2_conflict"][pg]["features"].setdefault(ik,set([])).update(features)
                            #self.no_edge.add(pg)
                        elif conflict:
                            conflict_assignments[i]["c1_conflict"].setdefault(pg, {'features':{}})
                            conflict_assignments[i]["c1_conflict"][pg]["features"].setdefault(ik,set([])).update(features)
                            num_c1+=1
                            #self.no_edge.add(pg)
                        else:
                            num_regular +=1

            i-=1            


    #considering available guides/pg-nodes and conflicts to determine round of assignments
    #put new node assignments under None
    def determineAssignments(self, feature_pile, pre_assignments, cur_node, conflicts=None):
        if conflicts == None:
            i=self.ksize-1
            while i >= 0:
                for ik, features in feature_pile[i]["by_instance"].items():
                    # 1. Separate guided features from unguided features
                    guided_assignments = collections.defaultdict(set)
                    unguided_features = set()
                    
                    max_keys = self.maxInstanceOverlap(ik, pre_assignments[i]["inst_key_map"].keys())
                    valid_pgs = set()
                    if max_keys[0] is not None:
                        for k in max_keys:
                            mapped_pgs = pre_assignments[i]["inst_key_map"][k]
                            if isinstance(mapped_pgs, set):
                                valid_pgs.update(mapped_pgs)
                            else:
                                valid_pgs.add(mapped_pgs)
                            
                    for f in features:
                        f_signed = self.sign_feature(f.new_feature, f.direction)
                        assigned_to_guide = False
                        
                        # Check if this specific feature is a known guide for any valid_pg
                        for pg in valid_pgs:
                            if pg in pre_assignments[i]["assignments"] and ik in pre_assignments[i]["assignments"][pg]["guides"]:
                                if f_signed in pre_assignments[i]["assignments"][pg]["guides"][ik]:
                                    guided_assignments[pg].add(f)
                                    assigned_to_guide = True
                                    break
                            if pg in pre_assignments[i]["new_nodes"] and ik in pre_assignments[i]["new_nodes"][pg]["guides"]:
                                if f_signed in pre_assignments[i]["new_nodes"][pg]["guides"][ik]:
                                    guided_assignments[pg].add(f)
                                    assigned_to_guide = True
                                    break
                        
                        if not assigned_to_guide:
                            unguided_features.add(f)
                            
                    # 2. Process guided features (they go exactly to their guided PGs)
                    for pg, sub_features in guided_assignments.items():
                        block = "assignments" if pg in pre_assignments[i]["assignments"] else "new_nodes"
                        cur_assignments = pre_assignments[i][block][pg]["features"].setdefault(ik, set())
                        cur_assignments.update(sub_features)
                        
                    # 3. Process unguided features (close SNP bubbles or emit new node)
                    if unguided_features:
                        if valid_pgs:
                            max_f = -1
                            best_pg = None
                            for pg in valid_pgs:
                                n_f = self.num_features_pg_node(pg) if type(pg) != str else 0
                                if n_f > max_f:
                                    max_f = n_f
                                    best_pg = pg
                            block = "assignments" if best_pg in pre_assignments[i]["assignments"] else "new_nodes"
                            if ik not in pre_assignments[i]["inst_key_map"]:
                                guide_feature = next(iter(unguided_features))
                                self.track_guide(pre_assignments, guide_feature.new_feature, i, negative=guide_feature.direction, target_pg=best_pg, new_node=(block == "new_nodes"))
                            cur_assignments = pre_assignments[i][block][best_pg]["features"].setdefault(ik, set())
                            cur_assignments.update(unguided_features)
                        else:
                            # No guides exist at all! Create a new node for them.
                            tmp_id = "nn"+str(len(pre_assignments[i]["new_nodes"].keys()))
                            guide_feature = next(iter(unguided_features))
                            self.track_guide(pre_assignments, guide_feature.new_feature, i, negative=guide_feature.direction, target_pg=tmp_id, new_node=True)
                            
                            cur_assignments = pre_assignments[i]["new_nodes"][tmp_id]["features"].setdefault(ik, set())
                            cur_assignments.update(unguided_features)
                i-=1
        else:
            i=self.ksize-1
            while i >= 0:
                #if there is a shift conflict create a new node. only need the pg_node and the offending pre_assignment instance key
                if "shift" in conflicts[i] and len(conflicts[i]["shift"]) > 0:
                    for pg in conflicts[i]["shift"]:
                        for ik in conflicts[i]["shift"][pg]["features"]:
                            self.alt_counter+=1
                            tmp_id = "nn"+str(len(pre_assignments[i]["new_nodes"].keys()))
                            pre_assignments[i]["new_nodes"].setdefault(tmp_id, {'guides':{},'features':{}})
                            pre_assignments[i]["new_nodes"][tmp_id]["features"][ik]=conflicts[i]["shift"][pg]["features"][ik]
                            #take out features of pre_assignments
                            pre_assignments[i]["assignments"][pg]["features"].pop(ik)
                i-=1

    
    #store alternate pg-nodes at a given k-position
    def storeAlternates(self, pg_ids):
        assert type(pg_ids == set)
        for i in pg_ids:
            self.pg_node_alt.setdefault(i, set([])).update(pg_ids)

    #make assignments that are in pre_assignments
    def makeAssignments(self, pre_assignments, conflicts):
        i=self.ksize-1
        while i >= 0:
            alt_nodes = set([])
            alt_nodes.update(pre_assignments[i]["assignments"].keys())
            for pg in pre_assignments[i]["assignments"]:
                conflict_status = []
                if ("shift" in conflicts[i] and pg in conflicts[i]["shift"]):
                    conflict_status.append("repeat_ambiguity")
                if ("c2_conflict" in conflicts[i] and pg in conflicts[i]["c2_conflict"]):
                    conflict_status.append("assembly_repeat_break")
                if ("c1_conflict" in conflicts[i] and pg in conflicts[i]["c1_conflict"]):
                    conflict_status.append("synteny_breakpoint")

                for ik, features in pre_assignments[i]["assignments"][pg]["features"].items():
                    for f in features:
                        raw_strand = self.feature_index[f.new_feature].strand
                        if f.direction == 1:
                            rel_strand = '+' if raw_strand == '-' else '-'
                        else:
                            rel_strand = raw_strand
                            
                        self.assign_pg_node(prev_feature=f.prev_feature, new_feature=f.new_feature, pg_node=pg, conflicts=conflict_status, relative_strand=rel_strand)
            for pg in pre_assignments[i]["new_nodes"]:
                new_pg = None
                for ik, features in pre_assignments[i]["new_nodes"][pg]["features"].items():
                    for f in features:
                        # Translate the raw strand based on the TFS traversal direction!
                        raw_strand = self.feature_index[f.new_feature].strand
                        if f.direction == 1:
                            rel_strand = '+' if raw_strand == '-' else '-'
                        else:
                            rel_strand = raw_strand
                            
                        new_pg=self.assign_pg_node(prev_feature=f.prev_feature, new_feature=f.new_feature, pg_node=new_pg, conflicts=[], relative_strand=rel_strand)
                        alt_nodes.add(new_pg)
            if len(alt_nodes) > 1: self.storeAlternates(alt_nodes)
            i-=1

    #walk down the kmer and fill out to_assign
    #process features remaining in cur_node.features
    #targets organized as targets["left" & "right" == 0 & 1][ "increasing" & "decreasing" == 0 & 1 ]
    #targets are all RHS representative
    def kfill_feature_pile(self, feature_pile, pre_assignments, cur_node, prev_node, temp_self_queue, temp_reg_queue, node_bundles, targets):
            i=0
            while i < self.ksize: #because any features remaining represent "new" threads need to assign the entire k-mer
                for direction in self.target_cat:
                    if targets != None:
                        feature_sets = [("old", self.ksize-1 ,targets[0][direction]), ("old", 0,targets[1][direction]), ("new",None,cur_node.features[direction])]
                    else: feature_sets = [("new",None,cur_node.features[direction])]
                    for status, skip_count, to_fill in feature_sets:
                        for rhs_feature in to_fill:
                            new_feature_adj= i
                            prev_feature_adj=i-1
                            #TODO Consider adjusting this based on signed features.
                            #TODO Consider what would be necessary to adapt this to consistency use for palindrome
                            if not direction:
                                new_feature_adj=new_feature_adj*-1
                                prev_feature_adj=prev_feature_adj*-1
                            new_feature=rhs_feature+new_feature_adj
                            instance_key = self.feature_index[new_feature].instance_key
                            prev_feature= rhs_feature+prev_feature_adj if i > 0 else None 
                            if status == "new":
                                #cur_node.features[direction].remove(rhs_feature)
                                cur_node.assigned_features[direction].add(rhs_feature) #only track rhs.
                                #there are two cases where the feature could be about to leave the kmer-frame. If they are on the left or right of the kmeri
                                if i == 0:
                                    #this feature is on the rhs of kmer
                                    self.queueFeature(cur_node, 1, direction, new_feature, temp_self_queue, temp_reg_queue, node_bundles, up_node=prev_node)
                                if i == self.ksize-1:
                                    #this feature is on the lhs of kmer
                                    self.queueFeature(cur_node, 0, direction, new_feature, temp_self_queue, temp_reg_queue, node_bundles, up_node=prev_node)
                            #only want to fill features/track guides that are from new (non-targets) and from the parts 
                            if status == "new" or i != skip_count:
                                if self.feature_index[new_feature].pg_assignment != None:
                                    self.track_guide(pre_assignments, new_feature, self.ksize-(i+1), negative=direction)
                                if status != "old":
                                    #track feature to "add" in case of anchor instance key expansion AND/OR pg-edge addition
                                    feature_info=self.featureContext(self.ksize-(i+1), direction, rhs_feature, prev_feature, new_feature, leaving_feature=None, conflict=False, split=False)
                                    self.track_feature(self.ksize-(i+1), feature_pile, new_feature, feature_info, instance_key)
                i+=1

    def getFeatureRecord(self, feature):
        return self.feature_index[abs(feature)]

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
                if self.getFeatureRecord(cur_sf).pg_assignment != None:
                    self.track_guide(pre_assignments, abs(cur_sf), count, negative=(cur_sf < 0)) 
            for sf in rhs:
                cur_sf = sf - count
                if self.getFeatureRecord(cur_sf).pg_assignment != None:
                    self.track_guide(pre_assignments, abs(cur_sf), (self.ksize-1)-count, negative=(cur_sf < 0))
            count+=1


    #keep track of which nodes are available for assignment in this 'column' of the kmer
    def track_guide(self, pre_assignments, unsigned_guide, col, negative, target_pg=None, new_node=False):
        new_guide = self.sign_feature(unsigned_guide, negative)
        guide_ikey=self.getFeatureRecord(new_guide).instance_key
        if target_pg == None:
            target_pg =self.getFeatureRecord(new_guide).pg_assignment
            if target_pg == None:
                raise Exception("LogicError")
        if not new_node:
            pre_assignments[col]["assignments"].setdefault(target_pg, {'guides':{},'features':{}})
            pre_assignments[col]["assignments"][target_pg]['guides'].setdefault(guide_ikey, set([])).add(new_guide)
            pre_assignments[col]["inst_key_map"].setdefault(guide_ikey, set()).add(target_pg)
        else:
            pre_assignments[col]["new_nodes"].setdefault(target_pg, {'guides':{},'features':{}})
            pre_assignments[col]["new_nodes"][target_pg]['guides'].setdefault(guide_ikey, set([])).add(new_guide)
            pre_assignments[col]["inst_key_map"].setdefault(guide_ikey, set()).add(target_pg)
            

    def track_feature(self, col, feature_pile, new_feature, feature_info, instance_key):
        if not new_feature in feature_pile[col]['all_features']:
            feature_pile[col]['by_instance'].setdefault(instance_key, set([])).add(feature_info)
            feature_pile[col]['all_features'].add(new_feature)

    #right now kmer_side is 0 to indicate left and 1 to indicate right
    def side_to_kcoord(self, kmer_side):
        answer = kmer_side*self.ksize
        if answer > 0:
            answer-=1
        return answer

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
            existing_label = self.rf_graph.nodes[cur_node.nodeID]["visit"]
            self.rf_graph.nodes[cur_node.nodeID]["visit"] =  ",".join([str(self.visit_number),existing_label]) if len(existing_label) else str(self.visit_number)
            logging.debug("visiting rf-"+str(cur_node.nodeID)+" number of pg-nodes is "+str(self.pg_graph.number_of_nodes()))
        #sys.stderr.write("number of pg-nodes is "+str(self.pg_graph.number_of_nodes())+"\n")
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
                rhs_adj_info=self.rhs_adj_table[new_guide_cat][new_guide_side]
                #the new guide is signed so can do simple increment/decrement depending on LHS or RHS. always take abs() to get actual ID
                new_guide= rhs_guide + rhs_adj_info['new_feature_adj'] #guides are aligned to the right side and walked left from there
                #REPLACED instance_key = self.feature_index[new_guide].instance_key
                self.track_guide(pre_assignments, new_guide, self.side_to_kcoord(new_guide_side), negative=new_guide_cat)
                # SEE IF YOU CAN GET AWAY WITHOUT DOING THIS TILL YOU NEED IT self.fill_guides(pre_assignments, rhs_guide, rhs_adj_info['new_feature_adj'])
                #REPLACED pg_set.add(self.feature_index[new_guide].pg_assignment)
                #REPLACED target_guides[instance_key]=[new_guide]

                        
        #whether this is an anchor or not there will be targets passed down if it is not the start of a traversal.
        temp_self_queue=[]
        temp_reg_queue=[]
        if (num_targets>0):
            # if there are targets then this isn't the first node visited
            # this means only one new column aka 'character' in the kmer needs to be expanded 
            # no matter which feature is being 'targeted' only the RHS of this kmer exists in the features structure.
            #targets organized as targets["left" & "right" == 0 & 1][ "increasing" & "decreasing" == 0 & 1 ]
            #REPLACED currently_seen=set([])

            #REPLACED available_nodes=set([])
            #PROCESS TARGETS
            for kmer_side in self.target_cat:
                for direction in self.target_cat:
                    rhs_adj_info=self.rhs_adj_table[direction][kmer_side]
                    #get the targets that have been assigned
                    for rhs_feature in targets[kmer_side][direction]:
                        new_feature=rhs_feature + rhs_adj_info['new_feature_adj']
                        prev_feature= rhs_feature + rhs_adj_info['prev_feature_adj']
                        leaving_feature=rhs_feature + rhs_adj_info['leaving_feature_adj']
                        cur_instance_key = self.feature_index[new_feature].instance_key
                        edge_only=False
                        if not rhs_feature in cur_node.features[direction]:
                            if not rhs_feature in cur_node.assigned_features[direction]:
                                logging.warning("missing projected "+str(rhs_feature)+" in "+str(cur_node.nodeID)+" from "+str(prev_node.nodeID))
                                raise Exception("LogicError")
                            elif self.feature_index[new_feature].pg_assignment != None:
                                #construct pg-edge to previously created node
                                edge_only=True
                                
                                raw_strand = self.feature_index[new_feature].strand
                                if direction == 1:
                                    rel_strand = '+' if raw_strand == '-' else '-'
                                else:
                                    rel_strand = raw_strand
                                    
                                self.construct_pg_edge(self.feature_index[prev_feature].pg_assignment, self.feature_index[new_feature].pg_assignment, self.feature_index[new_feature].genome_id, self.feature_index[new_feature].contig_id, strand=rel_strand)
                            else:
                                logging.warning("pre-processed, unassigned target "+str(new_feature)+" in "+str(cur_node.nodeID)+" from "+str(prev_node.nodeID))
                                raise Exception("LogicError")
                        else:
                            #this initial loop through the targets is really just to see if any have been assigned
                            if self.feature_index[new_feature].pg_assignment != None:
                                #REPLACED rhs_guide = rhs_feature
                                #REPLACED rhs_guide_cat=direction
                                #REPLACED rhs_guide_side=kmer_side #need this for palindromes
                                #REPLACED pg_set.add(self.feature_index[new_feature].pg_assignment)
                                #REPLACED target_guides.setdefault(instance_key,[]).append(new_feature)
                                self.track_guide(pre_assignments, new_feature, self.side_to_kcoord(kmer_side), negative=direction)
                            #else:
                            cur_node.features[direction].remove(rhs_feature)
                            cur_node.assigned_features[direction].add(rhs_feature) #only track rhs.
                            cur_info=self.featureContext(kmer_side, direction, rhs_feature, prev_feature, new_feature, leaving_feature, False, False)
                            self.track_feature(self.side_to_kcoord(kmer_side), feature_pile, new_feature, cur_info, cur_instance_key)
                        if up_targets:#if up_targets is true then these features were returned from a DFS exploration of an anchor node and passed here as target.
                            #when queueing based on return 'new' features make sure don't do a DFS "up"
                            q_rfid = self.queueFeature(cur_node, (not kmer_side), direction, leaving_feature, temp_self_queue, temp_reg_queue, node_bundles, up_node=prev_node) #no prevent_node
                        elif not edge_only:
                            q_rfid = self.queueFeature(cur_node, (not kmer_side), direction, leaving_feature, temp_self_queue, temp_reg_queue, node_bundles) #no prevent_node

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



        #the only place where all k-features are added is at an anchor node if there are features that are not currently part of a bundle
        #if a "feature thread" is already in the bundle then the previous k-1 features have been processed
        #in that case the previous k-1 features should have already been checked for anchor instance expansion
        if cur_node.anchorNode():
            #if this is an anchor node and a starting node then everything needs to expanded/assigned to a pg-node
            #if this is an anchor node and had targets incoming then everything remaining is new and needs to be fully expanded
            #at this point anything remaining is regarded as 'new' and can be passed as targets up or down !!!!!
            #rhs_guide is used to track a feature "thread" that has already been assigned so that current features can be assigned to the correct pg_node
            self.kfill_feature_pile(feature_pile, pre_assignments, cur_node, prev_node, temp_self_queue, temp_reg_queue, node_bundles, targets)
            self.fill_guides(pre_assignments)
            cur_node.features[0]=set([])#after assigning all features clear it out.
            cur_node.features[1]=set([])#after assigning all features clear it out.
	
        if self.traverse_priority == "area":
            #mod these to sort based on maximum area next
            temp_self_queue.sort(key=lambda x: self.rfNodeToMaxAlignArea(x[0],node_bundles), reverse=True)
            temp_reg_queue.sort(key=lambda x: self.rfNodeToMaxAlignArea(x[0],node_bundles), reverse=True)

        node_queue.extendleft(temp_self_queue)
        node_queue.extend(temp_reg_queue)


        #anchor_expanded=False
        #anchor_expanded = self.anchorInstanceExpansion(feature_pile)

        self.determineAssignments(feature_pile, pre_assignments, cur_node, conflicts=None)
        self.findKConflicts(pre_assignments, conflict_assignments)
        self.determineAssignments(feature_pile, pre_assignments, cur_node, conflict_assignments)
        self.makeAssignments(pre_assignments, conflict_assignments)
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
    #Because of this the new target will be on the opposte side than the previous targets
    def getTargetGuide(self,targets):
        kmer_side=0
        guide=None
        while kmer_side < len(targets):
            if guide != None:
                break
            direction=0
            while direction < len(targets[kmer_side]):
                if len(targets[kmer_side][direction]) > 0:
                    guide = (next(iter(targets[kmer_side][direction])), direction, not kmer_side)
                    break
                direction+=1
            kmer_side+=1
        if guide == None:
            raise Exception("LogicError")
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
                    self.expand_features(cv.prev_node, cv.cur_node, targets=cntargets, guide_array=new_guide, node_queue=cv.node_queue, node_bundles=cv.node_bundles, up_targets=True)
            pv=cv

class pFamGraph(nx.DiGraph):
    def __init__(self):
        nx.DiGraph.__init__(self)


def main():
    parser = argparse.ArgumentParser()

    parser.set_defaults(file_type="patric_tab")
    #parser.add_argument('--break_conflict', help='Uses methods for dealing with latent updating to APIs', required=False, default=False, action='store_true')
    parser.add_argument('--no_function', help='no functions as labels. keep file size smaller.', required=False, default=False, action='store_true')
    parser.add_argument('--order_contigs', choices=["none","area","tfs"], help='produce output that orders contigs for rectilinear layout', required=False, default="none")
    parser.add_argument('--contig_output', help='output file that orders contigs for rectilinear layout', required=False, default=None)
    parser.add_argument('--layout', help='run gephi layout code for gexf', required=False, default=False, action='store_true')
    parser.add_argument("--output", type=str, help="the path and base name give to the output files. if not given goes to stdout", required=False, default=sys.stdout)
    parser.add_argument("--rfgraph", type=str, help="create rf-graph gexf file at the following location", required=False, default=None)
    parser.add_argument("--diversity", type=str, help="calculate diversity quotient according to given taxa level", required=False, default="genus", choices=["genus","species"])
    parser.add_argument("--alpha", type=str, help="alphabet i.e. 'family type' to use", required=False, default="pgfam_id", choices=["figfam_id","plfam_id","pgfam_id"])
    input_type = parser.add_mutually_exclusive_group()
    input_type.add_argument("--patric", dest="file_type", help="table specifying the group, genome, contig, feature, and start in sorted order", action='store_const', const="patric_tab")
    input_type.add_argument("--patric_genomes", dest="file_type", action='store_const', const="patric_genomes", help="use the files listed in --feature_files as a comma or tab separated file specifying genome ids to pull from patric. automatically downloads and uses the data stream for those genome ids.")
    
    parser.add_argument("--feature_selector", type=str, required=False, default="genome_id", choices=["genome_id", "accession"], help="What key to use to filter the results from a patric stream. Only has an effect if the --patric_genomes flag is also used.")
    parser.add_argument("--context", type=str, required=False, default="genome", choices=["genome","contig","feature"], help="the synteny context")
    parser.add_argument("--ksize", type=int, default=3, required=False, choices=range(3,10), help="the size of the kmer to use in constructing synteny")
    parser.add_argument("--min", type=int, default=1, required=False, help="minimum required sequences aligned to be in the resulting graph")
    parser.add_argument("--gfa", type=str, help="Output the Pan-Synteny graph in GFA v1.1 format", required=False, default=None)
    parser.add_argument("--embed_info", help="MD5 and coordinates will be embedded at each node", required=False, default=False , action='store_true')
    parser.add_argument("--gexf_directed", action='store_true', help="Keep the GEXF graph directed (by default it is converted to an undirected graph for UI/Layout compatibility)")
    parser.add_argument("feature_files", type=str, nargs="*", default=["-"], help="Files of varying format specifing group, genome, contig, feature, and start in sorted order. stdin also accepted")


    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()
    pargs = parser.parse_args()
    if type(pargs.output) == str:
        logging.basicConfig(filename=pargs.output+'.log.txt',level=logging.INFO, filemode='w', format='%(message)s')
    root = logging.getLogger()
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.WARNING)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    handler.setFormatter(formatter)
    root.addHandler(handler)
    global embed_info
    embed_info=pargs.embed_info

    gmaker=GraphMaker(feature_files=pargs.feature_files, feature_selector=pargs.feature_selector, file_type=pargs.file_type, context=pargs.context, ksize=pargs.ksize, break_conflict=False, label_function= (not pargs.no_function),diversity=pargs.diversity, minSeq=pargs.min, alpha=pargs.alpha)
    if pargs.order_contigs !="none":
        if pargs.contig_output == None:
            sys.stderr.write("Need contig_ouptut parameter specified to output contig ordering\n")
            sys.exit()
        gmaker.order_contigs=True
        if pargs.order_contigs == "area":
            gmaker.traverse_priority = pargs.order_contigs
    gmaker.processFeatures()
    gmaker.RF_to_PG()
    if pargs.rfgraph != None:
        nx.readwrite.write_gexf(gmaker.rf_graph, pargs.rfgraph)
    gmaker.checkPGGraph()
    gmaker.checkRFGraph()
    gmaker.calcStatistics()
    gmaker.compute_graph_metrics()
    gmaker.detect_and_annotate_superbubbles()

    inversions_count, sv_edges_count = gmaker.tag_structural_variants()
    alt_paths_count, scaffolds_count = gmaker.flag_bridges()

    block_ids = gmaker.annotate_major_blocks_tfs(min_node_fraction=0.05) # Any component < 5% of nodes is a "Fragment"

    if pargs.gfa:
        gmaker.export_gfa(pargs.gfa)

    # Create the summary JSON blob and write TSV reports BEFORE serialize_graph_for_gexf mutates sets/dicts
    contig_map = {gen_id: list(contigs.keys()) for gen_id, contigs in gmaker.replicon_map.items()}
    total_contigs = sum(len(contigs) for contigs in contig_map.values())

    # Dynamically tally the biological ontology tags
    cnv_clusters = set()
    fragmented_cnv_clusters = set()
    superbubbles = set()
    synteny_breakpoints_count = 0
    assembly_repeat_breaks_count = 0
    
    for n, d in gmaker.pg_graph.nodes(data=True):
        # Unique Paralogs
        if d.get("cnv_cluster_id", 0) > 0:
            cnv_clusters.add(d["cnv_cluster_id"])
            if d.get("repeat_ambiguity_fragmentation"):
                fragmented_cnv_clusters.add(d["cnv_cluster_id"])
        
        # Unique Superbubbles
        sb_val = d.get("superbubble_id", set())
        if isinstance(sb_val, set):
            superbubbles.update(sb_val)
        elif sb_val:
            superbubbles.update(str(sb_val).split(","))
        
        # Unique Breakpoint Nodes
        if d.get("synteny_breakpoint"):
            synteny_breakpoints_count += 1
        if d.get("assembly_repeat_break"):
            assembly_repeat_breaks_count += 1

    # Compile the strict event counts
    summary_dict = {
        "total_genomes": len(gmaker.replicon_map),
        "total_contigs": total_contigs,
        "total_features": len(gmaker.feature_index),
        "total_nodes": gmaker.pg_graph.number_of_nodes(),
        "superbubbles": len(superbubbles),
        "copy_number_variants": len(cnv_clusters),
        "alternative_paths": alt_paths_count,                                    # 1 insertion = 1 event
        "structural_rearrangements": synteny_breakpoints_count + sv_edges_count, # Breakpoint nodes + Path Dropouts
        "inverted_blocks": inversions_count,                                     # 1 block = 1 event
        "assembly_breaks": assembly_repeat_breaks_count + scaffolds_count,       # Shattered repeat nodes + Scaffold bridges
        "fragmented_cnv": len(fragmented_cnv_clusters),                          # Shattered / fragmented paralog clusters
        "parameters": {"k": pargs.ksize, "min": pargs.min},
        "block_manifest": list(block_ids),
        "contig_map": contig_map
    }

    # Trigger the TSV exports alongside the main GEXF/GFA while pg_graph attributes are still native sets/dicts
    if isinstance(pargs.output, str):
        gmaker.write_reports(pargs.output)

    graph_to_write = gmaker.pg_graph
    if not pargs.gexf_directed:
        logging.warning("Converting to undirected graph for GEXF UI compatibility")
        graph_to_write = gmaker.get_undirected_pg_graph()
    gmaker.serialize_graph_for_gexf(graph_to_write)
    
    if pargs.layout:
        gexf_capture = io.BytesIO() 
        nx.readwrite.write_gexf(graph_to_write, gexf_capture) 
        # Decode the bytes into a string
        raw_gexf_str = gexf_capture.getvalue().decode('utf-8')
    
        # 1. Check if the user provided it via an environment variable (Great for HPC / module systems)
        layout_jar = os.environ.get('PANACONDA_LAYOUT_JAR')
        
        # 2. Check if it sits right next to the script in a flat /bin directory
        if not layout_jar:
            flat_bin_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "gexf_layout.jar")
            if os.path.exists(flat_bin_path):
                layout_jar = flat_bin_path
        layout_cmd=""        
        # 3. Fallback to the legacy submodule path
        if not layout_jar:
            layout_jar = os.path.join(os.path.dirname(os.path.realpath(__file__)), "layout", "pangenome_layout", "bin", "gexf_layout.jar")

        if not os.path.exists(layout_jar):
            logging.error(f"Layout failed: Cannot find gexf_layout.jar at {layout_jar}.")
            logging.error("Please set the PANACONDA_LAYOUT_JAR environment variable to the correct path.")
            # Gracefully save the un-laid-out graph instead of crashing
            raw_gexf_output = raw_gexf_str 
        else:
            logging.warning(f"Using layout jar at: {layout_jar}")
            layout_cmd = ["java", "-jar", layout_jar]
            logging.warning("laying out graph")
            logging.warning(" ".join(layout_cmd))
            
            # In Python 3, Popen with string input requires text=True
            p = Popen(layout_cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE, text=True)
            raw_gexf_output, err = p.communicate(input=raw_gexf_str)
            if err and err.strip():
                logging.warning("--- Layout Algorithm STDERR ---")
                logging.warning(err.strip())
                logging.warning("-------------------------------")
                
            if p.returncode != 0:
                logging.error(f"Layout algorithm failed with return code {p.returncode}")
                # Optional: Fallback to the un-laid-out graph if Java completely crashed
                if not raw_gexf_output:
                    raw_gexf_output = raw_gexf_str
                
    else:
        # If not using layout, intercept NetworkX bytes buffer to ensure clean JSON
        gexf_capture = io.BytesIO()
        nx.readwrite.write_gexf(graph_to_write, gexf_capture)
        
        raw_gexf_output = gexf_capture.getvalue().decode('utf-8')

    summary_json = json.dumps(summary_dict)
    summary_xml = f"    <summary>{summary_json}</summary>\n  "
    
    # Inject directly into the <meta> block (no regex cleaning needed on summary_xml)
    augmented_gexf = re.sub(r'</meta>', f'{summary_xml}</meta>', raw_gexf_output)  

    if isinstance(pargs.output, str):
        with open(pargs.output, 'w') as f:
            f.write(augmented_gexf)
    else:
        pargs.output.write(augmented_gexf)

    if pargs.order_contigs != "none":
        unsorted_file = pargs.contig_output+".unsorted"
        gmaker.write_contigs(pargs.contig_output, unsorted_file)

    
if __name__ == "__main__":
    start_time = time.time()
    main()
    sys.stderr.write(("--- %s seconds ---" % (time.time() - start_time))+"\n")
