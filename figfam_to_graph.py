#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
from networkx import Graph
from networkx import readwrite
import DOMLight, json
from collections import deque
##This script expects files which contain a list of consecutive protein family ID's for the replicons in an organism
## and some kind of summary information about where the kmers come from

##NCBI_TAX_ID     RANK    TAX_PATH
##590     genus   220341,90370,59201,28901,590

##NAME    NCBI_TAX_ID     ACCESSION       START_MIN
##FIG01045527     946034  AERV01000001    507

##Class for storing information about the origin of a Kmer
class figFamInfo():
	def __init__(self, acc, oid, pos):
		self.seq_accession=acc
		self.org_id=oid
		self.position=pos

##Class for storing all the figFamInfo in a particular node
class kmerNode():
	def __init__(self,nid):
		self.infoList={}
		self.nodeID=nid
		self.weightLabel=None
		self.weight=None
	#each cell in list stores info[x]=figFamInfo()
	def addInfo(self, cur_key, info):
		try: self.infoList[cur_key].append(info)
		except: self.infoList[cur_key]=[info]
	
##CALCULATE DIVERSITY QUOTIENT!!! GENUS/TOTAL GENOMES
##CALCULATE NORMALIZED NUMBER WEIGHT of NUMBER OF genomes in edge/ total number of genomes

##This class is for parsing figfam and summary file (table of taxonomy information) and storing it in a dictionary structure
##Parameters are filepaths and the size of kmer to use
class figFamStorage():
	def __init__(self, figfam_file, summary_file, ksize):
		print figfam_file
		print summary_file
		print str(ksize)
		self.kmerLookup={}#stores array for contig info and set for pointing to the next kmer
		self.figfamHash=[]#stores sets of coordinates for each figfam used to distinguish between paralogs/orthologs/distant orthologs
		self.summaryLookup={}
		self.summary_level=None#taxon level at which to summarize
		self.ksize=ksize #size of the kmer to store
		self.parseFigFam(figfam_file)
		self.parseSummary(summary_file)
		

	##This function checks whether the kmer is in the graph
	#and links kmer graph data structure appropriately
	#store kmers according to the combined protein family ids, and a set of IDs for which kmer comes next
	def addKmer(self, prev, kmer_key, fig_list):
		if not kmer_key in self.kmerLookup: 
			#pair: array for storing the contig info and set for storing the next kmer
			self.kmerLookup[kmer_key]=[kmerNode(kmer_key),set()]
		for fig_info in fig_list:
			self.kmerLookup[kmer_key][0].addInfo(fig_info[2], figFamInfo(fig_info[2], fig_info[1], fig_info[3]))#add information about kmer location
		if(prev!=None):
			#self.kkmerLookup[prev][1].add(kmer_key.split(',')[-1])#add the last figfam ID to the previous kmer so that it can find this kmer
			self.kmerLookup[prev][1].add(kmer_key)#add the last figfam ID to the previous kmer so that it can find this kmer
	##Create an ID for kmer
	##In case directionality is flipped for entire genome I am flipping each kmer
	##This shouldn't adversely affect inversions nor the overall result
	def makeID(self, k_info_list):
		k_list=[]
		for k in k_info_list:
			k_list.append(k[0])
		id_sep="|"
		result=None
		if(k_list[0]> k_list[-1]):
			k_list.reverse()
		result=id_sep.join(k_list)
		return result
	##get the id_set used for degree of similarity for protein family
	def getHashSet(self, fID, kmer):
		result=set()
		cur_info_list=[]
		try:
			cur_info_list=self.kmerLookup[kmer].infoList[fID]
		except:
			print "could not find "+fID+" in "+kmer
			pass
		for c in cur_info_list:
			result.add(c.position)
		return result
				
	##check to see if the figfam should be considered the same or not
	def checkIdentity(self, fID, id_set, threshold):
		matching_group=-1
		if fID in self.figfamHash:
			for idx, uid_set in enumerate(self.figfamHash[fID]):
				bigset= id_set if len(id_set) > len(uid_set) else uid_set
				score=len(id_set.intersection(uid_set))/float(len(bigset))
				if(score>=threshold):
					matching_group=idx
					#store the bigest id_set as the identifying one for this figfam
					if bigset != uid_set: self.figfamHash[fID][idx]=bigset.copy()
					break
			if matching_group == -1:
				self.figfamHash[fID].append(id_set)
				matching_group = len(self.figFamHash[fID])-1
		else:
			self.figfamHash[fID]=[id_set]
			matching_group=0
		return matching_group
					
	##Separate the kmer back into its parts
	def getParts(self, kmer):
		return(kmer.split('|'))
			
	def parseFigFam(self, figfam_file):
		
		num_fam=0
		inHandle=open(figfam_file, 'r')
		header=inHandle.readline()
		kmer_q=deque()
		prev_seq=""
		cur_seq=""
		prev_kmer=None
		#loop through figfams to create kmers
		for line in inHandle:
			num_fam+=1
			fig_info=line.strip().split("\t")
			cur_seq=fig_info[2]
			if(prev_seq != cur_seq and num_fam>1):
				kmer_q=deque()#clear kmer stack because switching replicons
				prev_kmer=None

			kmer_q.append(fig_info)#append the figfam ID
			if(len(kmer_q)>self.ksize):
				kmer_q.popleft()
				kmer=self.makeID(list(kmer_q))#put IDs together to make kmer ID
				self.addKmer(prev_kmer, kmer, list(kmer_q))
			elif(len(kmer_q)== self.ksize):
				kmer=self.makeID(list(kmer_q))#put IDs together to make kmer ID
				self.addKmer(prev_kmer, kmer, fig_info)#right now only passing in the last figfams information
			else:#kmer size is less than ksize
				kmer=None
			prev_seq=cur_seq # record which replicon we are on
			prev_kmer=kmer
		inHandle.close()
		
	#for a given kmer return a set of the organisms involved
	def getOrgSummary(self, kmer):
		result=set()
		if kmer in self.kmerLookup:
			for i in self.kmerLookup[kmer][0].infoList:
				result.add(i.org_id)
		return result
	
	def getTaxSummary(self,kmer):
		result=set()
		if kmer in self.kmerLookup:
			for i in self.kmerLookup[kmer][0].infoList:
				if(i.org_id in self.summaryLookup):
					result.add(self.summaryLookup[i.org_id])
		return result
	#for a given node return a set of the organisms involved
	def nodeOrgSummary(self,cnode):
		result=set()
		for i in cnode.infoList:
			result.add(i.org_id)
		return result
	
	def nodeTaxSummary(self,cnode):
		result=set()
		for i in cnode.infoList:
			if(i.org_id in self.summaryLookup):
				result.add(self.summaryLookup[i.org_id])
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
		header=inHandle.readline()
		for line in inHandle:
			summary_info=line.strip().split("\t")
			if(self.summary_level==None):
				self.summary_level=summary_info[1]
			summary_id=summary_info[0]
			ref_id=sumamry=summary_info[2].split(',')[0]
			self.summaryLookup[ref_id]=summary_id
		inHandle.close()

# undirected weighted
class pFamGraph(Graph):
	def __init__(self, storage, minOrg=2):
		#Graph.__init__(self, weighted=True)
		Graph.__init__(self)
		self.createGraph(storage, minOrg)
	def add_path_cumul_attr(self,nlist,**kwargs):
		edges=list(zip(nlist[:-1],nlist[1:]))#create list of edges
		for e in edges:
			if self.has_edge(*e):
				for k in kwargs:
					try: self.adj[e[0]][e[1]][k]=kwargs[k] | self.adj[e[0]][e[1]][k]# union of attribute
					except: 
						try: self.adj[e[0]][e[1]][k]=kwargs[k].copy()
						except: self.adj[e[0]][e[1]]=kwargs[k]
			else: self.add_edge(e[0],e[1],kwargs)
	def update_edge_weight(self, e_attr, divisor=1):
		for e in self.edges():
			#try: self.adj[e[0]][e[1]][e_attr]=list(self.adj[e[0]][e[1]][e_attr])
			#except: pass
			try: self.adj[e[0]][e[1]]['weight']=len(self.adj[e[0]][e[1]][e_attr])/float(divisor)
			except:
				try:self.adj[e[0]][e[1]]['weight']=0
				except: pass
			try: self.adj[e[0]][e[1]][e_attr]=",".join(self.adj[e[0]][e[1]][e_attr])
			except: pass
						
	##this function takes the storage class and constructs the graph from it
	def createGraph(self, storage, minOrg):
		num_orgs=len(storage.summaryLookup.keys())
		temp_size=len(storage.kmerLookup.keys())
		total_tax=len(storage.completeTaxSummary())
		print " ".join(["starting",str(temp_size),str(total_tax),str(num_orgs)])
		for kmer in storage.kmerLookup.keys():
			cur_node=storage.kmerLookup[kmer][0]
			org_part1=storage.nodeOrgSummary(cur_node)#a set of the organisms involved in this part of the graph
			tax_ids=storage.nodeTaxSummary(cur_node)#summarizes set of tax ids for this node
			cur_node.weightLabel="Percent genera"
			cur_node.weight=len(tax_ids)/float(total_tax)
			if(len(org_part1)>=minOrg):
				nodeSet1=storage.getParts(cur_node.nodeID)
				self.add_path_cumul_attr(nodeSet1, orgs=org_part1)
				for n in nodeSet1:
					self.node[n]['weight']=cur_node.weight
		self.update_edge_weight('orgs',divisor=float(total_tax))
			#edge_added=False
			#add all the edges to the graph
			#for next_kmer in storage.kmerLookup[kmer][1]:
			#	next_node=storage.kmerLookup[next_kmer][0]
			#	nxt_orgs=storage.nodeOrgSummary(next_node)
			#	orgs_in_edge=len(nxt_orgs.intersection(org_part1))
			#	if(orgs_in_edge>=2):
			#		nodeSet1=storage.getParts(cur_node.nodeID)
			#		nodeSet2=storage.getParts(next_node.nodeID)
			#		edge_added=True
			#		cur_weight=orgs_in_edge/float(num_orgs)
			#		self.add_edge(cur_node.nodeID, next_node.nodeID, weight=cur_weight)
					#need to fix this so that KmerNode is stored instead of string
					#work out write_graphml
			#		next_node.weight=len(storage.nodeTaxSummary(next_node))/float(total_tax)
			#		self.node[next_node.nodeID]['weight']=next_node.weight
			#if(edge_added):
			#	self.node[cur_node.nodeID]['weight']= cur_node.weight



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


def main(init_args):
	if(len(init_args)<4):
		sys.stderr.write("Usage: figfam_to_graph.py figfam_table summary_table output_folder k-size\n")
		sys.exit()
	k_size=int(init_args[3])
	fstorage=figFamStorage(init_args[0], init_args[1], k_size)
	out_basename=os.path.splitext(os.path.basename(init_args[0]))[0] #get basename of the file to name output
	out_folder=os.path.expanduser(init_args[2])
	out_file=os.path.join(out_folder,out_basename)
	pgraph=pFamGraph(fstorage)
	csize=pgraph.order()
	toGML(pgraph, out_file+".graphml")
	result_handle=open(out_file+".xgmml", 'w')
	pgraph.toXGMML(result_handle)
	result_handle.close()
	result_handle=open(out_file+".json", 'w')
	pgraph.toJSON(result_handle)
	result_handle.close()
	
if __name__ == "__main__":
	main(sys.argv[1:])
