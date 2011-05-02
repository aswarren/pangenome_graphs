#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
from networkx import Graph
from networkx import readwrite
import DOMLight
from collections import deque


##Class for storing information about the origin of a Kmer
class figFamInfo():
	def __init__(self, acc, oid, pos):
		self.seq_accession=acc
		self.org_id=oid
		self.position=pos

##Class for storing all the figFamInfo in a particular node
class kmerNode():
	def __init__(self,nid):
		self.infoList=[]
		self.nodeID=nid
		self.weightLabel=None
		self.weight=None
	def addInfo(self, info):
		self.infoList.append(info)
	
##CALCULATE DIVERSITY QUOTIENT!!! GENUS/TOTAL GENOMES
##CALCULATE NORMALIZED NUMBER WEIGHT of NUMBER OF genomes in edge/ total number of genomes

##This class is for parsing figfam and summary file (table of taxonomy information) and storing it in a dictionary structure
##Parameters are filepaths and the size of kmer to use
class figFamStorage():
	def __init__(self, figfam_file, summary_file, ksize):
		self.kmerLookup={}#stores array for contig info and set for pointing to the next kmer
		self.summaryLookup={}
		self.summary_level=None#taxon level at which to summarize
		self.ksize=ksize #size of the kmer to store
		self.parseFigFam(figfam_file)
		self.parseSummary(summary_file)
		

	##This function checks whether the kmer is in the graph
	#and links kmer graph data structure appropriately	
	def addKmer(self, prev, kmer_key, fig_info):
		if not kmer_key in self.kmerLookup: 
			#pair: array for storing the contig info and set for storing the next kmer
			self.kmerLookup[kmer_key]=[kmerNode(kmer_key),set()]
		self.kmerLookup[kmer_key][0].addInfo(figFamInfo(fig_info[2], fig_info[1], fig_info[3]))#add information about kmer location
		if(prev!=None):
			#self.kkmerLookup[prev][1].add(kmer_key.split(',')[-1])#add the last figfam ID to the previous kmer so that it can find this kmer
			self.kmerLookup[prev][1].add(kmer_key)#add the last figfam ID to the previous kmer so that it can find this kmer
	##Create an ID for kmer
	def makeID(self, k_list):
		id_sep="|"
		result=None
		if(k_list[0]> k_list[-1]):
			k_list.reverse()
		result=id_sep.join(k_list)
		return result
			
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

			kmer_q.append(fig_info[0])#append the figfam ID
							
			if(len(kmer_q)>self.ksize):
				kmer_q.popleft()
				kmer=self.makeID(list(kmer_q))#put IDs together to make kmer ID
				self.addKmer(prev_kmer, kmer, fig_info)
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
	def __init__(self, storage):
		#Graph.__init__(self, weighted=True)
		Graph.__init__(self)
		self.createGraph(storage)
		
						
	##this function takes the storage class and constructs the graph from it
	def createGraph(self, storage):
		num_orgs=len(storage.summaryLookup.keys())
		temp_size=len(storage.kmerLookup.keys())
		total_tax=len(storage.completeTaxSummary())
		for kmer in storage.kmerLookup.keys():
			cur_node=storage.kmerLookup[kmer][0]
			org_part1=storage.nodeOrgSummary(cur_node)#a set of the organisms involved in this part of the graph
			tax_ids=storage.nodeTaxSummary(cur_node)#summarizes set of tax ids for this node
			cur_node.weightLabel="Percent genera"
			cur_node.weight=len(tax_ids)/float(total_tax)
			#add all the edges to the graph
			for next_kmer in storage.kmerLookup[kmer][1]:
				next_node=storage.kmerLookup[next_kmer][0]
				cur_orgs=storage.nodeOrgSummary(next_node)
				orgs_in_edge=len(cur_orgs.intersection(org_part1))
				if(orgs_in_edge>2):
					cur_weight=orgs_in_edge/float(num_orgs)
					self.add_edge(cur_node.nodeID, next_node.nodeID, weight=cur_weight)
			if(len(storage.kmerLookup[kmer][1])>0):
				self.add_node(cur_node.nodeID, weight= cur_node.weight, ID=cur_node.nodeID)

	def toXGMML(self, fhandle):
		xml = DOMLight.XMLMaker()
		fhandle.write("""<?xml version="1.0" encoding="UTF-8"?>
		<graph xmlns="http://www.cs.rpi.edu/XGMML" directed="0" label="PFam assembly">
		""")
		cid = 0
		cur_ids = {}
		for cn in self.nodes_iter():
			cur_ids[cn] = cid
			fhandle.write(str(xml.cn({'id': cid, 'label': cn}, '<att type="real" name="weight" value="'+str(1)+'"/>')) + "\n")
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
	if(len(init_args)<3):
		sys.stderr.write("Usage: figfam_to_graph.py figfam_table summary_table output_folder\n")
		sys.exit()
	k_size=3
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
	
	
if __name__ == "__main__":
	main(sys.argv[1:])