#!/usr/bin/env python
import os, sys
from networkx import Graph

##Class for storing information about the origin of a Kmer
class figFamInfo():
	def __init__(self, acc, oid, pos):
		self.seq_accession=acc
		self.org_id=oid
		self.position=pos
		
##CALCULATE DIVERSITY QUOTIENT!!! GENUS/TOTAL GENOMES
##CALCULATE NORMALIZED NUMBER WEIGHT

##This class is for parsing figfam and summary file (table of taxonomy information) and storing it in a dictionary structure
##Parameters are filepaths and the size of kmer to use
class figFamStorage():
	def __init__(self, figfam_file, summary_file, ksize):
		kmerLookup={}
		summaryLookup={}
		self.parseFigFam(figfam_file)
		self.parseSummary(summary_file)
		self.summary_level=""#taxon level at which to summarize
		self.ksize=ksize #size of the kmer to store

	##This function checks whether the kmer is in the graph
	#and links kmer graph data structure appropriately	
	def addKmer(self, prev, kmer_key, fig_info):
		if not kmer_key in kmerLookup: 
			#pair: array for storing the contig info and set for storing the next kmer
			self.kmerLookup[kmer_key]=[[],set()]
		self.kmerLookup[kmer_key][0].append(figFamInfo(fig_info[2], fig_info[1], fig_info[3]))#add kmer to dictionary
		self.kkmerLookup[prev][1].add(kmer_key.split(',')[-1])#add the last figfam ID to the previous kmer so that it can find this kmer
	def parseFigFam(self, figfam_file):
		num_fam=0
		inHandle=open(figfam_file, 'r')
		header=inHandle.readline()
		kmer_stack=[]
		prev_seq=""
		cur_seq=""
		prev_kmer=None
		for line in inHandle:
			num_fam+=1
			fig_info=line.strip().split("\t")
			kmer_stack.append(fig_info[0])#append the figfam ID
			cur_seq=fig_info[2]
			if(prev_seq != cur_seq and num_fam>1):
				kmer_stack=[]#clear kmer stack
			if(len(kmer_stack)>self.ksize):
				kmer_stack.pop()
				kmer=",".join(kmer_stack)#put IDs together to make kmer ID
				addKmer(self, prev_kmer, kmer, fig_info)
			elif(len(kmer_stack)< self.ksize):
				continue
			else:#else exactly ksize
				kmer=",".join(kmer_stack)#put IDs together to make kmer ID
				addKmer(self, prev_kmer, kmer, fig_info)
			prev_seq=cur_seq # record which replicon we are on
			prev_kmer=kmer
		inHandle.close()


	#expects summary taxid, tax level, and the taxpath comma seperated values
	def parseSummary(self, summary_file):
		inHandle=open(summary_file, 'r')
		header=inHandle.readline()
		for line in inHandle:
			summary_info=line.split("\t")
			

# undirected weighted
class pfamGraph(DiGraph):
	def __init__(self, storage):
		DiGraph.__init__(self, weighted=True)
		self.createGraph(storage)
						
	##this function takes the storage class and constructs the graph from it
	def createGraph(self, storage):
		self.makeWeighted() #add weights to the graph
		
	def makeWeighted(self):
		wc = weightingClass(self.storage, species, aspect, directed=False)
		return wc.makeWeighted(self)
					

	def toXGMML(self, fhandle):
		xml = DOMLight.XMLMaker()
		fhandle.write("""<?xml version="1.0" encoding="UTF-8"?>
		<graph xmlns="http://www.cs.rpi.edu/XGMML" directed="1" label="GOGrapher Network">
		""")
		cid = 0
		goids = {}
		for node in self.nodes_iter():
			goids[node.goid] = cid
			fhandle.write(str(xml.node({'id': cid, 'label': node.goid})) + "\n")
			cid += 1
		count = 0
		for edge in self.edges_iter():
			fhandle.write(str(xml.edge({'source': goids[edge[0].goid], 'target': goids[edge[1].goid], 'weight': edge[2]['weight'], 'label': ""})) + "\n")
			if count == 1000:
				break
			count += 1
		fhandle.write("</graph>")
																													

	## Get weighted edgesD from this graph.
	def edges(self):
		# This is just the code from networkx.graph - except call our
		return list(self.edges_iter())
				

	## Overwrite Graph edges_iter method so we get weight information too
	def edges_iter(self, nbunch=None):
		for edge in DiGraph.edges_iter(self, nbunch, data=True):
			yield edge



def main(init_args):
	
	
if __name__ == "__main__":
	main(sys.argv[1:])