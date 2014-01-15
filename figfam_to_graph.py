#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
from networkx import Graph
from networkx import readwrite
import DOMLight, json
from collections import deque
from collections import OrderedDict
##This script expects files which contain a list of consecutive protein family ID's for the replicons in an organism
## and some kind of summary information about where the kmers come from

##NCBI_TAX_ID     RANK    TAX_PATH
##590     genus   220341,90370,59201,28901,590

##NAME    NCBI_TAX_ID     ACCESSION       START_MIN
##FIG01045527     946034  AERV01000001    507


def warning(*objs):
	for o in objs:
		print >> sys.stderr, o

##Class for storing information about the origin of a Kmer
class geneInfo():
	def __init__(self, acc, oid, pos, func):
		self.replicon_id=acc
		self.org_id=oid
		self.position=pos
		self.function=func
	def getLocation(self):
		return (self.replicon_id, self.position)
	def getReplicon(self):
		return self.replicon_id

##Class for storing all the geneInfo in a particular node
class kmerNode():
	def __init__(self,nid):
		self.infoList={}#stores information about location of the figfams that make up kmer as geneInfo objects infoList[figfam ID] = geneInfo()
		self.nodeID=nid
		self.weightLabel=None
		self.weight=None
	#each cell in list stores info[x]=geneInfo()
	def addInfo(self, cur_key, info):
		try: self.infoList[cur_key].append(info)
		except: self.infoList[cur_key]=[info]
	def getReplicons(self):
		result=set([])
		#all the replicons should be the same for each fam in this kmer
		fam=self.infoList.values()[0]
		for info in fam: 
			result.add(info.getReplicon())
		return result
	def testNode(self):
		#make sure that all the families in the kmer come from same replicons
		ref_set=set([x.getReplicon() for x in self.infoList.values()[0]])
		for fam in self.infoList.values():
			test_set=set([])
			for info in fam:
				test_set.add(info.getReplicon())
			if test_set != ref_set:
				warning("kmer "+self.nodeID+" has inconsistent replicons")
				



#provides a summary of where this family occurs
#a family may be differentiated into multiple version depending on its ocurrence in kmers
class famVersion():
	def __init__(self, id_set):
		self.instances=id_set.copy() #set of locations that identify this version of family
		self.organisms=set()
		self.tax_summary=set()
		self.replicons=set()
		self.locations=set()
		self.functions=set()
	#returns of summary items
	def get_summary(self):
		for i in self.instances:
			self.replicons.add(i.replicon_id)
			self.organisms.add(i.org_id)
			self.locations.add(':'.join(i.getLocation()))
			self.functions.add(i.function)
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
						self.versions[idx].instances=self.versions[idx].instances.union(v.instances)
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
				v.instances=v.instances.union(id_set)
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


		
		
		
##CALCULATE DIVERSITY QUOTIENT!!! GENUS/TOTAL GENOMES
##CALCULATE NORMALIZED NUMBER WEIGHT of NUMBER OF genomes in edge/ total number of genomes

##This class is for parsing figfam and summary file (table of taxonomy information) and storing it in a dictionary structure
##Parameters are filepaths and the size of kmer to use
class figFamStorage():
	def __init__(self, figfam_file, summary_file, ksize, ignore_fams=set([])):
		print figfam_file
		print summary_file
		print str(ksize)
		self.replicon_ids=set()
		self.ignore_fams=ignore_fams
		self.kmerLookup=OrderedDict()#Ordered so that results are reproducible. stores array for contig info and set for pointing to the next kmer
		self.figfamHash={}#stores sets of coordinates for each figfam used to distinguish between paralogs/orthologs/distant orthologs
		self.summaryLookup={}
		self.locationHash={}#stores the disambiguated 'version' of the protein family. hashed by (seq. accession, location)
		self.geneHash={} #storing information about the individual genes
		self.summary_level=None#taxon level at which to summarize
		self.ksize=ksize #size of the kmer to store
		self.parseFam(figfam_file)
		self.parseSummary(summary_file)
		
		

	##This function checks whether the kmer is in the graph
	#and links kmer graph data structure appropriately
	#store kmers according to the combined protein family ids, and a set of IDs for which kmer comes next
	def addKmer(self, prev, kmer_key, fig_list):
		if not kmer_key in self.kmerLookup: 
			#pair: array for storing the contig info and set for storing the next kmer
			self.kmerLookup[kmer_key]=[kmerNode(kmer_key),set()]
		for fig_info in fig_list:
			#print fig_info
			fID, replicon_id, organism_id, position, fam_function= fig_info[0], fig_info[3], fig_info[1], fig_info[4], fig_info[5]
			gene_lookup=(replicon_id, position)
			#key the gene lookup by replicon_id and position
			if not gene_lookup in self.geneHash:
				target=geneInfo(replicon_id, organism_id, position, fam_function)
				self.geneHash[gene_lookup]=target
			else:
				target=self.geneHash[gene_lookup]
			self.kmerLookup[kmer_key][0].addInfo(fID, target)#add information about kmer location
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
	#Parameters fID-figfamID; kmer-the kmer ID
	#Return a hash value (set) of all the positions of a figfam
	def getHashSet(self, fID, kmer):
		result=set()
		cur_info_list=[]
		try:
			cur_info_list=self.kmerLookup[kmer][0].infoList[fID]
		except:
			print "could not find "+fID+" in "+self.kmerLookup[kmer][0].infoList.keys()
			pass
		for c in cur_info_list:
			result.add(c)
		return result
		
				
	##check to see if the figfam should be considered the same or not
	#Parameters fID-figfam ID, id_set-
	#The graph is composed of figfams that represent many positions in many genomes
	#when expanding out from kmers to figfams the figfams need to be checked to see if they represent
	#equivalent positions
	#there is potentially a flaw if  the same figfam occurs at the same base in two different genomes
	def checkIdentity(self, fID, kmer, threshold):
		id_set = self.getHashSet(fID,kmer)				

		#check if the figfam is in storage
		if fID in self.figfamHash:
			#iterate through all the sets of coordinates for the figfam
			#to find if this figfam-node has already been encountered
			#if it has make sure to update the positions with the union
			#if it hasn't add these positions to the list.
			#each node (consists of a figfam that can occur in multiple organisms)
			self.figfamHash[fID].add_instance(id_set, threshold, self.locationHash)	
		else:	
			self.figfamHash[fID]=famInfo(fID)
			self.figfamHash[fID].add_instance(id_set, threshold, self.locationHash)
		#MODIFY Some kindof kmer data structure to keep track of which figfams are which version of themselves.
		return 0
					
	##Separate the kmer back into its parts
	def getParts(self, kmer):
		return(kmer.split('|'))
			
	def parseFam(self, figfam_file, ignore=True):
		
		num_fam=0
		inHandle=open(figfam_file, 'r')
		#header=inHandle.readline()
		kmer_q=deque()
		prev_seq=""
		cur_seq=""
		prev_kmer=None
		#loop through figfams to create kmers
		for line in inHandle:
			if line.startswith('#'):
				continue

			fig_info=line.strip().split("\t")

			if ignore and fig_info[0] in self.ignore_fams:
				continue
			num_fam+=1
			cur_seq=fig_info[3]
			self.replicon_ids.add(cur_seq)
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
				self.addKmer(prev_kmer, kmer, list(kmer_q))#right now only passing in the last figfams information
			else:#kmer size is less than ksize
				kmer=None
			prev_seq=cur_seq # record which replicon we are on
			prev_kmer=kmer
		inHandle.close()
		
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
					result.add(self.summaryLookup[i.org_id])
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
		for i in cnode.infoList.values()[0]:
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
		#header=inHandle.readline()
		for line in inHandle:
			if line.startswith('#'):
				continue
			summary_info=line.strip().split("\t")
			if(self.summary_level==None):
				self.summary_level=6
			ref_id=summary_info[1]
			summary_id=summary_info[-1].split(',')[self.summary_level]
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

	#update the edge weight based on a designated attribute
	#also flatten to a string since writing list objects isn't supported
	def update_edges(self, weight_attr, divisor=1, label_attr='replicons', remove_attrs=[]):
		for e in self.edges():
			#try: self.adj[e[0]][e[1]][e_attr]=list(self.adj[e[0]][e[1]][e_attr])
			#except: pass
			self.adj[e[0]][e[1]]['label']=''
			try: self.adj[e[0]][e[1]]['weight']=len(self.adj[e[0]][e[1]][weight_attr])/float(divisor)
			except:
				try:self.adj[e[0]][e[1]]['weight']=0
				except: pass
			try: self.adj[e[0]][e[1]][label_attr]=", ".join(self.adj[e[0]][e[1]][label_attr])
			except: pass
			for r in remove_attrs:
				try: self.adj[e[0]][e[1]].pop(r,None)
				except: pass
				
	def update_node_cumul_attr(self, n_id, **kwargs ):
		if n_id in self.node:
			for k in kwargs:
				try: self.node[n_id][k]=kwargs[k] | self.node[n_id][k]
				except:
					try: self.node[n_id][k]=kwargs[k].copy()
					except: print "cannot add attribute to node "+str(n_id)
	
	#calculate the node weight and change the set attributes to string
	#so that they can be written by graphml writer
	def update_node_attr_final(self, weight_attr, divisor=1):
		for n in self.nodes():
			try: self.node[n]['weight']=len(self.node[n][weight_attr])/float(divisor)
			except: pass
			for a in self.node[n]:
				if type(self.node[n][a])==set:
					self.node[n][a] = ','.join(self.node[n][a])

	## create a graph node from a kmer
	#Expands Each kmer to FIGFAM nodes ensuring that the kmer represents a certain number of unique positions/organisms before expanding
	#parameters: storage- the figfam storage used to look things up
	#kmer - the kmer being processed; total_tax- the total number of taxonomy groups represented in the graph; minOrg- the number of positions that need to be represented to qualify as a part of the summary graph
	def kmer_to_node(self, storage, kmer, total_tax, minOrg):
		cur_node=storage.kmerLookup[kmer][0]
		org_part1=storage.nodeOrgSummary(cur_node)#a set of the organisms involved in this part of the graph
		#print org_part1

		if(len(org_part1)>=minOrg):
			#get a list of individual figfams that make up the kmer
			pfamily_ids=storage.getParts(cur_node.nodeID)
			#for each figfam get list of locations and call checkIdentity
			for famID in pfamily_ids: storage.checkIdentity(famID, cur_node.nodeID,threshold=1)
				
	def kmer_to_node2(self, storage, kmer, total_tax, minOrg):
		cur_node=storage.kmerLookup[kmer][0]
		org_part1=storage.nodeOrgSummary(cur_node)#a set of the organisms involved in this part of the graph
		if(len(org_part1)>=minOrg):
			tax_ids=storage.nodeTaxSummary(cur_node)#summarizes set of tax ids for this node
			pfamily_ids=storage.getParts(cur_node.nodeID)
			nodeList2=[]
			summaryList=[]
			for memID in pfamily_ids:
				target_set=storage.getHashSet(memID,cur_node.nodeID)
				fam_version=storage.locationHash[list(target_set)[0]]
				nodeList2.append(fam_version[0]+':'+fam_version[1])
				summary_info=storage.figfamHash[fam_version[0]].versions[int(fam_version[1])].get_summary()
				summary_info['tax_summary']=tax_ids
				summaryList.append(summary_info)
				#if fam_version[0] =="FIG00229272" and fam_version[1]=="0":
				#	pass
			#for this kmer get all the replicons it is involved in
			replicons=cur_node.getReplicons()
			#cur_node.testNode()
			self.add_path_cumul_attr(nodeList2, orgs=org_part1, replicons=replicons) #this also creats nodes in the graph
			for idx, n in enumerate(nodeList2):#now that the nodesare created you can update them with attributes
				attr_list={}
				#HACK exclude replicons from list since locations already have them
				for k in summaryList[idx].keys():
					if k == 'functions':
						attr_list['label']=summaryList[idx][k] 
					if k != 'replicons':
						attr_list[k]=summaryList[idx][k]
				self.update_node_cumul_attr(n, **attr_list)
				
				
				
	##this function takes the storage class and constructs the graph from it
	def createGraph(self, storage, minOrg):
		num_orgs=len(storage.summaryLookup.keys())
		temp_size=len(storage.kmerLookup.keys())
		total_tax=len(storage.completeTaxSummary())
		print " ".join(["starting",str(temp_size),str(total_tax),str(num_orgs)])
		for kmer in storage.kmerLookup.keys():
			self.kmer_to_node(storage,kmer,total_tax,minOrg)
			#for some reason this second loop is necessary to explode certain kmers. see if fID =="FIG00229272": Brucella
			#for next_kmer in storage.kmerLookup[kmer][1]:
			#	self.kmer_to_node(storage, next_kmer,total_tax,minOrg)
		for kmer in storage.kmerLookup.keys():
			self.kmer_to_node2(storage,kmer,total_tax,minOrg)
		self.update_edges(weight_attr='orgs',divisor=float(num_orgs), label_attr="replicons", remove_attrs=['orgs'])
		self.update_node_attr_final('tax_summary', divisor=float(total_tax))
		node_handle=open('new_loop_node_list.txt','w')
		for n in self.nodes_iter():
			node_handle.write(n+"\n")
		node_handle.close()
		edge_handle=open('new_loop_edge_list.txt','w')
		for e in self.edges_iter():
			edge_handle.write(str(e)+"\n")
		edge_handle.close()
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
	if len(init_args)>=5:
		ignore_fams=init_args[4].replace(' ','').split(',')
	fstorage=figFamStorage(init_args[0], init_args[1], k_size, ignore_fams=set(['FIG00638284','FIG01306568']))
	out_basename=os.path.splitext(os.path.basename(init_args[0]))[0] #get basename of the file to name output
	out_folder=os.path.expanduser(init_args[2])
	out_file=os.path.join(out_folder,out_basename)
	pgraph=pFamGraph(fstorage)
	csize=pgraph.order()
	toGML(pgraph, out_file+".graphml")
	readwrite.write_gexf(pgraph, out_file+".gexf")
	result_handle=open(out_file+".xgmml", 'w')
	pgraph.toXGMML(result_handle)
	result_handle.close()
	result_handle=open(out_file+".json", 'w')
	pgraph.toJSON(result_handle)
	result_handle.close()
	
if __name__ == "__main__":
	main(sys.argv[1:])
