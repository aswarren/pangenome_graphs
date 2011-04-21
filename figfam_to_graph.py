#!/usr/bin/env python
import os, sys
from networkx import Graph


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