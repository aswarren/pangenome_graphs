#!/usr/bin/env python
# -*- coding: utf-8 -*-

#This script is for reading in a table specifying groups of genomes and then running ete2 to get distances from common ancestor
#Assumes the installation of ete2 package

import os, sys, json, getopt
from ete2 import Tree


def get_distances(input_dir, group, genomes):
	results = {}
	in_file=os.path.join(input_dir,group+".nwk")
	try:
		t = Tree(in_file)
		a=t.get_common_ancestor(*genomes)
	except Exception, e:
		sys.stderr.write("Problem with newick "+in_file+"\n")
		print "Unexpected error:", str(e)
		sys.exit()
	for leaf in genomes:
		results[leaf]=t.get_distance(a,leaf)
	return results

def usage_msg():
	print "Usage: tree_distances.py [-i newick dir] [-t groups table]"
	sys.exit()

def main(init_args):

	if(len(init_args)<4):
		usage_msg()

	optlist, args=getopt.getopt(init_args, 'i:t:')
	for o,v in optlist:
		if(o=='-i'):
			input_dir=v
		elif(o=="-t"):
			table_file=v
	table_handle=open(table_file, 'r')
	results={}
	for line in table_handle:
		cols=line.strip().split()
		group=cols[0]
		genomes=cols[1:]
		cur_dists=get_distances(input_dir, group, genomes)
		results[group]=cur_dists
		avg_dist = sum(cur_dists.values())/float(len(cur_dists))
		print "\t".join([group,str(max(cur_dists.values())),str(avg_dist)])
	table_handle.close()
	#print json.dumps(results)

if __name__ == "__main__":
    main(sys.argv[1:])
