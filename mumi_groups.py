#!/usr/bin/env python
# -*- coding: utf-8 -*-

#This script is for reading in a table specifying groups of genomes and then running mummer and MUMi on all pairs within a group
#Assumes the installation of mummer and /give_mumi2.pl script

import os, sys, textwrap, getopt
import subprocess

def GetMUMi(input_dir, g1, g2):
	#mummer -mum -b -c -l 19 family_rep_seq/NC_004741.fna family_rep_seq/NC_009085.fna >mumi_test/mum_out.txt
	mummer_com="mummer -mum -b -c -l 19 "
	#"./give_mumi2.pl mumi_test/mum_out.txt -l1 4599354 -l2 3976747"
	mumi_com="./give_mumi2.pl mum_out.tmp -l1 "
	fna1=os.path.join(input_dir,g1+".fna")#set genbank file name
	fna2=os.path.join(input_dir,g2+".fna")#set genbank file name
	mumi_dist="Error NoDistance"
	if(os.path.isfile(fna1) and os.path.isfile(fna2)):
		mummer_com+=" ".join([fna1, fna2, ">mum_out.tmp"])
		#Run mummer
		sys.stderr.write(mummer_com+"\n")
		returncode=subprocess.call(mummer_com, shell=(sys.platform!="win32"))	
		sys.stderr.write(str(returncode)+"\n")
		#Find the length of the sequences
		len_com="grep -v \">\" "+fna1+" | tr -d '\\012' | wc -m"
		child = subprocess.Popen(len_com, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
		seq_len1, error=child.communicate()
		sys.stderr.write(error+"\n")
		len_com="grep -v \">\" "+fna2+" | tr -d '\\012' | wc -m"
		child = subprocess.Popen(len_com, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
		seq_len2, error=child.communicate()
		sys.stderr.write(error+"\n")
		mumi_com+=seq_len1+" -l2 "+seq_len2
		child = subprocess.Popen(mumi_com, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
		mumi_dist, error=child.communicate()
	else:	
		sys.stderr.write("input sequences not found: "+fna1+ " or "+fna2+"\n")
	return float(mumi_dist.split()[-1])

def AvgMUMi(input_dir, genomes):
	sum_score=0
	i=0
	count=0
	while i< len(genomes):
		j=i+1
		while j< len(genomes):
			cur_score=GetMUMi(input_dir, genomes[i], genomes[j])
			sum_score+=cur_score
			count+=1
			print "\t".join([genomes[i], genomes[j], str(cur_score)])
			j+=1
		i+=1
	return float(sum_score)/float(count)

def UsageMess():
	print "Usage: mg_mumi_groups.py [-i fasta seq dir] [-t groups table]"
	sys.exit()

def main(init_args):

	if(len(init_args)<4):
		UsageMess()

	optlist, args=getopt.getopt(init_args, 'i:t:')
	for o,v in optlist:
		if(o=='-i'):
			input_dir=v
		elif(o=="-t"):
			table_file=v
	table_handle=open(table_file, 'r')
	for line in table_handle:
		cols=line.strip().split()
		genomes=cols[1:]
		avg_i=AvgMUMi(input_dir, genomes)
		print cols[0]+"\t"+str(avg_i)



if __name__ == "__main__":
    main(sys.argv[1:])
