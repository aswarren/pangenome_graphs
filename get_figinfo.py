#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Wed Dec 18 22:04:38 2013

@author: anwarren
"""

import requests, json, sys
import time
import pandas as pd

ts = time.time()


#brucella = 234
#abortus = 29461
taxURL="http://macleod.vbi.vt.edu:8080/solr/genomesummary/select/?"
taxQuery="q=rast_cds:[1+TO+*]+AND+"
taxidQuery="+taxon_lineage_ids:"
gidQuery="gid:"
taxFormat="&indent=on&wt=json&fl=genome_name,genome_info_id,ncbi_tax_id,taxon_lineage_ids"
feature_url="http://macleod.vbi.vt.edu:8080/solr/dnafeature/select/?q="
feature_conditions="+AND+figfam_id:[*+TO+*]&sort=sequence_info_id+asc,start_max+asc&fl=figfam_id,gid,ncbi_tax_id,sequence_info_id,start_max,end_min,na_feature_id"
format_string="indent=on&wt=json"
fam_url='http://macleod.vbi.vt.edu:8080/solr/figfam-dic/select/?q='
##get SOLR query results
def get_solr_result(query_url):
    results=[]
    start=0
    rows=10000
    currentQuery=query_url+"&start="+str(start)+"&rows=1"
    response=requests.get(currentQuery)
    if(not response.ok):
       print "PATRIC solr API not responding"
       return results
    result_summary=json.loads(response.content)
    try:
        numFound=result_summary['response']['numFound']
    except (KeyError):
        print "Error in getting results returned from PATRIC"
        return results
    while start < numFound:
        if (numFound-start)< rows:
            limit=numFound-start
        else:
            limit=rows
        try:
            currentQuery=query_url+"&start="+str(start)+"&rows="+str(rows)
            response=requests.get(currentQuery)
            assert response.ok
            json_result=json.loads(response.content)
            cur_rows=json_result['response']['docs']
        except (AssertionError, KeyError,IndexError):
            print "PATRIC API stopped responding"
            return results
        results=results+cur_rows    
        start=start+rows
    return results
				

def get_tax_info(tax_ids, gids, target_path):
    query_pieces=[]
    currentQuery=taxURL+taxQuery
    if len(tax_ids):
        query_pieces.append(taxidQuery+"("+" OR ".join([str(i) for i in tax_ids])+")")
    if len(gids):
        query_pieces.append(gidQuery+"("+" OR ".join([str(i) for i in gids])+")")
    currentQuery=currentQuery+"+OR+".join(query_pieces)+taxFormat
    print currentQuery
    

    #get replicon results
    taxonomy_results=get_solr_result(currentQuery)
    out_handle=open(target_path+str(ts)+".taxinfo.txt",'w')
    #create header
    if len(taxonomy_results) and 'genome_name' in taxonomy_results[0]:
        print "Creating tax table from SOLR for "+currentQuery

    gids=set()
    for r in taxonomy_results:
        try:
            out_handle.write("\t".join([r['genome_name'],str(r['genome_info_id']),str(r['ncbi_tax_id']),','.join(r['taxon_lineage_ids'])])+"\n")
        except:
            print "couldn't write line "+str(r)
        if 'genome_info_id' in r:
            gids.add(str(r["genome_info_id"]))
    out_handle.close()
    return gids

def get_patric_feature_info(gids, target_path):
    out_handle=open(target_path+str(ts)+".feature_info.txt",'w')
    figfams=set()
    for g in gids:
        gidQuery="gid:"+g
        currentQuery=feature_url+gidQuery+feature_conditions+"&"+format_string
        print currentQuery
        feature_results=get_solr_result(currentQuery)

        for f in feature_results:
            try:
                out_handle.write("\t".join([str(f["gid"]),str(f["sequence_info_id"]),str(f["na_feature_id"]),str(f["start_max"]),str(f["end_min"]),f["figfam_id"]])+"\n")
		figfams.add(f["figfam_id"])
            except:
                print "couldn't write line "+str(f)
    out_handle.close()
    return figfams

def get_figfam_info(figfams, target_path):
    out_handle=open(target_path+str(ts)+".family_info.txt",'w')
    fid_query="figfam_id:("+" OR ".join(figfams)+")&fl=figfam_id,figfam_product"
    currentQuery=fam_url+fid_query+"&"+format_string
    figfam_results=get_solr_result(currentQuery)
    for f in figfam_results: 
        try:
            out_handle.write("\t".join([f["figfam_id"],str(f["figfam_product"])])+"\n")
        except:
            print "couldn't write line "+str(f)
    out_handle.close()
    
           

def main(init_args):
    target_path='./suis110514/'
    gids=get_tax_info(tax_ids=[29461],gids=[],target_path=target_path)
    figfams=get_patric_feature_info(gids, target_path)
    #table=pd.read_table('./1410652018.19.feature_info.txt',names=['figfam_id','gid','ncbi_tax_id','sid','start','end'])
    #figfams=set(table['figfam_id'])
    get_figfam_info(figfams, target_path)

if __name__ == "__main__":
	main(sys.argv[1:])												
