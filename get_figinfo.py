#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Wed Dec 18 22:04:38 2013

@author: anwarren
"""

import requests, json, sys




taxURL="http://macleod.vbi.vt.edu:8080/solr/genomesummary/select/?"
taxQuery="q=rast_cds:[1+TO+*]+AND+genome_info_id:(38055+OR+135566+OR+86222)"
taxFormat="&indent=on&wt=json&fl=genome_name,genome_info_id,ncbi_tax_id,taxon_lineage_ids"

feature_url='http://macleod.vbi.vt.edu:8080/solr/figfam/select/?'

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
				

def patric_figinfo_from_solr(tax_id, target_tax_path, target_figfam_path):
    format_string="indent=on&wt=json"
    
    currentQuery=taxURL+taxQuery+taxFormat
    

    #get replicon results
    taxonomy_results=get_solr_result(currentQuery)
    out_handle=open(target_tax_path,'w')
    #create header
    if len(taxonomy_results) and 'genome_name' in taxonomy_results[0]:
        print "Creating tax table from SOLR for "+currentQuery

    gidQuery="gid:"
    gids=set()
    for r in taxonomy_results:
        try:
            out_handle.write("\t".join([r['genome_name'],str(r['genome_info_id']),str(r['ncbi_tax_id']),','.join(r['taxon_lineage_ids'])])+"\n")
        except:
            print "couldn't write line "+str(r)
        if 'genome_info_id' in r:
            gids.add(str(r["genome_info_id"]))
    out_handle.close()
    out_handle=open(target_figfam_path,'w')
    
    for g in gids:
        gidQuery="gid:"+g
        currentQuery=feature_url+"q="+gidQuery+"&"+"sort=sequence_info_id+asc,start_max+asc"+"&"+format_string
        print currentQuery
        feature_results=get_solr_result(currentQuery)

        for f in feature_results:
            try:
                out_handle.write("\t".join([f["figfam_id"],str(f["gid"]),str(f["ncbi_tax_id"]),str(f["sequence_info_id"]),str(f["start_max"]),f["figfam_product"]])+"\n")
            except:
                print "couldn't write line "+str(f)
            

def main(init_args):
    patric_figinfo_from_solr(29461, './abortus_brucella_taxon.txt', './abortus_brucella_figfam.txt')

if __name__ == "__main__":
	main(sys.argv[1:])												
