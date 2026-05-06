
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import time, json
import subprocess
import argparse
from Bio import SeqIO
from Bio import SeqRecord

def parse_contig_order(contig_file):
    contig_to_genome={}
    genome_to_contig={}
    with open(contig_file, 'r') as istream:
        for line in istream:
            parts = line.strip().split("\t")
            if len(parts) >1:
                for c in parts[1:]:
                    genome_to_contig.setdefault(parts[0],parts[1:])
                    contig_to_genome.setdefault(c,parts[0])
    return contig_to_genome, genome_to_contig

def partition_fasta(contigs_ordered, contigs_unordered, c2g_ordered, c2g_unordered, fasta_files):
    align_pairs={}#query key, subject value
    ordered_fasta={}#genome id key, ordered fasta value
    unordered_fasta={}#genome id key, ordered fasta value
    ref_opt="genome" #genome or ordered seqs
    contigs={}
    fasta_lookup={}
    for f in fasta_files:
        f_unorder=f+".tmp.unordered"
        with open(f,'r') as fh, open(f_unorder,'w') as fh_unorder:
            missing_contigs=[]
            genome_id=None
            for record in SeqIO.parse(fh, "fasta"):
                contigs.setdefault(record.id, record)
                genome_id = c2g_ordered.get(record.id, c2g_unordered.get(record.id, "missing"))
                if genome_id == "missing":
                    sys.stderr.write("Fasta record %s in fasta file %s is missing from contigs processed by panaconda" % (record.id, f)+"\n")
                    missing_contigs.append(record.id)
                else:
                    fasta_lookup.setdefault(genome_id, {}).setdefault(f,[]).append(record.id)
                if record.id in c2g_ordered and genome_id != "missing":
                #    SeqIO.write(record, fh_order, "fasta")
                    if ref_opt == "genome":
                        ordered_fasta.setdefault(c2g_ordered.get(record.id),f)
                    elif ref_opt == "ordered":
                        ordered_fasta.setdefault(c2g_ordered.get(record.id),f_order)
                elif record.id in c2g_unordered or genome_id == "missing":
                    SeqIO.write(record, fh_unorder, "fasta")
                    unordered_fasta.setdefault(c2g_unordered.get(record.id),f_unorder)
            if len(missing_contigs) > 0:
                if genome_id == None:
                    sys.stderr.write("Error: No parseable contigs from entire genome "+genome_id+"\n")
                    sys.exit()
                else:
                    contigs_unordered.get(genome_id).extend(missing_contigs)
    for uid, uf in unordered_fasta.iteritems():
        align_pairs.setdefault(uid,(uf,[of for oid, of in ordered_fasta.iteritems() if oid != uid]))
    for gid, cur_lookup in fasta_lookup.iteritems():
        for f in cur_lookup:
            f_order=f+".tmp.ordered"
            with open(f_order,'w') as fh_order:
                for c_id in contigs_ordered.get(gid,[]):
                    cur_rec=contigs.get(c_id,None)
                    if type(cur_rec) == SeqRecord.SeqRecord:
                        SeqIO.write(cur_rec, fh_order, "fasta")
                for c_id in contigs_unordered.get(gid,[]):
                    cur_rec=contigs.get(c_id,None)
                    if type(cur_rec) == SeqRecord.SeqRecord:
                        SeqIO.write(cur_rec, fh_order, "fasta")

    return align_pairs, fasta_lookup
        



def place_contigs(align_groupings, synteny_groupings, contigs_ordered):
    return

def create_pan_params(top_params):
    ofolder=top_params.get("output_folder")
    json_file=top_params.get("panaconda_settings")
    feature_files=top_params.get("feature_files")
    fasta_files=top_params.get("contig_files")
    if type(feature_files) != list:
        feature_files=feature_files.strip().split(",")
        top_params["feature_files"]=feature_files
    if type(fasta_files) != list:
        fasta_files=fasta_files.strip().split(",")
        top_params["contig_files"]=fasta_files
    contig_files=top_params.get("contig_files")
    new_params={}
    new_params["contig_output"]= os.path.join(ofolder, "area_contig_order.txt")
    new_params["output"]= os.path.join(ofolder,"pc_output")
    params=[]
    with open(json_file,'r') as fh:
        param_dict=json.load(fh)
        param_dict.update(new_params)
        for k,v in param_dict.iteritems():
            params.append("--"+k)
            if v != None:
                params.append(v)
    params.extend(feature_files)
    new_params["unordered_contigs"]=new_params["contig_output"]+".unsorted"
    top_params.update(new_params)
    return params

def run_mauve_sort(fasta_lookup):
    min_key, min_value = min(fasta_lookup.iteritems(), key = lambda x: len(x[1].values()[0]))
    min_file = min_value.keys()[0]
    for k, cur_lookup in fasta_lookup.iteritems():
        sort_file=cur_lookup.keys()[0]
        if sort_file != min_file:
            try:
                command=["contig_sort","--reference",min_file,"--other",sort_file]
                print " ".join(command)
                subprocess.check_call(["nucmer","-p","1094551.3_v_"+str(c)+"_nucmer", g, x])
            except Exception as err:
                sys.stderr.write("nucmer failed: %s %s\n" % (" ".join(["nucmer","-p","1094551.3"+"nucmer", g, x]), err))
                sys.exit(1)

    
def run_mummer():
    for x,y in [alignment_pairs.get('1094551.3')]:
        for c,g in enumerate(y):
            try:
                print " ".join(["nucmer","-p",os.path.join(ofolder, "1094551.3"+"_nucmer"), g, x])
                subprocess.check_call(["nucmer","-p","1094551.3_v_"+str(c)+"_nucmer", g, x])
            except Exception as err:
                sys.stderr.write("nucmer failed: %s %s\n" % (" ".join(["nucmer","-p","1094551.3"+"nucmer", g, x]), err))
                sys.exit(1)

def main(params):
    ofolder=params.get("output_folder")
    #prepare parameters
    pparams=create_pan_params(params)
    #invoke panaconda to order contigs, partition unordered
    print " ".join(["python","fam_to_graph.py"]+pparams)
    try:
        subprocess.check_call(["python","fam_to_graph.py"]+pparams)
    except Exception as err:
        sys.stderr.write("panaconda failed: %s %s\n" % (" ".join(pparams), err))
        sys.exit(1)
    #those not visited are partitioned from those that are
    contig_file=params.get("contig_output")
    unordered_file=params.get("unordered_contigs")
    c2g_ordered, contigs_ordered = parse_contig_order(contig_file)
    c2g_unordered, contigs_unordered = parse_contig_order(unordered_file)
    alignment_pairs, fasta_lookup= partition_fasta(contigs_ordered, contigs_unordered, c2g_ordered, c2g_unordered, params.get("contig_files"))
    run_mauve_sort(fasta_lookup)
    sys.exit()
    #align using mummer4, exclude own genome
    alignment_files = make_alignments(alignment_pairs)
    #for each alignment get the best
    align_groupings = parse_alignments(alignment_files)
    #need to know what contigs are paired 
    place_contigs(align_groupings, synteny_groupings, contigs_ordered)

if __name__ == "__main__":
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=str, dest="output_folder", default="", help="Output folder for the sorted fasta.")
    parser.add_argument("--feature_files", type=str, nargs='*', default="", help="CSV files of varying format specifing group, genome, contig, feature, and start in sorted order. stdin also accepted")
    parser.add_argument("--contig_files", type=str, nargs='*', default="", help="CSV files of sequences for contigs. Collated with feature_files. IDs of contigs should match.")
    parser.add_argument("--panaconda_settings", type=str, default="./default_settings.json", help="JSON formatted parameters for panaconda. null is given as value for key that should be flag.")
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()
    pargs = parser.parse_args()
    params=vars(pargs)
    main(params)
    sys.stderr.write(("--- %s seconds ---" % (time.time() - start_time))+"\n")
