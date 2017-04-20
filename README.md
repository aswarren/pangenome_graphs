usage: fam_to_graph.py [-h] [--no_function] [--layout] [--output OUTPUT]  
                       [--rfgraph RFGRAPH] [--diversity {genus,species}]  
                       [--patric_figfam | --patric_plfam | --patric_pgfam | --generic]  
                       [--context {genome,contig,feature}]  
                       [--ksize {3,4,5,6,7,8,9}]  
                       [feature_files [feature_files ...]]  
positional arguments:  
  feature_files         Files of varying format specifing group, genome,  
                        contig, feature, and start in sorted order. stdin also  
                        accepted  
optional arguments:  
  -h, --help            show this help message and exit  
  --no_function         No functions as labels. Keep file size smaller.  
  --layout              run gephi layout code for gexf  
  --output OUTPUT       the path and base name give to the output files. if  
                        not given goes to stdout  
  --rfgraph RFGRAPH     create rf-graph gexf file at the following location  
  --diversity {genus,species}  
                        calculate diversity quotient according to given taxa  
                        level  
  --patric_figfam       PATRIC feature file in tab format  
  --patric_plfam        PATRIC feature file in tab format  
  --patric_pgfam        PATRIC feature file in tab format. selecting pgfams  
  --generic             table specifying the group, genome, contig, feature,  
                        and start in sorted order  
  --context {genome,contig,feature}  
                        the synteny context  
  --ksize {3,4,5,6,7,8,9}  
                        the size of the kmer to use in constructing synteny  
