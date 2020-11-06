
# Panaconda


## Requirements

```
//core code
python 2.7 (for now)
networkx


// layout algorithm
java

// visualization
firefox for web based browsing
gephi for graph editing and manipulation
```


## Installation

```
git clone --recursive https://github.com/aswarren/pangenome_graphs.git
cd pangenome_graphs
pip install -r requirements.txt
```


## Running

```
usage: fam_to_graph.py [-h] [--no_function] [--order_contigs {none,area,tfs}]
                       [--contig_output CONTIG_OUTPUT] [--layout]
                       [--output OUTPUT] [--rfgraph RFGRAPH]
                       [--diversity {genus,species}]
                       [--alpha {figfam_id,plfam_id,pgfam_id}]
                       [--patric | --patric_genomes]
                       [--context {genome,contig,feature}]
                       [--ksize {3,4,5,6,7,8,9}] [--min MIN]
                       [feature_files [feature_files ...]]

positional arguments:
  feature_files         Files of varying format specifing group, genome,
                        contig, feature, and start in sorted order. stdin also
                        accepted

optional arguments:
  -h, --help            show this help message and exit
  --no_function         no functions as labels. keep file size smaller.
  --order_contigs {none,area,tfs}
                        produce output that orders contigs for rectilinear
                        layout
  --contig_output CONTIG_OUTPUT
                        output file that orders contigs for rectilinear layout
  --layout              run gephi layout code for gexf
  --output OUTPUT       the path and base name give to the output files. if
                        not given goes to stdout
  --rfgraph RFGRAPH     create rf-graph gexf file at the following location
  --diversity {genus,species}
                        calculate diversity quotient according to given taxa
                        level
  --alpha {figfam_id,plfam_id,pgfam_id}
                        alphabet i.e. 'family type' to use
  --patric              table specifying the group, genome, contig, feature,
                        and start in sorted order
  --patric_genomes      use the files listed in --feature_files as a comma or
                        tab separated file specifying genome ids to pull from
                        patric. automatically downloads and uses the data
                        stream for those genome ids.
  --context {genome,contig,feature}
                        the synteny context
  --ksize {3,4,5,6,7,8,9}
                        the size of the kmer to use in constructing synteny
  --min MIN             minimum required sequences aligned to be in the
                        resulting graph

```
#### Example run for creating a graph
python fam_to_graph.py --layout --output data/BrucellaInversion/test_psgraph.gexf --alpha pgfam_id ./data/BrucellaInversion/*.tab

#### Visualizing data
Resulting gexf files can be opened in Gephi or through the JS visualizer distributed with Panaconda.
Files can be loaded from local disk but Chrome currently restricts this. To run the javascript based visualizer
locally you can use python to host a webserver and use a URL to view the data.

To do this:

cd viewer/gexf-js/

Soft link or copy a gexf file you want to view into the gexf-js folder. e.g. ln -s ../../data/BrucellaInversion/psgraph.gexf ./brucellainversion.gexf

python -m SimpleHTTPServer 8080

In firefox navigate to http://localhost:8080/index.html#brucellainversion.gexf



## Data

Examples from the paper https://www.biorxiv.org/content/early/2017/11/08/215988

Zipped versions can be found in the data directory.

These data can also be found and manipulated at PATRIC BRC (currently requires free account) at the following 
https://patricbrc.org/workspace/public/panaconda@patricbrc.org/Panaconda/PanSyntenyExamples

Currently the most conveniently accessible supported format is PATRIC's feature tab format.
Groups can be downloaded from the "feature tab" in PATRIC.  


