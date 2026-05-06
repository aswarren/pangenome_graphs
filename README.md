# Panaconda

**Application of pan-synteny graph models to genome content analysis.**

Panaconda is a graph-based bioinformatics tool for comparative microbial genomics. It bridges the gap between high-level gene presence/absence analyses (like Roary) and ultra-high-resolution nucleotide graphs (like PGGB). 

By operating at the **"meso-scale"**—using an alphabet of annotated Gene Families (such as BV-BRC/PATRIC `pgfams`, `plfams`, or `figfams`) rather than nucleotides—Panaconda abstracts away the noise of point mutations and small indels. It produces clean, topologically sound **Pan-Synteny Graphs** that instantly highlight genome architecture, evolutionary hotspots, and structural variants across thousands of genomes.

**Read the original paper:**[bioRxiv: Panaconda (Warren et al., 2017)](https://doi.org/10.1101/215988)

## Key Features

*   **Meso-Scale Graph Topology:** Builds a de Bruijn graph of gene families and compresses it into a Pan-Synteny (PS) graph, perfectly preserving sequence continuity.
*   **Explicit Structural Variant (SV) Detection:** Topologically identifies and tags:
    *   **Inversions** (Path loop-backs)
    *   **Translocations** (Cross-contig paralogy and path divergence)
    *   **CNVs / Tandem Repeats** (Localized paralogous expansions via Alt-Node clusters)
*   **Modern Pangenome Ecosystem Integration:** Exports natively to **GFA v1.1**, making the graphs compatible with modern downstream tools like `odgi`, `Bandage`, and `vg`.
*   **Interactive Visualization:** Exports to **GEXF** for beautifully laid-out force-directed visualization in Gephi or browser-based Gexf-JS viewers.
*   **Direct API Integration:** Can dynamically stream feature annotations directly from the[BV-BRC (PATRIC) API](https://www.bv-brc.org/) without needing to download massive local FASTA/GFF files.

## Installation

Panaconda requires **Python 3.7+**. 

Method 1:

Clone the repository and install the required dependencies:
```bash
git clone https://github.com/aswarren/pangenome_graphs.git
cd pangenome_graphs
pip install networkx requests

Method 2:
pip install git+https://github.com/aswarren/pangenome_graphs

To use the --layout feature, download gexf_layout.jar from 
https://github.com/aswarren/pangenome_layout/releases 
and place it in the src directory if running from repo folder structure 
or set export PANACONDA_LAYOUT_JAR=/path/to/jar

Usage

Panaconda accepts flat tabular feature files or can pull directly from the
BV-BRC API.

Basic Example (Local Files)

Generate a GFA sequence graph and a visual GEXF graph from a set of BV-BRC
feature tables:

python fam_to_graph.py \
    --alpha pgfam_id \
    --gfa output.gfa \
    --output output.gexf \
    ./data/Brucella/*.tab

Direct API Example (Remote Fetching)

If you have a text file of Genome IDs (input_genomes.txt), you can fetch the
features dynamically:

python fam_to_graph.py \
    --alpha pgfam_id \
    --patric_genomes \
    --gfa output.gfa \
    --output output.gexf \
    input_genomes.txt

Advanced Visualization (--layout)

If you want Panaconda to pre-compute a force-directed layout for your GEXF file
(useful for immediate browser rendering), include the --layout flag (Requires
Java):

python fam_to_graph.py --layout --alpha pgfam_id --gfa out.gfa --output out.gexf ./data/*.tab

Output Formats

1. Graphical Fragment Assembly (GFA v1.1)

A standard for pangenomics.

  - Segments (S): Represent syntenic blocks (Gene Families). Tags include Family
    ID (fm), Function (fn), Diversity Quotient (dv), and visual length (LN).
  - Links (L): Represent syntenic connections. Tagged with exact weights (wc,
    gc) and SV breakpoints (iv:i:1 for inversions, tr:i:1 for translocations).
  - Walks (W): Perfectly preserves the full genomic architecture of every
    sequence, allowing downstream tools to query exactly which strains traverse
    which paths.

The .gfa file into Bandage for instant, interactive 2D viewing.

2. Graph Exchange XML Format (GEXF)

Optimized for the Gephi layout engine and frontend JS viewers. By default,
bidirectional paths (like inversions) are merged into undirected edges to
prevent UI layout engines from dropping edges. Contains a <meta> JSON block
with a graph manifest and global SV statistics.

How it Works

Panaconda uses a Thread First Search (TFS) algorithm over a reverse-complement
aware de Bruijn graph of gene families.

1.  It identifies Anchor Nodes (highly conserved, non-duplicated k-mers).
2.  It traverses the threads (genomes) connecting these anchors, dynamically
    collapsing homologous paths into syntenic blocks.
3.  When threads conflict (e.g., due to tandem repeats), it elegantly emits
    macro-node clusters to preserve sequence continuity without breaking the
    graph.
4.  It tags breakpoints where the sequence topology violates expected genomic
    contexts (detecting translocations and inversions).

License

This project is licensed under the MIT License - see the LICENSE file for
details.

Citation

If you use Panaconda in your research, please cite:

Warren, A. S., Davis, J. J., Wattam, A. R., Machi, D., Setubal, J. C., & Heath,
L. S. (2017). Panaconda: Application of pan-synteny graph models to genome
content analysis. bioRxiv, 215988.https://doi.org/10.1101/215988

