# ska-mtdna

The ska-mtdna.py script leverages SKA and FastTree to allow for rapid and automated haplotype analysis of mitochondrial genomes. 


This script is described in *Rapid and reproducible haplotyping of complete mitochondrial genomes using split k-mers* Stuehler Jr. et al.



Introduction
------------------
Ska-mtdna haplotype and phylogenetic assesments of mitogenomes utilizes boostrap resampling from FastTree and compositie scoring from adjustable alignment metrics. 
SKA, intoduced [here](https://github.com/simonrharris/SKA) is the main component under the hood of ska-mtdna which breaks down larger sequences into subsequences of a specified length.
Haplotyping mitochondrial genomes with split *k*-mers is a new method comparable to MLST or WGA methods, presenting major advantages such as detection of polymorphisms in intergenic sequence in a reproduclible manner.
The script is run in two modes *network* and *phylo* which prioritize haplotype network construction or phylogenetic tree construction - test datasets are available for each mode.



Setup
------------------
The ska-mtdna.py script is released as a Python script for download, tested with Python 3.10.14.

Before you begin, you will need to download and install the required software, which are installed on the command line most simply with the [Conda](https://docs.conda.io/en/latest/) package manager.
Software descriptions can be found at the following links:

-Python 3 
Preinstalled in most linux distributions or installed into Conda environment.

- SKA2 
https://github.com/bacpop/ska.rust

- FastTree
https://github.com/morgannprice/fasttree

To get ska-mtdna.py setup, setup [Conda](https://docs.conda.io/en/latest/), then run these commands in order:

1. Create and activate conda environment:

`conda create -n ska_mtdna python=3.10`

`conda activate ska_mtdna`


2. Install SKA2

`conda install bioconda::ska2`


3. Install FastTree

`conda install bioconda::fasttree`



Running ska-mtdna.py:
------------------

Make sure the necessary software is installed and accessible from the command-line, run `ska-mtdna.py --help` to print program information.


1. Network mode

python3 ska-mtdna.py --fasta-dir fasta_input_folder

2. Phylo mode

python3 ska-mtdna.py --fasta-dir fasta_input_folder --phylo



Options:
------------------

1. Mandatory arguments

```
--fasta-dir Directory containing *.fa or *.fasta files
```

2. Optional arguments
```
-h, --help            show this help message and exit
  --network             Network mode (default)
  --phylo               Phylogenetic mode
  -k K                  Comma-separated k-mer sizes: default=11,13,15,17,19,21,23,25,27,29
  -m M                  Comma-separated m values: 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1
  --output-dir OUTPUT_DIR Directory for output files
  --num_top_params NUM_TOP_PARAMS Number of top-scoring k,m parameter sets to examine with ska map (default: 5, network mode only)
  --wb WB               Weight for bootstrap score (network: 0.87, phylo: 0.87)
  --wh WH               Weight for haplotypes (network: 0.06, phylo: 0.06)
  --wn WN               Weight for valid nucleotide sites (network: 0.07, phylo: 0.07)
  --ws WS               Weight for segregating sites (network: 0.0, phylo: 0.0)
  --wg WG               Weight for gap proportion (network: 0.0, phylo: -0.1)
```

Output:
------------------

A successful run of ska-mtdna.py will result in a directory output_repeats_masked by default containing:

1. Files
`initial_results.csv`  Results of the *ska align* algorithm (network and phylo).

`final_results.csv`  Results of the *ska map* algorithm (network only). Number of paremeter sets reported by `--num_top_params`.

`bootstrap_score_heatmap.png`  Heatmap of boostrap results of the *ska align* algorithm (network and phylo).

`composite_score_heatmap.png`  Heatmap of composite score results from weights (mode specific).

`bootstrap_vs_haplotypes.png`  Scatter plot of boostrap values vs. number of haplotypes (network).

`*.skf` k-mer profiles built by ska build

`*.aln` unordered k-mer sequence alignments written by ska align for a given parameter set 

`*.distance` pairwise genetic distances for all samples for a given parameter set 

`*.tree` phyogenetic trees written by FastTree for a given parameter set 

`*map*.aln` ordered k-mer sequence alignments written by ska map for a given parameter set

`*.nex` ordered k-mer sequence alignments converted from .aln to .nex for a given parameter set



Citation:
------------------

Rapid and reproducible haplotyping of complete mitochondrial genomes using split k-mers
Douglas S. Stuehler Jr., Liliana M. Cano, Michelle Heck
bioRxiv 2025.03.23.644767; doi: https://doi.org/10.1101/2025.03.23.644767


