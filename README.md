Amplicon Clustering
===================

Performs *de novo* clustering and consensus generation from full-length amplicon sequence reads, including raw PacBio and Oxford Nanopore data


Installation
------------

Several common scientific computing packages for Python are required:

* matplotlib
* seaborn
* numpy
* scipy
* scikit-bio

The recommended method for building this environment is with Anaconda (or Miniconda):

    conda create -n amp python=3.6 anaconda
    conda activate amp
    conda install -c https://conda.anaconda.org/biocore scikit-bio


Usage
-----

    python denovo_amp.py -h

    usage: De novo clustering and consensus of long (full-length) amplicon sequences
    [-h] [--verbosity VERBOSITY] [--start START] [--end END] fa ref aln out

    positional arguments:
      fa                    Read FASTA file
      ref                   Initial consensus or ref FASTA
      aln                   Initial consensus alignment (m5)
      out                   Output prefix (incl. path)

    optional arguments:
      -h, --help            show this help message and exit
      --verbosity VERBOSITY
                            Level of output (0-3)
      --start START         Start position of ref to consider
      --end END             End position of ref to consider


Example
-------

    blasr reads.fasta ref.fasta  -m 5 -bestn 1 -out aln.m5
    python denovo_amp.py read.fasta ref.fasta aln.m5 denovo_clusters
