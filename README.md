# Turterra

Turterra is a portal for analysing protein families. It consists of two main parts: turterra, which runs a web portal from a folder tree, and turterra-build, which creates any files in the folder tree that may be missing from a .fasta file and a directory containing templates for homology modelling.

## Installation

Turterra and turterra-build are installed together as follows:

First, we recommend you install [Anaconda](https://www.anaconda.com/products/individual-b) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html). Then, create a new conda environment for turterra and activate it:

```sh
conda create -n turterra python=3.9
conda activate turterra
```

Next, clone the turterra repository into a location of your choice, navigate to the folder, and install turterra.

```sh
conda install -c bioconda epa-ng hmmer muscle fasttree
git clone https://github.com/TurtleTools/turterra.git
cd turterra
pip install .
```

The majority of turterra's dependencies are installed through the provided setup.py file. However, some dependencies will need to be installed through conda.

```sh
conda install -c bioconda epa-ng hmmer muscle fasttree
```

Congratulations! Turterra was installed!

## Turterra folder architecture

In order to run turterra with your own data, create a folder called 'data' in the top-level folder called turterra. This folder should contain the following files and folders:

```
turterra
    |--data
        |--data.txt
        |--sequences.fasta
        |--sequence_alignment.fasta
        |--smiles.tsv
        |--structure_alignment.fasta
        |--structures
            |--accession1_model.pdb
            |--accession2.pdb
            |--accession3_model.pdb
            |--...
        |--tree.txt
```

| file name | file contents |
| ------ | ------ |
| data.txt | tab-separated file, with categories in the first row and data for each sequence in the following rows. Any category can be defined. These are the categories that turterra will later be able to filter your data on. Currently, the categories 'Accession', 'Species' and 'Compounds' should always be present. |
| sequences.fasta | A .fasta file containing all the sequences in the analysis, with the accessions specified in data.txt as headers. |
| sequence_alignment.fasta | A .fasta file containing an alignment of all sequences in the analysis, with the accessions specified in data.txt as headers. |
| smiles.tsv | A tab-separated file, with as header 'Name\tSMILES', and in the remaining rows all compound names in the analysis and their corresponding structures in [SMILES format](http://opensmiles.org/opensmiles.html). |
| structure_alignment.fasta | A .fasta file containing a structure-based sequence alignment of all sequences in the analysis, with the accessions specified in data.txt as headers. We recommend you create this file with [caretta](https://github.com/TurtleTools/caretta) through turterra-build. |
| structures | Directory containing structural (homology) models for sequences in the analysis in .pdb format. File names should have the format 'accession_model.pdb' for homology-modelled structures, and 'accession.pdb' for crystal structures. Accessions should match the accessions in data.txt. |
| tree.txt | A phylogenetic tree in newick format. Leaf nodes should be labelled with the accessions specified in data.txt. |

All these files, with the exception of data.txt and smiles.tsv, can be automatically generated with turterra-build from a .fasta file containing the sequences you wish to analyse (use accessions of your choice as header), and a folder containing .pdb structures for homology modelling.

## Running turterra-build

Turterra-build uses [Muscle](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-113) (default settings) for sequence alignment, [FastTree](http://www.microbesonline.org/fasttree/) (default settings) for phylogenetic tree construction, [Modeller](https://salilab.org/modeller/) (250 models per sequence) for homology modelling, and Caretta for structure alignment. Turterra-build is run as follows:

```
Usage: turterra-build [OPTIONS]

Options:
  --fasta.                        [required] .fasta file containing all sequences
                                  to analyse, using accessions as headers.
  --build-models / --no-build-models
                                  Run modeller to construct homology models
                                  using MODELLER v10.0.  [default: False]

  --build-alignment / --no-build-alignment
                                  Build a multiple sequence alignment using
                                  Muscle v3.8.  [default: False]

  --alignment                     Multiple sequence alignment. Required if
                                  --build_alignment is false.  [default: '']

  --build-tree / --no-build-tree  Build phylogenetic tree using FastTree
                                  v2.1.10.  [default: False]

  --model-templates.              Directory containing templates for homology
                                  modelling in .pdb format. Required if
                                  --model is passed.  [default: '']

  --num-threads.                  Number of threads used for homology
                                  modelling.  [default: 1]

  --help                          Show this message and exit.
```

The result is a folder architecture that can be loaded into turterra by replacing the folder turterra/data with the generated folder, or replacing individual files or directories within turterra/data with files or directories generated by turterra-build.

## Running turterra

Once the correct folder architecture has been created, all the hard work is done! Now it is just a simple matter or running turterra from the command line. First, navigate to the top-level turterra folder. Then type:

```sh
bin/turterra
```



