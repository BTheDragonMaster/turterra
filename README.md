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
            |--sequence1_model.pdb
            |--sequence2.pdb
            |--sequence3_model.pdb
            |--...
        |--tree.txt
```

| file name | file contents |
| ------ | ------ |
| data.txt | tab-separated file, with categories in the first row and data for each sequence in the following rows. Any category can be defined. These are the categories that turterra will later be able to filter your data on. Currently, the categories 'Accession', 'Species' and 'Compounds' should always be present. |
| sequences.fasta | A .fasta file containing all the sequences in the analysis, with the accessions specified in data.txt as headers. |
| sequence_alignment.fasta | A .fasta file containing an alignment of all sequences in the analysis, with the accessions specified in data.txt as headers. |
| smiles.tsv | A tab-separated file, with as header 'Name\tSMILES', and all compounds names in the analysis and their corresponding structures in [SMILES format](http://opensmiles.org/opensmiles.html). |





