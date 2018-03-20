# RaGOO

### A tool to order and orient genome assembly contigs via minimap2 alignments to a reference genome.

## Description

RaGOO is a tool for coalescing genome assembly contigs into pseudochromosomes via minimap2 alignments to a closely related reference genome. The focus of this tool is on practicality and therefore has the following features:

1. Good performance. On a MacBook Pro using Arabidopsis data, pseudochromosome construction takes less than a minute and the whole pipeline with SV calling takes ~ 2 minutes.
2. In-tact ordering and orienting of contigs so as to retain original sequence content and any analysis that relies on it. 
3. Lift-over of gff files from the contig-level assembly to the pseudochromosomes.
4. Structural variant calling with [Assemblytics](http://assemblytics.com/).

RaGOO does not automatically create synteny plots/homology maps, but RaGOO output files can easily passed to visualization tools such as Assemblytics (see "Output Files" below).

## Installation

### Dependencies

RaGOO should install on OSX and most standard flavors of Linux. RaGOO depends on Python3 as well as the following packages:

1. [intervaltree](https://pypi.python.org/pypi/intervaltree)
2. numpy
3. [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
4. [minimap2](https://github.com/lh3/minimap2)

The first two packages will be installed automatically when installing RaGOO. The last two are straightforward to install following the instructions on their respective websites. RaGOO expects bedtools to be installed globally and for the bedtools executables to be in ones path, however, the location of the minimap2 executable can optionally passed to RaGOO if not installed globally (see "Usage" below).  

### Installation

Currently, the only way to install RaGOO is from source. Set up a virtualenv if desired, just be sure to make a python3 environment. Then, enter the following command to install RaGOO and its dependencies:

$python setup.py install

## Usage

```
usage: ragoo.py [-h] [-e <exclude.txt>] [-gff <annotations.gff>] [-m PATH]
                [-b] [-t 3] [-g 100] [-s]
                <contigs.fasta> <reference.fasta>

order and orient contigs according to minimap2 alignments to a reference

positional arguments:
  <contigs.fasta>       fasta file with contigs to be ordered and oriented
  <reference.fasta>     reference fasta file

optional arguments:
  -h, --help            show this help message and exit
  -e <exclude.txt>      single column text file of reference headers to ignore
  -gff <annotations.gff>
                        Annotations for lift over
  -m PATH               path to minimap2 executable
  -b                    Break chimeric contigs
  -t 3                  Number of threads when running minimap.
  -g 100                Gap size for padding in pseudomolecules.
  -s                    Call structural variants
```

Note that one can optionally break chimeric contigs and call structural variants. Turning off SV calling can be especially useful for large genomes such as Maize and Human for which minimap2 SAM alignments might take a very long time. 

RaGOO will try to be smart and not redo intermediate analysis already done in previous executions of the pipeline. For example, if the minimap2 alignment files are already present from previous runs, RaGOO will not recreate them. However, RaGOO is not that smart, so be sure to remove any files that you want to replace. To be safe, one can just remove the entire output directory if a new analysis is desired (see "Output Files" below).

### Output Files
