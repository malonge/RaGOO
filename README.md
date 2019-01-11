# RaGOO

### A tool to order and orient genome assembly contigs via Minimap2 alignments to a reference genome.

## Description

RaGOO is a tool for coalescing genome assembly contigs into pseudochromosomes via minimap2 alignments to a closely related reference genome. The focus of this tool is on practicality and therefore has the following features:

1. High performance. On a MacBook Pro using Arabidopsis data, pseudochromosome construction takes less than a minute and the whole pipeline with SV calling takes ~2 minutes.
2. Intact ordering and orienting of contigs. 
3. Chimeric contig correction.
4. Structural variant calling with [Assemblytics](http://assemblytics.com/).
5. Confidence scores associated with the grouping, localization, and orientation for each contig.

## Installation

### Dependencies

RaGOO should install on OSX and most standard flavors of Linux. RaGOO depends on Python3 as well as the following packages:

1. [intervaltree](https://pypi.python.org/pypi/intervaltree)
2. numpy
3. [Minimap2](https://github.com/lh3/minimap2)

The first two packages will be installed automatically when installing RaGOO. Minimap2 is straightforward to install following the instructions on its website. Place the minimap2 executable in your path, or specify its location with the -m paramter (see below).

### Installation

Currently, the only way to install RaGOO is from source. Set up a virtualenv if desired, just be sure to make a python3 environment. Then, enter the following command to install RaGOO and its dependencies:

```
$ python setup.py install
```

## Usage

```
usage: ragoo.py [-h] [-e <exclude.txt>] [-m PATH] [-b] [-t 3] [-g 100] [-s]
                <contigs.fasta> <reference.fasta>

order and orient contigs according to minimap2 alignments to a reference

positional arguments:
  <contigs.fasta>    fasta file with contigs to be ordered and oriented
  <reference.fasta>  reference fasta file

optional arguments:
  -h, --help         show this help message and exit
  -e <exclude.txt>   single column text file of reference headers to ignore
  -m PATH            path to minimap2 executable
  -b                 Break chimeric contigs
  -t 3               Number of threads when running minimap.
  -g 100             Gap size for padding in pseudomolecules.
  -s                 Call structural variants
``` 

RaGOO will try to be smart and not redo intermediate analysis already done in previous executions of the pipeline. For example, if the minimap2 alignment files are already present from previous runs, RaGOO will not recreate them. However, RaGOO is not that smart, so be sure to remove any files that you want to replace. To be safe, one can just remove the entire output directory if a new analysis is desired (see "Output Files" below).

### Example Run
Both the assembly and the reference must be in the current workding directory, so please either copy them or create a symbolic link. For example:

```
$ cd /path/to/current/working/directory
$ ln -s /path/to/contigs.fasta
$ ln -s /path/to/reference.fasta
$ ragoo.py contigs.fasta reference.fasta
```

### Output Files

All of the output will be in the "ragoo_output" directory. If breaking chimeric contigs and calling SVs, the contents of this output directory is as follows:

```
ragoo_output/
├── ragoo.fasta
├── chimera_break/
├── groupings/
├── orderings/
└── pm_alignments/
```

#### ragoo.fasta
The final pseudomolecules. Any unlocalized contigs are concatenated and placed in "Chr0_RaGOO".

#### chimera_break
This directory cotnains the results from chimeric contig breaking. It contains two alignment (PAF) files of contigs against the reference, one for intrachromosmal and interchromosmal chimeras. A new corrected fasta file is also associated with each alignment file.

The most notable file here is the **[prefix].intra.chimera.broken.fa**, as this is the final corrected assembly used for downstream scaffolding. All of the downstream information such as confidence scores refers to this assembly, not the orignal assembly.

#### groupings
There is one file per chromosome listing the contigs assigned to that chromosome, and their grouping confidence score. Please not that these contigs are not ordered. Also note that if chimeras were corrected, the headears in these files refer to the broekn assembly in "chimera_break", and not the original assembly.

#### orderings
There is one file per chromosome showing the ordering, orientation (second column), location confidence scores (third column), and orientation confidence scores (fourth column) for each chromosome.

#### pm_alignments
This directory contains all of the structural variant calling results. The final structural variants can be found in **assemblytics_out.Assemblytics_structural_variants.bed**. There is also the alignmens of the RaGOO pseudomolecules against the provided reference in both SAM and delta format (**pm_contigs_against_ref.sam** and **pm_contigs_against_ref.sam.delta**). 
