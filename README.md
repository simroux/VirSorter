# VirSorter

Source code of the VirSorter App, available on CyVerse (https://de.iplantcollaborative.org/de/)

# Publication

* VirSorter: mining viral signal from microbial genomic data
* https://peerj.com/articles/985/
* PubMed 26038737

# Result files
The main output files of VirSorter are:
- VIRSorter_global-phage-signal.csv: Comma-separated table listing the viral predictions from VirSorter (one row per prediction)
- Metrics_files/VIRSorter_affi-contigs.tab: Pipe("|")-delimited table listing the annotation of all predicted ORFs in all contigs. Lines starting with a ">" are "headers", i.e. information about the contig (contig name, number of genes, "c" for circular or "l" for linear). All other lines are information about the genes, with different columns as follows: Gene name, start, stop, length, strand, Hit in the virus protein cluster database, hit score, hit e-value, category of the virus protein cluster (see below), Hit in PFAM, hit score, hit e-value.
The categories of virus clusters represent the range of genomes in which this virus cluster was detected, i.e. 0: hallmark genes found in Caudovirales, 1: non-hallmark gene found in Caudovirales, 2: non-hallmarke gene found exclusively in virome(s), 3: hallmark gene not found in Caudovirales, 4: non-hallmark gene not found in Caudovirales
- Predicted_viral_sequences/: fasta and genbank files of predicted viral sequences
- Fasta_files/: intermediary files, including predicted proteins
- Tab_files/: intermediary files, including results of the search agasint PFAM and the virus database.

VirSorter results can be imported into [Anvi'o](http://merenlab.org/software/anvio/) by following [these instructions](http://merenlab.org/2018/02/08/importing-virsorter-annotations/).

# Using a conda virtual environment (tested on Ubuntu and CentOS)
* First install [Anaconda or Miniconda](https://conda.io/docs/user-guide/install/index.html)
* Download the databases required by VirSorter which have been converted to be used with HMMER version 3.1b2. Change to the directory where you want the databases be, and then run the following commands:
```
wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
md5sum virsorter-data-v2.tar.gz
#m5sum should return dd12af7d13da0a85df0a9106e9346b45
tar -xvzf virsorter-data-v2.tar.gz
```
* Create and install your conda virtual environment. Change to the directory where you want VirSorter to be installed and run the following commands:
```
conda create --name virsorter -c bioconda mcl=14.137 muscle blast perl-bioperl perl-file-which hmmer=3.1b2 perl-parallel-forkmanager perl-list-moreutils diamond
git clone https://github.com/simroux/VirSorter.git
cd VirSorter/Scripts
make
```
* To run VirSorter from any directory, you can make symbolic links to `VirSorter/wrapper_phage_contigs_sorter_iPlant.pl` and `VirSorter/Scripts` and place them in the `bin` folder for your "virsorter" conda environment. An example location of this `bin` folder is `~/miniconda/envs/virsorter/bin`. Substitute this path with the path to the `bin` folder for your newly created "virsorter" environment.
```
ln -s ~/Applications/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl ~/miniconda/envs/virsorter/bin
ln -s ~/Applications/VirSorter/Scripts ~/miniconda/envs/virsorter/bin
```
* Finally, you'll need to download MetaGeneAnnotator ([Noguchi et al, 2006](https://doi.org/10.1093/nar/gkl723)). You can put this directly in the "virsorter" environment's `bin` folder alongside the VirSorter symbolic links taht were just created.
```
cd ~/miniconda/envs/virsorter/bin
wget http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz
tar -xvzf metagene/mga_x86_64.tar.gz
```

To run VirSorter, type the following:

```
source activate virsorter
wrapper_phage_contigs_sorter_iPlant.pl -f assembly.fasta --db 1 --wdir output_directory --ncpu 4 --data-dir /path/to/virsorter-data
```

# Docker - from DockerHub

* Download the databases required by VirSorter, available as a tarball archive on iMicrobe: http://mirrors.iplantcollaborative.org/browse/iplant/home/shared/imicrobe/VirSorter/virsorter-data.tar.gz
or /iplant/home/shared/imicrobe/VirSorter/virsorter-data.tar.gz through iPlant Discovery Environment
* Untar this package in a directory, e.g. /host/path/to/virsorter-data
* Pull VirSorter from dockerhub: $ docker pull discoenv/virsorter:v1.0.3
* Create a working directory for VirSorter which includes the input fasta file, e.g. /host/path/to/virsorter-run
* Then run VirSorter from docker, mounting the data directory as data, and the run directory as wdir:

    $ docker run -v /host/path/to/virsorter-data:/data -v /host/path/to/virsorter-run:/wdir -w /wdir --rm discoenv/virsorter:v1.0.3 --db 2 --fna /wdir/Input_contigs.fna

After "virsorter:v1.0.3", the options correspond to the ones described in wrapper_phage_contigs_sorter_iPlant.pl (here selecting the database "Viromes" and pointing VirSorter to the file "Input_contigs.fna").


# Docker - building packages from scratch


## Dependencies

Install the following into a "bin" directory:

* HMMER (http://hmmer.janelia.org/)
* MCL (http://micans.org/mcl/)
* Metagene Annotator (http://metagene.nig.ac.jp/metagene/download_mga.html)
* MUSCLE (http://www.drive5.com/muscle/)
* BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* DIAMOND (https://github.com/bbuchfink/diamond)

## Data Container

The 12G of dependent data exists as a separate data container 
called "virsorter-data."

This is the Dockerfile for that:

    FROM perl:latest

    MAINTAINER Ken Youens-Clark <kyclark@email.arizona.edu>

    COPY Generic_ref_file.refs /data/

    COPY PFAM_27 /data/PFAM_27

    COPY Phage_gene_catalog /data/Phage_gene_catalog

    COPY Phage_gene_catalog_plus_viromes /data/Phage_gene_catalog_plus_viromes

    COPY VirSorter_Readme.txt /data

    COPY VirSorter_Readme_viromes.txt /data

    VOLUME ["/data"]
  
Then do:

    $ docker build -t kyclark/virsorter-data .
    $ docker create --name virsorter-data kyclark/virsorter-data /bin/true

## Build

    $ docker build -t kyclark/virsorter .

## Run

A sample "run" command to use the current working directory for input/output:

    $ docker run --rm --volumes-from virsorter-data -v $(pwd):/de-app-work \
    -w /de-app-work kyclark/virsorter --fna Mic_1.fna --db 1

# Authors

Simon Roux <sroux@lbl.gov> is the author of Virsorter

Ken Youens-Clark <kyclark@email.arizona.edu> packaged this for Docker/iPlant.

Bryan D Merrill <bmerrill@stanford.edu> provided the improvements and additions for v1.0.5
