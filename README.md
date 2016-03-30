# VirSorter

Source code of the VirSorter App, available on iPlant (https://de.iplantcollaborative.org/de/)

# Publication

* VirSorter: mining viral signal from microbial genomic data
* https://peerj.com/articles/985/
* PubMed 26038737

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
* BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/, not BLAST+)


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

    COPY SUP05_SAGs_with_viruses.fna /data/

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

Simon Roux <roux.8@osu.edu> is the author of Virsorter

Ken Youens-Clark <kyclark@email.arizona.edu> packaged this for Docker/iPlant.
