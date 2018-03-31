FROM ubuntu:14.04

LABEL maintainer="Ken Youens-Clark <kyclark@email.arizona.edu>, Simon Roux <siroux1@gmail.com>"

## Prepare the environment variables
ENV PATH=/miniconda/bin:${PATH} PERL5LIB=/miniconda/lib/perl5/site_perl/5.22.0/:${PERL5LIB}

## Copying the files
COPY wrapper_phage_contigs_sorter_iPlant.pl /usr/local/bin/
COPY Scripts /usr/local/bin/Scripts/

## Instal everything, including Conda
RUN apt-get update && apt-get install libdb-dev curl -y && \
	curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && \
	bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b && \
	rm Miniconda-latest-Linux-x86_64.sh && \
	curl -LO http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz && \
	tar -xvf mga_x86_64.tar.gz -C /usr/local/bin/ mga_linux_ia64 && \
	conda update -y conda && \
	conda install -y -c bioconda mcl=14.137 muscle blast perl-bioperl perl-file-which hmmer=3.1b2 perl-parallel-forkmanager perl-list-moreutils diamond && \
	conda clean --yes --tarballs --packages --source-cache && \
	apt-get purge -y --auto-remove curl ca-certificates && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENTRYPOINT ["wrapper_phage_contigs_sorter_iPlant.pl"]

CMD ["-h"]
