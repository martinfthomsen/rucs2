############################################################
# Dockerfile to build Primer Finder image
############################################################

# Set base image to Python Anaconda
FROM continuumio/anaconda

# File Author / Maintainer
MAINTAINER Martin Christen Frølund Thomsen

# Update the repository sources list
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils\
    ncbi-blast+ \
    perl \
    git \
    emacs \
    vim

# Install Additional Python Module Dependencies
RUN pip install primer3-py tabulate pysam

# Make tools directory
RUN mkdir /tools

# Install BLAST DB
RUN mkdir /tools/blastdb
ENV BLASTDB /tools/blastdb/
RUN perl update_blastdb.pl taxdb
RUN gunzip -cd refseq_protein.tar.gz | (cd $BLASTDB; tar xvf - )

# Install Samtools
RUN git clone https://github.com/samtools/samtools /tools/samtools
RUN git clone https://github.com/samtools/htslib /tools/htslib
RUN cd /tools/samtools; ./configure; make; make install

# Install BWA
RUN git clone --branch v0.7.15 https://github.com/lh3/bwa /tools/bwa
RUN cd /tools/bwa; make

# Add This repo to the image
RUN mkdir /repo
ADD . /repo

# Install Entry Points
WORKDIR /repo
RUN python setup.py install --force

# Set convenience aliases
RUN echo "alias edit='emacs'" >> /.bashrc
RUN echo "alias ls='ls -h --color=tty'" >> /.bashrc
RUN echo "alias ll='ls -lrt'" >> /.bashrc
RUN echo "alias l='less'" >> /.bashrc
RUN echo "alias du='du -hP --max-depth=1'" >> /.bashrc
RUN echo "alias cwd='readlink -f .'" >> /.bashrc
