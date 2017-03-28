############################################################
# Dockerfile to build Primer Finder image
############################################################

# Set base image to Python Anaconda
FROM continuumio/anaconda3

# File Author / Maintainer
MAINTAINER Martin Christen Frølund Thomsen

# Install apt-get essentials
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get install -y --no-install-recommends apt-utils debian-archive-keyring
RUN apt-get update && apt-key update && apt-get install -y build-essential

# Install general dependencies
RUN apt-get install -y --fix-missing --no-install-recommends \
    perl \
    git \
    emacs \
    wget \
    vim

# Make tools directory
RUN mkdir /tools

# Install Samtools dependencies
RUN apt-get install -y \
    autoconf \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev #libncursesw5-dev

# Install Tools
WORKDIR /tools

# Samtools
RUN git clone https://github.com/samtools/htslib /tools/htslib
RUN cd htslib; make
RUN git clone https://github.com/samtools/samtools /tools/samtools
RUN cd samtools; make; make install

# BWA
RUN git clone --branch v0.7.15 https://github.com/lh3/bwa /tools/bwa
RUN cd bwa; make
ENV PATH ${PATH}:/tools/bwa

# BLAST (newest blast is only available through ftp, sadly not via apt-get)
RUN mkdir ncbi-blast
RUN wget -r --no-parent -nv -A 'ncbi-blast-*+-x64-linux.tar.gz' -O ncbi-blast.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
RUN tar -xzf ncbi-blast.tar.gz -C ncbi-blast --strip-components=1 && rm ncbi-blast.tar.gz
ENV PATH $PATH:/tools/ncbi-blast/bin

# Install Python Module Dependencies
RUN pip install --upgrade pip primer3-py tabulate #pysam
ENV PATH $PATH:/opt/conda/bin/

# Set BLAST DB Environment
#ENV BLASTDB /blastdb
# NOTE: BLASTx does not work in Docker, keeps dying, best guess is due to high
# IO load, memory and CPU seems fine

# Add This repo to the image
RUN mkdir /repo
ADD . /repo

# Install Entry Points
ENV PATH $PATH:/repo

# Set convenience aliases
RUN echo "alias edit='emacs'" >> ~/.bashrc
RUN echo "alias ls='ls -h --color=tty'" >> ~/.bashrc
RUN echo "alias ll='ls -lrt'" >> ~/.bashrc
RUN echo "alias l='less'" >> ~/.bashrc
RUN echo "alias du='du -hP --max-depth=1'" >> ~/.bashrc
RUN echo "alias cwd='readlink -f .'" >> ~/.bashrc

# Set Entry Point
ENTRYPOINT ["primer_core_tools.py"]
CMD ["--help"]

WORKDIR /workdir
