############################################################
# Dockerfile to build image
############################################################

# Load base Docker image
FROM python:3.5.3-jessie

# Disable Debian frontend interaction promts
ENV DEBIAN_FRONTEND noninteractive

# Install dependencies with apt-get
RUN set -ex; \
    apt-get update -qq; \
    apt-get install -y -qq apt-utils; \
    apt-get install -y -qq \
        perl \
        git \
        wget \
        autoconf \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libz-dev \
    ; \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*;

# Re-enable Debian frontend interaction promts
ENV DEBIAN_FRONTEND Teletype

# Install Python dependencies
RUN pip install --upgrade pip cython tabulate primer3-py numpy

# Install tools
RUN mkdir -p /tools/test
ENV PATH $PATH:/tools/
WORKDIR /tools

# Copy repository files to /tools/
COPY ./primer_core_tools.py /tools/
COPY ./settings.default.cjson /tools/
COPY ./test/* /tools/test/

# Install BLAST (newest blast is only available through ftp, sadly not via apt-get)
ENV BLASTDB /blastdb
ENV PATH $PATH:/tools/ncbi-blast/bin
RUN mkdir ncbi-blast && \
    wget -r --no-parent -nv -A 'ncbi-blast-*+-x64-linux.tar.gz' -O ncbi-blast.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ && \
    tar -xzf ncbi-blast.tar.gz -C ncbi-blast --strip-components=1 && \
    rm ncbi-blast.tar.gz

# Install BWA
ENV PATH ${PATH}:/tools/bwa
RUN git clone --branch v0.7.15 https://github.com/lh3/bwa /tools/bwa && \
    cd bwa && \
    make

# Install Samtools
ENV HTSDIR /tools/htslib
RUN git clone https://github.com/samtools/htslib /tools/htslib && \
    git clone https://github.com/samtools/samtools /tools/samtools && \
    cd /tools/samtools && \
    make && \
    make install

# Set convenience aliases
RUN echo "alias ls='ls -h --color=tty'" >> ~/.bashrc && \
    echo "alias ll='ls -lrt'" >> ~/.bashrc && \
    echo "alias l='less'" >> ~/.bashrc && \
    echo "alias du='du -hP --max-depth=1'" >> ~/.bashrc && \
    echo "alias cwd='readlink -f .'" >> ~/.bashrc;

# Set entry point (This makes the container invoke the tool directly)
ENTRYPOINT ["primer_core_tools.py"]
CMD ["--help"]

WORKDIR /workdir
