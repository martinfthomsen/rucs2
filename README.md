# README #

## RUCS - Rapid Identification of PCR Primers Pairs for Unique Core Sequences ##
This repository contains the source code for a bioinformatics tool, which have
several usages:
1. Find sequences which are unique to a dataset of positive samples compared to a dataset of negative samples
2. Identify PCR primer pairs for a given set of sequences
3. Run PCR in silico, for a given set of primer pairs against a given set of references
4. Annotate a given set of sequences for protein annotations from the NCBI refseq database
5. Show PCR statistics for a given primer set to a given template
6. Combine all of the above functionalities in a pipeline to provide a tool for rapid identification of PCR primer Pairs for the unique target sequences of a positive dataset versus a negative dataset

Authors: 
   Martin Christen Fr�lund Thomsen,
   Henrik Hasman,
   Ole Lund


## Entry Point Descriptions ##

### full - Find PCR Primer Pairs for the Unique Core Sequences ###
This is the main method, which combines fucs and fppp into one serial execution.

Example of usage:
   docker run --rm -v `pwd`:/workdir rucs test \
          rucs full --positives positives/* --negatives negatives/*

### fucs - Find Unique Core Sequences ###
This method finds all the core sequences from the positive dataset, remove
any sequence found in the negative dataset and returns two fasta files:
contigs, containing the unique core sequences; dissected scaffolds,
containing the fragments of the scaffolds which are usable for primer design.

Example of usage:
   docker run --rm -v `pwd`:/workdir \
          rucs full --positives positives/* --negatives negatives/*

### fppp - Find PCR Primer Pairs ###
This method identifies primer pairs (and probe) using the Primer3 software.
The found pairs are then additionally tested for PCR suitability.
The eligible pairs are then clustered according to their position on the
template and sorted in the clusters according to their PCR suitability.
The best performing pairs are then BLASTed against the positive and negative
datasets and sorted according to their sensitivity, specificity, uniqueness,
noise, and PCR stats. The best candidate from each cluster gets their product
annotated with gene annotations and the list of candidates with all relevant
information is stored in a tab separated file.

Example of usage:
   docker run --rm -v `pwd`:/workdir \
          rucs fppp --template_file template.fa --positives positives/* --negatives negatives/*

### vpcr - Virtual PCR ###
Simulate PCR in silico for a list of primer pairs against a list of references

Example of usage:
   docker run --rm -v `pwd`:/workdir \
          rucs vpcr --pairs pair_file.tsv --references references/*

### anno - Annotate Sequences ###
Annotate provided sequence with BLAST refseq gene annotations

Example of usage: (No docker)
   primer_core_tools.py anno --template_file template.fa

### spst - Show PCR statistics ###
This method will annotate a PCR primer set with PCR statistics, such as
primer Tm, Hairpin Tm, primer-probe distance and much more.

Example of usage:
   docker run --rm -v `pwd`:/workdir \
          rucs spst --pairs pair_file.tsv --template_file template.fa

## How do I get set up? ##

1. Install docker (follow this guide: https://docs.docker.com/engine/getstarted/step_one/)
2. Start docker deamon
3. Clone this repository
4. Build docker image
5. (optional) Download BLAST annotation database
6. Run test and see if everything is ok
7. Ready to use

Commands for installation
   docker-machine start default # To start on MAC OS
   git clone https://bitbucket.org/genomicepidemiology/rucs.git
   cd rucs
   docker-compose build
   docker run --rm -v `pwd`:/workdir rucs test


Commands for downloading and preparing BLAST annotation DB:
   BLASTDB /blastdb
   mkdir $BLASTDB
   /usr/bin/update_blastdb --passive refseq_protein # taxdb
   gunzip -cd refseq_protein.tar.gz | (cd $BLASTDB; tar xvf - )

Downloading the refseq_protein database can take a while, since the database
is > 20 GB...
If you want to download all the databases in parallel for quicker download, feel
free to use the script install_db.sh.

OBS: If you don't need the annotation capabilities, you can opt out of this!

OBS: BLASTX does not work in docker, so to even use this feature, you need to
install the dependencies and scripts locally on your machine.

## Troubleshoot ##
### Cannot run build command on MAC OS ###
Try running this command to resolve any environment issues:
   eval "$(docker-machine env default)"

### MAC OS: Could not resolve 'http.debian.net' ###
To solve this, you need to set the DNS of your docker machine.
1. Edit ~/.docker/machine/machines/default/config.json
2. Locate "Dns" under "HostOptions" and "EngineOptions"
3. Add "8.8.8.8" to the list
   "Dns": ["8.8.8.8"],

Then restart the machine
 docker-machine restart default
 eval "$(docker-machine env default)"

### MAC OS: docker-compose build cannot download BLAST through FTP? ###
Try changing your docker preferences, from using system proxy to no proxy.

### MAC OS: Docker does not mount my directory! ###
Check if you have added the directory to the file sharing in the docker
preferences. If not add it, and if you have an old version of docker, consider
updating to a new docker version.

### MAC OS: Docker takes up too much space? ###
Check the disk space use on the machine:
1. ssh into the virtual machine
   docker-machine ssh default

2. Check the disk usage
   df -h

If you find that the machine is using too much space and you are not worried
about losing the data on the machine, you can delete the machine and recreate it.

Stop and delete the default machine. (WARNING, this will remove all your
containers, and all data is lost!)
   docker-machine stop default
   docker-machine rm default
   docker-machine create --driver virtualbox --virtualbox-disk-size "15360" default

--virtualbox-disk-size sets the disk size in MB. This example sets it to 15GB...
Now you can check the diskspace again, and it should be all good again.
Next step is now to reinstall all your images...


## Docker Cleanup Commands ##
   # Stop and remove all containers (instances of images)
   docker rm $(docker stop $(docker ps -aq))
   # Remove all exited containers
   docker rm -v $(docker ps -aq -f status=exited)
   # Remove all dangling images
   docker rmi $(docker images -qf "dangling=true")
   # Remove all dangling volumes
   docker volume rm $(docker volume ls -qf dangling=true)


## Who do I talk to? ##

* Repo owner or admin
