# README #

## RUCS - *R*apid Identification of PCR Primers Pairs for *U*nique *C*ore *S*equences ##
This repository contains the source code for a bioinformatics tool, which have
several usages:

1. Find sequences which are unique to a dataset of positive samples compared to a dataset of negative samples
2. Identify PCR primer pairs for a given set of sequences
3. Run PCR in silico, for a given set of primer pairs against a given set of references
4. Annotate a given set of sequences for protein annotations from the NCBI refseq database
5. Show PCR statistics for a given primer set to a given template
6. Combine all of the above functionalities in a pipeline to provide a tool for rapid identification of PCR primer Pairs for the unique target sequences of a positive dataset versus a negative dataset

Authors:
   Martin Christen Frølund Thomsen,
   Henrik Hasman,
   Ole Lund


## Entry Point Descriptions ##
The entrypoints provide a quick interface to running the services in this
software package.

To modify the settings and parameters for any entry point, modify or make a
modified copy of settings.default.cjson.

If you want more flexibility, the primer_core_tools.py module can be imported
into a python script, then you will have full access to all features. The
classes and functions should be well documented in the source code.

### full - Find PCR Primer Pairs for the Unique Core Sequences ###
This is the main method, which combines fucs and fppp into one serial execution.

Example of usage:
```
#!bash
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb \
       rucs full --positives positives/* other/positive.fa --negatives negatives/*
```
Notice how you can specify multiple paths with or without wild cards as input.
This is true for --positives and --negatives, and for vpcr's --references.

Note: To run the tool outside of docker, install the dependencies described in
the Dockerfile on your local machine, and make a symbolic link to
primer_core_tools.py in your /usr/local/bin called rucs. That way you can run
the tool directly using the second line of the command above.

Note: You can specify which entries in the template you want analysed. Default
for full is running against the first entry.

### fucs - Find Unique Core Sequences ###
This method finds all the core sequences from the positive dataset, remove
any sequence found in the negative dataset and returns two fasta files:
contigs, containing the unique core sequences; dissected scaffolds,
containing the fragments of the scaffolds which are usable for primer design.

Example of usage:
```
#!bash
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb \
       rucs fucs --positives positives/* --negatives negatives/*
```

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
```
#!bash
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb \
       rucs fppp --template template.fa --positives positives/* --negatives negatives/*
```

### vpcr - Virtual PCR ###
Simulate PCR in silico for a list of primer pairs against a list of references

Example of usage:
```
#!bash
docker run --rm -v `pwd`:/workdir \
       rucs vpcr --pairs pair_file.tsv --references references/*
```

### anno - Annotate Sequences ###
Annotate provided sequence with BLAST protein annotations.
A BLAST protein database must be installed. we recommend the swissprot database.
It is reccommended to use the environment variable $BLASTDB to point to the
BLAST repository, but a full path to the directory can also be provided.

Example of usage:
```
#!bash
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb \
       rucs anno --template template.fa
```   

### pcrs - Show PCR statistics ###
This method will annotate a PCR primer set with PCR statistics, such as
primer Tm, Hairpin Tm, primer-probe distance and much more.

Example of usage:
```
#!bash
docker run --rm -v `pwd`:/workdir \
       rucs spst --pairs pair_file.tsv --template template.fa
```

## Installation and setup? ##
1. Install and start docker
2. Install the RUCS program
3. Download BLAST annotation database

### Install the RUCS program ###
There are several ways to install the RUCS programs



#### Pull image from DockerHub (option 1) - The easiest ####
https://hub.docker.com/r/genomicepidemiology/rucs/

1. Pull image
2. Tag image
3. Run test and see if everything is ok and ready to use

Commands for installation:
```
#!bash
docker image pull genomicepidemiology/rucs
docker tag genomicepidemiology/rucs rucs
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb rucs test
```

#### Clone Git repository and build image (option 2) - the more challenging ####
1. Clone this repository
2. Build docker image
3. Run test and see if everything is ok and ready to use

Commands for installation:
```
#!bash
git clone https://github.com/martinfthomsen/rucs2
cd rucs2
docker-compose build
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb rucs test
```

#### Direct install (option 3) - the most challenging ####
**Not using docker?** You can check out what to install and configure on your
system from the Dockerfile.
After everything is installed and ready, Verify that all dependencies are in
your local PATH and the Python modules are properly installed:
```
#!bash
which python3 samtools bwa blastn blastx makeblastdb
python3 -c 'import gzip, json, types, shutil, glob, bisect, primer3, numpy, subprocess, difflib, tabulate'
```


### Download BLAST annotation database ###
A BLAST protein database must be installed separately if you want RUCS to
annotate the results with protein annotations. we recommend installing the
swissprot database. It is recommended to use the environment variable $BLASTDB
to point to the BLAST database.
```
#!bash
BLASTDB=/blastdb
mkdir $BLASTDB
docker run --rm --entrypoint update_blastdb.pl -v $BLASTDB:/workdir rucs --decompress swissprot
```

If you are using the the $BLASTDB environment variable, it is recommended to set
it in your bash_profile file:
```
#!bash
echo "BLASTDB=$BLASTDB" >> ~/.bash_profile
```

Alternatively, it is possible to run the local script to this repo
(install_db.sh), to download multiple databases in parallel:
```
#!bash
./install_db.sh $BLASTDB swissprot
 ```
The alternative script also downloads the file all_blast_db_files.txt, which
lists all available blast databases.
and a ...-metadata.json file with details of the database

You can also choose to download other databases such as refseq_protein.
Downloading the refseq_protein database can take a while, since the database
is > 20 GB...

**OBS: If you don't need the annotation capabilities, you can opt out of this!**


## Usage example? ##


### See help menu ###
To see the help page, run the following command:
```
#!bash
docker run --rm rucs --help
 ```

### Example of full run using default options ###
This example assumes that the $BLASTDB environment variable is set. Here are the
steps:
1. Create run directory
2. Setup positive dataset
   - download/copy the positive sequences
3. Setup negative dataset
   - download/copy the negative sequences
4. Run the full RUCS command
5. Inspect the results, etc.


```
#!bash
mkdir ~/my_first_rucs_analysis && cd $_
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb \
       rucs full --positives positives/* other/positive.fa --negatives negatives/*

 ```

### Example of running in interactive mode ###
This example assumes the same setup steps as the example of the full run
(step 1-3) has been performed.
Steps:
4. Run the docker image in interactive mode
5. Run the full analysis command directly on the main script
6. Inspect the results, etc.

```
#!bash
docker run -it --entrypoint /bin/bash --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb rucs
primer_core_tools.py full --positives positives/* other/positive.fa --negatives negatives/*

```
In the interactive mode, you are inside the RUCS container, and you will only be
able to interact with the container and the material you bring with you through
mounting, like the workdir and blastdb in the command above.
You can run any command available in there and manipulate what ever files you
wish. OBS! Only changes made in the mounted directories are kept after you exit
the container.


## Result Explanation ##

#### full run ####
The full run creates 8 note worthy files.
**debug.log** This file gives an overview and statistics on the progress.

The fucs part of the algorithm produces 3 important result files:
* core_sequences.contigs.fa
* unique_core_sequences.contigs.fa
* unique_core_sequences.disscafs.fa
These files contain the core sequences and unique core sequences respectively in
two formats. Contigs contains the clean sequences annotated with position in the
reference. disscafs contains the dissected scaffolds, where the sequences which
are close in proximity have been joined with stretches of n's to fill out the
gaps between the sequences.

The fppp part of the algorithm creates 3 important files:
**results_best.tsv** This file shows the best matches
**results.tsv** This shows the details for all the tested pairs
**products.tsv** This file shows the vpcr results for all the tested pairs

#### fucs run ####
See first part of the full run

#### fppp run ####
See last part of the full run

#### vpcr run ####
This entry point provides one important file
**products.tsv** This file shows the vpcr results for all the tested pairs

#### anno run ####
This entry point provides no result files, but instead shows the annotations
found directly on the screen in a json format.
EG.
```
#!json
{'0': [(1, 735, ['DIM/GIM/SIM family subclass B1 metallo-beta-lactamase'])]}
```
The key **'0'** is the contig name where the annotaion was discovered.

The value is a list of hits found on different places in the contig.

The first item in a hit **1** is the starting position of the hit.

The second item in a hit **735** is the end position of the hit.

The third item in a hit
**['DIM/GIM/SIM family subclass B1 metallo-beta-lactamase']** is a reduced list
of BLAST annotations matching the given hit.

#### pcrs run ####
This entry point provides no result files, but instead shows the the statistics
for the provided pairs directly on the screen.


## Troubleshoot ##
#### MAC OS: Cannot run build command ####
Try running this command to resolve any environment issues:
```
#!bash
eval "$(docker-machine env default)"
```

#### MAC OS: Could not resolve 'http.debian.net' ####
To solve this, you need to set the DNS of your docker machine.

1. Edit ~/.docker/machine/machines/default/config.json
2. Locate "Dns" under "HostOptions" and "EngineOptions"
3. Add "8.8.8.8" to the list
```
#!markup
    "Dns": ["8.8.8.8"],
```

Then restart the machine
```
#!bash
docker-machine restart default
eval "$(docker-machine env default)"
```

#### MAC OS: Docker does not mount my directory! ####
Check if you have added the directory to the file sharing in the docker
preferences. If not add it, and if you have an old version of docker, consider
updating to a new docker version.

#### MAC OS: BLAST annotation fails without errors! ####
This is probably caused by lack of RAM. For MACs and Windows, a virtual machine
is used to run a Linux environment where Docker can run. This machine has a
limit on how much of the host's RAM it may access. Check if you have enough RAM
allocated, or try to increase the amount of RAM for your machine. Otherwise,
consider using a smaller BLAST database.

The answer to "MAC OS: Docker takes up too much space?" presents a way to
remove the default machine and create a new with more CPU, RAM and disk space.

#### MAC OS: Docker takes up too much space? ####
Check the disk space use on the machine:
1. ssh into the virtual machine
```
#!bash
docker-machine ssh default
```
2. Check the disk usage
```
#!bash
df -h
```

If you find that the machine is using too much space and you are not worried
about losing the data on the machine, you can delete the machine and recreate it.

Stop and delete the default machine. (WARNING, this will remove all your
containers, and all data stored inside the containers is lost!)
```
#!bash
docker-machine stop default
docker-machine rm -f default
docker-machine create -d virtualbox --virtualbox-cpu-count=4 \
                                    --virtualbox-memory=12288 \
                                    --virtualbox-disk-size=15360 \
                                    default
```
--virtualbox-cpu-count sets the number of CPUs allowed to be used (4 CPUs in this example)

--virtualbox-memory sets the amount of RAM allowed to be used (12GB RAM in this example)

--virtualbox-disk-size sets the disk size. (15GB in this example)

Now you can check the diskspace again, and it should be all good again.
Next step is now to reinstall all your images...


## Docker Cleanup Commands ##
```
#!bash
# Stop and remove all containers (instances of images)
docker rm $(docker stop $(docker ps -aq))
# Remove all exited containers
docker rm -v $(docker ps -aq -f status=exited)
# Remove all dangling images
docker rmi $(docker images -qf "dangling=true")
# Remove all dangling volumes
docker volume rm $(docker volume ls -qf dangling=true)
```

## Who do I talk to? ##

* Repo owner or admin


License
=======

Copyright (c) 2017, Martin Christen Frølund Thomsen, Technical University of Denmark.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
