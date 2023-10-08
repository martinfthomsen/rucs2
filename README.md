# README #
This repo was forked from the original RUCS repo (2023): https://bitbucket.org/genomicepidemiology/rucs

This tool was created as part of a research project.

For publication of results made using this tool, please cite:
Martin Christen Frølund Thomsen, Henrik Hasman, Henrik Westh, Hülya Kaya, Ole Lund, RUCS: rapid identification of PCR primers for unique core sequences, Bioinformatics, Volume 33, Issue 24, December 2017, Pages 3917–3921, https://doi.org/10.1093/bioinformatics/btx526

A freely available online implementation of the original RUCS can be found here: https://cge.food.dtu.dk/services/RUCS/

## RUCS - *R*apid Identification of PCR Primers Pairs for *U*nique *C*ore *S*equences ##
This repository contains the source code for a bioinformatics tool, which have
several usage cases:

1. Find sequences which are unique to a dataset of positive samples compared to a dataset of negative samples
2. Identify PCR primer pairs for a given set of sequences
3. Run PCR in silico, for a given set of primer pairs against a given set of references
4. Annotate a given set of sequences for protein annotations from the NCBI refseq database
5. Show PCR statistics for a given primer set to a given template
6. Combine all of the above functionalities in a pipeline to provide a tool for rapid identification of PCR primer Pairs for the unique target sequences of a positive dataset versus a negative dataset

Authors:
   Martin Christen Frølund Thomsen

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
```bash
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb \
       rucs full -v --positives positives/* other/positive.fa --negatives negatives/*
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
```bash
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
```bash
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb \
       rucs fppp --template template.fa --positives positives/* --negatives negatives/*
```

### vpcr - Virtual PCR ###
Simulate PCR in silico for a list of primer pairs against a list of references

Example of usage:
```bash
docker run --rm -v `pwd`:/workdir \
       rucs vpcr --pairs pair_file.tsv --references references/*
```

### anno - Annotate Sequences ###
Annotate provided sequence with BLAST protein annotations.
A BLAST protein database must be installed. we recommend the swissprot database.
It is reccommended to use the environment variable $BLASTDB to point to the
BLAST repository, but a full path to the directory can also be provided.

Example of usage:
```bash
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb \
       rucs anno --template template.fa
```   

### pcrs - Show PCR statistics ###
This method will annotate a PCR primer set with PCR statistics, such as
primer Tm, Hairpin Tm, primer-probe distance and much more.

Example of usage:
```bash
docker run --rm -v `pwd`:/workdir \
       rucs pcrs --pairs pair_file.tsv --template template.fa
```

### expl - Explore the genomes ###
Explore the genomes of the positive dataset versus the negative dataset for
Over- and underrepresentation of k-mers.

This method computes the over- and underrepresentation of k-mers in the
positive genomes versus the negative genomes, and generates consensus sequences
from these.
It can be useful if you have two cluster of the same species that are pertain
different effects and you aren't sure of why. In such a case, this method should
provide you with results that can tell you what the differences are genetically.

**Output**
The method creates two output files in fasta format.
The first, containing the overrepresented sequences. The second, containing the
underrepresented sequences. An overrepresented sequence is a sequence of DNA
found significantly more often in the positive dataset compared to the negative
dataset. An underrepresented sequence is a sequence of DNA found significantly
more often in the negative dataset compared to the positive dataset.
The fasta description header will contain the genome names where the given sequence is found.

**Options**
There are four specific settings for this feature:
 - kmer_count_threshold [default 1] - kmers found below this limit within each file are ignored
 - sensitivity_threshold [default 0.6] - The sensitivity threshold defines how often a kmer must be found in the positive set to be included in the results. If provided as an integer the number is considered as a minimum count.
 - fall_out_threshold [default 0.2] - The fall-out threshold defines how often k-mers may be found in the negative set and still be included in the results. If provided as an integer the number is considered as a maximum count.
 - align_percent_threshold [default 0.05] - The alignment percent threshold defines the acceptable amount of kmers to not be aligned to a contig. These k-mers are lost from further analysis to speed up the process. Set to 0, if you want as much data as possible.

**Example of usage**
```bash
docker run --rm -v `pwd`:/workdir rucs expl -v --positives positives/* --negatives negatives/*
```

**Fails?**
If you experience this method failing due to being killed, you can try
increasing the allowed memory usage in the Docker setting. This method is quite
memory intensive, a dataset of ~30 bacteria genomes can easily use more than
10GB RAM (peak).

**Slow?**
This method is very cumbersome, which many heavy compute tasks, so the more
references are used and the bigger the references, this can really take a long
time to run. easpecially the step where the over- and underrepresented k-mers
are converted into contigs and scaffolds sequences.


## Other hidden features ##
The RUCS image comes with some hidden but convenient features/functionalities
that can make your experience with the tool easier.

### Downloading Genomes through Entrez Direct ###
The RUCS image contains the Entrez Direct cmd-line tools and a small wrapper to
make it easier for the user to download genomes when you know an accession ID
like CP000672.1 or ARBW00000000.

To use the download script, the entrypoint must be changed, and a directory
where the downloaded files should be stored must be mounted to the internal
/workdir. And last but not least, the accession IDs must be provided as
arguments.

The download script can be invoked like so:
```bash
docker run -it --entrypoint download_genomes.sh --rm -v `pwd`/positives:/workdir rucs CP063056 ARBW00000000 BBXJ00000000
docker run -it --entrypoint download_genomes.sh --rm -v `pwd`/negatives:/workdir rucs JWIZ01 JADGLC01 CP000672.1
```
Just switch out "`pwd`/positives" with your preferred download directory, and
"CP000672.1 ARBW00000000", with your list of space-separated accessions.

** Fails? **
If the download script fails, it could be due to bad internet connection to NCBI.
Try opening a browser and go to:
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=JWIZ01
If this site is able to load, the download script should also work.


## Installation and setup? ##
1. Install and start docker
2. Install the RUCS program
3. Download BLAST annotation database

### Install the RUCS program ###
There are several ways to install the RUCS programs

#### Pull latest image from DockerHub (option 1) - The easiest ####
https://hub.docker.com/repository/docker/mcft/rucs2

1. Pull image
2. Tag image
3. Run test and see if everything is ok and ready to use

Commands for installation:
```bash
docker image pull mcft/rucs2
docker tag mcft/rucs2 rucs
docker run --rm -v `pwd`/test:/workdir -v $BLASTDB:/blastdb rucs test
```

##### Pull the original version #####
```bash
docker image pull genomicepidemiology/rucs
docker tag genomicepidemiology/rucs rucs
docker run --rm -v `pwd`/test:/workdir -v $BLASTDB:/blastdb rucs test
```

#### Clone Git repository and build image (option 2) - the more challenging ####
1. Clone this repository
2. Build docker image
3. Run test and see if everything is ok and ready to use

Commands for installation:
```bash
git clone https://github.com/martinfthomsen/rucs2
cd rucs2
docker-compose build
docker run --rm -v `pwd`/test:/workdir -v $BLASTDB:/blastdb rucs test
```

#### Direct install (option 3) - the most challenging ####
**Not using docker?** You can check out what to install and configure on your
system from the Dockerfile.
After everything is installed and ready, Verify that all dependencies are in
your local PATH and the Python modules are properly installed:
```bash
which python3 samtools bwa blastn blastx makeblastdb
python3 -c 'import gzip, json, types, shutil, glob, bisect, primer3, numpy, subprocess, difflib, tabulate'
```


### Download BLAST annotation database ###
A BLAST protein database must be installed separately if you want RUCS to
annotate the results with protein annotations. we recommend installing the
swissprot database. It is recommended to use the environment variable $BLASTDB
to point to the BLAST database.
```bash
BLASTDB=/blastdb
mkdir $BLASTDB
docker run --rm --entrypoint update_blastdb.pl -v $BLASTDB:/workdir rucs --decompress swissprot
```

If you are using the the $BLASTDB environment variable, it is recommended to set
it in your bash_profile file:
```bash
echo "BLASTDB=$BLASTDB" >> ~/.bash_profile
```

Alternatively, it is possible to run the local script to this repo
(install_db.sh), to download multiple databases in parallel:
```bash
./install_db.sh $BLASTDB swissprot
 ```
The alternative script also downloads the file all_blast_db_files.txt, which
lists all available blast databases.
and a ...-metadata.json file with details of the database

You can also choose to download other databases such as refseq_protein.
Downloading the refseq_protein database can take a while, since the database
is > 20 GB...

***OBS: If you don't need the annotation capabilities, you can opt out of this!***


## Usage example? ##


### See help menu ###
To see the help page, run the following command:
```bash
docker run --rm rucs --help
 ```

### Example of full run using default options ###
This example below assumes that the $BLASTDB environment variable is set.

*Command:*
```bash
docker run --rm -v `pwd`/test:/workdir -v $BLASTDB:/blastdb rucs full --positives CP000672.1 ARBW00000000 --negatives JWIZ01
```
RUCS will create the run directory "test" in your current directory, download
the fasta files for the two provided accession IDs in the positives argument,
and the single accession ID in the negatives argument. RUCS will then run the
full analysis algorithm to identify the good PCR primer pairs that only produces
amplicons for the samples in the positives argument.


As the -v option is provided to RUCS, the screen should display the following
progress log information:
```bash
processing CP000672.1...
found GCA_000016485.1_ASM1648v1_genomic.fna.gz! Downloading...
processing ARBW00000000...
found GCA_000379905.1_ASM37990v1_genomic.fna.gz! Downloading...
processing JWIZ01...
found GCA_001038205.1_ASM103820v1_genomic.fna.gz! Downloading...
# Running RUCS...
#    Prepare reference: CP000672.1_GCA_000016485.1_ASM1648v1_genomic.fna.gz...
#    Prepare inputs: 2 positive and 1 negative genomes...
#    Finding unique core sequences...
#       Computing k-mer intersection...
#          Extracting k-mers from CP000672.1_GCA_000016485.1_ASM1648v1_genomic.fna.gz...
#          Extracting k-mers from ARBW00000000_GCA_000379905.1_ASM37990v1_genomic.fna.gz...
#       Aligning core_kmers.fq to reference.fa...
#       Computing k-mer contigs and scaffolds...
#       Computing k-mer complement to the negative genomes...
#          Extracting k-mers from JWIZ01_GCA_001038205.1_ASM103820v1_genomic.fna.gz...
#          Filtering the negative k-mers found in {os.path.basename(genome)}...
#       Computing k-mer contigs and scaffolds...
#    Find valid primer pairs for PCR...
#       Scanning contig 0_1375921_0 for primer pairs...
#       Aligning primers to positive references...
#          Aligning primers.fa to CP000672.1_GCA_000016485.1_ASM1648v1_genomic.fna.gz...
#          Aligning primers.fa to ARBW00000000_GCA_000379905.1_ASM37990v1_genomic.fna.gz...
#       Aligning primers to negative references...
#          Aligning primers.fa to JWIZ01_GCA_001038205.1_ASM103820v1_genomic.fna.gz...
#       Ranking primer pairs...
#       Sorting primer pairs...
#       Validating primer pairs...
#       Current good pp: 71...
#       Annotating PCR product environment...
372 sequence were ignored since they were not in the sequence selection!
```
Your new run directory should now hold the following files:

```bash
ARBW00000000_GCA_000379905.1_ASM37990v1_genomic.fna.gz
CP000672.1_GCA_000016485.1_ASM1648v1_genomic.fna.gz
JWIZ01_GCA_001038205.1_ASM103820v1_genomic.fna.gz
core_sequences.aux.tsv
core_sequences.contigs.fa
core_sequences.disscafs.fa
core_sequences.scaffolds.fa
pairs.json
pcr_products_with_skirts.blastx.err
pcr_products_with_skirts.blastx.tsv
pcr_products_with_skirts.fa
primers.fa
products.tsv
results.tsv
results_best.tsv
stats.log
unique_core_sequences.aux.tsv
unique_core_sequences.contigs.fa
unique_core_sequences.disscafs.fa
unique_core_sequences.scaffolds.fa
```

You will find a list of the best candidates in the "results_best.tsv" file.


### Example of running in interactive mode ###
This example assumes the same setup steps as the example of the full run
(step 1-3) has been performed.

Steps:

1. Run the docker image in interactive mode
2. Run the full analysis command directly on the main script
3. Inspect the results, etc.

```bash
docker run -it --entrypoint /bin/bash --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb rucs
primer_core_tools.py full --positives positives/* other/positive.fa --negatives negatives/*
```
In the interactive mode, you are inside the RUCS container, and you will only be
able to interact with the container and the material you bring with you through
mounting, like the workdir and blastdb in the command above.
You can run any command available in there and manipulate what ever files you
wish. OBS! Only changes made in the mounted directories are kept after you exit
the container.

#### Open interactive python terminal inside the container ####
You can also run python inside the container and import the primer core tools
from RUCS.

The following command will open an interactive terminal inside the RUCS container.
```bash
docker run -it --entrypoint python3 --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb rucs
```

#### Run full RUCS through interactive python terminal ####
The following commands will run the full RUCS algorithm using the files in the
positives directory as the positives references, and vice-versa for the
negatives references.
```python
>>> import sys
>>> sys.path.append('/tools/')
>>> from primer_core_tools import *
>>> load_global_settings('settings.default.cjson')
>>> positives = [x for x in glob.glob("positives/*") if check_file_type(x) == 'fasta']
>>> negatives = [x for x in glob.glob("negatives/*") if check_file_type(x) == 'fasta']
>>> reference = positives[0]
>>> main(positives, negatives, reference, quiet=False, clean_run=True, annotate=True)
...
```

#### Look for available files ####
To check which files are available to you in your directories you can run the
following command.

This example finds all files ending with ".fa" in your current working directory
```python
>>> find_files('*.fa')
['core_sequences.contigs.fa', 'core_sequences.disscafs.fa', ...]
```

#### Load fasta file ####
Sometimes when working with RUCS, it is nice to be able to interact with the
sequence data from the involved fasta files.

Below is an example of how to parse a fasta file and store the sequences in a list
```python
>>> my_seqs = [seq for seq, name, desc in seqs_from_file('bla.fa')]
>>> my_seqs[0]
'ATGCGACCATTTTT...'
```

#### Load results_best.tsv ####
The primary results from RUCS are stored in the results_best.tsv file. These are
the good primer pair candidates identified to spefically produce amplicons for
the positives references. The results are stored as a TSV file which can be
loadet and processed in the terminal.

Below is an example of how to load the identified primer pairs and statistics
from results_best.tsv and extract useful information.
```python
>>> results = parse_tsv('results_best.tsv')
>>> results[0]['sequence_id']
'0_1375921_0'
>>> results[0]['forward_primer']
'CACCCAGTAGAGCACACTTTG'
>>> results[0]['forward_position']
'2525'
```

#### Extract subsequence from sequence ####
In the results_best.tsv file you will find the primer/probe positions, but as
the positions here refer to the position in the dissected scaffold of the unique
core sequence, it rarely translates into the position in the reference sequence.
However, the results also provides the sequence id from the dissected scaffolds,
which provides the necessary information to compute the primer/probe position
relative to the reference sequence.

The first number of the sequence id refer to the index of the sequence in the
reference.
The second number of the sequence id refer to the position of the dissected
scaffold in the reference sequence.

Below is an example of how to use the sequence id and the forward_position to
extract the forward_primer sequence from the reference sequence.
```python
>>> sequence_id = "0_1375921_0"
>>> forward_primer = "CACCCAGTAGAGCACACTTTG"
>>> forward_position = "2525"
>>> sid = int(sequence_id.split('_')[0])
>>> spos = int(sequence_id.split('_')[1])
>>> ppos = int(forward_position)
>>> length = len(forward_primer)
>>> my_seqs[sid][spos+ppos:spos+ppos+length] == forward_primer
True
```


## Result Explanation ##

### full run ###
The full run creates 8 note worthy files.
**stats.log** This file gives an overview and statistics on the progress.

The fucs part of the algorithm produces 3 important result files:
* **core_sequences.contigs.fa** This file contain the core sequences, the sequences which are common to all the positive references .
* **unique_core_sequences.contigs.fa** This file contain the unique core sequences, the sequences which are common to all the positive references and not found in any of the negative references.
* **unique_core_sequences.disscafs.fa** This file contain the dissected scaffolds for the unique core sequences.

*contigs* contains the clean sequences annotated with position in the reference.

*disscafs* contains the dissected scaffolds, where the sequences which
are close in proximity have been joined with stretches of n's to fill out the
gaps between the sequences.
This makes it possible to identify primer pairs which are uniquely binding but where the amplicon contains non-unique stretches of DNA.

The fppp part of the algorithm creates 3 important files:
* **results_best.tsv** This file shows the best matches for primer pair candidates.
* **results.tsv** This file provides an exaustive list of all the tested primer pair candidates.
* **products.tsv** This file provides the results of the virtual PCR for all tested primer pair candidates.

### fucs run ###
See first part of the full run

### fppp run ###
See last part of the full run

### vpcr run ###
This entry point provides one important file
**products.tsv** This file shows the vpcr results for all the tested pairs

### anno run ###
This entry point provides no result files, but instead shows the annotations
found directly on the screen in a json format.
EG.
```json
{"0": [1, 735, ["DIM/GIM/SIM family subclass B1 metallo-beta-lactamase"]]}
```
The key **'0'** is the contig name where the annotaion was discovered.

The value is a list of hits found on different places in the contig.

The first item in a hit **1** is the starting position of the hit.

The second item in a hit **735** is the end position of the hit.

The third item in a hit
**['DIM/GIM/SIM family subclass B1 metallo-beta-lactamase']** is a reduced list
of BLAST annotations matching the given hit.

### pcrs run ###
This entry point provides no result files, but instead shows the the statistics
for the provided pairs directly on the screen.

### expl run ###
The Explore entrypoint is similar to the fucs entry point. But where fucs method is
focussed on identifying unique core sequences, the Explore method tries to find
over-represented (ors) and under-represented sequences (urs):
* **ors.contigs.fa** This file contain the contigs of the over-represented sequences.
* **ors.disscafs.fa** This file contain the dissected scaffolds of the over-represented sequences.
* **urs.contigs.fa** This file contain the contigs of the under-represented sequences.
* **urs.disscafs.fa** This file contain the dissected scaffolds of the under-represented sequences.

ORS can be used to identify primer pairs to target most of the positive references, and
URS can be used to identify primer pairs to target most of the negative references.

For instance, sometimes what strains have in common might not be what DNA they share, but ratcher what they are missing. The URS will provide this insight of any missing DNA.
The explore method stops at this stage. To get primer pair identification, Virtual PCR, etc. please run the FPPP method using either the ors.disscafs.fa or urs.disscafs.fa as the template option.


## Troubleshoot ##
### MAC OS: BLAST annotation fails without errors! ###
This is probably caused by lack of RAM. For MACs and Windows, a virtual machine
is used to run a Linux environment where Docker can run. This machine has a
limit on how much of the host's RAM it may access. Check if you have enough RAM
allocated, or try to increase the amount of RAM for your machine. Otherwise,
consider using a smaller BLAST database.

## Docker Cleanup Commands ##
```bash
# Stop and remove all containers (instances of images)
docker rm $(docker stop $(docker ps -aq))
# Remove all exited containers
docker rm -v $(docker ps -aq -f status=exited)
# Remove all dangling images
docker rmi $(docker images -qf "dangling=true")
# Remove all dangling volumes
docker volume rm $(docker volume ls -qf dangling=true)
```

## Contact ##
You can report your issues here: https://github.com/martinfthomsen/rucs2/issues


## License ##

**RUCS 2**: https://github.com/martinfthomsen/rucs2

Copyright ©️ 2023, Martin Christen Frølund Thomsen.


**Original RUCS**: https://bitbucket.org/genomicepidemiology/rucs

Copyright ©️ 2017, Martin Christen Frølund Thomsen, Technical University of Denmark.


*All rights reserved.*


Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0


Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

See the License for the specific language governing permissions and
limitations under the License.
