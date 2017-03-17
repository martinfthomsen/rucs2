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
   Martin Christen Frølund Thomsen
   Henrik Hasman
   Ole Lund


## Entry Point Descriptions ##

### fucs - Find Unique Core Sequences ###
This method finds all the core sequences from the positive dataset, remove
any sequence found in the negative dataset and returns two fasta files:
contigs, containing the unique core sequences; dissected scaffolds,
containing the fragments of the scaffolds which are usable for primer design.

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

### vpcr - Virtual PCR ###
# Simulate PCR in silico for a list of primer pairs against a list of references

### anno - Annotate Sequences ###
Annotate provided sequence with BLAST refseq gene annotations

### spst - Show PCR statistics ###
This method will annotate a PCR primer set with PCR statistics, such as
primer Tm, Hairpin Tm, primer-probe distance and much more.


## How do I get set up? ##

1. Install docker
2. Build docker image
3. Run tests
4. Ready to use


## Who do I talk to? ##

* Repo owner or admin
