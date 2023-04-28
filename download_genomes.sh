#!/bin/bash

# This script downloades genomes through Entrez Direct cmd-line tools based on a
# list of accessions through stdin. fx
# cat acc_list.txt | download_genomes.sh

# or as a string of space-separated arguments
# download_genomes.sh CP000672.1 ARBW00000000

if [ $# -eq 0 ]; then
    if ! tty -s; then
        acc_list=$(cat - | tr '\n' ' ')
    else
        echo "Error: no input provided" >&2
        exit 1
    fi
else
    acc_list=$@
fi

for accession in $acc_list; do
    echo "processing $accession..."
    esearch -db assembly -query "$accession" | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r line ;
        do
            fname=$(echo $line | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
            echo "found $fname..."
            wget -q --output-document $accession"_"$fname "$line/$fname";
        done
done
