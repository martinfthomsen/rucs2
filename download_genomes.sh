#!/bin/bash

# This script downloades genomes through Entrez Direct cmd-line tools based on a
# list of accessions through stdin. fx
# cat acc_list.txt | download_genomes.sh

# or as a string of space-separated arguments
# download_genomes.sh CP000672.1 ARBW00000000 NZ_KU341381.1
# download_genomes.sh CP000672.1 ARBW00000000 NZ_KU341381.1

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
    # Try download the accessions from the assembly db
    found=false
    echo "Checking against the assembly DB..."
    esearch -db assembly -query "$accession" | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r line ;
        do
            fname=$(echo $line | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
            found=true
            if test -f $accession"_"$fname || test -f $accession"_"${fname%.*}; then
                echo "found $fname! File previously downloaded! skipping..."
            else
                echo "found $fname! Downloading..."
                wget -q --output-document $accession"_"$fname "$line/$fname";
            fi
        done
    if [ "$found" = false ] ; then
        echo "Checking against the nuccore DB..."
        # Try download the accessions from the nuccore db
        AssemblyAcc=$(esearch -db nuccore -query "$accession" | esummary | xtract -pattern DocumentSummary -element AssemblyAcc)
        if [ ! -z "$AssemblyAcc" ] ;then
            echo "found $AssemblyAcc! Downloading..."
            efetch -db nuccore -id $AssemblyAcc -format fasta | gzip -c > $accession.fna.gz
        fi
    fi

done
