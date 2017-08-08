#!/bin/bash
### Install BLAST DB - Downloads and extractions run in parallel ###
# USAGE install_db.sh /blast/db swissprot
# DATABASES: swissprot refseq_protein
url='ftp://ftp.ncbi.nlm.nih.gov/blast/db/'

# OS X work-around for missing md5sum
if ! hash md5sum >/dev/null 2>&1 && hash md5 >/dev/null 2>&1; then
 # Set alias for md5sum to md5
 checkmd5sum="md5 -q"
else
 checkmd5sum="md5sum --status -c"
fi

# Set BLAST DB DIR
if [ -z "$1" ]; then
 if [ -z $BLASTDB ]; then
  echo 'No database location was specified, and the environment variable BLASTDB is not set!'
  exit 1
 fi
else
   export BLASTDB=$1
fi
# Set BLAST DB
if [ -z "$2" ]; then
 echo 'No database was specified, defaulting to refseq_protein!'
 export database='refseq_protein'
else
 export database=$2
fi

# Create BLAST directory if it does not exist
if [ ! -d $BLASTDB ]; then
 mkdir $BLASTDB
fi

# Switch to the BLAST directory
cd $BLASTDB

# Get list of files to download and install
curl $url > all_blast_db_files.txt
files=`egrep $database'.*' all_blast_db_files.txt | grep -v 'md5' | awk '{str = str" "$NF}END{print str}'`

fa=($files)
echo 'We found '${#fa[@]}' file entries!'

# Download all files
echo 'Downloading new database entries...'
for f in $files; do
 fbase=$(basename $f .tar.gz)
 curl -s -o $f'.md5' $url$f'.md5'
 # Verify that the database file does not exist before download
 if [ ! -f $fbase'.phr' ]; then
  curl -s -o $f $url$f &
  sleep 5
 else
  # Verify if the md5 sum has changed
  if [ ! -f $f'_old.md5' ] || ! cmp --silent $f'_old.md5' $f'.md5'; then
   # Remove old outdated files
   rm -f $fbase*.p*
   curl -s -o $f $url$f &
   sleep 5
  fi
 fi
done

# Wait for all downloads to finish
wait

# Verify if md5sum is ok, unpack the good files, and remove the unneeded archive files afterwards
echo 'Installing new database entries...'
nonew=true
for f in $files; do
 fbase=$(basename $f .tar.gz)
 # Verify that the database file does not exist before extraction
 if [ ! -f $fbase'.phr' ]; then
  # Verify that the tar.gz files exists before extraction
  if [ -f $f ]; then
   echo 'Installing new database file ('$fbase')'
   $checkmd5sum $f.md5 && tar zxpf $f && rm -f $f && mv -f $f'.md5' $f'_old.md5' &
   sleep 10
   nonew=false
  fi
 fi
done

# Wait for everything is unpacked
wait

if $nonew; then
 echo 'No new database entries were found!'
fi

# Log failed attempts
failed=false
for f in $files; do
 fbase=$(basename $f .tar.gz)
 # Verify that the database file does not exist before extraction
 if [ ! -f $fbase'.phr' ]; then
  # Verify that the tar.gz files exists before extraction
  if [ ! -f $f ]; then
   echo $f' was not downloaded!'
  else
   echo $f' was not extracted properly!'
  fi
  failed=true
 fi
done

if $failed; then
 echo 'Error not all databases was installed, try again, and see if the second time is a charm...'
fi
