#!/bin/bash
### Install BLAST DB (refseq_protein) ###
# USAGE install_db.sh /blast/db

# Get NCBI BLAST refseq protein database
# Set DB DIR
if [ -z "$1" ]; then
   export BLASTDB=`pwd`
elif [ -z "$BLASTDB" ]; then
   echo ""
else
   export BLASTDB=/blast/db
fi

mkdir $BLASTDB
cd $BLASTDB

# Download database MD5 sums
curl -s -o refseq_protein.00.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.00.tar.gz.md5 &
curl -s -o refseq_protein.01.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.01.tar.gz.md5 &
curl -s -o refseq_protein.02.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.02.tar.gz.md5 &
curl -s -o refseq_protein.03.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.03.tar.gz.md5 &
curl -s -o refseq_protein.04.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.04.tar.gz.md5 &
curl -s -o refseq_protein.05.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.05.tar.gz.md5 &
curl -s -o refseq_protein.06.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.06.tar.gz.md5 &
curl -s -o refseq_protein.07.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.07.tar.gz.md5 &
curl -s -o refseq_protein.08.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.08.tar.gz.md5 &
curl -s -o refseq_protein.09.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.09.tar.gz.md5 &
curl -s -o refseq_protein.10.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.10.tar.gz.md5 &
curl -s -o refseq_protein.11.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.11.tar.gz.md5 &
curl -s -o refseq_protein.12.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.12.tar.gz.md5 &
curl -s -o refseq_protein.13.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.13.tar.gz.md5 &
curl -s -o refseq_protein.14.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.14.tar.gz.md5 &
curl -s -o refseq_protein.15.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.15.tar.gz.md5 &
curl -s -o refseq_protein.16.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.16.tar.gz.md5 &
curl -s -o refseq_protein.17.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.17.tar.gz.md5 &
curl -s -o refseq_protein.18.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.18.tar.gz.md5 &
curl -s -o refseq_protein.19.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.19.tar.gz.md5 &
curl -s -o refseq_protein.20.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.20.tar.gz.md5 &
curl -s -o refseq_protein.21.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.21.tar.gz.md5 &
curl -s -o refseq_protein.22.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.22.tar.gz.md5 &
curl -s -o refseq_protein.23.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.23.tar.gz.md5 &
curl -s -o refseq_protein.24.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.24.tar.gz.md5 &
curl -s -o refseq_protein.25.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.25.tar.gz.md5 &
curl -s -o refseq_protein.26.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.26.tar.gz.md5 &
curl -s -o refseq_protein.27.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.27.tar.gz.md5 &
curl -s -o refseq_protein.28.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.28.tar.gz.md5 &
curl -s -o refseq_protein.29.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.29.tar.gz.md5 &
curl -s -o refseq_protein.30.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.30.tar.gz.md5 &

# Download database parts
curl -s -o refseq_protein.01.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.01.tar.gz &
curl -s -o refseq_protein.02.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.02.tar.gz &
curl -s -o refseq_protein.03.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.03.tar.gz &
curl -s -o refseq_protein.04.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.04.tar.gz &
curl -s -o refseq_protein.05.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.05.tar.gz &
curl -s -o refseq_protein.06.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.06.tar.gz &
curl -s -o refseq_protein.07.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.07.tar.gz &
curl -s -o refseq_protein.08.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.08.tar.gz &
curl -s -o refseq_protein.09.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.09.tar.gz &
curl -s -o refseq_protein.10.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.10.tar.gz &
curl -s -o refseq_protein.11.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.11.tar.gz &
curl -s -o refseq_protein.12.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.12.tar.gz &
curl -s -o refseq_protein.13.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.13.tar.gz &
curl -s -o refseq_protein.14.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.14.tar.gz &
curl -s -o refseq_protein.15.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.15.tar.gz &
curl -s -o refseq_protein.16.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.16.tar.gz &
curl -s -o refseq_protein.17.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.17.tar.gz &
curl -s -o refseq_protein.18.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.18.tar.gz &
curl -s -o refseq_protein.19.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.19.tar.gz &
curl -s -o refseq_protein.20.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.20.tar.gz &
curl -s -o refseq_protein.21.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.21.tar.gz &
curl -s -o refseq_protein.22.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.22.tar.gz &
curl -s -o refseq_protein.23.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.23.tar.gz &
curl -s -o refseq_protein.24.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.24.tar.gz &
curl -s -o refseq_protein.25.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.25.tar.gz &
curl -s -o refseq_protein.26.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.26.tar.gz &
curl -s -o refseq_protein.27.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.27.tar.gz &
curl -s -o refseq_protein.28.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.28.tar.gz &
curl -s -o refseq_protein.29.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.29.tar.gz &
curl -s -o refseq_protein.30.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.30.tar.gz &
curl -s -o refseq_protein.00.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.00.tar.gz 

# Sleep 10 minutes to make sure all database fragments are downloaded
sleep 600

# Check MD5 sums and unpack good downloads
for f in `ls refseq_protein.*.tar.gz`
do
 md51=`md5 -q $f`
 md52=`cut -d' ' -f1 $f.md5` 
 if [ "$md51" = "$md52" ]; then
  echo "$f - OK"
  tar zxpf $f
 else
  echo "$f - FAIL"
  echo "$md51 != $md52" 
 fi
done

# If all files were downloaded and unpacked correctly, do some clean up:
#rm -f refseq_protein.*.tar.gz*
