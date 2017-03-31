# Setup RUCS on COMPUTEROME
module load anaconda3/4.0.0 perl/5.24.0 samtools/0.1.18 bwa/0.7.10 ncbi-blast/2.6.0+
mkdir ~/envs
cd ~/envs
virtualenv rucs
source ~/envs/rucs/bin/activate
pip install primer3-py numpy tabulate
deactivate

mkdir ~/src ~/bin
cd ~/src
git clone https://bitbucket.org/genomicepidemiology/rucs.git
ln -s ~/src/rucs/primer_core_tools.py ~/bin/rucs
cat <<EOT >> ~/envs/rucs/bin/activate

# Load required modules
module load anaconda3/4.0.0 perl/5.24.0 samtools/0.1.18 bwa/0.7.10 ncbi-blast/2.6.0+

# Add BLAST DB to ENV to enable annotation
BLASTDB=/home/databases/ncbi-blast/refseq_protein
export BLASTDB

# Add $Home/bin to $PATH
PATH="\$PATH:~/bin"
export PATH
EOT

# Check installation
source ~/envs/rucs/bin/activate
which python3 samtools bwa blastn blastx makeblastdb
python3 -c 'import gzip, json, types, shutil, glob, bisect, primer3, numpy, subprocess, difflib, tabulate'
deactivate

# Test program
mkdir -p ~/test/rucs
cd ~/test/rucs
source ~/envs/rucs/bin/activate
rucs test
deactivate

# Show program help
source ~/envs/rucs/bin/activate
rucs --help
deactivate

## Run program
#source ~/envs/rucs/bin/activate
#rucs full --positives /path/to/my/positive_genomes/* --negatives /path/to/my/negative_genomes/*
#deactivate
