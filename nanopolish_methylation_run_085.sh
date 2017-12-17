#!/bin/bash
#PBS -q small
#PBS -m abe
#PBS -M gouil.q@wehi.edu.au

# Usage: qsub -F "path/to/genome.fasta absolute_or_relative_path/to/reads.fastq" nanopolish_methylation_run.sh
# point to the genome stored on torquelord. Give the full extension. The bwa mem index should be .fa.pac etc.
# the reads will be in fasta folder.


set -x
cd $PBS_O_WORKDIR
GENOME=$(realpath $1)
FASTA=$(realpath $2)

if [ $(echo $FASTA | grep -c -e "a$") -gt 0 ]; then 
  FMT="fasta"
elif [ $(echo $FASTA | grep -c -e "q$") -gt 0 ]; then 
  FMT="fastq"
else
  echo "ERROR: $FASTA format not recognised"
  exit 1
fi

N=10000
TMP_DIR="$PBS_O_HOME/tmp/$(dirname $FASTA)"
SCRIPTS_DIR=~/wehisan_home/scripts/quentin/nanopore-scripts

mkdir -p $TMP_DIR
cd $TMP_DIR
if [ ! -f $(basename $FASTA).1.$FMT ]; then
	python $SCRIPTS_DIR/split_fasta.py $FASTA $N
	echo "splitting reads for parallelisation"
fi
ARRAY_ID=$(qsub -F "$GENOME $FASTA $TMP_DIR" -t 1-$(ls -1 $TMP_DIR/$(basename $FASTA).*.$FMT | wc -l) $SCRIPTS_DIR/nanopolish_methylation_085.sh)
qsub -W "depend=afteranyarray:$ARRAY_ID" -F "$FASTA $TMP_DIR" $SCRIPTS_DIR/nanopolish_methylation_clean.sh

