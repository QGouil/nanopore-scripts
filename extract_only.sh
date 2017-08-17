#!/bin/bash
#PBS -q huge
#PBS -m abe
#PBS -M gouil.q@wehi.edu.au

# Usage: qsub -F "relative/path/to/basecalled_reads" ./extract_only.sh
# Only use phage_reads_subdir if a calibration strand has been added
# Run from AGRF_Data folder.

module load nanopolish
module load parallel


set -x
cd $PBS_O_WORKDIR
READS_DIR=$1
PARENT_DIR=$(dirname $READS_DIR)
OUT_NAME=$(basename $READS_DIR)
SCRIPTS DIR=~/wehisan_home/scripts/quentin/nanopore-scripts
NPROC=$(nproc)


# extract fastq
mkdir -p $PARENT_DIR/fasta
FASTQ=$PARENT_DIR/fasta/${OUT_NAME}.fastq
find $READS_DIR -name "*.fast5" | \
  parallel -X nanopolish extract -q {} > $FASTQ

