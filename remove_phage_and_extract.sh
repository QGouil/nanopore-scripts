#!/bin/bash
#PBS -q huge
#PBS -m abe
#PBS -M gouil.q@wehi.edu.au
#
# Usage: ./remove_phage_and_extract.sh /path/to/reads phage_reads_subdir
# Run from AGRF_Data folder.

module load bwa
module load samtools
module load nanopolish
module load parallel
module load bamtools

set -x
cd $PBS_O_WORKDIR
READS_DIR=$1
PARENT_DIR=$(dirname $READS_DIR)
OUT_NAME=$(basename $READS_DIR)
SCRIPTS_DIR=$PBS_O_HOME/wehisan_home/scripts/quentin/nanopore-scripts
NPROC=$(nproc)

PHAGE_DIR=$PARENT_DIR/$4
mkdir -p $PHAGE_DIR

# extract fasta
FASTA=$PHAGE_DIR/${OUT_NAME}.fasta
find $READS_DIR -name "*.fast5" | \
parallel -X nanopolish extract {} > $FASTA

# align phage
PHAGE_BAM=$PHAGE_DIR/${OUT_NAME}.sorted.bam
bwa mem -x ont2d -t $((NPROC-1)) ~/genomes/phage_lambda/phage_lambda_NC_001416.1.calibration_strand.fa $FASTA | samtools sort -@ 4 -T $PARENT_DIR/samtools.tmp -o $PHAGE_BAM
samtools index $PHAGE_BAM
  
# remove phage
bamtools split -in $PHAGE_BAM -mapped
python $SCRIPTS_DIR/remove_phage_reads.py $PHAGE_DIR/${OUT_NAME}.sorted.MAPPED.bam $FASTA $PHAGE_DIR

# extract fastq
mkdir -p $PARENT_DIR/fasta
FASTQ=$PARENT_DIR/fasta/${OUT_NAME}.fastq
find $READS_DIR -name "*.fast5" | \
  parallel -X nanopolish extract -q {} > $FASTQ