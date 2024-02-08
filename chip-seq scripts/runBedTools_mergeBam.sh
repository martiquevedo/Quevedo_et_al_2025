#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
## -A and --mail-user set in the submit job

## stop on error
set -ex

## we get one bed file with coordinates, multifiles (store in one variable) to calculate coverage and an output dir
usage(){
  echo >&2 \
  " Usage: $0  <outfile> <rep1_bam> <rep2_bam> <rep1_bai> <rep2_bai>
    Note: <outfile> is the final merged 
    <rep1_bam> BAM file
    <rep1_bai> index file
  "
  exit 1
}


outfile=$1
shift;

rep1_bam=$1
shift;

rep2_bam=$1
shift;

# rep1_bai=$1
# shift;
# 
# rep2_bai=$1
# shift;

## get the subtracted results
samtools merge $outfile $rep1_bam $rep2_bam #$rep1_bai $rep2_bai

