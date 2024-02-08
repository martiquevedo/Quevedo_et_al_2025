#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
## -A and --mail-user set in the submit job

## stop on error
set -ex

## we get one bed file with coordinates, multifiles (store in one variable) to calculate coverage and an output dir
usage(){
  echo >&2 \
  " Usage: $0  <a file> <out dir> <b files> [bed coverage option] 
    Note: 'a' is the 'subject' file, the file that contains the intervals
          of interest (e.g. the genes)
          'b' is the 'query' file, the file that contains the intervals to
          be counted/summarized (e.g. the reads)
  "
  exit 1
}

if [ "$#" -lt 3 ]; then
echo "This function requires 3 arguments"
usage;
fi

if [ ! -f $1 ]; then
echo "The first argument needs to be an existing bed file"    
usage;
fi
a=$1
shift;

outfile=$1
shift;

b=$1

## get the subtracted results
bedtools multicov  -bams $b $@ -bed $a > $outfile



