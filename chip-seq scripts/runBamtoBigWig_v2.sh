#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL

## stop on error
set -ex


## Create a bigWig file for visiualization
usage(){
    echo >&2 \
    " Usage: $0  <bam> <ext> <genome.sizes.file> <out>
    "
    exit 1
}


file=$1
shift;
ext=$1
shift;
genome_sizes=$1
shift;
out=$1
shift;

dir=$(dirname $file)"/"
input=$(basename $file .bam)
input_data=$dir$input

bamToBed -i $input_data".bam"  > $out"/"$input".bed" #BEDTools

#extend reads
slopBed -i $out"/"$input".bed" -g $genome_sizes -b $ext > $out"/"$input"_extended_"$ext"bp.bed"
rm $out"/"$input".bed"

genomeCoverageBed -i $out"/"$input"_extended_"$ext"bp.bed" -bg -g $genome_sizes > $out"/"$input"_extended_"$ext"bp.bedgraph"
rm $out"/"$input"_extended_"$ext"bp.bed"

bedGraphToBigWig $out"/"$input"_extended_"$ext"bp.bedgraph" $genome_sizes $out"/"$input"_nolowup_extended_"$ext"bp.bw" #kentUtils
rm $out"/"$input"_extended_"$ext"bp.bedgraph"
