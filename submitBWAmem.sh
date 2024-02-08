#!/bin/bash

## be verbose and print
set -ex

## source functions
source $UPSCb/src/bash/functions.sh

## tools
#module load bioinfo-tools bwa samtools


## variables
proj=u2017017
mail=marti.quevedo@umu.se
in=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Literature/trimmomatic/h3k9me2_h3k27me3
out=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Literature/bams
inx=/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/indices/bwa/genome.fa

## check vars
if [ -z $UPSCb ]; then
abort "The UPSCb var needs to be set."
fi

## create the out dir
if [ ! -d $out ]; then
mkdir -p $out
fi

## execute
for f in $(find $in -name "*.fq.gz"); do
fnam=$(basename ${f/.f*q.gz/})
sbatch -A $proj -t 00:40:00 -e $out/$fnam.err -o $out/$fnam.out \
-J BWA-$fnam -p core -n 12 $UPSCb/pipeline/runBwamem.sh -t 12 -s -f $f -i $inx -o $out 
done

