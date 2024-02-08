#!/bin/bash

## be verbose and print
set -ex

## variables
proj=u2017017
mail=marti.quevedo@umu.se

in=~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/profiles/beds/

## execute
for f in $(find $in -maxdepth 1 -name "*.bed"); do
  fnam=$(basename ${f/.bed/})
  dir=$(dirname $f)"/"
  input_data=$dir$fnam
  output=$input_data
  mkdir $output
  sbatch -A $proj -t 00:40:00 -e $in/$fnam.err \
  -J Hfm-$fnam -p core -n 12 $UPSCb/projects/arabidopsis-greening-ChIP-Seq/pipeline/runHomer_findMotifsGenome.sh $f $output
done 
