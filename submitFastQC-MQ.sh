#!/bin/bash 

set -ex

# functions
source $UPSCb/src/bash/functions.sh

## variables

proj=u2017017
mail=marti.quevedo@umu.se

in=/mnt/picea/projects/arabidopsis/aastrand/greening/phy/raw/far_red/combined_renamed/trimmed
out=/mnt/picea/projects/arabidopsis/aastrand/greening/phy/raw/far_red/combined_renamed/trimmed/fastqc

if [ ! -d $out ]; then
mkdir -p $out
fi

for f in $(find $in -name "*.f*q".gz); do 
    fnam=$(basename ${f/.f*q.gz/})
    sbatch -A $proj -t 00:15:00  \
    -e $out/$fnam.err -o $out/$fnam.out \
    -J FQC-$fnam -p core -n 12 \
    $UPSCb/pipeline/runFastQC.sh $out $f
done