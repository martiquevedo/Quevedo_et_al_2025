#!/bin/bash

set -ex

##Variable z to jump from folders day0,1,4,7 (both in and out)

## variables
proj=u2017017
mail=marti.quevedo@umu.se

## source functions
source $UPSCb/src/bash/functions.sh

### Tools
#module load bioinfo-tools macs2


for z in {0,1,4,7}; do

in=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/bams_bwa/16M/day$z
contr_dir=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/bams_bwa/16M/day$z/H3/
out=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/macs2/summits_intersects/vsH3_0.05

## create the out dir
if [ ! -d $out ]; then
mkdir -p $out
fi


    ## execute
    for f in $(find $in -maxdepth 1 -name "*.bam"); do
            for v in $(find $contr_dir -name "*.bam"); do
            contr=$v
            done
        fnam=$(basename ${f/_16M.sorted.bam/})
        sbatch -A $proj -t 00:45:00 -e $out/$fnam.err -o $out/$fnam.out \
        -J MACS-$fnam -p core -n 12 $UPSCb/pipeline/runMacs2_broad.sh -g 120000000  -n "$fnam""_0.05_vsH3" $f $contr $out
    done
done 






   