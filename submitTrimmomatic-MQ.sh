#!/bin/bash 

set -ex

#UPSCb=/mnt/picea/home/mquevedo/Git/UPSCb/
#TRIMMOMATIC_HOME=/mnt/picea/storage/Modules/apps/bioinfo/trimmomatic/0.36

# functions
source $UPSCb/src/bash/functions.sh


## variables

proj=u2017017
mail=marti.quevedo@umu.se

in=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Literature/raw
out=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Literature/trimmomatic/h3k9me2_h3k27me3

## tools ## DO NOT LOAD FOR ME WHEN IN SCRIPT!... I LOAD THEM MANUALLY
#module load bioinfo-tools 
#module load trimmomatic
# ADAPTED MINLEN:30 for fastQC lenght


for f in $(find $in -maxdepth 1 -name "*.f*q.gz"); do  
    fnam=$(basename ${f/.f*q.gz/}) 
    sbatch -A $proj -t 00:40:00  \
    -e $out/$fnam.err -o $out/$fnam.out \
    -J TRIM-$fnam -p core -n 20 \
    $UPSCb/pipeline/runTrimmomatic.sh -s  $f $out SLIDINGWINDOW:5:20 MINLEN:30
 
done


