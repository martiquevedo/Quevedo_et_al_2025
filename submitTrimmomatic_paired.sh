#!/bin/bash 

set -ex

# Temporary environment variables FROM VINCENT SCRIPT, DIDN´T WORK
#UPSCb=/mnt/picea/home/mquevedo/Git/UPSCb/
#TRIMMOMATIC_HOME=/mnt/picea/storage/Modules/apps/bioinfo/trimmomatic/0.36

# functions
source $UPSCb/src/bash/functions.sh


## variables

proj=u2017017
mail=marti.quevedo@umu.se

in=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Literature/raw/seedling_H3k4me3/deinterleaved
out=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Literature/trimmomatic/h3k4me3

## tools ## DO NOT LOAD FOR ME WHEN IN SCRIPT!... I LOAD THEM MANUALLY
#module load bioinfo-tools 
#module load trimmomatic


for f in $(find $in -name "*_1.fastq"); do  
    fnam=$(basename ${f/_1.fastq/}) 
    sbatch -A $proj -t 00:40:00  \
    -e $out/$fnam.err -o $out/$fnam.out \
    -J TRIM-$fnam -p core -n 20 \
    $UPSCb/pipeline/runTrimmomatic.sh $fnam"_1.fastq" $fnam"_2.fastq" $out
 
done


