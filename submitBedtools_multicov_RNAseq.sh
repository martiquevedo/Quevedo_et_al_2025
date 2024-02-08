#!/bin/bash

## be verbose and print
set -ex

proj=u2017017
mail=marti.quevedo@umu.se

#module load bioinfo-tools BEDTools

in=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/RNAseq/bams
out=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Sept2021_extra_analysis/profiles
outfile=$out/"k27ac_report_RNA_coverage.txt"
## create the out dir
    if [ ! -d $out ]; then
    mkdir -p $out
    fi
bed=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Sept2021_extra_analysis/profiles/27ac_report.bed
bams=`ls $in/*_STAR.bam`

fnam="RNAseq"
    sbatch -A $proj -t 04:00:00 -e $out/$fnam.err -o $out/$fnam.out \
    -J multicov-$fnam -p core -n 12 \
    /mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/pipeline/runBedToolsMultiCoverage.sh $bed $outfile $bams 

    
