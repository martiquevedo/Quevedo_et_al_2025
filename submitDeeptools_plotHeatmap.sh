#!/bin/bash

## be verbose and print
set -ex

proj=u2017017
mail=marti.quevedo@umu.se

# functions
# module load conda
# conda init zsh
# module unload conda



for histone in {"k27ac","k27me3","k4me3","k9me2"}; do

    bedfile='/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/heatmap/peaks_d7_bymark/d7_'$histone'_center.bed'
    bigwigs=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/heatmap/bigwigs
    outdir='/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/heatmap/'$histone
    outname=$histone
    ## create the out dir
    if [ ! -d $outdir ]; then
    mkdir -p $outdir
    fi
    
    fnam=$outname
    sbatch -A $proj -t 04:00:00 -e $outdir/$fnam.err -o $outdir/$fnam.out \
    -J matrix-$fnam -p core -n 12 \
    /mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/pipeline/runDeepTools_computematrix.R $bedfile $bigwigs $outdir $outname 

done


