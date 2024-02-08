#!/bin/bash

#### ADAPT VARIABLES TO SUBMITTER!!!!

set -ex

# functions
source $UPSCb/src/bash/functions.sh

# usage 
USAGETXT=\
"
Usage: runHomer_findMotifsGenome.sh <peak file> <outdir>
"

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/deepTools_v3.5.1.sif /conda/bin/computeMatrix
computeMatrix reference-point \
    --referencePoint center \
    -R '/mnt/e/bioinfo/side_quests/heatmap/peaks_d7_bymark/original_peaks/'$bedname'.bed' \
    -S $i \
    -b 3000 -a 3000 \
    -bs 10 -p "max" --missingDataAsZero \
    --skipZeros -o '/mnt/e/bioinfo/side_quests/heatmap/output/'$bedname'/'$fnam'_'$bedname'_peaks_multi.gz' \
    --outFileNameMatrix '/mnt/e/bioinfo/side_quests/heatmap/output/'$bedname'/'$fnam'_'$bedname'_peaks_multi.tab' \
    --outFileSortedRegions '/mnt/e/bioinfo/side_quests/heatmap/output/'$bedname'/'$fnam'_'$bedname'_peaks_multi.bed'  
    
annotatePeaks.pl $1  /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/tair10  > $2