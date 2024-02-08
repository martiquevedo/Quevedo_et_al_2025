#!/bin/bash

set -ex

# functions
# module load conda
# conda init zsh
# module unload conda
#module load bioinfo-tools deepTools
#conda activate deepTools

# usage
USAGETXT=\
"
Usage: computeMatrix  <bedfile> <bwfolder> <outdir> <outname> 
    With 1 peak coordinate bed file, computes matrix for all bigwigs in a folder and outputs into a directory
"
if [ "$#" -lt 4 ]; then
echo "This function requires 4 arguments"
usage;
fi


bed_file=$1
shift;

bw_folder=$1
shift;

outdir=$1
shift;

outname=$1
shift;


plotHeatmap -m ./k27me3/matrix_k27me3_MultiBW.gz -out k27me3_day7.png \
--heatmapHeight 15  --refPointLabel 'peak_center' \
--regionsLabel 'peaks' --plotTitle '27me3 signal' \
--zMin 0 0 0 0 --zMax 120 120 120 120 \
--colorList 'white, blue' 'white, darkorange' 'white, forestgreen' 'white, red'



