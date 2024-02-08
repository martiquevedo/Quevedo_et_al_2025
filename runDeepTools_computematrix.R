#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL


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

computeMatrix reference-point \
--referencePoint center \
-R $bed_file \
-S $bw_folder"/"*.bw \
-b 3000 -a 3000 \
-bs 10 -p "max" --missingDataAsZero \
--skipZeros -o $outdir"/"$outname".gz" \
--outFileNameMatrix $outdir"/"$outname".tab" \
--outFileSortedRegions $outdir"/"$outname".bed"  



