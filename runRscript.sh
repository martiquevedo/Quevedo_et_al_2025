#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH -A u2017017

module load R

Rscript /mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/pipeline/Rscripts/Sept2021_VST_deseq.R