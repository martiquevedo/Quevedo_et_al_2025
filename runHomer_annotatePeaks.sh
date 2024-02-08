#!/bin/bash

set -ex

# functions
source $UPSCb/src/bash/functions.sh

# usage 
USAGETXT=\
"
Usage: runHomer_findMotifsGenome.sh <peak file> <outdir>
"

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/nfcore-chipseq-20190626.simg \
annotatePeaks.pl $1  /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/tair10  > $2