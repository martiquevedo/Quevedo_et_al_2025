#!/bin/bash

## be verbose and print
set -ex

proj=u2017017
mail=marti.quevedo@umu.se

# module load bioinfo-tools bedops

in=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/macs2/replicate_merging/H3_Input_common/final_merge/
out=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/macs2/replicate_merging/H3_Input_common/final_merge/final

## execute
        for x in {"k27ac","k27me3","k4me3","k9me2"}; do
            
            rep1=$in"day0_"$x".bed"
            rep2=$in"day1_"$x".bed"
            rep3=$in"day4_"$x".bed"
            rep4=$in"day7_"$x".bed"
            
            fnam="all"$x

            #execute
        	sbatch -A $proj -t 00:20:00 --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
        	-J Join-$fnam -p core -n 12 \
        	/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/pipeline/runbedops_everything.sh $out $rep1 $rep2 $rep3 $rep4

           
            
        done
  