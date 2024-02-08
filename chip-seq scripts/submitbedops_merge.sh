#!/bin/bash

## be verbose and print
set -ex

proj=u2017017
mail=marti.quevedo@umu.se

# module load bioinfo-tools bedops

for z in {0,1,4,7}; do

in=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/macs2/initial_tests/vsH3/day$z
out=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/macs2/summits_intersects

## execute
        for x in {"k27ac","k27me3","k4me3","k9me2"}; do
            repu1=$in"/$x"_"day"$z"_rep1_vsH3_summits.bed"
            repu2=$in"/$x"_"day"$z"_rep2_vsH3_summits.bed"
            #cp $repu1 $in"/day"$z"_"$x"_rep1_16M_q0.05_vsInput_peaks.bed"
            #cp $repu2 $in"/day"$z"_"$x"_rep2_16M_q0.05_vsInput_peaks.bed"
            #rep1=$in"/day"$z"_"$x"_rep1_16M_q0.05_vsInput_peaks.bed"
            #rep2=$in"/day"$z"_"$x"_rep2_16M_q0.05_vsInput_peaks.bed"
            
            fnam=$x"_day"$z"_rep1_rep2_summits"

            #execute
        	sbatch -A $proj -t 00:20:00 --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
        	-J Join-$fnam -p core -n 12 \
        	/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/pipeline/runbedops_merge.sh $out $repu1 $repu2

           
            
        done
    
done

