#!/bin/bash -l

proj=u2017017
mail=marti.quevedo@umu.se

#module load bioinfo-tools 
#module load BEDTools kentUtils
#does not work on symbolic links

## variables
in=//mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/bams_bwa/16M/merged/
out=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/bigwigs/
genomesize=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/tair10/genome.sizes.athailana
##samtools view -H file.bam|grep @SQ|sed 's/@SQ\tSN:\|LN://g' > genome.sizes.athaliana ## makes the genomesize file from your own BAM if it has a header ##
ext=50 # bp to extend reads

for f in $(find $in -maxdepth 1 -name "*.sorted.bam" ); do
fnam=$(basename ${f/.sorted.bam/})


sbatch -A $proj -t 02:00:00 \
-e $out/$fnam.err -o $out/$fnam.out \
-J BtbW-$fnam -p core -n 20 ~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/pipeline/runBamtoBigWig_v2.sh $f $ext $genomesize $out
done

