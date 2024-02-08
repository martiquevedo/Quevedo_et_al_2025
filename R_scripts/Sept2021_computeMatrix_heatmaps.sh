for bed in $(find /mnt/e/bioinfo/side_quests/heatmap/peaks_d7_bymark/original_peaks/ -name '*.bed');do
  bedname=$(basename ${bed/.bed/})
  #mkdir '/mnt/e/bioinfo/side_quests/heatmap/output/'$bedname
  for i in $(find /mnt/e/bioinfo/side_quests/heatmap/bigwigs/ -name '*.bw');do
    fnam=$(basename ${i/.bw/})
    computeMatrix reference-point \
    --referencePoint center \
    -R '/mnt/e/bioinfo/side_quests/heatmap/peaks_d7_bymark/original_peaks/'$bedname'.bed' \
    -S $i \
    -b 3000 -a 3000 \
    -bs 10 -p "max" --missingDataAsZero \
    --skipZeros -o '/mnt/e/bioinfo/side_quests/heatmap/output/'$bedname'/'$fnam'_'$bedname'_peaks_multi.gz' \
    --outFileNameMatrix '/mnt/e/bioinfo/side_quests/heatmap/output/'$bedname'/'$fnam'_'$bedname'_peaks_multi.tab' \
    --outFileSortedRegions '/mnt/e/bioinfo/side_quests/heatmap/output/'$bedname'/'$fnam'_'$bedname'_peaks_multi.bed'  
    
    
  done
done

colors=( royalblue royalblue darkorange darkorange rebeccapurple rebeccapurple indianred indianred)
min_scale=( -100 -20 -100 -100 -100 -100 -50 -10 )
max_scale=( 300 70 300 300 350 300 150 20 )
i=-1
for matrix in $(find /mnt/e/bioinfo/side_quests/heatmap/output/day7_k9me2/ -maxdepth 1 -name '*.gz');do
outname=$(basename ${matrix/_multi.gz/})
i=$((i+1))
plotHeatmap -m $matrix \
-out '/mnt/e/bioinfo/side_quests/heatmap/graphs/k9me2_peaks/'$outname'.png' \
--heatmapHeight 5 \
--whatToShow 'heatmap and colorbar' \
--zMin ${min_scale[i]} --zMax ${max_scale[i]}   \
--colorList 'black, white, '${colors[i]}
done

set -ex

# functions
source $UPSCb/src/bash/functions.sh

# test
#isEnvVarSet $UPSCb

# UNION of all PEAKS, needs to be modified by n° of samples 

# usage
USAGETXT=\
"
Usage: runbedops.sh <output directory> <sample 1> <sample 2> <sample3> <sample4>
"


# test if $1 is -d

if [ ! -d $(dirname $1) ]; then
abort "The output directory $(dirname $1) does not exist."
fi

# test if $2 - $3 are -f

if [ ! -f $2 ]; then
abort "Sample 1 does not exist."
fi

if [ ! -f $3 ]; then
abort "Sample 2 does not exist."
fi


fnam=$(basename ${2/.bed/})

bedops -u $2 $3 $4 $5 > $1"/"$fnam".bed"
