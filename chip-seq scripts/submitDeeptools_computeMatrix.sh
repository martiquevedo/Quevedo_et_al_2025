#!/bin/bash

## be verbose and print
set -ex

proj=u2017017
mail=marti.quevedo@umu.se

# functions
# module load conda
# conda init zsh
# module unload conda

bed_files='/mnt/e/bioinfo/side_quests/heatmap/peaks_d7_bymark/original_peaks/' # folder where coordinates to compose matrix
out_dir='/mnt/e/bioinfo/side_quests/heatmap/output/' #general folder where to create subfolders of results
bigwigs='/mnt/e/bioinfo/side_quests/heatmap/bigwigs/' #folder of the bigwigs to compute

for bed in $(find $bed_files -name '*.bed');do #CHECK LOOP I CREATED AT HOME
  bedname=$(basename ${bed/.bed/})
  out=$outdir$bedname
  out_path_prefix=$out'/'$fnam'_'$bedname
  if [ ! -d $out ]; then
   mkdir -p $out
  fi
  
  
  for i in $(find /mnt/e/bioinfo/side_quests/heatmap/bigwigs/ -name '*.bw');do
    fnam=$(basename ${i/.bw/})
    computeMatrix reference-point \
    --referencePoint center \
    -R '/mnt/e/bioinfo/side_quests/heatmap/peaks_d7_bymark/original_peaks/'$bedname'.bed' \
    -S $i \
    -b 3000 -a 3000 \
    -bs 10 -p "max" --missingDataAsZero \
    --skipZeros -o $out_path_prefix'_peaks_multi.gz' \
    --outFileNameMatrix $out_path_prefix'_peaks_multi.tab' \
    --outFileSortedRegions $out_path_prefix'_peaks_multi.bed'  
    
    
  done
done
