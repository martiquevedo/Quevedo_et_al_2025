#!/bin/bash-l

# usage 
USAGETXT=\
"
Usage: perl <Program> <TAIR_GOs> <GO list> <Universe> <Input file> <output>
"

for f in $(find $1 -maxdepth 1 -name "*.txt"); do
    fnam=$(basename ${f/_geneIDs.txt/})
    out=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/profile_approach/GOs/processed/27ac/
    c1=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-mediator-stress/CDK8_analysis/scripts/GeneMerge1.4.pl
    c2=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-mediator-stress/CDK8_analysis/GO/ATH_GO_from_TAIR_2018-09-07_GeneMerge_format_OK_without_repeated_GOs_OK                   
    c3=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-mediator-stress/CDK8_analysis/GO/GO.XX.use                                                                                 
    c4=/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/GO_genemerge/Araport_universe.txt
    c5=$f
    c6=$out$fnam

perl $c1 $c2 $c3 $c4 $c5 $c6

done