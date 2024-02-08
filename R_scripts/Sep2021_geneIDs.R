
## make loop to intersect -wa -wb Araport_new with beds / shell script

module load bioinfo-tools BEDTools

in=~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/homer_beds/
out=~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/GOs_geneIDs_beds/27me3_27ac_overlap/

for f in $(find $in -maxdepth 1 -name "*.bed"); do
fnam=$(basename ${f/.bed/})
output=$out$fnam"_vs_araport.bed"

echo $f
echo $output

bedtools intersect -wa -wb -a /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/araport11_complete_NEW.bed \
-b $f > $output

done

## then extract IDs from file / Rscript

out<-"~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/GOs_geneIDs_beds/27me3_27ac_overlap/"
setwd(out)

intersected_IDs <-  dir(".", pattern="araport", recursive = F, full.names = TRUE) 

Gene_IDs <- lapply (intersected_IDs, function(x){
    table <- read.table(x, header = F, sep="\t")
    geneIDs <- table$V6
    geneIDs
})

names(Gene_IDs)<-intersected_IDs
names(Gene_IDs) <- strtrim(names(Gene_IDs), 25)
names(Gene_IDs) <- substring(names(Gene_IDs),3)

for(x in 1:length(Gene_IDs)) {
    elemento <- Gene_IDs[[x]]
    nombre <- names(Gene_IDs)[x]
    
    # Write table of gene IDs
    write.table(elemento, file = paste0(out,"gene_IDs/",nombre,"_GeneIDs.txt"), sep ="\t", row.names = F, col.names = F, quote=FALSE )
    
}

# dat = your list of genes of interest (a subset of a population)
dat <- scan("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/GOs_geneIDs_beds/27me3_27ac_overlap/gene_IDs/27ac_27me3_overlap_vs_a_GeneIDs.txt",what="character")
class(dat)
# bg = your population (think what defines it.)
bg <- read.table("/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/gopher/tair10_gene_to_go.tsv", sep = "\t")
bg <- bg$V1
class(bg)
# we silence warnings (not good practice, but we know what we're doing) - think to adjust the path
source("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/src/gopher.R")

# we just quantify the run time
# task has to be a list
# task can take any value from: go, kegg, mapman and pfam
enrichment <- gopher(genes=dat,task ="go",background = NULL, url="athaliana")

# enrichment will contain a list of tibbles (a "type of" data.frame), on per "task"

# for go, you can for example export the GO ID and the FDR to a file and then 
# upload that to REVIGO (http://revigo.irb.hr)

