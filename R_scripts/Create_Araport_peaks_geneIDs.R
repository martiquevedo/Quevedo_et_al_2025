#after bedtools intersect -wa -wb -a peak_file -b /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/araport11_complete_NEW.bed > peak_araport.txt

setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/summits_expanded/interesected_araport/")

intersected_IDs <-  dir(".", pattern="araport", recursive = F, full.names = TRUE) 

Gene_IDs <- lapply (intersected_IDs, function(x){
  table <- read.table(x, header = F, sep="\t")
  geneIDs <- table$V9
  geneIDs
  })

names(Gene_IDs)<-intersected_IDs
names(Gene_IDs) <- strtrim(names(Gene_IDs), 25)
names(Gene_IDs) <- substring(names(Gene_IDs),3)

for(x in 1:length(Gene_IDs)) {
  elemento <- Gene_IDs[[x]]
  nombre <- names(Gene_IDs)[x]
  
  # Write table of gene IDs
  write.table(elemento, file = paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/summits_expanded/interesected_araport/geneIDs/",nombre,"_GeneIDs.txt"), sep ="\t", row.names = F, col.names = F, quote=FALSE )
  
}


