setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/macs2/summits_intersects/vsInput_H3_expanded_merged/merged")

bed_list <-  dir(".", pattern="_final", recursive = F, full.names = TRUE) 


expanded_beds <- lapply (bed_list, function(x){
    y <- read.table(x, stringsAsFactors = F)
    y$V2 <- y$V2-25
    y$V3 <- y$V3+25
    y
    
})
    
names(expanded_beds)<-bed_list
names(expanded_beds) <- strtrim(names(expanded_beds), 15)
names(expanded_beds) <- substring(names(expanded_beds),3)

for(x in 1:length(expanded_beds)) {
    elemento <- expanded_beds[[x]]
    nombre <- names(expanded_beds)[x]
    
    # Write table of gene IDs
    write.table(elemento, file = paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/macs2/summits_intersects/vsInput_H3_expanded_merged/merged/",nombre,"_expanded.bed"), sep ="\t", row.names = F, col.names = F, quote=FALSE )
    
}
