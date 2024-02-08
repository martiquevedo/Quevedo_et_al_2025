#after bedtools intersect -wa -wb -a peak_file -b /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/araport11_complete_NEW.bed > peak_araport.txt
#then create a fuseTable of all HISTONEpeaks_vsAraport

#Convert a Chr START END file into 1st a FUSED REGION and then match to specfic ARAPORT FUSED_ID

setwd("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/homer_beds/k27ac/")
out="~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/GOs_geneIDs_beds/k27ac/gene_IDs"
i=27ac

intersected_IDs <-  dir(".", pattern="27ac", recursive = F, full.names = TRUE) 

nombres_lista <- as.list(intersected_IDs)
sub(".bed", "", nombres_lista)

x<-1
Gene_IDs <- lapply (1:length(intersected_IDs), function(x){
    elemento <- read.table(intersected_IDs[[x]])
    nombre <- sub("./", "", nombres_lista[[x]])
    nombre <- sub("./", "", nombre)
    
    elemento[,1] <- gsub('Chr', '', elemento[,1])
    elemento[,1] <- as.numeric(as.character(elemento[,1]))
    fused_region <- data.frame(region = paste(elemento$V1, elemento$V2, elemento$V3, sep = ":"))
    
    # Load araport annotation
    fused_reference <- read.csv(file = paste0("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/araport_annotation/k27me3_summits_araport_fused.txt"),sep = "\t", header = T)
    
    # Match 2 lists
    fused_region_merged <-fused_reference$gene_ID[fused_reference$region %in% fused_region$region]
    
    # Split multiple entries in gene_ID
    gene_IDs_list <- unlist(strsplit(as.character(fused_region_merged), ","))
    
    # Write table of gene IDs
    write.table(gene_IDs_list, file = paste0(out,nombre,"_geneIDs.txt"), sep ="", row.names = F, col.names = F, quote=FALSE )
    
})


   