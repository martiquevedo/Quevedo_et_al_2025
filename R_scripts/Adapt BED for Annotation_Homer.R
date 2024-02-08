
setwd("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/v.2.0/Diffbind_reports/forGOs_new/enhancers/")

   
addPEAK_IDs <-  dir(".", pattern=".bed", recursive = F, full.names = TRUE) 


PEAK_IDs <- lapply (1:length(addPEAK_IDs), function(x) {
    tablea <- addPEAK_IDs[[x]]
    elemento <- read.table(tablea,sep="\t")
    elemento$V4 <- seq(1:length(elemento$V3))
    elemento$V1 <- gsub('Chr', '', elemento$V1 )
    elemento$V5 <- rep("0",length(elemento$V3))
    elemento <- elemento[, c(1, 2, 3, 5)]
    elemento
})

names(PEAK_IDs)<-addPEAK_IDs
names(PEAK_IDs) <- strtrim(names(PEAK_IDs), 25)
names(PEAK_IDs) <- substring(names(PEAK_IDs),3)

for(x in 1:length(PEAK_IDs)) {
    elemento <- PEAK_IDs[[x]]
    nombre <- names(PEAK_IDs)[x]
    
    # Write table of gene IDs
    write.table(elemento, file = paste0("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/v.2.0/Diffbind_reports/forGOs_new/enhancers/",nombre,"_ffannotation.txt"), sep ="\t", row.names = F, col.names = F, quote=FALSE )
    
}