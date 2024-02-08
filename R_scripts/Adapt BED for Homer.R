setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/macs2/summits_intersects/vsInput_H3_expanded_merged/expanded")

elemento <- read.table("day7_k27ac_expanded.bed",sep="\t")
nombre <- names(profileList)[x]

# FOR HOMER
out1="/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/homer_motifs/all_peaks"
#REMOVE Chr and converting to numeric to avoid ""
elemento[,1] <- gsub('Chr', '', elemento[,1])
elemento[,1] <- as.numeric(as.character(elemento[,1]))
write.table(elemento[1:3], file = c(paste0(out1,"/k27ac_day7_allpeaks_forhomer.bed")), append = FALSE, sep="\t", dec = ".",
            row.names = FALSE, col.names = F)
