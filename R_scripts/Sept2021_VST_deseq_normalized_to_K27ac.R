
#Load tools
suppressPackageStartupMessages(library(DiffBind))
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library("flashClust"))
suppressPackageStartupMessages(library(ggplot2))

#Set 
setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq")

#load metadata file
meta <- read.csv("/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/Meta_histones.csv", header = T, sep = ";")
colnames(meta) <- c("X","SampleID","Histone","Time","Replicate","Group1")
meta.k4me3 <- meta[grep("h3k4me3",meta$Histone),]
meta.k4me3$Time <- as.integer(substr(meta.k4me3$Time, 4,4))



#Diffbind read,count and extract matrix
greening_4me3 <- dba(sampleSheet="/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/Sep2021_k4me3_sheet_diffbind_vsAll.csv")
greening_4me3 <- dba.count(greening_4me3, score=DBA_SCORE_READS)
k4me3_counts_matrix_raw <- dba.peakset(greening_4me3, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)


#REMOVE Chr and converting to numeric to avoid ""
k4me3_counts_matrix_raw [,1] <- gsub('Chr', '', k4me3_counts_matrix_raw [,1])
k4me3_counts_matrix_raw [,1] <- as.numeric(as.character(k4me3_counts_matrix_raw [,1]))

#Prepare matrix
k4me3_counts_matrix_raw_clean <- k4me3_counts_matrix_raw [,-c(1:3)] 
k4me3_counts_matrix_raw_coordinates <- paste(k4me3_counts_matrix_raw [,1], ':', k4me3_counts_matrix_raw [,2], ':', k4me3_counts_matrix_raw [,3], sep = "")
rownames(k4me3_counts_matrix_raw_clean)<-k4me3_counts_matrix_raw_coordinates

sapply(k4me3_counts_matrix_raw_clean, class)

#DESEQ
dds.time.4me3 <- DESeqDataSetFromMatrix(
    countData = k4me3_counts_matrix_raw_clean,
    colData = meta.k4me3,
    design = ~ Time)

#Blind = TRUE only for quality check
#VST model aware
vsd_dds.time.4me3 <- varianceStabilizingTransformation(dds.time.4me3, blind=F)
vst_dds.time.4me3 <- assay(vsd_dds.time.4me3)
vst_dds.time.4me3 <- vst_dds.time.4me3 - min(vst_dds.time.4me3)

#Save DBA.object
saveRDS(vst_dds.time.4me3, file = "/mnt/picea/home/mquevedo/R_objects/Sept2021_vst_4me3_vsAll.rds")






