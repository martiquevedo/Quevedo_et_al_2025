
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
meta <- read.csv("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/Meta_histones.csv", header = T, sep = ";")
colnames(meta) <- c("X","SampleID","Histone","Time","Replicate","Group1")
meta.k27ac <- meta[grep("h3k27ac",meta$Histone),]
meta.k27ac$Time <- as.integer(substr(meta.k27ac$Time, 4,4))



#Diffbind read,count and extract matrix
greening_27ac <- dba(sampleSheet="~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/Sep2021_k27ac_sheet_diffbind.csv")
greening_27ac <- dba.count(greening_27ac, score=DBA_SCORE_READS)
k27ac_counts_matrix_raw <- dba.peakset(greening_27ac, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)


#REMOVE Chr and converting to numeric to avoid ""
k27ac_counts_matrix_raw [,1] <- gsub('Chr', '', k27ac_counts_matrix_raw [,1])
k27ac_counts_matrix_raw [,1] <- as.numeric(as.character(k27ac_counts_matrix_raw [,1]))

#Prepare matrix
k27ac_counts_matrix_raw_clean <- k27ac_counts_matrix_raw [,-c(1:3)] 
k27ac_counts_matrix_raw_coordinates <- paste(k27ac_counts_matrix_raw [,1], ':', k27ac_counts_matrix_raw [,2], ':', k27ac_counts_matrix_raw [,3], sep = "")
rownames(k27ac_counts_matrix_raw_clean)<-k27ac_counts_matrix_raw_coordinates

sapply(k27ac_counts_matrix_raw_clean, class)

#DESEQ
dds.time.27ac <- DESeqDataSetFromMatrix(
    countData = k27ac_counts_matrix_raw_clean,
    colData = meta.k27ac,
    design = ~ Time)

#Blind = TRUE only for quality check
#VST model aware
vsd_dds.time.27ac <- varianceStabilizingTransformation(dds.time.27ac, blind=F)
vst_dds.time.27ac <- assay(vsd_dds.time.27ac)
vst_dds.time.27ac <- vst_dds.time.27ac - min(vst_dds.time.27ac)

#Save DBA.object
saveRDS(vst_dds.time.27ac, file = "~/R_objects/Sept2021_vst_27ac.rds")






