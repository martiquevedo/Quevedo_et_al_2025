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

setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Sept2021_extra_analysis/profiles")

#READ RNA coverage (from bedtools multicov, bams/bais)
rna_count_27ac <- read.table("k27ac_report_RNA_coverage.txt", stringsAsFactors = F)
rna_count_27ac$V17 <- c(rep("k27ac", length(rna_count_27ac$V1)))
head(rna_count_27ac)
# rna_count_27me3 <- read.table("k27me3_peaks_RNAseq_coverage.txt", stringsAsFactors = F)
# rna_count_27me3$V17 <- c(rep("k27me3", length(rna_count_27me3$V1)))
# 
# 
# rna_count_4me3 <- read.table("k4me3_peaks_RNA_coverage_NEW.txt", stringsAsFactors = F)
# rna_count_4me3$V17 <- c(rep("k4me3", length(rna_count_4me3$V1)))
# 
# 
# rna_count_9me2 <- read.table("k9me2_peaks_RNAseq_coverage.txt", stringsAsFactors = F)
# rna_count_9me2$V17 <- c(rep("k9me2", length(rna_count_9me2$V1)))
# 
# 
# all_raw_counts <- rbind(rna_count_27ac,rna_count_27me3,rna_count_4me3,rna_count_9me2)

#REMOVE Chr and converting to numeric to avoid ""
all_raw_counts<-rna_count_27ac
all_raw_counts [,1] <- gsub('Chr', '', all_raw_counts [,1])
all_raw_counts [,1] <- as.numeric(as.character(all_raw_counts [,1]))

all_raw_counts_clean <- all_raw_counts [,-c(1:4)] 
all_raw_counts_coordinates <- paste(all_raw_counts [,1], ':', all_raw_counts [,2], ':', all_raw_counts [,3], sep = "")
rownames(all_raw_counts_clean)<-all_raw_counts_coordinates

SampleID <- c("Day0_rep1","Day0_rep2","Day0_rep3","Day1_rep1","Day1_rep2","Day1_rep3","Day4_rep1","Day4_rep2","Day4_rep3","Day7_rep1","Day7_rep2","Day7_rep3")
Time <- as.integer(c(0,0,0,1,1,1,4,4,4,7,7,7))
Replicate <- as.integer(c(1,2,3,1,2,3,1,2,3,1,2,3))
Type <- c("RNAseq","RNAseq","RNAseq","RNAseq","RNAseq","RNAseq","RNAseq","RNAseq","RNAseq","RNAseq","RNAseq","RNAseq")    
meta <- data.frame(SampleID,Time,Replicate,Type)

colnames(all_raw_counts_clean) <- c(as.character(meta$SampleID))
nrow(all_raw_counts_clean)
head(all_raw_counts_clean)

# DESEQ 

dds.time.RNA.peaks <- DESeqDataSetFromMatrix(
    countData = all_raw_counts_clean[,1:12],
    colData = meta,
    design = ~ Time)


# VST normalization

#Blind = TRUE only for quality check
#VST model aware
vsd <- varianceStabilizingTransformation(dds.time.RNA.peaks, blind=F)
vst.time.RNA.peaks <- assay(vsd)
vst.time.RNA.peaks <- vst.time.RNA.peaks - min(vst.time.RNA.peaks)

vst_RNA<-as.data.frame(vst.time.RNA.peaks)
colnames(vst_RNA)
length(vst_RNA$Day0_rep1)
vst_RNA$Histone <- all_raw_counts_clean$Histone


#Save DBA.object
saveRDS(vst_RNA, file = "~/R_objects/vst_RNA_27ac_report_peaks_Sept2021.rds")

#VST preliminar evaulation
meanSdPlot(vst.time.RNA.peaks)
