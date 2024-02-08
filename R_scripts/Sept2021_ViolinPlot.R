library(tidyverse)
library(ggplot2)
vst_rna_4me3 <- readRDS(file = "~/R_objects/vst_RNA_4me3_Sept2021.rds")
vst_rna_27ac <- readRDS(file = "~/R_objects/vst_RNA_27ac_Sept2021.rds")
vst_rna_27me3 <- readRDS(file = "~/R_objects/vst_RNA_27me3_Sept2021.rds")
vst_rna_9me2 <- readRDS(file = "~/R_objects/vst_RNA_9me2_Sept2021.rds")

day7_4me3 <- data.frame(stringsAsFactors = FALSE, D7_4me3 = rowMeans(vst_rna_4me3[,10:12], na.rm=TRUE))
rownames(day7_4me3)<-NULL
day7_27ac <- data.frame(stringsAsFactors = FALSE, D7_27ac = rowMeans(vst_rna_27ac[,10:12], na.rm=TRUE))
rownames(day7_27ac)<-NULL
day7_27me3 <- data.frame(stringsAsFactors = FALSE, D7_27me3 = rowMeans(vst_rna_27me3[,10:12], na.rm=TRUE))
rownames(day7_27me3)<-NULL
day7_9me2 <- data.frame(stringsAsFactors = FALSE, D7_9me2 = rowMeans(vst_rna_9me2[,10:12], na.rm=TRUE))
rownames(day7_9me2)<-NULL


d7_df <- day7_27ac %>% rownames_to_column() %>% 
    left_join(day7_27me3 %>% rownames_to_column()) %>% 
    left_join(day7_4me3 %>% rownames_to_column()) %>% 
    left_join(day7_9me2 %>% rownames_to_column()) %>%
    select(-rowname)

# boxplot(d7_df,
#         main = "RNA-seq reads at Peaks",
#         names = c("H3K27ac", "H3K27me3", "H3K4me3", "H3K9me2"),
#         col =  c("#4169e1", "#ee7600", "#663399", "#cd0000"),
#         border = "black",
#         horizontal = F,
#         outline=FALSE,
#         notch = TRUE
# )

histone_colors <- c("#4169e1", "#ee7600", "#663399", "#cd0000")
#Convert to ggplot (remember that melt only accepts matrix)
d2.df <- reshape2::melt(as.matrix(d7_df), c("x", "y"), value.name = "z")
head(d2.df)
d2.df$y <- as.factor(d2.df$y) ## Convert the variable from a numeric to a factor variable

p <- ggplot(d2.df, aes(x=y, y=z, fill=y)) + 
    geom_violin(trim=FALSE)
   
class(d2.df$y)

p + geom_boxplot(outlier.shape = NA, width=0.1, fill="white") +
    #scale_fill_manual(values=c("#4169e1", "#663399", "#FF9900", "#CC0000")) +
    scale_fill_manual(values=histone_colors)+
    ylim(-2,12)+
    theme_classic()

pdf(file="/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/RNAseq/violin_plots/day7_RNA_zscore.pdf",
    useDingbats=FALSE)
p + geom_boxplot(outlier.shape = NA, width=0.1, fill="white") +
    #scale_fill_manual(values=c("#4169e1", "#663399", "#FF9900", "#CC0000")) +
    scale_fill_manual(values=histone_colors)+
    ylim(-2,12)+
    theme_classic()
dev.off()

