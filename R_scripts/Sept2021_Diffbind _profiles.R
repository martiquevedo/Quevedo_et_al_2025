
# Setup -------------------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))
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
suppressPackageStartupMessages(library(reshape2))

setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq")

##PlotEigene #

plotEigengene <- function(data, genes, condition, time, timeUnits = "Time",
                          inverse = F, title = "", noGrid = T, colors = c("deepskyblue", "darkorange2","darkgreen","red3"),
                          noLegend=T) {
  
  require(ggplot2)
  require(RColorBrewer)
  require(reshape2)
  
  #require(plotly)
  
  expr <- NA
  
  #d <- data[,which(colnames(data) %in% gene)]
  d <- data[, genes]
  
  if (length(genes) == 1 ) { #1 gene, directly the data
    expr <- scale(d)
  } else if (length(genes > 1)) { #cluster, eigengene
    pca <- prcomp((d))
    pc1 <- pca$x[, 1]
    
    #print(percents <- round(summary(pca)$importance[2,]*100))
    
    if(sum(sign(cor(pca$x[,1,drop = FALSE], d))) < 0) {
      pc1 <- pc1 * -1
    }
    if(inverse)
      pc1 <- pc1 * -1
    expr <- pc1
  } else { # 1 gene is needed
    stop("this function needs at least one gene")
  }
  
  myplot <- ggplot(data.frame(x = time, y = scale(expr), g = condition),
                   
                   aes(x = x, y = y, group = g)) +
    stat_summary(fun.data = mean_se, geom = "ribbon", fill = "lightgrey", alpha = 0.75) +
    stat_summary(fun.data = mean_se, geom = "line", aes(col = g), lwd = 2) + #      plot_output_list <- lapply(shiftedColors, function(color) {
    xlab(timeUnits) +
    ylab("z-score") +
    labs(color = "Condition") +
    ggtitle(title)
  
  if(noGrid) {
    myplot <- myplot +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
  }
  
  if(is.vector(colors) & length(colors) >= length(unique(condition))) {
    myplot <- myplot + scale_color_manual(values=colors)
  }
  
  if(noLegend) {
    myplot <- myplot +
      theme(legend.position = "none")
  }
  return (myplot)
}



#Multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# Approach ----------------------------------------------------------------

#1 Extract reports from Diffbind specifying ALL peaks, with th=1
#2 Construct supermatrix with Day0 to 1 to 4 to 7 transitions (FDR and Fold from reports) + coverage mean values of each day for plotting (apply H3)
#
#3 Set profiles as objects (test them in single plotting)
#4 Create function to feed the profiles and store regions and plots.
#5 Implemented the plotting to also plot the other histones DONE! 
#6 histone normalization by Day0


# Diffbind pre-analysis 
### DIFFBIND COUNTS AND REPORT
 meta <- read.csv("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/old/Meta_histones.csv", header = T, sep = ";")
 colnames(meta) <- c("X","SampleID","Histone","Time","Replicate","Group1")
 meta.k27ac <- meta[grep("h3k27ac",meta$Histone),]
 meta.k27ac$Time <- as.integer(substr(meta.k27ac$Time, 4,4))

# Using pre-saved diffbind object  
## IMPORTANT! SET WHICH HISTONE TO ANALYSE AND APPLY COLOR

 
 i="27ac"
 
 
    # load Diffbind object
     dba_object <- readRDS(file = paste0("~/R_objects/Sept2021_k",i,"_diffbind.rds"))
     
     stats_overview <- dba.show(dba_object, bContrast=T)
     stats_overview
 
    # th=1 gives back all peaks, bCounts the signal (apply mean, VST convertible?)
    report1 <- dba.report(dba_object,th=1, contrast=1,bCounts=T, bCalledDetail=T, DataType=DBA_DATA_FRAME)
    report2 <- dba.report(dba_object,th=1, contrast=4,bCounts=T, bCalledDetail=T, DataType=DBA_DATA_FRAME)
    report3 <- dba.report(dba_object,th=1, contrast=6,bCounts=T, bCalledDetail=T, DataType=DBA_DATA_FRAME) #in new implementation is changed in order, check later it makes sense in superMatrix
    report4 <- dba.report(dba_object,th=1, contrast=3,bCounts=T, bCalledDetail=T, DataType=DBA_DATA_FRAME)
    report5 <- dba.report(dba_object,th=1, contrast=5,bCounts=T, bCalledDetail=T, DataType=DBA_DATA_FRAME)
    report6 <- dba.report(dba_object,th=1, contrast=2,bCounts=T, bCalledDetail=T, DataType=DBA_DATA_FRAME)
    
    
    #MERGE TABLE by first 3 columns
    #Objective Day0-counts(mean) FDR0-1 Day1-counts(mean) FDR1-4 Day4-counts(mean) FDR4-7 Day7-counts(mean)
    # 1. sort, see if equal, just construct data.frame manually?
    report1 <- report1[order(report1$Chr,report1$Start),]
    report2 <- report2[order(report2$Chr,report2$Start),]
    report3 <- report3[order(report3$Chr,report3$Start),]
    report4 <- report4[order(report4$Chr,report4$Start),]
    report5 <- report5[order(report5$Chr,report5$Start),]
    report6 <- report5[order(report6$Chr,report6$Start),]
    
   
    ## EXPORT this coordinates for doing the bedtools_coverage at this sides (both RNA and chip BAMs)
    write.table(report1[1:3], file = "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Sept2021_extra_analysis/profiles/27ac_report.bed", append = FALSE, sep="\t", dec = ".",
                row.names = FALSE, col.names = F, quote = F)
    
    
    
    setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Sept2021_extra_analysis/profiles")
    
    #READ coverage (from bedtools multicov, bams/bais using report as coordinates)
    coverage_27ac_peaks <- read.table("k27ac_peaks_report_allhistones_coverage.txt", stringsAsFactors = F)
    coverage_names <- read.table("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/bams_bwa/16M/all_bams/all_links/bam_list.txt", stringsAsFactors = F)
    colnames(coverage_27ac_peaks) <- c("Chr","start","end","na",coverage_names$V1)
    
    #REMOVE Chr and converting to numeric to avoid ""
    coverage_27ac_peaks [,1] <- gsub('Chr', '', coverage_27ac_peaks [,1])
    coverage_27ac_peaks [,1] <- as.numeric(as.character(coverage_27ac_peaks [,1]))
    
    coverage_27ac_peaks_clean <- coverage_27ac_peaks [,-c(1:4)] 
    coverage_27ac_peaks_coordinates <- paste(coverage_27ac_peaks [,1], ':', coverage_27ac_peaks [,2], ':', coverage_27ac_peaks [,3], sep = "")
    rownames(coverage_27ac_peaks_clean)<-coverage_27ac_peaks_coordinates
    
    #construct table for all histones (loop possible by cycling through day, histone and rep)
    coverage_27ac_normalized <- coverage_27ac_peaks_clean[,FALSE]
    coverage_27ac_normalized$day7_h3k27ac_rep1 <- coverage_27ac_peaks_clean$day7_k27ac_rep1_16M.sorted.bam / coverage_27ac_peaks_clean$day7_H3_16M.sorted.bam
    coverage_27ac_normalized$day7_h3k27ac_rep2 <- coverage_27ac_peaks_clean$day7_k27ac_rep2_16M.sorted.bam / coverage_27ac_peaks_clean$day7_H3_16M.sorted.bam
    coverage_27ac_normalized$day7_h3k27me3_rep1 <- coverage_27ac_peaks_clean$day7_k27me3_rep1_16M.sorted.bam / coverage_27ac_peaks_clean$day7_H3_16M.sorted.bam
    coverage_27ac_normalized$day7_h3k27me3_rep2 <- coverage_27ac_peaks_clean$day7_k27me3_rep2_16M.sorted.bam / coverage_27ac_peaks_clean$day7_H3_16M.sorted.bam
    coverage_27ac_normalized$day7_h3k4me3_rep1 <- coverage_27ac_peaks_clean$day7_k4me3_rep1_16M.sorted.bam / coverage_27ac_peaks_clean$day7_H3_16M.sorted.bam
    coverage_27ac_normalized$day7_h3k4me3_rep2 <- coverage_27ac_peaks_clean$day7_k4me3_rep2_16M.sorted.bam / coverage_27ac_peaks_clean$day7_H3_16M.sorted.bam
    coverage_27ac_normalized$day7_h3k9me2_rep1 <- coverage_27ac_peaks_clean$day7_k9me2_rep1_16M.sorted.bam / coverage_27ac_peaks_clean$day7_H3_16M.sorted.bam
    coverage_27ac_normalized$day7_h3k9me2_rep2 <- coverage_27ac_peaks_clean$day7_k9me2_rep2_16M.sorted.bam / coverage_27ac_peaks_clean$day7_H3_16M.sorted.bam
    
    
    coverage_27ac_normalized_vst <- coverage_27ac_normalized
    coverage_27ac_normalized_vst[is.na(coverage_27ac_normalized_vst)] <- 0
    coverage_27ac_normalized_vst[coverage_27ac_normalized_vst == Inf] <- 0    
    coverage_27ac_normalized_vst <- sapply( coverage_27ac_normalized_vst, as.integer )
    coverage_27ac_normalized_vst <- as.matrix(coverage_27ac_normalized_vst)+1 #technical patch to avoid 0s in vsd, which provokes errors
    
    #DESEQ
    meta_all <- read.csv("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/Meta_histones.csv", header = T, sep = ";")
    colnames(meta_all) <- c("X","SampleID","Histone","Time","Replicate","Group1")
    
    dds.time.27ac <- DESeqDataSetFromMatrix(
      countData = round(coverage_27ac_normalized_vst),
      colData = meta_all,
      design = ~ Time)
    
    #Blind = TRUE only for quality check
    #VST model aware
    
    vsd_dds.time.27ac <- varianceStabilizingTransformation(dds.time.27ac, blind=F)
    vst_dds.time.27ac <- assay(vsd_dds.time.27ac)
    vst_dds.time.27ac <- vst_dds.time.27ac - min(vst_dds.time.27ac)
    vst_forplot <- vst_dds.time.27ac # all histone at 27ac report peaks
    rownames(vst_forplot)<-coverage_27ac_peaks_coordinates #to be used for ploteugenie
    
    
    # MADE a VST of 
    # read colnames from "/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/bams_bwa/16M/all_bams/all_links/bam_list.txt"
    vst_rna <- readRDS(file = paste0("~/R_objects/vst_RNA_27ac_report_peaks_Sept2021.rds"))
    
    
  
    # Combine everything, REMEMBER DAY4toDay7 is INVERSED... 
    supermatrix <- data.frame(Chr=report1$Chr,Start=report1$Start, End=report1$End, day0=rowMeans(vst_forplot[,1:2], na.rm=TRUE), day0to1.FDR=report1$FDR, day0to1.Fold=report1$Fold,
                                    day1=rowMeans(vst_forplot[,9:10], na.rm=TRUE), day1to4.FDR=report2$FDR, day1to4.Fold=report2$Fold,
                                    day4=rowMeans(vst_forplot[,17:18], na.rm=TRUE), day4to7.FDR=report3$FDR, day4to7.Fold=report3$Fold,
                                    day7=rowMeans(vst_forplot[,25:26], na.rm=TRUE), day0to7.FDR=report4$FDR, day0to7.Fold=report4$Fold, 
                                    day0to4.FDR=report5$FDR, day0to4.Fold=report5$Fold, day1to7.FDR=report5$FDR, day1to7.Fold=report5$Fold, region.lenght=(report1$End-report1$Start))
    
    rownames(supermatrix)<-coverage_27ac_peaks_coordinates
    rownames(supermatrix) 
    colnames(supermatrix) 
    

    
    #Create list of profiles
    #REMEMBER THAT FOLD IS CALCULATED FROM 1st to 2nd timepoint (negative!!! Fold)
    # FOLD is loq2
    # Remember that contrast 4to7 is INVERSED
    
    #UP-DOWN-DOWN-DOWN
    profile01 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold > 1 & between(supermatrix$day1to4.Fold,-0.5,0.5) & between(supermatrix$day1to7.Fold,-0.5,0.5) & between(supermatrix$day4to7.Fold,-0.5,0.5) 
    #DOWN-UP-UP-UP
    profile02 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold < -1 & between(supermatrix$day1to4.Fold,-0.5,0.5) & between(supermatrix$day1to7.Fold,-0.5,0.5)  & between(supermatrix$day4to7.Fold,-0.5,0.5) 
    #DOWN-UP-UP-DOWN
    profile03 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold < -1 &  between(supermatrix$day0to7.Fold,-0.5,0.5) &  between(supermatrix$day1to4.Fold,-0.5,0.5) 
    #UP-DOWN-DOWN-UP
    profile04 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold > 1 &  between(supermatrix$day0to7.Fold,-0.5,0.5) &  between(supermatrix$day1to4.Fold,-0.5,0.5) 
    #DOWN-DOWN-DOWN-UP
    profile05 <- supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold > 1 & between(supermatrix$day0to1.Fold,-0.5,0.5) &  between(supermatrix$day1to4.Fold,-0.5,0.5) & supermatrix$day0to7.FDR < 0.05 & supermatrix$day0to7.Fold < -1
    #UP-UP-UP-DOWN
    profile06 <- supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold < -1 &  supermatrix$day1to7.Fold > 1 & between(supermatrix$day0to1.Fold,-0.5,0.5) &  between(supermatrix$day1to4.Fold,-0.5,0.5) 
    #DOWN-DOWN-UP-DOWN
    profile07 <- supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold < -1 &  between(supermatrix$day0to1.Fold,-0.5,0.5) & supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold < -1 &  between(supermatrix$day0to7.Fold,-0.5,0.5)
    #UP-UP-DOWN-UP
    profile08 <- supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold > 1 &  between(supermatrix$day0to1.Fold,-0.5,0.5) & supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold > 1 &  between(supermatrix$day0to7.Fold,-0.5,0.5)
    #DOWN-UP-DOWN-DOWN
    profile09 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold < -1 &  supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold > 1 &  between(supermatrix$day0to7.Fold,-1,1) 
    #UP-DOWN-UP-UP
    profile10 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold > 1 &  supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold < -1 &  between(supermatrix$day0to7.Fold,-1,1) 
    # 
    # profile11 <- supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold < -1 &  supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold > 1 &  between(supermatrix$day0to7.Fold,-1,1) 
    # profile12 <- supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold > 1 &  supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold < -1 &  between(supermatrix$day0to7.Fold,-1,1) & supermatrix$day0to1.FDR > 0.05 &  between(supermatrix$day1to7.Fold,-1,1) &  between(supermatrix$day0to1.Fold,-1,1) 
    # profile13 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold < -1 &  supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold > 1 &  between(supermatrix$day0to4.Fold,-1,1) & supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold < -1
    # profile14 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold > 1 &  supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold < -1 &  between(supermatrix$day0to4.Fold,-1,1) & supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold > 1
    
    
    profileList <- list(profile01,profile02,profile03,profile04,profile05,profile06,profile07,profile08,profile09,profile10)
    names(profileList) <- sprintf("profile%02d", 1:10)
    
    #SUM(profile) and remove the ones with 0!!

    
    subseto_forplot <- melt(select(supermatrix[profile02,] , day0, day1, day4,day7))
    colnames(subseto_forplot)<- c("Days","Value" )
    subseto_forplot$Days<-as.integer(substr(subseto_forplot$Days, 4,4))
    
    
    
    subseto_plot <- ggplot(subseto_forplot, aes(x=Days,y=Value, group=1)) +
      stat_summary(fun.data = mean_se, geom = "ribbon", fill = "lightgrey", alpha = 0.75) +
      stat_summary(fun.data = mean_se, geom = "line", colour = "blue", lwd = 2)+
      ggtitle(paste("profile01","\n#",(length(subseto_forplot$Value)/4), sep = ""))
    subseto_plot + scale_x_continuous(name="Days of light", breaks=c(0,1, 4, 7), limits=c(0, 7))+
      theme(plot.title = element_text(size=10, hjust=0.5))+ scale_y_continuous(name="z score")
    
    
    
    #CREATE LOOP to process all profiles

    # Individual histone
    profileList <- list(profile01,profile02,profile03,profile04,profile05,profile06,profile07,profile08,profile09,profile10,profile11,profile12,profile13,profile14)
    names(profileList) <- sprintf("profile%02d", 1:14)
  
 
    profile_plotting <- lapply (1:length(profileList), function(x) {
        elemento <- profileList[[x]]
        nombre <- names(profileList)[x]
    
      #plot
      subseto_forplot <- melt(select(supermatrix[elemento,] , day0, day1, day4,day7))
      colnames(subseto_forplot)<- c("Days","Value" )
      subseto_forplot$Days<-as.integer(substr(subseto_forplot$Days, 4,4))



      subseto_plot <- ggplot(subseto_forplot, aes(x=Days,y=Value, group=1)) +
        stat_summary(fun.data = mean_se, geom = "ribbon", fill = "lightgrey", alpha = 0.75) +
        stat_summary(fun.data = mean_se, geom = "line", colour = "darkgreen", lwd = 2)+
        ggtitle(paste(nombre,"\n#",(length(subseto_forplot$Value)/4), sep = ""))
      subseto_plot + scale_x_continuous(name="Days of light", breaks=c(0,1, 4, 7), limits=c(0, 7))+
        theme(plot.title = element_text(size=10, hjust=0.5))+ scale_y_continuous(name="z score")
      })

    
    

    #RNA DATA
    supermatrix$day0.RNA <- rowMeans(vst_rna[,1:3], na.rm=TRUE)
    supermatrix$day1.RNA <- rowMeans(vst_rna[,4:6], na.rm=TRUE)
    supermatrix$day4.RNA <- rowMeans(vst_rna[,7:9], na.rm=TRUE)
    supermatrix$day7.RNA <- rowMeans(vst_rna[,10:12], na.rm=TRUE)
    
      profile_RNA_plotting <- lapply (1:length(profileList), function(x) {
        elemento <- profileList[[x]]
        nombre <- names(profileList)[x]
        #plot
        subseto_forplot <- melt(select(supermatrix[elemento,] , day0.RNA, day1.RNA, day4.RNA,day7.RNA))
        colnames(subseto_forplot)<- c("Days","Value" )
        subseto_forplot$Days<-as.integer(substr(subseto_forplot$Days, 4,4))
      
      
      
        subseto_plot <- ggplot(subseto_forplot, aes(x=Days,y=Value, group=1)) +
          stat_summary(fun.data = mean_se, geom = "ribbon", fill = "lightgrey", alpha = 0.75) +
          stat_summary(fun.data = mean_se, geom = "line", colour = "black", lwd = 2)+
          ggtitle(paste(" ",nombre, "\n RNAseq", sep = ""))
        subseto_plot + scale_x_continuous(name="days of light", breaks=c(0,1, 4, 7), limits=c(0, 7))+
          theme(plot.title = element_text(size=10, hjust=0.5), text=element_text(size=10))+ scale_y_continuous(name="z score")
        })
    
   
    
       paired_graphs_forplot <- list()
       for(x in 1:length(profileList)) {
         paired_graphs_forplot <- append(paired_graphs_forplot,profile_plotting[x], after = length(paired_graphs_forplot))
         paired_graphs_forplot <- append(paired_graphs_forplot,profile_RNA_plotting[x], after = length(paired_graphs_forplot))
      #   paired_graphs_forplot <- append(paired_graphs_forplot,profile_allmC_plotting[x], after = length(paired_graphs_forplot))
      #   paired_graphs_forplot
         }
        
     multiplot(plotlist = paired_graphs_forplot,  layout = matrix(1:28, nrow = 7, ncol = 4, byrow=TRUE))
    #   
    # png(paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/profile_approach/plots/",i,"_profiles_with_RNA_NEW.png"), width = 1800, height = 1800)
    # multiplot(plotlist = paired_graphs_forplot,  layout = matrix(1:28, nrow = 7, ncol = 4, byrow=TRUE))
    # dev.off()
    
      ############### 
      
      
      colnames(vst_forplot)
      vst_normalize <- vst_forplot
      
      rownames(vst_forplot)
      meta.histones <- read.csv("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/Meta_histones.csv", header = T, sep = ";")
      colnames(meta.histones) <- c("X","SampleID","Histone","Time","Replicate","Group1")
      meta.histones$Time <- as.integer(substr(meta.histones$Time, 4,4))
      
      #normalizing, I used DAY0 as normalizer
      
      for (z in c(0,2,4,6))  {
        vst_normalize[,c(1,2,9,10,17,18)+z] <- (vst_normalize[,c(1,2,9,10,17,18)+z] / rowMeans(vst_normalize[,c(1,2)+z]))
      }
      
############ PlotEigengene2 modified no scale and lines
      
      plotEigengene2 <- function(data, genes, condition, time, timeUnits = "Time",
                                 inverse = F, title = "", noGrid = T, colors = c("deepskyblue", "darkorange2","darkgreen","red3"),
                                 noLegend=T) {
        
        require(ggplot2)
        require(RColorBrewer)
        require(reshape2)
        
        #require(plotly)
        
        expr <- NA
        
        #d <- data[,which(colnames(data) %in% gene)]
        d <- data[, genes]
        
        if (length(genes) == 1 ) { #1 gene, directly the data
          expr <- scale(d)
        } else if (length(genes > 1)) { #cluster, eigengene
          pca <- prcomp((d))
          pc1 <- pca$x[, 1]
          
          #print(percents <- round(summary(pca)$importance[2,]*100))
          
          if(sum(sign(cor(pca$x[,1,drop = FALSE], d))) < 0) {
            pc1 <- pc1 * -1
          }
          if(inverse)
            pc1 <- pc1 * -1
          expr <- pc1
        } else { # 1 gene is needed
          stop("this function needs at least one gene")
        }
        
        myplot <- ggplot(data.frame(x = time, y = expr, g = condition),
                         
                         aes(x = x, y = y, group = g)) +
          stat_summary(fun.data = mean_se, geom = "ribbon", fill = "lightgrey", alpha = 0.75) +
          stat_summary(fun.data = mean_se, geom = "line", aes(col = g),linetype = c(5,5,5,5,5,5,5,5,1,1,1,1,5,5,5,5), lwd = 2) + #      plot_output_list <- lapply(shiftedColors, function(color) {
          xlab(timeUnits) +
          ylab("normalized signal") +
          labs(color = "Histone") +
          ggtitle(title)
        
        if(noGrid) {
          myplot <- myplot +
            theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank())
        }
        
        if(is.vector(colors) & length(colors) >= length(unique(condition))) {
          myplot <- myplot + scale_color_manual(values=colors)
        }
        
        if(noLegend) {
          myplot <- myplot +
            theme(legend.position = "none")
        }
        return (myplot)
      }
      
      ###########
      
      df2 <- t(vst_normalize)
 
      multigraphs_forplot <- lapply (1:length(profileList), function(x) {
        elemento <- profileList[[x]]
        nombre <- names(profileList)[x]
        regions <- colnames(df2)[profileList[[x]]]
        number_regions <- length(regions)
        ploty <- ggplot()
        
        if (number_regions == 0) {
          ploty <-ploty + geom_blank()
        } else {  
        title_graph <- paste(nombre,"\n ",i," #",number_regions, sep = "")
        ploty <- plotEigengene2(df2, regions, meta.histones$Histone, meta.histones$Time, title = title_graph, noLegend=F)
        ploty <- ploty + scale_x_continuous(name="days of light", breaks=c(0,1, 4, 7), limits=c(0, 7)) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        theme(legend.position = "none", plot.title = element_text(size=10, hjust=0.5),text=element_text(size=10))
        ploty
        
        }
        
      })
      

      paired_multigraphs <- list()
      
      for(x in 1:length(multigraphs_forplot)) {
        paired_multigraphs <- append(paired_multigraphs,multigraphs_forplot[x], after = length(paired_multigraphs))
        paired_multigraphs <- append(paired_multigraphs,profile_RNA_plotting[x], after = length(paired_multigraphs))
        
        paired_multigraphs
      }
      
      multiplot(plotlist = paired_multigraphs[1:20],  layout = matrix(1:20, nrow = 5, ncol = 4, byrow=TRUE))
      
      pdf("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/plots/profiles_27ac_multi_RNA.pdf")
      multiplot(plotlist = paired_multigraphs[1:20],  layout = matrix(1:20, nrow = 5, ncol = 4, byrow=TRUE))
      dev.off()

    ###########################
    
    
    # OUTPUT TABLES -----------------------------------------------------------
    
    
    profile_data <- lapply (1:length(profileList), function(x) {
                    elemento <- profileList[[x]]
                    nombre <- names(profileList)[x]
                    subset_table <- supermatrix[elemento,]
                    subset_table
                    
                })
    names(profile_data) <- sprintf("profile%02d", 1:10)
    
    for(x in 1:length(profile_data)) {
        elemento <- profile_data[[x]]
        nombre <- names(profileList)[x]
        
        # FOR HOMER
        out="~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/profiles/beds/"
        #REMOVE Chr and converting to numeric to avoid ""
        elemento[,1] <- gsub('Chr', '', elemento[,1])
        elemento[,1] <- as.numeric(as.character(elemento[,1]))
        write.table(elemento[1:3], file = c(paste0(out,nombre,".bed")), append = FALSE, sep="\t", dec = ".",
                    row.names = FALSE, col.names = F)
        
        
        
    }

    
    
    #### FOR GENE IDs
    
    ### SHELL SCRIPT
    
    module load bioinfo-tools BEDTools
    
    in=~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/profiles/beds/
      out=~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/profiles/GeneIDs/
      
      for f in $(find $in -maxdepth 1 -name "*.bed"); do
    fnam=$(basename ${f/.bed/})
    output=$out$fnam"_vs_araport.bed"
    
    echo $f
    echo $output
    
    bedtools intersect -wa -wb -a /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/araport11_complete_NEW.bed \
    -b $f > $output
    
    done
    
    ## then extract IDs from file / Rscript
    
    in_dir<-"~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/profiles/GeneIDs/"
    setwd("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Sept2021/profiles/GeneIDs/vsAraport/")
    
    intersected_IDs <-  dir(".", pattern="araport", recursive = F, full.names = TRUE) 
    
    Gene_IDs <- lapply (intersected_IDs, function(x){
      table <- read.table(x, header = F, sep="\t")
      geneIDs <- table$V6
      geneIDs
    })
    
    names(Gene_IDs)<-intersected_IDs
    names(Gene_IDs) <- strtrim(names(Gene_IDs), 11)
    names(Gene_IDs) <- substring(names(Gene_IDs),3)
    
    for(x in 1:length(Gene_IDs)) {
      elemento <- Gene_IDs[[x]]
      nombre <- names(Gene_IDs)[x]
      
      # Write table of gene IDs
      write.table(elemento, file = paste0(in_dir,nombre,"_GeneIDs.txt"), sep ="\t", row.names = F, col.names = F, quote=FALSE )
      
    }
    
    
   