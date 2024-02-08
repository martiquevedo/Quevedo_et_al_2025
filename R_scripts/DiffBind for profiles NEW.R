
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

#1 Extract reports from Diffbind especifying ALL peaks, with th=1
#2 Construct supermatrix with Day0 to 1 to 4 to 7 transitions (FDR and Fold from reports) + VST mean values of each day for plotting (vst.k27ac.time) 
## RUN VST prior this to have vst.k27ac.time object or load the object
#3 Set profiles as objects (biological relevants, or all?)
#4 Create function to feed the profiles and store regions and plots.
##5# Implement the plotting to also plot the other histones DONE! 
  # TODO, normalizing better? plotting RNA in the same graph?

#####REDO vst.all OBJECTS!!!


# Diffbind pre-analysis (not active) ---------------------------------------------------
### DIFFBIND COUNTS AND REPORT
# meta <- read.csv("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/Meta_histones.csv", header = T, sep = ";")
# colnames(meta) <- c("X","SampleID","Histone","Time","Replicate","Group1")
# meta.k27ac <- meta[grep("h3k27ac",meta$Histone),]
# meta.k27ac$Time <- as.integer(substr(meta.k27ac$Time, 4,4))
# 
# greening_27ac <- dba(sampleSheet="~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/k27ac_sheet_diffbind3.csv")
# greening_27ac <- dba.count(greening_27ac)
# greening_27ac <- dba.contrast(greening_27ac,categories = DBA_CONDITION, minMembers=2)
# greening_27ac <- dba.analyze(greening_27ac)





# Using pre-saved objects -------------------------------------------------

### IMPORTANT! SET WHICH HISTONE TO ANALYSE AND APPLY COLOR

i="27ac"


    #load Diffbind object
    dba_object <- readRDS(file = paste0("~/R_objects/greening_",i,"_NEW.rds"))
    
    
    
    # th=1 gives back all peaks, bCounts the signal (apply mean, VST convertible?)
    report1 <- dba.report(dba_object,th=1, contrast=1,bCounts=T, bCalledDetail=T, DataType=DBA_DATA_FRAME)
    report2 <- dba.report(dba_object,th=1, contrast=4,bCounts=T, bCalledDetail=T, DataType=DBA_DATA_FRAME)
    report3 <- dba.report(dba_object,th=1, contrast=6,bCounts=T, bCalledDetail=T, DataType=DBA_DATA_FRAME)
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
    
    # #load VST object
    vst_forplot <- readRDS(file = paste0("~/R_objects/vst_dds.time.",i,"_NEW.rds"))
    vst_rna <- readRDS(file = paste0("~/R_objects/vst_RNA_normalized_NEW_",i,".rds"))
    # vst_rna <- vst_rna[which(vst_rna$Histone==paste0("k",i)),]
    # 
    # 
    # vst_allmC <- readRDS(file = "~/R_objects/vst_allmC.rds")
    # vst_allmC <- vst_allmC[which(vst_allmC$Histone==paste0("k",i)),]
    
   
    
    

    #no need to transpose the saved objects
    
    #By luck when combining vst and the report the rownames assume the vst 
    supermatrix <- data.frame(Chr=report1$Chr,Start=report1$Start, End=report1$End, day0=rowMeans(vst_forplot[,1:2], na.rm=TRUE), day0to1.FDR=report1$FDR, day0to1.Fold=report1$Fold,
                                    day1=rowMeans(vst_forplot[,3:4], na.rm=TRUE), day1to4.FDR=report2$FDR, day1to4.Fold=report2$Fold,
                                    day4=rowMeans(vst_forplot[,5:6], na.rm=TRUE), day4to7.FDR=report3$FDR, day4to7.Fold=report3$Fold,
                                    day7=rowMeans(vst_forplot[,7:8], na.rm=TRUE), day0to7.FDR=report4$FDR, day0to7.Fold=report4$Fold, 
                                    day0to4.FDR=report5$FDR, day0to4.Fold=report5$Fold, day1to7.FDR=report5$FDR, day1to7.Fold=report5$Fold, region.lenght=(report1$End-report1$Start))
    
    
    #Create list of profiles
    #REMEMBER THAT FOLD IS CALCULATED FROM 1st to 2nd timepoint (negative!!! Fold)
    
    profile01 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold > 1.5 &  between(supermatrix$day4to7.Fold,-1.5,1.5) 
    profile02 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold < -1.5 &  supermatrix$day0to7.Fold < -1.5 &  between(supermatrix$day4to7.Fold,-1.5,1.5) &  between(supermatrix$day1to4.Fold,-1.5,1.5) 
    profile03 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold < -1.5 &  between(supermatrix$day0to7.Fold,-1.5,1.5) &  between(supermatrix$day1to4.Fold,-1,1) 
    profile04 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold > 1.5 &  between(supermatrix$day0to7.Fold,-1.5,1.5) &  between(supermatrix$day1to4.Fold,-1.5,1.5) 
    profile05 <- supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold < -1.5 &  supermatrix$day0to7.Fold < -1.5  &  between(supermatrix$day0to1.Fold,-1.5,1.5) &  between(supermatrix$day1to4.Fold,-1.5,1.5) 
    profile06 <- supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold > 1.5 &  supermatrix$day1to7.Fold > 1.5 &  supermatrix$day0to7.Fold > 1.5 & between(supermatrix$day0to1.Fold,-1.5,1.5) &  between(supermatrix$day1to4.Fold,-1.5,1.5) 
    profile07 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold > 1.5 &  between(supermatrix$day1to4.Fold,-1.5,1.5)  &  supermatrix$day1to7.FDR < 0.05 & supermatrix$day1to7.Fold > 1.5 &  supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold > 1.5 & between(supermatrix$day1to4.Fold,-1.5,1.5) 
    profile08 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold < -1.5 &  between(supermatrix$day1to4.Fold,-1.5,1.5)  &  supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold < -1.5 
    profile09 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold < -1.5 &  supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold > 1.5 &  between(supermatrix$day0to7.Fold,-1,1) 
    profile10 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold > 1.5 &  supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold < -1.5 &  between(supermatrix$day0to7.Fold,-1,1) 
    profile11 <- supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold < -1.5 &  supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold > 1.5 &  between(supermatrix$day0to7.Fold,-1,1) 
    profile12 <- supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold > 1.5 &  supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold < -1.5 &  between(supermatrix$day0to7.Fold,-1,1) & supermatrix$day0to1.FDR > 0.05 &  between(supermatrix$day1to7.Fold,-1,1) &  between(supermatrix$day0to1.Fold,-1,1) 
    profile13 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold < -1.5 &  supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold > 1.5 &  between(supermatrix$day0to4.Fold,-1.5,1.5) & supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold < -1.5
    profile14 <- supermatrix$day0to1.FDR < 0.05 & supermatrix$day0to1.Fold > 1.5 &  supermatrix$day1to4.FDR < 0.05 & supermatrix$day1to4.Fold < -1.5 &  between(supermatrix$day0to4.Fold,-1.5,1.5) & supermatrix$day4to7.FDR < 0.05 & supermatrix$day4to7.Fold > 1.5
    
    
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

    
    

    #for RNA
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
    
    
    # #for methylation
    # supermatrix$day0.allmC <- rep(0, length(supermatrix$day0))
    # supermatrix$day1.allmC <- rep(0, length(supermatrix$day0))
    # supermatrix$day4.allmC <- rowMeans(vst_allmC[,1:3], na.rm=TRUE)
    # supermatrix$day7.allmC <- rowMeans(vst_allmC[,4:6], na.rm=TRUE)
    #   
    # profile_allmC_plotting <- lapply (1:length(profileList), function(x) {
    #   elemento <- profileList[[x]]
    #   nombre <- names(profileList)[x]
    #   #plot
    #   subseto_forplot <- melt(select(supermatrix[elemento,] ,  day4.allmC,day7.allmC))
    #   colnames(subseto_forplot)<- c("Days","Value" )
    #   subseto_forplot$Days<-as.integer(substr(subseto_forplot$Days, 4,4))
    #   
    #   
    #   
    #   subseto_plot <- ggplot(subseto_forplot, aes(x=Days,y=Value, group=1)) +
    #     stat_summary(fun.data = mean_se, geom = "ribbon", fill = "lightgrey", alpha = 0.75) +
    #     stat_summary(fun.data = mean_se, geom = "line", colour = "purple", lwd = 2)+
    #     ggtitle(paste(" ",nombre, "\n CpG and CWG", sep = ""))
    #   subseto_plot + scale_x_continuous(name="days of light", breaks=c(0,1, 4, 7), limits=c(0, 7))+
    #     theme(plot.title = element_text(size=10, hjust=0.5), text=element_text(size=10))+ scale_y_continuous(name="z score")
    # })
    
     
    

    
       paired_graphs_forplot <- list()
       for(x in 1:length(profileList)) {
         paired_graphs_forplot <- append(paired_graphs_forplot,profile_plotting[x], after = length(paired_graphs_forplot))
         paired_graphs_forplot <- append(paired_graphs_forplot,profile_RNA_plotting[x], after = length(paired_graphs_forplot))
      #   paired_graphs_forplot <- append(paired_graphs_forplot,profile_allmC_plotting[x], after = length(paired_graphs_forplot))
      #   paired_graphs_forplot
         }
        
    # multiplot(plotlist = paired_graphs_forplot,  layout = matrix(1:28, nrow = 7, ncol = 4, byrow=TRUE))
    #   
    # png(paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/profile_approach/plots/",i,"_profiles_with_RNA_NEW.png"), width = 1800, height = 1800)
    # multiplot(plotlist = paired_graphs_forplot,  layout = matrix(1:28, nrow = 7, ncol = 4, byrow=TRUE))
    # dev.off()
    
      ############### REDO VST.ALLs with new csv_allpeaks!!!!!!!!!! ### APPLIED plotEigengene2 as no scale
      
      vst_multi_histones <- readRDS(file = paste0("~/R_objects/vst_dds.time.",i,"_vs_all_NEW.rds"))
      
      meta.histones <- read.csv("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/Meta_histones.csv", header = T, sep = ";")
      colnames(meta.histones) <- c("X","SampleID","Histone","Time","Replicate","Group1")
      meta.histones$Time <- as.integer(substr(meta.histones$Time, 4,4))
      
      #normalizing histones
      for (z in c(0,2,4,6))  {
        vst_multi_histones[,(c(1,2,9,10,17,18,25,26)+z)] <- vst_multi_histones[,(c(1,2,9,10,17,18,25,26)+z)] / max(melt(vst_multi_histones[,(c(1,2,9,10,17,18,25,26)+z)])[,3])
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
      
      df2 <- t(vst_multi_histones)
      
      multigraphs_forplot <- lapply (1:length(profileList), function(x) {
        elemento <- profileList[[x]]
        nombre <- names(profileList)[x]
        regions <- colnames(df2)[profileList[[x]]]
        number_regions <- length(regions)
        ploty <- ggplot()
        
        if (number_regions == 0) {
          ploty
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
        #paired_multigraphs <- append(paired_multigraphs,profile_allmC_plotting[x], after = length(paired_multigraphs))
        paired_multigraphs
      }
      
      multiplot(plotlist = paired_multigraphs[1:12],  layout = matrix(1:12, nrow = 3, ncol = 4, byrow=TRUE))
      
#       png(paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/profile_approach/plots/",i,"_MULTI_with_RNA_NEW.png"), width = 1800, height = 1800)
#       multiplot(plotlist = paired_multigraphs,  layout = matrix(1:28, nrow = 7, ncol = 4, byrow=TRUE))
#       dev.off()
#  
#       png(paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/profile_approach/plots/",i,"_MULTI_with_RNA_allmC_legend.png"), width = 200, height = 400)
#       plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#       legend("center", legend =c(levels(meta.histones$Histone),"MedSeq","RNAseq"), lty = c(5,5,5,1,1,1), bty='n',
#              col = c("deepskyblue", "darkorange2","darkgreen","red3","purple","black"), ncol=2)
#       dev.off()
#       
      
      
      
      
      
      # png(paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/profile_approach/plots/",i,"_MULTI_with_RNA.png"), width = 1800, height = 1800)
      # multiplot(plotlist = paired_multigraphs,  layout = matrix(1:28, nrow = 7, ncol = 4, byrow=TRUE))
      # dev.off()
      # 
    
    
    ###########################
    ###########################
    
    
    # OUTPUT TABLES -----------------------------------------------------------
    
    
    profile_data <- lapply (1:length(profileList), function(x) {
                    elemento <- profileList[[x]]
                    nombre <- names(profileList)[x]
                    subset_table <- supermatrix[elemento,]
                    subset_table
                    
                })
    names(profile_data) <- sprintf("profile%02d", 1:14)
    
    for(x in 1:length(profile_data)) {
        elemento <- profile_data[[x]]
        nombre <- names(profileList)[x]
        
        # FOR HOMER
        out1="/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/profile_approach/new_submits/homer_motifs/"
        #REMOVE Chr and converting to numeric to avoid ""
        elemento[,1] <- gsub('Chr', '', elemento[,1])
        elemento[,1] <- as.numeric(as.character(elemento[,1]))
        write.table(elemento[1:3], file = c(paste0(out1,i,"_",nombre,"_forhomer.bed")), append = FALSE, sep="\t", dec = ".",
                    row.names = FALSE, col.names = F)
        
        # FOR GOs
        # Write fused coordinates to match araport
        out2="/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/GOs/by_profile/"
        
        # Load fused coordinates
        fused_region <- read.table(text=rownames(elemento),col.names=c("region"))
        
        # Load araport annotation
        fused_reference <- read.csv(file = paste0("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/araport_annotation/k",i,"_araport_fused_NEW.txt"),sep = "\t", header = T)
        
        # Match 2 lists
        fused_region_merged <-fused_reference$gene_ID[fused_reference$region %in% fused_region$region]
        
        # Split multiple entries in gene_ID
        gene_IDs_list <- unlist(strsplit(as.character(fused_region_merged), ","))
        
        # Write table of gene IDs
        write.table(gene_IDs_list, file = paste0("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/profile_approach/new_submits/GOs/",i,"/",i,"_",nombre,"_geneIDs.txt"), sep ="", row.names = F, col.names = F, quote=FALSE )
        
    }


    
   