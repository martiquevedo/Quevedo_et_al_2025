

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


plotEigengene <- function(data, genes, condition, time, timeUnits = "Time",
                          inverse = F, title = "", noGrid = T, colors = c("deepskyblue", "#BE9230"),
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




setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/")


meta <- read.csv("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/RNA_samples_forDiffbind.csv", header = T, sep = ";")
indir <- "/mnt/picea/projects/arabidopsis/aastrand/greening/results/STAR"
list.files(indir)
filenames <- file.path(indir, paste0(meta$filenames))

library("Rsamtools")

bamfiles <- BamFileList(filenames, yieldSize = 2000000)
seqinfo(bamfiles[1])

library(Rsubread)
library(rtracklayer)

## import the bed file
bed.ranges <- import.bed("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/ALL_peaks.bed")
## export as a gff3 file
export.gff3(bed.ranges,"~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/ALL_peaks.gff3")
gff_allpeaks<-"~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Metadata/ALL_peaks.gff3"


featureCounts(files = filenames[1], annot.ext="annotation.gtf")
