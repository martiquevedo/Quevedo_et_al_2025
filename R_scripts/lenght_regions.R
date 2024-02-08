
#set multiplot
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

setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/macs2/replicate_merging/H3_Input_common/final_merge/final/")
  
  
i= c("k27ac")
file_lists <-  dir(".", pattern=i, recursive = F, full.names = TRUE) 
  
      table_raw <- read.table(file_lists, stringsAsFactors = F)
      table_raw <- table_raw[,2:3]
      colnames(table_raw)<- c("Start","End")
      table_raw$Lenght <- table_raw$End-table_raw$Start
      head(table_raw)
      forgg <- as.data.frame(table_raw$Lenght, )
      
      ploty<- ggplot(forgg, aes(x=table_raw$Lenght)) + geom_histogram(color="black", fill="deepskyblue", binwidth=250)+
        labs(title=paste0("Lenght distribution of H3",i," regions"), x ="Lenght (bp)", y = "Frequency")+      
        geom_vline(aes(xintercept=mean(table_raw$Lenght)), color="black", linetype="dashed", size=1)+
        scale_x_continuous(limit = c(0, 10000), breaks = c(0,750,1500,2500,5000,7500,10000))+
        geom_text(aes(x=mean(table_raw$Lenght) + 250), label=paste0("Mean= ",round(mean(table_raw$Lenght),digits=0),"bp"), y=5000, colour="Black", hjust = 0, fontface="plain")
      
      
      ppi<-600
      png(paste0("/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Plots/k27ac_lenght_regions.png"), width=6*ppi, height=5*ppi, res=ppi)
      ploty
      dev.off()
      
      

