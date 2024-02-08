### Shell bedtools intersects
#for known features
in=/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/2021_sept_allpeak/
bedtools intersect -wb -a $in'noChr/''day7_k9me2_noChr.bed' \
-b /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/araport11_complete_NEW.bed  > $in'all_peaks_k9me2_vsAraport.bed'

#for non feature peaks (aka distal_intergenic)
bedtools intersect -v -a $in'noChr/''day7_k9me2_noChr.bed' \
-b /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/araport11_complete_NEW.bed  > all_peaks_k9me2_NOoverlap_vsAraport.bed

# COUNT and record how many peaks: wc -l ./day7_k9me2_noChr.bed >> k9me2_peakcount.txt

dir_start <- '/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/2021_sept_allpeak/'
data <- read.table(paste0(dir_start,'all_peaks_k9me2_vsAraport.bed'), sep='\t', header=F)
library(plyr)
library(reshape)


# peaks to consider for distal
total_peaks <- system("awk 'END { print NR }'  /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/2021_sept_allpeak/noChr/day7_k9me2_noChr.bed", intern = TRUE)
total_peaks <- as.numeric(total_peaks)

# ADD intergenic by counting lines of file made by not intersecting peaks between araport and peak_file
intergenic <- system("awk 'END { print NR }'  /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/2021_sept_allpeak/all_peaks_k9me2_NOoverlap_vsAraport.bed", intern = TRUE)
intergenic <- as.numeric(intergenic)

percent_distal <- round((intergenic / total_peaks)*100,1) #this determine the remaining percentatge for araport features

temp_row <- data.frame(x='distal_intergenic', freq=intergenic,percent=percent_distal )



feature_counts <- count(data$V7)
feature_counts$percent <- round((feature_counts$freq / length(data$V7) * (100-percent_distal)), 1)

total_counts_partial <- rbind(feature_counts,temp_row)
total_counts_partial <- total_counts_partial[order(-total_counts_partial$percent),]
total_counts_partial <- total_counts_partial [c(1:5),]
other_percent <- 100-sum(total_counts_partial$percent)
temp_row2 <- data.frame(x='other', freq='na',percent=other_percent)
total_counts_partial <- rbind(total_counts_partial,temp_row2)

total_counts_frame <- data.frame(total_counts_partial$x)
rownames(total_counts_frame)<-total_counts_partial$x

total_counts_9me2 <- cbind(total_counts_frame[,-1], h3k9me2=total_counts_partial$percent)

#final table, manually shitty curation
total_counts <-
total_counts <- merge(total_counts_27ac, total_counts_27me3, by=0, all=TRUE) # by=0 (by rownames)
rownames(total_counts) <- total_counts$Row.names
total_counts <- merge(total_counts[,c(2:3)], total_counts_4me3, by=0, all=TRUE)
rownames(total_counts) <- total_counts$Row.names
total_counts <- merge(total_counts[,c(2:4)], total_counts_9me2, by=0, all=TRUE)
rownames(total_counts) <- total_counts$Row.names
total_counts <- total_counts[,c(2:5)]
total_counts <- total_counts[order(-total_counts$h3k27ac),]
total_counts[is.na(total_counts)] = 0
lnc_RNAs <- total_counts[4,]+total_counts[6,]
rownames(lnc_RNAs)<-'lncRNA'
total_counts <- rbind(total_counts[c(3,1,2),],lnc_RNAs[,],total_counts[c(7,5),])

#Transform in something ggplotable
total_plot <- data.matrix(total_counts)
melted_table <- melt(total_plot)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

safe_colorblind_palette[1:6]

p2<-ggplot(data=melted_table, aes(x=X2, y=value, fill=X1)) +
    geom_bar(stat="identity")+ scale_fill_brewer(palette = "Set3", direction = 1)+
    scale_y_continuous(breaks = c(0,5,10,25,50,75,100))+
    theme(panel.grid = element_blank(), axis.title = element_blank(), legend.position="right", legend.title = element_blank(), 
          axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))
p2

pdf(file="/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/annotation_byintersect/2021_sept_allpeak/day7_genomicfeatures.pdf",
    useDingbats=FALSE)
p2<-ggplot(data=melted_table, aes(x=X2, y=value, fill=X1)) +
    geom_bar(stat="identity")+ scale_fill_brewer(palette = "Set3", direction = 1)+
    scale_y_continuous(breaks = c(0,5,10,25,50,75,100))+
    theme(panel.grid = element_blank(), axis.title = element_blank(), legend.position="right", legend.title = element_blank(), 
          axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))
p2
dev.off()



