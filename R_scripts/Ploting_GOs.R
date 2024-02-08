#/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/GO_genemerge/genemerged/Bonferroni_

#Ploting GOs enrichments
#Alexander Vergara
#May 31 2019
#Here we will plot GOs enrcihments analysis 
# but we will plot also non signifant data in other color

#These analysis were performed using GOs annotatios 
# scripts at /mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-mediator-stress/Carmen_analysis/scripts

setwd("~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/GO_genemerge/genemerged")

library(ggplot2)
library(plyr)


###Bonferroni 0.05 ##################################################################################################################
##Up-regulated
k27ac_day0to1_up <- read.table("/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/GO_genemerge/genemerged/Bonferroni_0.05/k27ac_day0to1_up", sep= "\t", header = T, quote = "")
k27ac_day0to1_up$Comparison <- rep("k27ac_day0to1_up", times = length(k27ac_day0to1_up$Comparison)) 

k27ac_day1to4_up <- read.table("/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/GO_genemerge/genemerged/Bonferroni_0.05/k27ac_day1to4_up", sep= "\t", header = T, quote = "")
k27ac_day1to4_up$Comparison <- rep("k27ac_day1to4_up", times = length(k27ac_day1to4_up$Comparison)) 

k27ac_day4to7_up <- read.table("/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/GO_genemerge/genemerged/Bonferroni_0.05/k27ac_day4to7_up", sep= "\t", header = T, quote = "")
k27ac_day4to7_up$Comparison <- rep("k27ac_day4to7_up", times = length(k27ac_day4to7_up$Comparison)) 


#REVISAR SGTE LINEA, tuve que editar a mano la linea de GO:0102965 pues su descripcion no aparece en GO file 
# GO:0102965  alcohol-forming fatty acyl-CoA reductase activity
#WT_3h_vs_cdke1_3h_Up <- read.table("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-mediator-stress/Carmen_analysis/Enrichments_results/Bonferroni_0.05/WT_3h_vs_cdke1_3h_Up", sep= "\t", head=T)
#head(WT_3h_vs_cdke1_3h_Up)

#WT_72h_vs_cdke1_72h_Up<- read.table("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-mediator-stress/Carmen_analysis/Enrichments_results/Bonferroni_0.05/WT_72h_vs_cdke1_72h_Up", sep= "\t", head=T)
#head(WT_72h_vs_cdke1_72h_Up)

####################################
##Down-regulated 
k27ac_day0to1_down <- read.table("/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/GO_genemerge/genemerged/Bonferroni_0.05/k27ac_day0to1_down", sep= "\t", head=T, quote = "")
k27ac_day0to1_down$Comparison <- rep("k27ac_day0to1_down", times = length(k27ac_day0to1_down$Comparison)) 

k27ac_day1to4_down <- read.table("/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/GO_genemerge/genemerged/Bonferroni_0.05/k27ac_day1to4_down", sep= "\t", head=T, quote = "")
k27ac_day1to4_down$Comparison <- rep("k27ac_day1to4_down", times = length(k27ac_day1to4_down$Comparison)) 

k27ac_day4to7_down <- read.table("/mnt/picea/home/mquevedo/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/GO_genemerge/genemerged/Bonferroni_0.05/k27ac_day4to7_down", sep= "\t", head=T, quote = "")
k27ac_day4to7_down$Comparison <- rep("k27ac_day4to7_down", times = length(k27ac_day4to7_down$Comparison)) 

###################################
#########################################################
combined_up <- rbind.fill(k27ac_day0to1_up, k27ac_day1to4_up, k27ac_day4to7_up)
combined_down <- rbind.fill(k27ac_day0to1_down, k27ac_day1to4_down, k27ac_day4to7_down)

head(combined_up)
head(combined_down)

#Remove "roots tags" before plot (Biological process GO:0008150, Cellullar component GO:0005575, Molecular Function GO:0003674)

combined_up <- combined_up[!grepl("GO:0008150", combined_up$GO_ID),]
combined_up <- combined_up[!grepl("GO:0005575", combined_up$GO_ID),]
combined_up <- combined_up[!grepl("GO:0003674", combined_up$GO_ID),]

combined_down <- combined_down[!grepl("GO:0008150", combined_down$GO_ID),]
combined_down <- combined_down[!grepl("GO:0005575", combined_down$GO_ID),]
combined_down <- combined_down[!grepl("GO:0003674", combined_down$GO_ID),]

head(combined_up)
tail(combined_up)
#########################################################
###NOW PLOTING
#########################################################
#Up
ggplot(combined_up,aes(x=Comparison,y=Description))+ geom_point(aes(size=Percentage, color=Adjusted_Pvalue))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("GO enrichments Up-regulated genes \n")+scale_shape_discrete(solid=F)+
  scale_colour_gradientn(colours=c("royalblue4", "lightskyblue1")) 

#Down
ggplot(combined_down,aes(x=Comparison,y=Description))+ geom_point(aes(size=Percentage, color=Adjusted_Pvalue))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("GO enrichments Down-regulated genes \n")+scale_shape_discrete(solid=F)+
  scale_colour_gradientn(colours=c("royalblue4", "lightskyblue1")) 

#######################################################################################
############################################################

head(combined_up)
head(combined_down)

#DONE
###################################################################################
## Plotting Up and Down-regulated genes together 
##########################################################################################################
combined_up[,"Direction"]  <- "Up"
head(combined_up)
combined_down[,"Direction"]  <- "Down"
head(combined_down)


###########################################
###################

combined_all_up_and_down <- rbind(combined_up, combined_down)

head(combined_all_up_and_down)
tail(combined_all_up_and_down)
#######################################
# plot using a multi panel option of ggplot2 (facet_wrap).
# remove unnecesary information of Comparison column
combined_all_up_and_down$Comparison <- gsub('WT_vs_cdke1_Up','WT vs cdke1', combined_all_up_and_down$Comparison)
combined_all_up_and_down$Comparison <- gsub('WT_3h_vs_cdke1_3h_Up','WT 3h vs cdke1 3h', combined_all_up_and_down$Comparison)
combined_all_up_and_down$Comparison <- gsub('WT_72h_vs_cdke1_72h_Up','WT 72h vs cdke1 72h', combined_all_up_and_down$Comparison)

combined_all_up_and_down$Comparison <- gsub('WT_vs_cdke1_Down','WT vs cdke1', combined_all_up_and_down$Comparison)
combined_all_up_and_down$Comparison <- gsub('WT_3h_vs_cdke1_3h_Down','WT 3h vs cdke1 3h', combined_all_up_and_down$Comparison)
combined_all_up_and_down$Comparison <- gsub('WT_72h_vs_cdke1_72h_Down','WT 72h vs cdke1 72h', combined_all_up_and_down$Comparison)

head(combined_all_up_and_down)
tail(combined_all_up_and_down)

#Now we made the categories like factors to sort and we change the facet_wrap
head(combined_all_up_and_down) 
tail(combined_all_up_and_down)
combined_all_up_and_down$Direction_ok = factor(combined_all_up_and_down$Direction, levels=c('Up','Down'))
#Now we made the factors to sort time treatment
head(combined_all_up_and_down) 
combined_all_up_and_down$Comparison_ok = factor(combined_all_up_and_down$Comparison, levels=c('WT vs cdke1','WT 3h vs cdke1 3h','WT 72h vs cdke1 72h'))
####
#Up and Down together
myplot <-ggplot(combined_all_up_and_down,aes(x=Comparison_ok,y=Description))+ geom_point(aes(size=Counts, color=Adjusted_Pvalue))+ 
  facet_wrap(~Direction_ok, nrow=1) +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("GO enrichments Up and Down-regulated genes\n")+scale_shape_discrete(solid=F)+
  scale_colour_gradientn(colours=c("royalblue4", "lightskyblue1"))
myplot
myplot + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#Done
###HASTA ACA TODO OK 
#############################################################################################
#   END   *********************************************************************************
##*********************************************************************************
##*********************************************************************************
head(combined_all_up_and_down)

#we save the combined_all_up_and_down object as a file 
#to plot later these ones including in gray the non-sginficant counts values
write.table(combined_all_up_and_down,file="combined_all_up_and_down_sgnificant_GOs.txt", quote = F ,sep = "\t", row.names=FALSE)

bk = c(seq(0,0.01,length=100),seq(0.01,0.05,length=100),seq(0.05,1,length=100))
####Intento 1
myplot <-ggplot(combined_all_up_and_down,aes(x=Comparison_ok,y=Description))+ geom_point(aes(size=Counts, color=Adjusted_Pvalue))+ 
  facet_wrap(~Direction_ok, nrow=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("GO enrichments\n")+scale_shape_discrete(solid=F) + labs(title = "GO enrichments", x = " ", y = "GO")+
  scale_colour_gradientn(colours = c("royalblue4", "lightskyblue1"))
myplot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 360))

#######################
#END
###########
pdf("GO_enrichments.pdf", width=12, height=10)
myplot <-ggplot(combined_all_up_and_down,aes(x=Comparison_ok,y=Description))+ geom_point(aes(size=Counts, color=Adjusted_Pvalue))+ 
  facet_wrap(~Direction_ok, nrow=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("GO enrichments\n")+scale_shape_discrete(solid=F) + labs(title = "GO enrichments ", x = " ", y = "GO")+
  scale_colour_gradientn(colours = c("royalblue4", "lightskyblue1"))
myplot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(strip.text.x = element_text(size = 8, colour = "black", angle = 360))
dev.off()
######
myplot <-ggplot(combined_all_up_and_down,aes(x=Comparison_ok,y=Description))+ geom_point(aes(size=Counts, color=Adjusted_Pvalue))+ 
  facet_wrap(~Direction_ok, nrow=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #ggtitle("GO enrichments \n")
  scale_shape_discrete(solid=F)+
  scale_colour_gradientn(colours=c("royalblue4", "lightskyblue1"))
myplot
p4 <- myplot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ theme(strip.text.x = element_text(size = 8, colour = "black", angle = 360))
p4
p5 <- p4 +  theme(axis.text.y=element_text(angle=0, hjust=1, size=12)) +
  theme(axis.text.x=element_text(angle=45, hjust=1,size=12)) +
  labs(x="",y="GO")
print(p5) 
p6 <-p5+theme(strip.text.x = element_text(face = 'bold',size = 12, colour = "black", angle = 0))
p7 <- p6 + labs(color='Adjusted Pvalue') 
p7
p8 <- p7 + theme(axis.title.y = element_text(face = 'bold', angle = 90, size = 14))
p8
### Now we save this version in 600dpi for manuscript
####SAVE FIG 
# Calculate the height and width (in pixels) for a 10x10-inch image at 300 ppi
ppi <- 600
png("Fig_GOs_transcriptome_universe.png", width=12*ppi, height=10*ppi, res=ppi)
plot(p8)
dev.off()
##################################
###*******************************************************************************************************************

