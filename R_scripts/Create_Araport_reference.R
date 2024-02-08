

All_araport<-read.csv(file = "~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Annotation/Araport11_ALL.txt",sep = "\t", header = F)
araport <- All_araport
head(araport)

library(stringr)
araport$V9<-str_split_fixed(araport$V9, ",", 2)

IDs <- araport$V9[,1]

araport<- araport[,c(1,4,5,3,7)] 
araport$V6 <- IDs
araport <- araport[order(araport$V1, araport$V4),]


#REMOVE Chr and converting to numeric to avoid ""
araport[,1] <- gsub('Chr', '', araport[,1])
araport[,1] <- as.numeric(as.character(araport[,1]))
####

write.table(araport,file = "~/Git/UPSCb/projects/arabidopsis-greening-ChIP-Seq/docs/Annotation/Araport11_final_forintersect",sep = "\t",row.names = F, col.names = F, quote = F )
