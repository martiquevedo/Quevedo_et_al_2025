
setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Literature/comparison/h3k4me3/")

A <- read.table("day7_4me3_peaks", sep = "\t")
A <- A[!duplicated(A), ]
A <- length(A$V1)

B <- read.table("seedling_h3k4me3_peaks", sep = "\t")
B <- B[!duplicated(B), ]
B <- length(B$V1)

AB <- read.table("AB_intersected", sep = "\t")
AB <- AB[!duplicated(AB), ]
AB <- length(AB$V1)


library(eulerr)
pal2<-c("darkgreen","darkolivegreen3")
MyVenn <- euler(c(A = A, B = B,"A&B" = AB),input = "union")
plot(
    MyVenn,
    fills = pal2,
    edges = F,
    legend =F,
    labels = identical(legend, FALSE),
    quantities = T,
    strips = NULL,
    main = NULL,
    n = 200L,
    adjust_labels = TRUE)
