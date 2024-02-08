
setwd("/mnt/picea/projects/arabidopsis/aastrand/Arabidopsis-greening-chipseq/Literature/comparison")

library(VennDiagram)

A <- read.table("day7_k27ac_peaks", sep = "\t")
A <- A[!duplicated(A), ]
A <- length(A$V1)

B <- read.table("seedlings_H3K27ac_peaks", sep = "\t")
B <- B[!duplicated(B), ]
B <- length(B$V1)

C <- read.table("seedlings_H3K9ac_peaks", sep = "\t")
C <- C[!duplicated(C), ]
C <- length(C$V1)

AB <- read.table("AB_intersected", sep = "\t")
AB <- AB[!duplicated(AB), ]
AB <- length(AB$V1)

AC <- read.table("AC_intersected", sep = "\t")
AC <- AC[!duplicated(AC), ]
AC <- length(AC$V1)

BC <- read.table("BC_intersected", sep = "\t")
BC <- BC[!duplicated(BC), ]
BC <- length(BC$V1)

ABC <- read.table("ABC_intersected", sep = "\t")
ABC <- ABC[!duplicated(ABC), ]
ABC <- length(ABC$V1)

###### INSERT VALUES IN http://eulerr.co/

library(eulerr)
pal2<-c("deepskyblue","dodgerblue4","darkorchid")
MyVenn <- euler(c(A = A, B = B, C = C, "A&B" = AB, "A&C" = AC,"B&C" = BC ,"A&B&C" = ABC ),input = "union")
plot(
    MyVenn,
    fills = pal2,
    edges = F,
    legend =list(labels = c("Day7_H3K27ac", "Seedlings_H3K27ac","Seedlings_H3K9ac")),
    labels = identical(legend, FALSE),
    quantities = F,
    strips = NULL,
    main = NULL,
    n = 200L,
    adjust_labels = TRUE)

legend("right", legend=c("Test1", "Test2", "test3"),
       col=pal2, lty=1:2, cex=0.8)

