rm(list=ls())


#1. introduce threshold
BlastN <- read.delim("BlastN output.tsv", sep="\t")
str(BlastN)
BlastN_thr <- BlastN[BlastN$pident>95 &BlastN$qcovhsp>95,]
BlastN_thr <- BlastN_thr[-c(3:11)]
write.csv(BlastN_thr ,"BlastN 95 thr.csv",row.names = FALSE)


