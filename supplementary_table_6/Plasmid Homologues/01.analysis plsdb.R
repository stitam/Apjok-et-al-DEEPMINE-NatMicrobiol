library(plyr)
rm(list = ls())

annotation <- read.delim("BlastN output.tsv", header=FALSE)
names(annotation)[1:11] <- c(
  "Prodigal.Query.ID", "ACC_NUCCORE", "qstart", "qend","ssstart","send",
  "evalue","bitscore","pident","qcovs","qcovhsp"
)
annotation_1 <-ddply(annotation, .(Prodigal.Query.ID), function(x)x[which.min(x$evalue), ])
names(annotation_1)[2]<- "ACC_NUCCORE"
  
db <- read.delim("Plasmid Sequences included in PSDB database.tsv")
length(unique(db$ACC_NUCCORE))
db <- db[-c(1,3,4,5.67,10,11,12,13,14,15,16,18,19,22,40:49)]
m <- merge(db, annotation_1, by = "ACC_NUCCORE")

write.csv(m,"Plasmid Homologues.csv",row.names = FALSE)

