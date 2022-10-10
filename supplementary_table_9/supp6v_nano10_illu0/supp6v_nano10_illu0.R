library(stringr)
library(tidyverse)
library(dplyr)
library(plyr)
library(ggplot2)
library(grid)
rm(list = ls())

nano_thr <- 10 #nanopore threshold at least this many reads should be included
illu_thr <- 0 #illumina threshold at least this many reads should be included

print(nano_thr)
print(illu_thr)

##1. READ TRHOUGH THE FILE LIST OF NANOPORE BARCODES FROM FIRST RUN
dir_path_first_run = "../../run1/nanopore/"
filelist_first_run <- list.files(dir_path_first_run,pattern="*.fasta", full.names=TRUE)

#create a list
lst_first_run <- lapply(filelist_first_run, read.delim, header = TRUE, sep = "\n", stringsAsFactors=FALSE)

#convert list to dataframe (takes a couple of minutes to run)
nanopore_table_first_run <- lst_first_run %>%
  data.table::rbindlist(fill=T)%>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)%>%
  gather(variables, Contig.Sequence, -1, na.rm = TRUE) %>%
  select(-variables)

#break down information for each sequence (nanopores downstream and upstream are inverted)
nanopore_table_first_run$name <- gsub("\\.", "_", nanopore_table_first_run$name)
nanopore_table_first_run <- data.frame(nanopore_table_first_run,do.call(rbind,str_split(nanopore_table_first_run$name,"_")))
names(nanopore_table_first_run)[c(4,5,8,11)] <- c("Downstream.Barcode", "Upstream.Barcode", "Contig.Length", "Nanopore.Read.Count")
nanopore_first_run <- nanopore_table_first_run[, c(1,2,4,5,8,11)]
merge_nanopore_first_run <- nanopore_first_run[,c(1,3,4)]

#merge with annotations and mobile geen pool network
#A) BLASTX
blastx__first_run <- read.delim("../../run1/results/Nanopore sequences Run1 blastx.tsv", stringsAsFactors = FALSE)
blastx__first_run <- data.frame(blastx__first_run,do.call(rbind,str_split(blastx__first_run$qseqid,"_")))
names(blastx__first_run)[c(1,2,15:16)]<-c("Contig","ID","Downstream.Barcode","Upstream.Barcode")
blastx__first_run_name <- merge(blastx__first_run, merge_nanopore_first_run, by = c("Downstream.Barcode", "Upstream.Barcode"))
#B) PRODIGAL
prodigal_first_run <- read.delim("../../run1/results/Nanopore sequences Run1 prodigal.tsv")
prodigal_first_run <- data.frame(prodigal_first_run,do.call(rbind,str_split(prodigal_first_run$Contig,"_")))
names(prodigal_first_run)[c(5,19:20)]<-c("ID","Downstream.Barcode","Upstream.Barcode")
prodigal_first_run_name <- merge(prodigal_first_run, merge_nanopore_first_run, by = c("Downstream.Barcode", "Upstream.Barcode"))
#C) horizontal gene transfer annotation based on mobile gene pool network
hgt_first_run <- read.delim("../../run1/results/Nanopore sequences Run1 mobile gene pool.tsv")
hgt_first_run <- data.frame(hgt_first_run,do.call(rbind,str_split(hgt_first_run$Nanopore.contig,"_")))
names(hgt_first_run)[c(2,52,53)]<-c("Query.ID","Downstream.Barcode","Upstream.Barcode")
hgt_first_run_name <- merge(hgt_first_run, merge_nanopore_first_run, by = c("Downstream.Barcode","Upstream.Barcode"))




##2. READ THROUGH THE FILE LIST OF NANOPORE BARCODES FROM SECOND RUN
dir_path_second_run = "../../run2/nanopore/"
filelist_second_run <- list.files(dir_path_second_run,pattern="*.fasta", full.names=TRUE)

#create a list
lst_second_run <- lapply(filelist_second_run, read.delim, header = TRUE, sep = "\n", stringsAsFactors=FALSE)

#convert list to dataframe (takes a couple of minutes to run)
nanopore_table_second_run <- lst_second_run %>%
  data.table::rbindlist(fill=T)%>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)%>%
  gather(variables, Contig.Sequence, -1, na.rm = TRUE) %>%
  select(-variables)

#break down information for each sequence (nanopores downstream and upstream are inverted)
nanopore_table_second_run$name <- gsub("\\.", "_", nanopore_table_second_run$name)
nanopore_table_second_run <- data.frame(nanopore_table_second_run,do.call(rbind,str_split(nanopore_table_second_run$name,"_")))
names(nanopore_table_second_run)[c(4,5,8,11)] <- c("Downstream.Barcode", "Upstream.Barcode", "Contig.Length", "Nanopore.Read.Count")
nanopore_second_run <- nanopore_table_second_run[, c(1,2,4,5,8,11)]
merge_nanopore_second_run <- nanopore_second_run[,c(1,3,4)]

#merge with annotations and mobile geen pool network
#A) BLASTX
blastx__second_run <- read.delim("../../run2/results/Nanopore sequences Run2 blastx.tsv", stringsAsFactors = FALSE)
blastx__second_run <- data.frame(blastx__second_run,do.call(rbind,str_split(blastx__second_run$qseqid,"_")))
names(blastx__second_run)[c(1,2,15:16)]<-c("Contig","ID","Downstream.Barcode","Upstream.Barcode")
blastx__second_run_name <- merge(blastx__second_run, merge_nanopore_second_run, by = c("Downstream.Barcode", "Upstream.Barcode"))
#B) PRODIGAL
prodigal_second_run <- read.delim("../../run2/results/Nanopore sequences Run2 prodigal.tsv")
prodigal_second_run <- data.frame(prodigal_second_run,do.call(rbind,str_split(prodigal_second_run$Contig,"_")))
names(prodigal_second_run)[c(5,19:20)]<-c("ID","Downstream.Barcode","Upstream.Barcode")
prodigal_second_run_name <- merge(prodigal_second_run, merge_nanopore_second_run, by = c("Downstream.Barcode", "Upstream.Barcode"))
#C) horizontal gene transfer annotation based on mobile gene pool network
hgt_second_run <- read.delim("../../run2/results/Nanopore sequences Run2 mobile gene pool.tsv")
hgt_second_run <- data.frame(hgt_second_run,do.call(rbind,str_split(hgt_second_run$Nanopore.contig,"_")))
names(hgt_second_run)[c(2,52,53)]<-c("Query.ID","Downstream.Barcode","Upstream.Barcode")
hgt_second_run_name <- merge(hgt_second_run, merge_nanopore_second_run, by = c("Downstream.Barcode", "Upstream.Barcode"))




##3. MERGE THE NANOPORE FROM BOTH RUNS AND KEEP ONES WITH THE UNIQUE COMBINATION OF BARCODES HAVING THE HIGHEST RC
nanopore_merge <- rbind(nanopore_second_run, nanopore_first_run)
nanopore_merge <- transform(nanopore_merge, Contig = paste(Downstream.Barcode, Upstream.Barcode, sep = "_"))
str(nanopore_merge)
nanopore_merge$Nanopore.Read.Count <- as.numeric(nanopore_merge$Nanopore.Read.Count)
nanopore <-ddply(nanopore_merge, .(Contig), function(x)x[which.max(x$Nanopore.Read.Count), ])
nanopore_contigs <- unique(nanopore$Name)



##4. REMOVE DUBLICATES OF UPSTREAM AND DOWNSTEAM BARCODES TO REDUCE SEQUENCE REUNDANCY
nanopore_unique_down <-ddply(nanopore, .(Downstream.Barcode), function(x)x[which.max(x$Nanopore.Read.Count), ])
nanopore_unique_up_down <- ddply(nanopore_unique_down, .(Upstream.Barcode), function(x)x[which.max(x$Nanopore.Read.Count), ])
unique_nano_upstream_barcodes <- unique(nanopore_unique_up_down$Upstream.Barcode)
unique_nano_downstream_barcodes <- unique(nanopore_unique_up_down$Downstream.Barcode)



##5. KEEP THE ANNOTATION FOR THE CONIGS WIHT THE HIGHEST READ COUNTS FROM THE TWO RUNS FROM
nanopore_names <- nanopore_unique_up_down[,c(1,3,4)]

#A) BLASTX
blastx <- rbind(blastx__first_run_name, blastx__second_run_name)
blastx_corrected <- merge(blastx, nanopore_names, by=c("Downstream.Barcode", "Upstream.Barcode", "name"))
#B) PRODICAL
prodigal <- rbind(prodigal_first_run_name, prodigal_second_run_name)
prodigal_corrected <- merge(prodigal, nanopore_names, by=c("Downstream.Barcode", "Upstream.Barcode", "name"))
#C) MOBILE GENE POOL
hgt <- rbind(hgt_first_run_name, hgt_second_run_name)
hgt_corrected <- merge(hgt, nanopore_names, by=c("Downstream.Barcode", "Upstream.Barcode", "name"))




##6. PROCESS ILLUMINA DATA AND MERGE WITH NANOPORE
illumina <- read.csv("../../run3/results/HMT5GAFX3_bcnumbers.csv", sep = "\t")

# format columns
illumina <- illumina[,c(1,6,5,7)]
illumina <- dplyr::rename(
  illumina,
  Up.Down = "UP.DOWN",
  Barcode = "Mid.BC"
)
illumina$Barcode <- tolower(illumina$Barcode)

#illumina <- read.csv("Illumina barcodes Round1.csv", sep=",")
illumina <- data.frame(illumina,do.call(rbind,str_split(illumina$Sample,"_")))
#new line
illumina <- cbind(illumina[,1:5], data.frame(resistome = "combined"), illumina[,6])
names(illumina)[4:7] <- c("Illumina.Read.Count", "Host", "Resistome", "Antibiotic")
#job specific mod
illumina$Antibiotic <- ifelse(
  illumina$Antibiotic %in% illumina$Host, "none", illumina$Antibiotic
)

#A) Upstream barcodes
#exclude PXB antibiotic since not used in our study and Salmonella DOX
illumina_up <- illumina[!(grepl("PXB", illumina$Antibiotic)) & !(grepl("Salmonella_gut_DOX|Salmonella_soil_DOX|Salmonella_patho_DOX", illumina$Sample)) &illumina$Up.Down=="Upstream",]
illumina_up$Biological.Replicate <- 1

#include also the second run barcodes which is a replicate of Klebsiella to prove reproducibility of our pipeline
#second_illumina <- read.csv("Illumina barcodes Round2.csv", sep = "\t")
#replicate_Klebsiella <- subset(second_illumina, grepl("Klebsiella_K390", second_illumina$Sample))
#replicate_Klebsiella <- data.frame(replicate_Klebsiella,do.call(rbind,str_split(replicate_Klebsiella$Sample,"_")))
#names(replicate_Klebsiella)[c(4,5,7)] <- c("Illumina.Read.Count","Host", "Antibiotic")
#replicate_Klebsiella$Resistome <- "patho"
#replicate_Klebsiella <- replicate_Klebsiella[-c(1,6,8,9)]
#replicate_Klebsiella$Sample <-   paste(replicate_Klebsiella$Host, replicate_Klebsiella$Resistome, replicate_Klebsiella$Antibiotic, sep = "_")
#replicate_Klebsiella$Biological.Replicate <- 2

#merge the two replicates from upstream illumina barcodes indicating the biological replicates
#illumina_up_replicate <- rbind(illumina_up,replicate_Klebsiella)
#illumina_up_replicate_final1 <- illumina_up_replicate[c(1,3,8)]

#illumina_up_replicate_final2 <- aggregate(illumina_up_replicate_final1[3], illumina_up_replicate_final1[-3], FUN = function(X) paste(unique(X), collapse=", "))
#illumina_up_replicate_final3 <- illumina_up_replicate[-8]

#illumina_up_replicate_final4 <- left_join(illumina_up_replicate_final3, illumina_up_replicate_final2,by = c("Sample", "Barcode"))
#illumina_up_replicate_final4 <- illumina_up_replicate_final4[!duplicated(illumina_up_replicate_final4[c(1,3,8)]),]

#new line
illumina_up_replicate_final4 <- illumina_up

#decapilise the illumina so that you can do the mergering with nanopore
illumina_up_replicate_final4$Barcode=tolower(illumina_up_replicate_final4$Barcode) 
names(illumina_up_replicate_final4)[3] <- "Upstream.Barcode"

illumina_nano_up <- merge(nanopore_unique_up_down, illumina_up_replicate_final4, by = "Upstream.Barcode")

length(unique(nanopore_unique_up_down$Contig))
length(unique(illumina_nano_up$Contig))

#B) Downstream barcodes
#exclude PXB and Salmonella DOX, and keep downstream barcodes
illumina_down <- illumina[!(grepl("PXB", illumina$Antibiotic)) & !(grepl("Salmonella_gut_DOX|Salmonella_soil_DOX|Salmonella_patho_DOX", illumina$Sample)) &illumina$Up.Down=="Downstream",]

#decapilise the illumina so that you can do the merging with nanopore
illumina_down$Barcode=tolower(illumina_down$Barcode) 
illumina_down$Biological.Replicate <- 1
names(illumina_down)[3] <- "Downstream.Barcode"

illumina_nano_down <- merge(nanopore_unique_up_down, illumina_down, by = "Downstream.Barcode")
length(unique(illumina_nano_down$Contig))

#C) Merge upstream and downstream barcodes

illumina_nano_all <- rbind(illumina_nano_down, illumina_nano_up)
str(illumina_nano_all)
#if contig and sample is the same it means it aligns both to downstream and upstream so just leave one with highest read count
#2022-05-13 comment the line below, 

illumina_nano_all <- ddply(illumina_nano_all, .(Contig, Sample), function(x)x[which.max(x$Illumina.Read.Count), ])

#) put log values for Illumina RC so that you create scatter plot after
illumina_nano_all$log.nanopore <- log(illumina_nano_all$Nanopore.Read.Count)
illumina_nano_all$log.illumina <- log(illumina_nano_all$Illumina.Read.Count)






##7. APPLY TRESHOLDS FOR ILLUMINA AND NANANOPORE
all_thr_illumina_nano<- illumina_nano_all[illumina_nano_all$Illumina.Read.Count>=illu_thr & illumina_nano_all$Nanopore.Read.Count>=nano_thr,]
str(all_thr_illumina_nano)
all_thr_illumina_nano$Contig.Length <- as.numeric(all_thr_illumina_nano$Contig.Length)
all_thr_barcodes_illumina_nano <- unique(all_thr_illumina_nano$Contig)
all_thr_barcodes_illumina_nano_species <-ddply(all_thr_illumina_nano,~Host,summarise,number_of_unique_barcodes=length(unique(Contig)))
length(unique(all_thr_illumina_nano$Contig))






##8. EVALUATE CORRELATION BETWEEN ILLUMINA AND NANOPORE READ COUNT TO SUPPORT USE OF THRESHOLD
#no threshold
illumina_nano_all$Normalise.Nanopore <- illumina_nano_all$Nanopore.Read.Count/sum(illumina_nano_all$Nanopore.Read.Count)
illumina_nano_all$Normalise.Illumina <- illumina_nano_all$Illumina.Read.Count/sum(illumina_nano_all$Illumina.Read.Count)
illumina_nano_all$log.Normalise.Nanopore<- log(illumina_nano_all$Normalise.Nanopore)
illumina_nano_all$log.Normalise.Illumina <- log(illumina_nano_all$Normalise.Illumina)

res <- cor.test(illumina_nano_all$log.Normalise.Nanopore, illumina_nano_all$log.Normalise.Illumina, method = "spearman")
grob1 = grobTree(textGrob(paste("Before Thr: R = ", round(res$estimate, 2)), x = 0.67, y = 0.90, gp = gpar(fontsize = 18, fontface = "bold")))

plot1 <- ggplot(illumina_nano_all, aes(log.Normalise.Nanopore, log.Normalise.Illumina)) + geom_point(color="#CC3333",size = 1) + geom_smooth(method = "lm", se=FALSE, colour="#555555", size=2)
g <- plot1 +labs(x = "Nanopore Read Count (log scale)", y= "Illumina Read Count (log scale)") + theme_bw() +annotation_custom(grob1) +
  theme(axis.line = element_line(colour = "black"), text = element_text(size=22),
        panel.grid.major = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank() ,axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")))

#pdf(file = "rc_corrplot_before.pdf", width = 10, height = 8)
#g
#dev.off()

#with threshold
all_thr_illumina_nano$Normalise.Nanopore <- all_thr_illumina_nano$Nanopore.Read.Count/sum(all_thr_illumina_nano$Nanopore.Read.Count)
all_thr_illumina_nano$Normalise.Illumina <- all_thr_illumina_nano$Illumina.Read.Count/sum(all_thr_illumina_nano$Illumina.Read.Count)
all_thr_illumina_nano$log.Normalise.Nanopore <- log(all_thr_illumina_nano$Normalise.Nanopore)
all_thr_illumina_nano$log.Normalise.Illumina <-log(all_thr_illumina_nano$Normalise.Illumina)

res <- cor.test(all_thr_illumina_nano$log.Normalise.Nanopore, all_thr_illumina_nano$log.Normalise.Illumina, method = "spearman")
grob3 = grobTree(textGrob(paste("After Thr: R = ", round(res$estimate,3)), x = 0.67, y = 0.90, gp = gpar(fontsize = 18, fontface = "bold")))

g <- ggplot(all_thr_illumina_nano, aes(log.Normalise.Nanopore, log.Normalise.Illumina)) + geom_point(color="#CC3333",size = 1) + geom_smooth(method = "lm", se=FALSE, colour="#555555", size=2)+
  labs( x = "Nanopore Read Count (log scale)", y= "Illumina Read Count (log scale)") +annotation_custom(grob3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"), text = element_text(size=22),
        panel.grid.major = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank() ,axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")))        

fname <- paste0("rc_corrplot_after_thr - nano", nano_thr, " illu", illu_thr, ".pdf")
#pdf(file = fname, width = 10, height = 8)
#g
#dev.off()

all_thr_illumina_nano <- all_thr_illumina_nano[-c(17:20)]





##9. MERGE PRODIGAL ANNOTATION WITH THE NANOPORE AND ILLUMINA CONTIGS SURPASSING THE THRESHOLD
prodigal_illu_nano <- merge(prodigal_corrected, all_thr_illumina_nano, by = c("Downstream.Barcode","Upstream.Barcode", "Contig", "name"))

#A) create an extra column showing whether the barcode match to a CardRes ID or not
unique_barcodes_prodigal <- aggregate(prodigal_illu_nano$ID ~ prodigal_illu_nano$Contig, prodigal_illu_nano, paste, collapse = " ")
unique_barcodes_prodigal$Present.CardRes.ID <- ifelse(grepl("CardRes",unique_barcodes_prodigal$`prodigal_illu_nano$ID`), "Yes", "No")
names(unique_barcodes_prodigal)[1] <- "Contig"
unique_barcodes_prodigal <- unique_barcodes_prodigal[,c(1,3)]
prodigal_match_no_match <- merge(unique_barcodes_prodigal, prodigal_illu_nano, by = "Contig")

#B) the contigs that had no CardRes match keep them once for each unique sample
prodigal_no_match <- prodigal_match_no_match[prodigal_match_no_match$Present.CardRes.ID=="No",]

# EDIT TO ENSURE CGATCACCAA IS ROBUSTLY MAINTAINED (2022-06-13)
#prodigal_no_match_unique <- prodigal_no_match[!duplicated(prodigal_no_match[c(1,26)]),]
# NOTE THE CHANGE: DEDUPLICATION IS BASED ON QUERY ID AND NOT CONTIG ID
prodigal_no_match_unique <- prodigal_no_match[!duplicated(prodigal_no_match[c(6,26)]),]


#C) the contigs that have a match see whether there have ORFs matching to different resistance genes
prodigal_match <- prodigal_match_no_match[prodigal_match_no_match$Present.CardRes.ID=="Yes",]
#remove the fragmetns with no matches since there is already a match for these sequences
prodigal_match <- prodigal_match[!prodigal_match$ID=="No match",]

#D) keep contigs with unique CARD_Res_ID (if repeated keep the one with the minimal Expected value)
str(prodigal_match)
prodigal_match$Expected.value <- as.numeric(prodigal_match$Expected.value)
prodigal_match_higest_Evalue <- ddply(prodigal_match, .(Contig, ID, Sample), function(x)x[which.min(x$Expected.value), ])

#E) merge everything together
prodigal_final_annotation <- rbind(prodigal_no_match_unique, prodigal_match_higest_Evalue)
colnames(prodigal_final_annotation)[6:22] <- paste("Prodigal", colnames(prodigal_final_annotation[,c(6:22)]), sep = ".")





##10. USE BLASTX AS COMPLEMENARY TO PRODIGAL AND GIVE ANNOTAITON TO GENES USING CARDRES MERGED DATABASE
str(blastx_corrected)
blastx_higest_Evalue <- ddply(blastx_corrected, .(Contig), function(x)x[which.min(x$evalue), ])
colnames(blastx_higest_Evalue)[5:17] <- paste("Blastx", colnames(blastx_higest_Evalue[,c(5:17)]), sep = ".")

prodigal_blastx <- left_join(prodigal_final_annotation,blastx_higest_Evalue, by = "Contig", all=TRUE)
str(prodigal_blastx)
#use blastx as complemenrary method to prodigal
prodigal_blastx$Unify.Blastx.Prodigal.ID <- ifelse(!(prodigal_blastx$Prodigal.ID=="No match"), as.character(prodigal_blastx$Prodigal.ID), as.character(prodigal_blastx$Blastx.ID))

#read the database and create single column for a gene name 
database <- read.delim("../../databases/Cardres merged database.txt")
database[database==""]<- NA
database$Card.Res.Gene.Annotations <- ifelse(!is.na(database$CARD.Model.Name), as.character(database$CARD.Model.Name), as.character(database$RES.Gene_accession.no.))
database$Card.Res.Gene.Mechanism <- ifelse(!is.na(database$CARD.Resistance.Mechanism), as.character(database$CARD.Resistance.Mechanism), as.character(database$RES.Mechanism.of.resistance))
database$Card.Res.Gene.Classification <- ifelse(!is.na(database$CARD.AMR.Gene.Family), as.character(database$CARD.AMR.Gene.Family), as.character(database$RES.Notes))
database$Card.Res.Gene.Drug.Class <- ifelse(!is.na(database$CARD.Drug.Class), as.character(database$CARD.Drug.Class), as.character(database$RES.Class))

names(database)[1]<- "Unify.Blastx.Prodigal.ID"
prodigal_blastx_database <- left_join(prodigal_blastx,database, by = "Unify.Blastx.Prodigal.ID", all=TRUE)




##11. INFORMATION ON ANTIBIOTICS
capillary <- read.csv("../../supplementary_table_6/Information on Antibiotics and CFU values in vitro.csv")
capillary <- capillary[c(1,5,6)]
capillary <- capillary[!duplicated(capillary[c(1:3)]),]

# new line
capillary <- rbind(
  capillary,
  data.frame(Antibiotic = "none", Type = "none", Class = "none")
)

metagenome <- merge(capillary, prodigal_blastx_database, by = "Antibiotic")
length(unique(metagenome$Contig))
#total 768 have annotation wheres 35 don't
contigs_annotated <- aggregate(data=metagenome, Card.Res.Gene.Annotations ~Contig, function(x) length(unique(x)))



##12. FIX THE GENES MECHANISM AND CLASSIFICATION
genes<- ddply(metagenome, .(Card.Res.Gene.Annotations), summarise, mechanism = list(unique(na.omit(Card.Res.Gene.Mechanism))), classification = list(unique(na.omit(Card.Res.Gene.Classification))))
genes$mechanism <- vapply(genes$mechanism, paste, collapse = ", ", character(1L))
genes$classification <- vapply(genes$classification, paste, collapse = ", ", character(1L))
#write.table(genes, "Check gene mechanism and classification.tsv", row.names = FALSE, sep="\t")

#read the file you corrected and fix the metagenome data
genes_correct <- read.delim("../../supplementary_table_6/CORRECTED Check gene mechanism and classification.tsv")
metagenome_before_classification <- metagenome[,-c(14:19,22,24,35:36,42:47,50,52,73,74)]
metagenome_before_classification <- merge(metagenome_before_classification, genes_correct, by = "Card.Res.Gene.Annotations", all= TRUE)



##13. IF THE GENES ON THE SAME CONTIG ARE CLASSIFIED IN THE SAME GENE FAMILY LEAVE ONE OF GENES SINCE RESULT FRAGMENATION POSSIBLY
str(metagenome_before_classification)
metagenome_before_classification$Prodigal.Expected.value <- as.numeric(metagenome_before_classification$Prodigal.Expected.value)
metagenome_before_classification_no_match <- metagenome_before_classification[metagenome_before_classification$Prodigal.ID=="No match",]
metagenome_before_classification_match <- metagenome_before_classification[!(metagenome_before_classification$Prodigal.ID=="No match"),]
metagenome_hgt_plasmid_match<- ddply(metagenome_before_classification_match, .(Contig, Gene.AMR.Family, Host,Resistome,Antibiotic), function(x)x[which.min(x$Prodigal.Expected.value), ])
metagenome_after_classification <- rbind(metagenome_before_classification_no_match, metagenome_hgt_plasmid_match)




##14. CORRECT FOR GENES NOT MAKING SENSE TO PRESENT RESISTANCE AGAINST CERTAIN ANTIBIOTICS
metagenome_after_classification$type_resistance_antibiotic<- with(metagenome_after_classification, ifelse(grepl("CEF|CFD|FOX|MER|SUL", metagenome_after_classification$Antibiotic), "Beta-lactamase", ifelse(grepl("GEN|APS", metagenome_after_classification$Antibiotic),"Aminoglycoside", ifelse(grepl("MOX|DEL|GEP", metagenome_after_classification$Antibiotic),"Fluoroquinole",ifelse(grepl("DOX|ERA|OMA",metagenome_after_classification$Antibiotic),"Tetracycline", NA)))))

metagenome_after_classification$type_resistance_amino <- ifelse(grepl("minoglycoside", metagenome_after_classification$Card.Res.Gene.Drug.Class),"Aminoglycoside", NA)
metagenome_after_classification$type_resistance_beta_lact <- ifelse(grepl("eta-lact|cephalosporin|carbapene", metagenome_after_classification$Card.Res.Gene.Drug.Class)|grepl("eta-lact", metagenome_after_classification$Gene.AMR.Family),"Beta-lactamase", NA)
metagenome_after_classification$type_resistance_tetrac <- ifelse(grepl("etracyclin", metagenome_after_classification$Card.Res.Gene.Drug.Class),"Tetracycline", NA)
metagenome_after_classification$type_resistance_fluro <- ifelse(grepl("luoroqui", metagenome_after_classification$Card.Res.Gene.Drug.Class),"Fluoroquinole", NA)
metagenome_after_classification$type_resistance = paste(metagenome_after_classification$type_resistance_amino, metagenome_after_classification$type_resistance_beta_lact,metagenome_after_classification$type_resistance_tetrac,metagenome_after_classification$type_resistance_fluro, sep=", ")

metagenome_after_classification$resistance_true <-str_detect( metagenome_after_classification$type_resistance,metagenome_after_classification$type_resistance_antibiotic)

#if there is no gene annotation no need to correct
metagenome_after_classification_na <- metagenome_after_classification[is.na(metagenome_after_classification$Card.Res.Gene.Annotations),]

metagenome_after_classification_NO_na <- metagenome_after_classification[!is.na(metagenome_after_classification$Card.Res.Gene.Annotations),]
metagenome_after_classification_NO_na_1 <- metagenome_after_classification_NO_na[!(metagenome_after_classification_NO_na$resistance_true=="FALSE"),]

metagenome_1 <- rbind(metagenome_after_classification_na, metagenome_after_classification_NO_na_1)

# NEW LINE TO KEEP CONTIGS FOR SAMPLE WHERE NO ANTIBIOTICS WERE USED
index <- which(metagenome_after_classification$Antibiotic == "none")
metagenome_no_antibiotic <- metagenome_after_classification[index, ]
metagenome_1 <- rbind(metagenome_1, metagenome_no_antibiotic)

#check the contigs that have multiple genes if you miss anyone one of the of the same family on one contig
genes_contig <- ddply(metagenome_1,~Contig,summarise,distinct_res_genes=list(unique(na.omit(Universal.Gene.Name))))
metagenome_2 <-metagenome_1[!grepl("X58717",metagenome_1$Card.Res.Gene.Annotations),]
metagenome_3 <- metagenome_2[!(grepl("M21682",metagenome_2$Card.Res.Gene.Annotations) & grepl("gcaaactcgt_acatctaaga", metagenome_2$Contig)),]
length(unique(metagenome_1$Universal.Gene.Name))



##15. ADD THE IDs FROM CLUSTERING OF SEQUENCES
clusters <- read.csv("../../supplementary_table_6/Clustering of ARGs into Cluster IDs/Clustering of ARGs into Cluster IDs.csv",header=FALSE, sep="\t")
clusters$ID=as.numeric(factor(apply(clusters,1,function(x)toString(sort(x)))))
clusters_1 <- tidyr::separate_rows(clusters, V1,sep=",")
names(clusters_1)[1] <- "Prodigal.Query.ID"

cluster_ARGs <- left_join(metagenome_3,clusters_1, by = "Prodigal.Query.ID")


##16. CREATE A SINGLE COLUMN CONTANING SEQUENCE OF CONTIG IF ORF SEQUENCE IS MISSING
cluster_ARGs$Nucleotide_Sequence_Plasmid_Homologue <- ifelse(!(cluster_ARGs$Prodigal.ORF.sequence==""), cluster_ARGs$Prodigal.ORF.sequence, cluster_ARGs$Contig.Sequence)



##17. ANNOTATE IF THERE ARE HORIZONTALLY TRANSFER BASED ON THE MOBILE GENE POOL AND ALSO THE SPECIES THAT WERE INVOLED IN THIS TRANSFER
hgt_corrected$Nanopore_HGT_Blast.Expected.value <- as.numeric(hgt_corrected$Nanopore_HGT_Blast.Expected.value)
hgt_corrected$Nanopore_HGT_Blast.Percentage.of.identical.matches <- as.numeric((hgt_corrected$Nanopore_HGT_Blast.Percentage.of.identical.matches))
hgt_corrected$Nanopore_HGT_Blast.Query.coverage.per.HSP <- as.numeric((hgt_corrected$Nanopore_HGT_Blast.Query.coverage.per.HSP))

hgt_corrected <- hgt_corrected[,c(1,2,3,5,8,16,18,19)]
names(hgt_corrected)[c(1:4)] <- c("Downstream.Barcode.x","Upstream.Barcode.x","name.x","Prodigal.Query.ID")
metagenome_hgt <- right_join(hgt_corrected, cluster_ARGs, by = c("name.x", "Prodigal.Query.ID","Downstream.Barcode.x","Upstream.Barcode.x"))

hgt_IDs <- read.csv("../../databases/Mobile Gene Pool - Species per Event.tsv", sep="\t")

hgt_IDs$number_species_per_HGT <- str_count(hgt_IDs$Species, " ")
hgt_IDs <- hgt_IDs[c(1,2,5)]

names(hgt_IDs)[1] <- "Mobilepool.HGT"

metagenome_hgt_species <- left_join(metagenome_hgt, hgt_IDs, by = "Mobilepool.HGT")



##18. ANNONATE IF THEY HAVE PLASMID HOMOLOUES BASED ON THE PLDSB DATABASE
plasmid <- read.csv("../../supplementary_table_6/Plasmid Homologues/Plasmid Homologues.csv")
plasmid <- plasmid[c(1,13,28,33,35,37)]
names(plasmid)[c(1,2,4,5,6)] <- c("Taxon_Id no thr (plasmid database)", "taxon_species_name no thr (plasmid database)", "plasmid_database_evalue","plasmid_database_pident","plasmid_database_qcovhsp")
metagenome_hgt_plasmid <- left_join(metagenome_hgt_species, plasmid, by = "Prodigal.Query.ID")




##19. TAXONOMICALLY CLASSIFY THEM
taxonomy <- read.csv("../../supplementary_table_6/Taxonomy.csv")
taxonomy <- taxonomy[c(1,7,9,11,13:18)]
colnames(taxonomy)[2:4] <- paste("Taxonomy", colnames(taxonomy[,c(2:4)]), sep=".")

metagenome_final <- left_join(metagenome_hgt_plasmid, taxonomy, by = "Contig")

metagenome_final1 <- metagenome_final[!is.na(metagenome_final$ID),]

fname <- paste0("supp6v_nano", nano_thr, "_illu", illu_thr, ".tsv")
write.table(metagenome_final1, fname, row.names = FALSE, sep="\t")
