# This script is essentially the same as the script preparing Supplementary
# Table 7 from Supplementary Table 6. The difference is that the input table is
# filtered to rows where resistance_true == TRUE.

library(tidyverse)
#if("dplyr" %in% (.packages())){
#  detach("package:dplyr", unload=TRUE) 
#  detach("package:plyr", unload=TRUE) 
#} 
library(plyr)
library(dplyr)
library(data.table)
library(reshape2)
library(foreign)
library(MASS)
rm(list = ls())

# read supplementary table 6
metagenome <- read.delim("../../supplementary_table_6/Supplementary table 6.tsv")
# filter to resistance_true == TRUE
metagenome <- metagenome[which(metagenome$resistance_true == TRUE),]

metagenome$evalue <- ifelse(is.na(metagenome$Prodigal.Query.coverage.per.HSP),metagenome$Blastx.evalue, metagenome$Prodigal.Expected.value)
metagenome$evalue <- as.numeric(metagenome$evalue)
metagenome$align_length <- ifelse(is.na(metagenome$Prodigal.Query.coverage.per.HSP),metagenome$Blastx.length, metagenome$Prodigal.Alignment.length)
metagenome$align_length <- as.numeric(metagenome$align_length)
metagenome$identity <- ifelse(is.na(metagenome$Prodigal.Query.coverage.per.HSP),metagenome$Blastx.pident, metagenome$Prodigal.Percentage.of.identical.matches)
metagenome$coverage <- ifelse(is.na(metagenome$Prodigal.Query.coverage.per.HSP),metagenome$Blastx.qcovhsp, metagenome$Prodigal.Query.coverage.per.HSP)



#1) clusters of genes
Genes <- ddply(metagenome, ~ID,summarise, list_ORF_homologues=list(unique(Prodigal.Query.ID)), Genes = list(unique(Universal.Gene.Name)))
Genes$list_ORF_homologues <- vapply(Genes$list_ORF_homologues,paste, collapse=", ", character(1L))
Genes$Genes <- vapply(Genes$Genes,paste, collapse=", ", character(1L))



#2) Find the best annotation for each cluster
annotation <- metagenome[c(73,4,60:63, 74, 91:94)]
annotation <- annotation[!duplicated(annotation),]

annotation_max <- ddply(annotation, .(ID), function(x)x[which.min(x$evalue), ])
colnames(annotation_max)[2:11] <- paste("Cluster.Best.Annotation", colnames(annotation_max[,c(2:11)]), sep = ".")

variables <- merge(Genes, annotation_max, by = "ID",all=TRUE)



#3 ) Thresholds for HGT transfer
#mobile gene pool
metagenome$Hgt_transferred_mobile_gene_pool_90 <- ifelse(metagenome$Nanopore_HGT_Blast.Percentage.of.identical.matches>90 & metagenome$Nanopore_HGT_Blast.Query.coverage.per.HSP>90, "Yes","No")

hgt_tra_mobile_gp_90 <- ddply(metagenome, .(ID), summarise, mobile_gene_pool_thr_90=list(unique(na.omit(Hgt_transferred_mobile_gene_pool_90))))
hgt_tra_mobile_gp_90$mobile_gene_pool_thr_90 <- vapply(hgt_tra_mobile_gp_90$mobile_gene_pool_thr_90, paste,collapse=", ", character(1L))
hgt_tra_mobile_gp_90$mobile_gene_pool_thr_90[hgt_tra_mobile_gp_90$mobile_gene_pool_thr_90==""] <- "No"

#plasmid
metagenome$Hgt_transferred_plasmid_90 <- ifelse(metagenome$plasmid_database_pident>90 & metagenome$plasmid_database_qcovhsp>90, "Yes","No")
hgt_tra_plasmid_90 <- ddply(metagenome, .(ID), summarise, plasmid_thr_90=list(unique(na.omit(Hgt_transferred_plasmid_90))))
hgt_tra_plasmid_90$plasmid_thr_90 <- vapply(hgt_tra_plasmid_90$plasmid_thr_90, paste,collapse=", ", character(1L))
hgt_tra_plasmid_90$plasmid_thr_90[hgt_tra_plasmid_90$plasmid_thr_90==""] <- "No"


m1.1 <- merge(hgt_tra_mobile_gp_90, hgt_tra_plasmid_90, by = "ID", all=TRUE)

m1.1[m1.1=="Yes, No"|m1.1=="No, Yes"] <- "Yes"

m1.1$HGT <- ifelse(m1.1$mobile_gene_pool_thr_90=="Yes"|m1.1$plasmid_thr_90=="Yes", "HGT","Non-HGT")

variables_1 <- left_join(variables, m1.1, by ="ID")





#5 ) Information on the IDs 
unique_genes_per_host <- aggregate(metagenome$Host ~ metagenome$ID, metagenome, function(i) 
  paste(unique(i), collapse = ",")) #165 res genes in metagenome data
names(unique_genes_per_host)[c(1,2)] <- c("Card.Res.Genes", "Host")

unique_fre_genes_presented_host <- ddply(metagenome,~ID,summarise,number_host=length(unique(Host)))
names(unique_fre_genes_presented_host)[c(1,2)] <- c("Card.Res.Genes", "Number_Hosts")

unique_sample <- aggregate(metagenome$Sample ~ metagenome$ID, metagenome, function(i) 
  paste(unique(i), collapse = ","))
names(unique_sample)[c(1,2)] <- c("Card.Res.Genes", "Sample")

unique_class_antibiotic <- aggregate(metagenome$Class ~ metagenome$ID, metagenome, function(i) 
  paste(unique(i), collapse = ","))
names(unique_class_antibiotic)[c(1,2)] <- c("Card.Res.Genes", "Class Antibiotic")

unique_resistome <- aggregate(metagenome$Resistome ~ metagenome$ID, metagenome, function(i) 
  paste(unique(i), collapse = ","))
names(unique_resistome)[c(1,2)] <- c("Card.Res.Genes", "Resistome")

unique_antibiotics_per_host <- aggregate(metagenome$Antibiotic ~ metagenome$ID, metagenome, function(i) 
  paste(unique(i), collapse = ",")) #172 res genes in metagenome data
names(unique_antibiotics_per_host)[c(1,2)] <- c("Card.Res.Genes", "Antibiotc")

unique_type_per_host <- aggregate(metagenome$Type ~ metagenome$ID, metagenome, function(i) 
  paste(unique(i), collapse = ",")) #172 res genes in metagenome data
names(unique_type_per_host)[c(1,2)] <- c("Card.Res.Genes", "Type Antibiotc")


m1 <- merge(unique_genes_per_host, unique_type_per_host, by = "Card.Res.Genes", all=TRUE)
m5 <- left_join(m1, unique_fre_genes_presented_host, by= "Card.Res.Genes", all= TRUE)
m6 <- left_join(m5, unique_sample, by = "Card.Res.Genes", all=TRUE)
m7 <- left_join(m6, unique_class_antibiotic, by = "Card.Res.Genes", all=TRUE)
m8 <- left_join(m7, unique_resistome, by = "Card.Res.Genes", all=TRUE)
m9 <- left_join(m8, unique_antibiotics_per_host, by = "Card.Res.Genes", all=TRUE)
names(m9)[1] <- "ID"

variables_2 <- left_join(variables_1,m9, by = "ID")

variables_2$meet_thr_evalue <-ifelse(variables_2$Cluster.Best.Annotation.evalue<10^-5&variables_2$Cluster.Best.Annotation.align_length>50, "Yes","No")





#6 ) add the taxonomic classifcation
metagenome_taxo_na <- metagenome[!is.na(metagenome$Taxonomy.evalue),]
metagenome_taxo_na$taxonomic_classification <- ifelse(metagenome_taxo_na$Taxonomy.evalue<10^-10, "Yes","No")
#metagenome_taxo_na$taxonomic_classification <- ifelse(metagenome_taxo_na$Taxonomy.pident>90&metagenome_taxo_na$Taxonomy.qcovhsp>90, "Yes","No")
metagenome_taxo_1 <- metagenome_taxo_na[metagenome_taxo_na$taxonomic_classification =="Yes",]

taxonomic_class <- ddply(metagenome_taxo_1, .(ID,Prodigal.Query.ID), summarise,taxonomy_thr=list(unique(na.omit(taxonomic_classification))), taxonomy_phylum=list(unique(na.omit(phylum))), taxonomy_species=list(unique(na.omit(species))))
taxonomic_class$taxonomy_thr <- vapply(taxonomic_class$taxonomy_thr, paste,collapse=", ", character(1L))
taxonomic_class$taxonomy_thr[taxonomic_class$taxonomy_thr==""] <- "No"
taxonomic_class$taxonomy_species <- vapply(taxonomic_class$taxonomy_species, paste,collapse=", ", character(1L))

taxonomic_class$taxonomy_phylum <- vapply(taxonomic_class$taxonomy_phylum, paste,collapse=", ", character(1L))

taxonomic_class <- taxonomic_class[-3]
names(taxonomic_class)[2]<- "Cluster.Best.Annotation.Prodigal.Query.ID"

variables_3 <- left_join(variables_2, taxonomic_class, by = c("ID","Cluster.Best.Annotation.Prodigal.Query.ID"))

variables_3$`Class Antibiotic` <- gsub("Cephalosporin,Carbapenem","Carbapenem,Cephalosporin",variables_3$`Class Antibiotic`)
variables_3$`Class Antibiotic` <- gsub("Tetracycline,Fluoroquinolone,Cephalosporin","Tetracycline,Cephalosporin,Fluoroquinolone",variables_3$`Class Antibiotic`)



#7) Add homologue for HGT in plasmid and HGT species
metagenome_HGT <- metagenome[c(73,4,75,78)]
metagenome_HGT <- metagenome_HGT[!duplicated(metagenome_HGT),]
names(metagenome_HGT)[2:4]<- c("Cluster.Best.Annotation.Prodigal.Query.ID", "Plasmid.Species","Mobile.Gene.Pool.Species")

variables_4 <- left_join(variables_3, metagenome_HGT, by = c("ID","Cluster.Best.Annotation.Prodigal.Query.ID"))
variables_4$Plasmid.Species <- ifelse(variables_4$plasmid_thr_90=="No",NA,variables_4$Plasmid.Species)
variables_4$Mobile.Gene.Pool.Species <- ifelse(variables_4$mobile_gene_pool_thr_90=="No",NA,variables_4$Mobile.Gene.Pool.Species)

variables_annotated <- variables_4[!is.na(variables_3$Cluster.Best.Annotation.evalue)&variables_3$meet_thr_evalue=="Yes",]

variables_genes <- variables_3[!is.na(variables_3$Cluster.Best.Annotation.evalue)&variables_3$meet_thr_evalue=="Yes",]

variables_genes$HGT_numerical <- ifelse(variables_genes$HGT== "HGT",1,0)
variables_genes$HGT <- gsub("HGT", "Mobile",variables_genes$HGT)
variables_genes$HGT <- gsub("Non-HGT", "Non-Mobile",variables_genes$HGT)



#9 ) Numeber of hosts in the different resistomes
metagenome_no_genes_na <- metagenome[!is.na(metagenome$Universal.Gene.Name),]
metagenome_no_genes_na <- metagenome_no_genes_na[metagenome_no_genes_na$evalue<10^-5&metagenome_no_genes_na$align_length>50,]

metagenome_no_genes_na$Resistome_new <- ifelse(!metagenome_no_genes_na$Resistome=="soil","human-related","enviromental")

resistome_N_hosts <- ddply(metagenome_no_genes_na, .(Resistome_new,ID), summarise, N_hosts=length(unique(na.omit(Host))))
resistome_N_hosts <- resistome_N_hosts[!resistome_N_hosts$N_hosts==0,]
resistome_N_hosts_1 <- spread(resistome_N_hosts, key=Resistome_new, value=N_hosts)
resistome_N_hosts_1[is.na(resistome_N_hosts_1)]<-0
variables_5 <- left_join(variables_4, resistome_N_hosts_1, by = "ID")

variables_6 <- variables_5[variables_5$meet_thr_evalue=="Yes"&!is.na(variables_5$meet_thr_evalue),]
write.table(variables_6, "supp7f_nano10_illu9.tsv", row.names = FALSE, sep="\t")
