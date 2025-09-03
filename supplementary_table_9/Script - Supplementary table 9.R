library(dplyr)
library(tidyr)
rm(list=ls())

##### PART 1 #####

# Prepare a combined supplementary table 6 "like" table which contains:
# - contigs that came up against the same antibiotic in run1/run2 and also in run3
# - contigs that came up in run3 without any antibiotics 

# read supplementary table 6 - run1 and run2
supp6o <- read.csv(
  "../supplementary_table_6/Supplementary table 6.tsv", sep = "\t")

supp6o_cr <- supp6o[, which(names(supp6o) %in% c(
  "Downstream.Barcode.x", "Upstream.Barcode.x", "Antibiotic"
))] %>% distinct()

# read supplementary table 6 "like" table - run3
supp6n <- read.csv("supp6v_nano10_illu9/supp6v_nano10_illu9.tsv", sep = "\t")

# only keep rows that came up with both run1/run2 and run3
# this will remove rows in run3 that came up without any antibiotics
# these rows will be added in a subsequent step
supp6nf <- dplyr::inner_join(supp6n, supp6o_cr)

# only keep rows where the resistance gene matches the antibiotic it came up against
supp6nf2 <- supp6nf[which(supp6nf$resistance_true == TRUE),]

# read supplementary table 6 "like" table - run3, but with illumina threshold 0
supp6n_no_ab <- read.csv("supp6v_nano10_illu0/supp6v_nano10_illu0.tsv", sep = "\t")

# filter to rows where no antibiotics were used
supp6n_no_ab <- supp6n_no_ab[which(supp6n_no_ab$Antibiotic == "none"),]

# bind these rows to supp6nf
supp6nf2 <- dplyr::bind_rows(supp6nf2, supp6n_no_ab)

#export table
#write.table(
#  supp6nf2,
#  file = "supp6nf2.tsv",
#  sep = "\t",
#  row.names = FALSE,
#  quote = FALSE
#)

##### PART 2 #####

# Generate dummy variables for rows where there are antibiotics.

# Generate a Table 6 with nano >= 10, illu >= 0 for no AB data and create dummy variables.

# Full Join the two dummy tables and handle NAs.

prep_dummies <- function(df) {
  dfcr <- df[,which(names(df) %in% c("Host", "Antibiotic", "ID"))]
  dfcr <- distinct(dfcr)
  
  count_hosts <- dfcr %>% 
    group_by(Antibiotic, ID) %>% 
    summarise(hosts = length(unique(Host)),
              EC = "EC" %in% Host,
              KP = "KP" %in% Host,
              SE = "SE" %in% Host,
              SS = "SS" %in% Host)
  dff <- dplyr::left_join(dfcr, count_hosts)
  dff <- dff[,-which(names(dff) == "Host")]
  dff <- distinct(dff)
  return(dff)
}

# dummies for rows with ABs
df <-supp6nf2[which(supp6nf2$Antibiotic != "none"),]
dfcr <- df[,which(names(df) %in% c("Host", "Antibiotic", "ID"))]
dfcr <- distinct(dfcr)

count_hosts <- dfcr %>% 
  group_by(Antibiotic, ID) %>% 
  summarise(hosts = length(unique(Host)),
            EC = "EC" %in% Host,
            KP = "KP" %in% Host,
            SE = "SE" %in% Host,
            SS = "SS" %in% Host)
dff <- dplyr::left_join(dfcr, count_hosts)
dff <- dff[,-which(names(dff) %in% c("Antibiotic","Host"))]
dff <- distinct(dff)
for (i in unique(dff$ID)) {
  index <- which(dff$ID == i)
  dff$EC[index] <- "TRUE" %in% dff$EC[index]
  dff$KP[index] <- "TRUE" %in% dff$KP[index]
  dff$SE[index] <- "TRUE" %in% dff$SE[index]
  dff$SS[index] <- "TRUE" %in% dff$SS[index]
}
dff$hosts <- dff$EC + dff$KP + dff$SE + dff$SS
dff <- distinct(dff)

n_ab_dummies <- dff
n_ab_dummies <- dplyr::rename(n_ab_dummies, EC_AB = EC)
n_ab_dummies <- dplyr::rename(n_ab_dummies, KP_AB = KP)
n_ab_dummies <- dplyr::rename(n_ab_dummies, SE_AB = SE)
n_ab_dummies <- dplyr::rename(n_ab_dummies, SS_AB = SS)
n_ab_dummies <- dplyr::rename(n_ab_dummies, hosts_AB = hosts)

# dummies for rows without ABs

n_no_ab_dummies <- prep_dummies(supp6n_no_ab[which(supp6n_no_ab$Antibiotic == "none"),])
n_no_ab_dummies <- dplyr::rename(n_no_ab_dummies, EC_noAB = EC)
n_no_ab_dummies <- dplyr::rename(n_no_ab_dummies, KP_noAB = KP)
n_no_ab_dummies <- dplyr::rename(n_no_ab_dummies, SE_noAB = SE)
n_no_ab_dummies <- dplyr::rename(n_no_ab_dummies, SS_noAB = SS)
n_no_ab_dummies <- dplyr::rename(n_no_ab_dummies, hosts_noAB = hosts)
n_no_ab_dummies <- n_no_ab_dummies[-which(names(n_no_ab_dummies) == "Antibiotic")]

##### PART 3 #####

# Generate dummies for rows where there are ABs.
# For the "without AB" comparison, import data table from new experiment.
# Full Join the two dummy tables and handle NAs.

# dummies for rows with ABs

supp6o <- read.csv("../supplementary_table_6/Supplementary table 6.tsv", sep = "\t")
supp6of <- supp6o[which(supp6o$resistance_true == TRUE),]

df <-supp6of[which(supp6of$Antibiotic != "none"),]
dfcr <- df[,which(names(df) %in% c("Host", "Antibiotic", "ID"))]
dfcr <- distinct(dfcr)

count_hosts <- dfcr %>% 
  group_by(Antibiotic, ID) %>% 
  summarise(hosts = length(unique(Host)),
            EC = "Ecoli" %in% Host,
            KP = "Klebsiella" %in% Host,
            SE = "Salmonella" %in% Host,
            SS = "Shigella" %in% Host)
dff <- dplyr::left_join(dfcr, count_hosts)
dff <- dff[,-which(names(dff) %in% c("Antibiotic","Host"))]
dff <- distinct(dff)
for (i in unique(dff$ID)) {
  index <- which(dff$ID == i)
  dff$EC[index] <- "TRUE" %in% dff$EC[index]
  dff$KP[index] <- "TRUE" %in% dff$KP[index]
  dff$SE[index] <- "TRUE" %in% dff$SE[index]
  dff$SS[index] <- "TRUE" %in% dff$SS[index]
}
dff$hosts <- dff$EC + dff$KP + dff$SE + dff$SS
dff <- distinct(dff)

n_ab_dummies <- dff
n_ab_dummies <- dplyr::rename(n_ab_dummies, EC_AB = EC)
n_ab_dummies <- dplyr::rename(n_ab_dummies, KP_AB = KP)
n_ab_dummies <- dplyr::rename(n_ab_dummies, SE_AB = SE)
n_ab_dummies <- dplyr::rename(n_ab_dummies, SS_AB = SS)
n_ab_dummies <- dplyr::rename(n_ab_dummies, hosts_AB = hosts)

# dummies for rows without ABs
# object n_no_ab_dummies comes from validation results
# n_no_ab_dummies <- readRDS("n_no_ab_dummies.rds")

# full join rows
n_dummies <- dplyr::full_join(n_ab_dummies, n_no_ab_dummies, by = "ID")

# handle NAs
n_dummies$Antibiotic <- replace_na(n_dummies$Antibiotic, "none")
n_dummies$hosts_AB <- replace_na(n_dummies$hosts_AB, 0)
n_dummies$EC_AB <- replace_na(n_dummies$EC_AB, FALSE)
n_dummies$KP_AB <- replace_na(n_dummies$KP_AB, FALSE)
n_dummies$SE_AB <- replace_na(n_dummies$SE_AB, FALSE)
n_dummies$SS_AB <- replace_na(n_dummies$SS_AB, FALSE)
n_dummies$hosts_noAB <- replace_na(n_dummies$hosts_noAB, 0)
n_dummies$EC_noAB <- replace_na(n_dummies$EC_noAB, FALSE)
n_dummies$KP_noAB <- replace_na(n_dummies$KP_noAB, FALSE)
n_dummies$SE_noAB <- replace_na(n_dummies$SE_noAB, FALSE)
n_dummies$SS_noAB <- replace_na(n_dummies$SS_noAB, FALSE)

#7. Read Table 7 created from updated Table 7 with Chryso's pipeline.

#8. Left join Table 7 with the dummies

supp7o <- read.csv("supp7f_nano10_illu9/supp7f_nano10_illu9.tsv", sep = "\t")
#supp7o <- read.csv("supp7of - nano10 illu9.tsv", sep = "\t")

supp7oj <- dplyr::left_join(supp7o, n_dummies)

# filter to rows where hosts_noAB == 4

supp7oj <- supp7oj[which(supp7oj$hosts_noAB == 4),]

write.table(
  supp7oj,
  file = "Supplementary table 9.csv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
