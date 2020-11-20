# Network prep script

library(dplyr)
library(tidyverse)
library(readr)
library(tibble)
library(BiocGenerics)

# Read in data
df.fpkm <-read.csv(file = "files/FPKM_all.csv",  header= T)

datTraits <- df.fpkm %>%
  select(-Genes) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "sampleID") %>%
  select(sampleID) %>% 
  mutate(environment = if_else(grepl("GF", sampleID), "GF",
                               if_else(grepl("SPF", sampleID), "SPF", "ERR"))) %>%
  mutate(genotype = if_else(grepl("ASO", sampleID), "ASO",
                            if_else(grepl("WT", sampleID), "WT", "ERR"))) %>%
  mutate(group = paste(genotype, environment, sep = '_')) %>%
  dplyr::select(c(sampleID, environment, genotype, group)) 

df.output <- t(df.fpkm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleID") %>% 
  left_join(datTraits) %>% 
  column_to_rownames(var = "sampleID") %>% 
  t() %>% 
  as.data.frame() %>% 
  arrange(Genes)

df.output <- df.output[nrow(df.output):1, ]

colnames(df.output)[1] <- "#NAME"
df.output["group", 1] <- "#CLASS:group"
df.output["environment", 1] <- "#CLASS:environment"
df.output["genotype", 1] <- "#CLASS:genotype"
rownames(df.output) <- NULL

# write.table(df.output, file = "files/fpkm_NetworkAnalyst_formatted.txt",sep = "\t",
#             row.names = F, col.names = TRUE)



ASO_GF_vs_ASO_SPF <- read.delim(file = "files/DEG_ST_ASO_GF_ST_ASO_SPF.txt", header = TRUE, sep = "\t", dec = ".")
ASO_GF_vs_WT_GF <- read.delim(file = "files/DEG_ST_ASO_GF_ST_WT_GF.txt", header = TRUE, sep = "\t", dec = ".")
ASO_SPF_vs_WT_SPF <- read.delim(file = "files/DEG_ST_ASO_SPF_ST_WT_SPF.txt", header = TRUE, sep = "\t", dec = ".")

ASO_GF_vs_ASO_SPF.sig <- filter(ASO_GF_vs_ASO_SPF, significant == "yes")
# ASO_GF_vs_WT_GF.sig <- filter(ASO_GF_vs_WT_GF, significant == "yes")
ASO_SPF_vs_WT_SPF.sig <- filter(ASO_SPF_vs_WT_SPF, significant == "yes")

# Select features shared in (ASO SPF/GF) & (SPF ASO/WT) Comparisons
genes_of_interest <- ASO_GF_vs_ASO_SPF.sig %>% 
  filter(test_id %in% ASO_SPF_vs_WT_SPF.sig$test_id)



#### Only ASO_GF_vs_ASO_SPF.sig genes and Log2FC for network analysis
output2 <- ASO_GF_vs_ASO_SPF.sig %>% 
  select(test_id, log2.fold_change.) %>% 
  mutate(test_id = factor(test_id))
# write.csv(output2, file = "files/ASO_GF_vs_ASO_SPF.sig_NetworkAnalyst_formatted.csv",
#             row.names = F)

#### Only ASO_SPF_vs_WT_SPF.sig genes and Log2FC for network analysis
output3 <- ASO_SPF_vs_WT_SPF.sig %>% 
  select(test_id, log2.fold_change.) %>% 
  mutate(test_id = factor(test_id))
# write.csv(output3, file = "files/ASO_SPF_vs_WT_SPF.sig_NetworkAnalyst_formatted.csv",
#           row.names = F)

