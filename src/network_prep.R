# Network prep script

source("src/load_packages.R")
source("src/load_files.R")

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


# Select features shared in (ASO SPF/GF) & (SPF ASO/WT) Comparisons

genes_of_interest <- DEG_ST_ASO_GF_ST_ASO_SPF.SIG %>% 
  dplyr::filter(test_id %in% DEG_ST_ASO_SPF_ST_WT_SPF.SIG$test_id) %>% 
  dplyr::filter(test_id %ni% DEG_ST_WT_GF_ST_WT_SPF.SIG$test_id) %>% 
  dplyr::select(test_id, log2.fold_change.) %>% 
  dplyr::mutate(test_id = factor(test_id))

write.csv(genes_of_interest, file = "files/ASO_GF_vs_ASO_SPF_X_ASO_SPF_vs_WT_SPF.SIG_NetworkAnalyst_formatted.csv",
            row.names = F)


#### Only ASO_GF_vs_ASO_SPF.SIGgenes and Log2FC for network analysis
output2 <- DEG_ST_ASO_GF_ST_ASO_SPF.SIG %>% 
  select(test_id, log2.fold_change.) %>% 
  mutate(test_id = factor(test_id))
# write.csv(output2, file = "files/ASO_GF_vs_ASO_SPF.SIG_NetworkAnalyst_formatted.csv",
#             row.names = F)

#### Only ASO_SPF_vs_WT_SPF.SIGgenes and Log2FC for network analysis
output3 <- DEG_ST_ASO_SPF_ST_WT_SPF.SIG %>% 
  select(test_id, log2.fold_change.) %>% 
  mutate(test_id = factor(test_id))
# write.csv(output3, file = "files/DEG_ST_ASO_SPF_ST_WT_SPF.SIG_NetworkAnalyst_formatted.csv",
#           row.names = F)



