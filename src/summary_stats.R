# don't know what to call this yet

library(dplyr)
library(tidyverse)
library(readr)
library(tibble)
library(BiocGenerics)


## Plot gene expression correlation plots between groups of interest

# Read in data
df.fpkm <-read.csv(file = "files/FPKM_all.csv",  header= T)
# Make a gene key
df.fpkm2 <- df.fpkm %>%
  dplyr::mutate(GeneKey = paste0("Gene_", 1:nrow(df.fpkm)))
df.fpkm.key <- df.fpkm2 %>% 
  dplyr::select(Genes, GeneKey)
write.csv(df.fpkm.key, file = "files/group_medians_gene_key.csv",
          row.names = F)

df.fpkm3 <- df.fpkm2[1:500,]
rownames(df.fpkm3) <- NULL

df.plot <- df.fpkm2 %>%
  column_to_rownames(var = "GeneKey") %>%
  dplyr::select(-Genes) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleID") %>%
  dplyr::mutate(environment = if_else(grepl("GF", sampleID), "GF",
                               if_else(grepl("SPF", sampleID), "SPF", "ERR"))) %>%
  dplyr::mutate(genotype = if_else(grepl("ASO", sampleID), "ASO",
                            if_else(grepl("WT", sampleID), "WT", "ERR"))) %>%
  dplyr::mutate(group = paste(genotype, environment, sep = '_')) 

df.stats <- 
  df.plot %>% 
  group_by(group) %>% 
  dplyr::select(-c(sampleID, environment, genotype)) %>%
  dplyr::summarise(across(everything(), list(median))) %>% 
  pivot_longer(!group, names_to = "key", values_to = "count") %>% 
  dplyr::mutate(key =  str_sub(string = key, 1, -3))

write.csv(df.stats, file = "files/group_medians.csv",
          row.names = F)



# 
# ## dplyr::select only ASO_GF x ASO_SPF 
# ASO_GF_vs_ASO_SPF <- df.stats %>% 
#   filter(group == "ASO_GF" | group == "ASO_SPF") %>% 
#   pivot_wider(names_from = group, values_from = count)
# ## dplyr::select only ASO_SPF x WT_SPF 
# ASO_SPF_vs_WT_SPF  <- df.stats %>% 
#   filter(group == "ASO_SPF" | group == "WT_SPF") %>% 
#   pivot_wider(names_from = group, values_from = count)
# ## dplyr::select only ASO_GF x WT_GF 
# ASO_GF_vs_WT_GF <- df.stats %>% 
#   filter(group == "ASO_GF" | group == "WT_GF") %>% 
#   pivot_wider(names_from = group, values_from = count)
# ## dplyr::select only WT_GF x WT_SPF 
# WT_GF_vs_WT_SPF <- df.stats %>% 
#   filter(group == "WT_GF" | group == "WT_SPF") %>% 
#   pivot_wider(names_from = group, values_from = count)
# ## dplyr::select only WT_GF x ASO_SPF 
# WT_GF_vs_ASO_SPF <- df.stats %>% 
#   filter(group == "WT_GF" | group == "ASO_SPF") %>% 
#   pivot_wider(names_from = group, values_from = count)
# 
# 
# 
# #-----------------------------------------------------------------------------------
# # Plotting Function
# #-----------------------------------------------------------------------------------
# xysummary <- function(df, x, y, xgroup, ygroup, title){
#   
#   cols <- c("nodiff" = "#808080", "diff" = "#ff0000")
#   col.rims <- c("nodiff" = "#ffffff", "diff" = "#ff0000")
#   
#   ASO_GF_vs_ASO_SPF %>% 
#     dplyr::mutate(log2diff = log2((x + 1)/(y + 1))) %>% 
#     dplyr::mutate(goi_col = if_else(abs(log2diff) > 1, "diff", "nodiff")) %>% 
#     ggplot(aes(x=log2(x + 1), y=log2(y + 1))) +
#     geom_point(aes(fill = goi_col, color = goi_col),
#       shape=21, size=0.7, alpha = 0.8) +
#     theme_classic() +
#     geom_abline(intercept = 0, slope = 1) +
#     labs(x = paste0("expression (log2 FPKM + 1): ", xgroup),
#          y = paste0("expression (log2 FPKM + 1): ", ygroup),
#          title = title) +
#     scale_color_manual(values = col.rims, name ="Group") +
#     scale_fill_manual(values = cols, name ="Group") +
#     theme(legend.position = "none", 
#           plot.title = element_text(hjust = 0.5))
# }
# 
# #-----------------------------------------------------------------------------------
# 
# p1 <- xysummary(ASO_GF_vs_ASO_SPF, x = ASO_GF_vs_ASO_SPF$ASO_GF, y = ASO_GF_vs_ASO_SPF$ASO_SPF, 
#           xgroup = "ASO_GF", ygroup = "ASO_SPF", title = "ASO_GF vs ASO_SPF")
# 
# p2 <- xysummary(ASO_SPF_vs_WT_SPF, x = ASO_SPF_vs_WT_SPF$ASO_SPF, y = ASO_SPF_vs_WT_SPF$WT_SPF, 
#           xgroup = "ASO_SPF", ygroup = "WT_SPF", title = "ASO_SPF_vs_WT_SPF")
# 
# p3 <- xysummary(ASO_GF_vs_WT_GF, x = ASO_GF_vs_WT_GF$ASO_GF, y = ASO_GF_vs_WT_GF$WT_GF, 
#           xgroup = "ASO_GF", ygroup = "WT_GF", title = "ASO_GF_vs_WT_GF")
# 
# p4 <- xysummary(WT_GF_vs_WT_SPF, x = WT_GF_vs_WT_SPF$WT_GF, y = WT_GF_vs_WT_SPF$WT_SPF, 
#           xgroup = "WT_GF", ygroup = "WT_SPF", title = "WT_GF_vs_WT_SPF")
# 
# p5 <- xysummary(WT_GF_vs_ASO_SPF, x = WT_GF_vs_ASO_SPF$ASO_SPF, y = WT_GF_vs_ASO_SPF$WT_GF, 
#           xgroup = "ASO_SPF", ygroup = "WT_GF", title = "WT_GF_vs_ASO_SPF")
# 
# 
# 
# c1 <- cowplot::plot_grid(p1, p2, p3, p4, p5, nrow = 2, ncol = 3)
# ggsave(filename = "data/expression_comparisons_medianbased.png", 
#        height = 7, 
#        width = 11)


