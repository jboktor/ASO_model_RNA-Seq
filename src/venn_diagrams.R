# Venn Diagram scripts

library(eulerr)
library(dplyr)
library(tidyverse)
library(readr)
library(tibble)
library(BiocGenerics)

# Load significant genes in comparisons
DEG_ST_ASO_GF_ST_ASO_SPF <- 
  read.delim(file = "files/DEG_ST_ASO_GF_ST_ASO_SPF.txt", header = TRUE, sep = "\t", dec = ".") %>% 
  filter(significant == "yes")
DEG_ST_ASO_GF_ST_WT_GF <- 
  read.delim(file = "files/DEG_ST_ASO_GF_ST_WT_GF.txt", header = TRUE, sep = "\t", dec = ".") %>% 
  filter(significant == "yes")
DEG_ST_ASO_SPF_ST_WT_SPF <- 
  read.delim(file = "files/DEG_ST_ASO_SPF_ST_WT_SPF.txt", header = TRUE, sep = "\t", dec = ".") %>% 
  filter(significant == "yes")
DEG_ST_WT_GF_ST_WT_SPF <- 
  read.delim(file = "files/DEG_ST_WT_GF_ST_WT_SPF.txt", header = TRUE, sep = "\t", dec = ".") %>% 
  filter(significant == "yes")
DEG_ST_ASO_SPF_ST_WT_GF <- 
  read.delim(file = "files/DEG_ST_ASO_SPF_ST_WT_GF.txt", header = TRUE, sep = "\t", dec = ".") %>% 
  filter(significant == "yes")

#-------------------------------------------------------------------------------------------
# Venn Diagram for all Significant features in comparisons of interest
#-------------------------------------------------------------------------------------------

## Create joined matrix of features 
df.DEG_ST_ASO_GF_ST_ASO_SPF <- data.frame("features"=DEG_ST_ASO_GF_ST_ASO_SPF$gene_id, "ASO_GF_vs_ASO_SPF" = TRUE)
df.DEG_ST_ASO_GF_ST_WT_GF <- data.frame("features"=DEG_ST_ASO_GF_ST_WT_GF$gene_id, "ASO_GF_vs_WT_GF" = TRUE)
df.DEG_ST_ASO_SPF_ST_WT_SPF <- data.frame("features"=DEG_ST_ASO_SPF_ST_WT_SPF$gene_id, "ASO_SPF_vs_WT_SPF" = TRUE)
df.DEG_ST_WT_GF_ST_WT_SPF <- data.frame("features"=DEG_ST_WT_GF_ST_WT_SPF$gene_id, "WT_GF_vs_WT_SPF" = TRUE)
df.DEG_ST_ASO_SPF_ST_WT_GF <- data.frame("features"=DEG_ST_ASO_SPF_ST_WT_GF$gene_id, "ASO_SPF_vs_WT_GF" = TRUE)

all_sig <- full_join(df.DEG_ST_ASO_GF_ST_ASO_SPF, df.DEG_ST_ASO_GF_ST_WT_GF,  by="features") %>% 
  full_join(df.DEG_ST_ASO_SPF_ST_WT_SPF,  by="features") %>% 
  # full_join(df.DEG_ST_WT_GF_ST_WT_SPF,  by="features") %>% 
  full_join(df.DEG_ST_ASO_SPF_ST_WT_GF,  by="features")
all_sig[is.na(all_sig)] <-  F

venn_ALL <- plot(euler(all_sig[-1], shape = "circle"), quantities = TRUE)
pdf(file = "data/VennDiagram_all_significant.pdf",
    width = 7, 
    height = 5,
    pointsize = 12)
# units = "in", pointsize = 12, res=300)
plot(venn_ALL)
dev.off()

#-------------------------------------------------------------------------------------------
# Venn Diagram for Shared UP & DN Significant features in comparisons of interest
#-------------------------------------------------------------------------------------------

# Load directional significance 
DEG_ST_ASO_GF_ST_ASO_SPF_UP <- 
  read.delim(file = "files/DEG_ST_ASO_GF_ST_ASO_SPF_UP.txt", header = F, sep = "\t", dec = ".") 
DEG_ST_ASO_GF_ST_ASO_SPF_DN <- 
  read.delim(file = "files/DEG_ST_ASO_GF_ST_ASO_SPF_DN.txt", header = F, sep = "\t", dec = ".") 

DEG_ST_ASO_GF_ST_WT_GF_UP <- 
  read.delim(file = "files/DEG_ST_ASO_GF_ST_WT_GF_UP.txt", header = F, sep = "\t", dec = ".") 
DEG_ST_ASO_GF_ST_WT_GF_DN <- 
  read.delim(file = "files/DEG_ST_ASO_GF_ST_WT_GF_DN.txt", header = F, sep = "\t", dec = ".") 

DEG_ST_ASO_SPF_ST_WT_SPF_UP <- 
  read.delim(file = "files/DEG_ST_ASO_SPF_ST_WT_SPF_UP.txt", header = F, sep = "\t", dec = ".") 
DEG_ST_ASO_SPF_ST_WT_SPF_DN <- 
  read.delim(file = "files/DEG_ST_ASO_SPF_ST_WT_SPF_DN.txt", header = F, sep = "\t", dec = ".") 

DEG_ST_WT_GF_ST_WT_SPF_UP <- 
  read.delim(file = "files/DEG_ST_WT_GF_ST_WT_SPF_UP.txt", header = F, sep = "\t", dec = ".")
DEG_ST_WT_GF_ST_WT_SPF_DN <- 
  read.delim(file = "files/DEG_ST_WT_GF_ST_WT_SPF_DN.txt", header = F, sep = "\t", dec = ".")

DEG_ST_ASO_SPF_ST_WT_GF_UP <- 
  read.delim(file = "files/DEG_ST_ASO_SPF_ST_WT_GF_UP.txt", header = F, sep = "\t", dec = ".")
DEG_ST_ASO_SPF_ST_WT_GF_DN <- 
  read.delim(file = "files/DEG_ST_ASO_SPF_ST_WT_GF_DN.txt", header = F, sep = "\t", dec = ".")

#-------------------------------------------------------------------------------------------
# UP Regulated Genes
#-------------------------------------------------------------------------------------------

df.DEG_ST_ASO_GF_ST_ASO_SPF_UP <- df.DEG_ST_ASO_GF_ST_ASO_SPF %>% 
  filter(features %in% DEG_ST_ASO_GF_ST_ASO_SPF_UP$V1)
df.DEG_ST_ASO_GF_ST_WT_GF_UP <- df.DEG_ST_ASO_GF_ST_WT_GF %>% 
  filter(features %in% DEG_ST_ASO_GF_ST_WT_GF_UP$V1)
df.DEG_ST_ASO_SPF_ST_WT_SPF_UP <- df.DEG_ST_ASO_SPF_ST_WT_SPF %>% 
  filter(features %in% DEG_ST_ASO_SPF_ST_WT_SPF_UP$V1)
df.DEG_ST_WT_GF_ST_WT_SPF_UP <- df.DEG_ST_WT_GF_ST_WT_SPF %>% 
  filter(features %in% DEG_ST_WT_GF_ST_WT_SPF_UP$V1)
df.DEG_ST_ASO_SPF_ST_WT_GF_UP <- df.DEG_ST_ASO_SPF_ST_WT_GF %>% 
  filter(features %in% DEG_ST_ASO_SPF_ST_WT_GF_UP$V1)

UP_sig <- full_join(df.DEG_ST_ASO_GF_ST_ASO_SPF_UP, df.DEG_ST_ASO_GF_ST_WT_GF_UP,  by="features") %>% 
  full_join(df.DEG_ST_ASO_SPF_ST_WT_SPF_UP,  by="features") %>% 
  full_join(df.DEG_ST_ASO_SPF_ST_WT_GF_UP,  by="features")
UP_sig[is.na(UP_sig)] <-  F

venn_UP <- plot(euler(UP_sig[-1], shape = "circle"), quantities = TRUE)

pdf(file = "data/VennDiagram_UP_significant.pdf",
    width = 7, 
    height = 5,
    pointsize = 12)
plot(venn_UP)
dev.off()


#-------------------------------------------------------------------------------------------
# Down Regulated Genes
#-------------------------------------------------------------------------------------------

df.DEG_ST_ASO_GF_ST_ASO_SPF_DN <- df.DEG_ST_ASO_GF_ST_ASO_SPF %>% 
  filter(features %in% DEG_ST_ASO_GF_ST_ASO_SPF_DN$V1)
df.DEG_ST_ASO_GF_ST_WT_GF_DN <- df.DEG_ST_ASO_GF_ST_WT_GF %>% 
  filter(features %in% DEG_ST_ASO_GF_ST_WT_GF_DN$V1)
df.DEG_ST_ASO_SPF_ST_WT_SPF_DN <- df.DEG_ST_ASO_SPF_ST_WT_SPF %>% 
  filter(features %in% DEG_ST_ASO_SPF_ST_WT_SPF_DN$V1)
df.DEG_ST_WT_GF_ST_WT_SPF_DN <- df.DEG_ST_WT_GF_ST_WT_SPF %>% 
  filter(features %in% DEG_ST_WT_GF_ST_WT_SPF_DN$V1)
df.DEG_ST_ASO_SPF_ST_WT_GF_DN <- df.DEG_ST_ASO_SPF_ST_WT_GF %>% 
  filter(features %in% DEG_ST_ASO_SPF_ST_WT_GF_DN$V1)

DN_sig <- full_join(df.DEG_ST_ASO_GF_ST_ASO_SPF_DN, df.DEG_ST_ASO_GF_ST_WT_GF_DN,  by="features") %>% 
  full_join(df.DEG_ST_ASO_SPF_ST_WT_SPF_DN,  by="features") %>% 
  full_join(df.DEG_ST_ASO_SPF_ST_WT_GF_DN,  by="features")
DN_sig[is.na(DN_sig)] <-  F

venn_DN <- plot(euler(DN_sig[-1], shape = "circle"), quantities = TRUE)

pdf(file = "data/VennDiagram_DN_significant.pdf",
    width = 7, 
    height = 5,
    pointsize = 12)
plot(venn_DN)
dev.off()



