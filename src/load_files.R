# Heatmap of significant features in various comparisons

library(eulerr)
library(dplyr)
library(tidyverse)
library(readr)
library(tibble)
library(BiocGenerics)


# Load significant genes in comparisons
DEG_ST_ASO_GF_ST_ASO_SPF <- 
  read.delim(file = "files/DEG_ST_ASO_GF_ST_ASO_SPF.txt", header = TRUE, sep = "\t", dec = ".")
DEG_ST_ASO_GF_ST_WT_GF <- 
  read.delim(file = "files/DEG_ST_ASO_GF_ST_WT_GF.txt", header = TRUE, sep = "\t", dec = ".") 
DEG_ST_ASO_SPF_ST_WT_SPF <- 
  read.delim(file = "files/DEG_ST_ASO_SPF_ST_WT_SPF.txt", header = TRUE, sep = "\t", dec = ".")
DEG_ST_WT_GF_ST_WT_SPF <- 
  read.delim(file = "files/DEG_ST_WT_GF_ST_WT_SPF.txt", header = TRUE, sep = "\t", dec = ".")
DEG_ST_ASO_SPF_ST_WT_GF <- 
  read.delim(file = "files/DEG_ST_ASO_SPF_ST_WT_GF.txt", header = TRUE, sep = "\t", dec = ".")

# Load significant genes in comparisons
DEG_ST_ASO_GF_ST_ASO_SPF.SIG <- 
  DEG_ST_ASO_GF_ST_ASO_SPF %>% 
  filter(significant == "yes")
DEG_ST_ASO_GF_ST_WT_GF.SIG <- 
  DEG_ST_ASO_GF_ST_WT_GF %>% 
  filter(significant == "yes")
DEG_ST_ASO_SPF_ST_WT_SPF.SIG <- 
  DEG_ST_ASO_SPF_ST_WT_SPF %>% 
  filter(significant == "yes")
DEG_ST_WT_GF_ST_WT_SPF.SIG <- 
  DEG_ST_WT_GF_ST_WT_SPF %>% 
  filter(significant == "yes")
DEG_ST_ASO_SPF_ST_WT_GF.SIG <- 
  DEG_ST_ASO_SPF_ST_WT_GF %>% 
  filter(significant == "yes")

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

# Load TOPGENE Data
ASO_GF_ST_ASO_SPF_DN <- read.delim(file = "files/TOPPGENE_DEG_ST_ASO_GF_ST_ASO_SPF_DN.txt", 
                                   header = TRUE, sep = "\t", dec = ".")
ASO_GF_ST_ASO_SPF_UP <- read.delim(file = "files/TOPPGENE_DEG_ST_ASO_GF_ST_ASO_SPF_UP.txt", 
                                   header = TRUE, sep = "\t", dec = ".")
ASO_GF_ST_WT_GF_DN <- read.delim(file = "files/TOPPGENE_DEG_ST_ASO_GF_ST_WT_GF_DN.txt", 
                                   header = TRUE, sep = "\t", dec = ".")
ASO_GF_ST_WT_GF_UP <- read.delim(file = "files/TOPPGENE_DEG_ST_ASO_GF_ST_WT_GF_UP.txt", 
                                   header = TRUE, sep = "\t", dec = ".")
ASO_SPF_ST_WT_SPF_DN <- read.delim(file = "files/TOPPGENE_DEG_ST_ASO_SPF_ST_WT_SPF_DN.txt", 
                                   header = TRUE, sep = "\t", dec = ".")
ASO_SPF_ST_WT_SPF_UP <- read.delim(file = "files/TOPPGENE_DEG_ST_ASO_SPF_ST_WT_SPF_UP.txt", 
                                   header = TRUE, sep = "\t", dec = ".")

