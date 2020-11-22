# Plotting 

library(dplyr)
library(tidyverse)
library(readr)
library(tibble)
library(BiocGenerics)
library(ggplot2)

source("src/functions.R")

# Read in data
df.fpkm <-read.csv(file = "files/FPKM_all.csv",  header= T)
# Create Group Column
df.long <- df.fpkm %>%
  pivot_longer(!Genes, names_to = "sampleID", values_to = "count") %>% 
  mutate(environment = if_else(grepl("GF", sampleID), "GF",
                               if_else(grepl("SPF", sampleID), "SPF", "ERR"))) %>%
  mutate(genotype = if_else(grepl("ASO", sampleID), "ASO",
                            if_else(grepl("WT", sampleID), "WT", "ERR"))) %>%
  mutate(group = paste(genotype, environment, sep = '_')) 


#------------------------------------------------------------------------------


goi <- "Hsbp1"

boxplot_fpkm(goi)
boxplot_transformed(goi)

