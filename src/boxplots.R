# Plotting 

library(dplyr)
library(tidyverse)
library(readr)
library(tibble)
library(BiocGenerics)
library(ggplot2)



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
#  Functions
#------------------------------------------------------------------------------
boxplot <- function(goi, cols = group.cols, 
                    title = blank.title, ylabel = blank.ylabel){
  
  blank.title = " "; blank.ylabel = " "
  group.cols = c("ASO_GF"= "#FF6000", "ASO_SPF" = "#077E97", 
                 "WT_GF" = "#A00000", "WT_SPF" = "808080")
  
  set.seed(123)
  
  df.long %>% 
    filter(Genes == goi) %>% 
    ggplot(aes(x=group, y=count)) +
    geom_boxplot(aes(fill = group), alpha = 0.75, outlier.alpha = 0, width = 0.9) +
    geom_point(aes(fill = group), position = position_jitterdodge(jitter.width = 0.5), 
               shape=21, size=1.5, alpha = 1) +
    theme_classic() +
    ggtitle(title) +
    labs(y = ylabel) +
    scale_color_manual(values = cols, name ="Group") +
    scale_fill_manual(values = cols, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none")
  
}

#------------------------------------------------------------------------------


goi <- "Znhit2"
boxplot(goi)

