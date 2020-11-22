# Differential Expression & Functional Enrichment Plots

source("src/load_packages.R")
source("src/load_files.R")
source("src/functions.R")


df.stats <- read.csv(file = "files/group_medians.csv",  header= T)[-1]
gene.key <- read.csv(file = "files/group_medians_gene_key.csv",  header= T) %>% 
  dplyr::rename(key = GeneKey)

#----------------------------------------------------------------
# Functional Enrichment plots
#----------------------------------------------------------------


comparisons <- list(ASO_GF_ST_ASO_SPF_DN, ASO_GF_ST_ASO_SPF_UP, 
                 ASO_GF_ST_WT_GF_DN, ASO_GF_ST_WT_GF_UP,
                 ASO_SPF_ST_WT_SPF_DN, ASO_SPF_ST_WT_SPF_UP)
comp.names <- c("ASO_GF_ST_ASO_SPF_DN", "ASO_GF_ST_ASO_SPF_UP", 
  "ASO_GF_ST_WT_GF_DN", "ASO_GF_ST_WT_GF_UP",
  "ASO_SPF_ST_WT_SPF_DN", "ASO_SPF_ST_WT_SPF_UP")
names(comparisons) <- comp.names


# Loop to create Enrichment Plots
for (i in comp.names) {
  
  if (i == "ASO_GF_ST_ASO_SPF_DN"  | i == "ASO_GF_ST_ASO_SPF_UP") {
    fill = "#ff0000"
  } else if (i == "ASO_GF_ST_WT_GF_DN"  | i == "ASO_GF_ST_WT_GF_UP") {
    fill = "#7fb539"
  } else if (i == "ASO_SPF_ST_WT_SPF_DN"  | i == "ASO_SPF_ST_WT_SPF_UP") {
    fill = "#397fb5"
  }
    
  GOCC <- comparisons[[i]] %>%
    dplyr::filter(Category == "GO: Cellular Component")
  if (nrow(GOCC) > 0) {
    GOCC.plot <- GO_barplot(GOCC, fill = fill)
    ggsave(GOCC.plot, filename = paste0("data/Enrichment_Analysis/Cellular_Compartment/GO-CC_", i,".pdf"),
           height = 4, width = 7) 
  }

  GOBP <- comparisons[[i]] %>%
    dplyr::filter(Category == "GO: Biological Process")
  if (nrow(GOBP) > 0) {
    GOBP <- GO_barplot(GOBP, fill = fill)
    ggsave(GOBP, filename = paste0("data/Enrichment_Analysis/Biological_Processes/GO-BP_", i,".pdf"),
           height = 4, width = 7) 
  }
  
  Pathway <- comparisons[[i]] %>%
    dplyr::filter(Category == "Pathway")
  if (nrow(Pathway) > 0) {
    Pathway.plot <- GO_barplot(Pathway, fill = fill)
    ggsave(Pathway.plot, filename = paste0("data/Enrichment_Analysis/Pathways/GO-Pathway_", i,".pdf"),
           height = 4, width = 7) 
  }

  Disease <- comparisons[[i]] %>%
    dplyr::filter(Category == "Disease")
  if (nrow(Disease) > 0) {
    Disease.plot <- GO_barplot(Disease, fill = fill)
    ggsave(Disease.plot, filename = paste0("data/Enrichment_Analysis/Disease/GO-Disease_", i,".pdf"),
           height = 4, width = 7)
  }

  TFBS <- comparisons[[i]] %>%
    dplyr::filter(Category == "Transcription Factor Binding Site")
  if (nrow(TFBS) > 0) {
    TFBS.plot <- GO_barplot(TFBS, fill = fill)
    ggsave(TFBS.plot, filename = paste0("data/Enrichment_Analysis/Transcription_Factor_Binding_Site/GO-TFBS_", i,".pdf"),
           height = 4, width = 7) 
  }

  MP <- comparisons[[i]] %>%
    dplyr::filter(Category == "Mouse Phenotype")
  if (nrow(MP) > 0) {
    MP.plot <- GO_barplot(MP, fill = fill)
    ggsave(MP.plot, filename = paste0("data/Enrichment_Analysis/Mouse_Phenotype/GO-MP_", i,".pdf"),
           height = 4, width = 7) 
  }
}


#----------------------------------------------------------------
# Differential Expression plots
#----------------------------------------------------------------

## select only ASO_GF x ASO_SPF 
ASO_GF_vs_ASO_SPF <- df.stats %>% 
  dplyr::filter(group == "ASO_GF" | group == "ASO_SPF") %>% 
  pivot_wider(names_from = group, values_from = count)
## select only ASO_SPF x WT_SPF 
ASO_SPF_vs_WT_SPF  <- df.stats %>% 
  dplyr::filter(group == "ASO_SPF" | group == "WT_SPF") %>% 
  pivot_wider(names_from = group, values_from = count)
## select only ASO_GF x WT_GF 
ASO_GF_vs_WT_GF <- df.stats %>% 
  dplyr::filter(group == "ASO_GF" | group == "WT_GF") %>% 
  pivot_wider(names_from = group, values_from = count)
## select only WT_GF x WT_SPF 
WT_GF_vs_WT_SPF <- df.stats %>% 
  dplyr::filter(group == "WT_GF" | group == "WT_SPF") %>% 
  pivot_wider(names_from = group, values_from = count)
## select only WT_GF x ASO_SPF 
WT_GF_vs_ASO_SPF <- df.stats %>% 
  dplyr::filter(group == "WT_GF" | group == "ASO_SPF") %>% 
  pivot_wider(names_from = group, values_from = count)


p1 <- xysummary(df = DEG_ST_ASO_GF_ST_ASO_SPF, genekeydf = gene.key, siglist = DEG_ST_ASO_GF_ST_ASO_SPF.SIG,
                x = DEG_ST_ASO_GF_ST_ASO_SPF$value_1, y = DEG_ST_ASO_GF_ST_ASO_SPF$value_2, 
                xgroup = "ASO-GF", ygroup = "ASO-SPF", title = "ASO-GF vs. ASO-SPF", fill = "#ff0000", force = 1.5)
ggsave(p1, filename = "data/DEGs/expression_comparisons_ASO_GF_vs_ASO_SPF.pdf",
       height = 3.5,
       width = 3.66)

p2 <- xysummary(df = DEG_ST_ASO_SPF_ST_WT_SPF, genekeydf = gene.key, siglist = DEG_ST_ASO_SPF_ST_WT_SPF.SIG,
                x = DEG_ST_ASO_SPF_ST_WT_SPF$value_1, y = DEG_ST_ASO_SPF_ST_WT_SPF$value_2, 
                xgroup = "ASO-SPF", ygroup = "WT-SPF", title = "ASO-SPF vs. WT-SPF", fill = "#397fb5", force = 1)
ggsave(p2, filename = "data/DEGs/expression_comparisons_ASO_SPF_vs_WT_SPF.pdf", 
       height = 3.5, 
       width = 3.66)

p3 <- xysummary(df = DEG_ST_ASO_GF_ST_WT_GF, genekeydf = gene.key, siglist = DEG_ST_ASO_GF_ST_WT_GF.SIG,
                x = DEG_ST_ASO_GF_ST_WT_GF$value_1, y = DEG_ST_ASO_GF_ST_WT_GF$value_2, 
                xgroup = "ASO-GF", ygroup = "WT-GF", title = "ASO-GF vs. WT-GF", fill = "#7fb539", force = 1)
ggsave(p3, filename = "data/DEGs/expression_comparisons_ASO_GF_vs_WT_GF.pdf", 
       height = 3.5, 
       width = 3.66)

p4 <- xysummary(df = DEG_ST_WT_GF_ST_WT_SPF, genekeydf = gene.key, siglist = DEG_ST_WT_GF_ST_WT_SPF.SIG,
                x = DEG_ST_WT_GF_ST_WT_SPF$value_1, y = DEG_ST_WT_GF_ST_WT_SPF$value_2, 
                xgroup = "WT-GF", ygroup = "WT-SPF", title = "WT-GF vs. WT-SPF", fill = "#a9773d", force = 1)
ggsave(p4, filename = "data/DEGs/expression_comparisons_WT_GF_vs_WT_SPF.pdf", 
       height = 3.5, 
       width = 3.66)


p5 <- xysummary(df = DEG_ST_ASO_SPF_ST_WT_GF, genekeydf = gene.key, siglist = DEG_ST_ASO_SPF_ST_WT_GF.SIG,
                x = DEG_ST_ASO_SPF_ST_WT_GF$value_1, y = DEG_ST_ASO_SPF_ST_WT_GF$value_2, 
                xgroup = "ASO-SPF", ygroup = "WT-GF", title = "WT-GF vs. ASO-SPF", fill = "#8b64cc", force = 1.5)
ggsave(p5, filename = "data/DEGs/expression_comparisons_WT_GF_vs_ASO_SPF.pdf", 
       height = 3.5, 
       width = 3.66)


c1 <- cowplot::plot_grid(p1, p2, p3, p4, p5, nrow = 2, ncol = 3)
ggsave(c1, filename = "data/DEGs/expression_comparisons_medianbased.pdf", 
       height = 7, 
       width = 11)


