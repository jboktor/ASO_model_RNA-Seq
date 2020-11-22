# Heatmap of significant gene in various comparisons

source("src/load_files.R")
source("src/load_packages.R")

## Create joined matrix of gene 
df.DEG_ST_ASO_GF_ST_ASO_SPF <- data.frame("gene"=DEG_ST_ASO_GF_ST_ASO_SPF.SIG$gene_id, "ASO_GF_ASO_SPF)" = TRUE)
df.DEG_ST_ASO_GF_ST_WT_GF <- data.frame("gene"=DEG_ST_ASO_GF_ST_WT_GF.SIG$gene_id, "ASO_GF_WT_GF" = TRUE)
df.DEG_ST_ASO_SPF_ST_WT_SPF <- data.frame("gene"=DEG_ST_ASO_SPF_ST_WT_SPF.SIG$gene_id, "ASO_SPF_WT_SPF" = TRUE)
df.DEG_ST_WT_GF_ST_WT_SPF <- data.frame("gene"=DEG_ST_WT_GF_ST_WT_SPF.SIG$gene_id, "WT_GF_WT_SPF" = TRUE)
df.DEG_ST_ASO_SPF_ST_WT_GF <- data.frame("gene"=DEG_ST_ASO_SPF_ST_WT_GF.SIG$gene_id, "ASO_SPF_WT_GF" = TRUE)

all_sig <- full_join(df.DEG_ST_ASO_GF_ST_ASO_SPF, df.DEG_ST_ASO_GF_ST_WT_GF,  by="gene") %>% 
  full_join(df.DEG_ST_ASO_SPF_ST_WT_SPF,  by="gene") %>% 
  full_join(df.DEG_ST_WT_GF_ST_WT_SPF,  by="gene") %>%
  full_join(df.DEG_ST_ASO_SPF_ST_WT_GF,  by="gene")
all_sig[is.na(all_sig)] <-  F


# Join DEG datatables - genes and log2fc only
A <- DEG_ST_ASO_GF_ST_ASO_SPF %>% 
  mutate(log2.fold_change. = -log2.fold_change.) %>% 
  dplyr::rename(ASO_SPF_vs_ASO_GF =  log2.fold_change.) %>% 
  select(c(gene, ASO_SPF_vs_ASO_GF))
B <- DEG_ST_ASO_GF_ST_WT_GF %>% 
  dplyr::rename(ASO_GF_vs_WT_GF =  log2.fold_change.) %>% 
  select(c(gene, ASO_GF_vs_WT_GF))
C <- DEG_ST_ASO_SPF_ST_WT_SPF %>% 
  dplyr::rename(ASO_SPF_vs_WT_SPF =  log2.fold_change.) %>% 
  select(c(gene, ASO_SPF_vs_WT_SPF))
D <- DEG_ST_WT_GF_ST_WT_SPF %>% 
  dplyr::rename(WT_GF_vs_WT_SPF =  log2.fold_change.) %>% 
  select(c(gene, WT_GF_vs_WT_SPF))
E <- DEG_ST_ASO_SPF_ST_WT_GF %>% 
  dplyr::rename(ASO_SPF_vs_WT_GF =  log2.fold_change.) %>% 
  select(c(gene, ASO_SPF_vs_WT_GF))

joined.df <- full_join(A, B,  by="gene") %>% 
  full_join(C,  by="gene") %>% 
  full_join(D,  by="gene") %>%
  full_join(E,  by="gene")



# filter for genes significant in any comparison
joined.df.sig <- joined.df %>% 
  filter(gene %in% all_sig$gene)

gene.dendro <- as.dendrogram(hclust(d = dist(x = joined.df.sig[-1])))
comparison.dendro <- as.dendrogram(hclust(d = dist(x = t(joined.df.sig[-1]))))

gene.dendro.plot <- ggdendrogram(data = fc.dendro, rotate = TRUE)
gene.dendro.plot <- gene.dendro.plot + theme(axis.text.y = element_blank())
# gene.dendro.plot
comparison.dendro.plot <- ggdendrogram(data = comparison.dendro, rotate = TRUE)
# comparison.dendro.plot

### REORDERING HEATMAP ROWS BASED ON DENDOGRAM
gene.order <- order.dendrogram(gene.dendro)
joined.df.sig$gene <- factor(x=joined.df.sig$gene, 
                     levels = joined.df.sig$gene[gene.order], ordered = TRUE)
comparison.order <- order.dendrogram(comparison.dendro)

# convert df to longer format
df.plot <- 
  joined.df.sig %>% 
  pivot_longer(!gene,
               names_to = "comparison",
               values_to = "log2fc")

df.plot$comparison <- factor(x=df.plot$comparison, 
                             levels = df.plot$comparison[comparison.order], ordered = TRUE)

# Plot
heatmap <- 
  df.plot %>% 
  ggplot(aes(x = comparison, y = gene, fill = log2fc)) +
  geom_tile() +
  coord_flip() +
  labs(y = "genes", fill = expression(paste(log[2], "[FC]"))) +
  scale_fill_distiller(palette = "RdBu") +
  # Note - this following list was rearranged to align with dendro order
  scale_x_discrete(limit = c("WT_GF_vs_WT_SPF", "ASO_GF_vs_WT_GF",
                             "ASO_SPF_vs_WT_GF", "ASO_SPF_vs_WT_SPF", "ASO_SPF_vs_ASO_GF"),
                   labels = c("WT-GF vs WT-SPF", "ASO-GF vs WT-GF",
                             "ASO-SPF vs WT-GF", "ASO-SPF vs WT-SPF", "ASO-SPF vs ASO-GF")) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank())

ggsave(heatmap, 
       filename = "data/significant_gene_heatmap.pdf", height = 3.5, width = 8)
