# analysis functions


#------------------------------------------------------------------------------
#  Functions
#------------------------------------------------------------------------------
GO_bubble_plot <- function(df){
  
  df %>% 
    top_n(15, wt=-q.value.FDR.B.H) %>% 
    mutate(hitsPerc=Hit.Count.in.Query.List*100/Hit.Count.in.Genome) %>%
    arrange(hitsPerc) %>%
    mutate(Name=factor(Name, levels=Name)) %>%
    ggplot(aes(x=hitsPerc, 
               y=Name, 
               colour=q.value.FDR.B.H, 
               size=Hit.Count.in.Query.List)) +
    geom_point() +
    theme_bw()+
    scale_color_viridis_c(option = "cividis") +
    expand_limits(x=1) +
    labs(x="Hits (%)", y="", colour="FDR", size="Count")
}

#------------------------------------------------------------------------------

GO_barplot <- function(df, fill){
  
  df %>% 
    mutate(hitsPerc=Hit.Count.in.Query.List*100/Hit.Count.in.Genome) %>%
    arrange(q.value.FDR.B.H, desc(hitsPerc)) %>%
    slice_head(n = 15) %>% 
    mutate(yvar = paste(Name, ID)) %>% 
    ggplot(aes(x= -log10(q.value.FDR.B.H), 
               y=reorder(yvar, -q.value.FDR.B.H))) +
    geom_bar(stat = "identity", fill = fill, width = 0.5) + 
    labs(x = expression(paste(log[10], " (Adjusted P-value)")), 
         y = "") +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) 
}

#------------------------------------------------------------------------------


xysummary <- function(df, genekeydf, siglist, 
                      x, y, xgroup, ygroup, title, fill){
  
  cols <- c("nodiff" = "#808080", "diff" = fill)
  col.rims <- c("nodiff" = "#ffffff", "diff" = fill)
  
  # remerge genename
  df.merged <- left_join(df, genekeydf, by = "key")
  # select df of only significant genes
  sigfilter <- siglist %>% 
    top_n(20, wt = q_value)
  print(sigfilter$gene)
  
  df.merged %>% 
    dplyr::mutate(genelabels = if_else(Genes %in% sigfilter$gene, Genes, "")) %>% 
    dplyr::mutate(goi_col = if_else(Genes %in% DEG_ST_ASO_GF_ST_ASO_SPF.SIG$gene, "diff", "nodiff")) %>% 
    # mutate(goi_col = if_else(abs(log2diff) > 1, "diff", "nodiff")) %>% 
    mutate(log2diff = log2((x + 1)/(y + 1))) %>% 
    ggplot(aes(x=log2(x + 1), y=log2(y + 1))) +
    geom_point(aes(fill = goi_col, color = goi_col),
               shape=21, size=0.7, alpha = 0.8) +
    geom_text_repel(aes(label = genelabels), segment.alpha = 0.1, size = 2, force = 3) +
    theme_classic() +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = paste0("expression (log2 FPKM + 1): ", xgroup),
         y = paste0("expression (log2 FPKM + 1): ", ygroup),
         title = title) +
    scale_color_manual(values = col.rims, name ="Group") +
    scale_fill_manual(values = cols, name ="Group") +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5))
}


# xysummary <- function(df, x, y, xgroup, ygroup, title, fill){
#   
#   cols <- c("nodiff" = "#808080", "diff" = fill)
#   col.rims <- c("nodiff" = "#ffffff", "diff" = fill)
#   
#   df %>% 
#     mutate(log2diff = log2((x + 1)/(y + 1))) %>% 
#     mutate(goi_col = if_else(abs(log2diff) > 1, "diff", "nodiff")) %>% 
#     ggplot(aes(x=log2(x + 1), y=log2(y + 1))) +
#     geom_point(aes(fill = goi_col, color = goi_col),
#                shape=21, size=0.7, alpha = 0.8) +
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

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
