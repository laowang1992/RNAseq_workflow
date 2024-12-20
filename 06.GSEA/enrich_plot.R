mytheme <- theme(panel.grid = element_blank(),
                 axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 10,
                                            face = "plain", color = "black"),
                 axis.text.y = element_text(face = "plain", color = "black", size = 10),
                 axis.title = element_text(face = "bold", size = 12),
                 legend.text = element_text(size = 9, face = "bold"),
                 legend.title = element_text(size = 12),
                 strip.background = element_rect(fill = "#008B8B"))

plotGO <- function(df, prefix, headn = 15, width = 7, height = nrow(ego_df)*0.2+1.5){
  bp <- df %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::arrange(p.adjust) %>% head(n = headn)
  cc <- df %>% dplyr::filter(ONTOLOGY == "CC") %>% dplyr::arrange(p.adjust) %>% head(n = headn)
  mf <- df %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::arrange(p.adjust) %>% head(n = headn)
  ego_df <- rbind(bp, cc, mf)
  ego_df <- ego_df %>% 
    mutate(Description = if_else(str_length(Description)>50, paste(str_sub(Description, 1, 47), "...", sep = ""), Description))
  
  p <- ggplot(ego_df,
              aes(NES, fct_reorder(Description, NES), fill = -log10(pvalue))) +
    geom_col() +
    geom_segment(mapping = aes(x = min(NES, 0),
                               xend = ifelse(sign(NES) > 0, 0, NES),
                               yend = Description),
                 color = "gray", linetype = "dashed", linewidth = 1) +
    scale_x_continuous(expand = c(0.02 , 0.02)) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    cowplot::theme_half_open() +
    mytheme +
    labs(x = "Normalized Enrichment Score", y = NULL, fill = "-log10 P.adjust") +
    facet_grid(ONTOLOGY ~ ., scales = "free_y",space = "free_y")
  if (nrow(ego_df)*0.2+1.5 > 4) {
    GOheight <- nrow(ego_df)*0.2+1.5
  }else{
    GOheight <- 4
  }
  ggsave(p, filename = paste(prefix, "gseGO_enrichment.pdf", sep = "."), width = width, height = GOheight)
  ggsave(p, filename = paste(prefix, "gseGO_enrichment.png", sep = "."), width = width, height = GOheight, dpi = 500)
}

plotKEGG <- function(df, prefix, headn = 20, width = 7, height = nrow(ekp_df)*0.2+1.5){
  ekp_df <- df %>% dplyr::arrange(p.adjust) %>% head(n = headn)
  ekp_df <- ekp_df %>% 
    mutate(Description = if_else(str_length(Description)>50, paste(str_sub(Description, 1, 47), "...", sep = ""), Description))
  
  p <- ggplot(ekp_df,
              aes(NES, fct_reorder(Description, NES), fill = -log10(pvalue))) +
    geom_col() +
    geom_segment(mapping = aes(x = min(NES, 0),
                               xend = ifelse(sign(NES) > 0, 0, NES),
                               yend = Description),
                 color = "gray", linetype = "dashed", linewidth = 1) +
    scale_x_continuous(expand = c(0.02 , 0.02)) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    cowplot::theme_half_open() +
    mytheme +
    labs(x = "Normalized Enrichment Score", y = NULL, fill = "-log10 P.adjust")
  
  if (nrow(ekp_df)*0.2+1.5 > 4) {
    KEGGheight <- nrow(ekp_df)*0.2+1.5
  }else{
    KEGGheight <- 4
  }
  ggsave(p, filename = paste(prefix, "gseKEGG_enrichment.pdf", sep = "."), width = width, height = KEGGheight)
  ggsave(p, filename = paste(prefix, "gseKEGG_enrichment.png", sep = "."), width = width, height = KEGGheight, dpi = 500)
}

