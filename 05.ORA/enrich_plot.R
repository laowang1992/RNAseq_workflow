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
  ego_df <- ego_df %>% 
    separate(GeneRatio, c("term_deg", "ann_deg"), sep = "/", convert = TRUE) %>%
    mutate(GeneRatio = term_deg / ann_deg)
  p <- ggplot(ego_df, aes(x = reorder(Description, Count), y = GeneRatio)) +
    geom_point(aes(size = Count, fill = -log10(p.adjust)), color = "#005298") +
    geom_point(aes(size = Count, fill = -log10(p.adjust)), shape=21, color = "#DCDCDC", stroke = 1.5) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    scale_size_continuous(range = c(3,11)) +
    labs(x = NULL, y = "Gene Ratio", fill = "-log10 P.adjust") +
    facet_grid(ONTOLOGY ~ ., scales = "free_y",space = "free_y") +
    coord_flip() +
    cowplot::theme_half_open() +  
    mytheme
  if (nrow(ego_df)*0.2+1.5 > 4) {
    GOheight <- nrow(ego_df)*0.2+1.5
  }else{
    GOheight <- 4
  }
  ggsave(p, filename = paste(prefix, "GO_enrichment.pdf", sep = "."), width = width, height = GOheight)
  ggsave(p, filename = paste(prefix, "GO_enrichment.png", sep = "."), width = width, height = GOheight, dpi = 500)
}

plotKEGG <- function(df, prefix, headn = 20, width = 7, height = nrow(ekp_df)*0.2+1.5){
  ekp_df <- df %>% dplyr::arrange(p.adjust) %>% head(n = headn)
  ekp_df <- ekp_df %>% 
    mutate(Description = if_else(str_length(Description)>50, paste(str_sub(Description, 1, 47), "...", sep = ""), Description))
  ekp_df <- ekp_df %>% 
    separate(GeneRatio, c("term_deg", "ann_deg"), sep = "/", convert = TRUE) %>% 
    mutate(GeneRatio = term_deg / ann_deg)
  
  p <- ggplot(ekp_df, aes(x = reorder(Description, Count), y = GeneRatio)) +
    geom_point(aes(size = Count, fill = -log10(p.adjust)), color = "#005298") +
    geom_point(aes(size = Count, fill = -log10(p.adjust)), shape=21, color = "#DCDCDC", stroke = 1.5) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    scale_size_continuous(range=c(3,10)) +
    labs(x = NULL, y = "Gene Ratio", fill = "-log10 P.adjust") +
    coord_flip() +
    cowplot::theme_half_open() +
    mytheme
  if (nrow(ekp_df)*0.2+1.5 > 4) {
    KEGGheight <- nrow(ekp_df)*0.2+1.5
  }else{
    KEGGheight <- 4
  }
  ggsave(p, filename = paste(prefix, "KEGG_enrichment.pdf", sep = "."), width = width, height = KEGGheight)
  ggsave(p, filename = paste(prefix, "KEGG_enrichment.png", sep = "."), width = width, height = KEGGheight, dpi = 500)
}

