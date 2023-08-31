parse_kegg_json <- function(json = "ko00001.json", organism = "ko") {
  pathway2name <- tibble::tibble(Pathway = character(), Pathway_Name = character(), Pathway_Class = character())
  ko2pathway <- tibble::tibble(Ko = character(), Pathway = character())
  
  kegg <- jsonlite::fromJSON(json)
  
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]
      
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        
        if(stringr::str_detect(pathway_info, paste("PATH:", organism, "[0-9]{5}", sep = ""))){
          pathway_id <- str_match(pathway_info, paste(organism, "[0-9]{5}", sep = ""))[1]
          pathway_name <- stringr::str_replace(pathway_info, paste(" \\[PATH:", organism, "[0-9]{5}\\]", sep = ""), "") %>% 
            stringr::str_replace("[0-9]{5} ", "")
          pathway2name <- rbind(pathway2name, tibble::tibble(Pathway = pathway_id, Pathway_Name = pathway_name, Pathway_Class = stringr::str_sub(A,7), Pathway_Subclass = stringr::str_sub(B, 7)))
          
          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
          
          kos <- stringr::str_match(kos_info, "K[0-9]{5}")[,1]
          
          ko2pathway <- rbind(ko2pathway, tibble::tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
        }
        next
      }
    }
  }
  pathway2name <- pathway2name %>% na.omit() %>% dplyr::distinct()
  ko2pathway <- ko2pathway %>% na.omit() %>% dplyr::distinct()
  kegg_pathway <- ko2pathway %>% dplyr::left_join(pathway2name, by = "Pathway")
  return(kegg_pathway)
  #save(pathway2name, ko2pathway, file = paste(organism, "kegg_info.RData", sep = "."))
}

create_orgdb <- function(gene_info, gene2go, genus, species){
  my_orgdb <- AnnotationForge::makeOrgPackage(gene_info  = gene_info,
                                              go         = gene2go,
                                              maintainer = 'wangpf <wangpf0608@126.com>',
                                              author     = 'Wang Pengfei',
                                              outputDir  = getwd(),
                                              tax_id     = 0000,
                                              genus      = genus,
                                              species    = species,
                                              goTable    = "go",
                                              version    = "1.0")
  pkgbuild::build(my_orgdb, dest_path = "./")
}

annoStat <- function(all_gene, gene_info, gene2go, gene2cog, gene2pathway){
  # 
  total_gene = length(all_gene)
  eggnog_anno = length(gene_info$GID)
  go_anno = length(unique(gene2go$GID))
  cog_anno = length(unique(gene2cog$GID))
  pathway_anno = length(unique(gene2pathway$GID))
  
  anno_stat <- tibble::tibble(
    database = c("EggNOG", "GO", "COG/KOG", "KEGG Pathway"),
    number = prettyNum(c(eggnog_anno, go_anno, cog_anno, pathway_anno), big.mark=",", scientific=FALSE),
    percentage = sprintf("%1.0f%%", round(x = c(eggnog_anno, go_anno, cog_anno, pathway_anno)*100/total_gene, digits = 0))
  )
  write.table(anno_stat, "anno_stat.txt", quote = F, row.names = F, sep = "\t")
  
  # GO statistics and plot --------------------------------------------------
  go_bp <- clusterProfiler::groupGO(gene     = all_gene,
                                    OrgDb    = org.Bnapus.eg.db,
                                    keyType  = "GID",
                                    ont      = "BP",
                                    level    = 2,
                                    readable = FALSE)
  
  go_bp <- as.data.frame(go_bp)
  go_bp$GO_Class <- "Biological Process"
  
  go_cc <- clusterProfiler::groupGO(gene     = all_gene,
                                    OrgDb    = org.Bnapus.eg.db,
                                    keyType  = "GID",
                                    ont      = "CC",
                                    level    = 2,
                                    readable = FALSE)
  
  go_cc <- as.data.frame(go_cc)
  go_cc$GO_Class <- "Cellular Component"
  
  go_mf <- clusterProfiler::groupGO(gene     = all_gene,
                                    OrgDb    = org.Bnapus.eg.db,
                                    keyType  = "GID",
                                    ont      = "MF",
                                    level    = 2,
                                    readable = FALSE)
  go_mf <- as.data.frame(go_mf)
  go_mf$GO_Class <- "Molecular Function"
  
  go_all <- rbind(go_bp, go_cc, go_mf)
  write.table(go_all, "go.txt", sep = "\t", quote = F)
  p <- ggplot2::ggplot(go_all) + 
    ggplot2::geom_bar(aes(x = Description, 
                          y = Count,
                          fill = GO_Class),
                      stat = "identity") + ggplot2::facet_grid(.~GO_Class, scales = "free_x", space = "free_x") + 
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(title = "GO function classification", y = "Number of genes") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(color = "black"),
                   legend.position = "none")
  ggplot2::ggsave("go.pdf", p, width = 13, height = 5.5)
  
  # Pathway statistics and plot ---------------------------------------------
  pathway_stat <- gene2pathway %>% dplyr::select(GID, Pathway_Class, Pathway_Subclass) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(Pathway_Class, Pathway_Subclass) %>%
    dplyr::summarise(Count = n(), Percentage = sprintf("%1.0f%%", round(x = n()*100/pathway_anno, digits = 0)))
  
  pathway_stat$Pathway_Subclass <- ordered(pathway_stat$Pathway_Subclass, levels = pathway_stat$Pathway_Subclass) 
  
  p <- ggplot2::ggplot(pathway_stat, ggplot2::aes(x = Pathway_Subclass, y = Count/pathway_anno)) +
    ggplot2::geom_bar(ggplot2::aes(fill = Pathway_Class), stat = 'identity') +
    ggplot2::geom_text(aes(label = Count), nudge_y = 0.05, size = 3.5) +
    ggplot2::scale_y_continuous(labels=scales::label_percent()) + 
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(y = "Percent of genes(%)", x ="", fill = "Class") +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(color = "black"))
  ggplot2::ggsave("pathway.pdf", p, width = 10, height = 5)
  write.table(gene2pathway, file = "pathway.txt", sep = "\t", quote = F)
  write.table(pathway_stat, file = "pathway_stat.txt", sep = "\t", quote = F, row.names = F)
  
  # COG statistics and plot -------------------------------------------------
  gene2cog$COG_Name = paste("(", gene2cog$COG, ")", gene2cog$COG_Name, sep = " ")
  
  write.table(gene2cog, file = "cog.txt", sep = "\t", quote = F, row.names = F)
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
  p <- ggplot2::ggplot(data = gene2cog) + 
    ggplot2::geom_bar(aes(x = COG, 
                          fill = COG_Name)) +
    ggplot2::labs(title = "COG/KOG Function Classification ", 
                  x = "",
                  y = "Number of genes") +
    ggplot2::scale_fill_manual(values = getPalette(length(gene2cog %>% pull(COG) %>% unique()))) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text = ggplot2::element_text(color = "black"),
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = ggplot2::unit(1,"line"),
                   legend.text = ggplot2::element_text(size = 7.5)) +
    ggplot2::guides(fill=ggplot2::guide_legend(ncol=1))
  ggplot2::ggsave("cog.pdf", p, width = 12, height = 6)
}
