#' GetKEGGPathways
#'
#' @param df data frame containing differential gene expression results
#' @param species species to map . Available option: "human" and "mouse".
#'
#' @return Data frame containing KEGG pathways
#' @export
GetKEGGPathways <- function(df = NULL, species = "human"){

  up <- df %>% filter(avg_log2FC > 0) %>% pull(gene)
  down <- df %>% filter(avg_log2FC < 0) %>% pull(gene)

  if(species == "human"){
    if (length(up) > 0)
      up <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, up, "ENTREZID", "SYMBOL")
    if (length(down) > 0)
      down <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, down, "ENTREZID", "SYMBOL")

    genes <- list(Up = up, Down = down)
    pathways <- limma::kegga(genes, species = "Hs")

  }else if(species == "mouse"){
    if (length(up) > 0)
      up <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, up, "ENTREZID", "SYMBOL")
    if (length(down) > 0)
      down <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, down, "ENTREZID", "SYMBOL")

    genes <- list(Up = up, Down = down)
    pathways <- limma::kegga(genes, species = "Mm")
  }

  # adjust p values
  pathways <- pathways %>% mutate(adj.P.Up = p.adjust(pathways[["P.Up"]], method = "bonferroni"))
  pathways <- pathways %>% mutate(adj.P.Down = p.adjust(pathways[["P.Down"]], method = "bonferroni"))

  # signed log p-values
  pathways <- pathways %>% mutate(logP.Up = -log10(adj.P.Up)) # -log to obtain positive values for UP
  pathways <- pathways %>% mutate(logP.Down = log10(adj.P.Down))# log to obtain negative values for UP

  return(pathways)
}
