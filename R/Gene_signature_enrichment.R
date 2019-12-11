#' This is for geneset enrichment analysis
#' @param geneset user input
#' @param LFC CRISPR LFC value
#' @param LFC_name CRIPPR gene name relate with LFC
#' @param figure_title figure title
#' @return enrichment matrix and enrichment plot
#' @export
geneset_enrichment <- function(geneset, LFC,LFC_name,figure_title){
  #create pathways information
  names(LFC) <- LFC_name
  LFC <- LFC[!is.na(LFC)]
  examplePathways <- list(depletion <- geneset)
  names(examplePathways) <- c("depletion")
  require(fgsea)
  require(ggplot2)
  fgseaRes <- fgsea(pathways = examplePathways,
                    stats = LFC,
                    minSize=1,
                    maxSize=250000,
                    nperm=100000)
  # plot the most significantly enriched pathway
  if(nrow(fgseaRes) == 0){
    fgseaRes <- "This is no enrichment score"
    p <- 0
  }else{
    p <- plotEnrichment(examplePathways[[1]],LFC) +
      labs(title=figure_title) +
      theme(plot.title = element_text(hjust = 0.5,face = "bold"))
  }
  out = list(enrich_score = fgseaRes, enrich_plot = p)
  return(out)
}



