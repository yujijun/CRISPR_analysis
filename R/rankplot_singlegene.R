#' rankplot of single gene for mageck
#' @param gene_list gene signature by user(string)
#' @param rra_LFC LFC value of a RRA dataset
#' @param rra_genename gene name of a dataset
#' @param cutoff value for rank plot
#' @export
rankplot_singlegene <- function(gene_list, rra_LFC,rra_genename,cutoff = 1,figure_title){
  require(MAGeCKFlute)
  genelist <- rra_LFC
  names(genelist) <- rra_genename
  genelist <- sort(genelist,decreasing = T)
  names(genelist)[!names(genelist)==gene_list] <- NA
  top <- sum(genelist > 0)
  bottom <- sum(genelist <=0)
  p2 = RankView(genelist, top=top,bottom = bottom,cutoff = 1) +
    labs(x="LFC value",title = figure_title)

  plot(p2)
}
