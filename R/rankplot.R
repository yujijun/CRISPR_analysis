#' rankplot for mageck
#' @param gene_list gene signature by user
#' @param rra_LFC LFC value of a RRA dataset
#' @param rra_genename gene name of a dataset
#' @param cutoff value for rank plot
#' @export
rankplot <- function(gene_list, rra_LFC,rra_genename,cutoff = 1,figure_title){
  require(MAGeCKFlute)
  genelist <- rra_LFC
  names(genelist) <- rra_genename
  genelist <- sort(genelist,decreasing = T)
  genelist_new <- genelist
  names(genelist)[!(names(genelist) %in% gene_list)] <- NA
  genelist_new[!(names(genelist_new) %in% gene_list)] <- NA
  if (sum(genelist_new>cutoff,na.rm = T) > 10){
    genenames= names(rank(genelist_new,na.last = NA)[10])
    top = which(names(genelist_new) == genenames)
  }else{
    top = sum(genelist>cutoff,na.rm = T)
  }
  if(sum(genelist_new < -cutoff, na.rm = T) >10){
    genenames= names(rank(genelist_new,na.last = NA)[sum(!is.na(genelist_new)) - 10])
    bottom = which(names(genelist_new) == genenames)
    bottom = length(genelist_new) - bottom
  }else{
    bottom = sum(genelist < -cutoff, na.rm = T)
  }
  p2 = RankView(genelist, top=top,bottom = bottom,cutoff = 1) +
    labs(x="LFC value",title = figure_title)

  plot(p2)
}


