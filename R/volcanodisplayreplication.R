#' This is for volcano replication display
#' @param dd.rra input dataframe which you have calculated the mean value.
#' @param figure_title figure title
#' @param verified_gene the gene name vector which you would like to label.
#' @param output_path the output path which you would like to storage your figure.
#' @return A figure of volcano
#' @export
volcanodisplayrep <- function(dd.rra, figure_title, verified_gene, output_path){
  require(MAGeCKFlute)
  dd.rra$LFC = dd.rra$LFC
  dd.rra$color <- "non-verified"
  dd.rra$log_10 <- -log10(dd.rra$FDR)
  dd.rra$color[dd.rra$Official %in% verified_gene] <- "verified"
  subset = dd.rra[dd.rra$Official %in% verified_gene,]

  #draw the plot
  p<-Volcano(data=dd.rra, x="LFC",y="log_10",
             label_data = subset,fill="color",
             color="color",label = "Official",color_palette = c("grey","red"),title=figure_title)  #+ geom_text_repel(data = subset,aes(x = LFC, y = FDR),label = Official)
  #save the plot
  figure_name = paste(figure_title,".png",sep = "")
  png(paste(output_path,figure_name,sep = "/"),height = 12,width = 17,units ="cm",res=150)
  print(p)
  dev.off()
}
