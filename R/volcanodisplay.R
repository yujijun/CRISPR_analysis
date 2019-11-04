#' This is for volcanodisplay
#' @param basic_path your basic path for input file, string.
#' @param specific_path specific_path for your input now, string.
#' @param file_name file name you would like to input, string.
#' @param verified_gene the gene name vector which you would like to label.
#' @param output_path the output path which you would like to storage your figure.
#' @export
#' @return A figure of volcano

volcanodisplay <- function(basic_path, specific_path, file_name, verified_gene, output_path){
  source('./R/Volcano.R')
  require(MAGeCKFlute)
  input_file = paste(paste(basic_path,specific_path,sep = ""),file_name,sep = "/")
  #input data
  dd.rra = ReadRRA(input_file)
  head(dd.rra)
  dd.rra$LFC = dd.rra$LFC
  dd.rra$color <- "non-verified"
  dd.rra$log_10 <- -log10(dd.rra$FDR)
  dd.rra$color[dd.rra$Official %in% verified_gene] <- "verified"
  figure_title = gsub(".gene_summary.txt","",basename(file_name))
  first_author = strsplit(specific_path,split = "_")[[1]][1]
  figure_title = paste(first_author,figure_title,sep = "_")
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
