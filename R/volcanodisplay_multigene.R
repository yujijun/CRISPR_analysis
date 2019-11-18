#' This is for volcanodisplay of single gene

#' @param file_path input file with absolute path.
#' @param verfied_gene the gene name vector which you would like to label.
#' @param figure_name the figure title and output figure name.
#' @param output_path the output path which you would like to storage your figure.
#' @export
#' @return A figure of volcano

volcanodisplaymultigene <- function(file_path, verfied_gene, figure_name, output_path){
  source('./R/Volcano.R')
  require(MAGeCKFlute)
  #input data
  dd.rra = ReadRRA(file_path)
  head(dd.rra)
  dd.rra$Official = toupper(dd.rra$Official)
  if(sum(verfied_gene %in% dd.rra$Official) > 0){
    dd.rra$LFC = dd.rra$LFC
    dd.rra$color <- "non-verified"
    dd.rra$log_10 <- -log10(dd.rra$FDR)
    dd.rra$color[dd.rra$Official %in% verfied_gene] <- "verified"
    figure_title = figure_name
    subset = dd.rra[dd.rra$Official %in% verfied_gene,]
    #draw the plot
    p<-Volcano(data=dd.rra, x="LFC",y="log_10",
               label_data = subset,fill="color",
               color="color",label = "Official",color_palette = c("grey","red"),title=figure_title)  #+ geom_text_repel(data = subset,aes(x = LFC, y = FDR),label = Official)
    #save the plot
    figure_output = paste(figure_title,".png",sep = "")
    png(paste(output_path,figure_output,sep = "/"),height = 12,width = 17,units ="cm",res=150)
    print(p)
    dev.off()
  }else{
    print(paste("There isn't gene in ",figure_name))
  }
}
