library(crisprvarified)
library(dplyr)
library(ggplot2)
#### LFC after z-score####
#input file
input_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/data/"
meta = read.table(paste(input_path,"QC_meta.txt",sep = "/"),sep ="\t",header = T)
LFC = read.table(paste(input_path,"LFC_maxscale.txt",sep = "/"))
genename = "CARM1"
output_path <- "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/CRISPR_output_plot/Barplot"
#create barplot generation function for single gene
singlegenebarplot <- function(inputdata, x,y,fill,sort.by.groups,ylab,xlab,legend.title,title){
  require("ggpubr")
  p <- ggbarplot(inputdata, x = x, y = y,
                 fill = fill,           # change fill color by mpg_level
                 color = "white",            # Set bar border colors to white
                 palette = "jco",            # jco journal color palett. see ?ggpar
                 sort.val = "desc",          # Sort the value in descending order
                 sort.by.groups = sort.by.groups,     # Don't sort inside each group
                 x.text.angle = 90,
                 xlab = xlab,          # Rotate vertically x axis texts
                 ylab = ylab,
                 legend.title = legend.title,
                 rotate = TRUE,
                 ggtheme = theme_bw(),
                 title = title
  ) +
    theme(legend.title = element_text(face = "bold")) +
    theme(title = element_text(face = "bold",size = 13,hjust = 0.5)) +
    theme(axis.text = element_text(face = "bold")) +
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(face = "bold",size=8))
  plot(p)
}

singlegenebarplot_new <- function(inputdata, x,y,fill,sort.by.groups,ylab,xlab,legend.title,title){
  require("ggpubr")
  p <- ggplot(data = inputdata,aes(x=x,y=y)) +
    geom_bar(aes(fill=y < 0),stat="identity",
                 color = "white"            # Set bar border colors to white
  ) +
    scale_fill_manual(guide = FALSE, breaks=c(TRUE,FALSE),values=c("#B2182B","#2166AC")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(title = element_text(face = "bold",size = 20,hjust = 0.5)) +
    theme(axis.text = element_text(face = "bold")) +
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(face = "bold",size=15))
  plot(p)
}
Barplot_singlegene <- function(LFC,meta,genename,output_path,figure_title){
  # LFC is a dataframe include all samples' LFC in different genes;
  # meta is a dataframe include all meta information for all samples;
  # genename is a character of gene name which you woule like to perform.
  CARM1 <- LFC[which(rownames(LFC) == genename),] #you can change any of the gene here
  CARM1_matrix <- meta
  CARM1_matrix$carm1 = as.numeric(CARM1[1,])
  CARM1_matrix$file_treatment <- paste(meta$file_name_all,
                                       "_(",meta$Treatment_type,")",sep = "")
  CARM1_matrix_NONA = CARM1_matrix %>% filter(!is.na(carm1)) #filter out the NA value
  #draw figure and output
  p <- singlegenebarplot(CARM1_matrix_NONA,x = "file_treatment",y="carm1",fill = "KO_type",
                         sort.by.groups = T,ylab = "LFC Value",xlab="Cohort Name",legend.title = "KO_type",title = paste("Barplot for", genename, "gene",sep = " "))
  figure_name = paste(figure_title,".png",sep = "")
  png(paste(output_path,figure_name,sep = "/"),height = 12,width = 17,units ="cm",res=150)
  print(p)
  dev.off()
}

#barplot for immuno kill/ tumor kill coculture/tumor kill sorting
CARM1_immuno <- CARM1_matrix_NONA %>% filter(normal_type == "Immune cell KO")
p <- singlegenebarplot(CARM1_immuno,x = "file_treatment",y="carm1",fill = c(rep("#B2182B",3),rep("#2166AC",3)),
                       sort.by.groups = T,ylab = "LFC Value",xlab="Cohort Name",legend.title = "KO_type",title = paste("Barplot for", genename, "gene","in immunecell KO",sep = " "))

CARM1_tumor_coculture <- CARM1_matrix_NONA[grep("coculture", CARM1_matrix_NONA$Treatment_type),]
p <- singlegenebarplot(CARM1_tumor_coculture,x = "file_treatment",y="carm1",fill = c(rep("#B2182B",3),rep("#2166AC",4)),
                       sort.by.groups = T,ylab = "LFC Value",xlab="Cohort Name",legend.title = "KO_type",title = paste("Barplot for", genename, "gene","in TumorKO(co-culture)",sep = " "))

#CARM1_sorting <- CARM1_matrix_NONA[-grep("coculture", CARM1_matrix_NONA$Treatment_type),]
CARM1_sorting <- CARM1_matrix_NONA %>% filter(normal_type == "TumorKO_sorting")
p <- singlegenebarplot(CARM1_sorting,x = "file_treatment",y="carm1",fill = c(rep("#B2182B",1),rep("#2166AC",3)),
                       sort.by.groups = T,ylab = "LFC Value",xlab="Cohort Name",legend.title = "KO_type",title = paste("Barplot for", genename, "gene","in TumorKO(sorting)",sep = " "))

#barplot for all cohort
Barplot_singlegene(LFC,meta,genename,output_path,figure_title="Barplot for CARM1 in all cohort(maxscale)")

#barplot for Tcell coculture in type of tumor KO
meta_tcellco <- meta %>% filter(Treatment_type == "Tcell_coculture")
colnames(LFC) <- gsub("_LFC","",colnames(LFC))
LFC_tcellco <- LFC[,meta_tcellco$file_name_all]
figure_title <- "Barplot for CARM1 in Tcell coculture(maxscale)"
Barplot_singlegene(LFC_tcellco,meta_tcellco,genename,output_path,figure_title)

#barplot for Tcell in type of Immunecell KO
meta_tcell <- meta %>% filter(Treatment_type == "Tcell")
colnames(LFC) <- gsub("_LFC","",colnames(LFC))
LFC_tcell <- LFC[,meta_tcell$file_name_all]
figure_title <- "Barplot for CARM1 in Tcell(maxscale)"
Barplot_singlegene(LFC_tcell,meta_tcell,genename,output_path,figure_title)

####LFC before Z-score####
#input file
input_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/data/"
meta = read.table(paste(input_path,"QC_meta.txt",sep = "/"),header = T)
LFC_original = read.table(paste(input_path,"LFC_original.txt",sep = "/"))
genename = "CARM1"
Barplot_singlegene(LFC_original,meta,genename,output_path,figure_title="Barplot for CARM1 in all cohort(before z-score)")
#barplot for Tcell coculture in type of tumor KO(before z-score)
meta_tcellco <- meta %>% filter(Treatment_type == "Tcell_coculture")
colnames(LFC) <- gsub("_LFC","",colnames(LFC))
LFC_tcellco <- LFC_original[,meta_tcellco$file_name_all]
figure_title <- "Barplot for CARM1 in Tcell coculture(before z-score)"
Barplot_singlegene(LFC_tcellco,meta_tcellco,genename,output_path,figure_title)
#barplot for Tcell in type of Immunecell KO (before z-score)
meta_tcell <- meta %>% filter(Treatment_type == "Tcell")
colnames(LFC) <- gsub("_LFC","",colnames(LFC))
LFC_tcell <- LFC_original[,meta_tcell$file_name_all]
figure_title <- "Barplot for CARM1 in Tcell(before z-score)"
Barplot_singlegene(LFC_tcell,meta_tcell,genename,output_path,figure_title)

####LFC after Z-score but withour center####
#input file
input_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/data/"
meta = read.table(paste(input_path,"QC_meta.txt",sep = "/"),header = T)
LFC_scale = read.table(paste(input_path,"LFC_scale(nocenter).txt",sep = "/"))
genename = "CARM1"
Barplot_singlegene(LFC_scale,meta,genename,output_path,figure_title="Barplot for CARM1 in all cohort(just scale no center)")
#barplot for Tcell coculture in type of tumor KO(Just scale no center)
meta_tcellco <- meta %>% filter(Treatment_type == "Tcell_coculture")
colnames(LFC_scale) <- gsub("_LFC","",colnames(LFC_scale))
LFC_tcellco <- LFC_scale[,meta_tcellco$file_name_all]
figure_title <- "Barplot for CARM1 in Tcell coculture(Just scale no center)"
Barplot_singlegene(LFC_tcellco,meta_tcellco,genename,output_path,figure_title)

#barplot for Tcell in type of Immunecell KO (before z-score)
meta_tcell <- meta %>% filter(Treatment_type == "Tcell")
colnames(LFC_scale) <- gsub("_LFC","",colnames(LFC_scale))
LFC_tcell <- LFC_scale[,meta_tcell$file_name_all]
figure_title <- "Barplot for CARM1 in Tcell(Just scale no center)"
Barplot_singlegene(LFC_tcell,meta_tcell,genename,output_path,figure_title)


