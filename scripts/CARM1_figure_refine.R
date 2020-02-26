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
                 fill=fill,
                 #fill = c(rep("#B2182B",2),rep("#2166AC",5)),           # change fill color by mpg_level
                 color = "black",            # Set bar border colors to white
                 palette = "jco",            # jco journal color palett. see ?ggpar
                 sort.val = "desc",          # Sort the value in descending order
                 sort.by.groups = sort.by.groups,     # Don't sort inside each group
                 x.text.angle = 0,
                 #y.text.angle =90,
                 xlab = xlab,          # Rotate vertically x axis texts
                 ylab = ylab,
                 legend.title = legend.title,
                 rotate = TRUE,
                 ggtheme = theme_bw(),
                 width = 0.8
  ) + scale_y_continuous(limits=c(-2.05, 0.1),breaks = seq(-2,0.5)) +
    scale_x_discrete(breaks=CARM1_no4$file_treatment,labels=c("mNaiveTh GATA3 (Henriksson et al.) ","B16F10 OT-I (Kearney et al.)",
                                                              "MC38 OT-I aPD1 (Kearney et al.)","B16F10 Pmel-I IFNg (Pan et al.)",
                                                              "hCD8T CFSE (Shifrut et al.)","IFNGR1mut.SKCM MART-I (Vredevoogd et al.)"))+
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_blank()) +
    theme(title = element_blank()) +
    #theme(title = element_text(face = "bold",size = 13,hjust = 0.5)) +
    #theme(axis.text = element_text(face = "bold")) +
    theme(legend.position = "none") +
    #scale_fill_manual("legend", values = c("#2166AC" = "#2166AC", "#B2182B" = "#B2182B"))
    theme(axis.text.y = element_text(face = "bold",size=10))+
    theme(plot.margin = margin(5,0.5,5,0.5,"cm")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(aspect.ratio = 1/2) +
    theme(axis.text.x = element_text(angle = 90))
#, axis.line = element_line(colour = "black")
  #p <- p + coord_cartesian(xlim = c(-4, 4))
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
  CARM1_matrix_NONA_choose <- CARM1_matrix_NONA %>% filter(Treatment_type == "Tcell_coculture" | Treatment_type == "Tcell")
  #draw figure and output
  p <- singlegenebarplot(CARM1_matrix_NONA_choose,x = "file_treatment",y="carm1",fill= "color",
                         sort.by.groups = T,ylab = "LFC Value",xlab="Cohort Name",legend.title = "KO_type",title = paste("Barplot for", genename, "gene",sep = " "))

  figure_name = paste(figure_title,".png",sep = "")
  png(paste(output_path,figure_name,sep = "/"),height = 12,width = 17,units ="cm",res=150)
  print(p)
  dev.off()
}

#barplot for immuno kill/ tumor kill coculture/tumor kill sorting
#The first little figure
fill = c(rep("#B2182B",3),rep("#2166AC",3))
CARM1_immuno <- CARM1_matrix_NONA %>% filter(normal_type == "Immune cell KO")
CARM1_immuno_Tcell <- CARM1_immuno %>% filter(Treatment_type == "Tcell")
p <- singlegenebarplot(CARM1_immuno_Tcell,x = "file_treatment",y="carm1",fill=c(rep("#2166AC",2)),
                       sort.by.groups = T,ylab = "LFC Value",xlab="Cohort Name",legend.title = "KO_type",title = paste("Barplot for", genename, "gene","in immunecell KO",sep = " "))
show(p)
#the second little figure
CARM1_tumor_coculture <- CARM1_matrix_NONA[grep("coculture", CARM1_matrix_NONA$Treatment_type),]
CARM1_Tcellcoculture <- CARM1_tumor_coculture %>% filter(Treatment_type == "Tcell_coculture")
p <- singlegenebarplot(CARM1_Tcellcoculture,x = "file_treatment",y="carm1",fill = c(rep("#B2182B",2),rep("#2166AC",3)),
                       sort.by.groups = T,ylab = "LFC Value",xlab="Cohort Name",legend.title = "KO_type",title = paste("Barplot for", genename, "gene","in TumorKO(co-culture)",sep = " "))
show(p)

#The all figure
CARM1_no4 <- CARM1_matrix_NONA_choose[-4,]
p <- singlegenebarplot(CARM1_no4,x = "file_treatment",y="carm1",fill = c(rep("#B2182B",1),rep("#2166AC",5)),
                       sort.by.groups = T,ylab = "LFC Value",xlab="Cohort Name",legend.title = "KO_type",title = paste("Barplot for", genename, "gene","in TumorKO(co-culture)",sep = " "))
show(p)



#CARM1_sorting <- CARM1_matrix_NONA[-grep("coculture", CARM1_matrix_NONA$Treatment_type),]
CARM1_sorting <- CARM1_matrix_NONA %>% filter(normal_type == "TumorKO_sorting")
p <- singlegenebarplot(CARM1_sorting,x = "file_treatment",y="carm1",fill = c(rep("#B2182B",1),rep("#2166AC",3)),
                       sort.by.groups = T,ylab = "LFC Value",xlab="Cohort Name",legend.title = "KO_type",title = paste("Barplot for", genename, "gene","in TumorKO(sorting)",sep = " "))

Barplot_singlegene(LFC,meta,genename,output_path,figure_title="Barplot for CARM1 in all cohort(maxscale)")

ggplot(CARM1_matrix_NONA_choose, aes(x = "file_treatment",y="carm1",fill= "color")) +
  geom_bar(stat = "identity")
