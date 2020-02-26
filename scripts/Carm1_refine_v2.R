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
    scale_x_discrete(breaks=CARM1_no4$file_treatment,labels=labels_v1)+
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

CARM1 <- LFC[which(rownames(LFC) == genename),] #you can change any of the gene here
CARM1_matrix <- meta
CARM1_matrix$carm1 = as.numeric(CARM1[1,])
CARM1_matrix$file_treatment <- paste(meta$file_name_all,
                                     "_(",meta$Treatment_type,")",sep = "")
CARM1_matrix_NONA = CARM1_matrix %>% filter(!is.na(carm1)) #filter out the NA value
CARM1_matrix_NONA_choose <- CARM1_matrix_NONA %>% filter(Treatment_type == "Tcell_coculture" | Treatment_type == "Tcell")
#draw figure and output
#The all figure
CARM1_no4 <- CARM1_matrix_NONA_choose[-4,]
p <- singlegenebarplot(CARM1_no4,x = "file_treatment",y="carm1",fill = c(rep("#B2182B",1),rep("#2166AC",5)),
                       sort.by.groups = T,ylab = "LFC Value",xlab="Cohort Name",legend.title = "KO_type",title = paste("Barplot for", genename, "gene","in TumorKO(co-culture)",sep = " "))
show(p)


labels_v1=c("mNaiveTh(KO)_GATA3(sorting)","B16F10(KO)_OT-I(co_culture)",
         "MC38(KO)_OT-I(co_culture)","B16F10(KO)_Pmel-I(co_culture)",
         "hCD8T(KO)_CFSE(sorting)","IFNGR1mut.SKCM(KO)_MART-I(co_culture)")

labels = c("mNaiveTh(KO)_GATA3(sorting)","B16F10(KO)_OT-I(co_culture)",
           "MC38(KO)_OT-I(co_culture)","B16F10(KO)_Pmel-I(co_culture)",
           "hCD8T(KO)_CFSE(sorting)","IFNGR1mut.SKCM(KO)_MART-I(co_culture)")
