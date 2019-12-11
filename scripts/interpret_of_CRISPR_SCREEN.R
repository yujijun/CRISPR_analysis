# --------------
# Date:  2019-11-18 15:29:42
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: This is for heatmap
#
# """
# Unweighted signature:
# 1. selected the top and bottom (maybe 10) genes from
# each type of screens (LFC>1.5).
# 2. heatmap with clustering and label the
# screen type like the QC heatmap (If there are too many genes).
# 3. boxplot: each box represents one gene, and color the box by
# screen type. (limited genes).
# 4.volcano plot showing top 10 genes
# (from the list) with highest logFC
# """
library(dplyr)
library(stringi)
library(stringr)
library(tibble)
library(crisprvarified)
library(MAGeCKFlute)
#gene list
basic_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/data/"
gene_list = read.table(paste(basic_path,"DGElist_zexian.txt",sep = "/"))
gene_list = gene_list[!duplicated(gene_list$V1),]
input_meta = read.table(paste(basic_path,"QC_meta.txt",sep = "/"),header = T)
input_LFC = read.table(paste(basic_path,"LFC_maxscale.txt",sep = "/"))
input_LFC_genesignature <- input_LFC %>% filter(rownames(input_LFC) %in% as.vector(gene_list$V1))


#sort out all meta data
#change the file name all


####split the gene signature by tumor_KO and Immune_KO####
input_LFC_tumorKO = input_LFC_genesignature[,which(input_meta$KO_type == "Tumor_KO")]
input_LFC_ImmuneKO = input_LFC_genesignature[,-which(input_meta$KO_type == "Tumor_KO")]

input_LFC_tumorKO = input_LFC_tumorKO[rowMeans(is.na(input_LFC_tumorKO)) <= 0.5,]
input_LFC_ImmuneKO = input_LFC_ImmuneKO[rowMeans(is.na(input_LFC_ImmuneKO)) <= 0.5,]
tumorKO_mean <- rowMeans(input_LFC_tumorKO,na.rm = T)
immuneKO_mean <- rowMeans(input_LFC_ImmuneKO,na.rm = T)
#top 10 and bottom 10
cutoff_value = 10
Top10_gene_tumorKO <- names(tumorKO_mean)[order(tumorKO_mean,decreasing = T)[1:cutoff_value]]
Bottom10_gene_tumorKO <- names(tumorKO_mean)[order(tumorKO_mean,decreasing = F)[1:cutoff_value]]
Top10_gene_immuneKO <- names(immuneKO_mean)[order(immuneKO_mean,decreasing = T)[1:cutoff_value]]
Bottom10_gene_immuneKO <- names(immuneKO_mean)[order(immuneKO_mean,decreasing = F)[1:cutoff_value]]

cutoff_LFC <- 1.5
Top_LFC_tumorKO <- names(tumorKO_mean)[tumorKO_mean > cutoff_LFC]
bottom_LFC_tumorKO <- names(tumorKO_mean)[tumorKO_mean <- -cutoff_LFC]
Top_LFC_immuneKO <- names(immuneKO_mean)[immuneKO_mean > cutoff_LFC]
bottom_LFC_immuneKO <- names(immuneKO_mean)[immuneKO_mean <- -cutoff_LFC]

#draw heatmap:
input_LFC_tumorKO_top <- input_LFC_tumorKO[c(Top10_gene_tumorKO,Bottom10_gene_tumorKO),]
annotation_col <- data.frame(Genelist = c(rep("Top10",10),rep("Bottom10",10)))
rownames(annotation_col) <- c(Top10_gene_tumorKO,Bottom10_gene_tumorKO)

colnames(input_LFC_tumorKO_top) <- gsub("_LFC","",colnames(input_LFC_tumorKO_top))
colnames(input_LFC_tumorKO_top) <- gsub("[.]","_",colnames(input_LFC_tumorKO_top))
input_meta$file_name_all <- gsub("-","_",input_meta$file_name_all)
input_meta$file_name_all <- gsub("[.]","_",input_meta$file_name_all)
#colnames(input_LFC_tumorKO_top) <- gsub("-","_",colnames(input_LFC_tumorKO_top))
annotation_row <- data.frame(Treatment_type = input_meta$Treatment_type[input_meta$file_name_all %in% colnames(input_LFC_tumorKO_top)])
rownames(annotation_row) <- colnames(input_LFC_tumorKO_top)
HeatmapView(t(input_LFC_tumorKO_top), limit = c(-8,8), cluster_row = T, cluster_col = TRUE,
            filename = "/Users/yujijun/Desktop/QC_ImmuneScreens.png", width = 12, height = 8,
            annotation_col = annotation_col,annotation_row = annotation_row)



####not split and search for top10 and bottom10 ####
input_LFC_genesignature <- input_LFC_genesignature[rowMeans(is.na(input_LFC_genesignature)) <= 0.5,]
#make names consistent
input_meta <- input_meta[order(input_meta$KO_type,input_meta$Treatment_type),]
input_LFC_genesignature <- input_LFC_genesignature[,input_meta$file_name_all]
LFC_means <- rowMeans(input_LFC_genesignature,na.rm = T)
cutoff_value = 10
Top10 <- names(LFC_means)[order(LFC_means,decreasing = T)[1:cutoff_value]]
Bottom10 <- names(LFC_means)[order(LFC_means,decreasing = F)[1:cutoff_value]]


##heatmap for top10 and bottom10####
LFC_choose <- input_LFC_genesignature[c(Top10,Bottom10),]
annotation_col <- data.frame(Genelist = c(rep("Top10",10),rep("Bottom10",10)))
rownames(annotation_col) <- c(Top10,Bottom10)

annotation_row <- data.frame(Treatment_Type=input_meta$Treatment_type,KO_type = input_meta$KO_type)
rownames(annotation_row) <- colnames(LFC_choose)
LFC_choose <- LFC_choose[!rowSums(is.na(LFC_choose)) == dim(LFC_choose)[1],]
HeatmapView(t(LFC_choose), limit = c(-8,8), cluster_row = F, cluster_col = TRUE,
            filename = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/CRISPR_output_plot/heatmap/heatmap_for_top10_and_bottom.png", width = 12, height = 8,
            annotation_col = annotation_col,annotation_row = annotation_row)


####draw boxplot/violin, just for one gene, comparisons for Tumors and Immunecell_KO:####
library(tibble)
require(ggpubr)
rownames(LFC_choose)
LFC <- as.numeric(LFC_choose[1,])
df <- data.frame(LFC = LFC, KO_type =input_meta$KO_type)
my_comparisons <- list (c("Tumor_KO","Immunecell_KO"))
p <- ggviolin(df, x = "KO_type", y = "LFC", fill = "KO_type",
         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"),
         title = "STAT1")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") # Add significance levels

output_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/CRISPR_output_plot/Violinplot"
figure_name = paste("Violinplot_for_STAT1",".png",sep = "")
png(paste(output_path,figure_name,sep = "/"),height = 12,width = 17,units ="cm",res=150)
print(p)
dev.off()


####violin plot for all choose gene signature and all cohort ####
empty <- data.frame()  #concat all dataset together
datasets <- colnames(input_LFC_genesignature)
for(i in seq(1, dim(input_LFC_genesignature)[2])){
  Value = input_LFC_genesignature[,i]
  Dataset = rep(datasets[i],dim(input_LFC_genesignature)[1])
  dataframe_i <- data.frame(Value,Dataset)
  empty <- rbind(empty, dataframe_i)
  print(dim(empty))
}
empty <- empty[!is.na(empty$Value),]
empty$Dataset <- gsub("_LFC","",empty$Dataset)
p <- ggviolin(empty, x = "Dataset", y = "Value", fill = "Dataset",
              #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              add = "boxplot", add.params = list(fill = "white"),
              xlab = "Datasets",
              ylab = "LFC_value",
              title = "Gene signature in all datasets") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5,color = "#2F4F4F")) +
  theme(axis.text.x.bottom = element_text(size = 7,face = "bold")) +
  theme(axis.text.x.bottom = element_text(angle = 90,hjust=1)) +
  #scale_x_discrete(position = "") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 3,byrow = T)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 5,face = "bold")) +
  theme(axis.title.x = element_text(color = "#2F4F4F",face = "bold",size = 15)) +
  theme(axis.title.y = element_text(color = "#2F4F4F",face = "bold",size = 15)) +
  theme(axis.text.x = element_text(face = "bold",size = 20)) +
  theme(legend.position = "none")

show(p)
output_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/CRISPR_output_plot/Violinplot"
figure_name = paste("Violinplot_for_allgenes_allcohort",".png",sep = "")
png(paste(output_path,figure_name,sep = "/"),height = 15,width = 30,units ="cm",res=150)
print(p)
dev.off()



####draw volcano plot for top10 and bottom10 genes####
volcanodisplaymultigene <- function(file_path, top10_gene, bottom10_gene, figure_name, output_path){
  source('./R/Volcano.R')
  require(MAGeCKFlute)
  #input data
  dd.rra = ReadRRA(file_path)
  head(dd.rra)
  dd.rra$Official = toupper(dd.rra$Official)
  if(sum(top10_gene %in% dd.rra$Official) > 0){
    dd.rra$color <- "background"
    dd.rra$log_10 <- -log10(dd.rra$FDR)
    dd.rra$color[dd.rra$Official %in% top10_gene] <- "top10_gene"
    dd.rra$color[dd.rra$Official %in% bottom10_gene] <- "bottom10_gene"
    figure_title = figure_name
    subset = dd.rra[dd.rra$Official %in% c(top10_gene,bottom10_gene),]
    #draw the plot
    p<-Volcano(data=dd.rra, x="LFC",y="log_10",
               label_data = subset,fill="color",
               color="color",label = "Official",color_palette = c("#CFCFCF","#228B22","#FF0000"),title=figure_title)  #+ geom_text_repel(data = subset,aes(x = LFC, y = FDR),label = Official)
    #save the plot
    figure_output = paste(figure_title,".png",sep = "")
    png(paste(output_path,figure_output,sep = "/"),height = 12,width = 17,units ="cm",res=150)
    print(p)
    dev.off()
  }else{
    print(paste("There isn't gene in ",figure_name))
  }
}

basic_path = "/Users/yujijun/Dropbox (Partners HealthCare)/ImmuneScreens_verified"
all_file = list.files(basic_path,pattern = "*.gene_summary.txt",recursive = T)
all_file_list = paste(basic_path, all_file,sep = "/")
author_name = str_split_fixed(all_file,"_",n=2)[,1]
control_treatment = str_split_fixed(basename(all_file),"[.]",n=2)[,1]
figure_name = paste(author_name, control_treatment,sep = "_")
output_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/CRISPR_output_plot/Volcano"
for(i in seq(1,length(all_file_list))){
  volcanodisplaymultigene(all_file_list[i],Top10,Bottom10,figure_name[i],output_path)
}




####violin plot for postive and negtive genes####
gene_list$logical <- "postive weight"
gene_list$logical[gene_list$V2 <= 0] <- "negative weight"
rownames(input_LFC_genesignature) <- rownames(input_LFC)[rownames(input_LFC) %in% as.vector(gene_list$V1)]
input_LFC_genesignature <- add_column(input_LFC_genesignature,gene_name=rownames(input_LFC_genesignature),.before="Burr_PD_L1_CMTM6_INFr_LFC")
gene_signature <- merge(gene_list,input_LFC_genesignature,by.x = "V1",by.y = "gene_name")

empty <- data.frame()  #concat all dataset together
for(i in seq(4, dim(gene_signature)[2])){
  Value = gene_signature[,i]
  sample_name = rep(colnames(gene_signature)[i],dim(gene_signature)[1])
  logical_data = gene_signature$logical
  dataframe_i <- data.frame(Value,sample_name,logical_data)
  empty <- rbind(empty, dataframe_i)
  print(dim(empty))
}
#df <- data.frame(LFC = LFC, KO_type =input_meta$KO_type)
my_comparisons <- list (c("postive weight","negative weight"))
empty <- empty[!is.na(empty$Value),]
p <- ggviolin(empty, x = "sample_name", y = "Value", fill = "logical_data",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              add = "boxplot", add.params = list(fill = "logical_data"),
              title = "Positive and negative weight in all datasets",
              xlab = "Datasets", ylab = "LFC_Value")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +# Add significance levels
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5,color = "#2F4F4F")) +
  theme(axis.text.x.bottom = element_text(size = 7,face = "bold")) +
  theme(axis.text.x.bottom = element_text(angle = 90,hjust=1)) +
  #scale_x_discrete(position = "") +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 3,byrow = T)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 10,face = "bold")) +
  theme(axis.title.x = element_text(color = "#2F4F4F",face = "bold",size = 15)) +
  theme(axis.title.y = element_text(color = "#2F4F4F",face = "bold",size = 15))
output_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/CRISPR_output_plot/Violinplot"
figure_name = paste("Violinplot of Positive and negative weight in all cohort",".png",sep = "")
png(paste(output_path,figure_name,sep = "/"),height = 15,width = 30,units ="cm",res=150)
print(p)
dev.off()



####stouffer z-score####
intersection_gene <- intersect(gene_list$V1,rownames(input_LFC))
gene_list_interset <- gene_list %>% filter(V1 %in% intersection_gene)
input_LFC_interset <- input_LFC %>% filter(rownames(input_LFC) %in% intersection_gene)
not_all_na <- function(x) any(!is.na(x))
input_LFC_interset <- input_LFC_interset %>% select_if(not_all_na)
w <- gene_list_interset$V2
cohort <- colnames(input_LFC_interset)
Stouff_zscore <- c()

for (i in seq(1,ncol(input_LFC_interset))){
  print(paste(i,cohort[i],sep = " "))
  z <- input_LFC_interset[,i]
  Stouff_zscore_i <- Stoufferzscore(w,z)
  Stouff_zscore <- c(Stouff_zscore,Stouff_zscore_i)
}

Stouff_df <- data.frame(cohort, Stouff_zscore)
Stouff_df$cohort <- gsub("_LFC","",Stouff_df$cohort)
singlegenebarplot(Stouff_df,x="cohort",y ="Stouff_zscore",fill = "#2F4F4F",sort.by.groups = T,ylab = "Weighted sum score",xlab = "Dataset",title = "barplot of Weighted sum score",legend.title = NA )



