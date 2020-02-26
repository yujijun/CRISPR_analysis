# --------------
# Date:  2019-11-18 09:57:10
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: Add all gene signature together
#
require(dplyr)
require(stringi)
require(stringr)
require(MAGeCKFlute)
library(psycho)

#Input all file we would like to preprocess
basic_path = "~/Dropbox (Partners HealthCare)/ImmuneScreens_verified"
all_file = list.files(basic_path,pattern = "*.gene_summary.txt",recursive = T)
all_file_list = paste(basic_path, all_file,sep = "/")
author_name = str_split_fixed(all_file,"_",n=2)[,1]
file_name = basename(all_file_list)
file_name = gsub(".gene_summary.txt","",file_name)
file_name_all = paste(author_name,file_name,sep = "_")
file_name_all = gsub("Pan_Pan","Pan",file_name_all)
file_name_all = gsub("Shifrut_Shifrut","Shifrut",file_name_all)

# concatenate all data together
gg = data.frame(Human = "", stringsAsFactors = FALSE)
for(i in seq(1,length(all_file_list))){
  name_list = c("RRAscore","FDR","LFC")
  name_list = paste(file_name_all[i],name_list,sep = "_")
  tmp= read.table(all_file_list[i],header = T, sep = "\t")
  tmp_select = tmp %>% select(id,neg.score,neg.fdr,neg.lfc)
  tmp_select$id = toupper(tmp_select$id)
  if(any(grepl("RIK$", tmp_select$id))){
    tmp_select$Human = TransGeneID(tmp_select$id, "Symbol", "Symbol",
                            fromOrg = "mmu", toOrg = "hsa")
    #if(grepl("High|Low", all_file_list[i])) tmp_select$neg.fdr[tmp_select$neg.fdr<0] = 0
    tmp_select = tmp_select[,-1]
    colnames(tmp_select) = c(name_list,"Human")
  }else{
    tmp_select$Human = tmp_select$id
    tmp_select = tmp_select[,-1]
    #if(grepl("High|Low", all_file_list[i])) tmp_select$neg.fdr[tmp_select$neg.fdr<0] = 0
    colnames(tmp_select) = c(name_list,"Human")
  }
  tmp_select = tmp_select[order(-abs(tmp_select[,2])), ]
  tmp_select = tmp_select[!duplicated(tmp_select$Human), ]
  gg = merge.data.frame(gg, tmp_select, by = "Human", all = TRUE)
}
gg = gg[-1,]
gg = gg[!(is.na(gg$Human)|gg$Human==""), ]
rownames(gg) = gg$Human
gg = gg[,-1]
dim(gg)

#select
gg_RRAscore <- gg %>% select(contains("RRAscore"))
gg_LFC <- gg %>% select(contains("LFC"))

delete_colindex = c(2,3,7,9,11,20,22,23,25,27,29,30,32)
gg_RRAscore <- gg_RRAscore[,-delete_colindex]
gg_LFC  <- gg_LFC[,-delete_colindex]

#gg_LFC rnorm for unsorted situation
gg_LFC$`Burr_PD-L1_CMTM6_INFr_LFC`[is.na(gg_LFC$`Burr_PD-L1_CMTM6_INFr_LFC`)] <- 0
gg_LFC$`Burr_PD-L1_CMTM6_INFr_LFC`[gg_LFC$`Burr_PD-L1_CMTM6_INFr_LFC` < 0] <- rnorm(sum(gg_LFC$`Burr_PD-L1_CMTM6_INFr_LFC` < 0),0,0.1)
gg_LFC$Freeman_B16F10_MHCI_Low_LFC[is.na(gg_LFC$Freeman_B16F10_MHCI_Low_LFC)] <- 0
gg_LFC$Freeman_B16F10_MHCI_Low_LFC[gg_LFC$Freeman_B16F10_MHCI_Low_LFC < 0] <- rnorm(sum(gg_LFC$Freeman_B16F10_MHCI_Low_LFC < 0),0,0.1)
meltLFC <- melt(gg_LFC)
p <- ggplot(meltLFC, aes(x=variable, y=value)) +
  geom_violin(fill="variable") +
  theme(axis.text.x = element_text(angle = 90))
show(p)
#### max normalization of dataframe####
maxs <- apply(abs(gg_LFC),2,max,na.rm = TRUE)
scaledgg_LFC <- t(t(gg_LFC)/maxs*median(maxs))
meltLFC <- melt(scaledgg_LFC)
p <- ggplot(meltLFC, aes(x=X2, y=value)) +
  geom_violin(fill="variable") +
  theme(axis.text.x = element_text(angle = 90))
show(p)

#gg_matrix <- as.matrix(gg_RRAscore)
library(ggplot2)
p <- ggplot(meltData, aes(factor(variable), value))
p + geom_boxplot() + facet_wrap(~variable, scale="free")

data <- data.frame(v1=rnorm(100),v2=rnorm(100),v3=rnorm(100), v4=rnorm(100))
library(reshape)
meltData <- melt(data)
boxplot(data=meltData, value~variable)

meltRRA <- melt(gg_RRAscore)
meltRRA <- meltRRA[!is.na(meltRRA$value),]
meltRRA <- meltRRA %>% filter(value > -10 & value < 10)
boxplot(data=meltRRA, value~variable,)
meltLFC <- melt(gg_LFC)
p <- ggplot(meltLFC, aes(x=variable, y=value)) +
  geom_violin(fill="variable") +
  theme(axis.text.x = element_text(angle = 90))

show(p)

max_fun <- function(x){
  y <- max(x,na.rm = T)
  return(y)
}
nonzero_prop <- function(x){
  y <- (length(x)-sum(is.na(x)))/length(x)
}
zero <- as.data.frame(apply(scaledgg_LFC,2,nonzero_prop))
gg_matrix <- as.matrix(gg_LFC)
zero$name <- rownames(zero)
p <- ggplot(zero,aes(y=apply(scaledgg_LFC,2,nonzero_prop),x=name)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90))
show(p)
####zscore for matrix####
# gg_RRAscore_z <- gg_RRAscore %>%
#   psycho::standardize(center=FALSE,scale= TRUE)
# rownames(gg_RRAscore_z) <- rownames(gg)
# gg_LFC_z <- gg_LFC %>%
#   psycho::standardize()
# rownames(gg_LFC_z) <- rownames(gg)
# colnames(gg_RRAscore_z) <- gsub("[.]","_",colnames(gg_RRAscore_z))
# colnames(gg_RRAscore_z) <- gsub("-","_",colnames(gg_RRAscore_z))
# colnames(gg_LFC_z) <- gsub("[.]","_",colnames(gg_LFC_z))
# colnames(gg_LFC_z) <- gsub("-","_",colnames(gg_LFC_z))

#description file

KO_type = c("Tumor_KO","Tumor_KO","Tumor_KO","Immunecell_KO",
            "Immunecell_KO","Immunecell_KO","Tumor_KO","Tumor_KO","Tumor_KO","Tumor_KO",
            "Tumor_KO","Tumor_KO","Tumor_KO","Tumor_KO","Immunecell_KO",
            "Tumor_KO","Tumor_KO","Immunecell_KO","Immunecell_KO","Tumor_KO")
Treatment_type = c("IFNr_treatment","IFNr_treatment","NK_coculture","phagocytosis",
                  "Tcell","Bcell","Tcell_coculture","Tcell_coculture","Tcell_coculture","Dox_treatment",
                  "Tcell_coculture","Tcell_coculture","Tcell_coculture","Tcell_coculture","BMDCs",
                  "IFNr_treatment","NK_coculture","macrophage","Tcell","Tcell_coculture")
QC_meta = data.frame(columns_names,KO_type,Treatment_type)
colnames(meta)[1] <- "Sample"
meta$Final_type <- "Immune cell KO"
meta$Final_type[grep("_coculture",meta$Treatment_type)] <- rep("TumorKO_sorting",10)
sorting_index <- setdiff(seq(1,nrow(meta))[meta$KO_type == "Tumor_KO"],grep("_coculture",meta$Treatment_type))
meta$Final_type[grep("_coculture",meta$Treatment_type)] <- rep("TumorKO_coculture",10)
meta$Final_type[sorting_index] <- rep("TumorKO_sorting",4)
meta$Final_type <- gsub("Immune cell KO","Immuecell_KO",meta$Final_type)
write.table(meta,file = "./QC_meta.txt",sep = "\t",col.names = T,row.names = F, quote = F)

#QC_meta$file_name_all <- gsub("-","_",QC_meta$file_name_all)
#QC_meta$file_name_all <- gsub("[.]","_",QC_meta$file_name_all)


#output all file after Z-score
output = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/data"
output = "/Users/yujijun/Desktop"
write.table(gg_RRAscore,file=paste(output,"RRAscore.txt",sep = "/"),sep = "\t",row.names = T,col.names = T,quote = F)
write.table(gg_LFC,file=paste(output,"LFC.txt",sep = "/"),sep = "\t",row.names = T,col.names = T,quote = F)
write.table(QC_meta,file=paste(output,"QC_meta.txt",sep = "/"),sep = "\t",row.names = F,col.names = T,quote = F)

#output all file before z-score
rownames(gg_LFC) <- rownames(gg)
colnames(gg_RRAscore) <- gsub("[.]","_",colnames(gg_RRAscore))
colnames(gg_RRAscore) <- gsub("-","_",colnames(gg_RRAscore))
colnames(gg_LFC) <- gsub("[.]","_",colnames(gg_LFC))
colnames(gg_LFC) <- gsub("-","_",colnames(gg_LFC))
output = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/data"
write.table(gg_RRAscore,file=paste(output,"RRAscore_original.txt",sep = "/"),sep = "\t",row.names = T,col.names = T,quote = F)
write.table(gg_LFC,file=paste(output,"LFC_original.txt",sep = "/"),sep = "\t",row.names = T,col.names = T,quote = F)

# #### max normalization of dataframe####
maxs <- apply(abs(gg_LFC),2,max,na.rm = TRUE)
scaledgg_LFC <- t(t(gg_LFC)/maxs*median(maxs))
# rownames(scaledgg_LFC) <- rownames(gg)
# colnames(scaledgg_LFC) <- gsub("[.]","_",colnames(scaledgg_LFC))
# colnames(scaledgg_LFC) <- gsub("-","_",colnames(scaledgg_LFC))

# maxs <- apply(abs(gg_RRAscore),2,max,na.rm=TRUE)
# scaledgg_RRAscore <- t(t(gg_RRAscore)/maxs*median(maxs))
# rownames(scaledgg_RRAscore) <- rownames(gg)
# colnames(scaledgg_RRAscore) <- gsub("[.]","_",colnames(scaledgg_RRAscore))
# colnames(scaledgg_RRAscore) <- gsub("-","_",colnames(scaledgg_RRAscore))

columns_names <- colnames(scaledgg_LFC)
columns_names <- gsub("_LFC","",columns_names)
columns_names <- gsub("_r1","",columns_names)
columns_names <- gsub("_1r","",columns_names)
columns_names[19] <- "Shifrut_hCD8T_CFSE"
colnames(scaledgg_LFC) <- columns_names
output = "/Users/yujijun/Desktop"
#write.table(scaledgg_LFC,file=paste(output,"RRAscore_scale(nocenter).txt",sep = "/"),sep = "\t",row.names = T,col.names = T,quote = F)
write.table(scaledgg_LFC,file=paste(output,"LFC_maxscale.txt",sep = "/"),sep = "\t",row.names = T,col.names = T,quote = F)









