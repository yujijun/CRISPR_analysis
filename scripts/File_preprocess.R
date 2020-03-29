# --------------
# Date:  2020-03-03 22:44:05 
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project:CRISPR dataset preprocessing 
# 
library(dplyr)
library(tidyverse)

#data processing####
#all parameter to be used 
base.path <- "/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/"
rawdata.rowname12 <- c("sgRNA","Gene")
#all input and output files' name
Input.path <- "/Users/yujijun/Dropbox (Partners HealthCare)/ImmuneScreens_v2/"
Input.file1 <- "Jiang_2019_CellRep_31365872" #still need to change the two index file into one.
Input.file2 <- "Burr_2017_Nature_28813417"
Input.file3 <- "Burr_2019_CancerCell_31564637"
Input.file4 <- "Codina_2019_CellSyst_30797773"
Input.file5 <- "Dong_2019_Cell_31442407"
Input.file6 <- "Freeman_2019_CellRep_31509742"
Input.file7 <- "Haney_2018_NatGenet_30397336"
Input.file8 <- "Henriksson_2018_Cell_30639098"
Input.file9 <- "Kearney_2018_SciImmunol_29776993"
Input.file10 <- "Liu_2018_NatMed_30559422"
Input.file11 <- "Manguso_2017_Nature_28723893"
Input.file12 <- "Mezzadra_2017_Nature_28813410"
Input.file13 <- "Pan_2018_Science_29301958"
Input.file14 <- "Parnas_2015_Cell_26189680"
Input.file15 <- "Patel_2017_Nature_28783722"
Input.file16 <- "Pech_2019_Elife_31452512"
Input.file17 <- "Ritchie_2019_MolecularCell_31126740"
Input.file18 <- "Shifrut_2018_Cell_30449619"
Input.file19 <- "Vredevoogd_2019_Cell_31303383"
#file output files
Output.path <- paste0(base.path, "Output/")
Output.file1.1 <- "31365872_Benjamin E.Gewurz_CellRep_2019_1"
Output.file1.2 <- "31365872_Benjamin_EGewurz_CellRep_2019_2"
Output.file1.1 <- "31365872_Benjamin E.Gewurz_CellRep_2019_1"
Output.file1.2 <- "31365872_Benjamin_EGewurz_CellRep_2019_2"
Output.file2 <- "28813417_Mark A. Dawson_Nature_2017"
Output.file3 <- "31564637_Mark A. Dawson_CancerCell_2019"
Output.file4 <- "30797773_Sidi Chen_CellSyst_2019"
Output.file5.1 <- "31442407_Sidi Chen_Cell_2019_1"
Output.file5.2 <- "31442407_Sidi Chen_Cell_2019_2"
Output.file6.1 <- "31509742_Jane Oliaro_CellRep_2019_1"
Output.file6.2 <- "31509742_Jane Oliaro_CellRep_2019_2"
Output.file7 <- "30397336_Michael C. Bassik_NatGenet_2018"
Output.file8 <- "30639098_Sarah A. Teichmann_Cell_2018"
Output.file9.1 <- "29776993_Jane Oliaro_SciImmunol_2018_1"
Output.file9.2 <- "29776993_Jane Oliaro_SciImmunol_2018_2"
Output.file9.3 <- "29776993_Jane Oliaro_SciImmunol_2018_3"
Output.file10 <- "30559422_E. Robert McDonald III_NatMed_2018"
Output.file11.1 <- "28723893_W. Nicholas Haining_Nature_2017_1"
Output.file11.2 <- "28723893_W. Nicholas Haining_Nature_2017_2"
Output.file11.3 <- "28723893_W. Nicholas Haining_Nature_2017_3"
Output.file11.4 <- "28723893_W. Nicholas Haining_Nature_2017_4"
Output.file12 <- "28813410_Ton N.M. Schumacher_Nature_2017" #lib issue
Output.file13 <- "29301958_Kai W. Wucherpfennig_Science_2018"
Output.file14 <- "26189680_Aviv Regev_Cell_2015"
Output.file15 <- "28783722_Nicholas P. Restifo_Nature_2017"
Output.file16.1 <- "31452512_Jeffrey Settleman_Elife_2019_1"
Output.file16.2 <- "31452512_Jeffrey Settleman_Elife_2019_2"
Output.file17 <- "31126740_Lingyin Li_MolecularCell_2019"
Output.file18 <- "30449619_Alexander Marson_Cell_2018"
Output.file19 <- "31303383_Daniel S. Peeper_Cell_2019"
#### Data processing #### 
setwd(paste0(Input.path,Input.file10)) #para change first
#>>>>>>>>>>rawcount data####
#file two####
rawdata <- read.table("./count/PD-L1_CMTM6.count.txt",header = T,stringsAsFactors = F)
#file three####
#input dataset
rawdata <- read.table("./count/K562_HLA-ABC.count.txt",header = T)
rawdata.1 <- read.table("./count/K562_HLA-B.count.txt",header = T)
lib.1 <- readxl::read_xlsx("./lib/41467_2017_BFncomms15178_MOESM279_ESM.xlsx",sheet = "Guide sequences",col_names = F)
lib.2 <- readxl::read_xlsx("./lib/41467_2017_BFncomms15178_MOESM278_ESM.xlsx",sheet = "Guide sequences",col_names = F)
colnames(lib.1) <- c("sgRNA_ID","sgRNA_seq")
colnames(lib.2) <- c("sgRNA_ID","sgRNA_seq")

#file four:####
#input next:
rawdata <- read.table("./rawdata/count_final_v1.txt",header = T,stringsAsFactors = F)
rawdata <- rawdata[c(1001:nrow(rawdata)),]
lib <- read.delim("./lib/mbrie_lib.txt",header = T,stringsAsFactors = F)
rownames(rawdata) <- seq(1,nrow(rawdata),1)
lib <- rownames_to_column(lib,var = "sgRNA_ID")
rawdata <- rownames_to_column(rawdata,var = "sgRNA_ID")
merge.lib.rawdata <- merge(lib,rawdata,by.x = "sgRNA_ID",by.y = "sgRNA_ID")
merge.lib.rawdata <- merge.lib.rawdata[order(as.numeric(merge.lib.rawdata$sgRNA_ID)),]

lib <- merge.lib.rawdata[,c("sgRNA.Target.Sequence","Target.Gene.Symbol")]
rownames(lib) <- paste0("sgRNA",seq(1,nrow(lib),1))
lib <- rownames_to_column(lib,var = "sgRNA")
lib.colnames <- c("sgRNA_ID","sgRNA_seq","Gene")
colnames(lib) <- lib.colnames

rawdata <- merge.lib.rawdata[13:ncol(merge.lib.rawdata)]
rawdata[,1] <- paste0("sgRNA",seq(1,nrow(lib),1))

rownames(rawdata)[c(1,2)] <- rawdata.rowname12

#file five:####
#invitro analysis
rawdata <- read.delim(file = "./rawdata/Invitro_count.txt",sep = "\t",stringsAsFactors = F,header = T)
lib <- read.delim(file="./lib/mouse_geckov2_library_combine.csv",sep = ",",stringsAsFactors = F,header = F)
lib <- lib %>% rownames_to_column(var = "MKO_ID")
lib$MKO_ID <- seq(10000001,10000001+nrow(lib)-1,1)
lib$MKO_ID <- paste("MKO",lib$MKO_ID,sep = "")
merge.lib.rawdata <- merge(rawdata,lib, by.x = "sgRNA",by.y = "MKO_ID")
rownames(merge.lib.rawdata) <- paste0("sgRNA",seq(1, nrow(merge.lib.rawdata),1))
merge.lib.rawdata <- merge.lib.rawdata %>% rownames_to_column(var = "sgRNA_ID")

lib <- merge.lib.rawdata[c("sgRNA_ID","V2","gene")]
lib.colnames <- c("sgRNA_ID","sgRNA_seq","Gene")
colnames(lib) <- lib.colnames

rawdata <- merge.lib.rawdata[c(1:12)]
rawdata <- rawdata[-2]
rawdata.rowname12 <- c("sgRNA","Gene")
rownames(rawdata)[c(1,2)] <- rawdata.rowname12
#invivo analysis:

#file five.2：
rawdata <- read.delim(file = "./rawdata/Invivo_normcount.txt",sep = "\t",stringsAsFactors = F,header = T)
lib <- read.delim(file="./lib/mouse_geckov2_library_combine.csv",sep = ",",stringsAsFactors = F,header = F)
lib <- lib %>% rownames_to_column(var = "MKO_ID")
lib$MKO_ID <- seq(10000001,10000001+nrow(lib)-1,1)
lib$MKO_ID <- paste("MKO",lib$MKO_ID,sep = "")
merge.lib.rawdata <- merge(rawdata,lib, by.x = "sgRNA",by.y = "MKO_ID")

lib <- merge.lib.rawdata[c("sgRNA","V2","gene")]
lib.colnames <- c("sgRNA_ID","sgRNA_seq","Gene")
colnames(lib) <- lib.colnames

rawdata <- merge.lib.rawdata[c(1:15)]
rawdata.rowname12 <- c("sgRNA","Gene")
colnames(rawdata)[c(1,2)] <- rawdata.rowname12


#file six:####
rawdata.1 <- read.delim(file = "./rawdata/B16F10_NK.count.txt",sep = "\t",stringsAsFactors = F,header = T)
rawdata.2 <- read.delim(file = "./rawdata/MHCI_low.count.txt",sep = "\t",stringsAsFactors = F,header = T)
lib <- read.delim(file = "./lib/library_mBrie.txt",sep = "\t",stringsAsFactors = F,header = T)
lib <- lib[!duplicated(lib$sgRNA_ID),]
merge.lib.rawdata.1 <- merge(rawdata.1, lib,by.x = "sgRNA",by.y = "sgRNA_ID")
merge.lib.rawdata.2 <- merge(rawdata.2, lib,by.x = "sgRNA",by.y = "sgRNA_ID")
rawdata.1 <- merge.lib.rawdata.1[c(1:5)]
rawdata.1$sgRNA <- paste0("sgRNA",seq(1,nrow(rawdata.1),1))
rawdata.2 <- merge.lib.rawdata.2[c(1:4)]
rawdata.2$sgRNA <- paste0("sgRNA",seq(1,nrow(rawdata.2),1))
lib.1 <- merge.lib.rawdata.1[c("sgRNA","Gene")]
lib.2 <- merge.lib.rawdata.2[c("sgRNA","Gene")]
rownames(lib.1) <- paste0("sgRNA",seq(1,nrow(lib.1),1))
lib.1 <- lib.1 %>% rownames_to_column(var = "sgRNA_ID")
rownames(lib.2) <- paste0("sgRNA",seq(1,nrow(lib.2),1))
lib.2 <- lib.2 %>% rownames_to_column(var = "sgRNA_ID")
colnames(lib.1) <- c("sgRNA_ID","sgRNA_Seq","Gene")
colnames(lib.2) <- c("sgRNA_ID","sgRNA_Seq","Gene")

#file seven:####
library(dplyr)
rawdata <- read.delim(file = "./rawdata/WG_screen_count.txt",sep = "\t",stringsAsFactors = F,header = T)
rawdata.midbead <- rawdata[c(1,2,grep(pattern = "midbead",x = colnames(rawdata)))]
rawdata.smallbead <- rawdata[c(1,2,grep(pattern = "smallbead",x=colnames(rawdata)))]
rawdata.Zymosin <- rawdata[c(1,2,grep(pattern = "Zymosin",x=colnames(rawdata)))]
rawdata.PosBead <- rawdata[c(1,2,grep(pattern = "PosBead",x=colnames(rawdata)))]
rawdata.myelin <- rawdata[c(1,2,grep(pattern = "myelin",x=colnames(rawdata)))]
rawdata.iggRBC <- rawdata[c(1,2,grep(pattern = "iggRBC",x=colnames(rawdata)))]
rawdata.compRBC <- rawdata[c(1,2,grep(pattern = "compRBC",x=colnames(rawdata)))]
rawdata.bigbead <- rawdata[c(1,2,grep(pattern = "bigbead",x=colnames(rawdata)))]

#lib <- read.delim(file = "./rawdata/gRNA_library.txt",sep = "\t",stringsAsFactors = F,header = T)
#merge.lib.rawdata <- merge(rawdata,lib, by.x = "gR0",by.y = "id")
lib <- rawdata[c(1,2)]
rownames(lib) <- paste0("sgRNA",seq(1,nrow(lib),1))
lib <- lib %>% rownames_to_column(var = "sgRNA_ID")
lib.colnames <- c("sgRNA_ID","sgRNA_Seq","Gene")
colnames(lib) <- lib.colnames
rawdata.rowname12 <- c("sgRNA","Gene")

colnames(rawdata.midbead)[c(1,2)] <- rawdata.rowname12
colnames(rawdata.smallbead)[c(1,2)] <- rawdata.rowname12
colnames(rawdata.Zymosin)[c(1,2)] <- rawdata.rowname12
colnames(rawdata.PosBead)[c(1,2)] <- rawdata.rowname12
colnames(rawdata.myelin)[c(1,2)] <- rawdata.rowname12
colnames(rawdata.iggRBC)[c(1,2)] <- rawdata.rowname12
colnames(rawdata.compRBC)[c(1,2)] <- rawdata.rowname12
colnames(rawdata.bigbead)[c(1,2)] <- rawdata.rowname12

#file eight####
rawdata <- read.delim(file = "./rawdata/grnacnt.txt",sep = "\t",stringsAsFactors = F,header = T)
lib <- rawdata[c(1,2)]
rownames(lib) <- paste0("sgRNA",seq(1,nrow(lib),1))
lib <- lib %>% rownames_to_column(var="sgRNA_ID")
lib.colnames <- c("sgRNA_ID","sgRNA_seq","Gene")
colnames(lib) <- lib.colnames
rawdata.rowname12 <- c("sgRNA","Gene")
colnames(rawdata)[c(1,2)] <- rawdata.rowname12
rawdata$sgRNA <- paste0("sgRNA",seq(1,nrow(lib),1))

#file nine####
rawdata <- read.delim(file = "./rawdata/B16_screen.count.txt",sep = "\t",stringsAsFactors = F,header = T)
rawdata.1 <- read.delim(file = "./rawdata/MC38_screen.count.txt",sep = "\t",stringsAsFactors = F,header = T)
rawdata.2 <- read.delim(file = "./rawdata/MC38_screen2.count.txt",sep = "\t",stringsAsFactors = F,header = T)
rawdata.rowname12 <- c("sgRNA","Gene")
colnames(rawdata)[c(1,2)] <- rawdata.rowname12
colnames(rawdata.1)[c(1,2)] <- rawdata.rowname12
colnames(rawdata.2)[c(1,2)] <- rawdata.rowname12

merge.lib.rawdata <- merge(rawdata,rawdata.1,by.x = "sgRNA",by.y = "sgRNA")
merge.lib.rawdata <- merge(merge.lib.rawdata,rawdata.2,by.x = "sgRNA",by.y = "sgRNA")

lib <- rawdata[c(1,2)]
rownames(lib) <- paste0("sgRNA",seq(1,nrow(lib),1))
lib <- lib %>% rownames_to_column(var = "sgRNA_ID")

lib.1 <- rawdata.1[c(1,2)]
rownames(lib.1) <- paste0("sgRNA",seq(1,nrow(lib.1),1))
lib.1 <- lib.1 %>% rownames_to_column(var = "sgRNA_ID")

lib.2 <- rawdata.2[c(1,2)]
rownames(lib.2) <- paste0("sgRNA",seq(1,nrow(lib.2),1))
lib.2 <- lib.2 %>% rownames_to_column(var = "sgRNA_ID")

lib.colnames <- c("sgRNA_ID","sgRNA_seq","Gene")
colnames(lib) <- lib.colnames
colnames(lib.1) <- lib.colnames
colnames(lib.2) <- lib.colnames

#file ten####
library(dplyr)
library(tibble)
rawdata <- read.delim(file = "./rawdata/hsc4_ADAR_suppressor_count.txt",header = T,stringsAsFactors = F,sep = "\t")
lib <- rawdata[c(1,2)]
rownames(lib) <- paste0("sgRNA",seq(1,nrow(lib),1))
lib <- lib %>% rownames_to_column(var = "sgRNA_ID")
lib.colnames <- c("sgRNA_ID","sgRNA_seq","Gene")
colnames(lib) <- lib.colnames
rawdata.rowname12 <- c("sgRNA","Gene")
colnames(rawdata)[c(1,2)] <- rawdata.rowname12
#file eleven####
rawdata.1 <- read.delim(file = "./rawdata/Manguso_2017.pool1.count_normalized.txt",header = T,stringsAsFactors = F,sep = "\t")
rawdata.2 <- read.delim(file = "./rawdata/Manguso_2017.pool2.count_normalized.txt",header = T,stringsAsFactors = F,sep = "\t")
rawdata.3 <- read.delim(file = "./rawdata/Manguso_2017.pool3.count_normalized.txt",header = T,stringsAsFactors = F,sep = "\t")
rawdata.4 <- read.delim(file = "./rawdata/Manguso_2017.pool4.count_normalized.txt",header = T,stringsAsFactors = F,sep = "\t")
colnames(rawdata.1) <- paste0(colnames(rawdata.1),".pool1")
colnames(rawdata.2) <- paste0(colnames(rawdata.2),".pool2")
colnames(rawdata.3) <- paste0(colnames(rawdata.3),".pool3")
colnames(rawdata.4) <- paste0(colnames(rawdata.4),".pool4")
merge.all <- merge(merge(rawdata.1,rawdata.2,by.x = "Gene.pool1",by.y = "Gene.pool2"),
                   merge(rawdata.3,rawdata.4,by.x = "Gene.pool3",by.y = "Gene.pool4"),
                   by.x = "Gene.pool1",by.y = "Gene.pool3")
Annotation <- "/Users/yujijun/Dropbox (Partners HealthCare)/ImmuneScreens_paper/ImmuneScreens_Annotation.csv"
Annotation <- read.delim(file = Annotation,stringsAsFactors = F,sep = ",")
pool1_control <- paste0(stringi::stri_split_fixed(Annotation[43,"Control"],",")[[1]],".pool1")
pool2_control <- paste0(stringi::stri_split_fixed(Annotation[44,"Control"],",")[[1]],".pool2")
pool3_control <- paste0(stringi::stri_split_fixed(Annotation[45,"Control"],",")[[1]],".pool3")
pool4_control <- paste0(stringi::stri_split_fixed(Annotation[46,"Control"],",")[[1]],".pool4")
poolall_control <- c(pool1_control,pool2_control,pool3_control,pool4_control)

poolall_GVAX_aPD1_treatment <- c()
pool.vector <- c(".pool1",".pool2",".pool3",".pool4")
for (i in 43:46){
  pool_treatment <- paste0(stringi::stri_split_fixed(Annotation[i,"Treatment"],",")[[1]],pool.vector[(i-42)])
  poolall_GVAX_aPD1_treatment <- c(poolall_GVAX_aPD1_treatment,pool_treatment)
}
poolall_GVAX_treatment <- c()
pool.vector <- c(".pool1",".pool2",".pool3",".pool4")
for (i in 47:50){
  pool_treatment <- paste0(stringi::stri_split_fixed(Annotation[i,"Treatment"],",")[[1]],pool.vector[(i-46)])
  poolall_GVAX_treatment <- c(poolall_GVAX_treatment,pool_treatment)
}
poolall_TCRa_KO_treatment <- c()
pool.vector <- c(".pool1",".pool2",".pool3",".pool4")
for (i in 51:54){
  pool_treatment <- paste0(stringi::stri_split_fixed(Annotation[i,"Treatment"],",")[[1]],pool.vector[(i-50)])
  poolall_TCRa_KO_treatment <- c(poolall_TCRa_KO_treatment,pool_treatment)
}
poolall_In_Vitro_treatment <- c()
pool.vector <- c(".pool1",".pool2",".pool3",".pool4")
for (i in 55:58){
  pool_treatment <- paste0(stringi::stri_split_fixed(Annotation[i,"Treatment"],",")[[1]],pool.vector[(i-54)])
  poolall_In_Vitro_treatment <- c(poolall_In_Vitro_treatment,pool_treatment)
}
rawdata.GVAX_aPD1 <- merge.all[c(poolall_control,poolall_GVAX_aPD1_treatment)]
rawdata.GVAX <- merge.all[c(poolall_control,poolall_GVAX_treatment)]
rawdata.TCRa_KO <- merge.all[c(poolall_control,grep("TCRa_KO",colnames(merge.all),value = T))]
rawdata.In_Vitro <- merge.all[c(poolall_control,poolall_In_Vitro_treatment)]
lib <- merge.all[c(2,1)]
rownames(lib) <- paste0("sgRNA",seq(1,nrow(lib),1))
lib <- lib %>% rownames_to_column(var="sgRNA_ID")
colnames(lib) <- lib.colnames 

merge.all <- merge.all[-2]
rownames(merge.all) <- paste0("sgRNA",seq(1,nrow(lib),1))
merge.all <- merge.all %>% rownames_to_column(var="sgRNA")
rawdata.rowname12 <- c("sgRNA","Gene")
colnames(merge.all)[c(1,2)] <- rawdata.rowname12
rawdata.GVAX_aPD1 <- cbind(merge.all[c(1,2)],rawdata.GVAX_aPD1)
rawdata.GVAX <- cbind(merge.all[c(1,2)],rawdata.GVAX)
rawdata.TCRa_KO <- cbind(merge.all[c(1,2)],rawdata.TCRa_KO)
rawdata.In_Vitro <- cbind(merge.all[c(1,2)],rawdata.In_Vitro)
#file twelve####
rawdata <- read.delim(file = "./rawdata/all.count.txt",header = T,stringsAsFactors = F,sep="\t")
lib <- rawdata[c(1,2)]
rownames(lib) <- paste0("sgRNA",seq(1,nrow(lib),1))
lib <- lib %>% rownames_to_column(var = "sgRNA_ID")
colnames(lib) <- lib.colnames
colnames(rawdata) <- rawdata.rowname12

#file thirteen ####
rawdata <- read.delim(file ="./rawdata/Pan_2018_normcount.txt",header = T,sep = "\t",stringsAsFactors = F)
lib <- rawdata[c(1,2)]
rownames(lib) <- paste0("sgRNA",seq(1,nrow(lib),1))
lib <- lib %>% rownames_to_column(var = "sgRNA_ID")
colnames(lib) <- lib.colnames

#file fourteen####
rawdata.1 <- read.delim(file = "./rawdata/LPS_TNF_EX1.count.txt",header = T,stringsAsFactors = F,sep = "\t")
rawdata.2 <- read.delim(file = "./rawdata/LPS_TNF_EX2.count.txt",header = T,stringsAsFactors = F,sep = "\t")
rawdata.3 <- read.delim(file = "./rawdata/LPS_TNF_EX3.count.txt",header = T,stringsAsFactors = F,sep = "\t")
lib.name <- "/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/Reference/Library_original/mouse_geckov2_library_combine.csv"
lib <- read.delim(file = lib.name, header = F,sep = ",",stringsAsFactors = F)
merge.rawdata <- merge(rawdata.1,rawdata.2,by.x="sgRNA",by.y="sgRNA")
merge.rawdata <- merge(merge.rawdata,rawdata.3,by.x="sgRNA",by.y="sgRNA")
merge.rawdata.lib <- merge(merge.rawdata,lib,by.x = "sgRNA",by.y = "V1")
lib <- merge.rawdata.lib[c("sgRNA","V2","V3")]
colnames(lib) <- lib.colnames
lib$sgRNA_ID <- paste0("sgRNA",seq(1,nrow(lib),1))

rawdata <- merge.rawdata.lib[-c(9,14,19,20)]
colnames(rawdata)[c(1,2)] <- rawdata.rowname12
#file fifteen####
rawdata <- read.delim(file = "./rawdata/2CT_screens_ReadCountTable.txt",header = T,stringsAsFactors = F,sep = "\t")
rawdata <- rawdata[complete.cases(rawdata),]
lib <- rawdata[c(1,3,2)]
#colnames(lib) <- lib.colnames
#lib$sgRNA_ID <- paste0("sgRNA",seq(1,nrow(lib),1))
rawdata <- rawdata[-3]
#rawdata$UID <- paste0("sgRNA",seq(1,nrow(rawdata),1))
colnames(rawdata)[c(1,2)] <- rawdata.rowname12
#file sixteen####
rawdata.1 <- read.delim(file= "./rawdata/elife_47362_count.txt",sep = "\t",header = T,stringsAsFactors = F)
seq <- stringr::str_extract_all(rawdata.1$sgRNA,pattern = "[:upper:]{20}$",simplify = T)
locatnum <- stringr::str_locate(rawdata.1$sgRNA,pattern = "_[:upper:]")
locatnum <- locatnum[,1]
rawdata.1$sgRNA <- stringr::str_sub(rawdata.1$sgRNA,start = 1,end = locatnum-1)
rawdata.1$sgRNA_seq <- seq
lib.1 <- rawdata.1[c("sgRNA","sgRNA_seq","gene")]
rawdata.1 <- rawdata.1[-7]
rawdata.1$sgRNA <- paste0("sgRNA",seq(1,nrow(rawdata.1),1))
colnames(rawdata.1)[c(1,2)] <- rawdata.rowname12

rawdata.2 <- read.delim(file ="./rawdata/elife_47362_count_2.txt",sep = "\t",header = T,stringsAsFactors = F)
merge.rawdata.2 <- merge(rawdata.2,lib.1,by.x = "sgRNA",by.y = "sgRNA")
lib.2 <- merge.rawdata.2[c("sgRNA","gene.x","sgRNA_seq")]
rawdata.2 <- merge.rawdata.2[1:6]
colnames(rawdata.2)[c(1,2)] <- rawdata.rowname12

lib.1$sgRNA <- paste0("sgRNA",seq(1,nrow(lib.1),1))
colnames(lib.1) <- lib.colnames
colnames(lib.2) <- lib.colnames
lib.2$sgRNA_ID <- paste0("sgRNA",seq(1,nrow(lib.2),1))
#file seventeen####
rawdata <- read.delim(file ="./rawdata/Ritchie_count.txt",sep = "\t",header = T,stringsAsFactors = F)
lib <- read.delim(file = "./lib/Bassik_human_lib.tsv",sep = "\t",header = F,stringsAsFactors = F)
lib$V1 <- paste0(lib$V3,"_",lib$V1)
merge.lib.rawdata <- merge(rawdata, lib,by.x = "sgRNA",by.y = "V1")
lib <- merge.lib.rawdata[c("sgRNA","V2","gene.x.x")]
colnames(lib) <- lib.colnames
lib$sgRNA_ID <- paste0("sgRNA",seq(1,nrow(lib),1))
rawdata <- merge.lib.rawdata[1:6]
colnames(rawdata)[c(1,2)] <- rawdata.rowname12
rawdata$sgRNA <- paste0("sgRNA",seq(1,nrow(rawdata),1))
#file eighteen####
rawdata <- read.delim(file = "./rawdata/D1_div_nondiv.count_normalized.txt",sep = "\t",stringsAsFactors = F,header = T)
tmp_split <- stringi::stri_split_fixed(rawdata$sgRNA,pattern = "_",simplify = T)
lib <- as.data.frame(tmp_split)
rownames(lib) <- paste0("sgRNA",seq(1,nrow(lib),1))
lib <- lib %>% rownames_to_column(var = "sgRNA_ID")
colnames(lib) <- lib.colnames
rawdata$sgRNA <- paste0("sgRNA",seq(1,nrow(rawdata),1))

#file nineteen####

####Output####
setwd(Output.path)
dir.create(Output.file7)
setwd(Output.file7)
filename <- "Zymosin"
dir.create(filename)  #need to change each time
setwd(filename) #need to change each time
if(!file_test("-d","rawcount")){
  dir.create("rawcount")
}
if(!file_test("-d","lib")){
  dir.create("lib")
}
write.table(rawdata.Zymosin,"./rawcount/rawcount.txt",sep = "\t",quote = F,row.names = F) #need to change 
#write.table(l,file="./lib/library.csv",sep = "\t",quote = F,row.names = F)
setwd("..")
allvectors <- ls()
rm(list = allvectors[which(allvectors!="base.path" & allvectors !="Input.path" &
                             allvectors !="Output.path" & allvectors !="lib.colnames" &
                             allvectors != "rawdata.rowname12")])
file.copy("../../Reference/nonessential_ctrl_sgrna_list.txt","./lib/nonessential_ctrl_sgrna_list.txt")


#error checking: 
path <- "/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/Result/31365872_Benjamin E.Gewurz_CellRep_2019_1"
setwd(path)
rawcount <- read.table(file = "./rawcount/rawcount.txt",header = T,stringsAsFactors = F)
lib <- read.delim(file = "./lib/nonessential_ctrl_sgrna_list.txt",header = F,stringsAsFactors = F)
control_sgRNA <- rawcount[rawcount$Gene %in% lib$V1,]$sgRNA
#lib <- read.delim(file = "./lib/nonessential_ctrl_sgrna_list.txt",header = F,stringsAsFactors = F)
#lib <- lib[lib$V1 %in% rawcount[,2],]
write.table(control_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",row.names = F,col.names = F,quote = F)


#### reference:
# 1. config.yaml 中没有制定normalization 的方式，导致mageck 中没有进行normalization。 需要指定 control sgrna
# 做normalization。 默认使用350 个基因对应的sgrna做normalization，可以在github中进行做。
# 
# 2. sample,linrary, publication,cellline 信息需要严格按照格式要求放到TXT文件里面，TXT文件放到annotation 文件夹下面。
# 
# 3.默认不需要设置design matrix，只要设置day0.label就可以了。（所以不需要设置两个文件夹，需要一个文件夹就可以了）

#check the quantile of the lib and rawcount:
#this code is used to check if lib or rawcount is right.
lib.colnames <- c("sgRNA_ID","sgRNA_Seq","Gene")
for (i in as.vector(Output.filename)){
  print(i)
  setwd(i)
  lib <- read.table("./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
  rawcount <- read.table("./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
  colnames(lib) <- lib.colnames
  write.table(lib,file = "./lib/library.csv",sep = "\t",row.names = F,col.names = T,quote = F)
  #print(head(lib))
  #print(head(rawcount[c(1,2)]))
  setwd("..")
}

#refine lib and rawcount and generate nonessential_file:
#file1####
nonessential_gene <- read.delim("../../Reference/nonessential_ctrl_sgrna_list.txt",header = F,sep = "\t",stringsAsFactors = F)
Output.filename <- list.files()
file <- as.vector(Output.filename[2])
setwd(file)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
rawcount$sgRNA <- paste0("sgRNA",seq(1,nrow(rawcount),1))
write.table(rawcount, file = "./rawcount/rawcount.txt",sep = "\t",quote = F,row.names = F)
#nonessential_sgRNA
nonessential_sgRNA <- lib$sgRNA_ID[which(str_to_upper(lib$Gene) %in% nonessential_gene$V1)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)

#file2,3,4,5####
#file2 <- "28723893_W. Nicholas Haining_Nature_2017_1"
file2 <- as.vector(Output.filename[2])
setwd(file2)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.table(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
#rawcount$sgRNA <- paste0("sgRNA",seq(1,nrow(rawcount),1))
genes <- lib$Gene
#convert mouse gene symbol into human gene symbol
library(biomaRt)
human = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl")  
genes = getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",values = genes,
               mart = mouse,attributesL = c("hgnc_symbol"),martL = human , uniqueRows = T)
genes.intersect <- intersect(genes$HGNC.symbol,nonessential_gene$V1)
genes.original <- genes$MGI.symbol[which(genes$HGNC.symbol %in% genes.intersect)]
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% genes.original)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
colnames(lib)[2] <- "sgRNA_Seq"
setwd("..")

#file6:####
#delete rawcount with nan:
#28783722_Nicholas P. Restifo_Nature_2017
file6 <- Output.filename[6]
rawcount  <- rawcount[complete.cases(rawcount),]
rownames(lib) <- lib$sgRNA_ID
lib <- lib[rawcount$sgRNA,]
write.table(lib,file = "./lib/library.csv",quote = F,col.names = T,row.names = F)
write.table(rawcount,file = "./rawcount/rawcount.txt",quote = F,col.names = T,row.names = F)

#file7:####
Output.filename <- list.files()
file7 <- Output.filename[7]
setwd(file7)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
rawcount$sgRNA <- lib$sgRNA_ID
write.table(rawcount,file = "./rawcount/rawcount.txt",quote = F,col.names = T,row.names = F,sep = "\t")
nonessential_gene <- read.delim("../../Reference/nonessential_ctrl_sgrna_list.txt",header = F,sep = "\t",stringsAsFactors = F)
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% nonessential_gene$V1)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
setwd("..")

#file9####
Output.filename <- list.files()
file9 <- Output.filename[9]
setwd(file9)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
rawcount$sgRNA <- lib$sgRNA_ID
write.table(rawcount,file = "./rawcount/rawcount.txt",quote = F,col.names = T,row.names = F,sep = "\t")
library(biomaRt)
human = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl")  
genes = lib$Gene
genes = getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",values = genes,
               mart = mouse,attributesL = c("hgnc_symbol"),martL = human , uniqueRows = T)
genes.intersect <- intersect(genes$HGNC.symbol,nonessential_gene$V1)
genes.original <- genes$MGI.symbol[which(genes$HGNC.symbol %in% genes.intersect)]
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% genes.original)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
setwd("..")
#file10,11,12####
Output.filename <- list.files()
file10 <- Output.filename[12]
setwd(file10)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
rawcount$sgRNA <- lib$sgRNA_ID
genes = lib$Gene
genes = getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",values = genes,
               mart = mouse,attributesL = c("hgnc_symbol"),martL = human , uniqueRows = T)
genes.intersect <- intersect(genes$HGNC.symbol,nonessential_gene$V1)
genes.original <- genes$MGI.symbol[which(genes$HGNC.symbol %in% genes.intersect)]
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% genes.original)]

write.table(rawcount,file = "./rawcount/rawcount.txt",quote = F,col.names = T,row.names = F,sep = "\t")
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
setwd("..")

#file13####
Output.filename <- list.files()
file13 <- Output.filename[13]
setwd(file13)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
all(lib$Gene == rawcount$Gene) 
rawcount$sgRNA <- lib$sgRNA_ID
nonessential_gene <- read.delim("../../Reference/nonessential_ctrl_sgrna_list.txt",header = F,sep = "\t",stringsAsFactors = F)
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% nonessential_gene$V1)]
write.table(rawcount,file = "./rawcount/rawcount.txt",quote = F,col.names = T,row.names = F,sep = "\t")
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
setwd("..")
#file14####
Output.filename <- list.files()
file14 <- Output.filename[14]
setwd(file14)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
lib <- lib[c(1,3,2)]
colnames(lib) <- c("sgRNA_ID","sgRNA_Seq","Gene")
all(lib$Gene == rawcount$Gene) 
nonessential_gene <- read.delim("../../Reference/nonessential_ctrl_sgrna_list.txt",header = F,sep = "\t",stringsAsFactors = F)
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% nonessential_gene$V1)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(lib,file = "./lib/library.csv",sep = "\t",quote = F,row.names = F,col.names = T)
setwd("..")
#file15####
Output.filename <- list.files()
file15 <- Output.filename[15]
setwd(file15)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
all(lib$Gene == rawcount$Gene) 
rawcount$sgRNA = lib$sgRNA_ID
nonessential_gene <- read.delim("../../Reference/nonessential_ctrl_sgrna_list.txt",header = F,sep = "\t",stringsAsFactors = F)
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% nonessential_gene$V1)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(rawcount,file = "./rawcount/rawcount.txt",sep = "\t",quote = F,row.names = F,col.names = T)
setwd("..")
#file16####
Output.filename <- list.files()
file16 <- Output.filename[16]
setwd(file16)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
all(lib$Gene == rawcount$Gene) 
#rawcount$sgRNA = lib$sgRNA_ID
genes = lib$Gene
genes = getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",values = genes,
               mart = mouse,attributesL = c("hgnc_symbol"),martL = human , uniqueRows = T)
genes.intersect <- intersect(genes$HGNC.symbol,nonessential_gene$V1)
genes.original <- genes$MGI.symbol[which(genes$HGNC.symbol %in% genes.intersect)]
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% genes.original)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(rawcount,file = "./rawcount/rawcount.txt",sep = "\t",quote = F,row.names = F,col.names = T)
setwd("..")

#file17####
Output.filename <- list.files()
file17 <- Output.filename[17]
setwd(file17)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
all(lib$Gene == rawcount$Gene) 
#rawcount$sgRNA = lib$sgRNA_ID
genes = lib$Gene
genes = getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",values = genes,
               mart = mouse,attributesL = c("hgnc_symbol"),martL = human , uniqueRows = T)
genes.intersect <- intersect(genes$HGNC.symbol,nonessential_gene$V1)
genes.original <- genes$MGI.symbol[which(genes$HGNC.symbol %in% genes.intersect)]
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% genes.original)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(rawcount,file = "./rawcount/rawcount.txt",sep = "\t",quote = F,row.names = F,col.names = T)
setwd("..")

#file18####
Output.filename <- list.files()
file18 <- Output.filename[18]
setwd(file18)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
all(lib$Gene == rawcount$Gene) 
#rawcount$sgRNA = lib$sgRNA_ID
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% nonessential_gene$V1)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(rawcount,file = "./rawcount/rawcount.txt",sep = "\t",quote = F,row.names = F,col.names = T)
setwd("..")

#file19####
Output.filename <- list.files()
file19 <- Output.filename[19]
setwd(file19)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
all(lib$Gene == rawcount$Gene) 
#rawcount$sgRNA = lib$sgRNA_ID
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% nonessential_gene$V1)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(rawcount,file = "./rawcount/rawcount.txt",sep = "\t",quote = F,row.names = F,col.names = T)
setwd("..")

#file20####
Output.filename <- list.files()
file20 <- Output.filename[20]
setwd(file20)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp <- lib[rawcount$sgRNA,]
all(tmp$Gene == rawcount$Gene) 
setwd("..")
#file21,22####
Output.filename <- list.files()
file21 <- Output.filename[22]
setwd(file21)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
all(lib$Gene == rawcount$gene)
colnames(rawcount)[c(2)] <- "Gene"
genes = lib$Gene
genes = getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",values = genes,
               mart = mouse,attributesL = c("hgnc_symbol"),martL = human , uniqueRows = T)
genes.intersect <- intersect(genes$HGNC.symbol,nonessential_gene$V1)
genes.original <- genes$MGI.symbol[which(genes$HGNC.symbol %in% genes.intersect)]
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% genes.original)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
setwd("..")

#file23,24####
Output.filename <- list.files()
file23 <- Output.filename[23]
setwd(file23)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
all(lib$Gene == rawcount$gene)
#colnames(rawcount)[c(2)] <- "Gene"
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% nonessential_gene$V1)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
setwd("..")
#file24
Output.filename <- list.files()
file24 <- Output.filename[24]
setwd(file24)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
rawcount$sgRNA = paste0("sgRNA",seq(1,nrow(rawcount),1))
lib <- lib[c(1,3,2)]
colnames(lib) <- c("sgRNA_ID","sgRNA_Seq","Gene")
all(lib$Gene == rawcount$gene)
#colnames(rawcount)[c(2)] <- "Gene"
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% nonessential_gene$V1)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(rawcount,file = "./rawcount/rawcount.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(lib,file = "./lib/library.csv",sep = "\t",quote = F,row.names = F,col.names = T)
setwd("..")

#file25,26####
Output.filename <- list.files()
file25 <- Output.filename[26]
setwd(file25)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
#rawcount$sgRNA = paste0("sgRNA",seq(1,nrow(rawcount),1))
#lib <- lib[c(1,3,2)]
#colnames(lib) <- c("sgRNA_ID","sgRNA_Seq","Gene")
all(lib$Gene == rawcount$gene)
genes = lib$Gene
genes = getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",values = genes,
               mart = mouse,attributesL = c("hgnc_symbol"),martL = human , uniqueRows = T)
genes.intersect <- intersect(genes$HGNC.symbol,nonessential_gene$V1)
genes.original <- genes$MGI.symbol[which(genes$HGNC.symbol %in% genes.intersect)]
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% genes.original)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
setwd("..")

#file27####
Output.filename <- list.files()
file27 <- Output.filename[27]
setwd(file27)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
all(lib$Gene == rawcount$gene)
colnames(rawcount)[c(2,3)] <- c("Gene","Noselect")
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% nonessential_gene$V1)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(rawcount,file = "./rawcount/rawcount.txt",sep = "\t",quote = F,row.names = F,col.names = T)
#nonessential sgRNA human:
nonessential_gene <- read.delim("../../Reference/nonessential_ctrl_sgrna_list.txt",header = F,sep = "\t",stringsAsFactors = F)
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% nonessential_gene$V1)]
#nonessential sgRNA mouse:
library(biomaRt)
human = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl")  
genes = lib$Gene
genes = getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",values = genes,
               mart = mouse,attributesL = c("hgnc_symbol"),martL = human , uniqueRows = T)
genes.intersect <- intersect(genes$HGNC.symbol,nonessential_gene$V1)
genes.original <- genes$MGI.symbol[which(genes$HGNC.symbol %in% genes.intersect)]
nonessential_sgRNA <- lib$sgRNA_ID[which(lib$Gene %in% genes.original)]
write.table(nonessential_sgRNA,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,row.names = F,col.names = F)



#create annotation file:
#create annotation file####
setwd(base.path)
anno_file <- read.delim(file="./Reference/ImmuneScreens_Annotation.tsv",header = T,stringsAsFactors = F,sep = "\t")
#file tree:
##-------pulication.txt####
#### ID,PMID,Journal,Year,First_Author,Last_Authro,Title
colname_choose <- c("PMID","Journal","Year","FirstAuthor","Last_Author","Title")
publication <- anno_file[colname_choose]
write.table(publication,file = "./Reference/publication.txt",col.names = T,row.names = F,sep = "\t",quote = F)
#new output filename:
library(stringi)
Output.filename <- list.files(".",recursive = F)
pubmadid <- stri_split_fixed(Output.filename,pattern = "_",simplify = T)
pubmadid <- pubmadid[,1]
tmp_df <- data.frame()
for(pubmad in pubmadid){
  col.num <- which(publication$PMID == pubmad)[1]
  col.tmp <- publication[col.num,]
  tmp_df <- rbind(tmp_df,col.tmp)
}
Output.filename <- data.frame(Output.filename)
df <- cbind(Output.filename,tmp_df)
write.table(df,file = "../Reference/publication.new.txt",sep = "\t",col.names = T,row.names = F,quote = F)
publication <- read.delim(file = "./Reference/publication.new.txt",sep = "\t",header = T,stringsAsFactors = F)
setwd("./Output/")
Output.filename <- list.files()
Output.filename <- Output.filename[-c(14,29,30)]
library(tibble)
for(i in Output.filename){
  setwd(i)
  if(!file_test("-d","annotation")){
    dir.create("annotation")
  }
  tmp <- publication[which(publication$Output.filename == i),]
  tmp <- tmp[-1]
  tmp <- add_column(tmp,ID =1,.before = 1)
  #print(tmp)
  write.table(tmp,"./annotation/publication.txt",sep = "\t",col.names = T,row.names = F,quote = F)
  setwd("..")
}

##-------cellline.txt####
#### context: Cell_line, Cell_line_Search, Species,Tissue
anno_file$Condition[36:38] <- "NA"
cell_line_tmp <- data.frame()
for(pubmad in pubmadid){
  col.num <- which(anno_file$PMID == pubmad)[1]
  col.tmp <- anno_file[col.num,c(5,6)]
  cell_line_tmp <- rbind(cell_line_tmp,col.tmp)
}
cell_line_tmp <- cbind(Output.filename,cell_line_tmp)
write.table(cell_line_tmp,file = "../Reference/cellline.txt",sep = "\t",col.names = T,row.names = F,quote = F)
#split cell line into annotation folder:
cellline <- readxl::read_xlsx("../Reference/cell.xlsx")
cellline <- cellline[5:8]
id = 1 
for(i in Output.filename){
  setwd(i)
  if(!file_test("-d","annotation")){
    dir.create("annotation")
  }
  tmp <- cellline[id,]
  id = id +1 
  #print(tmp)
  write.table(tmp,"./annotation/cellline.txt",sep = "\t",col.names = T,row.names = F,quote = F)
  setwd("..")
}
##-------sample.txt####
#### Sample_ID, Dir_Name, Publication, Organism, Initial_Condition,Cell_Line, Type,Library
Sample_ID <- seq(1,27,1)
Dir_Name <- Output.filename
Publication <- rep("pmid",27)
Organism <- rep("NA",27)
Initial_Condition <- rep("HL60_initial,KBM7_initial",27)
Cell_Line <- rep("NA",27)
Type <- rep("ko",27)
Library <- rep(1,27)
Norm_method <- rep("350_noness",27)
Treatment <- rep("NA",27)
Other_Conditions <- rep("NA",27)
Culture_Days <- rep("NA",27)
Sample <- rep("NA",27)
Is_Initial_Condition <- rep("YES",27)
Source <- rep("counts",27)
Contact <- rep("NA",27)
Public_Status <- rep("public",27)
sample <- data.frame(Sample_ID,Dir_Name,Publication,Organism,Initial_Condition,Cell_Line,Type,
                     Library,Norm_method,Treatment,Other_Conditions,Culture_Days,Sample,Is_Initial_Condition,
                     Source,Contact,Public_Status)
write.table(sample,file = "../Reference/sample.txt",col.names = T,row.names = F,sep = "\t",quote = F)

sample <- readxl::read_xlsx("../Reference/sample.xlsx")
Output.filename <- sample$Dir_Name
for(i in seq(1,27,1)[-13]){
  setwd(Output.filename[i])
  print(Output.filename[i])
  rawcount <- read.delim(file = "./rawcount/rawcount.txt",header = T,sep = "\t",stringsAsFactors = F)
  rep_time <- ncol(rawcount)-2
  Sample_ID <- seq(1,rep_time,1)
  Dir_Name <- rep(Output.filename[i],rep_time)
  Publication <- rep("pmid",rep_time)
  Organism <- rep(sample$Organism[i],rep_time)
  Initial_Condition <- rep(sample$Initial_Condition[i],rep_time)
  Cell_Line <- rep(sample$Cell_Line[i],rep_time)
  Type <- rep("ko",rep_time)
  Library <- rep(1,rep_time)
  Norm_method <- rep("350_noness",rep_time)
  Treatment <- rep(sample$Treatment[i],rep_time)
  Other_Conditions <- rep(sample$Other_Conditions[i],rep_time)
  Culture_Days <- rep(sample$Culture_Days[i],rep_time)
  Sample <- colnames(rawcount)[3:ncol(rawcount)]
  Source <- rep("counts",rep_time)
  Contact <- rep("JijunYu",rep_time)
  Public_Status <- rep("public",rep_time)
  all_Initial <- strsplit(Initial_Condition[1],split = ",")[[1]]
  Is_Initial_Condition <- Sample  %in% all_Initial
  from <- c(FALSE,TRUE)
  to <- c("NO","YES")
  Is_Initial_Condition <- plyr::mapvalues(Is_Initial_Condition,from = from,to=to)
  Tmp.sample <- data.frame(Sample_ID,Dir_Name,Publication,Organism,Initial_Condition,Cell_Line,Type,
                       Library,Norm_method,Treatment,Other_Conditions,Culture_Days,Sample,Is_Initial_Condition,
                       Source,Contact,Public_Status)
  write.table(Tmp.sample,file = "./annotation/sample.txt",col.names = T,row.names = F,sep = "\t",quote = F)
  setwd("..")
}
#for specific 30397336:
setwd("./30397336_Michael C. Bassik_NatGenet_2018/")
file.list <- list.files()
for (i in file.list){
  setwd(i)
  rawcount <- read.delim(file = "./rawcount/rawcount.txt",header = T,sep = "\t",stringsAsFactors = F)
  rep_time <- ncol(rawcount)-2
  Sample_ID <- seq(1,rep_time,1)
  Dir_Name <- rep("30397336_Michael C. Bassik_NatGenet_2018",rep_time)
  Publication <- rep("pmid",rep_time)
  Organism <- rep(sample$Organism[13],rep_time)
  Initial_Condition <- grep("FT",colnames(rawcount),value = T)
  Initial_Condition <- paste0(Initial_Condition,collapse = ",")
  Cell_Line <- rep(sample$Cell_Line[13],rep_time)
  Type <- rep("ko",rep_time)
  Library <- rep(1,rep_time)
  Norm_method <- rep("350_noness",rep_time)
  Treatment <- rep(i,rep_time)
  Other_Conditions <- rep(sample$Other_Conditions[13],rep_time)
  Culture_Days <- rep(sample$Culture_Days[13],rep_time)
  Sample <- colnames(rawcount)[3:ncol(rawcount)]
  Source <- rep("counts",rep_time)
  Contact <- rep("JijunYu",rep_time)
  Public_Status <- rep("public",rep_time)
  all_Initial <- strsplit(Initial_Condition[1],split = ",")[[1]]
  Is_Initial_Condition <- Sample  %in% all_Initial
  from <- c(FALSE,TRUE)
  to <- c("NO","YES")
  Is_Initial_Condition <- plyr::mapvalues(Is_Initial_Condition,from = from,to=to)
  Tmp.sample <- data.frame(Sample_ID,Dir_Name,Publication,Organism,Initial_Condition,Cell_Line,Type,
                           Library,Norm_method,Treatment,Other_Conditions,Culture_Days,Sample,Is_Initial_Condition,
                           Source,Contact,Public_Status)
  write.table(Tmp.sample,file = "./annotation/sample.txt",col.names = T,row.names = F,sep = "\t",quote = F)
  setwd("..")
}
##-------library.txt####
#### context: ID, Library, Addgene_ID,Type, Species,Lentiviral Generation,PI, Genes_targeted, gRNAs_per_gene,Total_gRNAs
#reference:https://bitbucket.org/weililab/crisp-view/wiki/incorporatedata
ID <- seq(1,27,1)
Library <- rep("GeCKO",27)
AddGene_ID <- rep("NA",27)
Type <- rep("Knockout",27)
Species <- rep("Human",27)
Lentiviral_Generation <- rep("3rd",27)
PI <- rep("Zhang F",27)
Genes_Targeted <- rep("Varies",27)
gRNAs_per_gene <- rep(10,27)
Total_gRNAs <- rep("Varies",27)
lib <- data.frame(ID,Library,AddGene_ID,Type,Species,Lentiviral_Generation,
                  PI,Genes_Targeted,gRNAs_per_gene,Total_gRNAs)
write.table(lib,file = "../Reference/library.txt",sep = "\t",col.names = T,row.names = F,quote = F)

merge.lib.sample.celline <- cbind(lib,sample,cell_line_tmp)
merge.lib.sample.celline$Output.filename
merge.lib.sample.celline$PMID <- stri_split_fixed(as.vector(merge.lib.sample.celline$Output.filename),pattern = "_",simplify = T)[,1]

write.table(merge.lib.sample.celline,file = "../Reference/merge.lib.sample.cellline.txt",sep = "\t",
            col.names = T,row.names = F,quote = F)

#split the library.txt file into each folder:
#A <-iconv("./Reference/library.csv","WINDOWS-1252","UTF-8") 
lib <- readxl::read_xlsx("./Reference/library.xlsx")
Output.filename <- lib$Output.filename
lib <- lib[-2]
id = 1 
for(i in Output.filename){
  setwd(i)
  if(!file_test("-d","annotation")){
    dir.create("annotation")
  }
  tmp <- lib[which(lib$ID == id),]
  id = id +1 
  #print(tmp)
  write.table(tmp,"./annotation/library.txt",sep = "\t",col.names = T,row.names = F,quote = F)
  setwd("..")
}
#Sys.setlocale("LC_ALL","English")
#Sys.setlocale("LC_CTYPE","en_US.UTF-8")

setwd("..")
folder.list <- list.files()
nonessential_gene <- read.delim("../../Reference/nonessential_ctrl_sgrna_list.txt",sep = "\t",header = F)
for (i in folder.list){
  #i = folder.list[4]
  setwd(i)
  rawcount.tmp  <- read.delim("./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
  lib.tmp <- read.delim("./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
  print(all(rawcount.tmp$sgRNA == lib.tmp$sgRNA_ID))
  rawcount.tmp$sgRNA <- lib.bigbead$sgRNA_ID
  zero_row <- which(rowSums(rawcount.tmp[3:6]) == 0)
  rawcount.tmp <- rawcount.tmp[-zero_row,]
  print(dim(rawcount.tmp))
  write.table(rawcount.tmp,file = "./rawcount/rawcount.txt",sep = "\t",col.names = T,row.names = F,quote = F)
  setwd("..")
}
rawcount.bigbead  <- read.delim("./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
lib.bigbead <- read.delim("./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawcount.bigbead$sgRNA <- lib.bigbead$sgRNA_ID
nonessential_sgRNA.bigbead <- lib.bigbead$sgRNA_ID[which(lib.bigbead$Gene %in% nonessential_gene$V1)]
nonessential_sgRNA <- read.delim("./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t")


#-----refine the file 31452512 ####
setwd("/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/Output/31452512_Jeffrey\ Settleman_Elife_2019_2")
rawdata <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
sgRNA <- read.delim(file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",header = F,stringsAsFactors = F)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawdata <- rawdata[-grep(":",rawdata$Gene),]
rawdata <- rawdata[complete.cases(rawdata),]
zero_row <- which(rowSums(rawdata[3:6]) == 0)
rawdata <- rawdata[-zero_row,]
lib <- lib[which(lib$sgRNA_ID %in% rawdata$sgRNA),]
sgRNA.new <- intersect(sgRNA$V1,lib$sgRNA_ID)
write.table(rawdata,file = "./rawcount/rawcount.txt",sep = "\t",quote = F,col.names = T,row.names = F)
write.table(sgRNA.new,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,col.names = F,row.names = F)
write.table(lib,file = "./lib/library.csv",sep = "\t",quote = F,col.names = T,row.names = F)

#------refine the file  30559422 #####
setwd("/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/Output/30559422_E. Robert McDonald III_NatMed_2018/")
rawdata <- read.delim(file = "./rawcount/rawcount.txt",sep = "\t",header = T,stringsAsFactors = F)
sgRNA <- read.delim(file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",header = F,stringsAsFactors = F)
lib <- read.delim(file = "./lib/library.csv",sep = "\t",header = T,stringsAsFactors = F)
rawdata <- rawdata[-grep(":",rawdata$Gene),]
rawdata <- rawdata[complete.cases(rawdata),]
#zero_row <- which(rowSums(rawdata[3:6]) == 0)
#rawdata <- rawdata[-zero_row,]
lib <- lib[which(lib$sgRNA_ID %in% rawdata$sgRNA),]
sgRNA.new <- intersect(sgRNA$V1,lib$sgRNA_ID)
write.table(rawdata,file = "./rawcount/rawcount.txt",sep = "\t",quote = F,col.names = T,row.names = F)
write.table(sgRNA.new,file = "./lib/nonessential_ctrl_sgrna_list.txt",sep = "\t",quote = F,col.names = F,row.names = F)
write.table(lib,file = "./lib/library.csv",sep = "\t",quote = F,col.names = T,row.names = F)


#file copy annotation from Output folder into Result:
Output.folder <- "/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/Output/"
Result.folder <- "/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/Result/"
setwd(Output.folder)
filename <- list.files()
filename <- filename[-c(14,29,30)]
for(i in filename[-13]){
  print(i)
  from.folder <- paste0(Output.folder,i,"/annotation")
  to.folder <- paste0(Result.folder,i)
  file.copy(from=from.folder,to=to.folder,recursive = T)
  remove.file1 <- paste0(to.folder,"/config.yaml~")
  remove.file2 <- paste0(to.folder,"/config.yaml.bak")
  remove.file3 <- paste0(to.folder,"/#config.yaml#")
  remove.file4 <- paste0(to.folder,"/config.yaml.bak1")
  file.remove(remove.file1)
  file.remove(remove.file2)
  file.remove(remove.file3)
  file.remove(remove.file4)
}

#copy the annotation from output folder into result folder.
Output.folder <- "/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/Output/30397336_Michael C. Bassik_NatGenet_2018/"
Result.folder <- "/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/Result/30397336_Michael C. Bassik_NatGenet_2018/"
setwd(Output.folder)
filename <- list.files()
for(i in filename){
  print(i)
  from.folder <- paste0(Output.folder,i,"/annotation")
  to.folder <- paste0(Result.folder,i)
  file.copy(from=from.folder,to=to.folder,recursive = T)
  remove.file1 <- paste0(to.folder,"/config.yaml~")
  remove.file2 <- paste0(to.folder,"/config.yaml.bak")
  remove.file3 <- paste0(to.folder,"/#config.yaml#")
  remove.file4 <- paste0(to.folder,"/config.yaml.bak1")
  file.remove(remove.file1)
  file.remove(remove.file2)
  file.remove(remove.file3)
  file.remove(remove.file4)
}

#refine 30639098 
output.30639098 <- "/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/Result/30639098_Sarah\ A.\ Teichmann_Cell_2018"
setwd(output.30639098)
rawdata <- read.delim(file = "./rawcount/rawcount.txt",header = T,sep = "\t",stringsAsFactors = F)
rawdata.sc2a_irf4 <- rawdata[,grep("sc2a_irf4",colnames(rawdata))]
rawdata.sc2b_xbp1 <- rawdata[,grep("sc2b_xbp1",colnames(rawdata))]
rawdata.sc1_il13 <- rawdata[,grep("sc1_il13",colnames(rawdata))]
rawdata.s11_IL13 <- rawdata[,grep("s11_IL13",colnames(rawdata))]
rawdata.sx4_IRF4 <- rawdata[,grep("sx4_IRF4",colnames(rawdata))]
rawdata.sc3_gata3 <- rawdata[,grep("sc3_gata3",colnames(rawdata))]
#rawdata.sx2_STL <- rawdata[,grep("sx2_STL",colnames(rawdata))]
rawdata.s8a_STL <- rawdata[,grep("s8a_STL",colnames(rawdata))]
rawdata.s10_FOXP3 <- rawdata[,grep("s10_FOXP3",colnames(rawdata))]
rawdata.cd103 <- rawdata[,grep("cd103",colnames(rawdata))]
rawdata.s9_STG <- rawdata[,grep("s9_STG",colnames(rawdata))]
rawdata.s8b_XBP1 <- rawdata[,grep("s8b_XBP1",colnames(rawdata))]
rawdata.first_il4 <- rawdata[,grep("first_il4",colnames(rawdata))]

folder.create <- c("sc2a_irf4","sc2b_xbp1","sc1_il13","s11_IL13","sx4_IRF4",
                  "sc3_gata3","s8a_STL","s10_FOXP3","cd103","s9_STG",
                  "s8b_XBP1","first_il4")
for(i in folder.create){
  dir.create(i)
  setwd(i)
  #unlink("rawdata",recursive = T)
  dir.create("rawcount")
  setwd("rawcount")
  rawdata.tmp <- rawdata[,grep(i,colnames(rawdata))]
  rawdata.tmp <- cbind(rawdata[c("sgRNA","Gene")],rawdata.tmp)
  write.table(rawdata.tmp,file = "rawcount.txt",sep = "\t",quote = F,col.names = T,row.names = F)
  file.copy(from = "../../lib",to = "..",recursive = T)
  setwd("../..")
}

#annotation for 30639098
folder.create <- c("sc2a_irf4","sc2b_xbp1","sc1_il13","s11_IL13","sx4_IRF4",
                   "sc3_gata3","s8a_STL","s10_FOXP3","cd103","s9_STG",
                   "s8b_XBP1","first_il4")
basic.name <- "30639098_Sarah A. Teichmann_Cell_2018/"
for(i in folder.create){
  setwd(i)
  print(i)
  dir.create("annotation")
  rawcount <- read.delim(file = "./rawcount/rawcount.txt",header = T,sep = "\t",stringsAsFactors = F)
  rep_time <- ncol(rawcount)-2
  Sample_ID <- seq(1,rep_time,1)
  Dir.tmp <- paste0(basic.name,i)
  Dir_Name <- rep(Dir.tmp,rep_time)
  Publication <- rep("pmid",rep_time)
  Organism <- rep("Mouse",rep_time)
  if(i != "s11_IL13"){
    initial.tmp <- paste0(i,"_neg")
  }else{
    initial.tmp <- paste0(i,"_low")
  }
  Initial_Condition <- rep(initial.tmp,rep_time)
  Cell_Line <- rep("Th2",rep_time)
  Type <- rep("ko",rep_time)
  Library <- rep(1,rep_time)
  Norm_method <- rep("350_noness",rep_time)
  Treatment <- rep(i,rep_time)
  Other_Conditions <- rep(NA,rep_time)
  Culture_Days <- rep(7,rep_time)
  Sample <- colnames(rawcount)[3:ncol(rawcount)]
  Source <- rep("counts",rep_time)
  Contact <- rep("JijunYu",rep_time)
  Public_Status <- rep("public",rep_time)
  all_Initial <- strsplit(Initial_Condition[1],split = ",")[[1]]
  Is_Initial_Condition <- Sample  %in% all_Initial
  from <- c(FALSE,TRUE)
  to <- c("NO","YES")
  Is_Initial_Condition <- plyr::mapvalues(Is_Initial_Condition,from = from,to=to)
  Tmp.sample <- data.frame(Sample_ID,Dir_Name,Publication,Organism,Initial_Condition,Cell_Line,Type,
                           Library,Norm_method,Treatment,Other_Conditions,Culture_Days,Sample,Is_Initial_Condition,
                           Source,Contact,Public_Status)
  write.table(Tmp.sample,file = "./annotation/sample.txt",col.names = T,row.names = F,sep = "\t",quote = F)
  file.copy(from = "../annotation/cellline.txt",to = "./annotation/")
  file.copy(from = "../annotation/library.txt",to = "./annotation/")
  file.copy(from = "../annotation/publication.txt",to = "./annotation/")
  setwd("..")
}


#change the Dirname 
setwd("/Users/yujijun/Documents/01-Work/07-CRISPR_annotation/Result/30397336_Michael\ C.\ Bassik_NatGenet_2018")
folder.name <- list.files()
for(i in folder.name){
  setwd(i)
  setwd("annotation")
  sample <- read.delim("sample.txt",sep = "\t",stringsAsFactors = F,header = T)
  sample$Dir_Name <- paste0(sample$Dir_Name,"/",i)
  write.table(sample, file = "sample.txt",sep = "\t",row.names = F,quote = F,col.names = T)
  setwd("../..")
}


#delete all empty in folder name 
setwd("..")
all.folder <- list.files()
all.folder <- all.folder[-grep(".sh",all.folder)]
for(i in all.folder[14:27]){
  to <- stri_replace_all_fixed(i,".","")
  to <- stri_replace_all_fixed(to," ","")
  print(to)
  file.rename(from = i,to)
  setwd(to)
  sample <- read.delim("./annotation/sample.txt",sep = "\t",stringsAsFactors = F,header = T)
  sample$Dir_Name <- stri_replace_all_fixed(sample$Dir_Name,i,to)
  write.table(sample,file = "./annotation/sample.txt",sep = "\t",col.names = T,row.names = F,quote = F)
  setwd("..")
}

#propressing 30397336
setwd("./30397336_MichaelCBassik_NatGenet_2018/")
folder.3039 <- list.files()
for(i in folder.3039){
  setwd(i)
  sample <- read.delim("./annotation/sample.txt",sep = "\t",stringsAsFactors = F,header = T)
  to <- "30397336_MichaelCBassik_NatGenet_2018/"
  sample$Dir_Name <- paste0(to,i)
  write.table(sample,file = "./annotation/sample.txt",sep = "\t",col.names = T,row.names = F,quote = F)
  setwd("..")
}

#propressing 30639098
setwd("./30639098_SarahATeichmann_Cell_2018/")
folder.3063 <- list.files()
for(i in folder.3063){
  setwd(i)
  sample <- read.delim("./annotation/sample.txt",sep = "\t",stringsAsFactors = F,header = T)
  to <- "30639098_SarahATeichmann_Cell_2018/"
  sample$Dir_Name <- paste0(to,i)
  write.table(sample,file = "./annotation/sample.txt",sep = "\t",col.names = T,row.names = F,quote = F)
  setwd("..")
}

