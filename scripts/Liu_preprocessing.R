# --------------
# Date:  2019-11-07 16:32:11
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project:
#
require(tidyverse)
library("readxl")
require(dplyr)
library(tibble)

input_path = "/Users/yujijun/Documents/ImmuneScreens_v2/Ritchie_2019_MolecularCell_31126740/rawdata/1-s2.0-S1097276519303594-mmc2.xlsx"
#input_path = "/Users/yujijun/Documents/ImmuneScreens_v2/Liu_2018_NatMed_30559422/rawdata/hsc4_ADAR_suppressor_screen_NatMed_Liuetal.xlsx"
input_file1 = read_excel(input_path,sheet = "ctrl_1",col_names = F)
input_file2 = read_excel(input_path,sheet = "ctrl_2",col_names = F)
input_file3 = read_excel(input_path,sheet = "cGAMP_1",col_names = F)
input_file4 = read_excel(input_path,sheet = "cGAMP_2",col_names = F)

input_file1$gene= str_split_fixed(input_file1$...1, "_", 4)[,2]
input_file1$sgRNA = paste(str_split_fixed(input_file1$...1, "_", 4)[,3],str_split_fixed(input_file1$...1, "_", 4)[,4],sep = "_")

input_file2$gene= str_split_fixed(input_file2$...1, "_", 4)[,2]
input_file2$sgRNA = paste(str_split_fixed(input_file2$...1, "_", 4)[,3],str_split_fixed(input_file2$...1, "_", 4)[,4],sep = "_")

input_file3$gene= str_split_fixed(input_file3$...1, "_", 4)[,2]
input_file3$sgRNA = paste(str_split_fixed(input_file3$...1, "_", 4)[,3],str_split_fixed(input_file3$...1, "_", 4)[,4],sep = "_")

input_file4$gene= str_split_fixed(input_file4$...1, "_", 4)[,2]
input_file4$sgRNA = paste(str_split_fixed(input_file4$...1, "_", 4)[,3],str_split_fixed(input_file4$...1, "_", 4)[,4],sep = "_")

input_file1_new <- input_file1[,c("sgRNA","gene","...2")]
colnames(input_file1_new) <- c("sgRNA","gene","ctrl1")

input_file2_new <- input_file2[,c("sgRNA","gene","...2")]
colnames(input_file2_new) <- c("sgRNA","gene","ctrl2")

input_file3_new <- input_file3[,c("sgRNA","gene","...2")]
colnames(input_file3_new) <- c("sgRNA","gene","cGAMP_1")

input_file4_new <- input_file4[,c("sgRNA","gene","...2")]
colnames(input_file4_new) <- c("sgRNA","gene","cGAMP_2")

input_merge_all_1 = merge(input_file1_new,input_file2_new, by.x = "sgRNA",by.y = "sgRNA")
input_merge_all_2 = merge(input_file3_new,input_file4_new, by.x = "sgRNA",by.y = "sgRNA")
input_merge_all = merge(input_merge_all_1,input_merge_all_2,by.x = "sgRNA",by.y = "sgRNA")

input_merge_all = input_merge_all[,c("sgRNA","gene.x.x","ctrl1","ctrl2","cGAMP_1","cGAMP_2")]


write.table(input_merge_all,file = "/Users/yujijun/Documents/ImmuneScreens_v2/Ritchie_2019_MolecularCell_31126740/rawdata/Ritchie_count.txt",sep = "\t",col.names = T,row.names = F,quote = F)







