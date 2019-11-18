# --------------
# Date:  2019-11-12 09:11:53
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project:This project is for Carm1 analysis
#
require(stringi)
require(stringr)
library(crisprvarified)

#basic parameter
basic_path = "/Users/yujijun/Documents/ImmuneScreens_verified"
all_file = list.files(basic_path,pattern = "*.gene_summary.txt",recursive = T)
all_file_list = paste(basic_path, all_file,sep = "/")

author_name = str_split_fixed(all_file,"_",n=2)[,1]
control_treatment = str_split_fixed(basename(all_file),"[.]",n=2)[,1]
figure_name = paste(author_name, control_treatment,sep = "_")
# #verfied_gene = toupper(c('Impdh2','Rasgrp2','Pdlim1','Klrb1c','Atp1b1',
#                          'Eomes','Igfbp4','Pim2','Hmgn1','Id3','Rps2',
#                          'Satb1','Acp5','Tubb5','Dapl1','Fam101b',
#                          'S1pr1','Ly6c2','Nsg2','Ccr7','Ighm','Bcl2',
#                          'Fcer1g','Lef1','Sell','Tcf7','Klf2'))
#verfied_gene = toupper(c("Tcf7","Sell","Lef","Ccr7"))
verfied_gene = toupper(c("Lgals3","Rgs16","Cd8a","Tigit","Nkg7","S100a4","Klrc1","Cxcr6","Pdcd1","Id2","Capg","Gzmb","Lag3","Tnfrsf9","Lgals1","Gzma","Ccl3","S100a6","Ccl4"))
output_path = "/Users/yujijun/Documents/ImmuneScreens_QCoutput/liyejie_neglist"
#draw volcano plot for all verified
for(i in seq(1,length(all_file_list),1)){
  print(figure_name[i])
  volcanodisplaymultigene(all_file_list[i],verfied_gene,figure_name[i],output_path)
}




