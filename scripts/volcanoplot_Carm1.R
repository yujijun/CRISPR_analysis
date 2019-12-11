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
# qc_matrix = read_excel("/Users/yujijun/Documents/ImmuneScreens_v2/Quanlity_control.xlsx",sheet = "Sheet2")
# folder_name_choosed = unique(qc_matrix[which(qc_matrix$PassQC == "Y"),"Dir_name"])
# author_name_choosed = str_split_fixed(folder_name_choosed$Dir_name,"_",n=2)[,1]

all_file = list.files(basic_path,pattern = "*.gene_summary.txt",recursive = T)
all_file_list = paste(basic_path, all_file,sep = "/")

author_name = str_split_fixed(all_file,"_",n=2)[,1]
control_treatment = str_split_fixed(basename(all_file),"[.]",n=2)[,1]
figure_name = paste(author_name, control_treatment,sep = "_")
verfied_gene = "CARM1"
output_path = "/Users/yujijun/Documents/ImmuneScreens_QCoutput/Carm1_output"
#draw volcano plot for all verified
for(i in seq(1,length(all_file_list),1)){
  print(figure_name[i])
  volcanodisplaysinglegene(all_file_list[i],verfied_gene,figure_name[i],output_path)
}

