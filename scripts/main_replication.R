# --------------
# Date:  2019-11-03 22:19:54
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: this is the main function to draw volcano function.
library(crisprproject1)

####build the dd.rra mean file####
basic_path <- "/Users/yujijun/Documents/ImmuneScreens_v2/"
specific_path <- "Jiang_2019_CellRep_31365872/test"
output_path = "/Users/yujijun/Desktop/CRISPR/CRESPR_SCREEN_verified"
file_all = list.files(path = paste(basic_path,specific_path,sep = "/"))
file_choose = file_all[grepl(".gene_summary.txt", file_all)]
input_file_1 = paste(paste(basic_path,specific_path,sep = ""),file_choose[1],sep = "/")
input_file_2 = paste(paste(basic_path,specific_path,sep = ""),file_choose[2],sep = "/")
dd.rra_1 = ReadRRA(input_file_1)
dd.rra_2= ReadRRA(input_file_2)
dd.rra_merge = merge(dd.rra_1,dd.rra_2,by.x = "Official",by.y = "Official")
dd.rra_merge$LFC_mean = (dd.rra_merge$LFC.x + dd.rra_merge$LFC.y)/2
dd.rra_merge$FDR_mean = (dd.rra_merge$FDR.x + dd.rra_merge$FDR.y)/2
dd.rra = dd.rra_merge[,c("Official","LFC_mean","FDR_mean")]
colnames(dd.rra) = c("Official","LFC","FDR")

####set up other Parameters####
figure_title = "Jiang_Daudi_B_CD40_mean"
Jiang_Bcell_CD40 = c("FAS","IKBKB","FBXO11","CD40","CHUX","TRAF6","CELF1","TEAF5","SLC39A7")
output_path ="/Users/yujijun/Desktop/CRISPR/CRESPR_SCREEN_verified"

####display and save the plot with the function of volcanodisplayrep####
volcanodisplayrep(dd.rra,figure_title,Jiang_Bcell_CD40,output_path)
