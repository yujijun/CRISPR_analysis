# --------------
# Date:  2019-11-03 22:19:54
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: this is the main function to draw volcano function.
library(crisprproject1)

####build the dd.rra mean file####
basic_path <- "/Users/yujijun/Documents/ImmuneScreens_v2/"
specific_path <- "Parnas_2015_Cell_26189680/test"
output_path = "/Users/yujijun/Desktop/CRISPR/CRESPR_SCREEN_verified"
file_all = list.files(path = paste(basic_path,specific_path,sep = "/"))
file_choose = file_all[grepl(".gene_summary.txt", file_all)]
input_file_1 = paste(paste(basic_path,specific_path,sep = ""),file_choose[1],sep = "/")
input_file_2 = paste(paste(basic_path,specific_path,sep = ""),file_choose[2],sep = "/")
input_file_3 = paste(paste(basic_path,specific_path,sep = ""),file_choose[3],sep = "/")

dd.rra_1 = ReadRRA(input_file_1)
dd.rra_2= ReadRRA(input_file_2)
dd.rra_3= ReadRRA(input_file_3)
dd.rra_merge = merge(dd.rra_1,dd.rra_2,by.x = "Official",by.y = "Official")
df.rra_merge = merge(dd.rra_merge,dd.rra_3,by.x = "Official",by.y = "Official")
dd.rra_merge$LFC_mean = (dd.rra_merge$LFC.x + dd.rra_merge$LFC.y)/2
dd.rra_merge$FDR_mean = (dd.rra_merge$FDR.x + dd.rra_merge$FDR.y)/2
dd.rra = dd.rra_merge[,c("Official","LFC_mean","FDR_mean")]
colnames(dd.rra) = c("Official","LFC","FDR")

####set up other Parameters####
figure_title = "Parnas_Inf_highvslow_mean"
Parnas_Inf_zexian = c("Tnf", "Ly96", "Cd14", "Myd88", "Ticam2", "Dad1", "Ddost", "Ticam1", "Tirap", "Tlr4", "Zfp36", "Atxn7l3", "Rc3h1", "Yy1", "Dnttip1", "Rara", "Dusp1", "Eif5", "Stat5b", "Atp10b")

Jiang_Bcell_CD40 = c("FAS","IKBKB","FBXO11","CD40","CHUX","TRAF6","CELF1","TEAF5","SLC39A7","METTL14","KIAA1429","CHMP5","VPS25")
output_path ="/Users/yujijun/Desktop/CRISPR/CRESPR_SCREEN_verified"

####display and save the plot with the function of volcanodisplayrep####
volcanodisplayrep(dd.rra,figure_title,Parnas_Inf_zexian,output_path)




#Vredevoogd_2019_Cell_31303383 replication
####build the dd.rra mean file####
basic_path <- "/Users/yujijun/Documents/ImmuneScreens_v2/"
specific_path <- "Ritchie_2019_MolecularCell_31126740/test"
output_path = "/Users/yujijun/Desktop/CRISPR/CRESPR_SCREEN_verified"
file_all = list.files(path = paste(basic_path,specific_path,sep = "/"))
file_choose = file_all[grepl(".gene_summary.txt", file_all)]
input_file_1 = paste(paste(basic_path,specific_path,sep = ""),file_choose[1],sep = "/")
input_file_2 = paste(paste(basic_path,specific_path,sep = ""),file_choose[2],sep = "/")
#input_file_3 = paste(paste(basic_path,specific_path,sep = ""),file_choose[3],sep = "/")

dd.rra_1 = ReadRRA(input_file_1)
dd.rra_2= ReadRRA(input_file_2)
#dd.rra_3= ReadRRA(input_file_3)
dd.rra_merge = merge(dd.rra_1,dd.rra_2,by.x = "Official",by.y = "Official")
#df.rra_merge = merge(dd.rra_merge,dd.rra_3,by.x = "Official",by.y = "Official")
dd.rra_merge$LFC_mean = (dd.rra_merge$LFC.x + dd.rra_merge$LFC.y)/2
dd.rra_merge$FDR_mean = (dd.rra_merge$FDR.x + dd.rra_merge$FDR.y)/2
dd.rra = dd.rra_merge[,c("Official","LFC_mean","FDR_mean")]
colnames(dd.rra) = c("Official","LFC","FDR")

####set up other Parameters####
figure_title = "Ritchie_ctrlvscGAMP2_mean"
output_path ="/Users/yujijun/Desktop/CRISPR/CRESPR_SCREEN_verified"
Ritchie_U937_cGAMP = c("SLC19A1","TBK1","IRF3","TMEM173")
#Vredevoogd_IFNGR1mut.SKCM_MART1 = c("TRAF2","MAP3K7","BIRC2","TAP1","B2M","MLANA")

####display and save the plot with the function of volcanodisplayrep####
volcanodisplayrep(dd.rra,figure_title,Ritchie_U937_cGAMP,output_path)
