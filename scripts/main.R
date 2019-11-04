library(crisprproject1)
####genelist####
genelist_freemam = c("Jak1","Jak2","Ifngr2","Tap1","Tap2","B2m","Ifngr1","Stat1")
genelist_MHCI_sorted = c("IFNGR2","JAK1","STAT1","IRF1","B2M","JAK2","IRF2","IFNGR1")
genelist_pech_NK_coculture = c("DCAF15","PTPN2","PTPN1","ICAM1","STAT1","JAK1","IFNGR2","JAK2","IFNGR1","PVRL2")
genelist_Pan_Pmel1 = c("Jak1","Jak2","Ifngr2","Tap1","Tap2","B2m","Ifngr1","Stat1","Rela","Pbrm1","Cd274","Ikbkb")
genelist_Pan_OT1 = c("Jak1","Jak2","Ifngr2","Tap1","Tap2","B2m","Ifngr1","Stat1","Gale","Gne","Rraga","Rnf31")
genelist_Patel_Mel624_EffectT_Top20 = c("HLA_A","B2M","MLANA","SOX10","hsa-mir-101-2","CD58","TAPBP","TAP2","TAF3","SRP54","MAD2L1","hsa-mir-548s",
                                        "RNPS1","PSMB5","RPL10A","CRKL","APLNR","TAP1","COL17A1","RPL23")
genelist_Patel_Mel624_EffectT_common = c("RPL22L1","RPL8","CAPN8","JAK2","TM4SF1","SOX10","AURKA")
Vredevoogd_IFNGR1mut.SKCM_MART1 = c("TRAF2","MAP3K7","BIRC2","TAP1","B2M","MLANA")
Shifrut_Primary_Human_T = c("LCP2","VAV1","TRMT112","CD247","RPP21","EXOSC6","HAPA8","MYC","LCK",
                            "CD5","MEF2D","UBASH3A","DGKZ","DGKA","SOCS1","TCEB2","CDKN1B","TMEM222","CBLB")
Shifrut_pilots = c("CBLB","CD5","CD3D","LCP2")
Shifrut_Human_T_FigureB = c("LCP2","LCK","VAV1","CD3D","CD247","GRAP2","CBLB","CD8A","DGKZ","GGKZ","MEF2D","TMEM222")
Henriksson_Thelper_Gata3 = c("Sec62","Gata3","Stat6")
Henriksson_Thelper_xbp1 = c("Fam71b","Cd200","Ypel3","Rps14","Arhgef18","Phpt1","Hsd17b7","Sesn2","Apc")
Jiang_Bcell_CD40 = c("FAS","IKBKB","FBXO11","CD40","CHUX","TRAF6","CELF1","TEAF5","SLC39A7")

####folder name####
Jiang_2019_CellRep_31365872

####set up parameters####
basic_path <- "/Users/yujijun/Documents/ImmuneScreens_v2/"
specific_path <- "Jiang_2019_CellRep_31365872/test"
output_path = "/Users/yujijun/Desktop/CRISPR/CRESPR_SCREEN_verified"
file_all = list.files(path = paste(basic_path,specific_path,sep = "/"))
file_choose = file_all[grepl(".gene_summary.txt", file_all)]
file_name = file_choose[2] #choose the file you like
verified_gene=Jiang_Bcell_CD40  #choose a gene list you would like to display
output_path ="/Users/yujijun/Desktop/CRISPR/CRESPR_SCREEN_verified"

#####display and save the plot with the function of volcanodisplayrep####
volcanodisplay(basic_path, specific_path, file_name, verified_gene, output_path)
