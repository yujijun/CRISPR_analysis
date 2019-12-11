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
Jiang_Bcell_CD40 = c("FAS","IKBKB","FBXO11","CD40","CHUX","TRAF6","CELF1","TEAF5","SLC39A7","METTL14","KIAA1429","CHMP5","VPS25","NFKBIA")
Manguso_GVAX_PD1_TCRaKO = c("Stat1","Jak1","Ifnar1","Ifnar2","Eif2ak2","Ifngr1","Safb","Ripk1","Irf1","H2-T23")
Manguso_TCRaKO_invitro = c("Arf1","Adrbk1","Pdxk","Yipf3","Nf2","Chp1","Pon2","Kirrel")
Manguso_GVAX_TCRaKO = c("Ripk1","Ptpn2","Irf1","Calr","Stat1","Jak1","Eif2ak2")
Parnas_Inf_highvslow = c("Zfp36","Tlr4","Tnf","Ly96","Myd88","Cd14")
Parnas_Inf_zexian = c("Tnf", "Ly96", "Cd14", "Myd88", "Ticam2", "Dad1", "Ddost", "Ticam1", "Tirap", "Tlr4", "Zfp36", "Atxn7l3", "Rc3h1", "Yy1", "Dnttip1", "Rara", "Dusp1", "Eif5", "Stat5b", "Atp10b")
Haney_bigbead = c("NHLRC2", "TM2D2", "TLN1", "FERMT3", "LCMT1", "ITGAL", "ICAM1", "ITGB2", "SYS1-DBNDD2", "TLE3", "TM2D3", "ABI1", "PTPN7", "RPS6KA1", "STUB1", "NCKAP1L", "KAT6A", "PLEK", "MEF2D", "PRKCB")
Haney_igg = c("PRKCD","NCKAP1L","NANS","SLC35A2","SLC35A1","ICAM1","TSC2","PRKCB","GNE","ABI1","FNDC3B","UXS1","MAPK1","TSC1","PRKD2","DOCK2","RASA2")
Haney_comp = c("NCKAP1L","ICAM1","PRKCD","MAPK1","ABI1","RASA2","DOCK2","GNE","BRK1","CCNC","NANS","KCTD5","CREBBP","SLC35A2","CMAS","LZTR1","FBXO11","PTPN7","FNDC3B","PRKCB")
Haney_Myelin = c("MAPK1","ITGB2","NHLRC2","TLN1","PRKD2","PRKCB","PRKCD","SHOC2","FERMT3","SUCNR1","KAT6A","PLEK","RP11-45M22.4","LZTR1","NCKAP1L","PPP6R1","SRSF6","ARAF","FAM49B","KCTD5")
Haney_Myelin_r1 = c("MAPK1","ITGB2","NHLRC2","TLN1","PRKD2","PRKCD","FERMT3","NEDD8")

Haney_small = c("FERMT3","SLC35A2","NHLRC2","TLN1","PRKCD","TM2D2","PRKD2","TM2D3","NCKAP1L","ITGB2","PRKCB","ITGAL","PLEK","ICAM1","ARPC4-TTLL3","LCMT1","TM2D1","YPEL5","PPME1","RPS6KA1")
Hancy_pos = c("NHLRC2","ACTR2","NCKAP1L","SUPT20H","ACTR3","ABI1","ICAM1","DOCK2","CYFIP1","ARPC2","ARPC4-TTLL3","FLII","ARPC4","SPPL3","LCMT1","RAC1","ARPC3","MAPK1","AMBRA1","MYO9B")
Hancy_zymosan = c("MAPK1","PRKCD","NHLRC2","TLN1","PLEK","PRKD2","ITGB2","FERMT3","PRKCB","CERS2","NCKAP1L","SHOC2","SUCNR1","ZMYND8","PTPN7","ABI1","UBE2L3","JUNB","FADD","ELOVL1")
Hancy_mid = c("NHLRC2","NCKAP1L","ABI1","DOCK2","CD93","CYFIP1","LCMT1","ACTR2","BRK1","STUB1","ARPC2","RP11-45M22.4","TM2D2","PPME1","WASF2","LAMTOR4","ARPC3","ACTR3","CCNC","ARF1")
Burr_PDL1 = c("CMTM6","CD274","IFNGR1","JAK1","JAK2","STAT1","IFNGR2")
Kearney_B16F10_OT1ko = c("Samhd1","Isca2","Urm1","BC052040","Vmn2r95","Naa35","SGR000122364.1_XPR003.1","Rnf123","Cox4i1","Rtcb","Nrsn1","Pgm3","Ccni","Orc4","Pnpt1","Rlf")
Kearney_B16F10_WT = c("Ifngr1","Ifngr2","Jak1","Jak2","Stat1")
Kearney_MC38_NK = c("Ado","Casp8","Spata2","Gm13287","Tnfrsf1a","Capza2","Prss8","Smg9")
Kearney_MC38_OT1logG = c("Casp8","Jak1","Tnfrsf1a","Ado","B2m","Bmp3")
Kearney_MC38_OT1PDL1 = c("Casp8","Tnfrsf1a","Tap1","Tlx3","Ifngr1","Jak2")
Kearney_MC38_OT1Ko = c("Tnfrsf1a","Casp8","Tap1","B2m","Jak1","Ifngr1")
Ritchie_U937_cGAMP = c("SLC19A1","TBK1","IRF3","TMEM173")
Liu_HSC4_Dox = c("PKR","XPO5","IFNAR1","BDP1","BRF2","CHD8","GABPB1")
Liu_HSC4_DoxvsnoDox = c("RPK","MAVS","IFNAR1","STING")
Codina_P53 = c("Apc","Ap2s1","Ap2m1","Csnk1a1","Nf1","Nf2","Prkar1a","Pten","Tsc1","Tsc2","Vps45")

Done_cell_tumor = c("Cd247","Fam103a1","Havcr2","Tim3","Pdcd1","PD-1","Prkar1a","Stradb")
Done_invitro = c("Dhx37","Lyn","Odc1")
Mezzadra_HAP1_PDL1highvslow = c("CD274","JAK2","STAT1","JAK1","IFNGR2","IFNGR1","IRF1","CMTM6")
Burr_2019_K562_HLA_B = c("EED","SUZ12","MTF2","EZH2")
Burr_2019_K562_HLA_ABC = c("EED","SUZ12","MTF2","EZH2")
####folder name####
Freeman_2019_CellRep_31509742
Jiang_2019_CellRep_31365872
Manguso_2017_Nature_28723893
Parnas_2015_Cell_26189680
Vredevoogd_2019_Cell_31303383
Haney_2018_NatGenet_30397336
Pan_2018_Science_29301958
Burr_2017_Nature_28813417
Ritchie_2019_MolecularCell_31126740
Liu_2018_NatMed_30559422
Codina_2019_CellSyst_30797773
Dong_2019_Cell_31442407
Mezzadra_2017_Nature_28813410
Burr_2019_CancerCell_31564637
####set up parameters####

basic_path_1 <- "/Users/yujijun/Dropbox (Partners HealthCare)/ImmuneScreens"
basic_path <- "/Users/yujijun/Documents/ImmuneScreens_v2/"
specific_path <- "Freeman_2019_CellRep_31509742/test"
output_path = "/Users/yujijun/Desktop" #change the folder name
file_all = list.files(path = paste(basic_path,specific_path,sep = "/"))
file_choose = file_all[grepl("*.gene_summary.txt", file_all)]
file_name = file_choose[2] #choose the file you like
verified_gene= genelist_freemam#choose a gene list you would like to display

#input_file = paste(paste(basic_path,specific_path,sep = ""),file_name,sep = "/")
#dd.rra = ReadRRA(input_file)

####display and save the plot with the function of volcanodisplayrep####
volcanodisplay(basic_path, specific_path, file_name, verified_gene, output_path)


