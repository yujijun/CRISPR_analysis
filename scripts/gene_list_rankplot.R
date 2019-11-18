# library(MAGeCKFlute)
# input_path_sgrna <- "/Users/yujijun/Documents/ImmuneScreens_v2/Manguso_2017_Nature_28723893/test/TCRa_KO_vs_GVAX_PD-1.sgrna_summary.txt"
# input_path_gene <- "/Users/yujijun/Documents/ImmuneScreens_v2/Manguso_2017_Nature_28723893/test/TCRa_KO_vs_GVAX_PD-1.gene_summary.txt"
# dd.sgrna <- ReadsgRRA(input_path_sgrna)
# dd.rra <- ReadRRA(input_path_gene)
# dd.rra <- dd.rra[!grepl("SGR",dd.rra$Official),]
# geneList= dd.rra$LFC
# names(geneList) = dd.rra$Official
# p2 = RankView(geneList, top = 10, bottom = 10)
# print(p2)
#
# p2 = sgRankView(dd.sgrna, top = 20, bottom = 20)
# print(p2)

