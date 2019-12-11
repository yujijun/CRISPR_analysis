library(crisprvarified)
#gene list
basic_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/data/"
output_path = "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/CRISPR_output_plot/Enrichment_score"
gene_list = read.table(paste(basic_path,"DGElist_zexian.txt",sep = "/"))
input_LFC = read.table(paste(basic_path,"LFC_maxscale.txt",sep = "/"))

geneset <- gene_list$V1
LFC_name <- rownames(input_LFC)
cohort <- colnames(input_LFC)
cohort <- gsub("_LFC","",cohort)
enrichment_score <- c()
enrichment_cohort <- c()

for (i in seq(1,ncol(input_LFC))){
  print(i)
  LFC <- input_LFC[,i]
  output <- geneset_enrichment(geneset,LFC,LFC_name,cohort[i])
  if (is.character(output$enrich_score)){
    print(out$enrich_score)
  }else{
    enrichment_score <- c(enrichment_score,output$enrich_score$ES)
    enrichment_cohort <- c(enrichment_cohort,cohort[i])
    #save the plot
    p <-output$enrich_plot
    figure_output = paste(paste(cohort[i],"enrichment_score",sep = "_"),".png",sep = "")
    png(paste(output_path,figure_output,sep = "/"),height = 12,width = 17,units ="cm",res=150)
    print(p)
    dev.off()
  }
}

#cohort <- cohort[-19]
enrich_df <- data.frame(es = enrichment_score, cohort = enrichment_cohort)

output_path <- "/Users/yujijun/Documents/01-Work/05-CRESPR_SCREEN/crisprproject1/CRISPR_output_plot/Enrichment_score"
write.table(enrich_df,file=paste(output_path,"enriment_df.txt",sep = "/"),sep = "\t",quote = F,col.names = F,row.names = F)
p <- singlegenebarplot(enrich_df, x="cohort",y = "es",fill = "#2F4F4F",sort.by.groups = T,ylab = "cohort",xlab = "enrichment_score",title = "barplot of enrichment score",legend.title = NA )
figure_output = paste("barplot of enrichment score",".png",sep = "")
png(paste(output_path,figure_output,sep = "/"),height = 12,width = 17,units ="cm",res=150)
print(p)
dev.off()







#The warning produced indicates that there are few genes
#that have the same fold change and so are ranked equally.
#fgsea with arbitrarily order determine which comes first
#in the ranked list. As long as this number is small it
#shouldn’t significantly effect the results. If the number
#is large something is suspicious about the fold change results.
