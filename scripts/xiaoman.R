setwd("/Users/yujijun/Dropbox (Partners HealthCare)/ImmuneScreens_verified")
file <- list.files(path = ".",pattern = "*.sgrna_summary.txt",recursive = T)
library(MAGeCKFlute)
library(stringi)
library(stringr)
library(tibble)

ReadsgRRA <- function (sgRNA_summary)
{
  message(Sys.time(), " # Read sgRNA summary file ...")
  if (is.character(sgRNA_summary) && file.exists(sgRNA_summary)) {
    dd = read.table(file = sgRNA_summary, header = TRUE,
                    stringsAsFactors = FALSE,fill = T)
  }
  else if (is.data.frame(sgRNA_summary) && all(c("sgrna", "Gene",
                                                 "LFC", "FDR") %in% colnames(sgRNA_summary))) {
    dd = sgRNA_summary
  }
  else {
    stop("The parameter sgRNA_summary is below standard!")
  }
  dd = dd[, c("sgrna", "Gene", "LFC", "FDR")]
  return(dd)
}

dd.sgrna <- data.frame()
for( i in file.select){
  dd.sgrna_i <- ReadsgRRA(file[i])
  dd.sgrna_i$Gene <- str_to_upper(dd.sgrna_i$Gene)
  dd.sgrna_i_ESRRA <- dd.sgrna_i[dd.sgrna_i$Gene=="ESRRA",]
  print(dim(dd.sgrna_i_ESRRA))
  if(nrow(dd.sgrna_i_ESRRA) >=1){
    name <- basename(file[i])
    col_name <- str_split(name,pattern = "[.]")[[1]][1]
    col_name_all <- paste(col_name,seq(1,nrow(dd.sgrna_i_ESRRA)))
    dd.tmp <- dd.sgrna_i_ESRRA %>% add_column(sgdatabase = col_name_all,.before = "sgrna")
    dd.sgrna <- rbind(dd.sgrna,dd.tmp)
    print(dim(dd.sgrna))
  }
}
write.table(dd.sgrna,file = "/Users/yujijun/Desktop/ESRRA.sgRNA",sep = "\t",
            col.names = T,row.names = F,quote = F)

# dd.sgrna <- data.frame()
# for( i in file.select){
#   dd.sgrna_i <- ReadsgRRA(file[i])
#   dd.sgrna_i$Gene <- str_to_upper(dd.sgrna_i$Gene)
#   name <- basename(file[i])
#   col_name <- str_split(name,pattern = "[.]")[[1]][1]
#   dd.sgrna <- cbind()
# }
dd.sgrna %>%
ggplot(aes_string(x="expriment",y="scale_rank",color="type"))+
  geom_bar(stat = "identity",position = position_fill())+
  geom_hline(yintercept=c(0.8),size=.3,color="black",linetype="dashed")+
  theme_pubr()+labs(x="",y="Scaled Rank")

#how to exact the file of sgRNA
file.select.id<- c(4,24,25,19,20,26,27,5)
file.select <- file[file.select.id]
list_all <- list()
for(i in file.select){
  name = str_split_fixed(basename(i),pattern = "[.]",n=3)[,1]
  list_all[[name]] =ReadsgRRA(i)
}

