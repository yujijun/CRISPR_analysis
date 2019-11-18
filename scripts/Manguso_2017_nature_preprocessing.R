# --------------
# Date:  2019-11-04 13:11:04
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
# About project: This is for data preprocessing.
#
require(tidyverse)
library("readxl")
require(dplyr)
library(tibble)
#input all excel data
basic_input = "/Users/yujijun/Documents/ImmuneScreens_v2/Manguso_2017_Nature_28723893/rawdata"
pool1 = read_excel(paste(basic_input,"pool_1.xlsx",sep = "/"),sheet = "Totals")
pool2 = read_excel(paste(basic_input,"pool_2.xlsx",sep = "/"),sheet = "Totals")
pool3 = read_excel(paste(basic_input,"pool_3.xlsx",sep = "/"),sheet = "Totals")
pool4 = read_excel(paste(basic_input,"pool_4.xlsx",sep = "/"),sheet = "averages")

pool1_choose <- pool1 %>% select(c("Construct IDs","In vitro","GVAX + PD-1 AVERAGE",
                                   "GVAX only AVERAGE","TCRa KO AVERAGE"))
pool2_choose <- pool2 %>% select(c("Construct IDs","in vitro","GVAX + PD1 AVERAGE",
                                   "GVAX only AVERAGE","TCRa KO AVERAGE"))
pool3_choose <- pool3 %>% select(c("Construct IDs","in vitro","GVAX + PD1 AVERAGE",
                                   "GVAX only AVERAGE","TCRa KO AVERAGE"))
pool4_choose <- pool4 %>% select(c("Construct IDs","in vitro","GVAX + PD-1 AVERAGE",
                                   "GVAX only AVERAGE","TCRa KO AVERAGE"))

colnames(pool2_choose) <- c("Construct IDs","In vitro","GVAX + PD-1 AVERAGE",
                            "GVAX only AVERAGE","TCRa KO AVERAGE")
colnames(pool3_choose) <- c("Construct IDs","In vitro","GVAX + PD-1 AVERAGE",
                            "GVAX only AVERAGE","TCRa KO AVERAGE")
colnames(pool4_choose) <- c("Construct IDs","In vitro","GVAX + PD-1 AVERAGE",
                            "GVAX only AVERAGE","TCRa KO AVERAGE")
pool_all <- bind_rows(pool1_choose,pool2_choose,pool3_choose,pool4_choose)
#tibble::separate(data,"Construct IDs",sep = ";",into=c("gRNA","gene"))
pool_all$sgRNA <- gsub(";\\w*","",pool_all$`Construct IDs`)
pool_all$gene <- gsub("\\w*;","",pool_all$`Construct IDs`)
pool_all <- pool_all %>% select(c("sgRNA","gene","In vitro","GVAX + PD-1 AVERAGE",
                                  "GVAX only AVERAGE","TCRa KO AVERAGE"))
colnames(pool_all) <- c("sgRNA","gene","In_vitro","GVAX_PD-1_AVERAGE",
                        "GVAX_only_AVERAGE","TCRa_KO_AVERAGE")

cal2 <- function(x){
  y = 2^x
  return(y)
}
pool_all[,3:6] <- pool_all[,3:6] %>% mutate_all(cal2)
write.table(pool_all,file=paste(basic_input,"pool_all_count.txt",sep = "/"),sep = "\t",quote = F, row.names = F)

