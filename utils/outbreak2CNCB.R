# map Outbreak.info table of mutations to CNCB

setwd("/Users/sal/projects/BalchLab/Covid-19/VSPsnap/autoEWAD")

library(xlsx)

# read Outbreak.info file

mut_lambda<-read.delim("../data/Lambda_charMut_nextstrain.txt")
mut_mu<-read.delim("../data/Mu_charMut_nextstrain.txt")
mut_omicron<-read.delim("../data/Omicron_charMut_nextstrain.txt")

rewrite_mut<-function(mut){
  
  if(grepl("del",mut)==F){
  first<-substr(mut,1,1)
  last<-substr(mut,nchar(mut),nchar(mut))
  pos<-substr(mut,2,nchar(mut)-1)
  mut_out<-paste0(pos,first,">",last)
  }else{(mut_out<-mut)}
  return(mut_out)
}

mut_lambda$mut<-sapply(mut_lambda$aa,rewrite_mut)
mut_mu$mut<-sapply(mut_mu$aa,rewrite_mut)
mut_omicron$mut<-sapply(mut_omicron$aa,rewrite_mut)

# load CNCB mutations - latest csv files with more mutations

genetab_all<-read.csv("08_22_21.csv")

mut_lambda_table<-genetab_all[which(genetab_all$AA%in%mut_lambda$mut),]
write.xlsx(mut_lambda_table,file="lambda.xlsx")

mut_mu_table<-genetab_all[which(genetab_all$AA%in%mut_mu$mut),]
write.xlsx(mut_mu_table,file="mu.xlsx")

mut_omicron_table<-genetab_all[which(genetab_all$AA%in%mut_omicron$mut),]
write.xlsx(mut_omicron_table,file="omicron.xlsx")

# BA1 and BA2

mut_BA1<-read.delim("../data/BA1_nextstrain.txt")
mut_BA2<-read.delim("../data/BA2_nextstrain.txt")

mut_BA1$mut<-sapply(mut_BA1$aa,rewrite_mut)
mut_BA2$mut<-sapply(mut_BA2$aa,rewrite_mut)

# load CNCB mutations - latest csv files with more mutations

genetab_all<-read.csv("../data/cumulative_daily_AF_v4_lag_081421_011622/01_15_22.csv")

mut_BA1_table<-genetab_all[which(genetab_all$AA%in%mut_BA1$mut),]
write.xlsx(mut_BA1_table,file="BA1.xlsx")

mut_BA2_table<-genetab_all[which(genetab_all$AA%in%mut_BA2$mut),]
write.xlsx(mut_BA2_table,file="BA2.xlsx")







