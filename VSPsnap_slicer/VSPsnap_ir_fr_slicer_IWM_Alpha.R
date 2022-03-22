#! /usr/bin/Rscript

options(stringsAsFactors = F)

# kriging (VSPipe)
library(ggplot2)
library(xlsx)
library(sp)
library(gstat)
library(grid)
library(RColorBrewer)
library(gridExtra)
#library(ggrepel)

source("dict_to_df.R")
source("VSPsnap_accessory_functions2.R")


# load data
inFiles<-read.delim("../VSPsnap/data/dates_8_22_21.tsv",header=F)$V1


for(i in 416:431){
  print(i)

dt<-strsplit(inFiles[i],"[.]")[[1]][1]

  print(dt)


genetab_all<-read.csv(paste0("../VSPsnap/data/cumulative_daily_AF_v4_lag_5_31/",inFiles[i]))
genetab_all$AF<-genetab_all$counts/sum(genetab_all$counts)

# filter < 3 counts
genetab_all=genetab_all[-which(genetab_all$counts<3),]

# filter >=3 countries
genetab=genetab_all[-which(genetab_all$countries<3),]

# remove FR=0 if any
if(length(which(genetab$fatality_rate==0))>0){
  genetab=genetab[-which(genetab$fatality_rate==0),]
}


# log transform IR and FR 
  genetab$IR<-log10(genetab$infection_rate)
  genetab$FR<-log10(genetab$fatality_rate)

# sequence
seqX<-genetab$Start/29903

# get IR axis 
IR<-genetab$IR

# FR axis
FR<-genetab$FR

krige_tab<-as.data.frame(cbind(seqX,IR,FR))
krige_tab=krige_tab[order(krige_tab$seqX),]


rownames(krige_tab)<-seq(1:nrow(krige_tab))

write.table(krige_tab,quote=F,row.names=F,sep="\t",file =paste0("krige_tab_genome_log_",dt,".tsv"))

# max euclidean distance between xy points - to compute width and cutoff
maxd<-max(dist(krige_tab[,1:2])) 

coordinates(krige_tab)=~seqX+IR

# generate variogram

lzn.vgm<-variogram(FR~1,krige_tab,cutoff=maxd*0.75)


# debug for "zero distance semivariance" issue

if(any(lzn.vgm$dist==0)){
  lzn.vgm=lzn.vgm[-which(lzn.vgm$dist==0),]
}

# fit model to best vgm
if(i>20){
 lzn.fit = fit.variogram(lzn.vgm, model =vgm(c("Exp","Sph","Gau","Exc","Mat","Ste","Cir","Lin","Bes","Pen","Per","Wav","Hol","Log","Pow","Spl","Leg")))#,debug.level = 2)
}else{
lzn.fit = fit.variogram(lzn.vgm, model =vgm("Sph"))
}

# bug with power model - range cannot exceed 2 - fallback to Sph model

if(lzn.fit$range[2]>2.0){
  lzn.fit = fit.variogram(lzn.vgm, model =vgm("Sph"))
}


# expand grid to nucleotide level on x (29903 nc), 100 pts on y (total 3M)

krige_grid<-expand.grid(seqX=seq(0,1,by=1/29903)[2:29904],IR=seq(min(genetab$IR),max(genetab$IR),by=0.01))  # 3M grid for compatibility with first pass


coordinates(krige_grid)=~seqX+IR
gridded(krige_grid)=TRUE

if(nrow(zerodist(krige_tab))!=0){
  krige_tab<-krige_tab[-zerodist(krige_tab)[,1],]
}

lzn.kriged = krige(FR~1, krige_tab,krige_grid, model = lzn.fit,nmin=5,nmax=50)

out<-krige.cv(FR~1,krige_tab,lzn.fit,nmin=5,nmax=50,nfold=5) # k-fold cv instead of LOO cuts down time drastically

out_v<-c(dt,mean(out$residual),mean(out$residual^2),mean(out$zscore^2),cor(out$observed, out$observed - out$residual),cor(out$observed - out$residual, out$residual))
names(out_v)<-c("Item","ME","MSPE","MSNE","cor_obs_pred","cor_pred_res")

#out_v=c(out_v,vgm_par)


print(round(as.numeric(out_v[5]),3))

lzn.kriged_df<-as.data.frame(lzn.kriged)
lzn.kriged_df$seqX=sapply(lzn.kriged_df$seqX,function(x) round(x,9))
write.table(lzn.kriged_df,quote=F,row.names=F,sep="\t",file =paste0("kriged_df_genome_cmt_ls_3M_",dt,".tsv"))

cov_proteins<-read.xlsx2("../VSPsnap/data/unipCov2Chain_UCSC_Wuh1_edited_CNBC_4.xlsx",sheetIndex = 1)
cov_proteins[,2:3]<-apply(cov_proteins[,2:3],2,as.numeric)

cols_out<-brewer.pal(n=11,name="Spectral")
pal <- colorRamp(rev(cols_out))

krige_df<-as.data.frame(krige_tab)

AF_tab<-as.data.frame(cbind(seqX,IR,FR,genetab$AF))
colnames(AF_tab)[4]<-"AF"


for(i in 1:nrow(cov_proteins)){

 getIWMZ(cov_proteins$chromStart[i],cov_proteins$chromEnd[i],cov_proteins$name[i],cov_proteins$L[i],29903,lzn.kriged_df,dt)
 getMapSlice(cov_proteins$chromStart[i],cov_proteins$chromEnd[i],cov_proteins$name[i],cov_proteins$L[i],29903,lzn.kriged_df,krige_df,genetab,AF_tab,dt,muts=T,minVar=F,IWM=T)

}


}












