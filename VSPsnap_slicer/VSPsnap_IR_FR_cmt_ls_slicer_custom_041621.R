
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

source("dict_to_df.R")
source("VSPsnap_log_accessory_functions.R")

# initialize out_tab
#out_v_tab<-vector()

# load data
inFiles<-read.delim("dates.tsv",header=F)$V1

# get mid-month dates
dates<-grep("15",inFiles)

# just march, may, sept, feb

for(i in dates[c(8,13)]){

  print(i)

  dt<-strsplit(inFiles[i],"[.]")[[1]][1]

  print(dt)

genetab_all<-read.csv(paste0("cumulative_daily_AF_v4_lag_2_25/",inFiles[i]))
genetab_all$AF<-genetab_all$counts/sum(genetab_all$counts)


# filter < 3 counts
  genetab_all=genetab_all[-which(genetab_all$counts<3),]

# filter >=3 countries
genetab=genetab_all[-which(genetab_all$countries<3),]

# log transform IR and FR (no 0-1 scaling)

  genetab$IR<-log10(genetab$infection_rate)
  genetab$FR<-log10(genetab$fatality_rate)

# scale IR and FR - local
#genetab$IR<-(genetab$infection_rate-min(genetab$infection_rate))/(max(genetab$infection_rate)-min(genetab$infection_rate))
#genetab$FR<-(genetab$fatality_rate-min(genetab$fatality_rate))/(max(genetab$fatality_rate)-min(genetab$fatality_rate))

# sequence
# Genome length (Wuhan-1): 1-29903

#seqX<-(genetab$X0-min(genetab$X0)+1)/(max(genetab$X0)-min(genetab$X0)+1)
seqX<-genetab$Start/29903

# get IR axis
IR<-genetab$IR

# FR axis
FR<-genetab$FR

krige_tab<-as.data.frame(cbind(seqX,IR,FR))
krige_tab=krige_tab[order(krige_tab$seqX),]

#if(outputTab==T){
#  write.table(krige_tab,quote=F,row.names=F,sep="\t",file = "cov2_pos_IR_FR_filt_3countries_fullgenome.tsv")
#}

rownames(krige_tab)<-seq(1:nrow(krige_tab))

write.table(krige_tab,quote=F,row.names=F,sep="\t",file =paste0("krige_tab_genome_cmt_log_",dt,".tsv"))

# max euclidean distance between xy points - to compute width and cutoff
maxd<-max(dist(krige_tab[,1:2]))

coordinates(krige_tab)=~seqX+IR

# generate variogram
#lzn.vgm<-variogram(FR~1,krige_tab,cutoff=0.1, width=0.001)#width=maxd/50,cutoff=maxd/3)
#vm <- vgm(psill=0.02928, model="Sph",range=0.006, nugget=0)

 lzn.vgm<-variogram(FR~1,krige_tab,cutoff=maxd*0.75)


# debug for "zero distance semivariance" issue

if(any(lzn.vgm$dist==0)){
  lzn.vgm=lzn.vgm[-which(lzn.vgm$dist==0),]
}

# fit model to vgm
  lzn.fit = fit.variogram(lzn.vgm, model =vgm(c("Exp","Sph","Gau","Exc","Mat","Ste","Cir","Lin","Bes","Pen","Per","Wav","Hol","Log","Pow","Spl","Leg")),debug.level = 2)



# get vgm parameters to be saved
# SSErr is a (weighted) sum of squared errors of the fitted model
#vgm_par<-c(as.character(lzn.fit$model[2]),lzn.fit$psill[1],lzn.fit$psill[2],lzn.fit$range[2],attr(lzn.fit,"SSErr"))
#names(vgm_par)<-c("vgm_model","Nug","psill","range","SSErr")

# print variogram

#  preds = variogramLine(lzn.fit, maxdist = max(lzn.vgm$dist))
#  vgm.plot<-ggplot(lzn.vgm,aes(x=dist,y=gamma))+geom_point() + geom_line(data=preds)
#  ggsave(vgm.plot,filename = paste0("Vgm_CoV2_genome_cmt_ls_",dt,".png"),dpi=100)


#IR_bottom<-signif(range(IR)[1],3)-0.01
#IR_top<-signif(range(IR)[2],3)+0.01

# expand grid to nucleotide level on x (29903 nc), 100 pts on y (total 3M)

krige_grid<-expand.grid(seqX=seq(0,1,by=1/29903)[2:29904],IR=seq(min(genetab$IR),max(genetab$IR),by=0.01))  # 3M grid for compatibility with first pass



coordinates(krige_grid)=~seqX+IR
gridded(krige_grid)=TRUE

if(nrow(zerodist(krige_tab))!=0){
  krige_tab<-krige_tab[-zerodist(krige_tab)[,1],]
}

lzn.kriged = krige(FR~1, krige_tab,krige_grid, model = lzn.fit,nmin=5,nmax=50)

# c("red","orange","yellow","green") brewer.pal(n=11,name="Spectral")

out<-krige.cv(FR~1,krige_tab,lzn.fit,nmin=5,nmax=50)

#out_v<-c(dt,mean(out$residual),mean(out$residual^2),mean(out$zscore^2),cor(out$observed, out$observed - out$residual),cor(out$observed - out$residual, out$residual))
#names(out_v)<-c("Item","ME","MSPE","MSNE","cor_obs_pred","cor_pred_res")

#out_v=c(out_v,vgm_par)


#print(round(as.numeric(out_v[5]),3)) # nmax=20: 0.423 /nmax=30:  0.434 /nmax=100:  0.443

lzn.kriged_df<-as.data.frame(lzn.kriged)
lzn.kriged_df$seqX=sapply(lzn.kriged_df$seqX,function(x) round(x,9))
write.table(lzn.kriged_df,quote=F,row.names=F,sep="\t",file =paste0("kriged_df_genome_cmt_ls_3M_",dt,".tsv"))

cov_proteins<-read.xlsx2("unipCov2Chain_UCSC_Wuh1_structures.xlsx",sheetIndex = 1)
cov_proteins[,2:3]<-apply(cov_proteins[,2:3],2,as.numeric)

cols_out<-brewer.pal(n=11,name="Spectral")
pal <- colorRamp(rev(cols_out))

krige_df<-as.data.frame(krige_tab)

AF_tab<-as.data.frame(cbind(seqX,IR,FR,genetab$AF))
colnames(AF_tab)[4]<-"AF"

# get Spike and NC only
for(i in 1:nrow(cov_proteins)){

 getMinVarZ(cov_proteins$chromStart[i],cov_proteins$chromEnd[i],cov_proteins$name[i],cov_proteins$L[i],29903,lzn.kriged_df,dt)
 getMapSlice(cov_proteins$chromStart[i],cov_proteins$chromEnd[i],cov_proteins$name[i],cov_proteins$L[i],29903,lzn.kriged_df,krige_df,genetab,AF_tab,dt,muts=T,minVar=T)

}


}

