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

# initialize out_tab
out_v_tab<-vector()

# load data
inFiles<-read.delim("dates.tsv",header=F)$V1

for(i in %i:%i){
  print(i)

dt<-strsplit(inFiles[i],"[.]")[[1]][1]

  print(dt)

genetab_all<-read.csv(paste0("cumulative_daily1_7/",inFiles[i]))
genetab_all$AF<-genetab_all$counts/sum(genetab_all$counts)

# filter >=3 countries

genetab=genetab_all[-which(genetab_all$countries<3),]

# normalize IR and FR

genetab$IR<-(genetab$infection_rate-min(genetab$infection_rate))/(max(genetab$infection_rate)-min(genetab$infection_rate))
genetab$FR<-(genetab$fatality_rate-min(genetab$fatality_rate))/(max(genetab$fatality_rate)-min(genetab$fatality_rate))

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

write.table(krige_tab,quote=F,row.names=F,sep="\t",file =paste0("krige_tab_genome_cmt_ls_",dt,".tsv"))

# max euclidean distance between xy points - to compute width and cutoff
maxd<-max(dist(krige_tab[,1:2])) # 1.34

coordinates(krige_tab)=~seqX+IR

# generate variogram
#lzn.vgm<-variogram(FR~1,krige_tab,cutoff=0.1, width=0.001)#width=maxd/50,cutoff=maxd/3)
#vm <- vgm(psill=0.02928, model="Sph",range=0.006, nugget=0)

lzn.vgm<-variogram(FR~1,krige_tab) #cutoff=0.1, width=0.001)


# debug for "zero distance semivariance" issue

if(any(lzn.vgm$dist==0)){
  lzn.vgm=lzn.vgm[-which(lzn.vgm$dist==0),]
}

# fit model to best vgm among three model options
#lzn.fit = fit.variogram(lzn.vgm, model = vgm("Sph"))
lzn.fit = fit.variogram(lzn.vgm, model = vgm("Exp"))

# get vgm parameters to be saved
# SSErr is a (weighted) sum of squared errors of the fitted model
vgm_par<-c(as.character(lzn.fit$model[2]),lzn.fit$psill[1],lzn.fit$psill[2],lzn.fit$range[2],attr(lzn.fit,"SSErr"))
names(vgm_par)<-c("vgm_model","Nug","psill","range","SSErr")

# print variogram

  preds = variogramLine(lzn.fit, maxdist = max(lzn.vgm$dist))
  vgm.plot<-ggplot(lzn.vgm,aes(x=dist,y=gamma))+geom_point() + geom_line(data=preds)
  ggsave(vgm.plot,filename = paste0("Vgm_CoV2_genome_cmt_ls_",dt,".png"),dpi=100)


#IR_bottom<-signif(range(IR)[1],3)-0.01
#IR_top<-signif(range(IR)[2],3)+0.01

# expand grid to nucleotide level on x (29903 nc), 100 pts on y (total 3M)
# nc level not required for visualization / movie. Go for lower res

#krige_grid<-expand.grid(seqX=seq(0,1,by=1/29903)[2:29904],IR=seq(0,1,by=0.01))  # 3M grid for compatibility with first pass
krige_grid<-expand.grid(seqX=seq(0,1,by=0.001),IR=seq(0,1,by=0.01)) 


coordinates(krige_grid)=~seqX+IR
gridded(krige_grid)=TRUE

if(nrow(zerodist(krige_tab))!=0){
  krige_tab<-krige_tab[-zerodist(krige_tab)[,1],]
}

lzn.kriged = krige(FR~1, krige_tab,krige_grid, model = lzn.fit,nmin=5,nmax=30)

# c("red","orange","yellow","green") brewer.pal(n=11,name="Spectral")

out<-krige.cv(FR~1,krige_tab,lzn.fit,nmin=5,nmax=30)

out_v<-c(dt,mean(out$residual),mean(out$residual^2),mean(out$zscore^2),cor(out$observed, out$observed - out$residual),cor(out$observed - out$residual, out$residual))
names(out_v)<-c("Item","ME","MSPE","MSNE","cor_obs_pred","cor_pred_res")

out_v=c(out_v,vgm_par)


print(round(as.numeric(out_v[5]),3)) # nmax=20: 0.423 /nmax=30:  0.434 /nmax=100:  0.443

lzn.kriged_df<-as.data.frame(lzn.kriged)
#lzn.kriged_df$seqX=sapply(lzn.kriged_df$seqX,function(x) round(x,9))
write.table(lzn.kriged_df,quote=F,row.names=F,sep="\t",file =paste0("kriged_df_genome_cmt_ls_",dt,".tsv"))

# plot the genome map

#lowlim<-ifelse(min(lzn.kriged_df$var1.pred)<0,min(lzn.kriged_df$var1.pred),0)
#highlim<-ifelse(max(lzn.kriged_df$var1.pred)>1,max(lzn.kriged_df$var1.pred),1)

cols_out<-brewer.pal(n=11,name="Spectral")


pred.plot <- ggplot(aes(x = seqX, y = IR), data = lzn.kriged_df)+xlim(0,1)+ylim(0,1)#+theme(legend.key.size = unit(0.2, "cm"))
pred.plot <- pred.plot + geom_tile(aes(fill = var1.pred))

pred.plot_simple <- pred.plot + scale_fill_gradientn(colors=rev(cols_out),limits=c(0,1),guide=F)

pred.plot <- pred.plot + scale_fill_gradientn(colors=rev(cols_out),limits=c(0,1)) + labs(fill = "FR") + ggtitle(dt)

#pred.plot<- pred.plot+geom_point(data=as.data.frame(krige_tab),size=0.05)

pred.plot_simple<- pred.plot_simple+theme_void()+geom_point(data=as.data.frame(krige_tab),size=0.05)

# dot sizes as AF

dots_tab<-as.data.frame(cbind(seqX,IR,FR,genetab$AF))
colnames(dots_tab)[4]<-"AF"

pred.plot<- pred.plot+geom_point(data=dots_tab,aes(x=seqX,y=IR,size=AF))+scale_size(range = c(0.01,2))




contour=F
  
  #var_sim<-runif(nrow(lzn.kriged_df), min = 0, max = 1)
  
  if(contour==T){
    # contours drawn at specific quantile thresholds of computed variance
    #  breaks<-sapply(c(0.01,0.05,0.1,0.25),function(x) quantile(var_sim,x))
    breaks<-sapply(c(0.01,0.05,0.1,0.25),function(x) quantile(lzn.kriged_df$var1.var,x))
    pred.plot_contour=pred.plot+geom_contour(aes(z=var1.var),breaks=breaks,colour="black",size=0.1)
  }
  
  # extract top 10 countries from dataset

countries_df<-dict_to_df(genetab$counted_countries)
countries_vec<-apply(countries_df,2,function(x) sum(as.numeric(x),na.rm=T))
countries_vec=countries_vec[order(countries_vec,decreasing = T)]

# arrange plot and table

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.7)),
  colhead = list(fg_params=list(cex = 0.7)),
  rowhead = list(fg_params=list(cex = 0.7)))

tbl <- tableGrob(t(countries_vec[1:10]),cols=names(countries_vec[1:10]),theme=mytheme)

p3<-grid.arrange(pred.plot,tbl, nrow = 2,heights=c(6,1))

  ggsave(p3,filename = paste0("CoV2_genome_IR_FR_cmt_ls_",dt,".png"),width = 7.34, height = 5.5)
  #ggsave(pred.plot_contour,filename = paste0("CoV2_contour_genome_IR_FR_",dt,".png"),width = 7.34, height = 4.5)
  ggsave(pred.plot_simple,filename = paste0("raw_CoV2_genome_IR_FR_cmt_ls_",dt,".png"),width = 7.02, height = 6)
  
  
  out_v_tab<-rbind(out_v_tab,out_v)
  
}

write.table(out_v_tab,quote=F,row.names=F,sep="\t",file =paste0("params_tab_genome_cmt_ls_",dt,".tsv"))
