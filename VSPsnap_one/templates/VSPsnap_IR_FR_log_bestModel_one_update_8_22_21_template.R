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
library(ggrepel)

#source("dict_to_df.R")

# initialize out_tab
out_v_tab<-vector()

# load data
inFiles<-read.delim("../VSPsnap/data/dates_update_8_22_21.tsv",header=F)$V1

# load mutations for specific strains. Load CoV2 proteins.

B117<-read.xlsx("../VSPsnap/data/cov2_strains_mutations_5_25_21.xlsx",sheetName = "B117")
rownames(B117)<-B117$descriptor

P1<-read.xlsx("../VSPsnap/data/cov2_strains_mutations_5_25_21.xlsx",sheetName = "P1")
rownames(P1)<-P1$descriptor

B1351<-read.xlsx("../VSPsnap/data/cov2_strains_mutations_5_25_21.xlsx",sheetName = "B1351")
rownames(B1351)<-B1351$descriptor

B1617<-read.xlsx("../VSPsnap/data/cov2_strains_mutations_5_25_21.xlsx",sheetName = "B1617")
rownames(B1617)<-B1617$descriptor

cov_proteins<-read.xlsx2("../VSPsnap/data/unipCov2Chain_UCSC_Wuh1_edited_CNBC_3.xlsx",sheetIndex = 1)
cov_proteins[,2:3]<-apply(cov_proteins[,2:3],2,as.numeric)

for(i in %i:%i){
  print(i)

dt<-strsplit(inFiles[i],"[.]")[[1]][1]

  print(dt)


genetab_all<-read.csv(paste0("../VSPsnap/data/cumulative_daily_AF_v4_lag_8_22/",inFiles[i]))
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

# get vgm parameters to be saved
vgm_par<-c(as.character(lzn.fit$model[2]),lzn.fit$psill[1],lzn.fit$psill[2],lzn.fit$range[2],attr(lzn.fit,"SSErr"))
names(vgm_par)<-c("vgm_model","Nug","psill","range","SSErr")

# print variogram

  preds = variogramLine(lzn.fit, maxdist = max(lzn.vgm$dist))
  vgm.plot<-ggplot(lzn.vgm,aes(x=dist,y=gamma))+geom_point() + geom_line(data=preds)
  ggsave(vgm.plot,filename = paste0("Vgm_cutoff_maxd075_log_bestModel_",dt,".png"),dpi=100)


# nc level not required for visualization / movie. Go for lower res

krige_grid<-expand.grid(seqX=seq(0,1,by=0.001),IR=seq(min(genetab$IR),max(genetab$IR),by=0.01)) 


coordinates(krige_grid)=~seqX+IR
gridded(krige_grid)=TRUE

if(nrow(zerodist(krige_tab))!=0){
  krige_tab<-krige_tab[-zerodist(krige_tab)[,1],]
}

lzn.kriged = krige(FR~1, krige_tab,krige_grid, model = lzn.fit,nmin=5,nmax=50)


out<-krige.cv(FR~1,krige_tab,lzn.fit,nmin=5,nmax=50)

out_v<-c(dt,mean(out$residual),mean(out$residual^2),mean(out$zscore^2),cor(out$observed, out$observed - out$residual),cor(out$observed - out$residual, out$residual))
names(out_v)<-c("Item","ME","MSPE","MSNE","cor_obs_pred","cor_pred_res")

out_v=c(out_v,vgm_par)


print(round(as.numeric(out_v[5]),3))

lzn.kriged_df<-as.data.frame(lzn.kriged)

write.table(lzn.kriged_df,quote=F,row.names=F,sep="\t",file =paste0("kriged_df_maxd075_log_bestModel_",dt,".tsv"))

# plot the genome map

cols_out<-brewer.pal(n=11,name="Spectral")


pred.plot <- ggplot(aes(x = seqX, y = IR), data = lzn.kriged_df)+ylim(min(genetab$IR),max(genetab$IR)+0.15)

pred.plot <- pred.plot + geom_tile(aes(fill = var1.pred))

#pred.plot_simple <- pred.plot + scale_fill_gradientn(colors=rev(cols_out),limits=c(min(lzn.kriged_df$var1.pred),max(lzn.kriged_df$var1.pred)),guide=F)

pred.plot <- pred.plot + scale_fill_gradientn(colors=rev(cols_out),limits=c(min(lzn.kriged_df$var1.pred),max(lzn.kriged_df$var1.pred))) + labs(fill = "FR") + ggtitle(dt)


# pred.plot_simple<- pred.plot_simple+theme_void()#+geom_point(data=as.data.frame(krige_tab),size=0.05)

# dot sizes as AF - alpha shaded

dots_tab<-as.data.frame(cbind(seqX,IR,FR,genetab$AF))
colnames(dots_tab)[4]<-"AF"

# plot only top 10 percent AF dots (smaller png size)
dots_tab_filt=dots_tab[-which(dots_tab$AF<quantile(dots_tab$AF,0.90)),]

pred.plot<- pred.plot+geom_point(data=dots_tab_filt,aes(x=seqX,y=IR,size=AF),alpha=0.1)+scale_size(range = c(0.01,2))

contour=T

if(contour==T){
    # contours drawn at specific quantile thresholds of computed variance
    breaks4<-sapply(0.25,function(x) quantile(lzn.kriged_df$var1.var,x))
    breaks3<-sapply(0.1,function(x) quantile(lzn.kriged_df$var1.var,x))

    pred.plot_contour3=pred.plot+geom_contour(aes(z=var1.var),breaks=breaks3,colour="black",size=0.1)
    pred.plot_contour4=pred.plot+geom_contour(aes(z=var1.var),breaks=breaks4,colour="black",size=0.1)
    
  }
  
  # add vertical lines delimiting proteins
  pred.plot_contour3=pred.plot_contour3+geom_vline(xintercept = c(cov_proteins$chromStart/29903,cov_proteins$chromEnd[24]/29903), linetype="dotted",color="blue", size=0.3)
  pred.plot_contour4=pred.plot_contour4+geom_vline(xintercept = c(cov_proteins$chromStart/29903,cov_proteins$chromEnd[24]/29903), linetype="dotted",color="blue", size=0.3)
  
  # add protein names on top
  
  pred.plot_contour3=pred.plot_contour3+annotate(geom = "text", x =(cov_proteins$chromStart+(cov_proteins$chromEnd-cov_proteins$chromStart)/2)/29903,y =max(genetab$IR), label = cov_proteins$name, color = "blue",angle = 45,vjust=0,hjust=0)
  pred.plot_contour4=pred.plot_contour4+annotate(geom = "text", x =(cov_proteins$chromStart+(cov_proteins$chromEnd-cov_proteins$chromStart)/2)/29903,y =max(genetab$IR), label = cov_proteins$name, color = "blue",angle = 45,vjust=0,hjust=0)
  

  # add mutations for B117, P1,B1351,B1617

  if(length(which(match("27972_27972_REF=C_ALT=T",genetab$descriptor,nomatch=0)>0))>0){
    genetab$AA[which(match(genetab$descriptor,"27972_27972_REF=C_ALT=T",nomatch=0)>0)]<-"Q27stop"
  }
  
  
  # B117
  dots_tab_B117<-dots_tab[which(match(genetab$descriptor,B117$descriptor,nomatch=0)>0),]
  dots_tab_B117$name<-B117[genetab$descriptor[which(match(genetab$descriptor,B117$descriptor,nomatch=0)>0)],10]
  if(nrow(dots_tab_B117)>0) dots_tab_B117$cl<-"1"
  

  dots_tab_P1<-dots_tab[which(match(genetab$descriptor,P1$descriptor,nomatch=0)>0),]
  dots_tab_P1$name<-P1[genetab$descriptor[which(match(genetab$descriptor,P1$descriptor,nomatch=0)>0)],10]
  if(nrow(dots_tab_P1)>0) dots_tab_P1$cl<-"2"
  
  dots_tab_B1351<-dots_tab[which(match(genetab$descriptor,B1351$descriptor,nomatch=0)>0),]
  dots_tab_B1351$name<-B1351[genetab$descriptor[which(match(genetab$descriptor,B1351$descriptor,nomatch=0)>0)],10] 
  if(nrow(dots_tab_B1351)>0) dots_tab_B1351$cl<-"3"
  
  dots_tab_B1617<-dots_tab[which(match(genetab$descriptor,B1617$descriptor,nomatch=0)>0),]
  dots_tab_B1617$name<-B1617[genetab$descriptor[which(match(genetab$descriptor,B1617$descriptor,nomatch=0)>0)],10] 
  if(nrow(dots_tab_B1617)>0) dots_tab_B1617$cl<-"4"
  
  # fixed order of legends (AF top, FR bottom)
  pred.plot_contour3=pred.plot_contour3+guides(fill = guide_colourbar(order = 2),size = guide_legend(order = 1,override.aes=list(color="black",alpha=1)))  # overrides previous; add override.aes here
  pred.plot_contour4=pred.plot_contour4+guides(fill = guide_colourbar(order = 2),size = guide_legend(order = 1,override.aes=list(color="black",alpha=1)))  # overrides previous; add override.aes here



  # all annotations
  dots_tab_combo<-rbind(dots_tab_B117,dots_tab_P1,dots_tab_B1351,dots_tab_B1617)
  dots_tab_combo=dots_tab_combo[!duplicated(dots_tab_combo$name),]


  if(nrow(dots_tab_combo)>0){
  colsL <- c("1" = "black", "2" = "green4", "3" = "blue3", "4" = "brown4")
  pred.plot_contour3_all<-pred.plot_contour3+geom_text_repel(data=dots_tab_combo,aes(x=seqX,y=IR,label=name,colour=cl),size=4,max.overlaps=Inf)+scale_colour_manual(values=colsL)
  pred.plot_contour3_all=pred.plot_contour3_all+geom_point(data=dots_tab_combo,aes(x=seqX,y=IR,colour=cl),size=1,inherit.aes=F)+guides(colour=FALSE)
  
  pred.plot_contour4_all<-pred.plot_contour4+geom_text_repel(data=dots_tab_combo,aes(x=seqX,y=IR,label=name,colour=cl),size=4,max.overlaps=Inf)+scale_colour_manual(values=colsL)
  pred.plot_contour4_all=pred.plot_contour4_all+geom_point(data=dots_tab_combo,aes(x=seqX,y=IR,colour=cl),size=1,inherit.aes=F)+guides(colour=FALSE)
  }else{
    pred.plot_contour3_all<-pred.plot_contour3
    pred.plot_contour4_all<-pred.plot_contour4
  }
  
  ggsave(pred.plot_contour3_all,filename = sprintf("ctTen_all_%03d.png",i),width = 10, height = 7)
  ggsave(pred.plot_contour4_all,filename = sprintf("ctTwoFive_all_%03d.png",i),width = 10, height = 7)

  
 # B117
  if(nrow(dots_tab_B117)>0){
  colsL <- c("1" = "black", "2" = "green4", "3" = "blue3", "4" = "brown4")
  pred.plot_contour3_B117<-pred.plot_contour3+geom_text_repel(data=dots_tab_B117,aes(x=seqX,y=IR,label=name,colour=cl),size=4,max.overlaps=Inf)+scale_colour_manual(values=colsL)
  pred.plot_contour3_B117=pred.plot_contour3_B117+geom_point(data=dots_tab_B117,aes(x=seqX,y=IR,colour=cl),size=1,inherit.aes=F)+guides(colour=FALSE)
  
  pred.plot_contour4_B117<-pred.plot_contour4+geom_text_repel(data=dots_tab_B117,aes(x=seqX,y=IR,label=name,colour=cl),size=4,max.overlaps=Inf)+scale_colour_manual(values=colsL)
  pred.plot_contour4_B117=pred.plot_contour4_B117+geom_point(data=dots_tab_B117,aes(x=seqX,y=IR,colour=cl),size=1,inherit.aes=F)+guides(colour=FALSE)
  }else{
    pred.plot_contour3_B117<-pred.plot_contour3
    pred.plot_contour4_B117<-pred.plot_contour4
  }
  
  ggsave(pred.plot_contour3_B117,filename = sprintf("ctTen_B117_%03d.png",i),width = 10, height = 7)
  ggsave(pred.plot_contour4_B117,filename = sprintf("ctTwoFive_B117_%03d.png",i),width = 10, height = 7)


# P1
  if(nrow(dots_tab_P1)>0){
  colsL <- c("1" = "black", "2" = "green4", "3" = "blue3", "4" = "brown4")
  pred.plot_contour3_P1<-pred.plot_contour3+geom_text_repel(data=dots_tab_P1,aes(x=seqX,y=IR,label=name,colour=cl),size=4,max.overlaps=Inf)+scale_colour_manual(values=colsL)
  pred.plot_contour3_P1=pred.plot_contour3_P1+geom_point(data=dots_tab_P1,aes(x=seqX,y=IR,colour=cl),size=1,inherit.aes=F)+guides(colour=FALSE)
  
  pred.plot_contour4_P1<-pred.plot_contour4+geom_text_repel(data=dots_tab_P1,aes(x=seqX,y=IR,label=name,colour=cl),size=4,max.overlaps=Inf)+scale_colour_manual(values=colsL)
  pred.plot_contour4_P1=pred.plot_contour4_P1+geom_point(data=dots_tab_P1,aes(x=seqX,y=IR,colour=cl),size=1,inherit.aes=F)+guides(colour=FALSE)
  }else{
    pred.plot_contour3_P1<-pred.plot_contour3
    pred.plot_contour4_P1<-pred.plot_contour4
  }
  
  ggsave(pred.plot_contour3_P1,filename = sprintf("ctTen_P1_%03d.png",i),width = 10, height = 7)
  ggsave(pred.plot_contour4_P1,filename = sprintf("ctTwoFive_P1_%03d.png",i),width = 10, height = 7)

  
  # B1351
  if(nrow(dots_tab_B1351)>0){
  colsL <- c("1" = "black", "2" = "green4", "3" = "blue3", "4" = "brown4")
  pred.plot_contour3_B1351<-pred.plot_contour3+geom_text_repel(data=dots_tab_B1351,aes(x=seqX,y=IR,label=name,colour=cl),size=4,max.overlaps=Inf)+scale_colour_manual(values=colsL)
  pred.plot_contour3_B1351=pred.plot_contour3_B1351+geom_point(data=dots_tab_B1351,aes(x=seqX,y=IR,colour=cl),size=1,inherit.aes=F)+guides(colour=FALSE)
  
  pred.plot_contour4_B1351<-pred.plot_contour4+geom_text_repel(data=dots_tab_B1351,aes(x=seqX,y=IR,label=name,colour=cl),size=4,max.overlaps=Inf)+scale_colour_manual(values=colsL)
  pred.plot_contour4_B1351=pred.plot_contour4_B1351+geom_point(data=dots_tab_B1351,aes(x=seqX,y=IR,colour=cl),size=1,inherit.aes=F)+guides(colour=FALSE)
  }else{
    pred.plot_contour3_B1351<-pred.plot_contour3
    pred.plot_contour4_B1351<-pred.plot_contour4
  }
  
  ggsave(pred.plot_contour3_B1351,filename = sprintf("ctTen_B1351_%03d.png",i),width = 10, height = 7)
  ggsave(pred.plot_contour4_B1351,filename = sprintf("ctTwoFive_B1351_%03d.png",i),width = 10, height = 7)



  # B1617
  if(nrow(dots_tab_B1617)>0){
  colsL <- c("1" = "black", "2" = "green4", "3" = "blue3", "4" = "brown4")
  pred.plot_contour3_B1617<-pred.plot_contour3+geom_text_repel(data=dots_tab_B1617,aes(x=seqX,y=IR,label=name,colour=cl),size=4,max.overlaps=Inf)+scale_colour_manual(values=colsL)
  pred.plot_contour3_B1617=pred.plot_contour3_B1617+geom_point(data=dots_tab_B1617,aes(x=seqX,y=IR,colour=cl),size=1,inherit.aes=F)+guides(colour=FALSE)
  
  pred.plot_contour4_B1617<-pred.plot_contour4+geom_text_repel(data=dots_tab_B1617,aes(x=seqX,y=IR,label=name,colour=cl),size=4,max.overlaps=Inf)+scale_colour_manual(values=colsL)
  pred.plot_contour4_B1617=pred.plot_contour4_B1617+geom_point(data=dots_tab_B1617,aes(x=seqX,y=IR,colour=cl),size=1,inherit.aes=F)+guides(colour=FALSE)
  }else{
    pred.plot_contour3_B1617<-pred.plot_contour3
    pred.plot_contour4_B1617<-pred.plot_contour4
  }
  
  ggsave(pred.plot_contour3_B1617,filename = sprintf("ctTen_B1617_%03d.png",i),width = 10, height = 7)
  ggsave(pred.plot_contour4_B1617,filename = sprintf("ctTwoFive_B1617_%03d.png",i),width = 10, height = 7)



  
  out_v_tab<-rbind(out_v_tab,out_v)

}

write.table(out_v_tab,quote=F,row.names=F,sep="\t",file =paste0("params_tab_maxd075_log_bestModel_",dt,".tsv"))

















