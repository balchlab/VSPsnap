# get a genome-wide map, grid on x at nucleotide resolution
# modify seqX so values are relative to nucleotide positions
# slice the map for single proteins
# extract z values with minimal variance for each aa, to be painted on structures


# kriging map for Infectivity / fatality rate vs genome position - genome wide

library(ggplot2)
library(xlsx)
library(sp)
library(gstat)
library(grid)
library(RColorBrewer)

options(stringsAsFactors=F)

# load IR/FR data
setwd("/Users/sal/projects/BalchLab/Covid-19/2020_10_05")
genetab_all<-read.csv("../2020_09_21/unique_mutations_country9_10_20.csv")

# filter >=3 countries

genetab=genetab_all[-which(genetab_all$countries<3),]

# normalize IR and FR

genetab$IR<-(genetab$infection_rate-min(genetab$infection_rate))/(max(genetab$infection_rate)-min(genetab$infection_rate))
genetab$FR<-(genetab$fatality_rate-min(genetab$fatality_rate))/(max(genetab$fatality_rate)-min(genetab$fatality_rate))

# sequence
# Genome length (Wuhan-1): 1-29903

#seqX<-(genetab$X0-min(genetab$X0)+1)/(max(genetab$X0)-min(genetab$X0)+1)
seqX<-genetab$X0/29903

# get IR axis 
IR<-genetab$IR

# FR axis
FR<-genetab$FR

krige_tab<-as.data.frame(cbind(seqX,IR,FR))
krige_tab=krige_tab[order(krige_tab$seqX),]

if(outputTab==T){
  write.table(krige_tab,quote=F,row.names=F,sep="\t",file = "cov2_pos_IR_FR_filt_3countries_fullgenome.tsv")
}

rownames(krige_tab)<-seq(1:nrow(krige_tab))

# max euclidean distance between xy points - to compute width and cutoff
maxd<-max(dist(krige_tab[,1:2])) # 1.34

coordinates(krige_tab)=~seqX+IR

# generate variogram
#lzn.vgm<-variogram(FR~1,krige_tab,cutoff=0.1, width=0.001)#width=maxd/50,cutoff=maxd/3)
#vm <- vgm(psill=0.02928, model="Sph",range=0.006, nugget=0)

lzn.vgm<-variogram(FR~1,krige_tab,cutoff=0.2, width=0.001)
vm <- vgm(psill=0.02882, model="Exp",range=0.003, nugget=0.00462)


# debug for "zero distance semivariance" issue

if(any(lzn.vgm$dist==0)){
  lzn.vgm=lzn.vgm[-which(lzn.vgm$dist==0),]
}

# fit model to best vgm among three model options
#lzn.fit = fit.variogram(lzn.vgm, model = vgm("Sph"))
lzn.fit = fit.variogram(lzn.vgm, model =vm)

# get vgm parameters to be saved
# SSErr is a (weighted) sum of squared errors of the fitted model
vgm_par<-c(as.character(lzn.fit$model[2]),lzn.fit$psill[1],lzn.fit$psill[2],lzn.fit$range[2],attr(lzn.fit,"SSErr"))
names(vgm_par)<-c("vgm_model","Nug","psill","range","SSErr")

# print variogram option
if(vg==T){
  preds = variogramLine(lzn.fit, maxdist = max(lzn.vgm$dist))
  vgm.plot<-ggplot(lzn.vgm,aes(x=dist,y=gamma))+geom_point() + geom_line(data=preds)
  ggsave(vgm.plot,filename = paste0("CoV2_IR_FR_10_6_vgm.png"),dpi=100)
}

#IR_bottom<-signif(range(IR)[1],3)-0.01
#IR_top<-signif(range(IR)[2],3)+0.01

# expand grid to nucleotide level on x (29903 nc), 100 pts on y (total 3M)

krige_grid<-expand.grid(seqX=seq(0,1,by=1/29903)[2:29904],IR=seq(0,1,by=0.001))


coordinates(krige_grid)=~seqX+IR
gridded(krige_grid)=TRUE

if(nrow(zerodist(krige_tab))!=0){
  krige_tab<-krige_tab[-zerodist(krige_tab)[,1],]
}

lzn.kriged = krige(FR~1, krige_tab,krige_grid, model = lzn.fit,nmin=5,nmax=30)

# c("red","orange","yellow","green") brewer.pal(n=11,name="Spectral")

out<-krige.cv(FR~1,krige_tab,lzn.fit,nmin=5,nmax=30)

out_v<-c("Covid-19",mean(out$residual),mean(out$residual^2),mean(out$zscore^2),cor(out$observed, out$observed - out$residual),cor(out$observed - out$residual, out$residual))
names(out_v)<-c("Item","ME","MSPE","MSNE","cor_obs_pred","cor_pred_res")

out_v=c(out_v,vgm_par)

print(round(as.numeric(out_v[5]),3)) # nmax=20: 0.423 /nmax=30:  0.434 /nmax=100:  0.443

lzn.kriged_df<-as.data.frame(lzn.kriged)
lzn.kriged_df$seqX=sapply(lzn.kriged_df$seqX,function(x) round(x,9))
save(lzn.kriged_df,file="lzn.kriged_df_30M_10_8_20.RData")

# plot the 30M genome map

cols_out<-brewer.pal(n=11,name="Spectral")

pred.plot <- ggplot(aes(x = seqX, y = IR), data = lzn.kriged_df)#+theme(legend.key.size = unit(0.2, "cm"))
pred.plot <- pred.plot + geom_tile(aes(fill = var1.pred))
pred.plot <- pred.plot + scale_fill_gradientn(colors=rev(cols_out),limits=c(0,1))

pred.plot<- pred.plot+geom_point(data=as.data.frame(krige_tab),size=0.05)


cov_proteins<-read.xlsx2("unipCov2Chain_UCSC_Wuh1_structures.xlsx",sheetIndex = 1)
cov_proteins[,2:3]<-apply(cov_proteins[,2:3],2,as.numeric)

# assemble minvar for all proteins

minVar_all<-vector()
for(i in 1:24){
  minVar<-read.csv(paste0(cov_proteins$name[i],"_minVarZ_30M_rgb_pymol.csv"))
  minVar_all=rbind(minVar_all,minVar)
}
minVar_all=minVar_all[-7174,]
# add minVar dots and line on genome-wide map

pred.plot<-pred.plot+geom_point(data=minVar_all,aes(x = seqX, y = IR,colour=var1.var),size=0.05) + scale_colour_gradient(low="black", high="grey")
pred.plot<- pred.plot+geom_line(data=minVar_all,linetype="dotted",aes(x = seqX, y = IR))

ggsave(pred.plot,filename = "CoV2_genome_IR_FR_30M_MinVar.png",width = 24.34, height = 12.5)

# get colors - same as in the kriging map - remember to reverse!

cols_out<-brewer.pal(n=11,name="Spectral")
pal <- colorRamp(rev(cols_out))

# plot slices corrisponding to the proteins

krige_df<-as.data.frame(krige_tab)
i=24
start=cov_proteins$chromStart[i]
stop<-cov_proteins$chromEnd[i]
kriged_df<-lzn.kriged_df
L<-cov_proteins$L[i]
prot<-cov_proteins$name[i]

getMapSlice<-function(start,stop,prot,L,genomeLen,kriged_df,krige_df,muts=F,minVar=T){
  
  print(prot)
  
  L=as.numeric(L)
  
  startN<-round(start/genomeLen,9)
  stopN<-round(stop/genomeLen,9)
  
  aa<-round((seq(start+1,stop,by=3)/genomeLen)[1:L],9)
  
  df_prot<-subset(kriged_df,seqX>=startN&seqX<=stopN)
  krige_prot<-subset(krige_df,seqX>=startN&seqX<=stopN)
  
  pred.plot <- ggplot(aes(x = seqX, y = IR), data = df_prot) #+theme(legend.key.size = unit(0.2, "cm"))
  pred.plot <- pred.plot + geom_tile(aes(fill = var1.pred))
  pred.plot <- pred.plot + scale_fill_gradientn(colors=rev(cols_out),limits=c(0,1))

  if(muts==T){
   pred.plot<- pred.plot+geom_point(data=krige_prot,size=0.05)
  }
  
  if(minVar==T){
    minVar<-read.csv(paste0(prot,"_minVarZ_30M_rgb_pymol.csv"))
    pred.plot<-pred.plot+geom_point(data=minVar,aes(x = seqX, y = IR,colour=var1.var)) + scale_colour_gradient(low="black", high="grey")
    pred.plot<- pred.plot+geom_line(data=minVar,linetype="dotted",aes(x = seqX, y = IR))
  }
  
  dv<-ifelse(L<180,10,as.numeric(L)/15)
  rnd<-round(seq(1,L,by=dv),-1)
  rnd[1]<-1
  aa_pos<-aa[rnd]
  names(aa_pos)<-rnd
  aa_pos_df<-data.frame(aa_pos)
  
  pred.plot<-pred.plot + geom_text(data = aa_pos_df, aes(x = aa_pos, label = rownames(aa_pos_df), y = -0.02),size = 3, col = 'black')
  
  pred.plot=pred.plot+ggtitle(prot)
  
  W<-8.34
  H<-4.5
  if(L>300){
    W<-L/36
    H<-L/66.7
  }
  
  ggsave(pred.plot,filename = paste0(prot,"_minVar.png"),width = W, height = H,limitsize=F)
 
}

for(i in 1:nrow(cov_proteins)){
  
  print(i)
  getMapSlice(cov_proteins$chromStart[i],cov_proteins$chromEnd[i],cov_proteins$name[i],cov_proteins$L[i],29903,lzn.kriged_df,krige_df)
  
}


#### add min uncertainty line
i=1
start<-cov_proteins$chromStart[i]
stop<-cov_proteins$chromEnd[i]
prot<-cov_proteins$name[i]
L<-cov_proteins$L[i]


startN<-round(start/genomeLen,9)
stopN<-round(stop/genomeLen,9)

aa<-round((seq(start+1,stop,by=3)/genomeLen)[1:L],9)

df_prot<-subset(lzn.kriged_df,seqX>=startN&seqX<=stopN)
krige_prot<-subset(krige_df,seqX>=startN&seqX<=stopN)

pred.plot <- ggplot(aes(x = seqX, y = IR), data = df_prot) #+theme(legend.key.size = unit(0.2, "cm"))
pred.plot <- pred.plot + geom_tile(aes(fill = var1.pred))
pred.plot <- pred.plot + scale_fill_gradientn(colors=rev(cols_out),limits=c(0,1))

minVar<-read.csv("nsp1_minVarZ_30M_rgb_pymol.csv")

pred.plot<-pred.plot+geom_line(data=minVar,aes(x = seqX, y = IR,colour=var1.var)) + scale_colour_gradient(low="black", high="grey")

#pred.plot<-pred.plot+geom_point(data=minVar,aes(x = seqX, y = IR,colour=var1.var)) + scale_colour_gradient(low="black", high="grey")
pred.plot<-pred.plot+geom_point(data=minVar,aes(x = seqX, y = IR,colour=var1.var)) + scale_colour_gradient(low="black", high="grey")

pred.plot+geom_line(data=minVar,aes(x = seqX, y = IR,colour=var1.var))

pred.plot<- pred.plot+geom_line(data=minVar,linetype="dotted",aes(x = seqX, y = IR))

pred.plot<- pred.plot+geom_point(data=krige_prot,size=0.05)

###################



# get interval corresponding to selected protein

# get prediction value for each residue with minimal variance
getMinVar<-function(input,df_prot,pal){
  block<-df_prot[which(df_prot[,1]==input),]
  res<-block[which.min(block[,4]),]
  res=cbind(res,pal(res[1,3])/255)
  colnames(res)[5:7]<-c("R","G","B")
  return(res)
}

# given a protein, get a table of MinVarZ - one per aa 

getMinVarZ<-function(start,stop,prot,L,genomeLen,kriged_df,writeOut=T){
 
 print(prot)
  
 startN<-round(start/genomeLen,9)
 stopN<-round(stop/genomeLen,9)

 df_prot<-subset(kriged_df,seqX>=startN&seqX<=stopN)

 aa<-round((seq(start+1,stop,by=3)/genomeLen)[1:L],9)

 out<-do.call("rbind",lapply(aa,function(x) getMinVar(x,df_prot,pal))) 
 out=cbind(1:nrow(out),out)
 colnames(out)[1]<-"residue"

 if(writeOut==T){
   write.csv(out,file=paste0(prot,"_minVarZ_30M_rgb_fixed_pymol.csv"),quote=F,row.names=F)
 }
 return(out)
}

a<-getMinVarZ(cov_proteins$chromStart[1],cov_proteins$chromEnd[1],cov_proteins$name[1],cov_proteins$L[1],29903,lzn.kriged_df)

# loop on all proteins available in cov_proteins

for(i in 1:nrow(cov_proteins)){
  
  print(i)
  getMinVarZ(cov_proteins$chromStart[i],cov_proteins$chromEnd[i],cov_proteins$name[i],cov_proteins$L[i],29903,lzn.kriged_df)
  
}


# get tables for single proteins - add mutation type and normalized  genomic position

# get aa mutations
genetab_aa<-read.csv("unique_mutations_country9_10_20_aa.csv")
 
load("../2020_09_28/CoV2_genetab_4663mut_VEP.Rdata")

length(which(genetab$descriptor%in%genetab_aa$descriptor)) # 4663, all

genetab_aa=genetab_aa[!duplicated(genetab_aa$descriptor),]

genetab_annot<-merge(genetab,genetab_aa[,c(1,7:10)],by="descriptor",all.x=T)

genetab_annot=genetab_annot[order(genetab_annot$X0),]

genetab_annot$counted_countries=gsub("\\{|\\}|\\'","",genetab_annot$counted_countries)

genetab_annot=genetab_annot[,-15]

genetab_annot$pos<-genetab_annot$X0/29903

genetab_annot=genetab_annot[,c(1:3,6:10,15:18,11:12,19,13,14)]

colnames(genetab_annot)[1:3]<-c("Mutation","Start","Stop")
  
  
for(i in 1:nrow(cov_proteins)){
  
 genetab_s<-subset(genetab_annot,Start>=cov_proteins$chromStart[i]&Stop<=cov_proteins$chromEnd[i])

 aa_pos<-read.csv(paste0(cov_proteins$name[i],"_minVarZ_30M_rgb_pymol.csv"))

 genetab_s$pos_AA<-sapply(genetab_s$pos,function(x) which.min(abs(x-aa_pos$seqX)))
 
 genetab_s=genetab_s[,c(1:9,18,10,11,13:17)]
 
write.xlsx2(genetab_s,sheetName="all_mut",row.names=F,file=paste0(cov_proteins$name[i],"_mut_table.xlsx"))
write.xlsx2(genetab_s[-which(genetab_s$FR<0.7),],row.names=F,sheetName="FR_07",append=T,file=paste0(cov_proteins$name[i],"_mut_table.xlsx"))

}




