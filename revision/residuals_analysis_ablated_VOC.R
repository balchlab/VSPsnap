# focused analysis on VoC from GPR maps - ablation experiments for Cell Patterns review (4/2023)

options(stringsAsFactors = F)

library(xlsx)
library(ggplot2)
library(Rmisc)

### Load functions

# function to locate closest element on the grid
getKrigedCoords<-function(inp,kriged){
  #  print(inp[1])
  inp[1]=as.numeric(inp[1])
  inp[1]
  inp[2]=as.numeric(inp[2])
  x_block<-kriged[which(abs(kriged$seqX-as.numeric(inp[1]))==min(abs(kriged$seqX-as.numeric(inp[1])))),]
  out<-x_block[which(abs(x_block$IR-as.numeric(inp[2]))==min(abs(x_block$IR-as.numeric(inp[2])))),]
  
  return(out)
}

# function to compute residuals and plot them
getResiduals<-function(dots_tab,genetab,VOC,kriged,lab,rnd=F,rndsize=10,plot=F){
  
  rownames(VOC)<-VOC$descriptor
  
  if(rnd==F){
    dots_tab_VOC<-dots_tab[which(genetab$descriptor%in%VOC$descriptor),]
    dots_tab_VOC$name<-VOC[genetab$descriptor[which(genetab$descriptor%in%VOC$descriptor)],10]
  }else{
    rnd_ind<-sample(1:dim(dots_tab)[1],rndsize)
    dots_tab_VOC<-dots_tab[rnd_ind,]
    dots_tab_VOC$name<-genetab$AA[rnd_ind]
    lab=paste0(lab,sample(0:1000,1))
    
  }
  
  pred<-do.call("rbind",apply(dots_tab_VOC,1,function(x) getKrigedCoords(x,kriged)))
  
  out_p<-cbind(dots_tab_VOC[,c(5,1:4)],pred[,3:4])
  
  #res_p<-out_p$var1.pred-out_p$FR
  res_p<-out_p$FR-out_p$var1.pred
  names(res_p)<-out_p$name
  
  out_p=cbind(out_p,res_p)
  
  if(lab=="B117"){
    y_lim<-c(-0.3,1.2)
  }else if(lab=="P1"){
    y_lim<-c(-1.2,1.1)
  }else if(lab=="B1351"){
    y_lim<-c(-1,1)
  }else if(lab=="B1617"){
    y_lim<-c(-1,0.5)
  }
  
  if(plot==T){
    png(paste0(lab,"_resid2_",dt,".png"),width=850,height=500)
    plt<-barplot(res_p,xaxt="n",main=paste0(lab," - ",dt),cex.main=1.3,cex.axis =1.2,ylim=y_lim)
    text(plt, par("usr")[3], labels = names(res_p), srt = 60, adj = c(1.1,1.1), xpd = TRUE, cex=1.3)
    dev.off()
  }
  
  return(out_p)
}

# function to plot mean residuals

plotMeanRes<-function(res_l,dt_vec,VOC){
  res_avg_alpha<-unlist(lapply(res_l,function(x) mean(x$res_p)))
  res_sd_alpha<-unlist(lapply(res_l,function(x) sd(x$res_p)))
  res_CI_up_alpha<-unlist(lapply(res_l,function(x) CI(x$res_p)[1]))
  res_CI_low_alpha<-unlist(lapply(res_l,function(x) CI(x$res_p)[3]))
  
  res_tab_alpha<-cbind(dt_vec,res_avg_alpha,res_CI_up_alpha,res_CI_low_alpha)
  colnames(res_tab_alpha)<-c("date","mean_res","CI_up","CI_low")
  
  res_tab_alpha=as.data.frame(res_tab_alpha)
  res_tab_alpha[,2:4]=apply(res_tab_alpha[,2:4],2,as.numeric)
  
  res_tab_alpha$date <- factor(res_tab_alpha$date, levels = res_tab_alpha$date)
  
  mean_res_plot<-ggplot(res_tab_alpha, aes(date, mean_res)) +
    geom_line(color="blue",group=1) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_up),group=1,alpha=0.2)+
    geom_hline(yintercept=0, linetype='dashed', col = 'darkred')+
    theme(axis.text.x = element_text(angle = 45, hjust=1,size=14,color=c("black","transparent","transparent","transparent")),axis.text.y=element_text(size=14))
  
  ggsave(mean_res_plot,filename = sprintf("mean_res_CI_%s.png",VOC))
}

###

# 1. read in VoC mutation data 

setwd("/Users/sal/projects/BalchLab/Covid-19/2023_04_27-revision/")

# load dates
inFiles<-read.delim("../VSPsnap/data/dates_5_31_21.tsv",header=F)$V1

# weekly sampling starting Sept 20

dates<-seq(214,486,by=7)

# load mutations for VOCs

B117<-read.xlsx("../VSPsnap/data/cov2_strains_mutations_5_25_21.xlsx",sheetName = "B117")
P1<-read.xlsx("../VSPsnap/data/cov2_strains_mutations_5_25_21.xlsx",sheetName = "P1")
B1351<-read.xlsx("../VSPsnap/data/cov2_strains_mutations_5_25_21.xlsx",sheetName = "B1351")
#CAL20<-read.xlsx("../VSPsnap/data/cov2_strains_mutations.xlsx",sheetName = "CAL.20C")
B1617<-read.xlsx("../VSPsnap/data/cov2_strains_mutations_5_25_21.xlsx",sheetName = "B1617")

cov_proteins<-read.xlsx2("../VSPsnap/data/unipCov2Chain_UCSC_Wuh1_edited_CNBC_3.xlsx",sheetIndex = 1)
cov_proteins[,2:3]<-apply(cov_proteins[,2:3],2,as.numeric)


#getResiduals(dots_tab,genetab,B117,kriged,"rnd",rnd=T)

ind=0
res_l_B117<-list()
res_l_P1<-list()
res_l_B1351<-list()
res_l_B1617<-list()

# I. ablated alpha

for(i in dates){
  
  ind=ind+1
  
  dt<-strsplit(inFiles[i],"[.]")[[1]][1]
  print(dt)
  
  genetab_all<-read.csv(paste0("../VSPsnap_current/data/cumulative_daily_AF_v4_lag_5_31/",inFiles[i]))
  genetab_all$AF<-genetab_all$counts/sum(genetab_all$counts)
  
  # filter < 3 counts
  genetab_all=genetab_all[-which(genetab_all$counts<3),]
  
  # filter >=3 countries
  genetab=genetab_all[-which(genetab_all$countries<3),]
  
  # remove FR=0 if any
  
  if(length(which(genetab$fatality_rate==0))>0){
    genetab=genetab[-which(genetab$fatality_rate==0),]
  }
  
  # log transform IR and FR (no 0-1 scaling)
  genetab$IR<-log10(genetab$infection_rate)
  genetab$FR<-log10(genetab$fatality_rate)
  
  
  if(length(which("27972_27972_REF=C_ALT=T"%in%genetab$descriptor))>0){
    genetab$AA[which(genetab$descriptor%in%"27972_27972_REF=C_ALT=T")]<-"Q27stop"
  }
  
  # sequence
  # Genome length (Wuhan-1): 1-29903
  
  seqX<-genetab$Start/29903
  
  # get IR axis 
  IR<-genetab$IR
  
  # FR axis
  FR<-genetab$FR
  
  dots_tab<-as.data.frame(cbind(seqX,IR,FR,genetab$AF))
  colnames(dots_tab)[4]<-"AF"
  
  # get predicted - alpha ablation
  kriged<-read.delim(paste0("alpha/kriged_df_maxd075_log_bestModel_",dt,".tsv"),sep="\t")


  res_l_B117[[ind]]<-getResiduals(dots_tab,genetab,B117,kriged,"B117")
  #res_l_P1[[ind]]<-getResiduals(dots_tab,genetab,P1,kriged,"P1")
  #res_l_B1351[[ind]]<-getResiduals(dots_tab,genetab,B1351,kriged,"B1351")
  #res_l_B1617[[ind]]<-getResiduals(dots_tab,genetab,B1617,kriged,"B1617")

}

dt_vec<-vector()
for(i in dates) dt_vec=c(dt_vec,strsplit(inFiles[i],"[.]")[[1]][1])


plotMeanRes(res_l_B117,dt_vec,"Alpha")
#plotMeanRes(res_l_B1351,dt_vec,"Beta")
#plotMeanRes(res_l_P1,dt_vec,"Gamma")
#plotMeanRes(res_l_B1617,dt_vec,"Delta")

# I. ablated beta

ind=0

for(i in dates){
  
  ind=ind+1
  
  dt<-strsplit(inFiles[i],"[.]")[[1]][1]
  print(dt)
  
  genetab_all<-read.csv(paste0("../VSPsnap_current/data/cumulative_daily_AF_v4_lag_5_31/",inFiles[i]))
  genetab_all$AF<-genetab_all$counts/sum(genetab_all$counts)
  
  # filter < 3 counts
  genetab_all=genetab_all[-which(genetab_all$counts<3),]
  
  # filter >=3 countries
  genetab=genetab_all[-which(genetab_all$countries<3),]
  
  # remove FR=0 if any
  
  if(length(which(genetab$fatality_rate==0))>0){
    genetab=genetab[-which(genetab$fatality_rate==0),]
  }
  
  # log transform IR and FR (no 0-1 scaling)
  genetab$IR<-log10(genetab$infection_rate)
  genetab$FR<-log10(genetab$fatality_rate)
  
  
  if(length(which("27972_27972_REF=C_ALT=T"%in%genetab$descriptor))>0){
    genetab$AA[which(genetab$descriptor%in%"27972_27972_REF=C_ALT=T")]<-"Q27stop"
  }
  
  # sequence
  # Genome length (Wuhan-1): 1-29903
  
  seqX<-genetab$Start/29903
  
  # get IR axis 
  IR<-genetab$IR
  
  # FR axis
  FR<-genetab$FR
  
  dots_tab<-as.data.frame(cbind(seqX,IR,FR,genetab$AF))
  colnames(dots_tab)[4]<-"AF"
  
  # get predicted - alpha ablation
  kriged<-read.delim(paste0("beta/kriged_df_maxd075_log_bestModel_",dt,".tsv"),sep="\t")
  
  
  
  
  res_l_B1351[[ind]]<-getResiduals(dots_tab,genetab,B1351,kriged,"B1351")
 
  
}

plotMeanRes(res_l_B1351,dt_vec,"Beta")


# I. ablated gamma

ind=0

for(i in dates){
  
  ind=ind+1
  
  dt<-strsplit(inFiles[i],"[.]")[[1]][1]
  print(dt)
  
  genetab_all<-read.csv(paste0("../VSPsnap_current/data/cumulative_daily_AF_v4_lag_5_31/",inFiles[i]))
  genetab_all$AF<-genetab_all$counts/sum(genetab_all$counts)
  
  # filter < 3 counts
  genetab_all=genetab_all[-which(genetab_all$counts<3),]
  
  # filter >=3 countries
  genetab=genetab_all[-which(genetab_all$countries<3),]
  
  # remove FR=0 if any
  
  if(length(which(genetab$fatality_rate==0))>0){
    genetab=genetab[-which(genetab$fatality_rate==0),]
  }
  
  # log transform IR and FR (no 0-1 scaling)
  genetab$IR<-log10(genetab$infection_rate)
  genetab$FR<-log10(genetab$fatality_rate)
  
  
  if(length(which("27972_27972_REF=C_ALT=T"%in%genetab$descriptor))>0){
    genetab$AA[which(genetab$descriptor%in%"27972_27972_REF=C_ALT=T")]<-"Q27stop"
  }
  
  # sequence
  # Genome length (Wuhan-1): 1-29903
  
  seqX<-genetab$Start/29903
  
  # get IR axis 
  IR<-genetab$IR
  
  # FR axis
  FR<-genetab$FR
  
  dots_tab<-as.data.frame(cbind(seqX,IR,FR,genetab$AF))
  colnames(dots_tab)[4]<-"AF"
  
  # get predicted - alpha ablation
  kriged<-read.delim(paste0("gamma/kriged_df_maxd075_log_bestModel_",dt,".tsv"),sep="\t")
  
  
  
  
   res_l_P1[[ind]]<-getResiduals(dots_tab,genetab,P1,kriged,"P1")
  
  
}

plotMeanRes(res_l_P1,dt_vec,"Gamma")

# I. ablated delta

ind=0

for(i in dates){
  
  ind=ind+1
  
  dt<-strsplit(inFiles[i],"[.]")[[1]][1]
  print(dt)
  
  genetab_all<-read.csv(paste0("../VSPsnap_current/data/cumulative_daily_AF_v4_lag_5_31/",inFiles[i]))
  genetab_all$AF<-genetab_all$counts/sum(genetab_all$counts)
  
  # filter < 3 counts
  genetab_all=genetab_all[-which(genetab_all$counts<3),]
  
  # filter >=3 countries
  genetab=genetab_all[-which(genetab_all$countries<3),]
  
  # remove FR=0 if any
  
  if(length(which(genetab$fatality_rate==0))>0){
    genetab=genetab[-which(genetab$fatality_rate==0),]
  }
  
  # log transform IR and FR (no 0-1 scaling)
  genetab$IR<-log10(genetab$infection_rate)
  genetab$FR<-log10(genetab$fatality_rate)
  
  
  if(length(which("27972_27972_REF=C_ALT=T"%in%genetab$descriptor))>0){
    genetab$AA[which(genetab$descriptor%in%"27972_27972_REF=C_ALT=T")]<-"Q27stop"
  }
  
  # sequence
  # Genome length (Wuhan-1): 1-29903
  
  seqX<-genetab$Start/29903
  
  # get IR axis 
  IR<-genetab$IR
  
  # FR axis
  FR<-genetab$FR
  
  dots_tab<-as.data.frame(cbind(seqX,IR,FR,genetab$AF))
  colnames(dots_tab)[4]<-"AF"
  
  # get predicted - alpha ablation
  kriged<-read.delim(paste0("delta/kriged_df_maxd075_log_bestModel_",dt,".tsv"),sep="\t")
  
  
  
  
  res_l_B1617[[ind]]<-getResiduals(dots_tab,genetab,B1617,kriged,"B1617")
  
  
}

plotMeanRes(res_l_B1617,dt_vec,"Delta")








# Omicron

# 1. read in VoC mutation data 

setwd("/Users/sal/projects/BalchLab/Covid-19/2022_02_11/")


# load dates
inFiles<-read.delim("../VSPsnap/data/dates_01_16_22.tsv",header=F)$V1

# weekly sampling starting Sept 20

dates<-seq(561,716,by=7)

# load mutations for VOCs

Omicron<-read.xlsx("../VSPsnap/data/omicron_8_22_21.xlsx",sheetIndex=1)
#rownames(Omicron)<-Omicron$descriptor

cov_proteins<-read.xlsx2("../VSPsnap/data/unipCov2Chain_UCSC_Wuh1_edited_CNBC_3.xlsx",sheetIndex = 1)
cov_proteins[,2:3]<-apply(cov_proteins[,2:3],2,as.numeric)

ind=0
res_l_Omicron<-list()


for(i in dates){
  
  ind=ind+1
  
  dt<-strsplit(inFiles[i],"[.]")[[1]][1]
  print(dt)
  
  genetab_all<-read.csv(paste0("../VSPsnap/data/cumulative_daily_AF_v4_lag_081421_011622/",inFiles[i]))
  genetab_all$AF<-genetab_all$counts/sum(genetab_all$counts)
  
  # filter < 3 counts
  genetab_all=genetab_all[-which(genetab_all$counts<3),]
  
  # filter >=3 countries
  genetab=genetab_all[-which(genetab_all$countries<3),]
  
  # remove FR=0 if any
  
  if(length(which(genetab$fatality_rate==0))>0){
    genetab=genetab[-which(genetab$fatality_rate==0),]
  }
  
  # log transform IR and FR (no 0-1 scaling)
  genetab$IR<-log10(genetab$infection_rate)
  genetab$FR<-log10(genetab$fatality_rate)
  
  
  if(length(which("27972_27972_REF=C_ALT=T"%in%genetab$descriptor))>0){
    genetab$AA[which(genetab$descriptor%in%"27972_27972_REF=C_ALT=T")]<-"Q27stop"
  }
  
  # sequence
  # Genome length (Wuhan-1): 1-29903
  
  seqX<-genetab$Start/29903
  
  # get IR axis 
  IR<-genetab$IR
  
  # FR axis
  FR<-genetab$FR
  
  dots_tab<-as.data.frame(cbind(seqX,IR,FR,genetab$AF))
  colnames(dots_tab)[4]<-"AF"
  
  # get predicted
  kriged<-read.delim(paste0("../VSPsnap/data/krige_output/kriged_df_maxd075_log_bestModel_",dt,".tsv"),sep="\t")
  
  
  res_l_Omicron[[ind]]<-getResiduals(dots_tab,genetab,Omicron,kriged,"Omicron")
  
}

dt_vec<-vector()
for(i in dates) dt_vec=c(dt_vec,strsplit(inFiles[i],"[.]")[[1]][1])

plotMeanRes(res_l_Omicron,dt_vec,"Omicron")

