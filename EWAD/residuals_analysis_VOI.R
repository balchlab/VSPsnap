# focused analysis on VoC from GPR maps

library(xlsx)
library(ggplot2)

# 1. read in Variant mutation data 

setwd("/Users/sal/projects/BalchLab/Covid-19/VSPsnap/autoEWAD")

options(stringsAsFactors = F)

# load dates
inFiles<-read.delim("../data/dates_5_31_21.tsv",header=F)$V1

# weekly sampling starting Sept 20

dates<-seq(214,486,by=7)

# load mutations for VOI

mu<-read.xlsx("../data/VOI.xlsx",sheetName = "mu")
lambda<-read.xlsx("../data/VOI.xlsx",sheetName = "lambda")

cov_proteins<-read.xlsx2("../data/unipCov2Chain_UCSC_Wuh1_edited_CNBC_3.xlsx",sheetIndex = 1)
cov_proteins[,2:3]<-apply(cov_proteins[,2:3],2,as.numeric)

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
getResiduals<-function(dots_tab,genetab,VOC,kriged,lab,rnd=F,rndsize=10){
    
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
    
    res_p<-out_p$var1.pred-out_p$FR
    names(res_p)<-out_p$name
    
    out_p=cbind(out_p,res_p)
    
#    if(lab=="B117"){
    y_lim<-c(-1,1)
#    }else if(lab=="P1"){
#      y_lim<-c(-1.2,1.1)
#    }else if(lab=="B1351"){
#      y_lim<-c(-1,1)
#    }else if(lab=="B1617"){
#      y_lim<-c(-1,0.5)
#    }
    
    png(paste0(lab,"_resid_",dt,".png"),width=850,height=500)
    plt<-barplot(res_p,xaxt="n",main=paste0(lab," - ",dt),cex.main=1.3,cex.axis =1.2,ylim=y_lim)
    text(plt, par("usr")[3], labels = names(res_p), srt = 60, adj = c(1.1,1.1), xpd = TRUE, cex=1.3)
    dev.off()
    
    return(out_p)
  }

#getResiduals(dots_tab,genetab,B117,kriged,"rnd",rnd=T)

ind=0
res_l_mu<-list()
res_l_lambda<-list()




for(i in dates){
  
  ind=ind+1
  
  dt<-strsplit(inFiles[i],"[.]")[[1]][1]
  print(dt)
  
  genetab_all<-read.csv(paste0("../data/cumulative_daily_AF_v4_lag_5_31/",inFiles[i]))
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
  kriged<-read.delim(paste0("../data/krige_output/kriged_df_maxd075_log_bestModel_",dt,".tsv"),sep="\t")


  res_l_mu[[ind]]<-getResiduals(dots_tab,genetab,mu,kriged,"Mu")
  res_l_lambda[[ind]]<-getResiduals(dots_tab,genetab,lambda,kriged,"Lambda")
  
}


# random selection of mutations 

VOCrnd1<-genetab[sample(1:dim(genetab)[1],10),]
VOCrnd2<-genetab[sample(1:dim(genetab)[1],10),]
VOCrnd3<-genetab[sample(1:dim(genetab)[1],10),]
VOCrnd4<-genetab[sample(1:dim(genetab)[1],15),]
VOCrnd5<-genetab[sample(1:dim(genetab)[1],15),]


ind=0
res_l_rnd1<-list()
res_l_rnd2<-list()
res_l_rnd3<-list()
res_l_rnd4<-list()
res_l_rnd5<-list()


for(i in dates){
  
  ind=ind+1
  
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
  
  
  res_l_rnd1[[ind]]<-getResiduals(dots_tab,genetab,VOCrnd1,kriged,"rnd1")
  res_l_rnd2[[ind]]<-getResiduals(dots_tab,genetab,VOCrnd2,kriged,"rnd2")
  res_l_rnd3[[ind]]<-getResiduals(dots_tab,genetab,VOCrnd3,kriged,"rnd3")
  res_l_rnd4[[ind]]<-getResiduals(dots_tab,genetab,VOCrnd4,kriged,"rnd4")
  res_l_rnd5[[ind]]<-getResiduals(dots_tab,genetab,VOCrnd5,kriged,"rnd5")
  
}  
  
  













































