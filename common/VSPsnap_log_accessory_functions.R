
getMapSlice<-function(start,stop,prot,L,genomeLen,kriged_df,krige_df,genetab,AF_tab,dt,muts=F,minVar=F){
  
  print(prot)
  
  L=as.numeric(L)
  
  startN<-round(start/genomeLen,9)
  stopN<-round(stop/genomeLen,9)
  
  aa<-round((seq(start+1,stop,by=3)/genomeLen)[1:L],9)
  
  df_prot<-subset(kriged_df,seqX>=startN&seqX<=stopN)
  krige_prot<-subset(krige_df,seqX>=startN&seqX<=stopN)
  genetab_prot<-subset(genetab,Start>=start&End<=stop)
  AF_prot<-subset(AF_tab,seqX>=startN&seqX<=stopN)

  pred.plot <- ggplot(aes(x = seqX, y = IR), data = df_prot) #+theme(legend.key.size = unit(0.2, "cm"))
  pred.plot <- pred.plot + geom_tile(aes(fill = var1.pred))
  pred.plot <- pred.plot + scale_fill_gradientn(colors=rev(cols_out),limits=c(min(kriged_df$var1.pred),max(kriged_df$var1.pred))) + labs(fill = "FR") + ggtitle(paste0(prot,"-",dt))
  
  if(muts==T){
  #  pred.plot<- pred.plot+geom_point(data=krige_prot,size=0.05)
    pred.plot<- pred.plot+geom_point(data=AF_prot,aes(x=seqX,y=IR,size=AF))+scale_size(range = c(0.01,2))

  }
  
  if(minVar==T){
    minVar<-read.csv(paste0(prot,"_minVarZ_3M_rgb_pymol_",dt,".csv"))
    pred.plot<-pred.plot+geom_point(data=minVar,aes(x = seqX, y = IR,colour=var1.var)) + scale_colour_gradient(low="black", high="grey",guide=F)
    pred.plot<- pred.plot+geom_line(data=minVar,linetype="dotted",aes(x = seqX, y = IR))
  }
  
  dv<-ifelse(L<180,10,as.numeric(L)/15)
  rnd<-round(seq(1,L,by=dv),-1)
  rnd[1]<-1
  aa_pos<-aa[rnd]
  names(aa_pos)<-rnd
  aa_pos_df<-data.frame(aa_pos)
  
  #pred.plot<-pred.plot + geom_text(data = aa_pos_df, aes(x = aa_pos, label = rownames(aa_pos_df), y = -0.02),size = 3, col = 'black')
  
  W<-8.34
  H<-7
 # if(L>300){
 #   W<-L/36
 #   H<-L/30
 # }
  
# extract top 10 countries from dataset

if(nrow(genetab_prot)>0){
countries_df<-dict_to_df(genetab_prot$counted_countries)
countries_vec<-apply(countries_df,2,function(x) sum(as.numeric(x),na.rm=T))
countries_vec=countries_vec[order(countries_vec,decreasing = T)]

# arrange plot and table

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.7)),
  colhead = list(fg_params=list(cex = 0.7)),
  rowhead = list(fg_params=list(cex = 0.7)))

tbl <- tableGrob(t(countries_vec[1:10]),cols=names(countries_vec[1:10]),theme=mytheme)

p3<-grid.arrange(pred.plot,tbl, nrow = 2,heights=c(6,1))


  ggsave(p3,filename = paste0(prot,"_",dt,".png"),width = W, height = H,limitsize=F)
 }else{
  ggsave(pred.plot,filename = paste0(prot,"_",dt,".png"),width = W, height = H,limitsize=F)
 }

}


# get prediction value for each residue with minimal variance
getMinVar<-function(input,df_prot,pal){
  block<-df_prot[which(df_prot[,1]==input),]
  res<-block[which.min(block[,4]),]
  res=cbind(res,pal(res[1,3])/255)
  colnames(res)[5:7]<-c("R","G","B")
  return(res)
}

# given a protein, get a table of MinVarZ - one per aa 

getMinVarZ<-function(start,stop,prot,L,genomeLen,kriged_df,dt,writeOut=T){
  
  
  startN<-round(start/genomeLen,9)
  stopN<-round(stop/genomeLen,9)
  
  df_prot<-subset(kriged_df,seqX>=startN&seqX<=stopN)
  
  aa<-round((seq(start+1,stop,by=3)/genomeLen)[1:L],9)
  
  out<-do.call("rbind",lapply(aa,function(x) getMinVar(x,df_prot,pal))) 
  out=cbind(1:nrow(out),out)
  colnames(out)[1]<-"residue"
  
  if(writeOut==T){
    write.csv(out,file=paste0(prot,"_minVarZ_3M_rgb_pymol_",dt,".csv"),quote=F,row.names=F)
  }
  return(prot)
}