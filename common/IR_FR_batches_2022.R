# generate dates


setwd("../2020_11_11/")

inFiles<-read.delim("dates.tsv",header=F)$V1

dt<-strsplit(strsplit(inFiles[i],"window21_")[[1]][2],"[.]")[[1]][1]

a<-read.csv(paste0("window_21/",inFiles[i]))

# determine min and max for IR, FR across all time tables

options(stringsAsFactors = F)

inFiles<-read.delim("dates.tsv",header=F)$V1

minIR<-vector()
maxIR<-vector()
minFR<-vector()
maxFR<-vector()

allIR<-vector()
allFR<-vector()
allIR_3<-vector()
allFR_3<-vector()

for(i in 1:length(inFiles)){
  
  a<-read.csv(paste0("window_21/",inFiles[i]))
  a=a[-which(a$countries<3),]
#  minIR=c(minIR,min(a$infection_rate))
#  maxIR=c(maxIR,max(a$infection_rate))
  
#  minFR=c(minFR,min(a$fatality_rate))
#  maxFR=c(maxFR,max(a$fatality_rate)

allIR_3<-c(allIR_3,a$infection_rate)
allFR_3<-c(allFR_3,a$fatality_rate)
    
  
  gc()
}

# max(maxFR) 217.7066
# min(minFR) 0
# min(minIR) 7.292497e-05   
# max(maxIR) 40.90636

# repeat for cumulative

allIR_3<-vector()
allFR_3<-vector()

for(i in 1:length(inFiles)){
  
  a<-read.csv(paste0("cumulative_daily/",inFiles[i]))
  a=a[-which(a$countries<3),]
  #  minIR=c(minIR,min(a$infection_rate))
  #  maxIR=c(maxIR,max(a$infection_rate))
  
  #  minFR=c(minFR,min(a$fatality_rate))
  #  maxFR=c(maxFR,max(a$fatality_rate)
  
  allIR_3<-c(allIR_3,a$infection_rate)
  allFR_3<-c(allFR_3,a$fatality_rate)
  
  
  gc()
}

quantile(allIR_3,0.99) #62235.49
quantile(allFR_3,0.99) #422.4026

# summary(allFR_3)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000    7.989   14.871   41.061   31.297 5529.640
# summary(allIR_3)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0     708.8    1570.6    5265.0    3708.0 1244444.7


# meta spawner

library("readtext")

batches<-read.table("batches.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnap",i,".R"))
  
}

# assemble the spawner line

spawner_cmd<-paste(c("./foo_spawner.sh","run_Rscript.sh",dir(pattern = "VSPsnap[^(_)]")),sep="",collapse=" ") # regexp is to exclude the VSPsnap_template.R based on the underscore.


# no 0-1 scaling on IR and FR
tmp<-readtext("VSPsnap_IR_FR_noscaled_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnapNS",i,".R"))
  
}

spawner_cmd_NS<-paste(c("./foo_spawner.sh","run_Rscript.sh",dir(pattern = "VSPsnapNS[^(_)]")),sep="",collapse=" ")


# cumulative

options(stringsAsFactors = F)
library("readtext")

batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_cmt_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnapCMT",i,".R"))
  
}

spawner_cmd_CMT<-paste(c("./foo_spawner.sh","run_Rscript.sh",dir(pattern = "VSPsnapCMT[^(_)]")),sep="",collapse=" ")

# cumulative - local scaling

options(stringsAsFactors = F)
library("readtext")

batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_cmt_ls_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnapCMT_LS",i,".R"))
  
}

spawner_cmd_CMT_LS<-paste(c("./foo_spawner.sh","run_Rscript.sh",dir(pattern = "VSPsnapCMT_LS[^(_)]")),sep="",collapse=" ")

# cumulative - slicer

batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_cmt_ls_slicer_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnapCMT_LS_slicer",i,".R"))
  
}

spawner_cmd_CMT_LS_slicer<-paste(c("./foo_spawner.sh","run_Rscript.sh",dir(pattern = "VSPsnapCMT_LS_slicer[^(_)]")),sep="",collapse=" ")

# cumulative - slicer minVar

setwd("/Users/sal/projects/BalchLab/Covid-19/2020_12_04")
batches<-read.table("../2020_11_11/batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_cmtls_minvar_slicer_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnap_minVar_slicer",i,".R"))
  
}

spawner_cmd_minvar_slicer<-paste(c("./foo_spawner.sh","run_Rscript.sh",dir(pattern = "VSPsnap_minVar_slicer[^(_)]")),sep="",collapse=" ")

# cumulative - local scaling - new scheduler (slurm) and updated dataset

library("readtext")

setwd("/Users/sal/projects/BalchLab/Covid-19/2021_1_6/")
batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_cmt_ls_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnapCMT_LS",i,".R"))
  
}

spawner_cmd_CMT_LS<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnapCMT_LS[^(_)]")),sep="",collapse=" ")

# cumulative - local scaling - new scheduler (slurm) and updated datase - nodots

library("readtext")

setwd("/Users/sal/projects/BalchLab/Covid-19/2021_1_6/")
batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_cmt_ls_nodots_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnapCMT_LS_nodots",i,".R"))
  
}

spawner_cmd_CMT_LS_nodots<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnapCMT_LS_nodots[^(_)]")),sep="",collapse=" ")


# cumulative - slicer - new dataset - slurm

batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_cmt_ls_slicer_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnapCMT_LS_slicer",i,".R"))
  
}

spawner_cmd_CMT_LS_slicer_slurm<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnapCMT_LS_slicer[^(_)]")),sep="",collapse=" ")


# cumulative local scaling + B117

batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_cmt_ls_B117_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnapCMT_LS_B117_",i,".R"))
  
}

spawner_cmd_CMT_LS_B117<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnapCMT_LS_B117_[^(_)]")),sep="",collapse=" ")

# cumulative local scaling + B117 only

batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_cmt_ls_B117only_template.R")

for(i in 1:nrow(batches)){
  
  out<-sprintf(tmp$text,batches[i,1],batches[i,2])
  cat(out,file=paste0("VSPsnapCMT_LS_B117only_",i,".R"))
  
}

spawner_cmd_CMT_LS_B117only<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnapCMT_LS_B117only_[^(_)]")),sep="",collapse=" ")


# cumulative local scaling / new dataset, norm v4 lag - log transformed
# sprintf accepts max 8192 bytes array lengths

setwd("/Users/sal/projects/BalchLab/Covid-19/2021_03_31/IR_FR_cmt_ls_4_21")

batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_cmt_ls_log_template2.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_cmt_ls_log_",i,".R"))
  
}

spawner_cmd_cmt_ls_log<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_cmt_ls_log_[^(_)]")),sep="",collapse=" ")


# bestmodel

setwd("/Users/sal/projects/BalchLab/Covid-19/2021_04_22/timelapse_log_bestModel/")

batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_log_bestModel_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_IR_FR_log_bestModel_",i,".R"))
  
}

spawner_cmd_log_bestModel<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_IR_FR_log_bestModel_[^(_)]")),sep="",collapse=" ")

# b117 annot only
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_04_22/timelapse_log_bestModel_B117")

batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_log_B117_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_B117_",i,".R"))
  
}

spawner_cmd_log_bestModel<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_B117_[^(_)]")),sep="",collapse=" ")


# b1351 annot only
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_04_22/timelapse_log_bestModel_B1351")

#batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_log_B1351_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_B1351_",i,".R"))
  
}

spawner_cmd_log_bestModel<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_B1351_[^(_)]")),sep="",collapse=" ")


# P1 annot only
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_04_22/timelapse_log_bestModel_P1")

#batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_log_P1_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_P1_",i,".R"))
  
}

spawner_cmd_log_bestModel<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_P1_[^(_)]")),sep="",collapse=" ")

# CAL20 annot only
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_04_22/timelapse_log_bestModel_CAL20")

#batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_log_CAL20_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_CAL20_",i,".R"))
  
}

spawner_cmd_log_bestModel<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_CAL20_[^(_)]")),sep="",collapse=" ")

# B1617 annot only
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_04_22/timelapse_log_bestModel_B1617")

#batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("VSPsnap_IR_FR_log_B1617_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_B1617_",i,".R"))
  
}

spawner_cmd_log_bestModel<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_B1617_[^(_)]")),sep="",collapse=" ")

# same as above (all labels, single VOC)for dataset update 5/31/21.

# bestmodel

batches<-read.table("../VSPsnap/data/batchesCMT.tsv",header=F,quote="",sep=" ")

setwd("/Users/sal/projects/BalchLab/Covid-19/2021_06_15/timelapse_log_bestModel_allLabels/")


tmp<-readtext("../../VSPsnap/VSPsnap_IR_FR_log_bestModel_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_IR_FR_log_bestModel_",i,".R"))
  
}

spawner_cmd_log_bestModel<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_IR_FR_log_bestModel_[^(_)]")),sep="",collapse=" ")

# b117 annot only
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_06_15/timelapse_log_bestModel_B117/")


tmp<-readtext("../../VSPsnap/VSPsnap_IR_FR_log_B117_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_B117_",i,".R"))
  
}

spawner_cmd_log_b117<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_B117_[^(_)]")),sep="",collapse=" ")


# b1351 annot only
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_06_15/timelapse_log_bestModel_B1351")

#batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("../../VSPsnap/VSPsnap_IR_FR_log_B1351_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_B1351_",i,".R"))
  
}

spawner_cmd_log_b1351<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_B1351_[^(_)]")),sep="",collapse=" ")


# P1 annot only
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_06_15/timelapse_log_bestModel_P1")

#batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("../../VSPsnap/VSPsnap_IR_FR_log_P1_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_P1_",i,".R"))
  
}

spawner_cmd_log_P1<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_P1_[^(_)]")),sep="",collapse=" ")


# B1617 annot only
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_06_15/timelapse_log_bestModel_B1617")

#batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("../../VSPsnap/VSPsnap_IR_FR_log_B1617_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_B1617_",i,".R"))
  
}

spawner_cmd_log_B1617<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_B1617_[^(_)]")),sep="",collapse=" ")


# One script to rule them all
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_06_21/timelapse_log_bestModel_one/")

#batches<-read.table("batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("../../VSPsnap/VSPsnap_IR_FR_log_bestModel_one_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_one_",i,".R"))
  
}

spawner_cmd_log_one<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_one_[^(_)]")),sep="",collapse=" ")


# modify the above only annotating alpha and delta
setwd("/Users/sal/projects/BalchLab/Covid-19/2021_08_04/")

batches<-read.table("/Users/sal/projects/BalchLab/Covid-19/VSPsnap/data/batchesCMT.tsv",header=F,quote="",sep=" ")

tmp<-readtext("../VSPsnap/VSPsnap_IR_FR_log_bestModel_alphaDelta_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_alphaDelta_",i,".R"))
  
}

spawner_cmd_log_alphaDelta<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_alphaDelta_[^(_)]")),sep="",collapse=" ")

# update 8/22/2021

library("readtext")

setwd("/Users/sal/projects/BalchLab/Covid-19/2021_08_31/")

batches<-read.table("/Users/sal/projects/BalchLab/Covid-19/VSPsnap/data/batchesCMT_update_8_22_21.tsv",header=F,quote="",sep=" ")

tmp<-readtext("../VSPsnap/VSPsnap_IR_FR_log_bestModel_one_update_8_22_21_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_one_",i,".R"))
  
}

spawner_cmd_log_one<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_one_[^(_)]")),sep="",collapse=" ")



# update 12/15/2021

library("readtext")

setwd("/Users/sal/projects/BalchLab/Covid-19/update_12_15_21")

batches<-read.table("/Users/sal/projects/BalchLab/Covid-19/VSPsnap/data/batchesCMT_update_12_15_21.tsv",header=F,quote="",sep=" ")

tmp<-readtext("../VSPsnap/VSPsnap_IR_FR_log_bestModel_one_update_12_15_21_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_one_121521_",i,".R"))
  
}

spawner_cmd_log_one<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_one_121521_[^(_)]")),sep="",collapse=" ")



# update 01/16/2022

library("readtext")

setwd("/Users/sal/projects/BalchLab/Covid-19/update_01_16_22/")

batches<-read.table("/Users/sal/projects/BalchLab/Covid-19/VSPsnap/data/batchesCMT_update_01_16_22.tsv",header=F,quote="",sep=" ")

tmp<-readtext("../VSPsnap/VSPsnap_IR_FR_log_bestModel_one_update_01_16_22_template.R")

tmp0<-substr(tmp$text,1,5000)
tmp1<-substr(tmp$text,5001,nchar(tmp$text))

for(i in 1:nrow(batches)){
  
  out0<-sprintf(tmp0,batches[i,1],batches[i,2])
  
  out<-paste0(out0,tmp1)
  
  cat(out,file=paste0("VSPsnap_ir_fr_log_one_011622_",i,".R"))
  
}

spawner_cmd_log_one<-paste(c("./foo_spawner_slurm.sh","run_Rscript_slurm.sh",dir(pattern = "VSPsnap_ir_fr_log_one_011622_[^(_)]")),sep="",collapse=" ")






################

# assemble animations for maps and painted structures
cd /Users/sal/projects/BalchLab/Covid-19/2020_11_23/colored_prots/nsp7
convert -antialias -delay 50 *.png nsp7_IR_FR_cumulative_ls_temporal_seq_fast.gif
# same for other proteins

# nsp4
resid<-c(4,5,7,9,11,16,17,18,19,20,21,24,28,37,38,39,40,41,50,51,62,70,90,91,94,99,102,110,113,143,163,175,180,181,183,184,218,231,235,246,250,261,264,265,267,272,275,280,284,290,293,297,303,312,317,319,328,334,345,353,363,366,367,370,372,379,380,382,383,387,400,403,408,433,434,458,483)
nsp4<-read.csv("../2020_10_26/nsp4_minVarZ_3M_rgb_pymol.csv")
