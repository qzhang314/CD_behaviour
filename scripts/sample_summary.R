library(ggplot2)
library(reshape2)
library(DiagrammeR)
library(ggsignif)
library(kableExtra)
library(survival)
library(knitr)
library(scales)
library(gridExtra)
library(grid)
library(ggpubr)
library(forestplot)
library(MASS)

source("/nfs/users/nfs_q/qz2/R_functions/func.R")

wkdir<- "~/CD_behaviour"
setwd(wkdir)

dataCD <- read.table("~/CD_behaviour/pheno/ibdbr_pheno.txt",header=T)

pheno_baseline <- readRDS("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/PhenoQC/2025-09-29/ibdbr_baseline.rds")

data_sim <- pheno_baseline[,c("study_participant_id_1","smoke_diag","smoking_history","long_quit_smok")]

dataCD <- merge(dataCD,data_sim,by.x="SP0011_id_1",by.y="study_participant_id_1")


drugTimeSinceDiagAll <- readRDS("result/drugTimeSinceDiagAll.rds")

#summary table show the number of patients in each category
demoTable<-c()

for(i in 1:3){
  dataCDpart<- dataCD[dataCD$behaviour==i, ]
#  drugworkCDpart<- drugworkAll[dataCD$behaviour==i,]
  drugTimeSinceDiagAllpart<- drugTimeSinceDiagAll[dataCD$behaviour==i,]
  for(j in 1:4){
    drugTimeSinceDiagAllpart[,j] <- as.double(drugTimeSinceDiagAllpart[,j]<2)
  }
  tmp<-c(median(dataCDpart$follow_up,na.rm=T),
         median(dataCDpart$age_diag,na.rm=T),sum(dataCDpart$age_diag_section=="A1",na.rm=T),sum(dataCDpart$age_diag_section=="A2",na.rm=T),sum(dataCDpart$age_diag_section=="A3",na.rm=T),
         sum(dataCDpart$sex==2,na.rm=T),sum(dataCDpart$sex==1,na.rm=T),
         sum(dataCDpart$smoking_history==1,na.rm=T),
         sum(dataCDpart$smoking_history==2,na.rm=T),
         sum(dataCDpart$smoking_history==2 & dataCDpart$smoke_diag==1,na.rm=T),sum(dataCDpart$smoking_history==2 & dataCDpart$smoke_diag==2,na.rm=T),sum(dataCDpart$smoking_history==2 & dataCDpart$smoke_diag==4,na.rm=T),
         sum(dataCDpart$smoking_history==3,na.rm=T),
         sum(dataCDpart$smoking_history==3 & dataCDpart$long_quit_smok==1,na.rm=T),sum(dataCDpart$smoking_history==3 & dataCDpart$long_quit_smok==2,na.rm=T),sum(dataCDpart$smoking_history==3 & dataCDpart$long_quit_smok==3,na.rm=T),
         sum(dataCDpart$diag_year<2002,na.rm=T),sum(dataCDpart$diag_year>=2002 &dataCDpart$diag_year<2012 ,na.rm=T),sum(dataCDpart$diag_year>=2012,na.rm=T),
         sum(dataCDpart$location=="L1",na.rm=T),sum(dataCDpart$location=="L2",na.rm=T),sum(dataCDpart$location=="L3",na.rm=T),
         sum(dataCDpart$location1=="L4",na.rm=T),sum(dataCDpart$location1=="noL4",na.rm=T),
         sum(drugTimeSinceDiagAllpart$Thiopurine==1,na.rm=T),sum(drugTimeSinceDiagAllpart$Thiopurine==0,na.rm=T),
         sum(drugTimeSinceDiagAllpart$Adalimumab==1,na.rm=T),sum(drugTimeSinceDiagAllpart$Adalimumab==0,na.rm=T),
         sum(drugTimeSinceDiagAllpart$Infliximab==1,na.rm=T),sum(drugTimeSinceDiagAllpart$Infliximab==0,na.rm=T),
         sum(drugTimeSinceDiagAllpart$Methotrexate==1,na.rm=T),sum(drugTimeSinceDiagAllpart$Methotrexate==0,na.rm=T))
  demoTable<-cbind(demoTable,round(tmp))
}

for(i in 1:ncol(demoTable)){
  demoTable[,i] <- comma(demoTable[,i])
}

rownames(demoTable)<-c("Median time since diagnosis (year)","Median Age at diagnosis (year)","A1","A2","A3","Female","Male","Non-Smoker","Smoker","Less than 5 cigarettes per day","More than 5 cigarettes per day","Pipe only","Ex-Smoker","Quit smoke less than 1 month","Quit smoke 1-6 months","Quit smoke more than 6 months","< 2002",">= 2002 & < 2012",">= 2012","L1","L2","L3","UGI involvement (Yes)","UGI involvement (No)","Thioprine early prescription","Thioprine late prescription","Adalimumab early prescription","Adalimumab late prescription","Infliximab early prescription","Infliximab late prescription","Methotrexate early prescription","Methotrexate late prescription")

colnames(demoTable)<- c(paste0("B1 (N = ",sum(dataCD$behaviour==1),")"),paste0("B2 (N = ",sum(dataCD$behaviour==2),")"),paste0("B3 (N = ",sum(dataCD$behaviour==3),")"))

for(i in 1:ncol(demoTable)){
  demoTable[,i] <- comma(demoTable[,i])
}

write.csv(demoTable,file="result/demoTable.csv")




