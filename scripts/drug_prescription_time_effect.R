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

pheno_baseline <- readRDS("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/PhenoQC/2025-09-29/ibdbr_baseline.rds")

dataCD <- read.table("~/CD_behaviour/pheno/ibdbr_pheno.txt",header=T)

dataCD <- dataCD[!is.na(dataCD$timeSinceDiag_5yrs) & !is.na(dataCD$age_diag_section) & !is.na(dataCD$smoker_type) & !is.na(dataCD$year_diag_section) &  !is.na(dataCD$location) &  !is.na(dataCD$location1) &  !is.na(dataCD$sex),]

startIDpart<- c("inflix_start_date","ada_start_date","aza_start_date","metho_start_date")
drugworkIDpart<- c("inflix_work","ada_work","aza_work","metho_work")
drugnamespart<-c("Infliximab","Adalimumab","Thiopurine","Methotrexate")

data_sim <- pheno_baseline[,c("study_participant_id_1",startIDpart,drugworkIDpart,"merc_start_date","merc_work","clinic_review_date")]

dataCD <- merge(dataCD,data_sim,by.x="SP0011_id_1",by.y="study_participant_id_1")


# effect of drug prescription time on Crohn's disease progression (focused on patients with drug taken)
# definition of early drug prescription
# time since diagnosis, focusing on patients with drug taken only: early diagnosis: drug start date - year of diagnosis < 2, late diagnosis: drug start date - year of diagnosis >= 2

drugTimeSinceDiagAll<-c()
drugworkAll<-c()

for(i in startIDpart){
  dataCD[,i] <- sapply(dataCD[,i],getyear_new)
}

dataCD$merc_start_date <- sapply(dataCD$merc_start_date,getyear_new)
dataCD$clinic_review_date <- sapply(dataCD$clinic_review_date,getyear_new)

for(j in 1:length(startIDpart)){
  drugTimeSinceDiag<-c()
  #merge Mercaptopurine and Azathioprine into Thioprine
  if(j==3){
    Thioprine_startID<-c("aza_start_date","merc_start_date")
    Thioprine_workID<- c("aza_work","merc_work")
    for(i in 1:nrow(dataCD)){
      k<-1
      if(sum(is.na(dataCD[i,Thioprine_startID]))==2){
        drugTimeSinceDiag<-c(drugTimeSinceDiag,NA)
        # drugworktmp<-c(drugworktmp,NA)
        next
      }else if(sum(is.na(dataCD[i,Thioprine_startID]))==1){
        if(is.na(dataCD[i,"aza_start_date"])){
          k<-2
        }
      }else{
        if(dataCD[i,"aza_start_date"]>dataCD[i,"merc_start_date"]){
          k<-2
        }
      }
      
      if((dataCD[i,Thioprine_startID[k]] > dataCD[i,"clinic_review_date"]) &!is.na(dataCD[i,"clinic_review_date"])){
        drugTimeSinceDiag<-c(drugTimeSinceDiag,NA)
      }else{
        drugTime2Diag <- trunc(dataCD[i,Thioprine_startID[k]]) - dataCD[i,"diag_year"]
        drugTime2Diag[!is.na(drugTime2Diag) & (drugTime2Diag<0)]<- NA
        if(sum(!is.na(drugTime2Diag))==0){
          drugTimeSinceDiag<-c(drugTimeSinceDiag,NA)
        }else{
          drugTimeSinceDiag<-c(drugTimeSinceDiag,drugTime2Diag)
        }
      }
    }
  }else{
    for(i in 1:nrow(dataCD)){
      if(is.na(dataCD[i,startIDpart[j]])){
        drugTimeSinceDiag<-c(drugTimeSinceDiag,NA)
      }else if((dataCD[i,startIDpart[j]] > dataCD[i,"clinic_review_date"]) & !is.na(dataCD[i,"clinic_review_date"])){
        drugTimeSinceDiag<-c(drugTimeSinceDiag,NA)
      }
      else{
        drugTime2Diag <- trunc(dataCD[i,startIDpart[j]]) - dataCD[i,"diag_year"]
        drugTime2Diag[!is.na(drugTime2Diag) & (drugTime2Diag<0)]<- NA
        if(sum(!is.na(drugTime2Diag))==0){
          drugTimeSinceDiag<-c(drugTimeSinceDiag,NA)
        }else{
          drugTimeSinceDiag<-c(drugTimeSinceDiag,drugTime2Diag)
        }
      }
    }
  }
  
  dataCD$drugTimeSinceDiag<- drugTimeSinceDiag
  drugTimeSinceDiagAll <- cbind(drugTimeSinceDiagAll,drugTimeSinceDiag)
}

drugTimeSinceDiagAll<- as.data.frame(drugTimeSinceDiagAll)
colnames(drugTimeSinceDiagAll) <- drugnamespart

smy_drug<-c()
for(i in 1:4){
  smy<-summary(glm(bev ~ as.double(drugTimeSinceDiagAll[,i]<=2) + timeSinceDiag_5yrs  + as.factor(smoker_type) + year_diag_section + location + location1,data=dataCD,family = "binomial"))
  smy_drug<- rbind(smy_drug,smy$coefficients[2,])
}


smy_drug<- as.data.frame(smy_drug)
colnames(smy_drug)<-c("beta","se","z","pval")
smy_drug$drugname <- drugnamespart
smy_drug$variable<-"Early prescription v.s. Late prescription"

smy_drug$or <- round(exp(smy_drug$beta),2)
smy_drug$ci <- paste0(round(exp(smy_drug$beta - smy_drug$se*1.96),2),"-",round(exp(smy_drug$beta + smy_drug$se*1.96),2))

smy_drug$type <- c("Infliximab","Adalimumab","Thiopurine","Methotrexate")

write.csv(smy_drug,file="result/drug_presb_effect.csv",row.names=F,quote=F)

saveRDS(drugTimeSinceDiagAll,file="result/drugTimeSinceDiagAll.rds")




# explore the change of prescription time across eras
drugTimeSinceDiagAllnew <- drugTimeSinceDiagAll

for(i in 1:4){
  drugTimeSinceDiagAllnew[,i] <- as.double(drugTimeSinceDiagAll[,i]>=3)+1
}
dataCD23<- cbind(dataCD,drugTimeSinceDiagAllnew)
dataCD23_2002<- dataCD23[dataCD23$year_diag_section=="pre-2002",]
dataCD23_2012<- dataCD23[dataCD23$year_diag_section=="2002-2011",]
dataCD23_2020<- dataCD23[dataCD23$year_diag_section=="post-2011",]

dataplot<-data.frame(value=c(sum(dataCD23_2002$Adalimumab==1,na.rm=T),sum(dataCD23_2002$Adalimumab==2,na.rm=T),sum(dataCD23_2002$Infliximab==1,na.rm=T),sum(dataCD23_2002$Infliximab==2,na.rm=T),sum(dataCD23_2002$Thiopurine==1,na.rm=T),sum(dataCD23_2002$Thiopurine==2,na.rm=T),sum(dataCD23_2002$Methotrexate==1,na.rm=T),sum(dataCD23_2002$Methotrexate==2,na.rm=T),sum(dataCD23_2012$Adalimumab==1,na.rm=T),sum(dataCD23_2012$Adalimumab==2,na.rm=T),sum(dataCD23_2012$Infliximab==1,na.rm=T),sum(dataCD23_2012$Infliximab==2,na.rm=T),sum(dataCD23_2012$Thiopurine==1,na.rm=T),sum(dataCD23_2012$Thiopurine==2,na.rm=T),sum(dataCD23_2012$Methotrexate==1,na.rm=T),sum(dataCD23_2012$Methotrexate==2,na.rm=T),sum(dataCD23_2020$Adalimumab==1,na.rm=T),sum(dataCD23_2020$Adalimumab==2,na.rm=T),sum(dataCD23_2020$Infliximab==1,na.rm=T),sum(dataCD23_2020$Infliximab==2,na.rm=T),sum(dataCD23_2020$Thiopurine==1,na.rm=T),sum(dataCD23_2020$Thiopurine==2,na.rm=T),sum(dataCD23_2020$Methotrexate==1,na.rm=T),sum(dataCD23_2020$Methotrexate==2,na.rm=T)),drug=c(rep(drugnamespart,each=2),rep(drugnamespart,each=2),rep(drugnamespart,each=2)),type=rep(c("Early prescription","Late prescription"),12),dateDig=rep(c("pre-2002","2002-2011","post-2011"),each=8))

dat_text<- data.frame(
  drug=rep(drugnamespart,3),
  y=c(sum(dataCD23_2002$Adalimumab==1,na.rm=T)+sum(dataCD23_2002$Adalimumab==2,na.rm=T),sum(dataCD23_2002$Infliximab==1,na.rm=T)+sum(dataCD23_2002$Infliximab==2,na.rm=T),sum(dataCD23_2002$Thiopurine==1,na.rm=T)+sum(dataCD23_2002$Thiopurine==2,na.rm=T),sum(dataCD23_2002$Methotrexate==1,na.rm=T)+sum(dataCD23_2002$Methotrexate==2,na.rm=T),sum(dataCD23_2012$Adalimumab==1,na.rm=T)+sum(dataCD23_2012$Adalimumab==2,na.rm=T),sum(dataCD23_2012$Infliximab==1,na.rm=T)+sum(dataCD23_2012$Infliximab==2,na.rm=T),sum(dataCD23_2012$Thiopurine==1,na.rm=T)+sum(dataCD23_2012$Thiopurine==2,na.rm=T),sum(dataCD23_2012$Methotrexate==1,na.rm=T)+sum(dataCD23_2012$Methotrexate==2,na.rm=T),sum(dataCD23_2020$Adalimumab==1,na.rm=T)+sum(dataCD23_2020$Adalimumab==2,na.rm=T),sum(dataCD23_2020$Infliximab==1,na.rm=T)+sum(dataCD23_2020$Infliximab==2,na.rm=T),sum(dataCD23_2020$Thiopurine==1,na.rm=T)+sum(dataCD23_2020$Thiopurine==2,na.rm=T),sum(dataCD23_2020$Methotrexate==1,na.rm=T)+sum(dataCD23_2020$Methotrexate==2,na.rm=T)),
  label=c(percent(sum(dataCD23_2002$Adalimumab==1,na.rm=T)/sum(!is.na(dataCD23_2002$Adalimumab))),percent(sum(dataCD23_2002$Infliximab==1,na.rm=T)/sum(!is.na(dataCD23_2002$Infliximab))),percent(sum(dataCD23_2002$Thiopurine==1,na.rm=T)/sum(!is.na(dataCD23_2002$Thiopurine))),percent(sum(dataCD23_2002$Methotrexate==1,na.rm=T)/sum(!is.na(dataCD23_2002$Methotrexate))),percent(sum(dataCD23_2012$Adalimumab==1,na.rm=T)/sum(!is.na(dataCD23_2012$Adalimumab))),percent(sum(dataCD23_2012$Infliximab==1,na.rm=T)/sum(!is.na(dataCD23_2012$Infliximab))),percent(sum(dataCD23_2012$Thiopurine==1,na.rm=T)/sum(!is.na(dataCD23_2012$Thiopurine))),percent(sum(dataCD23_2012$Methotrexate==1,na.rm=T)/sum(!is.na(dataCD23_2012$Methotrexate))),percent(sum(dataCD23_2020$Adalimumab==1,na.rm=T)/sum(!is.na(dataCD23_2020$Adalimumab))),percent(sum(dataCD23_2020$Infliximab==1,na.rm=T)/sum(!is.na(dataCD23_2020$Infliximab))),percent(sum(dataCD23_2020$Thiopurine==1,na.rm=T)/sum(!is.na(dataCD23_2020$Thiopurine))),percent(sum(dataCD23_2020$Methotrexate==1,na.rm=T)/sum(!is.na(dataCD23_2020$Methotrexate)))),
  dateDig=rep(c("pre-2002","2002-2011","post-2011"),each=4),
  type="Early prescription"
)

dataplot$dateDig<- factor(dataplot$dateDig,levels=unique(dataplot$dateDig))

p <- ggplot(dataplot, aes(x=dateDig,y=value,fill=type)) + geom_bar(col="black",alpha=0.5,position="stack", stat="identity")+ylab("Number of recipients")+xlab("Date of diagnosis")+theme_bw()+theme(axis.text= element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),strip.text = element_text(size=14,face="bold"),legend.text=element_text(face="bold",size=14),legend.title=element_text(face="bold",size=14),legend.position = "bottom",panel.grid.major.x= element_blank())+scale_fill_manual(name="Prescription",labels=c("Early (within two years of diangosis)","Late"), values=c("#F8766D","#00BFC4"))+geom_text(data = dat_text,aes(x=dateDig,y=y,label=label),fontface="bold",size=4,vjust=-0.3)+facet_wrap(~drug)

ggsave("figures/prescription_per_era.png",width = 12, height = 8,dpi=300)

