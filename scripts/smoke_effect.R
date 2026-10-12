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

data_sim <- pheno_baseline[,c("study_participant_id_1","smoke_diag","smoking_history","long_quit_smok")]

dataCD <- merge(dataCD,data_sim,by.x="SP0011_id_1",by.y="study_participant_id_1")


#effect of number of cigarettes on Crohn's disease progression
smy_part1<-summary(glm(as.double(behaviour!=1) ~ as.factor(smoke_diag) , data=dataCD[dataCD$smoking_history==2 ,]  ,family = "binomial"))

smy_part1coef<- as.data.frame(smy_part1$coefficients[2:3,])
colnames(smy_part1coef)<- c("beta","se","z","pval")

#effect of time to quit smoking on Crohn's disease progression
smy_part2<-summary(glm(as.double(behaviour!=1) ~ as.factor(long_quit_smok),data=dataCD[dataCD$smoking_history==3,],family = "binomial"))

smy_part2coef<- as.data.frame(smy_part2$coefficients[2:3,])
colnames(smy_part2coef)<- c("beta","se","z","pval")

smytotal <- rbind(smy_part1coef,smy_part2coef)
smytotal$or <- round(exp(smytotal$beta),2)
smytotal$ci <- paste0(round(exp(smytotal$beta - smytotal$se*1.96),2),"-",round(exp(smytotal$beta + smytotal$se*1.96),2))
smytotal$orci <- paste0(smytotal$or,"(",smytotal$ci,")")
smytotal$variable <- c("5+ v.s. 5- cigarettes per day","pipe only v.s. 5- cigarettes per day","Quit: 1-6 months v.s. < 1 month","Quit: > 6 months v.s. < 1 month")

write.csv(smytotal,file="result/smoking_detailed_effect.csv",quote=F,row.names=F)


