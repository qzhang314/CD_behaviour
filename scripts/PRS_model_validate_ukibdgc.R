library(pROC)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggcorrplot)
library(dplyr)


wkdir<- "/nfs/users/nfs_q/qz2/CD_behaviour/"
setwd(wkdir)

predmodel <- readRDS("phenoEff/pheno_model.rds")

#Apply the model to samples from GWAS1 and GWAS3

data <- read.table("~/CD_behaviour/pheno/ukibdgc_pheno.txt",header=T)

prs1 <- read.csv("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/PRS/score/UKIBDGC_PRS_EUR.csv")
prs2 <- read.csv("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/score/UKIBDGC_PRS.csv")

identical(prs1$SPID,prs2$SPID)
# [1] TRUE

prs <- cbind(prs1,prs2$CD_bev_PRS)
colnames(prs)[ncol(prs)] <- "CD_bev_PRS"

data <- merge(data,prs,by.x="FID",by.y="SPID")

data_noNA <- data[!is.na(data$era) & !is.na(data$location) & !is.na(data$location1) & !is.na(data$smoke_ed) & !is.na(data$follow_up),]



cov <- read.table("UKIBDGC/pheno/behav_covar.regenie",header=T,stringsAsFactors = F)
prs <- merge(prs,cov,by="FID")

data13 <- prs

data13 <- data13[,c("FID","bev","timesincediag","loc","loc1","smoke","era","bev_raw","PRS")]
colnames(data13) <- c("FID","bev","timesincediag","location","location1","smoker_type","year_diag_section","bev_raw","bev_PRS")

data13[ data13$smoker_type=="non-Smoker" &!is.na(data13$smoker_type),"smoker_type"]<- "Never-smoker"
data13[ data13$smoker_type=="Smoker" &!is.na(data13$smoker_type),"smoker_type"]<- "Current-smoker"

data13$timeSinceDiag_5yrs <- data13$timesincediag/5

other_prs1 <- read.csv("UKIBDGC/score/GWAS1_PRS.csv")
other_prs3 <- read.csv("UKIBDGC/score/GWAS3_PRS.csv")

other_prs <- rbind(other_prs1,other_prs3)

data13 <- merge(data13,other_prs,by.x="FID",by.y="SPID")

data13$pred <- predict(step.model,data13,type = "response")
data13$pred_lm <- predict(step.model,data13,type = "link")

data13 <- data13[!is.na(data13$pred),]
#estimate the genetic correlation between triats

data13_pgs <- data13[,c("IBD_PRS","CD_PRS","UC_PRS","CDvsUC_PRS","bev_PRS")]
colnames(data13_pgs) <- c("IBD risk PGS","CD risk PGS","UC risk PGS","IBD subtype PGS","CD behaviour PGS")
cormat <- cor(data13_pgs)

p1 <- ggcorrplot(cormat, hc.order = TRUE, type = "lower",
           outline.col = "black",colors = c("#6D9EC1", "white", "#E46726"),
           lab = TRUE,lab_size=5)+
  theme(axis.text=element_text(size=10,face="bold"),
        legend.title = element_text(size=10,face="bold"),
        legend.text = element_text(size=10,face="bold"))


#fit risk score as covariates
smy <- c()
tmp <- summary(glm(bev ~ scale(IBD_PRS) + pred_lm,family = "binomial",data=data13))
smy <- rbind(smy,tmp$coefficients[2,])
tmp <- summary(glm(bev ~ scale(CD_PRS)+ pred_lm,family = "binomial",data=data13))
smy <- rbind(smy,tmp$coefficients[2,])
tmp <- summary(glm(bev ~ scale(UC_PRS)+ pred_lm,family = "binomial",data=data13))
smy <- rbind(smy,tmp$coefficients[2,])
tmp <- summary(glm(bev ~ scale(CDvsUC_PRS)+ pred_lm,family = "binomial",data=data13))
smy <- rbind(smy,tmp$coefficients[2,])
tmp <- summary(glm(bev ~ scale(bev_PRS)+ pred_lm,family = "binomial",data=data13))
smy <- rbind(smy,tmp$coefficients[2,])

smy <- as.data.frame(smy)
colnames(smy)<- c("est","se","zscore","p")
smy$prs <- c("IBD risk PGS","CD risk PGS","UC risk PGS","IBD subtype PGS","CD bebaviour PGS")

smy$prs <- factor(smy$prs,levels=c("IBD risk PGS","CD risk PGS","UC risk PGS","IBD subtype PGS","CD bebaviour PGS"))

p2 <- ggplot(data=smy,aes(prs,est,col=prs))+geom_point(position = position_dodge(width = 0.2),size=6)+
  geom_errorbar(aes(ymin=(est-1.96*se),ymax=(est+1.96*se)),width=0.1,position = position_dodge(width = 0.2))+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=10,face="bold"),
        axis.title.y=element_text(size=10,face="bold"),
        axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10,face="bold"),
        legend.position="top")+
  guides(col=guide_legend(nrow=2,byrow=TRUE))+
  ylab("standardized effect size (s.e.)")+
  scale_color_brewer(palette = "Set2")+
  geom_hline(yintercept = 0,lty=2)



png(file="figures/ukibdgc_prs_auc.png",units="in", width=8, height=8, res=300)
rocobj <- roc(data13$bev, data13$bev_PRS,ci=T,stratified=FALSE,plot=T,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
# Call:
#   roc.default(response = data13$bev, predictor = data13$bev_PRS,     ci = T, plot = T, stratified = FALSE, print.auc = F, show.thres = TRUE,     font = 2, font.lab = 2, cex = 1.3, col = "lightblue")
# 
# Data: data13$bev_PRS in 1613 controls (data13$bev 0) < 1177 cases (data13$bev 1).
# Area under the curve: 0.5621
# 95% CI: 0.5405-0.5836 (DeLong)

sens.ci <- ci.se(rocobj)
plot(sens.ci, type="shape", col="lightblue",font=1.5,font.lab=2)
text(0.85, 0.05, pos = 4,paste0("CDb PGS: AUC = ", round(rocobj$auc, 3),", 95% CI: ",round(rocobj$ci[1], 3),"-",round(rocobj$ci[3], 3)),cex=1.2,font=2)
dev.off()

summary(aov(lm(bev ~ bev_PRS + pred, data=data13)))
# Df Sum Sq Mean Sq F value   Pr(>F)    
# bev_PRS        1   10.1   10.15   51.88 7.58e-13 ***
#   pred           1  125.2  125.23  640.28  < 2e-16 ***
#   Residuals   2787  545.1    0.20      

10.1/(10.1 + 125.2 + 545.1)


set.seed(314)

idx <- 1:nrow(data13)
datanew <- c()
#use 10-fold cross validation to explore the benifit of adding PGS in the model

num <- 10
for(i in 1:num){
  if(i!=num){
    idx_sel <- sample(idx,round(nrow(data13)/num),replace=F)
  }else{
    idx_sel <- idx
  }
  idx <- idx[ !(idx %in% idx_sel)]
  datatrain <- data13[-idx_sel,]
  datatest <- data13[idx_sel,]
  newmodel <- glm(bev ~ bev_PRS + pred_lm, data=datatrain,family="binomial")
  datatest$totalscore <- predict(newmodel,newdata=datatest)
  datanew <- rbind(datanew,datatest)
}


png(file="figures/baseline_auc.png",units="in", width=8, height=8, res=300)
rocobj1<-roc(datanew$bev, datanew$pred,ci=T,stratified=FALSE,plot=T,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
sens.ci <- ci.se(rocobj1)
plot(sens.ci, type="shape", col="lightblue",font=1.5,font.lab=2)
text(0.85, 0.05, pos = 4,paste0("Baseline model: AUC = ", round(rocobj1$auc, 3),", 95% CI: ",round(rocobj1$ci[1], 3),"-",round(rocobj1$ci[3], 3)),cex=1.2,font=2)
dev.off()

png(file="figures/baseline_prs_auc.png",units="in", width=8, height=8, res=300)
rocobj2<-roc(datanew$bev, datanew$totalscore,ci=T,stratified=FALSE,plot=T,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
sens.ci <- ci.se(rocobj2)
plot(sens.ci, type="shape", col="lightblue",font=1.5,font.lab=2)
text(0.91, 0.05, pos = 4,paste0("Baseline model + CDb PGS: AUC = ", round(rocobj2$auc, 3),", 95% CI: ",round(rocobj2$ci[1], 3),"-",round(rocobj2$ci[3], 3)),cex=1.2,font=2)
dev.off()


roc.test(rocobj2,rocobj1,method="delong")

# DeLong's test for two correlated ROC curves
# 
# data:  rocobj2 and rocobj1
# Z = 1.0829, p-value = 0.2788
# alternative hypothesis: true difference in AUC is not equal to 0
# 95 percent confidence interval:
#  -0.001637289  0.005680479
# sample estimates:
# AUC of roc1 AUC of roc2 
#   0.7639469   0.7619253 



data13$tile1 <- ntile(data13$bev_PRS, 100)

smys <- c()
for(i in 50:99){
  tmp <- summary(glm(bev ~ as.double(tile1 > i) + pred,family="binomial",data=data13))
  smys <- rbind(smys,tmp$coefficients[2,])
}

smys <- as.data.frame(smys)
colnames(smys) <- c("beta","se","z","pval")
smys$idx <- 50:99

p3 <- ggplot(data=smys,aes(idx,exp(beta)))+geom_point()+
  scale_y_log10()+theme_bw()+
  geom_line(aes(y=exp(beta))) + 
  geom_ribbon(aes(ymin=exp(beta - 1.96*se), ymax=exp(beta + 1.96*se)),alpha=0.2)+
  xlab("Percentile of CD behaviour PGS")+
  ylab("Odds Ratio (95% CI) of\ndeveloping B2/B3")+
  theme(axis.text = element_text(size=10,face="bold"),
        axis.title = element_text(size=10,face="bold"))

p <- ggarrange(p1,ggarrange(p2,p3,nrow=2,labels = c("B","C")),labels = c("A",""),ncol=2,widths = c(0.5,0.8))

ggsave(file="figures/PRS_pred.png",width=12,height=6)



data13_noNA$tile2 <- ntile(data13_noNA$phenoscore, 100)

smys <- c()
for(i in 50:99){
  tmp <- summary(glm(bev ~ as.double(tile2 > i) + phenoscore,family="binomial",data=data13_noNA))
  smys <- rbind(smys,tmp$coefficients[3,])
}

smys <- as.data.frame(smys)
colnames(smys) <- c("beta","se","z","pval")
smys$idx <- 50:99

p4 <- ggplot(data=smys,aes(idx,exp(beta)))+geom_point()+
  scale_y_log10()+theme_bw()+
  geom_line(aes(y=exp(beta))) +
  geom_ribbon(aes(ymin=exp(beta - 1.96*se), ymax=exp(beta + 1.96*se)),alpha=0.2)+
  xlab("Percentile of CD behaviour PGS")+
  ylab("Odds Ratio (95% CI) of\ndeveloping B2/B3")+
  theme(axis.text = element_text(size=10,face="bold"),
        axis.title = element_text(size=10,face="bold"))


datanew$totalscore_prob <- exp(datanew$totalscore)/(1+exp(datanew$totalscore))

datanew$bev_anno <- "B2/B3"
datanew$bev_anno[datanew$bev==0] <- "B1"

p4 <- ggplot(data=datanew,aes(bev_score,phenoscore_prob,col=bev_anno))+geom_point()+
  theme_bw()+
  theme(axis.title = element_text(face="bold",size=14),
        axis.text = element_text(face="bold",size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=14,face="bold"))+
  xlab("CD behaviour PGS")+
  ylab("Baseline model predicted absolute risk of CD B2/B3")+
  scale_color_manual(values=c( "#E69F00", "#56B4E9"))

ggsave(file="figures/PRS_vs_baseline_model.png",width=10,height=8)



