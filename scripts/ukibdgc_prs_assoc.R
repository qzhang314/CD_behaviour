library(dplyr)

setwd("/nfs/users/nfs_q/qz2/CD_behaviour")

data <- read.table("~/CD_behaviour/pheno/ukibdgc_pheno.txt",header=T)

prs1 <- read.csv("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/PRS/score/UKIBDGC_PRS_EUR.csv")
prs2 <- read.csv("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/score/UKIBDGC_PRS.csv")

identical(prs1$SPID,prs2$SPID)
# [1] TRUE

prs <- cbind(prs1,prs2$CD_bev_PRS)
colnames(prs)[ncol(prs)] <- "CD_bev_PRS"

data <- merge(data,prs,by.x="FID",by.y="SPID")

smy_pc <- c()

#only PCs as covariates
coefs <- summary(glm(bev ~ scale(IBD_PRS) + PC1 + PC2 + PC3 + PC4,data=data,family="binomial"))$coefficients
smy_pc <- rbind(smy_pc,coefs[2,])
coefs <- summary(glm(bev ~ scale(CD_PRS) + PC1 + PC2 + PC3+ PC4,data=data,family="binomial"))$coefficients
smy_pc <- rbind(smy_pc,coefs[2,])
coefs <- summary(glm(bev ~ scale(UC_PRS) + PC1 + PC2 + PC3+ PC4,data=data,family="binomial"))$coefficients
smy_pc <- rbind(smy_pc,coefs[2,])
coefs <- summary(glm(bev ~ scale(CDvsUC_PRS) + PC1 + PC2 + PC3+ PC4,data=data,family="binomial"))$coefficients
smy_pc <- rbind(smy_pc,coefs[2,])
coefs <- summary(glm(bev ~ scale(CD_shaPRS_PRS) + PC1 + PC2 + PC3+ PC4,data=data,family="binomial"))$coefficients
smy_pc <- rbind(smy_pc,coefs[2,])
coefs <- summary(glm(bev ~ scale(UC_shaPRS_PRS) + PC1 + PC2 + PC3+ PC4,data=data,family="binomial"))$coefficients
smy_pc <- rbind(smy_pc,coefs[2,])
coefs <- summary(glm(bev ~ scale(CD_bev_PRS) + PC1 + PC2 + PC3+ PC4,data=data,family="binomial"))$coefficients
smy_pc <- rbind(smy_pc,coefs[2,])
smy_pc <- as.data.frame(smy_pc)
colnames(smy_pc) <- c("beta","se","z","pval")
smy_pc$var <- c("IBD_PRS","CD_PRS","UC_PRS","CDvsUC_PRS","CD_shaPRS_PRS","UC_shaPRS_PRS","CD_bev_PRS")



smy_cov <- c()
coefs <- summary(glm(bev ~ scale(IBD_PRS) + PC1 + PC2 + PC3 + PC4 + smoke_ed + location + location1 + era + follow_up,data=data,family="binomial"))$coefficients
smy_cov <- rbind(smy_cov,coefs[2,])
coefs <- summary(glm(bev ~ scale(CD_PRS) + PC1 + PC2 + PC3+ PC4 + smoke_ed + location + location1 + era + follow_up,data=data,family="binomial"))$coefficients
smy_cov <- rbind(smy_cov,coefs[2,])
coefs <- summary(glm(bev ~ scale(UC_PRS) + PC1 + PC2 + PC3+ PC4 + smoke_ed + location + location1 + era + follow_up,data=data,family="binomial"))$coefficients
smy_cov <- rbind(smy_cov,coefs[2,])
coefs <- summary(glm(bev ~ scale(CDvsUC_PRS) + PC1 + PC2 + PC3+ PC4 + smoke_ed + location + location1 + era + follow_up,data=data,family="binomial"))$coefficients
smy_cov <- rbind(smy_cov,coefs[2,])
coefs <- summary(glm(bev ~ scale(CD_shaPRS_PRS) + PC1 + PC2 + PC3+ PC4 + smoke_ed + location + location1 + era + follow_up,data=data,family="binomial"))$coefficients
smy_cov <- rbind(smy_cov,coefs[2,])
coefs <- summary(glm(bev ~ scale(UC_shaPRS_PRS) + PC1 + PC2 + PC3+ PC4 + smoke_ed + location + location1 + era + follow_up,data=data,family="binomial"))$coefficients
smy_cov <- rbind(smy_cov,coefs[2,])
coefs <- summary(glm(bev ~ scale(CD_bev_PRS) + PC1 + PC2 + PC3+ PC4 + smoke_ed + location + location1 + era + follow_up,data=data,family="binomial"))$coefficients
smy_cov <- rbind(smy_cov,coefs[2,])
smy_cov <- as.data.frame(smy_cov)
colnames(smy_cov) <- c("beta","se","z","pval")
smy_cov$var <- c("IBD_PRS","CD_PRS","UC_PRS","CDvsUC_PRS","CD_shaPRS_PRS","UC_shaPRS_PRS","CD_bev_PRS")

smy_pc$cov <- "PCs"
smy_cov$cov <- "all covs (including PCs)"

smy <- rbind(smy_pc,smy_cov)

write.csv(smy,file="result/ukibdgc_prs_assoc.csv",row.names = F,quote=F)


# # top quantile v..s remaining
# 
# data_noNA <- data[!is.na(data$era) & !is.na(data$location) & !is.na(data$location1) & !is.na(data$smoke_ed) & !is.na(data$follow_up),]
# 
# data_noNA$tile1 <- ntile(data_noNA$CD_bev_PRS, 100)
# 
# smy <- c()
# 
# for(i in c(50:99)){
#   coefs <- summary(glm(bev ~ as.double(tile1 > i) + PC1 + PC2 + PC3 + PC4 + smoke_ed + location + location1 + era + follow_up,data=data_noNA,family="binomial"))$coefficients
#   smy <- rbind(smy,coefs[2,])
# }


data_noNA <- data[!is.na(data$era) & !is.na(data$location) & !is.na(data$location1) & !is.na(data$smoke_ed) & !is.na(data$follow_up),]
data_noNA$tile1 <- ntile(data_noNA$CD_bev_PRS, 50)

smy <- c()

for(i in c(1:49)){
  coefs <- summary(glm(bev ~ as.double(tile1 > i) + PC1 + PC2 + PC3 + PC4 ,data=data_noNA,family="binomial"))$coefficients
  smy <- rbind(smy,coefs[2,])
}

smy <- as.data.frame(smy)
colnames(smy) <- c("beta","se","z","pval")
smy$thresh <- 1:49

write.csv(smy,file="result/ukibdgc_prs_top_vs_remaining.csv",row.names = F,quote=F)


