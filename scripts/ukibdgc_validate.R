library(dplyr)
library(pROC)
library(boot)

setwd("/nfs/users/nfs_q/qz2/CD_behaviour")

data <- read.table("~/CD_behaviour/pheno/ukibdgc_pheno.txt",header=T)

prs1 <- read.csv("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/PRS/score/UKIBDGC_PRS_EUR.csv")
prs2 <- read.csv("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/score/UKIBDGC_PRS.csv")

identical(prs1$SPID,prs2$SPID)
# [1] TRUE

prs <- cbind(prs1,prs2$CD_bev_PRS)
colnames(prs)[ncol(prs)] <- "CD_bev_PRS"

data <- merge(data,prs,by.x="FID",by.y="SPID")

data_noNA <- data[!is.na(data$era) & !is.na(data$location) & !is.na(data$location1) & !is.na(data$smoke_ed) & !is.na(data$follow_up),]

#AUC of PRS
rocobj <- roc(data_noNA$bev, data_noNA$CD_bev_PRS,ci=T,stratified=FALSE,plot=T,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
# Call:
#   roc.default(response = data_noNA$bev, predictor = data_noNA$CD_bev_PRS,     ci = T, plot = T, stratified = FALSE, print.auc = F, show.thres = TRUE,     font = 2, font.lab = 2, cex = 1.3, col = "lightblue")
# 
# Data: data_noNA$CD_bev_PRS in 1622 controls (data_noNA$bev 0) < 1189 cases (data_noNA$bev 1).
# Area under the curve: 0.5587
# 95% CI: 0.5372-0.5802 (DeLong)

#variance explained
summary(aov(lm(bev ~ CD_bev_PRS, data=data_noNA)))
# Df Sum Sq Mean Sq F value   Pr(>F)    
# CD_bev_PRS     1   10.2  10.236   42.55 8.16e-11 ***
#   Residuals   2809  675.8   0.241          

10.2/(675.8 + 10.2)
# [1] 0.0148688

# K: Population prevalence of the disease
# P: Sample case proportion (case / total in your sample)

R2_obs2lib<-function(k,p,R2obs){
  t<-qnorm(1-k)
  z<-dnorm(t)
  C <- (k*(1-k))^2/(z^2 * p*(1-p))
  theta <- z*(p-k)/(k*(1-k))*(z*(p-k)/(k*(1-k))-t)
  R2lib <- R2obs*C/(1 + R2obs*C*theta)
  return(R2lib)
}

R2_obs2lib(0.5,mean(data_noNA$bev==1),10.2/(675.8 + 10.2))
 # [1] 0.02391487


  
data_noNA$smoker_type <- data_noNA$smoke_ed
data_noNA$smoker_type[which(data_noNA$smoker_type=="non-Smoker")] <- "Never-smoker"
data_noNA$smoker_type[which(data_noNA$smoker_type=="Smoker")] <- "Current-smoker"

data_noNA$smoker_type <- factor(data_noNA$smoker_type,levels=unique(data_noNA$smoker_type))

# data_noNA$location <- data_noNA$loc
# data_noNA$location1 <- data_noNA$loc1

data_noNA$year_diag_section <- data_noNA$era

data_noNA$timeSinceDiag_5yrs <- data_noNA$follow_up/5


#load prediction model
model_fit <- readRDS("result/pheno_model.rds")

data_noNA$pred_score <- predict(model_fit,newdata = data_noNA)

data_noNA$pred_score_prob <- predict(model_fit,newdata = data_noNA,type="response")

roc_pheno <- roc(bev ~ pred_score, data=data_noNA, ci = TRUE)

# Call:
#   roc.formula(formula = bev ~ pred_score, data = data_noNA, ci = TRUE)
# 
# Data: pred_score in 1622 controls (bev 0) < 1189 cases (bev 1).
# Area under the curve: 0.7607
# 95% CI: 0.743-0.7784 (DeLong)



# combine with PRS, use 10-fold CV
set.seed(314)

k <- 10
folds <- sample(rep(1:k, length.out = nrow(data_noNA)))
preds <- rep(NA, nrow(data_noNA))

for (i in 1:k) {
  train <- data_noNA[folds != i, ]
  test  <- data_noNA[folds == i, ]
  
  fit <- glm(bev ~ pred_score + CD_bev_PRS, data = train, family = binomial)
  preds[folds == i] <- predict(fit, newdata = test,type = "response")
}

data_noNA$total_score <- preds

roc_full <- roc(bev ~ total_score, data=data_noNA,ci=T)

# Call:
#   roc.formula(formula = bev ~ total_score, data = data_noNA, ci = T)
# 
# Data: total_score in 1622 controls (bev 0) < 1189 cases (bev 1).
# Area under the curve: 0.7638
# 95% CI: 0.7462-0.7813 (DeLong)



#use different method to compare the AUC
roc.test(roc_full,roc_pheno, method = "delong")
# DeLong's test for two correlated ROC curves
# 
# data:  roc_full and roc_pheno
# Z = 1.3662, p-value = 0.1719
# alternative hypothesis: true difference in AUC is not equal to 0
# 95 percent confidence interval:
#  -0.001328348  0.007441724
# sample estimates:
# AUC of roc1 AUC of roc2 
#   0.7637618   0.7607051 


roc.test(roc_full,roc_pheno, method = "bootstrap", boot.n = 2000)
# Bootstrap test for two correlated ROC curves
# 
# data:  roc_full and roc_pheno
# D = 1.3416, boot.n = 2000, boot.stratified = 1, p-value = 0.1797
# alternative hypothesis: true difference in AUC is not equal to 0
# sample estimates:
#   AUC of roc1 AUC of roc2 
# 0.7637618   0.7607051 

roc.test( roc_full,roc_pheno, method = "venkatraman")
# Venkatraman's test for two paired ROC curves
# 
# data:  roc_full and roc_pheno
# E = 21762, boot.n = 2000, p-value = 0.4475
# alternative hypothesis: true difference in at least one ROC operating point is not equal to 0

write.csv(data_noNA,file="result/ukibdgc_fit_score.csv",row.names=F,quote=F)






#how above excluding follow-up time?
data_noNA$timeSinceDiag_5yrs <- 0
data_noNA$pred_score <- predict(model_fit,newdata = data_noNA)

roc_pheno <- roc(bev ~ pred_score, data=data_noNA, ci = TRUE)
roc_pheno
# Call:
#   roc.formula(formula = bev ~ pred_score, data = data_noNA, ci = TRUE)
# 
# Data: pred_score in 1622 controls (bev 0) < 1189 cases (bev 1).
# Area under the curve: 0.7389
# 95% CI: 0.7207-0.7571 (DeLong)








#plot ROC

png(file="figures/ukibdgc_auc.png",units="in", width=8, height=8, res=300)

library(pROC)

par(
  mar = c(4.5, 4.5, 2, 1),
  las = 1,
  bty = "l"
)

plot(
  rocobj,
  col = "#4C78A8",
  lwd = 3,
  legacy.axes = TRUE,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  cex.lab = 1.2,
  cex.axis = 1.1
)

plot(roc_pheno, col = "#F58518", lwd = 3, add = TRUE, legacy.axes = TRUE)
plot(roc_full,  col = "#54A24B", lwd = 3, add = TRUE, legacy.axes = TRUE)

# optional reference line
abline(0, 1, lty = 2, col = "grey80", lwd = 1.5)

ci1 <- ci.auc(rocobj)
ci2 <- ci.auc(roc_pheno)
ci3 <- ci.auc(roc_full)

text(
  0.52, 0.30,
  labels = sprintf(
    "behaviour PGS: AUC = %.2f (95%% CI %.2f-%.2f)",
    auc(rocobj), ci1[1], ci1[3]
  ),
  col = "#4C78A8",
  adj = 0,
  cex = 0.95
)

text(
  0.52, 0.22,
  labels = sprintf(
    "CD-ProBe: AUC = %.2f (95%% CI %.2f-%.2f)",
    auc(roc_pheno), ci2[1], ci2[3]
  ),
  col = "#F58518",
  adj = 0,
  cex = 0.95
)

text(
  0.52, 0.14,
  labels = sprintf(
    "CD-ProBe + behaviour PGS: AUC = %.2f (95%% CI %.2f-%.2f)",
    auc(roc_full), ci3[1], ci3[3]
  ),
  col = "#54A24B",
  adj = 0,
  cex = 0.95
)
dev.off()

