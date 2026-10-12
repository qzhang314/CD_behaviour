library(pROC)
library(ggplot2)
library(ggpubr)

wkdir<- "~/CD_behaviour/"
setwd(wkdir)

predmodel <- readRDS("result/pheno_model.rds")

# Coefficients:
#   (Intercept)                  timeSinceDiag_5yrs  as.factor(smoker_type)Never-smoker     as.factor(smoker_type)Ex-smoker  
# 0.01895                             0.08091                            -0.17290                            -0.23858  
# year_diag_section2002-2011          year_diag_sectionpost-2011                          locationL2                          locationL3  
# -0.28417                            -0.69363                            -1.23875                            -0.13312  
# location1L4  
# 0.22386  

data <- read.csv("result/ukibdgc_fit_score.csv")

data$pred_loc <- NA
data[!is.na(data$location),"pred_loc"] <- 0.01895
data[which(data$location=="L2"),"pred_loc"] <- data[which(data$location=="L2"),"pred_loc"] + -1.23875 
data[which(data$location=="L3"),"pred_loc"] <- data[which(data$location=="L3"),"pred_loc"] + -0.13312

data$pred_loc1 <- NA
data[!is.na(data$location1),"pred_loc1"] <- 0.01895
data[which(data$location1=="L4"),"pred_loc1"] <- data[which(data$location1=="L4"),"pred_loc1"] + 0.22386

data$pred_timesincediag <- data$timeSinceDiag_5yrs

data$pred_smoke <- NA
data[!is.na(data$smoker_type),"pred_smoke"] <- 0.01895
data[which(data$smoker_type=="Never-smoker"),"pred_smoke"] <- data[which(data$smoker_type=="Never-smoker"),"pred_smoke"] + -0.17290 
data[which(data$smoker_type=="Ex-smoker"),"pred_smoke"] <- data[which(data$smoker_type=="Ex-smoker"),"pred_smoke"] + -0.23858

data$pred_era <- NA
data[!is.na(data$year_diag_section),"pred_era"] <- 0.01895
data[which(data$year_diag_section=="2002-2011"),"pred_era"] <- data[which(data$year_diag_section=="2002-2011"),"pred_era"] + -0.28417
data[which(data$year_diag_section=="post-2011"),"pred_era"] <- data[which(data$year_diag_section=="post-2011"),"pred_era"] + -0.69363



aucs <- c()
rocobj<-roc(data$bev, data$pred_score,ci=T,stratified=FALSE,plot=F,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
tmp <- c(rocobj$auc,rocobj$ci[1],rocobj$ci[3])
aucs <- rbind(aucs,tmp)
rocobj<-roc(data$bev, data$pred_loc,ci=T,stratified=FALSE,plot=F,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
tmp <- c(rocobj$auc,rocobj$ci[1],rocobj$ci[3])
aucs <- rbind(aucs,tmp)
rocobj<-roc(data$bev, data$pred_loc1,ci=T,stratified=FALSE,plot=F,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
tmp <- c(rocobj$auc,rocobj$ci[1],rocobj$ci[3])
aucs <- rbind(aucs,tmp)
rocobj<-roc(data$bev, data$pred_timesincediag,ci=T,stratified=FALSE,plot=F,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
tmp <- c(rocobj$auc,rocobj$ci[1],rocobj$ci[3])
aucs <- rbind(aucs,tmp)
rocobj<-roc(data$bev, data$pred_era,ci=T,stratified=FALSE,plot=F,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
tmp <- c(rocobj$auc,rocobj$ci[1],rocobj$ci[3])
aucs <- rbind(aucs,tmp)
rocobj<-roc(data$bev, data$pred_smoke,ci=T,stratified=FALSE,plot=F,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
tmp <- c(rocobj$auc,rocobj$ci[1],rocobj$ci[3])
aucs <- rbind(aucs,tmp)
aucs <- as.data.frame(aucs)

colnames(aucs)<-c("auc","downci","upci")
aucs$var <- c("Full model","Location","Upper GI involvement","Time since diagnosis","Diagnosis era","Smoking at diagnosis")
aucs$var <- factor(aucs$var,levels=c("Full model","Location","Upper GI involvement","Time since diagnosis","Diagnosis era","Smoking at diagnosis"))

p1 <- ggplot(data=aucs,aes(var,auc,fill=var))+geom_col(col="black")+
  geom_errorbar(aes(ymin=downci,ymax=upci),width=0.2)+
  theme_bw()+
  theme(axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_blank(),
        axis.text = element_text(size=12,face="bold"),
        legend.position = "none",
        strip.text = element_text(size=12,face="bold"))+
  scale_fill_brewer(palette = "Set2")+
  coord_flip()+
  geom_hline(yintercept = 0.5,lty=2)+
  ylab("AUC (95% CI) in predicting CD behaviour")


data$pred_score_risk <- 1 / (1 + exp(-data$pred_score))

props1 <- c()
props2 <- c()
for(i in seq(0.10,0.78,0.05)){
  thresh1 <- i
  thresh2 <- i+0.05
  props1 <- c(props1,mean(data[ data$pred_score_risk > thresh1 & data$pred_score_risk <= thresh2 & !is.na(data$pred_score_risk), "bev"]))
  props2 <- c(props2,mean(data[ data$pred_score_risk > thresh1 & data$pred_score_risk <= thresh2 & !is.na(data$pred_score_risk), "pred_score_risk"],na.rm=T))
}


dataplot <- data.frame(props=props1,idx=props2)

p2 <- ggplot(data=dataplot,aes(idx,props))+
  geom_point(size=3,col="violetred")+
  xlab("Mean predicted absolute risk of CD B2/B3")+
  ylab("Observed CD B2/B3 risk")+
  theme_bw()+
  theme(axis.title = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        legend.position = c(0.2,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        strip.text = element_text(size=12,face="bold"))+
  scale_color_brewer(palette = "Set2")+
  geom_abline(intercept = 0,slope = 1,lty=2)+
  geom_smooth(method="lm",se=F,lty=2,col="violetred")


data_L1 <- data[ data$location=="L1",]
rocobj<-roc(data_L1$bev, data_L1$pred_score,ci=T,stratified=FALSE,plot=F,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")

# Call:
#   roc.default(response = data_L1$bev, predictor = data_L1$pred_score,     ci = T, plot = F, stratified = FALSE, print.auc = F, show.thres = TRUE,     font = 2, font.lab = 2, cex = 1.3, col = "lightblue")
# 
# Data: data_L1$pred_score in 404 controls (data_L1$bev 0) < 497 cases (data_L1$bev 1).
# Area under the curve: 0.6589
# 95% CI: 0.6231-0.6948 (DeLong)


data_L1$bev_cate <- "B1"
data_L1$bev_cate[data_L1$bev==1] <- "B2/B3"

p3 <- ggplot(data=data_L1,aes(pred_score_risk,fill=factor(bev_cate)))+theme_bw()+geom_histogram(position = "dodge")+
  theme(axis.title=element_text(size=12,face="bold"),
        axis.text=element_text(size=12,face="bold"),
        legend.position = c(0.8,0.8),
        legend.title=element_blank(),
        legend.text=element_text(size=12,face="bold"))+
  xlab("Predicted absolute risk of CD B2/B3")+
  ylab("Number of patients")+
  scale_fill_manual(values=c( "#E69F00", "#56B4E9"))

p <- ggarrange(p1,p2,p3,ncol=3,widths = c(0.4,0.3,0.3),labels=c("A","B","C"))

ggsave(file="figures/ukibdgc_all_factor_combined_figure.pdf",dpi=300,width=15,height=4)

ggsave(file="figures/ukibdgc_all_factor_combined_figure.png",dpi=300,width=15,height=4)



