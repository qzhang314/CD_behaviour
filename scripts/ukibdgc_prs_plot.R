library(dplyr)
library(pROC)
library(boot)
library(ggcorrplot)
library(Hmisc)


setwd("/nfs/users/nfs_q/qz2/CD_behaviour")

data_noNA <- read.csv("result/ukibdgc_fit_score.csv")

#estimate the genetic correlation between triats
data_noNA_pgs <- data_noNA[,c("IBD_PRS","CD_PRS","UC_PRS","CDvsUC_PRS","CD_bev_PRS")]
colnames(data_noNA_pgs) <- c("IBD risk PGS","CD risk PGS","UC risk PGS","IBD subtype PGS","CD behaviour PGS")
cormat <- cor(data_noNA_pgs)

p1 <- ggcorrplot(cormat, hc.order = F, type = "lower",
                 outline.col = "black",colors = c("#6D9EC1", "white", "#E46726"),
                 lab = TRUE,lab_size=5)+
  theme(axis.text=element_text(size=10,face="bold"),
        legend.title = element_text(size=10,face="bold"),
        legend.text = element_text(size=10,face="bold"),
        legend.position = "top")


smy <- read.csv("result/ukibdgc_prs_assoc.csv")
smy <- smy[ smy$cov %in% c("all covs (including PCs)"),]
smy <- smy[smy$var %in% c("IBD_PRS","CD_PRS","UC_PRS","CDvsUC_PRS","CD_bev_PRS"),]
smy$prs <- c("IBD risk PGS","CD risk PGS","UC risk PGS","IBD subtype PGS","CD bebaviour PGS")

smy$prs <- factor(smy$prs,levels=c("IBD risk PGS","CD risk PGS","UC risk PGS","IBD subtype PGS","CD bebaviour PGS"))

p2 <- ggplot(data=smy,aes(prs,beta,col=prs))+geom_point(position = position_dodge(width = 0.2),size=4)+
  geom_errorbar(aes(ymin=(beta-1.96*se),ymax=(beta+1.96*se)),width=0.1,position = position_dodge(width = 0.2))+
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


data_noNA$tile1 <- ntile(data_noNA$CD_bev_PRS, 50)

data_noNA$tile2 <- ntile(data_noNA$pred_score, 50)

prop1s <- c()
for(i in 1:50){
  prop1s <- c(prop1s,mean(data_noNA[data_noNA$tile1==i,"bev"]))
}

prop2s <- c()
for(i in 1:50){
  prop2s <- c(prop2s,mean(data_noNA[data_noNA$tile2==i,"bev"]))
}

dataplot <- data.frame(props=c(prop1s,prop2s),idx=c(1:50,1:50+0.1),type=rep(c("CD behaviour PGS","CD-ProBe"),each=50))

p3 <- ggplot(data=dataplot[ dataplot$type=="CD behaviour PGS",],aes(idx,props))+geom_point()+theme_bw()+
  geom_smooth()+theme(axis.text = element_text(size=10,face="bold"),
                      axis.title = element_text(size=10,face="bold"))+
  scale_y_continuous(breaks=seq(0.2,0.8,0.1))+
  ylab("Proportion of B2/B3\nwithin each bin")+
  xlab("CD behaviour genetics score bin (rank-based, 50 groups)")


# smys1 <- read.csv("result/ukibdgc_prs_top_vs_remaining.csv")
# smys1$type <- "CD behaviour PGS"
# 
# p3 <- ggplot(data=smys1,aes(thresh,exp(beta),group=type))+geom_point()+
#   scale_y_log10()+theme_bw()+
#   geom_line() +
#   geom_ribbon(aes(ymin=exp(beta - 1.96*se), ymax=exp(beta + 1.96*se)),alpha=0.1)+
#   xlab("Percentile of CD behaviour PGS")+
#   ylab("Odds Ratio (95% CI) of \ndeveloping B2/B3")+
#   theme(axis.text = element_text(size=10,face="bold"),
#         axis.title = element_text(size=10,face="bold"),
#         legend.title = element_blank(),
#         legend.position = "top",
#         legend.text = element_text(size=10,face="bold"))+
#   geom_hline(yintercept = 1,linetype=2)
# 
# 
p <- ggarrange(p1,ggarrange(p2,p3,nrow=2,labels = c("B","C")),labels = c("A",""),ncol=2,widths = c(0.5,0.6))

ggsave(file="figures/ukibdgc_PRS_plot.pdf",width=10,height=6)

ggsave(file="figures/ukibdgc_PRS_plot.png",dpi=300,width=10,height=6)


data_noNA$tile1 <- ntile(data_noNA$CD_bev_PRS, 100)
data_noNA$pred_score_risk <- 1/(1+exp(-data_noNA$pred_score))

sum(data_noNA[ data_noNA$tile1==100 & data_noNA$bev==1,"pred_score_risk"] < 0.5)
# [1] 8

data_noNA$bev_anno <- "B2/B3"
data_noNA$bev_anno[data_noNA$bev==0] <- "B1"

p4 <- ggplot(data=data_noNA,aes(CD_bev_PRS,pred_score_risk,col=bev_anno))+geom_point()+
  theme_bw()+
  theme(axis.title = element_text(face="bold",size=14),
        axis.text = element_text(face="bold",size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=14,face="bold"))+
  xlab("CD behaviour PGS")+
  ylab("CD-ProBe predicted absolute risk of developing B2/B3")+
  scale_color_manual(values=c( "#E69F00", "#56B4E9"))+
  geom_hline(yintercept = 0.5,linetype=2)+
  geom_vline(xintercept = min(data_noNA[data_noNA$tile1==100,"CD_bev_PRS"]),linetype=2)

ggsave(file="figures/PRS_vs_baseline_model.png",width=10,height=8)




# 
# 
# 
# data13_noNA<-data13[!is.na(data13$phenoscore_notime),]
# 
# set.seed(314)
# 
# idx <- 1:nrow(data13_noNA)
# datanew <- c()
# #use 10-fold cross validation to explore the benifit of adding PGS in the model
# 
# num <- 10
# for(i in 1:num){
#   if(i!=num){
#     idx_sel <- sample(idx,round(nrow(data13_noNA)/num),replace=F)
#   }else{
#     idx_sel <- idx
#   }
#   idx <- idx[ !(idx %in% idx_sel)]
#   datatrain <- data13_noNA[-idx_sel,]
#   datatest <- data13_noNA[idx_sel,]
#   newmodel <- glm(bev ~ bev_score  + phenoscore_notime, data=datatrain,family="binomial")
#   datatest$totalscore <- predict(newmodel,newdata=datatest)
#   datanew <- rbind(datanew,datatest)
# }
# 
# rocobj1<-roc(datanew$bev, datanew$phenoscore_notime,ci=T,stratified=FALSE,plot=T,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
# 
# rocobj2<-roc(datanew$bev, datanew$totalscore,ci=T,stratified=FALSE,plot=T,print.auc=F, show.thres=TRUE,font=2,font.lab=2,cex=1.3,col="lightblue")
# 
# roc.test(rocobj2,rocobj1,method="delong")
# 
# 
