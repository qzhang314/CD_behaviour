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

wkdir<- "/nfs/users/nfs_q/qz2/CD_behaviour"

setwd(wkdir)

dataCD <- read.table("~/CD_behaviour/pheno/ibdbr_pheno.txt",header=T)

#dataCD_nomiss <- dataCD[!is.na(dataCD$timeSinceDiag_5yrs) & !is.na(dataCD$age_diag_section) & !is.na(dataCD$smoker_type) & !is.na(dataCD$year_diag_section) &  !is.na(dataCD$location) &  !is.na(dataCD$location1) &  !is.na(dataCD$sex),]

nrow(dataCD) 
# [1] 13952

nrow(dataCD) - nrow(dataCD_nomiss)
#[1] 1388

nrow(dataCD_nomiss)
#[1] 12564

dataCD_nomiss$smoker_type <- factor(dataCD_nomiss$smoker_type,levels = c("Current-smoker","Never-smoker","Ex-smoker"))
dataCD_nomiss$year_diag_section <- factor(dataCD_nomiss$year_diag_section,levels = c("pre-2002","2002-2011","post-2011"))
dataCD_nomiss$location1 <- factor(dataCD_nomiss$location1,levels = c("noL4","L4"))

full.model<-glm(bev ~ timeSinceDiag_5yrs + age_diag_section + as.factor(smoker_type) + year_diag_section + location + location1+sex,data=dataCD_nomiss,family = "binomial")

step.model <- stepAIC(full.model, direction = "both",trace = FALSE)

saveRDS(step.model,file="result/pheno_model.rds")


smy <- summary(step.model)$coefficients[-1,]
smy <- as.data.frame(smy)
colnames(smy)<-c("beta","se","z","pval")

smy$or <- round(exp(smy$beta),2)
smy$ci <- paste0(round(exp(smy$beta - smy$se*1.96),2),"-",round(exp(smy$beta + smy$se*1.96),2))
smy$orci <- paste0(smy$or,"(",smy$ci,")")
smy$type<-"B2/B3 v.s B1"

smy$variable <- c("Time since diagnosis (5 years)","Never v.s. Current Smoker","Ex v.s. Current Smoker","2002-2011 v.s. pre-2002","post-2011 v.s. pre-2002","Location (L2 v.s. L1)","Location (L3 v.s. L1)","UGI involvement (Yes v.s. No)")

smy$variable <- factor(smy$variable,levels=c("Time since diagnosis (5 years)","Never v.s. Current Smoker","Ex v.s. Current Smoker","2002-2011 v.s. pre-2002","post-2011 v.s. pre-2002","Location (L2 v.s. L1)","Location (L3 v.s. L1)","UGI involvement (Yes v.s. No)"))

write.csv(smy,file="result/phenoeff.csv",quote=F,row.names=F)


p1 <- ggplot(data = smy,aes(variable,exp(beta)))+geom_point()+geom_errorbar(aes(ymin=exp(beta - 1.96*se),ymax=exp(beta + 1.96*se)),width=0.1)+
  ylab("Odds ratio (95% CI)") +theme_bw()+
  theme(legend.title=element_blank(),
                   legend.text=element_text(size=14,face="bold"),
                   axis.title.x=element_text(size=14,face="bold"),
                    axis.title.y=element_blank(),
                   axis.text=element_text(size=14,face="bold"),
                   legend.position = "none")+coord_flip()+geom_hline(yintercept = 1,lty=2)

ggsave(file="figures/pheno_effect_b23_vs_b1_stepwise.png",width=8,height=6,dpi=300)


smyfull <- summary(full.model)$coefficients[-1,]
smyfull <- as.data.frame(smyfull)
colnames(smyfull)<-c("beta","se","z","pval")

smyfull$or <- round(exp(smyfull$beta),2)
smyfull$ci <- paste0(round(exp(smyfull$beta - smyfull$se*1.96),2),"-",round(exp(smyfull$beta + smyfull$se*1.96),2))
smyfull$orci <- paste0(smyfull$or,"(",smyfull$ci,")")
smyfull$type<-"B2/B3 v.s B1"
smyfull$variable <- c("Time since diagnosis (5 years)","Age at diagnosis (A2 v.s. A1)","Age at diagnosis (A3 v.s. A1)","Never v.s. Current Smoker","Ex v.s. Current Smoker","2002-2011 v.s. pre-2002","post-2011 v.s. pre-2002","Location (L2 v.s. L1)","Location (L3 v.s. L1)","UGI involvement (Yes v.s. No)","Male v.s. Female")

smyfull_part <- smyfull[smyfull$variable %in% c("Age at diagnosis (A2 v.s. A1)","Age at diagnosis (A3 v.s. A1)","Male v.s. Female"),]

p1_2 <- ggplot(data = smyfull_part,aes(variable,exp(beta)))+geom_point()+geom_errorbar(aes(ymin=exp(beta - 1.96*se),ymax=exp(beta + 1.96*se)),width=0.1)+
  ylab("Odds ratio (95% CI)") +theme_bw()+
  theme(legend.title=element_blank(),
                   legend.text=element_text(size=14,face="bold"),
                   axis.title.x=element_text(size=14,face="bold"),
                    axis.title.y=element_blank(),
                   axis.text=element_text(size=14,face="bold"),
                   legend.position = "none")+coord_flip()+geom_hline(yintercept = 1,lty=2)

ggsave(file="figures/pheno_effect_b23_vs_b1_age_sex.png",width=8,height=3,dpi=300)


#variance explained by all associated variables, on "observed scale" (treating case/control as 0/1)
var_exp <- anova(lm(as.double(behaviour!=1) ~ timeSinceDiag_5yrs + as.factor(smoker_type) + year_diag_section + location + location1,data=dataCD))
var_exp <- as.data.frame(var_exp)
var_exp$var_pop <- var_exp[,2]/sum(var_exp[,2])
var_exp$variable <- c("Time since diagnosis","Smoking","Diagnosis era","Location","UGI involvement","others")
var_exp <- rbind(var_exp,c("",sum(var_exp[-6,2]),0,"","",sum(var_exp[-6,"var_pop"]),"Total"))
var_exp$var_pop <- as.double(var_exp$var_pop)
var_exp <- var_exp[order(var_exp$var_pop),]
var_exp$variable <- factor(var_exp$variable,levels = rev(var_exp$variable))



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

var_exp$var_pop_lib <- R2_obs2lib(0.5,mean(dataCD_nomiss$behaviour!=1),var_exp$var_pop)

p2 <- ggplot(data=var_exp[-7,],aes(variable,var_pop_lib,fill=variable))+
  geom_bar(stat = 'identity',col="black",alpha=0.6)+
  theme_bw()+
  theme(axis.text= element_text(face="bold",size=14),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(face="bold",size=14),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=9,face="bold"))+
  ylab("Proportion of CD behaviour variance explained by each factor")+
scale_fill_brewer(palette="Set2")+scale_y_continuous(expand = c(0, 0), limits = c(0, 0.18))+coord_flip()


ggsave(file="figures/pheno_effect_var_explained.png",width=9,height=6,dpi=300)

save(smy,var_exp,file="result/model_fit_plot.RObject")

write.csv(var_exp,file="result/var_exp.csv",quote=F,row.names=F)

