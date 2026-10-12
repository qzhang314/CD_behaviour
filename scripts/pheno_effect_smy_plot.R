library(ggplot2)
library(ggpubr)

source("/nfs/users/nfs_q/qz2/R_functions/func.R")

wkdir<- "~/CD_behaviour"
setwd(wkdir)

#load variables generated in model_fit.Rmd

load("result/model_fit_plot.RObject")

smy$variable_fmt <- c("Time since diagnosis\n(5 years)","Smoking at diagnosis\n(Never v.s. Current Smoker)","Smoking at diagnosis\n(Ex v.s. Current Smoker)","Diagnosis era\n(2002-2011 v.s. pre-2002)","Diagnosis era\n(post-2011 v.s. pre-2002)","Location\n(L2 v.s. L1)","Location\n(L3 v.s. L1)","UGI involvement\n(Yes v.s. No)")

smy$variable_fmt <- factor(smy$variable_fmt ,levels = c("Time since diagnosis\n(5 years)","Smoking at diagnosis\n(Never v.s. Current Smoker)","Smoking at diagnosis\n(Ex v.s. Current Smoker)","Diagnosis era\n(2002-2011 v.s. pre-2002)","Diagnosis era\n(post-2011 v.s. pre-2002)","Location\n(L2 v.s. L1)","Location\n(L3 v.s. L1)","UGI involvement\n(Yes v.s. No)"))


p1 <- ggplot(data = smy,aes(variable_fmt,exp(beta)))+geom_point()+geom_errorbar(aes(ymin=exp(beta - 1.96*se),ymax=exp(beta + 1.96*se)),width=0.1)+
  ylab("Odds ratio (95% CI)") +theme_bw()+
  theme(legend.title=element_blank(),
                   legend.text=element_text(size=10,face="bold"),
                   axis.title.x=element_text(size=10,face="bold"),
                    axis.title.y=element_blank(),
                   axis.text=element_text(size=10,face="bold"),
                   legend.position = "none")+coord_flip()+geom_hline(yintercept = 1,lty=2)+
  scale_y_log10()

var_exp$variable <- as.character(var_exp$variable)
var_exp$variable[3] <- "Smoking at diagnosis" 
var_exp$variable <- factor(var_exp$variable,levels=rev(var_exp$variable))

p2 <- ggplot(data=var_exp[-7,],aes(variable,var_pop_lib,fill=variable))+
  geom_bar(stat = 'identity',col="black",alpha=0.6)+
  theme_bw()+
  theme(axis.text= element_text(face="bold",size=10),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(face="bold",size=10),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=9,face="bold"))+
  ylab("Proportion of CD behaviour variance explained by each factor")+
scale_fill_brewer(palette="Set2")+scale_y_continuous(expand = c(0, 0), limits = c(0, 0.18))+coord_flip()


#delayed diagnosis

dataplotm <- readRDS("result/delayed_diagnosis_plot.rds")

p3<-ggplot(data=dataplotm,aes(cate,value,fill=variable))+geom_rect(aes(ymax=ymax,ymin=ymin, xmax=xmax,xmin=xmin),col="black")+theme_bw()+scale_fill_brewer(palette = 2,direction=-1, na.translate = F)+
  theme(axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size=10,face="bold"),
        panel.border = element_blank())+
  scale_y_continuous(labels=scales::percent) +
  ylab("Proportion")+
  xlab("Time interval (years) between date of first symptom and date of diagnosis")+
  coord_polar(theta = "y")+
  theme(panel.grid=element_blank())+
  annotate(geom = 'text', x = 0.5, y = 0, label = "Delayed diagnosis > 1.5 \nyears v.s. others (inner)",size=3.5)


# effect of drug prescription time on Crohn's disease progression (focused on patients with drug taken)
smy_drug <- read.csv("result/drug_presb_effect.csv")

p4 <- ggplot(data = smy_drug,aes(type,exp(beta)))+geom_point()+geom_errorbar(aes(ymin=exp(beta - 1.96*se),ymax=exp(beta + 1.96*se)),width=0.1)+
  ylab("Odds ratio (95% CI) of Early v.s. Late prescription") +theme_bw()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=10,face="bold"),
        axis.title.y=element_blank(),
        axis.text=element_text(size=10,face="bold"))+
  coord_flip()+geom_hline(yintercept = 1,lty=2) + theme(plot.margin = margin(0,0,0,10))



#combine p1, p2, p3 and p4 as figure 2 in the manuscript
p <- ggarrange(p1,p3,p2,p4,nrow=2,ncol=2,labels = c("A","C","B","D"),widths = c(0.6,0.45)) + theme(plot.margin = margin(0, 0, 0, 15))

ggsave(file="figures/pheno_effect_b23_vs_b1_full.png",dpi=300,width=15,height=9,bg = "transparent")

ggsave(file="figures/pheno_effect_b23_vs_b1_full.pdf",width=15,height=9,bg = "transparent")






