
source("/nfs/users/nfs_q/qz2/R_functions/func.R")

wkdir<- "~/CD_behaviour"
setwd(wkdir)

dataCD <- read.table("~/CD_behaviour/pheno/ibdbr_pheno.txt",header=T)

dataCD <- dataCD[!is.na(dataCD$timeSinceDiag_5yrs) & !is.na(dataCD$age_diag_section) & !is.na(dataCD$smoker_type) & !is.na(dataCD$year_diag_section) &  !is.na(dataCD$location) &  !is.na(dataCD$location1) &  !is.na(dataCD$sex),]


# Sree suggested that CD patients with older age at diangosis were less likely to have disease progression than patients with younger age at diagnosis (Charpentier et al., GUT, 2013)
# so maybe we should try to explore the interaction effect between age at diagnosis and PRS on disease progression

# this part is to show Sree that the age at diagnosis do have an effect on disease progression in univariate analysis, but the effect diaappares in multivariate analysis.
dataCD$age_diag_section1 <- NA
dataCD[dataCD$age_diag<17 & !is.na(dataCD$age_diag), "age_diag_section1"] <- "< 17"
dataCD[dataCD$age_diag>=17 & dataCD$age_diag<=40 & !is.na(dataCD$age_diag), "age_diag_section1"] <- "17 - 40"
dataCD[dataCD$age_diag>40 & dataCD$age_diag<=60 & !is.na(dataCD$age_diag), "age_diag_section1"] <- "41 - 60"
dataCD[dataCD$age_diag>60 & !is.na(dataCD$age_diag), "age_diag_section1"] <- "> 60"
dataCDnoNA <- dataCD[!is.na(dataCD$age_diag_section1),]
dataCDnoNA$age_diag_section1 <- factor(dataCDnoNA$age_diag_section1,levels=c("< 17","17 - 40","41 - 60","> 60"))

p1 <- ggplot(data=dataCDnoNA,aes(age_diag_section1,fill=as.factor(behaviour)))+geom_bar(position="dodge",alpha=0.8,col="black") + theme_bw()+theme(axis.text = element_text(face="bold",size=12),axis.title=element_text(face="bold",size=12),legend.title=element_text(face="bold",size=12),legend.text=element_text(face="bold",size=12),legend.position="bottom")+scale_fill_brewer("Behaviour",labels=c("B1","B2","B3"),palette = 2)+ylab("Number of patients")+xlab("Age at diagnosis")

v1 <- as.double(table(dataCDnoNA[dataCDnoNA$age_diag_section1=="< 17","behaviour"])/sum(dataCDnoNA$age_diag_section1=="< 17"))
v2 <- as.double(table(dataCDnoNA[dataCDnoNA$age_diag_section1=="17 - 40","behaviour"])/sum(dataCDnoNA$age_diag_section1=="17 - 40"))
v3 <- as.double(table(dataCDnoNA[dataCDnoNA$age_diag_section1=="41 - 60","behaviour"])/sum(dataCDnoNA$age_diag_section1=="41 - 60"))
v4 <- as.double(table(dataCDnoNA[dataCDnoNA$age_diag_section1=="> 60","behaviour"])/sum(dataCDnoNA$age_diag_section1=="> 60"))

dataplot<- data.frame(value=c(v1,v2,v3,v4),behaviour=rep(c(1,2,3),4),age_diag_section1=rep(c("< 17","17 - 40","41 - 60","> 60"),each=3))
dataplot$age_diag_section1 <- factor(dataplot$age_diag_section1,levels=c("< 17","17 - 40","41 - 60","> 60"))

p2 <- ggplot(data=dataplot,aes(x=age_diag_section1,y=value,fill=as.factor(behaviour)))+geom_bar(position="fill", stat="identity",alpha=0.8,col="black") + theme_bw()+theme(axis.text = element_text(face="bold",size=12),axis.title=element_text(face="bold",size=12),legend.title=element_text(face="bold",size=12),legend.text=element_text(face="bold",size=12),legend.position="bottom")+scale_fill_brewer("Behaviour",labels=c("B1","B2","B3"),palette = 2)+ylab("Proportion of patients")+xlab("Age at diagnosis")+scale_y_continuous(labels = scales::percent)

p <- ggarrange(p1,p2,ncol=2,nrow=1)

ggsave("figures/behaviour_aad.png",width=9,height=6,dpi=300)

coefs1 <- summary(glm(as.double(behaviour!=1) ~  age_diag_section1,data=dataCD,family = "binomial"))$coefficients[2:4,]

# Estimate Std. Error   z value     Pr(>|z|)
# age_diag_section1> 60    -0.6217259 0.09055773 -6.865520 6.624959e-12
# age_diag_section117 - 40 -0.1658339 0.05418504 -3.060511 2.209597e-03
# age_diag_section141 - 60 -0.5156542 0.06603720 -7.808541 5.785356e-15

#after correcting for following-up time, disease location and diagosis ear, the age at diagnosis showed no effect on disease progression
coefs2 <- summary(glm(as.double(behaviour!=1) ~  age_diag_section1 + timeSinceDiag_5yrs  + year_diag_section + location + location1 ,data=dataCD,family = "binomial"))
coefs2$coefficients

# Estimate Std. Error      z value      Pr(>|z|)
# (Intercept)                -0.15318762 0.11038469  -1.38776157  1.652097e-01
# age_diag_section1> 60       0.05388218 0.09887866   0.54493233  5.858001e-01
# age_diag_section117 - 40   -0.00412795 0.05795649  -0.07122499  9.432187e-01
# age_diag_section141 - 60   -0.09567959 0.07195717  -1.32967411  1.836257e-01
# timeSinceDiag_5yrs          0.08075750 0.01637852   4.93069498  8.193759e-07
# year_diag_sectionpost-2011 -0.42265289 0.05686906  -7.43203576  1.069389e-13
# year_diag_sectionpre-2002   0.29138348 0.07249202   4.01952507  5.831558e-05
# locationL2                 -1.25405785 0.05260982 -23.83695259 1.382926e-125
# locationL3                 -0.14066234 0.04473678  -3.14422134  1.665294e-03
# location1noL4              -0.20757127 0.08213458  -2.52720939  1.149729e-02

coefs2 <- coefs2$coefficients[2:4,]

coefs1 <- as.data.frame(coefs1)
coefs2 <- as.data.frame(coefs2)
colnames(coefs1) <- c("beta","se","z","pval")
colnames(coefs2) <- c("beta","se","z","pval")
coefs1$type <- "no phenotypic factor as covariates"
coefs2$type <- "phenotypic factor as covariates"
coefs1$variable <- rownames(coefs1)
coefs2$variable <- rownames(coefs2)
coefs <- rbind(coefs1,coefs2)
coefs$or <- exp(coefs$beta)
coefs$ci <- paste0(round(exp(coefs$beta - coefs$se*1.96),2),"-",round(exp(coefs$beta + coefs$se*1.96),2))

write.csv(coefs,file="result/old_age_effect.csv",row.names=F,quote=F)




