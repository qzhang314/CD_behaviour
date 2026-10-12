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
library(dplyr)
library(stringr)
library(lubridate)



source("/nfs/users/nfs_q/qz2/R_functions/func.R")

wkdir<- "~/CD_behaviour"
setwd(wkdir)

pheno_baseline <- readRDS("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/PhenoQC/2025-09-29/ibdbr_baseline.rds")

dataCD <- read.table("~/CD_behaviour/pheno/ibdbr_pheno.txt",header=T)

dataCD <- dataCD[!is.na(dataCD$timeSinceDiag_5yrs) & !is.na(dataCD$age_diag_section) & !is.na(dataCD$smoker_type) & !is.na(dataCD$year_diag_section) &  !is.na(dataCD$location) &  !is.na(dataCD$location1) &  !is.na(dataCD$sex),]

data_sim <- pheno_baseline[,c("study_participant_id_1","smoke_diag","smoking_history","long_quit_smok")]

dataCD <- merge(dataCD,data_sim,by.x="SP0011_id_1",by.y="study_participant_id_1")


hlq <- readRDS("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/PhenoQC/2026-02-12/ibdbr_hlq.rds")
hlq <- as.data.frame(hlq)

impute_mid_date <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "NA", "Na", "N/A", "NULL")] <- NA_character_
  
  out <- rep(as.Date(NA), length(x))
  
  # only test non-missing values
  not_na <- !is.na(x)
  
  # YYYY
  idx_year <- not_na & str_detect(x, "^\\d{4}$")
  out[idx_year] <- as.Date(paste0(x[idx_year], "-07-01"))
  
  # YYYY-M or YYYY-MM
  idx_ym <- not_na & str_detect(x, "^\\d{4}-\\d{1,2}$")
  if (any(idx_ym)) {
    yr <- str_sub(x[idx_ym], 1, 4)
    mo <- str_extract(x[idx_ym], "(?<=-)\\d{1,2}$")
    out[idx_ym] <- as.Date(sprintf("%s-%02d-15", yr, as.integer(mo)))
  }
  
  # YYYY-M-D or YYYY-MM-DD
  idx_ymd <- not_na & str_detect(x, "^\\d{4}-\\d{1,2}-\\d{1,2}$")
  if (any(idx_ymd)) {
    parts <- str_match(x[idx_ymd], "^(\\d{4})-(\\d{1,2})-(\\d{1,2})$")
    out[idx_ymd] <- as.Date(sprintf(
      "%s-%02d-%02d",
      parts[, 2],
      as.integer(parts[, 3]),
      as.integer(parts[, 4])
    ))
  }
  
  out
}

hlq <- hlq %>%
  mutate(
    date1_imp = impute_mid_date(ibd_diag_yr_month),
    date2_imp = impute_mid_date(ibd_symptoms_yr_month),
    delay_years = as.numeric(date1_imp - date2_imp) / 365.25
  )

#treating delay_years = -1 as 0, treating delay_years< -1 as NA
hlq$delay_years[which(hlq$delay_years >= -1 & hlq$delay_years<0)] <- 0
hlq$delay_years[which(hlq$delay_years < -1)] <- NA

dataCD <- merge(dataCD,hlq[,c("SP0011_id_1","delay_years")],all.x=T)

#check if patients with longer delayed diagnosis are more likely to have progressed disease behaivour. (it turns out patients with at least 2 years late diagnosis were more likely to develop B2/B3 than those with short diagnosis delay, OR = 1.246140)

dataplot <- c()
for(i in seq(0.5,4.5,0.5)){
  tmp <- table(dataCD[ dataCD$delay_years>=i & dataCD$delay_years<i+0.5 & !is.na(dataCD$delay_years), "behaviour"])/sum(dataCD$delay_years>=i & dataCD$delay_years<i+0.5 & !is.na(dataCD$delay_years))
  if(length(tmp)<3){
    tmp <- c(tmp["1"],tmp["2"],tmp["3"])
    names(tmp) <- c(1,2,3)
    tmp[is.na(tmp)]<- 0
  }
  dataplot <- rbind(dataplot,tmp)  
}

tmp <- table(dataCD[ dataCD$delay_years>=5 & !is.na(dataCD$delay_years), "behaviour"])/sum(dataCD$delay_years>=5 & !is.na(dataCD$delay_years))
dataplot <- rbind(dataplot,tmp)  
  
dataplot <- as.data.frame(dataplot)
colnames(dataplot) <- c("B1","B2","B3")
dataplot$cate <- c("[0-0.5)","[0.5-1)","[1-1.5)","[1.5-2)","[2.5-3)","[3-3.5)","[3.5-4)","[4-4.5)","[4.5-5)",">=5")

dataplotm <- reshape2::melt(dataplot)
dataplotm$variable <- factor(dataplotm$variable,levels=c("B3","B2","B1"))

p<-ggplot(data=dataplotm,aes(cate,value,fill=variable))+geom_col()+theme_bw()+scale_fill_brewer(palette = 2,direction=-1)+
  theme(axis.text = element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        legend.title=element_blank(),
        legend.text = element_text(size=12,face="bold"))+
  scale_y_continuous(labels=scales::percent) +
  ylab("Proportion")+
  xlab("Time interval (years) between date of first symptom and date of diagnosis")

ggsave(file="figures/delay_diag.png",width=9,height=6,dpi=300)


#plot 0 - 1.5 and > 1.5 only
dataplot <- c()
tmp <- table(dataCD[ dataCD$delay_years>=0 & dataCD$delay_years<1.5 & !is.na(dataCD$delay_years), "behaviour"])/sum(dataCD$delay_years>=0 & dataCD$delay_years<1.5 & !is.na(dataCD$delay_years))
dataplot <- rbind(dataplot,tmp)
tmp <- table(dataCD[ dataCD$delay_years>=1.5 & !is.na(dataCD$delay_years), "behaviour"])/sum(dataCD$delay_years>=1.5 & !is.na(dataCD$delay_years))
dataplot <- rbind(dataplot,tmp)


dataplot <- as.data.frame(dataplot)
colnames(dataplot) <- c("B1","B2","B3")
dataplot$cate <- c("[0-1.5)",">=1.5")

dataplotm <- reshape2::melt(dataplot)
dataplotm$variable <- factor(dataplotm$variable,levels=c("B3","B2","B1"))
dataplotm <- dataplotm %>%
  dplyr::group_by(cate) %>%
  dplyr::mutate(
    ymax = cumsum(value),
    ymin = dplyr::lag(ymax, default = 0)
  ) %>%
  dplyr::ungroup()

# Add x limits
baseNum <- 3
#numCat <- length(unique(dat$ring))
dataplotm$xmax <- c(1,2,1,2,1,2) + baseNum
dataplotm$xmin = dataplotm$xmax - 1

saveRDS(dataplotm,file="result/delayed_diagnosis_plot.rds")


coefs <-c()
for(i in seq(0.5,5,0.5)){
  dataCD1 <- dataCD[ which(dataCD$delay_years>=i |  dataCD$delay_years==0),]
  smy<-summary(glm(as.double(behaviour!=1) ~ as.double(delay_years>=i) + timeSinceDiag_5yrs  + as.factor(smoker_type) + year_diag_section + location + location1,data=dataCD1,family = "binomial"))
  coefs <- rbind(coefs,smy$coefficients[2,])
}

coefs <- as.data.frame(coefs)
colnames(coefs) <- c("beta","se","z","pval")
coefs$idx <- seq(0.5,5,0.5)
coefs$group<- 1

p <- ggplot(data=coefs,aes(idx,exp(beta),group=group))+
  geom_point(size=2)+
  geom_line(aes(y=exp(beta))) +
  geom_ribbon(aes(ymin=exp(beta - 1.96*se), ymax=exp(beta + 1.96*se)),alpha=0.2)+theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text=element_text(size=12,face="bold"))+
  ylab("Odds ratio (95% CI)")+xlab("Years of delayed diagnosis")+
  geom_abline(intercept = 1,slope=0,lty=2)+
  ylim(c(0.9,1.7))

ggsave(file="figures/delay_diag_CD_behaviour.png",width=8,height=6,dpi=300)

