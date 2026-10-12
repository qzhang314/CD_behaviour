library(dplyr)
library(stringr)
library(lubridate)

source("/nfs/users/nfs_q/qz2/R_functions/func.R")

wkdir<- "~/CD_behaviour"
setwd(wkdir)

dataCD <- read.table("~/CD_behaviour/pheno/ibdbr_pheno.txt",header=T)

dataCD <- dataCD[!is.na(dataCD$timeSinceDiag_5yrs) & !is.na(dataCD$age_diag_section) & !is.na(dataCD$smoker_type) & !is.na(dataCD$year_diag_section) &  !is.na(dataCD$location) &  !is.na(dataCD$location1) &  !is.na(dataCD$sex),]

p <- ggplot(dataCD,aes(follow_up,fill=year_diag_section))+geom_histogram(col="black",bins = 67)+
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 67, by = 1))+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        legend.position="top")+
  scale_fill_brewer(palette = "Set2")


props1 <- c()
props2 <- c()

for(i in 6:13){
  props1 <- c(props1,mean(dataCD[ which(dataCD$year_diag_section=="2002-2011" & dataCD$follow_up==i) ,"bev"]))
  props2 <- c(props2,mean(dataCD[ which(dataCD$year_diag_section=="post-2011" & dataCD$follow_up==i) ,"bev"]))
}


props1 <- c()
props2 <- c()

for(i in 17:23){
  props1 <- c(props1,mean(dataCD[ which(dataCD$year_diag_section=="2002-2011" & dataCD$follow_up==i) ,"bev"]))
  props2 <- c(props2,mean(dataCD[ which(dataCD$year_diag_section=="pre-2002" & dataCD$follow_up==i) ,"bev"]))
}



dat <- dataCD %>%
  mutate(
    event = ifelse(behaviour %in% c(2,3), 1, 0),
    left  = ifelse(event == 1, 0, follow_up),
    right = ifelse(event == 1, follow_up, Inf)
  )
