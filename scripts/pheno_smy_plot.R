library(ggplot2)

source("/nfs/users/nfs_q/qz2/R_functions/func.R")

wkdir<- "~/CD_behaviour"
setwd(wkdir)

dataCD <- read.table("~/CD_behaviour/pheno/ibdbr_pheno.txt",header=T)

# dataCD <- dataCD[!is.na(dataCD$timeSinceDiag_5yrs) & !is.na(dataCD$age_diag_section) & !is.na(dataCD$smoker_type) & !is.na(dataCD$year_diag_section) &  !is.na(dataCD$location) &  !is.na(dataCD$location1) &  !is.na(dataCD$sex),]

count <- table(dataCD$behaviour)/nrow(dataCD)
dataplot<-c()
dataplot$ymax <- cumsum(count)
dataplot$ymin <- c(0, head(dataplot$ymax, n=-1))
dataplot<- as.data.frame(dataplot)
dataplot$labelPosition <- (dataplot$ymax + dataplot$ymin) / 2

dataplot$Behaviour<-c("B1","B2","B3")
dataplot$label <- paste0(dataplot$Behaviour, ": ", table(dataCD$behaviour))

dataplot$Behaviour<- factor(dataplot$Behaviour,levels=c("B1","B2","B3"))

p1<-ggplot(dataplot, aes(xmin=3,xmax=4,ymin=ymin,ymax=ymax,fill=Behaviour)) + geom_rect() + coord_polar(theta="y")+ geom_text( x=3, aes(y=labelPosition, label=label), size=6) + theme_bw() + theme(axis.text = element_blank(),strip.text=element_text(face="bold",size=14),axis.title = element_blank(),legend.title=element_text(face="bold",size=14, margin = margin(r = 16)),axis.ticks=element_blank(),legend.text=element_text(face="bold",size=14),panel.grid.major = element_blank(),panel.grid.minor  = element_blank(),panel.border = element_blank())+xlim(c(2, 4))+scale_fill_brewer(palette=2)

ggsave("figures/pheno_samplesize.png",width = 9, height = 6,dpi=300)

dataplot<- dataCD
dataplot$behaviour <- factor(dataplot$behaviour,levels=c(1,2,3))

p2<-ggplot(dataplot, aes(x=behaviour,y=age_diag,fill=as.factor(behaviour))) + geom_violin() + theme_bw() + ylab("Age at diagnosis (year)")+ theme(axis.text.y = element_text(face="bold",size=14),axis.text.x = element_blank(),strip.text=element_text(face="bold",size=14),axis.title.y = element_text(face="bold",size=14),axis.title.x=element_blank(),legend.title=element_blank(),axis.ticks.x = element_blank(),legend.text=element_text(face="bold",size=14))+scale_fill_brewer(labels=c("B1","B2","B3"),palette=2)+geom_boxplot(width=0.1, color="grey", alpha=0.2)


p3<-ggplot(dataplot, aes(x=behaviour,y=follow_up,fill=as.factor(behaviour))) + geom_violin() + theme_bw() + ylab("Time since diagnosis (year)")+ theme(axis.text.y = element_text(face="bold",size=14),axis.text.x = element_blank(),strip.text=element_text(face="bold",size=14),axis.title.y = element_text(face="bold",size=14),axis.title.x=element_blank(),legend.title=element_blank(),axis.ticks.x = element_blank(),legend.text=element_text(face="bold",size=14))+scale_fill_brewer(labels=c("B1","B2","B3"),palette=2)+geom_boxplot(width=0.1, color="grey", alpha=0.2)



props<-c()
for(i in 1:50){
  tmp<- table(dataCD[ceiling(dataCD$follow_up)==i,"behaviour"])/sum(ceiling(dataCD$follow_up)==i,na.rm=T)
  props<- rbind(props,tmp)
}

# tmp<- table(dataCD[ceiling(dataCD$follow_up)>50,"behaviour"])/sum(ceiling(dataCD$follow_up)>50,na.rm=T)
# props<- rbind(props,tmp)
props<- as.data.frame(props)
# props$date<- c(1:51)
props$date<- c(1:50)
colnames(props)<-c("B1","B2","B3","date")
propsm<- melt(props,id.vars = "date")
propsm$date<- factor(propsm$date,levels=unique(propsm$date))
colnames(propsm)[2]<-"behaviour"

propsm$behaviour<- factor(propsm$behaviour,levels=c("B3","B2","B1"))


p4 <- ggplot(data=propsm,aes(x=as.double(date),y=value*100,fill=behaviour))+geom_bar(stat="identity",position = "stack",width = 1)+xlab("Time since diagnosis (year)")+ylab("Proportion of CD behaviors (%)")+theme_bw()+theme(axis.text.x = element_text(face="bold",size=14),axis.text.y = element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),strip.text = element_text(size=14,face="bold"),legend.text=element_text(face="bold",size=14),legend.position = "top",legend.title=element_text(face="bold",size=14))+scale_x_continuous(limits = c(1,51), expand = c(0, 0))+scale_y_continuous(limits = c(0,100.1), expand = c(0, 0))+scale_fill_brewer(labels=c("B1","B2","B3"),breaks = c("B1","B2","B3"),palette=2,direction = -1)


# p4 <- ggplot(data=propsm,aes(x=as.double(date),y=value*100,fill=behaviour))+geom_area(stat="identity",position = "stack")+xlab("Time since diagnosis (year)")+ylab("Proportion of CD behaviors (%)")+theme_bw()+theme(axis.text.x = element_text(face="bold",size=14),axis.text.y = element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),strip.text = element_text(size=14,face="bold"),legend.text=element_text(face="bold",size=14),legend.position = "top",legend.title=element_text(face="bold",size=14))+scale_x_continuous(limits = c(1,51), expand = c(0, 0))+scale_y_continuous(limits = c(0,100.1), expand = c(0, 0))+scale_fill_brewer(labels=c("B1","B2","B3"),breaks = c("B1","B2","B3"),palette=2,direction = -1)


# ggsave("figures/bev_prop.png",width = 9, height = 6,dpi=300)

# p<-ggarrange(p1,p2,p3,p4, nrow = 2,ncol=2,labels = c("A","B","C","D"),common.legend = TRUE,widths = c(0.4,0.6))

p<-ggarrange(p1,p4, nrow = 2,ncol=1,labels = c("A","B"),common.legend = TRUE,widths = c(0.4,0.6))

ggsave("figures/pheno_smy.png",width = 12, height = 9,dpi=300,bg = "white")

ggsave("figures/pheno_smy.pdf",width = 12, height = 9,bg = "white")

