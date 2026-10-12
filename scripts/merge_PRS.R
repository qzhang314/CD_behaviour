setwd("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour")

cohorts <- c("all_hce","gwas1","gwas2")

CD_bev_total <- c()
for(cohort in cohorts){
	CD_bev<- read.table(paste0("score/",cohort,"/chr1.sscore"),header=F,stringsAsFactors=F)
	CD_bev <- CD_bev[,c("V1","V6")]
	for(i in 2:22){
        	score <- read.table(paste0("score/",cohort,"/chr",i,".sscore"),header=F,stringsAsFactors=F)
        	CD_bev$V6 <- CD_bev$V6 + score$V6
	}
	CD_bev_total <- rbind(CD_bev_total,CD_bev)
}

prs <- data.frame(SPID = CD_bev_total$V1, CD_bev_PRS= CD_bev_total$V6)

write.csv(prs,"score/UKIBDGC_PRS.csv",quote=F,row.names=F)
