jobindex <- Sys.getenv('LSB_JOBINDEX')
chrnum<- as.double(jobindex)

library(data.table)
library(dplyr)

setwd("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour")

path<-"/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/"

info<-fread(paste(path,"imputed_data/release_interval_ibdbioresource_Aug24/combined_hrc_uk10k/",chrnum,"_combined_hrc_uk10k_eur_with_consensus_ids_imputation_info_hwe_ctr_marker_names_b37_b38_imputation_info_dbsnp154.tsv.gz",sep=""),head=T)

info<-as.data.frame(info)
info<-info[which(is.na(info$low_conconrdance_freq_hrc_uk10k)),]
info<-info[which(info$HW_lrt_p_value_ctr>1E-5),]
info<-info[which( (is.na(info$unexpected_a1freq_discrepancy_cd_ibdbr_vs_allhce) | info$unexpected_a1freq_discrepancy_cd_ibdbr_vs_allhce==0) &
             (is.na(info$unexpected_a1freq_discrepancy_uc_ibdbr_vs_allhce) | info$unexpected_a1freq_discrepancy_uc_ibdbr_vs_allhce==0) &
             (is.na(info$unexpected_a1freq_discrepancy_ctr_ibdbr_vs_allhce) | info$unexpected_a1freq_discrepancy_ctr_ibdbr_vs_allhce==0)),]
info<-info[which( (is.na(info$unexpected_a1freq_discrepancy_cd_ctr_ibdbr) | info$unexpected_a1freq_discrepancy_cd_ctr_ibdbr==0) &
                (is.na(info$unexpected_a1freq_discrepancy_uc_ctr_ibdbr) | info$unexpected_a1freq_discrepancy_uc_ctr_ibdbr==0)),]
info <- info[info$dbsnp154!="",]

eff <- read.table(paste0("PRScs/chr",chrnum,"_pst_eff_a1_b0.5_phiauto_chr",chrnum,".txt"),header=F,stringsAsFactor=F)
colnames(eff)<- c("CHR","rsID","BP_b37","A1","A2","weight_A1")

info <- info[,c("dbsnp154","id_b38","id_b37","impute_info_cases")]
info$id_b38 <- gsub("_",":",info$id_b38)
info$id_b38 <- paste0("chr",info$id_b38)

info <- info[ info$dbsnp154 %in% eff$rsID,]

info <- info %>%
	group_by(dbsnp154) %>%
	slice_max(impute_info_cases) %>%
	ungroup()

eff <- merge(eff,info,by.y="dbsnp154",by.x="rsID")
eff <- eff[,c("CHR","rsID","BP_b37","A1","A2","weight_A1","id_b37","id_b38")]

write.table(eff, file=paste0("PRScs/weights/chr",chrnum,".coef"),quote=F,row.names=F,col.names=T)
