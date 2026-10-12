library(data.table)

setwd("~/CD_behaviour/")

read.csv("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/phenotype_data/phenotype_data_freezes/20221101/phenotype_data_UKIBDGCGWAS_20221101.tsv.gz",sep="\t",header=T,stringsAsFactors=F)->data

# keep IBD patients
data <- data[ data$diag %in% c("Crohn's Disease","Indeterminate","Ulcerative Colitis"),]

data$orig_sample_id_new <- gsub("urn:wtsi:","",data$orig_sample_id)

gwas3 <- read.table("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/genotyped_data/QC/all_hce/all_hce_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam")
gwas2 <- read.table("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/genotyped_data/QC/gwas2/gwas2_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam")
gwas1 <- read.table("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/genotyped_data/QC/gwas1/gwas1_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam")
gwas1$batch <- "GWAS1"
gwas2$batch <- "GWAS2"
gwas3$batch <- "GWAS3"

gwas <- rbind(gwas1,gwas2,gwas3)
gwas$ID_new <- gsub("urn:wtsi:","", gwas$V1)

#keeping case
gwas <- gwas[ gwas$V6==2,]

data <- merge(data[,-2],gwas[,c("ID_new","V1","batch")],by.x="orig_sample_id_new",by.y="ID_new")


#remove duplicated, inconsistent
dup1 <- read.table("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/phenotype_data/ukiibdgc_list_duplicated_samples_and_pheno_inconsistencies_to_remove.tsv.gz")


#remove duplicated samples

ID_removed <- unique(dup1$V1)

sum(data$orig_sample_id %in% ID_removed)
# [1] 144

data <- data[!(data$orig_sample_id %in% ID_removed), ]


#remove duplicates and 1st relatives between IBDBR and UKIBDGC
dup1 <- read.table("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/cohort_overlaps/duplicated_pairs_interval_ibdbioresource_ukibdgc_toremove.tsv",header=T,stringsAsFactors = F)

dup2 <- read.table("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/cohort_overlaps/first_degree_related_pairs_interval_ibdbioresource_ukibdgc_toremove.tsv",header=T,stringsAsFactors = F)

dup <- rbind(dup1,dup2[,colnames(dup1)])

ID1_1 <- dup[ dup$cohort_FID1%in%c("all_hce","gwas1","gwas2") & dup$cohort_FID2=="interval_ibdbioresource","FID1"]
ID2_1 <- dup[ dup$cohort_FID1%in%c("all_hce","gwas1","gwas2") & dup$cohort_FID2=="interval_ibdbioresource","FID2"]

ID1_2 <- dup[ dup$cohort_FID2%in%c("all_hce","gwas1","gwas2") & dup$cohort_FID1=="interval_ibdbioresource","FID2"]
ID2_2 <- dup[ dup$cohort_FID2%in%c("all_hce","gwas1","gwas2") & dup$cohort_FID1=="interval_ibdbioresource","FID1"]

idmap <- data.frame(ID1=c(ID1_1,ID1_2),ID2=c(ID2_1,ID2_2))

# remove the overlapped samples between UKIBDGC and IBDBR
sum(data$orig_sample_id %in%idmap$ID1)
# [1] 1814

data <- data[!(data$orig_sample_id %in%idmap$ID1),]


#keep Crohn's disease
data <- data[ which(data$diag=="Crohn's Disease"),]


#define behaviour
data <- data[ data$dis_behavior!="Unknown",] 
data$bev <- as.double(data$dis_behavior %in% c("B2","B3"))


#define age at diagnosis
data$age_diag <- data$diag_yr - data$yob
data$age_diag_section <- NA
data$age_diag_section[which(data$age_diag<17)] <- "A1"
data$age_diag_section[which(data$age_diag>=17 & data$age_diag<=40)] <- "A2"
data$age_diag_section[which(data$age_diag>40)] <- "A3"

#define follow-up time
data$follow_up <- data$review_yr - data$diag_yr

#define CD location
data$location <- NA
data[data$dis_loc_ileal=="Yes" & data$dis_loc_colorectal=="No","location"] <- "L1"
data[data$dis_loc_ileal=="No" & data$dis_loc_colorectal=="Yes","location"] <- "L2"
data[data$dis_loc_ileal=="Yes" & data$dis_loc_colorectal=="Yes","location"] <- "L3"

#define Upper GI involvement
data$location1 <- NA
data[data$dis_loc_gi=="Yes" | data$dis_loc_jejunal=="Yes","location1"] <- "L4"
data[data$dis_loc_gi=="No" & data$dis_loc_jejunal=="No","location1"] <- "noL4"

#define smoking
data$smoke_ed <- NA
data$smoke_ed[which(data$smoking=="Yes")] <- "Current-smoker"
data$smoke_ed[which(data$smoking=="Ex-smoker")] <- "Ex-smoker"
data$smoke_ed[which(data$smoking=="No")] <- "Never-smoker"


datam <- data

datam <- datam[,c("V1","V1","bev","age_diag","age_diag_section","sex","follow_up","location","location1","smoke_ed","batch","diag_yr","dis_behavior")]

datam$era <- NA
datam$era[which(datam$diag_yr<2002)] <- "pre-2002"
datam$era[which(datam$diag_yr>=2002 & datam$diag_yr<=2011)] <- "2002-2011"
datam$era[which(datam$diag_yr>2011)] <- "post-2011"



colnames(datam)[1:2]<- c("FID","IID")

pcs <- fread("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/phenotype_data/ukbibdgc_all_40pcs.tsv.gz",header=T)

datam <- merge(datam,pcs[,c("FID","PC1","PC2","PC3","PC4")],by="FID")

#datam <- datam[!is.na(datam$age_diag) & !is.na(datam$follow_up) & !is.na(datam$location) & 
#		!is.na(datam$location1) & !is.na(datam$smoke_ed) & !is.na(datam$diag_yr),]
#
write.table(datam,file="~/CD_behaviour/pheno/ukibdgc_pheno.txt",row.names=F,quote=F,col.names=T)
