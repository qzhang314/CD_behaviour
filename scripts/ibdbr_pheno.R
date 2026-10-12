setwd("~/CD_behaviour/")
source("~/R_functions/func.R")

pheno <- readRDS("/lustre/scratch125/humgen/projects_v2/ibd_bioresource_az/phenotype_data/release_20250912_AZ_project/derived_phenotypes/ibdbr_2025-09-29.rds")

#more detailed CD surgery
#Surgery for colonic disease

## Colectomy and ileostomy
pheno$cd_surgery_colonic_colectomy_ileostomy <- rowSums(pheno[,c("op_1___1","op_2___1","op_3___1","op_4___1","op_5___1")],na.rm=T)
pheno$cd_surgery_colonic_colectomy_ileostomy[which(pheno$diagnosis!=1)]<- NA

pheno$cd_surgery_colonic_colectomy_ileostomy_date<- NA
pheno$cd_surgery_colonic_colectomy_ileostomy_date[which(pheno$op_1___1==1)] <- pheno$surg_date_1[which(pheno$op_1___1==1)]
pheno$cd_surgery_colonic_colectomy_ileostomy_date[which(pheno$op_2___1==1)] <- pheno$surg_date_1[which(pheno$op_2___1==1)]
pheno$cd_surgery_colonic_colectomy_ileostomy_date[which(pheno$op_3___1==1)] <- pheno$surg_date_1[which(pheno$op_3___1==1)]
pheno$cd_surgery_colonic_colectomy_ileostomy_date[which(pheno$op_4___1==1)] <- pheno$surg_date_1[which(pheno$op_4___1==1)]
pheno$cd_surgery_colonic_colectomy_ileostomy_date[which(pheno$op_5___1==1)] <- pheno$surg_date_1[which(pheno$op_5___1==1)]
pheno$cd_surgery_colonic_colectomy_ileostomy_date[is.na(pheno$cd_surgery_colonic_colectomy_ileostomy)] <- NA
pheno$cd_surgery_colonic_colectomy_ileostomy_date <- as.double(pheno$cd_surgery_colonic_colectomy_ileostomy_date)

## Colectomy and ileo-anal pouch
pheno$cd_surgery_colonic_colectomy_ileoanal_pouch <- rowSums(pheno[,c("op_1___2","op_2___2","op_3___2","op_4___2","op_5___2")],na.rm=T)
pheno$cd_surgery_colonic_colectomy_ileoanal_pouch[which(pheno$diagnosis!=1)]<- NA

pheno$cd_surgery_colonic_colectomy_ileoanal_pouch_date<- NA
pheno$cd_surgery_colonic_colectomy_ileoanal_pouch_date[which(pheno$op_1___2==1)] <- pheno$surg_date_1[which(pheno$op_1___2==1)]
pheno$cd_surgery_colonic_colectomy_ileoanal_pouch_date[which(pheno$op_2___2==1)] <- pheno$surg_date_1[which(pheno$op_2___2==1)]
pheno$cd_surgery_colonic_colectomy_ileoanal_pouch_date[which(pheno$op_3___2==1)] <- pheno$surg_date_1[which(pheno$op_3___2==1)]
pheno$cd_surgery_colonic_colectomy_ileoanal_pouch_date[which(pheno$op_4___2==1)] <- pheno$surg_date_1[which(pheno$op_4___2==1)]
pheno$cd_surgery_colonic_colectomy_ileoanal_pouch_date[which(pheno$op_5___2==1)] <- pheno$surg_date_1[which(pheno$op_5___2==1)]
pheno$cd_surgery_colonic_colectomy_ileoanal_pouch_date[is.na(pheno$cd_surgery_colonic_colectomy_ileoanal_pouch)] <- NA
pheno$cd_surgery_colonic_colectomy_ileoanal_pouch_date <- as.double(pheno$cd_surgery_colonic_colectomy_ileoanal_pouch_date)

## Defunctioning ileostomy/colostomy
pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy <- rowSums(pheno[,c("op_1___3","op_2___3","op_3___3","op_4___3","op_5___3")],na.rm=T)
pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy[which(pheno$diagnosis!=1)]<- NA

pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy_date<- NA
pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy_date[which(pheno$op_1___3==1)] <- pheno$surg_date_1[which(pheno$op_1___3==1)]
pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy_date[which(pheno$op_2___3==1)] <- pheno$surg_date_1[which(pheno$op_2___3==1)]
pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy_date[which(pheno$op_3___3==1)] <- pheno$surg_date_1[which(pheno$op_3___3==1)]
pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy_date[which(pheno$op_4___3==1)] <- pheno$surg_date_1[which(pheno$op_4___3==1)]
pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy_date[which(pheno$op_5___3==1)] <- pheno$surg_date_1[which(pheno$op_5___3==1)]
pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy_date[is.na(pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy)] <- NA
pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy_date <- as.double(pheno$cd_surgery_colonic_defunctioning_ileostomy_colostomy_date)

## Partial colectomy
pheno$cd_surgery_colonic_partial_colectomy <- rowSums(pheno[,c("op_1___9","op_2___9","op_3___9","op_4___9","op_5___9")],na.rm=T)
pheno$cd_surgery_colonic_partial_colectomy[which(pheno$diagnosis!=1)]<- NA

### no L1
pheno$cd_surgery_colonic_partial_colectomy_noL1 <- pheno$cd_surgery_colonic_partial_colectomy
pheno$cd_surgery_colonic_partial_colectomy_noL1[which(pheno$location=="L1")]<- NA


pheno$cd_surgery_colonic_partial_colectomy_noL1_date<- NA
pheno$cd_surgery_colonic_partial_colectomy_noL1_date[which(pheno$op_1___3==1)] <- pheno$surg_date_1[which(pheno$op_1___3==1)]
pheno$cd_surgery_colonic_partial_colectomy_noL1_date[which(pheno$op_2___3==1)] <- pheno$surg_date_1[which(pheno$op_2___3==1)]
pheno$cd_surgery_colonic_partial_colectomy_noL1_date[which(pheno$op_3___3==1)] <- pheno$surg_date_1[which(pheno$op_3___3==1)]
pheno$cd_surgery_colonic_partial_colectomy_noL1_date[which(pheno$op_4___3==1)] <- pheno$surg_date_1[which(pheno$op_4___3==1)]
pheno$cd_surgery_colonic_partial_colectomy_noL1_date[which(pheno$op_5___3==1)] <- pheno$surg_date_1[which(pheno$op_5___3==1)]
pheno$cd_surgery_colonic_partial_colectomy_noL1_date[is.na(pheno$cd_surgery_colonic_partial_colectomy_noL1)] <- NA
pheno$cd_surgery_colonic_partial_colectomy_noL1_date <- as.double(pheno$cd_surgery_colonic_partial_colectomy_noL1_date)

## Proctectomy
pheno$cd_surgery_colonic_proctectomy <- rowSums(pheno[,c("op_1___10","op_2___10","op_3___10","op_4___10","op_5___10")],na.rm=T)
pheno$cd_surgery_colonic_proctectomy[which(pheno$diagnosis!=1)]<- NA

## Proctectomy, no pCD+
pheno$cd_surgery_colonic_proctectomy_noPerianal <- pheno$cd_surgery_colonic_proctectomy
pheno$cd_surgery_colonic_proctectomy_noPerianal[which(pheno$pCD==1 & pheno$cd_surgery_colonic_proctectomy>=1)] <- NA


#do not include Proctectomy and Partial colectomy with L1
pheno$cd_surgery_colonic <- rowSums(pheno[,c("cd_surgery_colonic_colectomy_ileostomy","cd_surgery_colonic_colectomy_ileoanal_pouch","cd_surgery_colonic_defunctioning_ileostomy_colostomy","cd_surgery_colonic_partial_colectomy_noL1")],na.rm=T)






# #reclassify CD location of L1 -> L3 is the patients is CD & L1 & cd_surgery_colonic > 0.
# pheno[which(pheno$location=="L1" & pheno$cd_surgery_colonic>0 & pheno$diagnosis==1),"location"] <- "L3"


#Surgery for small bowel disease
## Drainage of intra-abdominal abscess
pheno$cd_surgery_small_bowel_intra_abdomial_abscess <- rowSums(pheno[,c("op_1___4","op_2___4","op_3___4","op_4___4","op_5___4")],na.rm=T)
pheno$cd_surgery_small_bowel_intra_abdomial_abscess[which(pheno$diagnosis!=1)]<- NA


## Ileal / jejunal resection
pheno$cd_surgery_small_bowel_ileal_jejunal_resection <- rowSums(pheno[,c("op_1___6","op_2___6","op_3___6","op_4___6","op_5___6")],na.rm=T)
pheno$cd_surgery_small_bowel_ileal_jejunal_resection[which(pheno$diagnosis!=1)]<- NA


## Ileal / jejunal resection, no B1
pheno$cd_surgery_small_bowel_ileal_jejunal_resection_noB1 <- pheno$cd_surgery_small_bowel_ileal_jejunal_resection

pheno$cd_surgery_small_bowel_ileal_jejunal_resection_noB1[which(pheno$behaviour==1 & pheno$cd_surgery_small_bowel_ileal_jejunal_resection>=1)] <- NA


## Ileal / jejunal stricturoplasty
pheno$cd_surgery_small_bowel_ileal_jejunal_stricturoplasty<- rowSums(pheno[,c("op_1___7","op_2___7","op_3___7","op_4___7","op_5___7")],na.rm=T)
pheno$cd_surgery_small_bowel_ileal_jejunal_stricturoplasty[which(pheno$diagnosis!=1)]<- NA


## Ileo-caecal resection (Right hemicolectomy)
pheno$cd_surgery_small_bowel_ileal_caecal_resection<- rowSums(pheno[,c("op_1___8","op_2___8","op_3___8","op_4___8","op_5___8")],na.rm=T)
pheno$cd_surgery_small_bowel_ileal_caecal_resection[which(pheno$diagnosis!=1)]<- NA

## Ileo-caecal resection (Right hemicolectomy), no B1
pheno$cd_surgery_small_bowel_ileal_caecal_resection_noB1<- pheno$cd_surgery_small_bowel_ileal_caecal_resection

pheno$cd_surgery_small_bowel_ileal_caecal_resection_noB1[which(pheno$behaviour==1 & pheno$cd_surgery_small_bowel_ileal_caecal_resection>=1)] <- NA

## Strictureplasty
pheno$cd_surgery_small_bowel_strictureplasty<- rowSums(pheno[,c("op_1___11","op_2___11","op_3___11","op_4___11","op_5___11")],na.rm=T)
pheno$cd_surgery_small_bowel_strictureplasty[which(pheno$diagnosis!=1)]<- NA


pheno$cd_surgery_small_bowel <- rowSums(pheno[,c("cd_surgery_small_bowel_intra_abdomial_abscess","cd_surgery_small_bowel_ileal_jejunal_resection_noB1","cd_surgery_small_bowel_ileal_jejunal_stricturoplasty","cd_surgery_small_bowel_ileal_caecal_resection_noB1","cd_surgery_small_bowel_strictureplasty")],na.rm=T)


data <- pheno[which(pheno$diagnosis==1),c("ID_postqc_data","ID_postqc_data","behaviour_baseline","location","location1","PC1","PC2","PC3","PC4","age_diag","sex","smoker_type","IBD_PRS","CD_PRS","UC_PRS","CDvsUC_PRS","sub_certain_ibd","certain_diagnosis","follow_up","diag_year","SP0011_id_1","surgery_baseline","cd_surgery_small_bowel","cd_surgery_colonic")]

colnames(data)[3] <- "behaviour"

data <- data[!is.na(data$behaviour),]
data$bev <- as.double(data$behaviour>1)
data$bev3to1 <- as.double(data$behaviour==3)
data$bev3to1[which(data$behaviour==2)] <- NA
data$bev2to1 <- as.double(data$behaviour==2)
data$bev2to1[which(data$behaviour==3)] <- NA
data$bev3to12 <- as.double(data$behaviour==3)
data$bev2to13 <- as.double(data$behaviour==2)


colnames(data)[1:2] <- c("FID","IID")

data$timeSinceDiag_5yrs <- data$follow_up/5
data$age_diag_section <- NA
data$age_diag_section[which(data$age_diag<=16)] <- "A1"
data$age_diag_section[which(data$age_diag>16 & data$age_diag<=40)] <- "A2"
data$age_diag_section[which(data$age_diag>40)] <- "A3"

data$year_diag_section <- NA
data$year_diag_section[which(data$diag_year<2002)] <- "pre-2002"
data$year_diag_section[which(data$diag_year>=2002 & data$diag_year<=2011)] <- "2002-2011"
data$year_diag_section[which(data$diag_year>2011)] <- "post-2011"

write.table(data,file="~/CD_behaviour/pheno/ibdbr_pheno.txt",quote=F,row.names=F)


