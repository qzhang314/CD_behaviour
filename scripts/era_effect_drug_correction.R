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

drugTimeSinceDiagAll <- readRDS("result/drugTimeSinceDiagAll.rds")


dataCD$year_diag_section <- factor(dataCD$year_diag_section,levels = c("pre-2002","2002-2011","post-2011"))

# explore the era effect when condictioning on anti-TNF prescription (whether or not), thiopurine prescription (whether or not) and methotreate prescription (whether or not)
era_eff <- summary(glm(bev ~  timeSinceDiag_5yrs  + as.factor(smoker_type) + year_diag_section + location + location1 + as.double(is.na(drugTimeSinceDiagAll[,1])&is.na(drugTimeSinceDiagAll[,2])) + as.double(is.na(drugTimeSinceDiagAll[,3])) + as.double(is.na(drugTimeSinceDiagAll[,4]))+ as.double(delay_years >= 1.5),data=dataCD,family = "binomial"))

era_eff_coef <- era_eff$coefficients[5:6,]
era_eff_coef <- as.data.frame(era_eff_coef)

colnames(era_eff_coef)<-c("beta","se","z","pval")
era_eff_coef$variable<-c("2002-2011 v.s. pre-2002","post-2011 v.s. pre-2002")
era_eff_coef$or <- exp(era_eff_coef$beta)
era_eff_coef$ci <- paste0(round(exp(era_eff_coef$beta - era_eff_coef$se*1.96),2),"-",round(exp(era_eff_coef$beta + era_eff_coef$se*1.96),2))

write.csv(era_eff_coef,file="result/era_effect_drug_correction.csv",quote=F,row.names=F)
