setwd("~/CD_behaviour/")

combine_two_fixed <- function(y1, se1, y2, se2) {
  v1 <- se1^2
  v2 <- se2^2

  
  # denominator for GLS weights
  den <- v1 + v2 

  # fixed-effect pooled estimate
  theta <- (v2 * y1 + v1 * y2) / den
  var_theta <- (v1 * v2) / den
  se_theta <- sqrt(var_theta)
  p_meta <- pchisq((theta/se_theta)^2, df=1, lower.tail = FALSE)
  
  
  # --- Heterogeneity (Q-test) ---
  Q <- (y1 - y2)^2 / (v1 + v2)
  df <- 1
  p_hetero <- pchisq(Q, df, lower.tail = FALSE)
  I2 <- max(0, 100 * (Q - df) / Q)
  
  c(
    beta1=y1,
    beta2=y2,
    se1=se1,
    se2=se2,
    p_meta = p_meta,
    beta_meta = theta,
    se_meta = se_theta,
    Q = Q,
    df = df,
    p_heterogeneity = p_hetero,
    I2 = I2
  )
}

ibdbr <- read.csv("~/CD_behaviour/result/ibdbr_hla_aa_4PC.csv")
ukibdgc <- read.csv("~/CD_behaviour/result/ukibdgc_hla_aa_4PC.csv")

ibdbr_ukibdgc <- merge(ibdbr,ukibdgc,by="ID",suffix=c("_ibdbr","_ukibdgc"))

ibdbr_ukibdgc_meta <- mapply(combine_two_fixed, ibdbr_ukibdgc$beta_ibdbr, ibdbr_ukibdgc$se_ibdbr, ibdbr_ukibdgc$beta_ukibdgc, ibdbr_ukibdgc$se_ukibdgc)
ibdbr_ukibdgc_meta <- as.data.frame(t(ibdbr_ukibdgc_meta))

ibdbr_ukibdgc$p_meta <- ibdbr_ukibdgc_meta$p_meta
ibdbr_ukibdgc$beta_meta <- ibdbr_ukibdgc_meta$beta_meta
ibdbr_ukibdgc$se_meta <- ibdbr_ukibdgc_meta$se_meta
ibdbr_ukibdgc$p_heterogeneity <- ibdbr_ukibdgc_meta$p_heterogeneity

write.csv(ibdbr_ukibdgc,"~/CD_behaviour/result/hla_aa_ibdbr_ukibdgc_4PC_meta.csv",quote=F,row.names=F)




ibdbr <- read.csv("~/CD_behaviour/result/ibdbr_hla_aa_4PC_loc.csv")
ukibdgc <- read.csv("~/CD_behaviour/result/ukibdgc_hla_aa_4PC_loc.csv")

ibdbr_ukibdgc <- merge(ibdbr,ukibdgc,by="ID",suffix=c("_ibdbr","_ukibdgc"))

ibdbr_ukibdgc_meta <- mapply(combine_two_fixed, ibdbr_ukibdgc$beta_ibdbr, ibdbr_ukibdgc$se_ibdbr, ibdbr_ukibdgc$beta_ukibdgc, ibdbr_ukibdgc$se_ukibdgc)
ibdbr_ukibdgc_meta <- as.data.frame(t(ibdbr_ukibdgc_meta))

ibdbr_ukibdgc$p_meta <- ibdbr_ukibdgc_meta$p_meta
ibdbr_ukibdgc$beta_meta <- ibdbr_ukibdgc_meta$beta_meta
ibdbr_ukibdgc$se_meta <- ibdbr_ukibdgc_meta$se_meta
ibdbr_ukibdgc$p_heterogeneity <- ibdbr_ukibdgc_meta$p_heterogeneity

write.csv(ibdbr_ukibdgc,"~/CD_behaviour/result/hla_aa_ibdbr_ukibdgc_4PC_loc_meta.csv",quote=F,row.names=F)


ibdbr <- read.csv("~/CD_behaviour/result/ibdbr_hla_aa_4PC_loc_DQA1_175.csv")
ukibdgc <- read.csv("~/CD_behaviour/result/ukibdgc_hla_aa_4PC_loc_DQA1_175.csv")

ibdbr_ukibdgc <- merge(ibdbr,ukibdgc,by="ID",suffix=c("_ibdbr","_ukibdgc"))

ibdbr_ukibdgc_meta <- mapply(combine_two_fixed, ibdbr_ukibdgc$beta_ibdbr, ibdbr_ukibdgc$se_ibdbr, ibdbr_ukibdgc$beta_ukibdgc, ibdbr_ukibdgc$se_ukibdgc)
ibdbr_ukibdgc_meta <- as.data.frame(t(ibdbr_ukibdgc_meta))

ibdbr_ukibdgc$p_meta <- ibdbr_ukibdgc_meta$p_meta
ibdbr_ukibdgc$beta_meta <- ibdbr_ukibdgc_meta$beta_meta
ibdbr_ukibdgc$se_meta <- ibdbr_ukibdgc_meta$se_meta
ibdbr_ukibdgc$p_heterogeneity <- ibdbr_ukibdgc_meta$p_heterogeneity

write.csv(ibdbr_ukibdgc,"~/CD_behaviour/result/hla_aa_ibdbr_ukibdgc_4PC_loc_DQA1_175_meta.csv",quote=F,row.names=F)


