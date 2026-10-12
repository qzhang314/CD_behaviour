library(data.table)

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







setwd("~/CD_behaviour/")

prs <- read.table("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step1/pc_only_1.loco",header=T)

prs <- t(prs)
prs <- prs[-1,]
prs <- as.data.frame(prs)
prs$ID <- rownames(prs)
colnames(prs)[1:22] <- paste0("CHR",1:22)


aa <- fread("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/HLA/Aug24/SNP2HLA_AA.csv")
aa <- as.data.frame(aa)

hla <- fread("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/HLA/Aug24/Broad_4D.csv")
hla <- as.data.frame(hla)


pheno <- read.table("pheno/ibdbr_pheno.txt",header=T)

data <- merge(pheno,hla,by.x="FID",by.y="sample_id")
data <- merge(data,aa,by.x="FID",by.y="sample_id")

data$ID <- paste0(data$FID,"_",data$FID)
data <- merge(data,prs,by="ID")


smy <- c()

coefs <- summary(glm(bev ~ AA_DQA1_175_32610016_exon3_E + PC1 + PC2 + PC3 + PC4 + CHR6, data=data[data$location=="L1",],family="binomial"))$coefficients
smy <-rbind(smy,coefs[2,])

coefs <- summary(glm(bev ~ AA_DQA1_175_32610016_exon3_E + PC1 + PC2 + PC3 + PC4 + CHR6, data=data[data$location=="L2",],family="binomial"))$coefficients
smy <-rbind(smy,coefs[2,])

coefs <- summary(glm(bev ~ AA_DQA1_175_32610016_exon3_E + PC1 + PC2 + PC3 + PC4 + CHR6, data=data[data$location=="L3",],family="binomial"))$coefficients
smy <-rbind(smy,coefs[2,])

smy <- as.data.frame(smy)
colnames(smy) <- c("beta","se","z","pval")

smy_ibdbr <- smy



aa <- fread("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/HLA/HLA_Qian/result/SNP2HLA_AA.csv")
aa <- as.data.frame(aa)

hla <- fread("/lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/HLA/HLA_Qian/result/Broad_4D.csv")
hla <- as.data.frame(hla)


pheno <- read.table("pheno/ukibdgc_pheno.txt",header=T)

data <- merge(pheno,hla,by.x="FID",by.y="sample_id")
data <- merge(data,aa,by.x="FID",by.y="sample_id")

smy <- c()

coefs <- summary(glm(bev ~ AA_DQA1_175_32610016_exon3_E + PC1 + PC2 + PC3 + PC4  + batch ,data=data[data$location=="L1",],family="binomial"))$coefficients
smy <-rbind(smy,coefs[2,])

coefs <- summary(glm(bev ~ AA_DQA1_175_32610016_exon3_E + PC1 + PC2 + PC3 + PC4  + batch ,data=data[data$location=="L2",],family="binomial"))$coefficients
smy <-rbind(smy,coefs[2,])

coefs <- summary(glm(bev ~ AA_DQA1_175_32610016_exon3_E + PC1 + PC2 + PC3 + PC4  + batch ,data=data[data$location=="L3",],family="binomial"))$coefficients
smy <-rbind(smy,coefs[2,])

smy <- as.data.frame(smy)
colnames(smy) <- c("beta","se","z","pval")

smy_ukibdgc <- smy

smy_ukibdgc$location <- c("L1","L2","L3")
smy_ibdbr$location <- c("L1","L2","L3")

ibdbr_ukibdgc <- merge(smy_ibdbr,smy_ukibdgc,by="location",suffix=c("_ibdbr","_ukibdgc"))

ibdbr_ukibdgc_meta <- mapply(combine_two_fixed, ibdbr_ukibdgc$beta_ibdbr, ibdbr_ukibdgc$se_ibdbr, ibdbr_ukibdgc$beta_ukibdgc, ibdbr_ukibdgc$se_ukibdgc)
ibdbr_ukibdgc_meta <- as.data.frame(t(ibdbr_ukibdgc_meta))

ibdbr_ukibdgc$p_meta <- ibdbr_ukibdgc_meta$p_meta
ibdbr_ukibdgc$beta_meta <- ibdbr_ukibdgc_meta$beta_meta
ibdbr_ukibdgc$se_meta <- ibdbr_ukibdgc_meta$se_meta
ibdbr_ukibdgc$p_heterogeneity <- ibdbr_ukibdgc_meta$p_heterogeneity


