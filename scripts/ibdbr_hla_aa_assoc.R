library(data.table)

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


idxs <- grep("AA_|HLA",colnames(data))

smy <- c()
maf <- c()
for(i in idxs){
  coefs <- summary(glm(data[,"bev"] ~ data[,i] + data$PC1 + data$PC2 + data$PC3 + data$PC4 + data$CHR6 ,family="binomial"))$coefficients
  smy <-rbind(smy,coefs[2,])
  maf <- c(maf,mean(data[,i])/2)
}

smy <- as.data.frame(smy)
colnames(smy) <- c("beta","se","z","pval")
smy$maf <- maf
smy$ID <- colnames(data)[idxs]

smy <- smy[smy$maf > 0.001,]

write.csv(smy,file="result/ibdbr_hla_aa_4PC.csv",quote=F,row.names=F)



smy <- c()
maf <- c()
for(i in idxs){
  coefs <- summary(glm(data[,"bev"] ~ data[,i] + data$PC1 + data$PC2 + data$PC3 + data$PC4 + data$CHR6 + data$location,family="binomial"))$coefficients
  smy <-rbind(smy,coefs[2,])
  maf <- c(maf,mean(data[,i])/2)
}

smy <- as.data.frame(smy)
colnames(smy) <- c("beta","se","z","pval")
smy$maf <- maf
smy$ID <- colnames(data)[idxs]

smy <- smy[smy$maf > 0.001,]

write.csv(smy,file="result/ibdbr_hla_aa_4PC_loc.csv",quote=F,row.names=F)



idx1 <- which(colnames(data)=="AA_DQA1_175_32610016_exon3_E")

smy <- c()
maf <- c()
for(i in idxs){
  coefs <- summary(glm(data[,"bev"] ~ data[,i] + data$PC1 + data$PC2 + data$PC3 + data$PC4 + data$CHR6 + data$location + data[,idx1],family="binomial"))$coefficients
  smy <-rbind(smy,coefs[2,])
  maf <- c(maf,mean(data[,i])/2)
}

smy <- as.data.frame(smy)
colnames(smy) <- c("beta","se","z","pval")
smy$maf <- maf
smy$ID <- colnames(data)[idxs]

smy <- smy[smy$maf > 0.001,]

write.csv(smy,file="result/ibdbr_hla_aa_4PC_loc_DQA1_175.csv",quote=F,row.names=F)

