library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)

setwd("~/CD_behaviour/")

pheno <- read.table("pheno/ibdbr_pheno.txt",header=T)


idmap <- fread("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/genotyped_data/QC/release_interval_ibdbioresource_Aug24/bioresource_hg19_noind_posstr_nodup_edited_id_map_file_with_batches_updated_version_20251128.tsv.gz",sep=",")

hla <- fread("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/HLA/Aug24/Broad_4D.csv")
hla <- hla[ hla$ancestry=="eur",]
hla <- as.data.frame(hla)
idx <- which(colnames(hla)=="HLA-4D-DRB1*01:03")
hla_sel <- hla[,c(1,idx)]
colnames(hla_sel) <- c("FID","DRB1_01_03")

hla_sel$disease <- "Control"
hla_sel$disease[ hla_sel$FID %in% idmap$ID_postqc_data]<- "IBD"



aa <- fread("/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/HLA/Aug24/SNP2HLA_AA.csv")
aa <- as.data.frame(aa)
idx <- which(colnames(aa)=="AA_DQA1_175_32610016_exon3_E")
aa_sel <- aa[,c(1,idx)]
colnames(aa_sel) <- c("FID","DQA1_Gln175Glu")

identical(hla_sel$FID,aa_sel$FID)
# [1] TRUE

hla_sel$DQA1_Gln175Glu <- aa_sel$DQA1_Gln175Glu

data <- merge(hla_sel,pheno,by="FID",all.x=T)

data <- data[which(data$disease=="Control" | !is.na(data$behaviour)),]
data[ data$disease=="IBD","disease"] <- NA
data[ which(data$bev==0 & data$location=="L2"),"disease"] <- "B1, Colonic CD"
data[ which(data$bev>0 & data$location=="L2"),"disease"] <- "B2/B3, Colonic CD"
data[ which(data$bev==0 & data$location=="L1"),"disease"] <- "B1, Ileal CD"
data[ which(data$bev>0 & data$location=="L1"),"disease"] <- "B2/B3, Ileal CD"
data[ which(data$bev==0 & data$location=="L3"),"disease"] <- "B1, Ileocolonic CD"
data[ which(data$bev>0 & data$location=="L3"),"disease"] <- "B2/B3, Ileocolonic CD"
data1 <- data[!is.na(data$disease),]
# data[ which(data$behaviour==3),"disease"] <- "B3"
# data[ which(data$behaviour==4),"disease"] <- "B2B3"



# dataplot_mean1 <- data %>%
#   dplyr::group_by(disease) %>%
#   dplyr::summarise(
#     mean_value = mean(DRB1_01_03, na.rm = TRUE)/2,
#     n = n(),                       # optional: sample size
#     .groups = "drop"
#   )


dataplot_mean2 <- data1 %>%
  dplyr::group_by(disease) %>%
  dplyr::summarise(
    mean_value = mean(DQA1_Gln175Glu, na.rm = TRUE)/2,
    n = n(),                       # optional: sample size
    .groups = "drop"
  )


set2_cols <- brewer.pal(8, "Set2")

# dataplot_mean1$variant <- "HLA-DRB1*01:03 frequency"
dataplot_mean2$variant <- "HLA-DQA1:p.Gln175Glu frequency"

# dataplot_mean <- rbind(dataplot_mean1,dataplot_mean2)

dataplot_mean2$disease <- factor(dataplot_mean2$disease,levels=c("Control","B1, Colonic CD","B2/B3, Colonic CD","B1, Ileocolonic CD","B2/B3, Ileocolonic CD","B1, Ileal CD","B2/B3, Ileal CD"))

set2_cols <- brewer.pal(8, "Set2")

p <- ggplot(data=dataplot_mean2,aes(disease,mean_value))+facet_wrap(~variant,ncol=2,scale="free_x")+
  geom_col(width=0.4,fill=set2_cols[2])+
  theme_bw()+
  theme_classic()+
  theme(axis.text = element_text(size=14,face="bold"),
        axis.title = element_blank(),
        strip.text = element_text(size=14,face="bold"),
        legend.position = "none")+
  geom_text(
    aes(label = percent(mean_value, accuracy = 0.1)),
    position = position_dodge2(width = 0.8, preserve = "single"),hjust = -0.1,size = 6)+
  scale_y_continuous(labels = percent_format(accuracy = 1),expand = expansion(mult = c(0, 0.2)))+coord_flip()


ggsave(file="~/CD_behaviour/figures/freq.png",dpi=300,width=8,height=8)

ggsave(file="~/CD_behaviour/figures/freq.pdf",width=8,height=8)
