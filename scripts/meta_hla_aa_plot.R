library(ggplot2)

setwd("~/CD_behaviour/")

hla <- read.csv("result/hla_aa_ibdbr_ukibdgc_4PC_meta.csv")
hla4D <- hla[ substr(hla$ID,0,3)=="HLA",]

geneid <- unlist(strsplit(hla4D$ID,split="\\*|-"))

hla4D$gene <- geneid[seq(3,length(geneid),4)]
hla4D$allele <- geneid[seq(4,length(geneid),4)]

hla4D$gene_allele <- paste0(hla4D$gene,"*",hla4D$allele)

hla4D$label <- hla4D$gene_allele
hla4D$label[ hla4D$p_meta>=5e-8] <- NA

hla4D$gene <- factor(hla4D$gene, levels=c("A","C","B","DRB1","DQA1","DQB1","DPB1"))

aa <- hla[substr(hla$ID,0,2)=="AA",]
geneid <- sapply(strsplit(aa$ID,split="_"), `[`, 2)

aa$gene <- geneid
aa$allele <- aa$ID
aa$gene_allele <- aa$ID

aa$label <- NA
aa$label[aa$ID=="AA_DQA1_175_32610016_exon3_E"] <- "DQA1:p.Gln175Glu"

aa <- aa[aa$gene!="DPA1",]

aa$gene <- factor(aa$gene, levels=c("A","C","B","DRB1","DQA1","DQB1","DPB1"))


dataplot <- rbind(hla4D,aa)

p1 <- ggplot(data=dataplot,aes(ID,-log10(p_meta),label=label,col=gene))+geom_point(size=3,alpha = 0.7)+facet_wrap(~gene,ncol=7,scale="free_x",strip.position = "bottom") +
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.y=element_text(size = 12,face = "bold"),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size = 12,face = "bold"),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        strip.placement = "outside",
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=12,face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=12,face="bold"))+
  scale_fill_brewer(palette = "Set2")+
  geom_hline(yintercept = -log10(5e-8),lty=1,color="violetred")+
  geom_hline(yintercept = -log10(5e-6),lty=2)+
  geom_text_repel(nudge_y = ifelse(-log10(hla4D$p_meta) < 2, 2 + log10(hla4D$p_meta) , 0.5),col="black")+
  ylab("-log10(P)")


ggsave(file="~/CD_behaviour/figures/meta_HLA_4D.png",dpi=300,height=4.5,width=12)



hla <- read.csv("result/hla_aa_ibdbr_ukibdgc_4PC_loc_meta.csv")
hla4D <- hla[ substr(hla$ID,0,3)=="HLA",]

geneid <- unlist(strsplit(hla4D$ID,split="\\*|-"))

hla4D$gene <- geneid[seq(3,length(geneid),4)]
hla4D$allele <- geneid[seq(4,length(geneid),4)]

hla4D$gene_allele <- paste0(hla4D$gene,"*",hla4D$allele)

hla4D$label <- hla4D$gene_allele
hla4D$label[ hla4D$p_meta>=5e-8] <- NA

hla4D$gene <- factor(hla4D$gene, levels=c("A","C","B","DRB1","DQA1","DQB1","DPB1"))


aa <- hla[substr(hla$ID,0,2)=="AA",]
geneid <- sapply(strsplit(aa$ID,split="_"), `[`, 2)

aa$gene <- geneid
aa$allele <- aa$ID
aa$gene_allele <- aa$ID

aa$label <- NA
aa$label[aa$ID=="AA_DQA1_175_32610016_exon3_E"] <- "DQA1:p.Gln175Glu"

aa <- aa[aa$gene!="DPA1",]

aa$gene <- factor(aa$gene, levels=c("A","C","B","DRB1","DQA1","DQB1","DPB1"))

dataplot <- rbind(hla4D,aa)

p1 <- ggplot(data=dataplot,aes(ID,-log10(p_meta),label=label,col=gene))+geom_point(size=3,alpha = 0.7)+facet_wrap(~gene,ncol=7,scale="free_x",strip.position = "bottom") +
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.y=element_text(size = 12,face = "bold"),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size = 12,face = "bold"),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        strip.placement = "outside",
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=12,face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=12,face="bold"))+
  scale_fill_brewer(palette = "Set2")+
  geom_hline(yintercept = -log10(5e-8),lty=1,color="violetred")+
  geom_hline(yintercept = -log10(5e-6),lty=2)+
  geom_text_repel(nudge_y = ifelse(-log10(hla4D$p_meta) < 2, 2 + log10(hla4D$p_meta) , 0.5),col="black")+
  ylab("-log10(P)")


ggsave(file="~/CD_behaviour/figures/meta_HLA_4D_loc.png",dpi=300,height=4.5,width=12)


hla <- read.csv("result/hla_aa_ibdbr_ukibdgc_4PC_loc_DQA1_175_meta.csv")
hla4D <- hla[ substr(hla$ID,0,3)=="HLA",]

geneid <- unlist(strsplit(hla4D$ID,split="\\*|-"))

hla4D$gene <- geneid[seq(3,length(geneid),4)]
hla4D$allele <- geneid[seq(4,length(geneid),4)]

hla4D$gene_allele <- paste0(hla4D$gene,"*",hla4D$allele)

hla4D$label <- hla4D$gene_allele
hla4D$label[ hla4D$p_meta>=5e-8] <- NA

hla4D$gene <- factor(hla4D$gene, levels=c("A","C","B","DRB1","DQA1","DQB1","DPB1"))

aa <- hla[substr(hla$ID,0,2)=="AA",]
geneid <- sapply(strsplit(aa$ID,split="_"), `[`, 2)

aa$gene <- geneid
aa$allele <- aa$ID
aa$gene_allele <- aa$ID

aa$label <- NA
aa$label[aa$ID=="AA_DQA1_175_32610016_exon3_E"] <- "DQA1:p.Gln175Glu"

aa <- aa[aa$gene!="DPA1",]

aa$gene <- factor(aa$gene, levels=c("A","C","B","DRB1","DQA1","DQB1","DPB1"))

dataplot <- rbind(hla4D,aa)

p1 <- ggplot(data=dataplot,aes(ID,-log10(p_meta),label=label,col=gene))+geom_point(size=3,alpha = 0.7)+facet_wrap(~gene,ncol=7,scale="free_x",strip.position = "bottom") +
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.y=element_text(size = 12,face = "bold"),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size = 12,face = "bold"),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        strip.placement = "outside",
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=12,face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=12,face="bold"))+
  scale_fill_brewer(palette = "Set2")+
  geom_hline(yintercept = -log10(5e-8),lty=1,color="violetred")+
  geom_hline(yintercept = -log10(5e-6),lty=2)+
  geom_text_repel(nudge_y = ifelse(-log10(hla4D$p_meta) < 2, 2 + log10(hla4D$p_meta) , 0.5),col="black")+
  ylab("-log10(P)")


ggsave(file="~/CD_behaviour/figures/meta_HLA_4D_loc_DQA1_175.png",dpi=300,height=4.5,width=12)
