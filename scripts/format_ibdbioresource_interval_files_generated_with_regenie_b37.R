# MEM=10000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G ukbb_team152 /software/bin/R-3.6.1

library(data.table)
library(ggplot2)

path<-"/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/"
path_gwas<-"/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step2/"

args <- commandArgs(T)
cov_type <- as.character(args[1])

chr <- as.double(Sys.getenv("LSB_JOBINDEX"))

if(chr==23){
        chridx <- "X"
}else{
        chridx <- chr
}

print(chridx)
print(cov_type)

info<-fread(paste(path,"imputed_data/release_interval_ibdbioresource_Aug24/combined_hrc_uk10k/",chridx,"_combined_hrc_uk10k_eur_with_consensus_ids_imputation_info_hwe_ctr_marker_names_b37_b38_imputation_info_dbsnp154.tsv.gz",sep=""),head=T)

info<-as.data.frame(info)

dups<-info$id_b38[duplicated(info$id_b38)]
dups<-dups[!duplicated(dups)]

for (i in 1:length(dups)) {
  tmp<-info[which(info$id_b38==dups[i]),]
  tmp<-tmp[which(tmp$impute_info_cases==min(tmp$impute_info_cases)),]

  if(i==1) {
    to_exclude<-as.character(tmp$id_b37)
  } else {
    to_exclude<-c(to_exclude,as.character(tmp$id_b37))
  }
}

table(length(to_exclude)==length(dups))
info<-info[which(!info$id_b37 %in% to_exclude),]

  
  rg<-fread(paste0(path_gwas,cov_type,"/chr",chridx,"_bev.regenie"),head=T)
    
  print(paste("N variants analysed: ",nrow(rg),sep=""))
  
#  rg$ID<-paste(rg$CHROM,rg$GENPOS,rg$ALLELE0,rg$ALLELE1,sep="_")
  
  print(paste("N variants liftover: ",nrow(rg[which(rg$ID %in% info$id_b37),])))
  
  rg<-merge(rg,info,by.x="ID",by.y="id_b37",all.x=T)
  
  rg$BETA_ed<-NA
  rg$A1FREQ_ed<-NA
  rg$A1FREQ_CASES_ed<-NA
  rg$A1FREQ_CONTROLS_ed<-NA
  
  # VARIANTS WITH SAME REF AND ALT ALLELES - same effect
  # print(nrow(rg[which( (rg$ALLELE1==rg$Alt_b38) & (rg$ALLELE0==rg$Ref_b38)),])/nrow(rg)) [1] 0.9997235 of variants
  
  rg$BETA_ed[which( (rg$ALLELE1==rg$Alt_b38) & (rg$ALLELE0==rg$Ref_b38))]<-rg$BETA[which( (rg$ALLELE1==rg$Alt_b38) & (rg$ALLELE0==rg$Ref_b38))]
  rg$A1FREQ_ed[which( (rg$ALLELE1==rg$Alt_b38) & (rg$ALLELE0==rg$Ref_b38))]<-rg$A1FREQ[which( (rg$ALLELE1==rg$Alt_b38) & (rg$ALLELE0==rg$Ref_b38))]
  rg$A1FREQ_CASES_ed[which( (rg$ALLELE1==rg$Alt_b38) & (rg$ALLELE0==rg$Ref_b38))]<-rg$A1FREQ_CASES[which( (rg$ALLELE1==rg$Alt_b38) & (rg$ALLELE0==rg$Ref_b38))]
  rg$A1FREQ_CONTROLS_ed[which( (rg$ALLELE1==rg$Alt_b38) & (rg$ALLELE0==rg$Ref_b38))]<-rg$A1FREQ_CONTROLS[which( (rg$ALLELE1==rg$Alt_b38) & (rg$ALLELE0==rg$Ref_b38))]
  
  
  # VARIANTS WHERE REF AND ALT HAVE BEEN FLIPPED - effect*-1:
  # print(nrow(rg[which( (rg$ALLELE0==rg$Alt_b38) & (rg$ALLELE1==rg$Ref_b38)),])/nrow(rg)) # 0.0002355134 of variants
  
  rg$BETA_ed[which( (rg$ALLELE0==rg$Alt_b38) & (rg$ALLELE1==rg$Ref_b38))]<-rg$BETA[which( (rg$ALLELE0==rg$Alt_b38) & (rg$ALLELE1==rg$Ref_b38))]*(-1)
  rg$A1FREQ_ed[which( (rg$ALLELE0==rg$Alt_b38) & (rg$ALLELE1==rg$Ref_b38))]<-1-rg$A1FREQ[which( (rg$ALLELE0==rg$Alt_b38) & (rg$ALLELE1==rg$Ref_b38))]
  rg$A1FREQ_CASES_ed[which( (rg$ALLELE0==rg$Alt_b38) & (rg$ALLELE1==rg$Ref_b38))]<-1-rg$A1FREQ_CASES[which( (rg$ALLELE0==rg$Alt_b38) & (rg$ALLELE1==rg$Ref_b38))]
  rg$A1FREQ_CONTROLS_ed[which( (rg$ALLELE0==rg$Alt_b38) & (rg$ALLELE1==rg$Ref_b38))]<-1-rg$A1FREQ_CONTROLS[which( (rg$ALLELE0==rg$Alt_b38) & (rg$ALLELE1==rg$Ref_b38))]
  
  
  # VARIANTS WHERE REF AND ALT HAVE BEEN SWAPPED -  same effect:
  
  print(nrow(rg[which( is.na(rg$BETA_ed) & !(is.na(rg$position_b38))),])/nrow(rg)) # 5.119856e-06 of variants
  
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]
  
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]<-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]<-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]<-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]<-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]<-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]<-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]<-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]<-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]
  
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]
  
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="T") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="G") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="T") & (rg$Ref_b38=="A") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="C") & (rg$ALLELE0=="A") & (rg$Ref_b38=="T") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="C") & (rg$Ref_b38=="G") & !(is.na(rg$position_b38)))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]<-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="A") & (rg$ALLELE0=="G") & (rg$Ref_b38=="C") & !(is.na(rg$position_b38)))]
  
  
  # VARIANTS WHERE REF AND ALT HAVE BEEN SWAPPED AND FLIPPED - effect*-1:
  print(nrow(rg[which( is.na(rg$BETA_ed) & !(is.na(rg$position_b38))),])/nrow(rg)) # [1] 3.89109e-05
  
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="T"))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="T"))]*(-1)
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="T"))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="T"))]*(-1)
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="G"))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="G"))]*(-1)
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="G"))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="G"))]*(-1)
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="C"))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="C"))]*(-1)
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="C"))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="C"))]*(-1)
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="A"))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="A"))]*(-1)
  rg$BETA_ed[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="A"))]<-rg$BETA[which( is.na(rg$BETA_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="A"))]*(-1)
  
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="T"))]<-1-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="T"))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="T"))]<-1-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="T"))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="G"))]<-1-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="G"))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="G"))]<-1-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="G"))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="C"))]<-1-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="C"))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="C"))]<-1-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="C"))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="A"))]<-1-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="A"))]
  rg$A1FREQ_ed[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="A"))]<-1-rg$A1FREQ[which( is.na(rg$A1FREQ_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="A"))]
  
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="T"))]<-1-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="T"))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="T"))]<-1-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="T"))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="G"))]<-1-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="G"))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="G"))]<-1-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="G"))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="C"))]<-1-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="C"))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="C"))]<-1-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="C"))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="A"))]<-1-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="A"))]
  rg$A1FREQ_CASES_ed[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="A"))]<-1-rg$A1FREQ_CASES[which( is.na(rg$A1FREQ_CASES_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="A"))]
  
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="T"))]<-1-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="T"))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="T"))]<-1-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="A") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="T"))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="G"))]<-1-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="G"))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="G"))]<-1-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="C") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="G"))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="C"))]<-1-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="A") & (rg$ALLELE0=="T") & (rg$Ref_b38=="C"))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="C"))]<-1-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="G") & (rg$Alt_b38=="T") & (rg$ALLELE0=="A") & (rg$Ref_b38=="C"))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="A"))]<-1-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="G") & (rg$ALLELE0=="C") & (rg$Ref_b38=="A"))]
  rg$A1FREQ_CONTROLS_ed[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="A"))]<-1-rg$A1FREQ_CONTROLS[which( is.na(rg$A1FREQ_CONTROLS_ed) & (rg$ALLELE1=="T") & (rg$Alt_b38=="C") & (rg$ALLELE0=="G") & (rg$Ref_b38=="A"))]
  
  
  # OTHER - A/T C/G variants:
  print(nrow(rg[which( is.na(rg$BETA_ed) & !(is.na(rg$position_b38))),])/nrow(rg))  # 0 of variants
  
  flip<-rg[which( ((rg$A1FREQ_CONTROLS_ed<0.45 & rg$AF_gnomad_nfe>0.55) | (rg$A1FREQ_CONTROLS_ed>0.55 & rg$AF_gnomad_nfe<0.45)) &
                     ( (rg$Alt_b38=="G" & rg$Ref_b38=="C") | (rg$Alt_b38=="C" & rg$Ref_b38=="G") | 
                         (rg$Alt_b38=="A" & rg$Ref_b38=="T") | (rg$Alt_b38=="T" & rg$Ref_b38=="A") ) ),]
  
  # CHANGE BETA, AND FREQUENCIES FOR THOSE 
  if(nrow(flip)>0) {
    rg$BETA_ed[which(rg$id_b38 %in% flip$id_b38)]<-rg$BETA_ed[which(rg$id_b38 %in% flip$id_b38)]*(-1)
    rg$A1FREQ_ed[which(rg$id_b38 %in% flip$id_b38)]<-1-rg$A1FREQ_ed[which(rg$id_b38 %in% flip$id_b38)]
    rg$A1FREQ_CASES_ed[which(rg$id_b38 %in% flip$id_b38)]<-1-rg$A1FREQ_CASES_ed[which(rg$id_b38 %in% flip$id_b38)]
    rg$A1FREQ_CONTROLS_ed[which(rg$id_b38 %in% flip$id_b38)]<-1-rg$A1FREQ_CONTROLS_ed[which(rg$id_b38 %in% flip$id_b38)]
  }
  
  
  # EXCLUDE VARIANTS NOT LIFTOVER
  print(paste("N variants not liftover:",nrow(rg[which(is.na(rg$id_b38)),])))
  rg<-rg[which(!is.na(rg$id_b38)),]
  
  # EXCLUDE VARIANTS WITH LARGE DISCREPANCIES IN FREQ BETWEEN HRC AND UK10K
  print(paste("N variants low concond frq hrc uk10k:",nrow(rg[which(rg$low_conconrdance_freq_hrc_uk10k==1),])))
  rg<-rg[which(is.na(rg$low_conconrdance_freq_hrc_uk10k)),]

  # EXCLUDE VARIANTS BY HWE:
  print(paste("N variants hwe:",nrow(rg[which(rg$HW_lrt_p_value_ctr<1E-4),])))

  if(chr==23) {
    rg<-rg[which( (rg$HW_exact_p_value_all>1E-4) | is.na(rg$HW_exact_p_value_all)),]
  } else {
    rg<-rg[which(rg$HW_lrt_p_value_ctr>1E-4),]
  }

  
  # EXCLUDE VARIANTS BY ALLELE FREQ DISCREPANCY:
  print(paste("N variants freq:",nrow(rg[which( rg$unexpected_a1freq_discrepancy_cd_ibdbr_vs_allhce==1 |
                                                  rg$unexpected_a1freq_discrepancy_uc_ibdbr_vs_allhce==1 | 
                                                  rg$unexpected_a1freq_discrepancy_ctr_ibdbr_vs_allhce==1),])))
  
  rg<-rg[which( (is.na(rg$unexpected_a1freq_discrepancy_cd_ibdbr_vs_allhce) | rg$unexpected_a1freq_discrepancy_cd_ibdbr_vs_allhce==0) & 
               (is.na(rg$unexpected_a1freq_discrepancy_uc_ibdbr_vs_allhce) | rg$unexpected_a1freq_discrepancy_uc_ibdbr_vs_allhce==0) &
               (is.na(rg$unexpected_a1freq_discrepancy_ctr_ibdbr_vs_allhce) | rg$unexpected_a1freq_discrepancy_ctr_ibdbr_vs_allhce==0)),]

  
  # EXCLUDE VARIANTS BY UNEXPECTED FREQ DISCREPANCY BETWEEN CASES AND CONTROLS:
  print(paste("N variants freq CA/CTR:",nrow(rg[which( rg$unexpected_a1freq_discrepancy_cd_ctr_ibdbr==1 |
                                                  rg$unexpected_a1freq_discrepancy_uc_ctr_ibdbr==1),])))
  
  rg<-rg[which( (is.na(rg$unexpected_a1freq_discrepancy_cd_ctr_ibdbr) | rg$unexpected_a1freq_discrepancy_cd_ctr_ibdbr==0) & 
                  (is.na(rg$unexpected_a1freq_discrepancy_uc_ctr_ibdbr) | rg$unexpected_a1freq_discrepancy_uc_ctr_ibdbr==0)),]
  
  
  rg<-rg[order(rg$position_b38,decreasing=F),]
  rg$ID<-gsub("_",":",rg$id_b38)
  rg$ID<-paste("chr",rg$ID,sep="")
  
  
  rg<-rg[,c("CHROM","position_b38","ID","Ref_b38","Alt_b38","A1FREQ_ed","A1FREQ_CASES_ed","A1FREQ_CONTROLS_ed",
            "INFO","N","N_CASES","N_CONTROLS","TEST",
              "BETA_ed","SE","CHISQ","LOG10P","EXTRA","dbsnp154")]
  colnames(rg)<-c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","A1FREQ_CASES","A1FREQ_CONTROLS","INFO","N","N_CASES","N_CONTROLS","TEST",
                   "BETA","SE","CHISQ","LOG10P","EXTRA","rsID")
  
  print(paste("N variants for meta-analysis: ",nrow(rg),sep=""))
  
  if(chr==23) {
    rg$CHROM<-as.character(rg$CHROM)
    rg$CHROM<-"X"
    rg$ID<-gsub("chr23","chrX",rg$ID)
  }
  

write.table(rg,paste0(path_gwas,"/../QC/",cov_type,"/chr",chr,"_bev.regenie"),
           col.names=T,row.names=F,quote=F,na = "NA")
  
q("no")



