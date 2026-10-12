pheno=$1
cov_type=$2

if [[ ${cov_type} == "pc_only" ]]
then
	regenie \
		--step 1 \
		--bed /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/genotyped_data/QC/release_interval_ibdbioresource_Aug24/interval_ibdbioresource_postqc_hg19_nomonom_RefAlt_flip2_pruned_nolongLD_no_freq_discrep \
		--extract /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step1/comm_var_maf01.snplist \
		--covarFile ~/CD_behaviour/pheno/ibdbr_pheno.txt \
		--phenoFile ~/CD_behaviour/pheno/ibdbr_pheno.txt \
		--phenoCol ${pheno} \
		--covarColList PC1,PC2,PC3,PC4 \
		--bsize 1000 \
		--bt --lowmem --loocv \
		--out /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step1/${pheno}/${cov_type}

else
	regenie \
                --step 1 \
                --bed /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/genotyped_data/QC/release_interval_ibdbioresource_Aug24/interval_ibdbioresource_postqc_hg19_nomonom_RefAlt_flip2_pruned_nolongLD_no_freq_discrep \
		--extract /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step1/comm_var_maf01.snplist \
		--covarFile ~/CD_behaviour/pheno/ibdbr_pheno.txt \
                --phenoFile ~/CD_behaviour/pheno/ibdbr_pheno.txt \
                --phenoCol ${pheno} \
                --covarColList PC1,PC2,PC3,PC4,timeSinceDiag_5yrs \
                --catCovarList smoker_type,location,location1,year_diag_section \
                --bsize 1000 \
                --bt --lowmem --loocv \
                --out /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step1/${pheno}/${cov_type}

fi

