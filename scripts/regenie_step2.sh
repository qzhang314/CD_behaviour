pheno=$1
cov_type=$2
chrnum=${LSB_JOBINDEX}

if [ ${chrnum} == 23 ];then
        chrnum="X"
fi


if [[ ${cov_type} == "pc_only" ]]
then
	regenie \
		--step 2 \
		--bgen /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/imputed_data/release_interval_ibdbioresource_Aug24/combined_hrc_uk10k/${chrnum}_combined_hrc_uk10k_eur_with_consensus_ids.bgen \
		--sample /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/imputed_data/release_interval_ibdbioresource_Aug24/combined_hrc_uk10k/${chrnum}_combined_hrc_uk10k_eur_with_consensus_ids.sample \
		--covarFile ~/CD_behaviour/pheno/ibdbr_pheno.txt \
                --phenoFile ~/CD_behaviour/pheno/ibdbr_pheno.txt \
                --phenoCol ${pheno} \
                --covarColList PC1,PC2,PC3,PC4 \
		--bsize 1000 \
		--bt --lowmem --af-cc \
		--firth --approx --pThresh 0.01 \
		--minINFO 0.6 --minMAC 50 \
		--pred /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step1/${pheno}/pc_only_pred.list \
		--out /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step2/pc_only/${pheno}/chr${chrnum}

else
	regenie \
                --step 2 \
                --bgen /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/imputed_data/release_interval_ibdbioresource_Aug24/combined_hrc_uk10k/${chrnum}_combined_hrc_uk10k_eur_with_consensus_ids.bgen \
		--sample /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/imputed_data/release_interval_ibdbioresource_Aug24/combined_hrc_uk10k/${chrnum}_combined_hrc_uk10k_eur_with_consensus_ids.sample \
		--covarFile ~/CD_behaviour/pheno/ibdbr_pheno.txt \
                --phenoFile ~/CD_behaviour/pheno/ibdbr_pheno.txt \
                --phenoCol ${pheno} \
                --covarColList PC1,PC2,PC3,PC4,timeSinceDiag_5yrs \
                --catCovarList smoker_type,location,location1,year_diag_section \
                --bsize 1000 \
                --bt --lowmem --af-cc \
                --firth --approx --pThresh 0.01 \
                --minINFO 0.6 --minMAC 50 \
                --pred /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step1/${pheno}/full_cov_pred.list \
                --out /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step2/full_cov/${pheno}/chr${chrnum}

fi

