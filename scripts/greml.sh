datadir="/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour"

outdir="/nfs/users/nfs_q/qz2/CD_behaviour"

pheno_idx=${LSB_JOBINDEX}
cov_type=$1

if [[ $1 == "pc_only" ]]
then
	gcta64 --reml --grm ${datadir}/GRM/ibdbr_eur_unrel --pheno ${outdir}/pheno/greml_pheno.txt --mpheno ${pheno_idx} --qcovar ${outdir}/pheno/greml_pcs.txt --out ${outdir}/GREML/CD_${cov_type}${pheno_idx} --prevalence 0.5
else
	gcta64 --reml --grm ${datadir}/GRM/ibdbr_eur_unrel --pheno ${outdir}/pheno/greml_pheno.txt --mpheno ${pheno_idx} --qcovar ${outdir}/pheno/greml_qcov.txt --covar ${outdir}/pheno/greml_ccov.txt --out ${outdir}/GREML/CD_${cov_type}${pheno_idx} --prevalence 0.5
fi

