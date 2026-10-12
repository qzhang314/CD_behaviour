#!/bin/bash

# analysis pipeline for CD behaviour project (B2/B3 v.s. B1)

wkdir="/nfs/users/nfs_q/qz2/CD_behaviour"

cd ${wkdir}

mkdir figures

mkdir pheno

mkdir result

# 1 phenotype definition in IBDBR
Rscript scripts/ibdbr_pheno.R

# 2 Heritability analysis (i.e. based on high quality common variants MAF > 1%) to determine the best control definition

## 2.1 build GRM 
### 2.1.1 extract common variants

mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/GRM
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/GRM/data
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log

module load HGI/softpack/groups/trynka/popgen/1.0

bsub -J GRM_data -q normal -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/GRM_data.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/GRM_data.err -W 10:00 -M 100000 -R'select[mem>100000] rusage[mem=100000]' plink --vcf /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/genotyped_data/imputation_ready/release_interval_ibdbioresource_Aug24/interval_ibdbioresource_eur_hg19.vcf.gz --double-id --keep /nfs/users/nfs_q/qz2/CD_behaviour/pheno/ibdbr_pheno.txt --maf 0.01 --make-bed --out /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/GRM/data/ibdbr_eur

### 2.1.2 calculate GRM & get genetically unrelated samples (genetic relatedness < 0.05)

module load HGI/softpack/groups/team152/genetics_tools/1

bsub -J GRM_calc -w "done(GRM_data)" -q normal -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/GRM_calc.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/GRM_calc.err -W 10:00 -M 100000 -R'select[mem>100000] rusage[mem=100000]' bash scripts/grm.sh

## 2.2 estimate heritability
mkdir GREML

awk '{print $1,$2,$22}' pheno/ibdbr_pheno.txt > pheno/greml_pheno.txt
awk '{print $1,$2,$6,$7,$8,$9}' pheno/ibdbr_pheno.txt > pheno/greml_pcs.txt
awk '{print $1,$2,$6,$7,$8,$9,$23}' pheno/ibdbr_pheno.txt > pheno/greml_qcov.txt
awk '{print $1,$2,$4,$5,$12,$25}' pheno/ibdbr_pheno.txt > pheno/greml_ccov.txt

### 2.2.1 with 4 PCs as covarites only
mkdir GREML/log

cov_type="pc_only"
bsub -J GREML_${cov_type} -q normal -o GREML/log/GREML_${cov_type}.log -e GREML/log/GREML_${cov_type}.err -W 10:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' bash scripts/greml.sh ${cov_type}

### 2.2.2 with 4 PCs and all other covariates
cov_type="full_cov"

bsub -J GREML_${cov_type} -q normal -o GREML/log/GREML_${cov_type}.log -e GREML/log/GREML_${cov_type}.err -W 10:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' bash scripts/greml.sh ${cov_type}


# based on the above result, I choose to not restricting the following-up time
# baseline model fitting based on the phenotype data
mkdir result

#phenotypic model development
Rscript scripts/model_fit.R


#effect of smoking dose and time of quitting smoking on behaviour
Rscript scripts/smoke_effect.R

#effect of delayed diagnosis on behaviour
Rscript scripts/delayed_diag_effect.R

#effect of drug prescription time on behaviour
Rscript scripts/drug_prescription_time_effect.R

#effect of era on behaviour after correcting for the drug use
Rscript scripts/era_effect_drug_correction.R

#effect of older age at diangnosis on behaviour
Rscript scripts/old_age_at_dianogsis_effect.R

#make a summary figure for all the dataset (Figure 1 in the manuscript)
Rscript scripts/pheno_smy_plot.R

#make a summary figure for all the effects (Figure 2 in the manuscript)
Rscript scripts/pheno_effect_smy_plot.R

#summary table of samples
Rscript scripts/sample_summary.R



# 3 perform GWAS by using regenie
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie

## 3.1 step 1
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step1

### get variant list with MAF > 1% in CD only data (if a variant with MAF = 0 in the selected samples, regenie will report error)

plink2 --bfile /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/genotyped_data/QC/release_interval_ibdbioresource_Aug24/interval_ibdbioresource_postqc_hg19_nomonom_RefAlt_flip2_pruned_nolongLD_no_freq_discrep --keep ~/CD_behaviour/pheno/ibdbr_pheno.txt --maf 0.01 --write-snplist --out /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step1/comm_var_maf01

### 3.1.1 with 4 PCs
bsub -J regenie_step1_pc_only -q normal -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/regenie_step1_pc_only.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/regenie_step1_pc_only.err -W 10:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' bash scripts/regenie_step1.sh "pc_only"


### 3.1.2 with 4 PCs and other covariates
bsub -J regenie_step1_full_cov -q normal -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/regenie_step1_full_cov.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/regenie_step1_full_cov.err -W 10:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' bash scripts/regenie_step1.sh "full_cov"

## 3.2 step 2, high quality common variants only (MAF > 1%, info > 0.3)
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step2
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step2/log

### 3.2.1 with 4 PCs

cov_type="pc_only"
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step2/${cov_type}

bsub -J regenie_step2_${cov_type}[1-23] -w "done(regenie_step1_pc_only)" -q long -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/regenie_step2_${cov_type}_%I.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/regenie_step2_${cov_type}_%I.err -W 40:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' bash scripts/regenie_step2.sh ${cov_type}


### 3.2.2 with 4 PCs and other covariates

cov_type="full_cov"
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/step2/${cov_type}

bsub -J regenie_step2_${cov_type}[1-23] -w "done(regenie_step1_full_cov)" -q long -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/regenie_step2_${cov_type}_%I.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/regenie_step2_${cov_type}_%I.err -W 40:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' bash scripts/regenie_step2.sh ${cov_type}


### liftover to b38, remove bad quality variant
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/QC

cov_types=("pc_only" "full_cov")

for cov_type in ${cov_types[@]}
do
	mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/QC/${cov_type}
	bsub -J smy_fmt_QC_${cov_type}[1-23] -q normal -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/smy_fmt_QC_${cov_type}_%I.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/smy_fmt_QC_${cov_type}_%I.err -W 1:00 -M 8000 -R'select[mem>8000] rusage[mem=8000]' Rscript scripts/format_ibdbioresource_interval_files_generated_with_regenie_b37.R ${cov_type}
done


### merge summary statistics
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics

head -n 1 /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/QC/pc_only/chr1_bev.regenie > /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_pc_only.regenie

tail -n +2 -q /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/QC/pc_only/chr{1..23}_bev.regenie >> /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_pc_only.regenie

head -n 1 /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/QC/full_cov/chr1_bev.regenie > /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_full_cov.regenie

tail -n +2 -q /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/regenie/QC/full_cov/chr{1..23}_bev.regenie >> /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_full_cov.regenie


### keep results of common variants only (MAF in total > 1%)
awk 'NR == 1 { print; next } {if($6 > 0.01 && $6 < 0.99) print $0}' /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_pc_only.regenie > /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_pc_only_maf01.regenie

awk 'NR == 1 { print; next } {if($6 > 0.01 && $6 < 0.99) print $0}' /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_full_cov.regenie > /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_full_cov_maf01.regenie


### manhattan plot and Q-Q plot

cov_types=("pc_only" "full_cov")

for cov_type in ${cov_types[@]}
do
        bsub -J gwas_smy_manh_plot_${cov_type} -q normal -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/gwas_smy_manh_plot_${cov_type}.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/gwas_smy_manh_plot_${cov_type}.err -W 1:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' Rscript scripts/manh_plot.R ${cov_type}
	bsub -J gwas_smy_qq_plot_${cov_type} -q normal -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/gwas_smy_qq_plot_${cov_type}.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/log/gwas_smy_qq_plot_${cov_type}.err -W 1:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' Rscript scripts/qq_plot.R ${cov_type}
done


# 4 HLA associations

Rscript scripts/ibdbr_hla.R


# 5 calculate CD behaviour PRS, using PC_only summary statistics

# 5.1 format summary statistics

awk '{if(NR!=1) print $19,$5,$4,$14,$15}' /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_pc_only_maf01.regenie|grep -v NA > /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_pc_only_maf01_PRScs.regenie
sed -i '1s/^/SNP A1 A2 BETA SE\n/' /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_pc_only_maf01_PRScs.regenie

## since PRS-CS use a common sample size for all SNPs, here I choose to use a median sample size across all SNPs, as Schultz et al. https://www.sciencedirect.com/science/article/pii/S2666247722000070

mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/PRScs
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/PRScs/log

module load HGI/softpack/users/cc53/base/7

# 5.2 re-estimate per SNP weight based on PRS-cs
bsub -J PRScs[1-22] -q normal -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/PRScs/log/PRScs_%I.log \
-e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/PRScs/log/PRScs_%I.err -W 10:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' \
python ~/software/program/PRScs/PRScs.py \
--ref_dir="/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/PRS/ref/ldblk_1kg_eur" \
--bim_prefix="/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/PRS/PRScs/bim/chr\${LSB_JOBINDEX}" \
--sst_file="/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/gwas_summary_statistics/CD_bev_pc_only_maf01_PRScs.regenie" --n_gwas=13952 \
--out_dir="/lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/PRScs/chr\${LSB_JOBINDEX}" --a=1 --b=0.5 --n_iter=1000 --n_burnin=500 \
--thin=5 --chrom=\${LSB_JOBINDEX} --seed=314


# 5.3 add ID to the summary statistics

mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/PRScs/weights

mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/PRScs/weights/log

bsub -J add_ID[1-22] -q normal -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/PRScs/weights/log/CD_bev_%I.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/PRScs/weights/log/CD_bev_%I.err -W 10:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' Rscript scripts/PRS_add_ID.R


# 6 UKIBDGC phenotype definition (excluding overlapped patients between IBDBR and UKIBDGC)

Rscript scripts/ukibdgc_pheno.R

Rscript scripts/ukibdgc_cov.R

# 7 calculate CD behaviour PRS on UKIBDGC

module load HGI/softpack/users/bh18/bcftools/1.0

mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/score
mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/score/log

cohorts=("all_hce" "gwas1" "gwas2")
for cohort in ${cohorts[@]}
do
	mkdir /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/score/${cohort}      
        bsub -J ${cohort}_PRS[1-22] -q normal -o /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/score/log/PRS_${cohort}_%I.log -e /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/score/log/PRS_${cohort}_%I.err -W 10:00 -M 20000 -R'select[mem>20000] rusage[mem=20000]' plink2 --bgen /lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/imputed_data/${cohort}/eur/chr\${LSB_JOBINDEX}.dose.hweremoved.bgen ref-first --sample /lustre/scratch125/humgen/projects_v2/ibdgwas_ukibdgc/imputed_data/${cohort}/eur/chr\${LSB_JOBINDEX}.dose.hweremoved.bgen.sample --score /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/PRScs/weights/chr\${LSB_JOBINDEX}.coef 8 4 6 header list-variants cols=+scoresums --out /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/score/${cohort}/chr\${LSB_JOBINDEX}
done

Rscript scripts/merge_PRS.R

# check the effect of PRS on disease progression
Rscript scripts/ukibdgc_prs_assoc.R


# 8 model validation in UKIBDGC
Rscript scripts/ukibdgc_validate.R


# show the effect in ukibdgc (Figure 3 in the manuscript)
Rscript scripts/ukibdgc_all_factor_effect.R


# 7 GWAS on UKIBDGC

## 7.1 step 1

#MAF > 1%, missingness < 10%, HWE P > 1e−15 and linkage disequilibrium (LD) pruning (1,000 variant windows, 100 variant sliding windows and r2 <0.9). Exclude from step 1 any SNPs in high LD region.

### 6.1.1 with 4 PCs

### 6.1.2 with 4 PCs and covariates

## 6.2 step 2, high quality common variants only (MAF > 1%, info > 0.3)

### 6.2.1 with 4 PCs

### 6.2.2 with 4 PCs and covariates

# 7 GWAS meta-analysis with UKIBDGC

# 8 prediction model

## 8.1 build phenotype-based prediction model based on IBDBR samples

## 8.2 validate model in UKIBDGC

## 8.3 combine phenotype-based prediction model with PRS in UKIBDGC (10-cross validation) 

