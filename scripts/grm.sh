#!/bin/bash

gcta64 --bfile /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/GRM/data/ibdbr_eur --make-grm --out /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/GRM/ibdbr_eur

gcta64 --grm /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/GRM/ibdbr_eur --grm-cutoff 0.05 --make-grm --out /lustre/scratch125/humgen/projects_v2/ibdgwas_bioresource/qz2/CD_behaviour/GRM/ibdbr_eur_unrel

