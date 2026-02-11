# genetic associations analyses where run for phenotypes on UK Biobank (UKB) whole exome sequencing data via the UKB Research Analysis Platform (RAP) hosted by DNAnexus
# the following prvoides the framework run for sex combined analyses of quantitative traits
#this script outlines analyses that were applied to the first set core phenotypes, that were also later applied to a wider range of quant traits
# sex specific analyses were additionally run using the additional regenie flags --sex-specific male and --sex-specific female 
# regenie mask files contained variants that passed QC parameters so only high quality variants were included in each burden test (see methods)

module load anaconda/3.2019-10
conda activate dx
############################################
######### run regenie step 1 UKB ###########
############################################
#swiss army knige 4.11.1 - contains regenie version: 3.1.1
run_regenie_step1="regenie --step 1\
 --lowmem --out core_pheno --bed ukb22418_c1_22_v2_merged\
 --phenoFile core_phenos_visit1_participant_processed_regenie.pheno --covarFile GWAS_covariates_REGENIE.txt  --keep KEEPFOR_White_Euro.txt\
 --extract GWAS_array_snps_qc_pass.snplist \
--phenoCol ivn_Body_mass_index_BMI --phenoCol ivn_Hip_circumference --phenoCol ivn_Standing_height --phenoCol ivn_Seated_height --phenoCol ivn_Sitting_height --phenoCol ivn_Body_fat_percentage --phenoCol ivn_Heel_bone_BMD --phenoCol ivn_Heel_BMD_Tscore --phenoCol ivn_diastolic_bp --phenoCol ivn_pulse_rate --phenoCol ivn_systolic_bp --phenoCol ivn_WHR --phenoCol ivn_WHRadjBMI --phenoCol ivn_mean_hand_grip_strength \
--covarCol age --covarCol age2 --covarCol sex --covarCol PC1 --covarCol PC2 --covarCol PC3 --covarCol PC4 --covarCol PC4 --covarCol PC5 --covarCol PC6 --covarCol PC7 --covarCol PC8 --covarCol PC9 --covarCol PC10 --covarCol genotyping_batch\
 --bsize 1000 --use-relative-path --gz --threads 16"

#need file paths here but not in the run_regenie_step1 command as these files will be downloaded to the cloud to run
dx run swiss-army-knife -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bed" \
-iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bim" \
-iin="${data_file_dir}/ukb22418_c1_22_v2_merged.fam"\
-iin="${data_file_dir}/core_phenos_visit1_participant_processed_regenie.pheno" \
-iin="${project}:/resources/GWAS_covariates_REGENIE.txt" \
-iin="${project}:/resources/KEEPFOR_White_Euro.txt" \
-iin="${data_file_dir}/genotype_array_QC/GWAS_array_snps_qc_pass.snplist" \
-icmd="${run_regenie_step1}" --tag="TNAP_step1" --instance-type "mem2_ssd1_v2_x16"\
--destination="${project}:/Data/Regenie_step1" --brief --yes

############################################
######### run regenie step 2 UKB ###########
############################################
#Regenie Step 2 - Association Analyses

#core phenotypes 
run_regenie_cmd="regenie --step 2 --out core_TNAP_additional_domain_assoc \
--bgen ukb23159_c1_b0_v1.bgen \
--sample ukb23159_c1_b0_v1.sample \
    --qt \
--ref-first \
 --phenoFile core_phenos_visit1_participant_processed_regenie.pheno --covarFile GWAS_covariates_REGENIE.txt  --keep KEEPFOR_White_Euro.txt\
 --covarCol age --covarCol age2 --covarCol sex --covarCol PC1 --covarCol PC2 --covarCol PC3 --covarCol PC4 --covarCol PC4 --covarCol PC5 --covarCol PC6 --covarCol PC7 --covarCol PC8 --covarCol PC9 --covarCol PC10 \
--pred core_pheno_pred.list --bsize 1000\
 --set-list TNAP_ALPL_set_file_additional_domains_incl_overlapping.txt \
--anno-file TNAP_ALPL_annotation_file_additional_domains_as_genes_incl_overlapping.txt \
--mask-def TNAP_ALPL_mask_definition_combined_incl_overlapping.txt \
--aaf-bins 0.0001,0.001,0.01,0.05,0.1,0.5,1 \
--write-mask-snplist \
--vc-MACthr 100 \
--vc-tests skato,acato-full \
--minMAC 1"

dx run swiss-army-knife -iin="${project}:/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release/ukb23159_c1_b0_v1.bgen" \
-iin="${project}:/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release/ukb23159_c1_b0_v1.bgen.bgi" \
-iin="${project}:/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release/ukb23159_c1_b0_v1.sample" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_pred.list" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_1.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_2.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_3.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_4.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_5.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_6.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_7.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_8.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_9.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_10.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_11.loco.gz" \
-iin="${project}:/Data/Regenie_step1/qc_core_pheno/core_pheno_12.loco.gz" \
-iin="/Data/core_phenos_visit1_participant_processed_regenie.pheno" \
-iin="/resources/GWAS_covariates_REGENIE.txt" \
-iin="/resources/KEEPFOR_White_Euro.txt" \
-iin="${project}:/TNAP_association/additional_domains/input/TNAP_ALPL_annotation_file_additional_domains_as_genes_incl_overlapping.txt" \
-iin="${project}:/TNAP_association/additional_domains/input/TNAP_ALPL_mask_definition_combined_incl_overlapping.txt" \
-iin="${project}:/TNAP_association/additional_domains/input/TNAP_ALPL_set_file_additional_domains_incl_overlapping.txt" \
-icmd="${run_regenie_cmd}" --tag="TNAP_Step2" --instance-type "mem1_ssd1_v2_x16"\
--destination="${project}:/TNAP_association/domains/results/" --brief --yes


