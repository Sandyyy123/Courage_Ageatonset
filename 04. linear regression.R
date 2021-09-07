

plink --bfile PLINK_R08_BARDIEN_CAUCASIAN_reference_updated --keep bardien_sample_list.txt --keep-allele-order --make-bed --out bardien_ageatonset

plink --bfile cohort_name_ageatonset --recode vcf --keep-allele-order --out cohort_name_ageatonset --threads 3 --memory 30000 


# Step : Linear regression (Cohort specific); Tool: rvtest on Linux terminal 
rvtest --inVcf cohort_name_ageatonset.vcf --pheno cohort_name_pheno.ped --pheno-name y1 --covar cohort_name_covar.covar --covar-name sex,pc1,pc2,pc3,pc4,pc5 --out cohort_name --single wald,score

# Step: Gender specific linear regression (Cohort specific); Tool: rvtest on Linux terminal 
rvtest --inVcf cohort_name_ageatonset.vcf --pheno cohort_name_pheno.ped --pheno-name y1 --covar cohort_name_covar.covar --covar-name pc1,pc2,pc3,pc4,pc5 --peopleIncludeFile cohort_name_males_sample_list.txt --out cohort_male --single wald,score

# Step: Formating the regression output (Cohort specific); Tool: Standard R package on Windows/Terminal  
fread("cohort_name.SingleScore.assoc") %>% mutate(SNP = paste(CHROM, POS, sep = ":")) %>% rename(N = N_INFORMATIVE) %>% fwrite("cohort_name_ageatonset_summstat.txt")

# Step: Meta-analysis; Tool: Metal on Terminal
metal metal_script_ageatonset_courage_250221.txt
