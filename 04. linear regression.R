
## A. Create cohort specific sample list from finalized phenotypic sheet for filtering samples in GWAS files; Tool: Command-line on linux terminal
study_list<-unique(clinical_cases_sporadic_pcs$study_recoded)
for(study in study_list)
{
  file=paste0(study, "_sample_list.txt")
  clinical_cases_sporadic_pcs %>% filter(study_recoded == study) %>% select("FID", "IID") %>%  write.table(file, row.names = F, quote = FALSE)
}


## B. Create GWAS data for age at onset GWAS from in-house PD GWAS with finalized for age at onset; Tool: Plink on linux terminal
### Repeat the process for each study cohort
plink --bfile study_gwasdata --keep study_sample_list.txt --keep-allele-order --make-bed --out study_ageatonset

## C. Convert plink files to vcf files; Tool: Plink on linux terminal
### Repeat the process for each study cohort
plink --bfile study_ageatonset --recode vcf --keep-allele-order --out study_ageatonset --threads 3 --memory 30000 

## D. Run Linear regression (Cohort specific); Tool: rvtest on Linux terminal 
### Repeat the process for each study cohort and for each subgroup (e.g.males as shown below in gender specific regression command)
rvtest --inVcf study_ageatonset.vcf --pheno study_pheno.ped --pheno-name y1 --covar study_covar.covar --covar-name sex,pc1,pc2,pc3,pc4,pc5 --out study --single wald,score
### Gender specific linear regression (Cohort specific); Tool: rvtest on Linux terminal 
rvtest --inVcf study_ageatonset.vcf --pheno study_pheno.ped --pheno-name y1 --covar study_covar.covar --covar-name pc1,pc2,pc3,pc4,pc5 --peopleIncludeFile study_males_sample_list.txt --out study_male --single wald,score

## E. Format the regression output (Cohort specific); Tool: Standard R package on Windows/Terminal  
### Repeat the process for each study cohort
fread("study.SingleScore.assoc") %>% mutate(SNP = paste(CHROM, POS, sep = ":")) %>% rename(N = N_INFORMATIVE) %>% fwrite("study_ageatonset_summstat.txt")

## F. Meta-analysis; Tool: Metal on Terminal
metal metal_script_ageatonset_consortium_date.txt

