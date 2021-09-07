# PRS SNP extraction from case control analysis
## A. Load relevant libraries
library(tidyverse)
library(data.table)
library(readxl)

## Load top hits from case control analysis
df <- read_excel("courage+ipdgc_010721.xlsx")
mergedageatonset_bim <- fread("mergedageatonset_pooled.bim")




write.table(prs_bim, "snpforprs.bim", row.names = F, quote = F, sep = "\t")

patterns <- paste0("^", df$`CHR:BP`)
result <- filter(mergedageatonset_bim, grepl(paste(patterns, collapse="|"), mergedageatonset_bim$V2))
write.table(result$V2, "snplist.txt", col.names = F, row.names = F, quote = F)


filter(mergedageatonset_bim, grepl("^1:88362", mergedageatonset_bim$V2))

plink --bfile mergedageatonset_pooled --extract snplist.txt  --recode --out snpforprs
prs_ped <- fread("snpforprs.ped")
prs_map <- fread("snpforprs.map")
prs_ped <- prs_ped %>% unite("1:152192927",  V7:V8, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("1:155033317",  V9:V10, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("1:155135036",  V11:V12, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("1:156007988",  V13:V14, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("1:205663478",  V15:V16, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("1:205744546",  V17:V18, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("2:135537119",  V19:V20, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("2:169110394",  V21:V22, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("3:18199602",  V23:V24, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("3:28705690",  V25:V26, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("3:182760073",  V27:V28, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:893687",  V29:V30, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:941518",  V31:V32, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:951947",  V33:V34, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:1030374",  V35:V36, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:15723514",  V37:V38, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:15737348",  V39:V40, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:77139510",  V41:V42, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:90546486",  V43:V44, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:90666041",  V45:V46, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:90744993",  V47:V48, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:90757294",  V49:V50, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:90949825",  V51:V52, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:91211524",  V53:V54, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:94151811",  V55:V56, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:95194740",  V57:V58, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("4:103037478",  V59:V60, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("5:60137959",  V61:V62, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("5:60247815",  V63:V64, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("5:60330171",  V65:V66, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("5:60395289",  V67:V68, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("6:32213052",  V69:V70, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("6:32582650",  V71:V72, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("6:133118216",  V73:V74, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("7:23391509",  V75:V76, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("8:16701281",  V77:V78, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("9:17579763",  V79:V80, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("11:14913645",  V81:V82, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("11:133774863",  V83:V84, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:40120100",  V85:V86, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:40388109",  V87:V88, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:40620808",  V89:V90, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:40734202",  V91:V92, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:40773684",  V93:V94, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:40787742",  V95:V96, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:41303140",  V97:V98, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:41506535",  V99:V100, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:41589415",  V101:V102, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:123093541",  V103:V104, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("12:123326598",  V105:V106, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("15:61993702",  V107:V108, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("16:20242614",  V109:V110, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("16:30666367",  V111:V112, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("16:30923602",  V113:V114, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("17:43573419",  V115:V116, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("17:43832337",  V117:V118, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("17:44121579",  V119:V120, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("17:44222460",  V121:V122, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("17:44293020",  V123:V124, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("17:44801340",  V125:V126, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("18:40672964",  V127:V128, sep = "",  remove = T)
prs_ped <- prs_ped %>% unite("22:41755105",  V129:V130, sep = "",  remove = T)
write.csv(prs_ped, "prs_input_aao.csv")

prs_ped <- prs_ped %>% rename(FID = V1)
prs_input_cov <- left_join(clinical_cases_sporadic_pcs, prs_ped)
prs_input_cov <- prs_input_cov %>% rename(SEX = gender) %>% select(FID, IID, SEX, PC1:PC5) %>% write.table("prs_cov_aao.txt", row.names = F, quote = F, sep = "\t")
prs_input_cov <- left_join(clinical_cases_sporadic_pcs, prs_ped)

prs_input_cov <-  prs_input_cov %>% rename(PHENOTYPE = age_onset) %>% select(FID, IID, PHENOTYPE) %>% write.table("prs_pheno_aao.txt", row.names = F, quote = F, sep = "\t")
df <- read_excel("courage+ipdgc_010721.xlsx")

df <- df %>% rename(snp = MarkerName, A1 = Allele1, A2= Allele2, beta = Effect, pvalue = 'P-value of Courage + IPDGC meta-analysis') %>% separate('CHR:BP', into = c("chr", "pos"), sep = ":")
df <- df %>% select(snp, chr, pos, A1, A2, beta, pvalue)
write.table(df, "prs_base_ipdgc_courage.assoc", row.names = F, quote = F, sep = "\t")
prs_ped <- fread("snpforprs.ped")
prs_ped %>% 
  prs_bim <- fread("snpforprs.bim")
prs_bim$V2 <- str_sub(prs_bim$V2, end=-5)
write.table(prs_bim, "snpforprs.bim", row.names = F, quote = F, sep = "\t")



# extract SNPdata
plink --bfile mergedageatonset_pooled --extract snplist.txt --make-bed --out snpforprs
#update rsid
plink --bfile snpforprs --update-name updatersid.txt --make-bed --out snpforprs



setwd("D:/courage_data/data/metanalysis_05_05_2021")
Rscript PRSice.R --dir . \
--prsice ./PRSice_linux \
--base prs_base_ipdgc_courage.assoc \
--ignore-fid \
--snp snp --chr chr --bp bp --A1 A1 --A2 A2 --stat beta --pvalue pvalue \
--target snpforprs \
--pheno prs_pheno_aao.txt \
--pheno-col PHENOTYPE \
--cov prs_cov_aao.txt \
--cov-col SEX,PC1,PC2,PC3,PC4,PC5 \
--score sum \
--stat beta \
--beta \
--binary-target F \
--perm 10000 \
--cov-factor SEX \


prsice_data <- fread("PRSice.summary")
