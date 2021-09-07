# Merge Principal compnents in the phenotypic data file
# Folder with files with cohort names
all.files <- list.files(path = ".", pattern = ".txt")
l <- lapply(all.files, data.table::fread)
dt <- data.table::rbindlist(l)
dt <-
  dt %>% rename(
    FID = V1,
    IID = V2,
    PC1 = V3,
    PC2 = V4,
    PC3 = V5,
    PC4 = V6,
    PC5 = V7
  )
dt <-  dt %>% select(-IID)
clinical_cases_sporadic_pcs <-
  clinical_cases_sporadic %>% left_join(dt)
clinical_cases_sporadic_pcs %>% select(
  FID,
  IID,
  study_recoded,
  gender,
  age_study,
  age_diag,
  age_onset,
  pd_duration,
  PC1,
  PC2,
  PC3,
  PC4,
  PC5
) %>% data.table::fwrite("phenotype_data_ageatonset.txt")
