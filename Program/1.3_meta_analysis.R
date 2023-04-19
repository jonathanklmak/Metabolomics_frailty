#==============================================================================
# FILENAME: 1.3_meta_analysis.R
# PROJECT: 	Metabolomics_frailty
# PURPOSE:  To perform meta-analysis for observational results
# AUTHOR:   Jonathan Mak
# CREATED:	2023-02-02
# UPDATED: 	2023-02-21
# R VERSION: 4.1.3
#==============================================================================

### Required packages
library(dplyr)
library(metafor)

#=================== IMPORT PREVIOUSLY SAVED RESULTS ==========================

### Results from UK Biobank
ukb_lm_univariate_fi_age_sex_168metabolomics <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_age_sex_168metabolomics.txt", sep=" ")
ukb_lm_univariate_fi_age_sex_32biomarkers <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_age_sex_32biomarkers.txt", sep=" ")
ukb_lm_univariate_fi_allcov_168metabolomics <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics.txt", sep=" ")
ukb_lm_univariate_fi_allcov_32biomarkers <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers.txt", sep=" ")
ukb_coef_lasso_fi_249metabolomics <- read.delim("Output/Observational_results/ukb_coef_lasso_fi_249metabolomics.txt", sep=" ")
ukb_coef_lasso_fi_32biomarker <- read.delim("Output/Observational_results/ukb_coef_lasso_fi_32biomarker.txt", sep=" ")

### Results from TwinGene
twingene_gee_univariate_fi_age_sex_biomarkers <- read.delim("Output/Observational_results/twingene_gee_univariate_fi_age_sex_biomarkers.txt", sep=" ")
twingene_gee_univariate_fi_allcov_biomarkers <- read.delim("Output/Observational_results/twingene_gee_univariate_fi_allcov_biomarkers.txt", sep=" ")

### Results from Health 2000
results_health2000_met_fi <- read.csv("Output/Observational_results/output12_FI_ass_w_48biomarkers_01022023_n6073.csv", sep=";")

# Change variable names
results_health2000_met_fi <- results_health2000_met_fi %>% mutate(variable=T2000_Output1_01022023.w.6073.individuals..FIx100.lm.biomarker.age.sex.)
results_health2000_met_fi[results_health2000_met_fi$variable=="GP",]$variable <- "GlycA"
results_health2000_met_fi[results_health2000_met_fi$variable=="VLDL_D",]$variable <- "VLDL_size"
results_health2000_met_fi[results_health2000_met_fi$variable=="LDL_D",]$variable <- "LDL_size"
results_health2000_met_fi[results_health2000_met_fi$variable=="HDL_D",]$variable <- "HDL_size"
results_health2000_met_fi[results_health2000_met_fi$variable=="PC",]$variable <- "Phosphatidylc"
results_health2000_met_fi[results_health2000_met_fi$variable=="SM",]$variable <- "Sphingomyelins"
results_health2000_met_fi[results_health2000_met_fi$variable=="UnSat",]$variable <- "Unsaturation"
results_health2000_met_fi[results_health2000_met_fi$variable=="FAw6",]$variable <- "Omega_6"
results_health2000_met_fi[results_health2000_met_fi$variable=="Glc",]$variable <- "Glucose"
results_health2000_met_fi[results_health2000_met_fi$variable=="Cit",]$variable <- "Citrate"
results_health2000_met_fi[results_health2000_met_fi$variable=="Ace",]$variable <- "Acetate"
results_health2000_met_fi[results_health2000_met_fi$variable=="bOHBut",]$variable <- "bOHbutyrate"
results_health2000_met_fi[results_health2000_met_fi$variable=="Crea",]$variable <- "Creatinine"
results_health2000_met_fi[results_health2000_met_fi$variable=="Alb",]$variable <- "Albumin"
results_health2000_met_fi[results_health2000_met_fi$variable=="GGT_129",]$variable <- "GGT_serum"
results_health2000_met_fi[results_health2000_met_fi$variable=="KOL_114",]$variable <- "CHOL_serum"
results_health2000_met_fi[results_health2000_met_fi$variable=="KOL_LDL_L_116",]$variable <- "LDL_serum"
results_health2000_met_fi[results_health2000_met_fi$variable=="TRIGLY_124",]$variable <- "TRIG_serum"
results_health2000_met_fi[results_health2000_met_fi$variable=="LIPOB_122",]$variable <- "APOB_serum"
results_health2000_met_fi[results_health2000_met_fi$variable=="B_GHb_A1C",]$variable <- "HBA1C_RBC"
results_health2000_met_fi[results_health2000_met_fi$variable=="D_25_136_ST",]$variable <- "VITD_serum"

health2000_lm_univariate_fi_age_sex_biomarkers <- results_health2000_met_fi[,c(14,4:7)]
colnames(health2000_lm_univariate_fi_age_sex_biomarkers) <- c("variable","beta","lb","ub","p")
health2000_lm_univariate_fi_age_sex_biomarkers$se <- (health2000_lm_univariate_fi_age_sex_biomarkers$ub - health2000_lm_univariate_fi_age_sex_biomarkers$beta)/1.96
health2000_lm_univariate_fi_age_sex_biomarkers$n <- 6073
health2000_lm_univariate_fi_age_sex_biomarkers$study <- "Health 2000"
health2000_lm_univariate_fi_age_sex_biomarkers <- health2000_lm_univariate_fi_age_sex_biomarkers[,c("variable","beta","se","p","n","study")]

health2000_lm_univariate_fi_allcov_biomarkers <- results_health2000_met_fi[,c(14,10:13)]
colnames(health2000_lm_univariate_fi_allcov_biomarkers) <- c("variable","beta","lb","ub","p")
health2000_lm_univariate_fi_allcov_biomarkers$se <- (health2000_lm_univariate_fi_allcov_biomarkers$ub - health2000_lm_univariate_fi_allcov_biomarkers$beta)/1.96
health2000_lm_univariate_fi_allcov_biomarkers$n <- 6073
health2000_lm_univariate_fi_allcov_biomarkers$study <- "Health 2000"
health2000_lm_univariate_fi_allcov_biomarkers <- health2000_lm_univariate_fi_allcov_biomarkers[,c("variable","beta","se","p","n","study")]


#======================== PREPARATION OF DATASET =============================

### Select only the identified biomarkers for further analyses
selected_metabolomics <- merge(ukb_lm_univariate_fi_allcov_168metabolomics, # 1. significant in fully-adjusted models
                               ukb_coef_lasso_fi_249metabolomics, by="variable") %>% # 2. retained from LASSO
  filter(bon.p.fi.allcov & beta.LASSO.fi!=0) %>% 
  select(variable) %>% unlist() %>% as.character() # 41 NMR metabolomic biomarkers selected 

selected_clinical_biomarkers <- merge(ukb_lm_univariate_fi_allcov_32biomarkers, # 1. significant in fully-adjusted models
                                      ukb_coef_lasso_fi_32biomarker, by="variable") %>% # 2. retained from LASSO
  filter(bon.p.fi.allcov & beta.LASSO.fi!=0) %>% # Those missing in replication cohort are also selected
  select(variable) %>% unlist() %>% as.character() # 18 clinical biomarkers selected 

### Combine datasets
combined <- rbind(ukb_lm_univariate_fi_allcov_168metabolomics[ukb_lm_univariate_fi_allcov_168metabolomics$variable %in% selected_metabolomics,] %>%
                    rename(beta=beta.fi.allcov, se=se.fi.allcov, p=pvalue.fi.allcov) %>%
                    mutate(n=90573, study="UKB") %>%
                    select(variable, beta, se, p, n, study),
                  
                  ukb_lm_univariate_fi_allcov_32biomarkers[ukb_lm_univariate_fi_allcov_32biomarkers$variable %in% selected_clinical_biomarkers,] %>%
                    rename(beta=beta.fi.allcov, se=se.fi.allcov, p=pvalue.fi.allcov) %>%
                    mutate(n=67488, study="UKB") %>%
                    select(variable, beta, se, p, n, study),
                  
                  twingene_gee_univariate_fi_allcov_biomarkers[twingene_gee_univariate_fi_allcov_biomarkers$variable %in% c(selected_metabolomics,selected_clinical_biomarkers),] %>%
                    rename(beta=beta.allcov.twingene, se=se.allcov.twingene, p=pvalue.allcov.twingene) %>%
                    mutate(n=11025, study="TwinGene") %>%
                    select(variable, beta, se, p, n, study),
                  
                  health2000_lm_univariate_fi_allcov_biomarkers) %>% arrange(variable)

write.table(combined, "Output/Observational_results/combined_allcov_UKB_TwinGene_Health2000.txt", row.names = F)


#========================= PERFORM META-ANALYSES =============================

### Fully adjusted models for TwinGene and Health 2000
meta_TwinGene_Health2000_allcov <- data.frame(variable=c(selected_metabolomics,selected_clinical_biomarkers[c(2,3,8,11,17)]),
                                                 beta=NA, se=NA, lb=NA, ub=NA, p=NA, I2=NA)

for (i in 1:length(c(selected_metabolomics,selected_clinical_biomarkers[c(2,3,8,11,17)]))) {
  
  # Random-effects meta-analysis
  m <- rma(yi = beta,
           sei = se,
           data = combined %>% 
             filter(study!="UKB" & variable==c(selected_metabolomics,selected_clinical_biomarkers[c(2,3,8,11,17)])[i]),
           method = "DL")
  
  meta_TwinGene_Health2000_allcov[i,"beta"] <- m$beta
  meta_TwinGene_Health2000_allcov[i,"se"] <- m$se
  meta_TwinGene_Health2000_allcov[i,"lb"] <- m$ci.lb
  meta_TwinGene_Health2000_allcov[i,"ub"] <- m$ci.ub
  meta_TwinGene_Health2000_allcov[i,"p"] <- m$pval
  meta_TwinGene_Health2000_allcov[i,"I2"] <- m$I2
}
rm(list=c("i","m"))
write.table(meta_TwinGene_Health2000_allcov, "Output/Observational_results/meta_TwinGene_Health2000_allcov.txt")


### Fully adjusted models for UKB, TwinGene and Health 2000
meta_UKB_TwinGene_Health2000_allcov <- data.frame(variable=c(selected_metabolomics,selected_clinical_biomarkers[c(2,3,5,7,8,11,17,18)]),
                                              beta=NA, se=NA, lb=NA, ub=NA, p=NA, I2=NA)

for (i in 1:length(c(selected_metabolomics,selected_clinical_biomarkers[c(2,3,5,7,8,11,17,18)]))) {
  
  # Random-effects meta-analysis
  m <- rma(yi = beta,
           sei = se,
           data = combined %>% 
             filter(variable==c(selected_metabolomics,selected_clinical_biomarkers[c(2,3,5,7,8,11,17,18)])[i]),
           method = "DL")
  
  meta_UKB_TwinGene_Health2000_allcov[i,"beta"] <- m$beta
  meta_UKB_TwinGene_Health2000_allcov[i,"se"] <- m$se
  meta_UKB_TwinGene_Health2000_allcov[i,"lb"] <- m$ci.lb
  meta_UKB_TwinGene_Health2000_allcov[i,"ub"] <- m$ci.ub
  meta_UKB_TwinGene_Health2000_allcov[i,"p"] <- m$pval
  meta_UKB_TwinGene_Health2000_allcov[i,"I2"] <- m$I2
}
rm(list=c("i","m"))
write.table(meta_UKB_TwinGene_Health2000_allcov, "Output/Observational_results/meta_UKB_TwinGene_Health2000_allcov.txt")

# =============================== END OF FILE  ===============================