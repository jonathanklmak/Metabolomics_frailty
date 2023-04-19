#==============================================================================
# FILENAME: 2.2_MR_frailty_sensitivity.R
# PROJECT: 	Metabolomics_frailty
# PURPOSE:  To perform MR sensitivity analyses for metabolites and frailty
# AUTHOR:   Jonathan Mak
# CREATED:	2022-09-14
# UPDATED: 	2023-02-21
# R VERSION: 4.1.3
#==============================================================================

### Required packages
library(dplyr)
library(ieugwasr)
library(ggforestplot) # For forest plots of metabolites
library(ggplot2); library(patchwork); library(ggforce) # Plotting graphs
library(TwoSampleMR) # To perform MR analysis
library(MRPRESSO) # MR PRESSO

### Load previously saved MR results
load("Data/R_data/Metabolites_frailty_MR.Rdata")


#============== SENSITIVITY ANALYSIS: REMOVE PLEIOTROPIC SNPs ================

metabolomics_iv_nonoverlap <- metabolomics_iv[!(duplicated(metabolomics_iv$SNP) | duplicated(metabolomics_iv$SNP, fromLast = TRUE)), ]

fi_nometabolomics_gwas <- read.delim("Data/GWAS_summary_statistics_data/Frailty/fi_nometabolomics.fastGWA", header=T)
fp_nometabolomics_gwas <- read.delim("Data/GWAS_summary_statistics_data/Frailty/fp_nometabolomics.fastGWA", header=T)
mr_outcome_nometabolomics_nonoverlap <- format_data(rbind(fi_nometabolomics_gwas %>% mutate(outcome="FI"),
                                                          fp_nometabolomics_gwas %>% mutate(outcome="FP")),
                                                    type = "outcome",
                                                    snps = metabolomics_iv_nonoverlap$SNP,
                                                    phenotype_col = "outcome",
                                                    snp_col = "SNP",
                                                    beta_col = "BETA",
                                                    se_col = "SE",
                                                    eaf_col = "AF1",
                                                    effect_allele_col = "A1",
                                                    other_allele_col = "A2",
                                                    pval_col = "P",
                                                    pos_col = "POS",
                                                    samplesize_col = "N")

mr_data_metabolomics_nonoverlap <- harmonise_data(metabolomics_iv[!(duplicated(metabolomics_iv$SNP) | duplicated(metabolomics_iv$SNP, fromLast = TRUE)), ], mr_outcome_nometabolomics_nonoverlap, action = 2)
mr_metabolomics_results_nonoverlap <- mr(mr_data_metabolomics_nonoverlap, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
mr_metabolomics_het_nonoverlap <- mr_heterogeneity(mr_data_metabolomics_nonoverlap)
mr_metabolomics_plt_nonoverlap <- mr_pleiotropy_test(mr_data_metabolomics_nonoverlap)

rm(list=c("fi_nometabolomics_gwas","fp_nometabolomics_gwas"))


#====== SENSITIVITY ANALYSIS: REMOVE ITEMS FROM EACH CATEGORY OF FI ===========

### 11 modified FIs in UKB, removing items from each category
mod_fi_list <- c("nosensory", # 1. Sensory (3 items): Glaucoma, Cataracts, Hearing difficulty
                 "nocranial", # 2. Cranial (2 items): Migraine, Dental problems
                 "nomental",  # 3. Mental wellbeing (8 items): Self-rated health, Fatigue, Sleep, Depressed feelings, Self-described nervous personality, Severe anxiety/ panic attacks, Common to feel loneliness, Sense of misery
                 "noinfirmity", # 4. Infirmity (3 items): Infirmity, Falls in last year, Fractures/broken bones in last five years
                 "nocardio", # 5. Cardiometabolic (8 items): diabetes, myocardial infarction, angina, stroke, high blood pressure, hypothyroidism, deep-vein thrombosis, high cholesterol
                 "norespiratory", # 6. Respiratory (4 items): Breathing: wheeze in last year, Pneumonia, Chronic bronchitis/emphysema, Asthma
                 "nomusculo", # 7. Musculoskeletal (4 items): Rheumatoid arthritis, Osteoarthritis, Gout, Osteoporosis
                 "noimmune", # 8. Immunological (2 items): Hayfever, Psoriasis 
                 "nocancer", # 9. Cancer (2 items): Any cancer diagnosis, Multiple cancers diagnosed (number reported)
                 "nopain", # 10. Pain (9 items): Chest pain, Head and/or neck pain, Back pain, Stomach/abdominal pain, Hip pain, Knee pain, Whole-body pain, Facial pain, Sciatica
                 "nogastro" # 11. Gastrointestinal (4 items): Gastric reflux, Hiatus hernia, Gall stones, Diverticulitis
)


### Repeat MR analysis for each modified FI
for (i in 1:length(mod_fi_list)) {
  
  # Import GWAS summary statistics for the modified FI
  assign(paste0("fi_",mod_fi_list[i],"_gwas"), read.delim(paste0("Data/GWAS_summary_statistics_data/Frailty/fi_",mod_fi_list[i],".fastGWA"), header=T))
  assign(paste0("fi_",mod_fi_list[i],"_nometabolomics_gwas"), read.delim(paste0("Data/GWAS_summary_statistics_data/Frailty/fi_",mod_fi_list[i],"_no_metabolomics.fastGWA"), header=T))
  assign(paste0("fi_",mod_fi_list[i],"_random50_gwas"), read.delim(paste0("Data/GWAS_summary_statistics_data/Frailty/fi_",mod_fi_list[i],"_random50.fastGWA"), header=T))
  
  # Format outcome data for MR sensitivity analysis
  mr_metabolomics_sensitivity <- format_data(get(paste0("fi_",mod_fi_list[i],"_nometabolomics_gwas")) %>% mutate(outcome=paste0("FI, ",mod_fi_list[i])),
                                             type = "outcome",
                                             snps = metabolomics_iv[metabolomics_iv$exposure %in% c("Creatinine","GlycA","HDL_size","IDL_CE","IDL_FC","IDL_P","M_LDL_CE","M_LDL_P","MUFA","Omega_6","S_LDL_CE","S_LDL_L","S_LDL_P","S_LDL_PL","XS_VLDL_CE"),]$SNP,
                                             phenotype_col = "outcome",
                                             snp_col = "SNP",
                                             beta_col = "BETA",
                                             se_col = "SE",
                                             eaf_col = "AF1",
                                             effect_allele_col = "A1",
                                             other_allele_col = "A2",
                                             pval_col = "P",
                                             pos_col = "POS",
                                             samplesize_col = "N")
  mr_clinical_biomarkers_sensitivity <-  format_data(get(paste0("fi_",mod_fi_list[i],"_gwas")) %>% mutate(outcome=paste0("FI, ",mod_fi_list[i])),
                                                     type = "outcome",
                                                     snps = clinical_biomarkers_iv[clinical_biomarkers_iv$exposure %in% c("Total cholesterol","LDL-C","Triglycerides"),]$SNP,
                                                     phenotype_col = "outcome",
                                                     snp_col = "SNP",
                                                     beta_col = "BETA",
                                                     se_col = "SE",
                                                     eaf_col = "AF1",
                                                     effect_allele_col = "A1",
                                                     other_allele_col = "A2",
                                                     pval_col = "P",
                                                     pos_col = "POS",
                                                     samplesize_col = "N")
  mr_clinical_biomarkers_sensitivity_2 <-  format_data(get(paste0("fi_",mod_fi_list[i],"_random50_gwas")) %>% mutate(outcome=paste0("FI, ",mod_fi_list[i])),
                                                       type = "outcome",
                                                       snps = clinical_biomarkers_iv[clinical_biomarkers_iv$exposure %in% c("ApoB"),]$SNP,
                                                       phenotype_col = "outcome",
                                                       snp_col = "SNP",
                                                       beta_col = "BETA",
                                                       se_col = "SE",
                                                       eaf_col = "AF1",
                                                       effect_allele_col = "A1",
                                                       other_allele_col = "A2",
                                                       pval_col = "P",
                                                       pos_col = "POS",
                                                       samplesize_col = "N")
  
  ### Data Harmonization
  assign(paste0("mr_data_metabolomics_",mod_fi_list[i],"_sensitivity"),
         harmonise_data(metabolomics_iv[metabolomics_iv$exposure %in% c("Creatinine","GlycA","HDL_size","IDL_CE","IDL_FC","IDL_P","M_LDL_CE","M_LDL_P","MUFA","Omega_6","S_LDL_CE","S_LDL_L","S_LDL_P","S_LDL_PL","XS_VLDL_CE"),], mr_metabolomics_sensitivity, action = 2))
  
  assign(paste0("mr_data_clinical_biomarkers_",mod_fi_list[i],"_sensitivity"),
         rbind(harmonise_data(clinical_biomarkers_iv[clinical_biomarkers_iv$exposure %in% c("Total cholesterol","LDL-C","Triglycerides"),], mr_clinical_biomarkers_sensitivity, action = 2),
               harmonise_data(clinical_biomarkers_iv[clinical_biomarkers_iv$exposure %in% c("ApoB"),], mr_clinical_biomarkers_sensitivity_2, action = 2)))
  
  ### MR analysis
  assign(paste0("mr_results_metabolomics_",mod_fi_list[i],"_sensitivity"),
         mr(get(paste0("mr_data_metabolomics_",mod_fi_list[i],"_sensitivity")), method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode")))
  assign(paste0("mr_results_metabolomics_",mod_fi_list[i],"_sensitivity_het"),
         mr_heterogeneity(get(paste0("mr_data_metabolomics_",mod_fi_list[i],"_sensitivity"))))
  assign(paste0("mr_results_metabolomics_",mod_fi_list[i],"_sensitivity_plt"),
         mr_pleiotropy_test(get(paste0("mr_data_metabolomics_",mod_fi_list[i],"_sensitivity"))))
  
  assign(paste0("mr_results_clinical_biomarkers_",mod_fi_list[i],"_sensitivity"),
         mr(get(paste0("mr_data_clinical_biomarkers_",mod_fi_list[i],"_sensitivity")), method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode")))
  assign(paste0("mr_results_clinical_biomarkers_",mod_fi_list[i],"_sensitivity_het"),
         mr_heterogeneity(get(paste0("mr_data_clinical_biomarkers_",mod_fi_list[i],"_sensitivity"))))
  assign(paste0("mr_results_clinical_biomarkers_",mod_fi_list[i],"_sensitivity_plt"),
         mr_pleiotropy_test(get(paste0("mr_data_clinical_biomarkers_",mod_fi_list[i],"_sensitivity"))))
  
  ### Remove unnecessary objects
  rm(list=c(paste0("fi_",mod_fi_list[i],"_gwas"),paste0("fi_",mod_fi_list[i],"_nometabolomics_gwas"),paste0("fi_",mod_fi_list[i],"_random50_gwas"),
            "mr_metabolomics_sensitivity","mr_clinical_biomarkers_sensitivity","mr_clinical_biomarkers_sensitivity_2"))
}


### Summary of results

mr_summary_mod_fi_sensitivity <- c()
for (i in 1:length(mod_fi_list)) {
  for (j in c("Creatinine","GlycA","HDL_size","IDL_CE","IDL_FC","IDL_P","M_LDL_CE","M_LDL_P","MUFA","Omega_6","S_LDL_CE","S_LDL_L","S_LDL_P","S_LDL_PL","XS_VLDL_CE")) {
    mr_summary_mod_fi_sensitivity <-
      rbind(
        mr_summary_mod_fi_sensitivity,
        # Results for IVW, MR Egger, Weighted median, Weighted mode
        data.frame(get(paste0("mr_results_metabolomics_",mod_fi_list[i],"_sensitivity")) %>% filter(exposure==j) %>% select(outcome,exposure,method,nsnp,b,se,pval),
                   p_het=rbind(get(paste0("mr_results_metabolomics_",mod_fi_list[i],"_sensitivity_het")) %>% filter(exposure==j & method=="Inverse variance weighted") %>% select(Q_pval),get(paste0("mr_results_metabolomics_",mod_fi_list[i],"_sensitivity_het")) %>% filter(exposure==j & method=="MR Egger") %>% select(Q_pval),NA,NA),
                   p_plt=rbind(NA,get(paste0("mr_results_metabolomics_",mod_fi_list[i],"_sensitivity_plt")) %>% filter(exposure==j) %>% select(pval),NA,NA))
      )
  }
    for (j in c("Total cholesterol","LDL-C","Triglycerides","ApoB")) {
      mr_summary_mod_fi_sensitivity <-
        rbind(
          mr_summary_mod_fi_sensitivity,
          # Results for IVW, MR Egger, Weighted median, Weighted mode
          data.frame(get(paste0("mr_results_clinical_biomarkers_",mod_fi_list[i],"_sensitivity")) %>% filter(exposure==j) %>% select(outcome,exposure,method,nsnp,b,se,pval),
                     p_het=rbind(get(paste0("mr_results_clinical_biomarkers_",mod_fi_list[i],"_sensitivity_het")) %>% filter(exposure==j & method=="Inverse variance weighted") %>% select(Q_pval),get(paste0("mr_results_clinical_biomarkers_",mod_fi_list[i],"_sensitivity_het")) %>% filter(exposure==j & method=="MR Egger") %>% select(Q_pval),NA,NA),
                     p_plt=rbind(NA,get(paste0("mr_results_clinical_biomarkers_",mod_fi_list[i],"_sensitivity_plt")) %>% filter(exposure==j) %>% select(pval),NA,NA))
        )
    }
  }

clipr::write_clip(mr_summary_mod_fi_sensitivity) # Copy results
write.table(mr_summary_mod_fi_sensitivity, "Output/MR_results/mr_summary_mod_fi_sensitivity.txt", row.names = F)


#==================== SAVE DATA FOR FURTHER ANALYSIS =========================

save.image("Data/R_data/Sensitivity_frailty_MR.Rdata")

# =============================== END OF FILE  ===============================