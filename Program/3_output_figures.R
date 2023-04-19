#==============================================================================
# FILENAME: 3_output_figures.R
# PROJECT: 	Metabolomics_frailty
# PURPOSE:  To summarize and produce figures for observational and MR results
# AUTHOR:   Jonathan Mak
# CREATED:	2022-02-18
# UPDATED: 	2023-02-21
# R VERSION: 4.1.3
#==============================================================================

setwd("P:/Frailty_omics/Frailty_omics_Research/Jonathan/Metabolomics_frailty")

### Required packages
library(dplyr)
library(TwoSampleMR)
library(ggplot2); library(patchwork); library(ggforce); library(cowplot) # Plotting graphs
#BiocManager::install("ComplexHeatmap")
#devtools::install_github("mattlee821/EpiViz/R_package")
library(EpiViz) # For circos plots of metabolites (see https://github.com/mattlee821/EpiViz)
library(ggforestplot) # For forest plots of metabolites


#============= CIRCOS PLOTS FOR METABOLITE-FRAILTY ASSOCIATIONS ===============

### Import UKB linear regression results
ukb_lm_univariate_fi_age_sex_168metabolomics <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_age_sex_168metabolomics.txt", sep=" ")
ukb_lm_univariate_fi_age_sex_32biomarkers <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_age_sex_32biomarkers.txt", sep=" ")
ukb_lm_univariate_fi_allcov_168metabolomics <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics.txt", sep=" ")
ukb_lm_univariate_fi_allcov_32biomarkers <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers.txt", sep=" ")
ukb_lm_univariate_fp_age_sex_168metabolomics <- read.delim("Output/Observational_results/ukb_lm_univariate_fp_age_sex_168metabolomics.txt", sep=" ")
ukb_lm_univariate_fp_age_sex_32biomarkers <- read.delim("Output/Observational_results/ukb_lm_univariate_fp_age_sex_32biomarkers.txt", sep=" ")
ukb_lm_univariate_fp_allcov_168metabolomics <- read.delim("Output/Observational_results/ukb_lm_univariate_fp_allcov_168metabolomics.txt", sep=" ")
ukb_lm_univariate_fp_allcov_32biomarkers <- read.delim("Output/Observational_results/ukb_lm_univariate_fp_allcov_32biomarkers.txt", sep=" ")
ukb_coef_lasso_fi_249metabolomics <- read.delim("Output/Observational_results/ukb_coef_lasso_fi_249metabolomics.txt", sep=" ")
ukb_coef_lasso_fi_32biomarker <- read.delim("Output/Observational_results/ukb_coef_lasso_fi_32biomarker.txt", sep=" ")

### Import labels of biomarkers from a prepared Excel file (information was based on description in UKB website)
ukb_labels <- readxl::read_excel("Documents/Biomarker_list.xlsx", sheet = "UKB")
ukb_labels <- ukb_labels[ukb_labels$Biomarker!="Oestradiol" & ukb_labels$Biomarker!="Rheumatoid factor",] # Remove oestradiol and rheumatoid factor (high missingness)

### Order for plotting metabolite categories
category_order <- data.frame(metabolite = c("Apolipoproteins","Amino acids","Branched-chain amino acids",
                                            "Aromatic amino acids","Fluid balance","Inflammation","Fatty acids", 
                                            "Glycolysis related metabolites","Ketone bodies","Total lipids",
                                            "Cholesterol","Free cholesterol","Cholesteryl esters","Phospholipids",
                                            "Triglycerides","Other lipids","Lipoprotein particle sizes",
                                            "Lipoprotein particle concentrations","Chylomicrons and extremely large VLDL",
                                            "Very large VLDL","Large VLDL","Medium VLDL","Small VLDL","Very small VLDL",
                                            "IDL","Large LDL","Medium LDL","Small LDL","Very large HDL","Large HDL",
                                            "Medium HDL", "Small HDL",
                                            "Bone and Joint","Cardiovascular","Diabetes","Hormone","Liver","Renal"),
                             order=c(1:38))

### Plots for age and sex-adjusted models of FI and FP in UKB
# Prepare dataset for plotting of FI models
circos_ukb_fi_age_sex <- rbind(ukb_lm_univariate_fi_age_sex_168metabolomics,ukb_lm_univariate_fi_age_sex_32biomarkers) %>% # FI models
  select(variable, beta.fi.age.sex, se.fi.age.sex, pvalue.fi.age.sex)  %>% 
  rename(beta=beta.fi.age.sex, se=se.fi.age.sex, p=pvalue.fi.age.sex)
row.names(circos_ukb_fi_age_sex) <- circos_ukb_fi_age_sex$variable
circos_ukb_fi_age_sex <- circos_ukb_fi_age_sex[ukb_labels[["Variable name"]][ukb_labels[["Variable name"]] %in% circos_ukb_fi_age_sex$variable],] # Re-order according to the labels
circos_ukb_fi_age_sex$Metabolite <- ukb_labels$Biomarker[ukb_labels$`Variable name` %in% circos_ukb_fi_age_sex$variable]
circos_ukb_fi_age_sex$Category <- ukb_labels$Category[ukb_labels$`Variable name` %in% circos_ukb_fi_age_sex$variable]
row.names(circos_ukb_fi_age_sex) <- NULL
circos_ukb_fi_age_sex$lower.bound <- circos_ukb_fi_age_sex$beta - 1.96*circos_ukb_fi_age_sex$se
circos_ukb_fi_age_sex$upper.bound <- circos_ukb_fi_age_sex$beta + 1.96*circos_ukb_fi_age_sex$se
circos_ukb_fi_age_sex$order <- category_order[match(circos_ukb_fi_age_sex$Category, category_order$metabolite),"order"]

# Prepare dataset for plotting of FP models
circos_ukb_fp_age_sex <- rbind(ukb_lm_univariate_fp_age_sex_168metabolomics,ukb_lm_univariate_fp_age_sex_32biomarkers) %>% # FP models
  select(variable, beta.fp.age.sex, se.fp.age.sex, pvalue.fp.age.sex)  %>% 
  rename(beta=beta.fp.age.sex, se=se.fp.age.sex, p=pvalue.fp.age.sex)
row.names(circos_ukb_fp_age_sex) <- circos_ukb_fp_age_sex$variable
circos_ukb_fp_age_sex <- circos_ukb_fp_age_sex[ukb_labels[["Variable name"]][ukb_labels[["Variable name"]] %in% circos_ukb_fp_age_sex$variable],] # Re-order according to the labels
circos_ukb_fp_age_sex$Metabolite <- ukb_labels$Biomarker[ukb_labels$`Variable name` %in% circos_ukb_fp_age_sex$variable]
circos_ukb_fp_age_sex$Category <- ukb_labels$Category[ukb_labels$`Variable name` %in% circos_ukb_fp_age_sex$variable]
row.names(circos_ukb_fp_age_sex) <- NULL
circos_ukb_fp_age_sex$lower.bound <- circos_ukb_fp_age_sex$beta - 1.96*circos_ukb_fp_age_sex$se
circos_ukb_fp_age_sex$upper.bound <- circos_ukb_fp_age_sex$beta + 1.96*circos_ukb_fp_age_sex$se
circos_ukb_fp_age_sex$order <- category_order[match(circos_ukb_fp_age_sex$Category, category_order$metabolite),"order"]

# Circos plots for NMR metabolomic biomarkers
#trace("circos_plot", edit=TRUE) 
pdf("Output/Observational_results/circos_plot_ukb_metabolomics_age_sex.pdf",
    width = 30, height = 30, pointsize = 40)
circos_plot(track_number = 2, # how many tracks do you want to plot
            track1_data = circos_ukb_fi_age_sex[1:168,] %>% mutate(Category=factor(Category)), # what is the dataframe for your first track
            track2_data = circos_ukb_fp_age_sex[1:168,] %>% mutate(Category=factor(Category)), # what is the dataframe for your second track
            track1_type = "points", # how do you want to plot your first track
            track2_type = "points", # how do you want to plot your second track
            label_column = 5, # whats is the column of your labels
            section_column = 6, # what is the column of your sections
            estimate_column = 2, # what is the column of your estimate (beta, OR etc.)
            pvalue_column = 4, # what is the column of your p-value
            pvalue_adjustment = .05/200, # what do you want your p-value adjustment to be = 0.05/X
            lower_ci = 7, # what is the column of your lower confidence interval
            upper_ci = 8, # what is the column of your upper confidence interval
            legend = F,
            track1_label = "Frailty index",
            track2_label = "Frailty phenotype",
            pvalue_label = "p>0.05/200",
            order = F,
            order_column = 9,
            circle_size = 15,
            colours = c("#E64B35FF","#3C5488FF","#00A087FF"))
dev.off()

# Circos plots for clinical biomarkers
#trace("circos_plot", edit=TRUE)
pdf("Output/Observational_results/circos_plot_ukb_clinical_biomarkers_age_sex.pdf",
    width = 20, height = 20, pointsize = 40)
circos_plot(track_number = 2, # how many tracks do you want to plot
            track1_data = circos_ukb_fi_age_sex[169:200,] %>% mutate(Category=factor(Category)), # what is the dataframe for your first track
            track2_data = circos_ukb_fp_age_sex[169:200,] %>% mutate(Category=factor(Category)), # what is the dataframe for your second track
            track1_type = "points", # how do you want to plot your first track
            track2_type = "points", # how do you want to plot your second track
            label_column = 5, # whats is the column of your labels
            section_column = 6, # what is the column of your sections
            estimate_column = 2, # what is the column of your estimate (beta, OR etc.)
            pvalue_column = 4, # what is the column of your p-value
            pvalue_adjustment = .05/200, # what do you want your p-value adjustment to be = 0.05/X
            lower_ci = 7, # what is the column of your lower confidence interval
            upper_ci = 8, # what is the column of your upper confidence interval
            legend = F,
            track1_label = "Frailty index",
            track2_label = "Frailty phenotype",
            pvalue_label = "p>0.05/200",
            order = F,
            order_column = 9,
            circle_size = 10,
            colours = c("#E64B35FF","#3C5488FF","#00A087FF"))
dev.off()


### Plots for fully-adjusted models of FI and FP in UKB
# Prepare dataset for plotting of FI models
circos_ukb_fi_allcov <- rbind(ukb_lm_univariate_fi_allcov_168metabolomics,ukb_lm_univariate_fi_allcov_32biomarkers) %>% # FI models
  select(variable, beta.fi.allcov, se.fi.allcov, pvalue.fi.allcov)  %>% 
  rename(beta=beta.fi.allcov, se=se.fi.allcov, p=pvalue.fi.allcov)
row.names(circos_ukb_fi_allcov) <- circos_ukb_fi_allcov$variable
circos_ukb_fi_allcov <- circos_ukb_fi_allcov[ukb_labels[["Variable name"]][ukb_labels[["Variable name"]] %in% circos_ukb_fi_allcov$variable],] # Re-order according to the labels
circos_ukb_fi_allcov$Metabolite <- ukb_labels$Biomarker[ukb_labels$`Variable name` %in% circos_ukb_fi_allcov$variable]
circos_ukb_fi_allcov$Category <- ukb_labels$Category[ukb_labels$`Variable name` %in% circos_ukb_fi_allcov$variable]
row.names(circos_ukb_fi_allcov) <- NULL
circos_ukb_fi_allcov$lower.bound <- circos_ukb_fi_allcov$beta - 1.96*circos_ukb_fi_allcov$se
circos_ukb_fi_allcov$upper.bound <- circos_ukb_fi_allcov$beta + 1.96*circos_ukb_fi_allcov$se
circos_ukb_fi_allcov$order <- category_order[match(circos_ukb_fi_allcov$Category, category_order$metabolite),"order"]

# Prepare dataset for plotting of FP models
circos_ukb_fp_allcov <- rbind(ukb_lm_univariate_fp_allcov_168metabolomics,ukb_lm_univariate_fp_allcov_32biomarkers) %>% # FP models
  select(variable, beta.fp.allcov, se.fp.allcov, pvalue.fp.allcov)  %>% 
  rename(beta=beta.fp.allcov, se=se.fp.allcov, p=pvalue.fp.allcov)
row.names(circos_ukb_fp_allcov) <- circos_ukb_fp_allcov$variable
circos_ukb_fp_allcov <- circos_ukb_fp_allcov[ukb_labels[["Variable name"]][ukb_labels[["Variable name"]] %in% circos_ukb_fp_allcov$variable],] # Re-order according to the labels
circos_ukb_fp_allcov$Metabolite <- ukb_labels$Biomarker[ukb_labels$`Variable name` %in% circos_ukb_fp_allcov$variable]
circos_ukb_fp_allcov$Category <- ukb_labels$Category[ukb_labels$`Variable name` %in% circos_ukb_fp_allcov$variable]
row.names(circos_ukb_fp_allcov) <- NULL
circos_ukb_fp_allcov$lower.bound <- circos_ukb_fp_allcov$beta - 1.96*circos_ukb_fp_allcov$se
circos_ukb_fp_allcov$upper.bound <- circos_ukb_fp_allcov$beta + 1.96*circos_ukb_fp_allcov$se
circos_ukb_fp_allcov$order <- category_order[match(circos_ukb_fp_allcov$Category, category_order$metabolite),"order"]

# Circos plots for NMR metabolomic biomarkers
#trace("circos_plot", edit=TRUE)
pdf("Output/Observational_results/circos_plot_ukb_metabolomics_allcov.pdf",
    width = 30, height = 30, pointsize = 40)
circos_plot(track_number = 2, # how many tracks do you want to plot
            track1_data = circos_ukb_fi_allcov[1:168,] %>% mutate(Category=factor(Category)), # what is the dataframe for your first track
            track2_data = circos_ukb_fp_allcov[1:168,] %>% mutate(Category=factor(Category)), # what is the dataframe for your second track
            track1_type = "points", # how do you want to plot your first track
            track2_type = "points", # how do you want to plot your second track
            label_column = 5, # whats is the column of your labels
            section_column = 6, # what is the column of your sections
            estimate_column = 2, # what is the column of your estimate (beta, OR etc.)
            pvalue_column = 4, # what is the column of your p-value
            pvalue_adjustment = .05/200, # what do you want your p-value adjustment to be = 0.05/X
            lower_ci = 7, # what is the column of your lower confidence interval
            upper_ci = 8, # what is the column of your upper confidence interval
            legend = F,
            track1_label = "Frailty index",
            track2_label = "Frailty phenotype",
            pvalue_label = "p>0.05/200",
            order = F,
            order_column = 9,
            circle_size = 15)
dev.off()

# Circos plots for clinical biomarkers
#trace("circos_plot", edit=TRUE) 
pdf("Output/Observational_results/circos_plot_ukb_clinical_biomarkers_allcov.pdf",
    width = 20, height = 20, pointsize = 40)
circos_plot(track_number = 2, # how many tracks do you want to plot
            track1_data = circos_ukb_fi_allcov[169:200,] %>% mutate(Category=factor(Category)), # what is the dataframe for your first track
            track2_data = circos_ukb_fp_allcov[169:200,] %>% mutate(Category=factor(Category)), # what is the dataframe for your second track
            track1_type = "points", # how do you want to plot your first track
            track2_type = "points", # how do you want to plot your second track
            label_column = 5, # whats is the column of your labels
            section_column = 6, # what is the column of your sections
            estimate_column = 2, # what is the column of your estimate (beta, OR etc.)
            pvalue_column = 4, # what is the column of your p-value
            pvalue_adjustment = .05/200, # what do you want your p-value adjustment to be = 0.05/X
            lower_ci = 7, # what is the column of your lower confidence interval
            upper_ci = 8, # what is the column of your upper confidence interval
            legend = F,
            track1_label = "Frailty index",
            track2_label = "Frailty phenotype",
            pvalue_label = "p>0.05/200",
            order = F,
            order_column = 9,
            circle_size = 10)
dev.off()




#============= FOREST PLOTS FOR METABOLITE-FRAILTY ASSOCIATIONS ==============

### Observatinal estimates from UKB, TwinGene, and Health 2000
combined_allcov <- read.delim("Output/Observational_results/combined_allcov_UKB_TwinGene_Health2000.txt", sep=" ")
ukb_lm_univariate_fp_allcov_168metabolomics <- read.delim("Output/Observational_results/ukb_lm_univariate_fp_allcov_168metabolomics.txt", sep=" ")
ukb_lm_univariate_fp_allcov_32biomarkers <- read.delim("Output/Observational_results/ukb_lm_univariate_fp_allcov_32biomarkers.txt", sep=" ")
ukb_lm_univariate_fi_allcov_168metabolomics_age60minus <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_age60minus.txt", sep=" ")
ukb_lm_univariate_fi_allcov_168metabolomics_age60plus <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_age60plus.txt", sep=" ")
ukb_lm_univariate_fi_allcov_168metabolomics_men <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_men.txt", sep=" ")
ukb_lm_univariate_fi_allcov_168metabolomics_women <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_women.txt", sep=" ")
ukb_lm_univariate_fi_allcov_168metabolomics_nonwhite <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_nonwhite.txt", sep=" ")
ukb_lm_univariate_fi_allcov_32biomarkers_age60minus <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_age60minus.txt", sep=" ")
ukb_lm_univariate_fi_allcov_32biomarkers_age60plus <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_age60plus.txt", sep=" ")
ukb_lm_univariate_fi_allcov_32biomarkers_men <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_men.txt", sep=" ")
ukb_lm_univariate_fi_allcov_32biomarkers_women <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_women.txt", sep=" ")
ukb_lm_univariate_fi_allcov_32biomarkers_nonwhite <- read.delim("Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_nonwhite.txt", sep=" ")

### Meta-analysis results
meta_TwinGene_Health2000_allcov  <- read.delim("Output/Observational_results/meta_TwinGene_Health2000_allcov.txt", sep=" ")
meta_UKB_TwinGene_Health2000_allcov  <- read.delim("Output/Observational_results/meta_UKB_TwinGene_Health2000_allcov.txt", sep=" ")

### MR results
mr_results_metabolomics <- read.delim("Output/MR_results/mr_results_metabolomics.txt", sep=" ")
mr_results_clinical_biomarkers <- read.delim("Output/MR_results/mr_results_clinical_biomarkers.txt", sep=" ")

### Import labels of biomarkers from a prepared Excel file (information was based on description in UKB website)
ukb_labels <- readxl::read_excel("Documents/Biomarker_list.xlsx", sheet = "UKB")
ukb_labels <- ukb_labels[ukb_labels$Biomarker!="Oestradiol" & ukb_labels$Biomarker!="Rheumatoid factor",] # Remove oestradiol and rheumatoid factor (high missingness)

### Biomarkers selected for MR analysis
selected_metabolomics <- c("Ala","Creatinine","Glucose","GlycA","HDL_size","IDL_CE",
                           "IDL_FC","IDL_P","LA","M_HDL_TG","M_LDL_CE","M_LDL_P",
                           "MUFA","Omega_6","Phe","Phosphatidylc","S_LDL_CE",
                           "S_LDL_L","S_LDL_P","S_LDL_PL","S_VLDL_TG","Sphingomyelins",
                           "VLDL_size","XL_HDL_FC","XS_VLDL_CE","XXL_VLDL_TG")
selected_clinical_biomarkers <- c("ALP_serum","APOB_serum","CHOL_serum","CRE_serum","CRP_serum",
                                  "CYS_serum","GGT_serum","HBA1C_RBC","IGF1_serum","K_urine",
                                  "LDL_serum","MALB_urine","PHOS_serum","SHBG_serum",
                                  "TBIL_serum","TES_serum","TRIG_serum","VITD_serum")


##### Forest plot combining observational & MR results ##### 

# FI models for metabolomics
results_metabolomics_fi_combined <- rbind(
  combined_allcov[combined_allcov$variable %in% selected_metabolomics,c("variable","beta","se","p","study")] %>%
    rename(b=beta, group=study),
    mr_results_metabolomics[mr_results_metabolomics$exposure %in% selected_metabolomics & mr_results_metabolomics$outcome=="FI" & mr_results_metabolomics$method=="Inverse variance weighted",c("exposure","b","se","pval","method")] %>%
    rename(variable=exposure,p=pval,group=method) %>%
    mutate(group=case_when(group=="Inverse variance weighted" ~ "IVW-MR")))
results_metabolomics_fi_combined$group <- factor(results_metabolomics_fi_combined$group, levels=c("IVW-MR","Health 2000","TwinGene","UKB"), labels=c("IVW-MR","Health 2000","TwinGene","UK Biobank"))
results_metabolomics_fi_combined$metabolite <- as.data.frame(ukb_labels)[match(results_metabolomics_fi_combined$variable, ukb_labels[["Variable name"]]),"Biomarker"] 
results_metabolomics_fi_combined$category <- as.data.frame(ukb_labels)[match(results_metabolomics_fi_combined$variable, ukb_labels[["Variable name"]]),"Category"] 
results_metabolomics_fi_combined$category <- factor(results_metabolomics_fi_combined$category, levels = c("Apolipoproteins","Amino acids","Branched-chain amino acids","Aromatic amino acids","Fluid balance","Inflammation","Fatty acids", "Glycolysis related metabolites","Ketone bodies","Total lipids","Cholesterol","Free cholesterol","Cholesteryl esters","Phospholipids","Triglycerides","Other lipids","Lipoprotein particle sizes","Lipoprotein particle concentrations","Chylomicrons and extremely large VLDL","Very large VLDL","Large VLDL","Medium VLDL","Small VLDL","Very small VLDL","IDL","Large LDL","Medium LDL","Small LDL","Very large HDL","Large HDL","Medium HDL", "Small HDL"))

# FP models for metabolomics
results_metabolomics_fp_combined <- rbind(
  ukb_lm_univariate_fp_allcov_168metabolomics[ukb_lm_univariate_fp_allcov_168metabolomics$variable %in% selected_metabolomics,c("variable","beta.fp.allcov","se.fp.allcov","pvalue.fp.allcov")] %>%
    rename(b=beta.fp.allcov,se=se.fp.allcov,p=pvalue.fp.allcov) %>%
    mutate(group="UKB"),
  mr_results_metabolomics[mr_results_metabolomics$exposure %in% selected_metabolomics & mr_results_metabolomics$outcome=="FP" & mr_results_metabolomics$method=="Inverse variance weighted",c("exposure","b","se","pval","method")] %>%
    rename(variable=exposure,p=pval,group=method) %>%
    mutate(group=case_when(group=="Inverse variance weighted" ~ "IVW-MR")))
results_metabolomics_fp_combined$group <- factor(results_metabolomics_fp_combined$group, levels=c("IVW-MR","UKB"), labels=c("IVW-MR","UK Biobank"))
results_metabolomics_fp_combined$metabolite <- as.data.frame(ukb_labels)[match(results_metabolomics_fp_combined$variable, ukb_labels[["Variable name"]]),"Biomarker"] 
results_metabolomics_fp_combined$category <- as.data.frame(ukb_labels)[match(results_metabolomics_fp_combined$variable, ukb_labels[["Variable name"]]),"Category"] 
results_metabolomics_fp_combined$category <- factor(results_metabolomics_fp_combined$category, levels = c("Apolipoproteins","Amino acids","Branched-chain amino acids","Aromatic amino acids","Fluid balance","Inflammation","Fatty acids", "Glycolysis related metabolites","Ketone bodies","Total lipids","Cholesterol","Free cholesterol","Cholesteryl esters","Phospholipids","Triglycerides","Other lipids","Lipoprotein particle sizes","Lipoprotein particle concentrations","Chylomicrons and extremely large VLDL","Very large VLDL","Large VLDL","Medium VLDL","Small VLDL","Very small VLDL","IDL","Large LDL","Medium LDL","Small LDL","Very large HDL","Large HDL","Medium HDL", "Small HDL"))

# FI models for clinical biomarkers
results_clinical_biomarker_fi_combined <- rbind(
  combined_allcov[combined_allcov$variable %in% selected_clinical_biomarkers,c("variable","beta","se","p","study")] %>%
    rename(b=beta, group=study),
  mr_results_clinical_biomarkers[mr_results_clinical_biomarkers$exposure %in% as.data.frame(ukb_labels)[match(selected_clinical_biomarkers, ukb_labels[["Variable name"]]),"Biomarker"]  & mr_results_clinical_biomarkers$outcome=="FI" & mr_results_clinical_biomarkers$method=="Inverse variance weighted",c("exposure","b","se","pval","method")] %>%
    rename(variable=exposure,p=pval,group=method) %>%
    mutate(group=case_when(group=="Inverse variance weighted" ~ "IVW-MR")))
results_clinical_biomarker_fi_combined$group <- factor(results_clinical_biomarker_fi_combined$group, levels=c("IVW-MR","Health 2000","TwinGene","UKB"), labels=c("IVW-MR","Health 2000","TwinGene","UK Biobank"))
results_clinical_biomarker_fi_combined$metabolite <- as.data.frame(ukb_labels[250:281,])[match(results_clinical_biomarker_fi_combined$variable, ukb_labels[250:281,][["Variable name"]]),"Biomarker"] 
results_clinical_biomarker_fi_combined$metabolite <- ifelse(is.na(results_clinical_biomarker_fi_combined$metabolite), as.vector(results_clinical_biomarker_fi_combined$variable), results_clinical_biomarker_fi_combined$metabolite)
results_clinical_biomarker_fi_combined$category <- as.data.frame(ukb_labels[250:281,])[match(results_clinical_biomarker_fi_combined$metabolite, ukb_labels[250:281,][["Biomarker"]]),"Category"] 
results_clinical_biomarker_fi_combined$category <- factor(results_clinical_biomarker_fi_combined$category)

# FP models for clinical biomarkers
results_clinical_biomarker_fp_combined <- rbind(
  ukb_lm_univariate_fp_allcov_32biomarkers[ukb_lm_univariate_fp_allcov_32biomarkers$variable %in% selected_clinical_biomarkers,c("variable","beta.fp.allcov","se.fp.allcov","pvalue.fp.allcov")] %>%
    rename(b=beta.fp.allcov,se=se.fp.allcov,p=pvalue.fp.allcov) %>%
    mutate(group="UKB"),
  mr_results_clinical_biomarkers[mr_results_clinical_biomarkers$exposure %in% as.data.frame(ukb_labels)[match(selected_clinical_biomarkers, ukb_labels[["Variable name"]]),"Biomarker"]  & mr_results_clinical_biomarkers$outcome=="FP" & mr_results_clinical_biomarkers$method=="Inverse variance weighted",c("exposure","b","se","pval","method")] %>%
    rename(variable=exposure,p=pval,group=method) %>%
    mutate(group=case_when(group=="Inverse variance weighted" ~ "IVW-MR")))
results_clinical_biomarker_fp_combined$group <- factor(results_clinical_biomarker_fp_combined$group, levels=c("IVW-MR","UKB"), labels=c("IVW-MR","UK Biobank"))
results_clinical_biomarker_fp_combined$metabolite <- as.data.frame(ukb_labels[250:281,])[match(results_clinical_biomarker_fp_combined$variable, ukb_labels[250:281,][["Variable name"]]),"Biomarker"] 
results_clinical_biomarker_fp_combined$metabolite <- ifelse(is.na(results_clinical_biomarker_fp_combined$metabolite), as.vector(results_clinical_biomarker_fp_combined$variable), results_clinical_biomarker_fp_combined$metabolite)
results_clinical_biomarker_fp_combined$category <- as.data.frame(ukb_labels[250:281,])[match(results_clinical_biomarker_fp_combined$metabolite, ukb_labels[250:281,][["Biomarker"]]),"Category"] 
results_clinical_biomarker_fp_combined$category <- factor(results_clinical_biomarker_fp_combined$category)

### Forest plot for FI models
cairo_pdf("Output/Forestplot_FI_selected_biomarkers.pdf", width = 14, height = 10, pointsize = 12)

(ggforestplot::forestplot(
  df = results_metabolomics_fi_combined[results_metabolomics_fi_combined$category %in% c("Apolipoproteins","Amino acids","Branched-chain amino acids", "Aromatic amino acids","Fluid balance","Inflammation","Fatty acids", "Glycolysis related metabolites","Ketone bodies","Total lipids","Cholesterol","Free cholesterol","Cholesteryl esters","Phospholipids","Triglycerides","Other lipids","Lipoprotein particle sizes","Lipoprotein particle concentrations"),],
  name = metabolite,
  estimate = b,
  se = se,
  pvalue = p,
  xlab = "Difference in FI (%) per SD increase",
  colour = group,
  shape = group,
  psignif = .011,
  xtickbreaks = c(-2,-1,0,1,2),
  xlim = c(-2,2)) +
    ggtitle("NMR metabolomic biomarkers") + 
    scale_color_manual(values = c("#925E9F","#DF8F44","#3C5488","#DC0000")) +
    scale_shape_manual(values = c(24,19,19,19)) +
    guides(colour = guide_legend(override.aes = list(size=2.5), reverse = T)) +
    labs(colour="Estimate", shape="Estimate") + 
    theme(legend.position="bottom",
          legend.title = element_text(size = 15, face = "bold"),
          legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=0.5),
          legend.text = element_text(size = 14 )) +
    ggforce::facet_col(facets = ~ category, scales = "free_y", space = "free" )) +
  
  (ggforestplot::forestplot(
    df = results_metabolomics_fi_combined[results_metabolomics_fi_combined$category %in% c("Chylomicrons and extremely large VLDL","Very large VLDL","Large VLDL","Medium VLDL","Small VLDL","Very small VLDL","IDL","Large LDL","Medium LDL","Small LDL","Very large HDL","Large HDL","Medium HDL", "Small HDL"),],
    name = metabolite,
    estimate = b,
    se = se,
    pvalue = p,
    xlab = "Difference in FI (%) per SD increase",
    colour = group,
    shape = group,
    psignif = .011,
    xtickbreaks = c(-2,-1,0,1,2),
    xlim = c(-2,2)) +
     ggtitle("") + 
     scale_color_manual(values = c("#925E9F","#DF8F44","#3C5488","#DC0000")) +
     scale_shape_manual(values = c(24,19,19,19)) +
     guides(colour = guide_legend(override.aes = list(size=2.5), reverse = T)) +
     labs(colour="Estimate", shape="Estimate") + 
     theme(legend.position="bottom",
           legend.title = element_text(size = 15, face = "bold"),
           legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=0.5),
           legend.text = element_text(size = 14)) +
     ggforce::facet_col(facets = ~ category, scales = "free_y", space = "free" )) + 
  
  (ggforestplot::forestplot(
    df = results_clinical_biomarker_fi_combined,
    name = metabolite,
    estimate = b,
    se = se,
    pvalue = p,
    xlab = "Difference in FI (%) per SD increase",
    colour = group,
    shape = group,
    psignif = .011,
    xtickbreaks = c(-2,-1,0,1,2),
    xlim = c(-2,2)) +
     ggtitle("Clinical biomarkers") + 
     scale_color_manual(values = c("#925E9F","#DF8F44","#3C5488","#DC0000")) +
     scale_shape_manual(values = c(24,19,19,19)) +
     guides(colour = guide_legend(override.aes = list(size=2.5), reverse = T)) +
     labs(colour="Estimate", shape="Estimate") + 
     theme(legend.position="bottom",
           legend.title = element_text(size = 15, face = "bold"),
           legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=0.5),
           legend.text = element_text(size = 14)) +
     ggforce::facet_col(facets = ~ category, scales = "free_y", space = "free" )) +
  
  plot_layout(guides = "collect") &  theme(legend.position="bottom")

dev.off()


### Forest plot for FP models
cairo_pdf("Output/Forestplot_FP_selected_biomarkers.pdf", width = 14, height = 10, pointsize = 12)

(ggforestplot::forestplot(
  df = results_metabolomics_fp_combined[results_metabolomics_fp_combined$category %in% c("Apolipoproteins","Amino acids","Branched-chain amino acids", "Aromatic amino acids","Fluid balance","Inflammation","Fatty acids", "Glycolysis related metabolites","Ketone bodies","Total lipids","Cholesterol","Free cholesterol","Cholesteryl esters","Phospholipids","Triglycerides","Other lipids","Lipoprotein particle sizes","Lipoprotein particle concentrations"),],
  name = metabolite,
  estimate = b,
  se = se,
  pvalue = p,
  xlab = "Difference in FP score per SD increase",
  colour = group,
  shape = group,
  psignif = .011,
  xlim = c(-.2,.2)) +
    ggtitle("NMR metabolomic biomarkers") + 
    scale_color_manual(values = c("#925E9F","#DC0000")) +
    scale_shape_manual(values = c(24,19)) +
    guides(colour = guide_legend(override.aes = list(size=2.5), reverse = T)) +
    labs(colour="Estimate", shape="Estimate") + 
    theme(legend.position="bottom",
          legend.title = element_text(size = 15, face = "bold"),
          legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=0.5),
          legend.text = element_text(size = 14)) +
    ggforce::facet_col(facets = ~ category, scales = "free_y", space = "free" )) +
  
  (ggforestplot::forestplot(
    df = results_metabolomics_fp_combined[results_metabolomics_fp_combined$category %in% c("Chylomicrons and extremely large VLDL","Very large VLDL","Large VLDL","Medium VLDL","Small VLDL","Very small VLDL","IDL","Large LDL","Medium LDL","Small LDL","Very large HDL","Large HDL","Medium HDL", "Small HDL"),],
    name = metabolite,
    estimate = b,
    se = se,
    pvalue = p,
    xlab = "Difference in FP score per SD increase",
    colour = group,
    shape = group,
    psignif = .011,
    xlim = c(-.2,.2)) +
     ggtitle("") + 
     scale_color_manual(values = c("#925E9F","#DC0000")) +
     scale_shape_manual(values = c(24,19)) +
     guides(colour = guide_legend(override.aes = list(size=2.5), reverse = T)) +
     labs(colour="Estimate", shape="Estimate") + 
     theme(legend.position="bottom",
           legend.title = element_text(size = 15, face = "bold"),
           legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=0.5),
           legend.text = element_text(size = 14)) +
     ggforce::facet_col(facets = ~ category, scales = "free_y", space = "free" )) + 
  
  (ggforestplot::forestplot(
    df = results_clinical_biomarker_fp_combined,
    name = metabolite,
    estimate = b,
    se = se,
    pvalue = p,
    xlab = "Difference in FP score per SD increase",
    colour = group,
    shape = group,
    psignif = .011,
    xlim = c(-.2,.2)) +
     ggtitle("Clinical biomarkers") + 
     scale_color_manual(values = c("#925E9F","#DC0000")) +
     scale_shape_manual(values = c(24,19)) +
     guides(colour = guide_legend(override.aes = list(size=2.5), reverse = T)) +
     labs(colour="Estimate", shape="Estimate") + 
     theme(legend.position="bottom",
           legend.title = element_text(size = 15, face = "bold"),
           legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=0.5),
           legend.text = element_text(size = 14)) +
     ggforce::facet_col(facets = ~ category, scales = "free_y", space = "free" )) +
  
  plot_layout(guides = "collect") &  theme(legend.position="bottom")

dev.off()


##### Forest plot for subgroup analysis #####

# FI models for metabolomics
results_metabolomics_fi_subgroup <- rbind(
  ukb_lm_univariate_fi_allcov_168metabolomics_age60minus[ukb_lm_univariate_fi_allcov_168metabolomics_age60minus$variable %in% selected_metabolomics,c("variable","beta.fi.allcov.age60minus","se.fi.allcov.age60minus","pvalue.fi.allcov.age60minus")] %>%
    rename(b=beta.fi.allcov.age60minus,se=se.fi.allcov.age60minus,p=pvalue.fi.allcov.age60minus) %>%
    mutate(group="Age <60 years"),
  ukb_lm_univariate_fi_allcov_168metabolomics_age60plus[ukb_lm_univariate_fi_allcov_168metabolomics_age60plus$variable %in% selected_metabolomics,c("variable","beta.fi.allcov.age60plus","se.fi.allcov.age60plus","pvalue.fi.allcov.age60plus")] %>%
    rename(b=beta.fi.allcov.age60plus,se=se.fi.allcov.age60plus,p=pvalue.fi.allcov.age60plus) %>%
    mutate(group="Age ≥60 years"),
  ukb_lm_univariate_fi_allcov_168metabolomics_women[ukb_lm_univariate_fi_allcov_168metabolomics_women$variable %in% selected_metabolomics,c("variable","beta.fi.allcov.women","se.fi.allcov.women","pvalue.fi.allcov.women")] %>%
    rename(b=beta.fi.allcov.women,se=se.fi.allcov.women,p=pvalue.fi.allcov.women) %>%
    mutate(group="Women"),
  ukb_lm_univariate_fi_allcov_168metabolomics_men[ukb_lm_univariate_fi_allcov_168metabolomics_men$variable %in% selected_metabolomics,c("variable","beta.fi.allcov.men","se.fi.allcov.men","pvalue.fi.allcov.men")] %>%
    rename(b=beta.fi.allcov.men,se=se.fi.allcov.men,p=pvalue.fi.allcov.men) %>%
    mutate(group="Men"),
  ukb_lm_univariate_fi_allcov_168metabolomics_nonwhite[ukb_lm_univariate_fi_allcov_168metabolomics_nonwhite$variable %in% selected_metabolomics,c("variable","beta.fi.allcov.nonwhite","se.fi.allcov.nonwhite","pvalue.fi.allcov.nonwhite")] %>%
    rename(b=beta.fi.allcov.nonwhite,se=se.fi.allcov.nonwhite,p=pvalue.fi.allcov.nonwhite) %>%
    mutate(group="Non-white ethnicity"))
results_metabolomics_fi_subgroup$group <- factor(results_metabolomics_fi_subgroup$group, levels=c("Non-white ethnicity","Women","Men","Age ≥60 years","Age <60 years"))
results_metabolomics_fi_subgroup$metabolite <- as.data.frame(ukb_labels)[match(results_metabolomics_fi_subgroup$variable, ukb_labels[["Variable name"]]),"Biomarker"] 
results_metabolomics_fi_subgroup$category <- as.data.frame(ukb_labels)[match(results_metabolomics_fi_subgroup$variable, ukb_labels[["Variable name"]]),"Category"] 
results_metabolomics_fi_subgroup$category <- factor(results_metabolomics_fi_subgroup$category, levels = c("Apolipoproteins","Amino acids","Branched-chain amino acids","Aromatic amino acids","Fluid balance","Inflammation","Fatty acids", "Glycolysis related metabolites","Ketone bodies","Total lipids","Cholesterol","Free cholesterol","Cholesteryl esters","Phospholipids","Triglycerides","Other lipids","Lipoprotein particle sizes","Lipoprotein particle concentrations","Chylomicrons and extremely large VLDL","Very large VLDL","Large VLDL","Medium VLDL","Small VLDL","Very small VLDL","IDL","Large LDL","Medium LDL","Small LDL","Very large HDL","Large HDL","Medium HDL", "Small HDL"))

# FI models for clinical_biomarkers
results_clinical_biomarker_fi_subgroup <- rbind(
  ukb_lm_univariate_fi_allcov_32biomarkers_age60minus[ukb_lm_univariate_fi_allcov_32biomarkers_age60minus$variable %in% selected_clinical_biomarkers,c("variable","beta.fi.allcov.age60minus","se.fi.allcov.age60minus","pvalue.fi.allcov.age60minus")] %>%
    rename(b=beta.fi.allcov.age60minus,se=se.fi.allcov.age60minus,p=pvalue.fi.allcov.age60minus) %>%
    mutate(group="Age <60 years"),
  ukb_lm_univariate_fi_allcov_32biomarkers_age60plus[ukb_lm_univariate_fi_allcov_32biomarkers_age60plus$variable %in% selected_clinical_biomarkers,c("variable","beta.fi.allcov.age60plus","se.fi.allcov.age60plus","pvalue.fi.allcov.age60plus")] %>%
    rename(b=beta.fi.allcov.age60plus,se=se.fi.allcov.age60plus,p=pvalue.fi.allcov.age60plus) %>%
    mutate(group="Age ≥60 years"),
  ukb_lm_univariate_fi_allcov_32biomarkers_women[ukb_lm_univariate_fi_allcov_32biomarkers_women$variable %in% selected_clinical_biomarkers,c("variable","beta.fi.allcov.women","se.fi.allcov.women","pvalue.fi.allcov.women")] %>%
    rename(b=beta.fi.allcov.women,se=se.fi.allcov.women,p=pvalue.fi.allcov.women) %>%
    mutate(group="Women"),
  ukb_lm_univariate_fi_allcov_32biomarkers_men[ukb_lm_univariate_fi_allcov_32biomarkers_men$variable %in% selected_clinical_biomarkers,c("variable","beta.fi.allcov.men","se.fi.allcov.men","pvalue.fi.allcov.men")] %>%
    rename(b=beta.fi.allcov.men,se=se.fi.allcov.men,p=pvalue.fi.allcov.men) %>%
    mutate(group="Men"),
  ukb_lm_univariate_fi_allcov_32biomarkers_nonwhite[ukb_lm_univariate_fi_allcov_32biomarkers_nonwhite$variable %in% selected_clinical_biomarkers,c("variable","beta.fi.allcov.nonwhite","se.fi.allcov.nonwhite","pvalue.fi.allcov.nonwhite")] %>%
    rename(b=beta.fi.allcov.nonwhite,se=se.fi.allcov.nonwhite,p=pvalue.fi.allcov.nonwhite) %>%
    mutate(group="Non-white ethnicity"))
results_clinical_biomarker_fi_subgroup$group <- factor(results_clinical_biomarker_fi_subgroup$group, levels=c("Non-white ethnicity","Women","Men","Age ≥60 years","Age <60 years"))
results_clinical_biomarker_fi_subgroup$metabolite <- as.data.frame(ukb_labels[250:281,])[match(results_clinical_biomarker_fi_subgroup$variable, ukb_labels[250:281,][["Variable name"]]),"Biomarker"] 
results_clinical_biomarker_fi_subgroup$metabolite <- ifelse(is.na(results_clinical_biomarker_fi_subgroup$metabolite), as.vector(results_clinical_biomarker_fi_subgroup$variable), results_clinical_biomarker_fi_subgroup$metabolite)
results_clinical_biomarker_fi_subgroup$category <- as.data.frame(ukb_labels[250:281,])[match(results_clinical_biomarker_fi_subgroup$metabolite, ukb_labels[250:281,][["Biomarker"]]),"Category"] 
results_clinical_biomarker_fi_subgroup$category <- factor(results_clinical_biomarker_fi_subgroup$category)

### Forest plot for FI models
cairo_pdf("Output/Forestplot_FI_subgroup.pdf", width = 14, height = 10, pointsize = 12)

(ggforestplot::forestplot(
  df = results_metabolomics_fi_subgroup[results_metabolomics_fi_subgroup$category %in% c("Apolipoproteins","Amino acids","Branched-chain amino acids", "Aromatic amino acids","Fluid balance","Inflammation","Fatty acids", "Glycolysis related metabolites","Ketone bodies","Total lipids","Cholesterol","Free cholesterol","Cholesteryl esters","Phospholipids","Triglycerides","Other lipids","Lipoprotein particle sizes","Lipoprotein particle concentrations"),],
  name = metabolite,
  estimate = b,
  se = se,
  pvalue = p,
  xlab = "Difference in FI (%) per SD increase",
  colour = group,
  shape = group,
  psignif = .05/200,
  xtickbreaks = c(-2,-1,0,1,2),
  xlim = c(-2,2)) +
    ggtitle("NMR metabolomic biomarkers") + 
    #ggsci::scale_color_nejm(palette = c("default"), alpha = 1) +
    scale_color_manual(values = ggpubfigs::friendly_pal("vibrant_seven")) +
    scale_shape_manual(values = c(23,22,22,21,21)) +
    guides(colour = guide_legend(override.aes = list(size=2.5), reverse = T)) +
    labs(colour="UKB subgroup", shape="UKB subgroup") + 
    theme(legend.position="bottom",
          legend.title = element_text(size = 15, face = "bold"),
          legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=0.5),
          legend.text = element_text(size = 14 )) +
    ggforce::facet_col(facets = ~ category, scales = "free_y", space = "free" )) +
  
  (ggforestplot::forestplot(
    df = results_metabolomics_fi_subgroup[results_metabolomics_fi_subgroup$category %in% c("Chylomicrons and extremely large VLDL","Very large VLDL","Large VLDL","Medium VLDL","Small VLDL","Very small VLDL","IDL","Large LDL","Medium LDL","Small LDL","Very large HDL","Large HDL","Medium HDL", "Small HDL"),],
    name = metabolite,
    estimate = b,
    se = se,
    pvalue = p,
    xlab = "Difference in FI (%) per SD increase",
    colour = group,
    shape = group,
    psignif = .05/200,
    xtickbreaks = c(-2,-1,0,1,2),
    xlim = c(-2,2)) +
     ggtitle("") + 
     #ggsci::scale_color_nejm(palette = c("default"), alpha = 1) +
     scale_color_manual(values = ggpubfigs::friendly_pal("vibrant_seven")) +
     scale_shape_manual(values = c(23,22,22,21,21)) +
     guides(colour = guide_legend(override.aes = list(size=2.5), reverse = T)) +
     labs(colour="UKB subgroup", shape="UKB subgroup") + 
     theme(legend.position="bottom",
           legend.title = element_text(size = 15, face = "bold"),
           legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=0.5),
           legend.text = element_text(size = 14)) +
     ggforce::facet_col(facets = ~ category, scales = "free_y", space = "free" )) + 
  
  (ggforestplot::forestplot(
    df = results_clinical_biomarker_fi_subgroup,
    name = metabolite,
    estimate = b,
    se = se,
    pvalue = p,
    xlab = "Difference in FI (%) per SD increase",
    colour = group,
    shape = group,
    psignif = .05/200,
    xtickbreaks = c(-3,-2,-1,0,1,2),
    xlim = c(-3.5,2)) +
     ggtitle("Clinical biomarkers") + 
     #ggsci::scale_color_nejm(palette = c("default"), alpha = 1) +
     scale_color_manual(values = ggpubfigs::friendly_pal("vibrant_seven")) +
     scale_shape_manual(values = c(23,22,22,21,21)) +
     guides(colour = guide_legend(override.aes = list(size=2.5), reverse = T)) +
     labs(colour="UKB subgroup", shape="UKB subgroup") + 
     theme(legend.position="bottom",
           legend.title = element_text(size = 15, face = "bold"),
           legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=0.5),
           legend.text = element_text(size = 14)) +
     ggforce::facet_col(facets = ~ category, scales = "free_y", space = "free" )) +
  
  plot_layout(guides = "collect", widths = c(1,1,1.1)) &  theme(legend.position="bottom")

dev.off()



##### Forest plots for additional analysis for GlycA-FI association ##### 

# Subgroup analysis by LDL and CRP
ukb_lm_GlycA_fi_subgroup <- rbind(
  ukb_lm_GlycA_fi_ldl %>% mutate(subgroup="Categories of LDL cholesterol") %>% rename(b="Estimate",se="Std. Error",p="Pr(>|t|)"),
  ukb_lm_GlycA_fi_crp %>% mutate(subgroup="Categories of C-reactive protein") %>% rename(b="Estimate",se="Std. Error",p="Pr(>|t|)")
)
plot_ukb_lm_GlycA_fi_subgroup <- ggforestplot::forestplot(
  df = ukb_lm_GlycA_fi_subgroup %>%
    mutate( model=case_when(
      model=="<2.59 mmol/L" ~ "<2.59 mmol/L\n(n = 5,783)",
      model=="2.60–3.36 mmol/L" ~ "2.60–3.36 mmol/L\n(n = 22,380)",
      model=="3.37–4.14 mmol/L" ~ "3.37–4.14 mmol/L\n(n = 31,600)",
      model=="4.15–4.92 mmol/L" ~ "4.15–4.92 mmol/L\n(n = 18,672)",
      model=="≥4.92 mmol/L" ~ "≥4.92 mmol/L\n(n = 7,627)",
      model=="<1 mg/L" ~ "<1 mg/L\n(n = 34,845)",
      model=="1–3 mg/L" ~ "1–2.9 mg/L\n(n = 32,911)",
      model=="3–10 mg/L" ~ "3–9.9 mg/L\n(n = 15,473)",
      model=="≥10 mg/L" ~ "≥10 mg/L\n(n = 3,249)")),
  name = model,
  estimate = b,
  se = se,
  pvalue = p,
  xlab = "Difference in FI (%) per SD increase in GlycA",
  psignif = .05/200) +
  ggforce::facet_col(facets = ~ subgroup, scales = "free_y", space = "free" ) +
  labs(title="Subgroup analysis for GlycA and FI in UK Biobank")  +
  theme(plot.title = element_text(size = 9.5, face = "bold", hjust = .6),
        plot.title.position = "plot",
        axis.text=element_text(size=7, color="black"),
        axis.title=element_text(size=7, color="black"),
        strip.text.x = element_text(size = 7.5, face = "plain"))
cairo_pdf("Output/Observational_results/Subgroup_ukb_GlycA_FI.pdf", width = 4.5, height = 4, pointsize = 10)
plot_ukb_lm_GlycA_fi_subgroup
dev.off()

# Co-twin control
plot_cotwin_GlycA_fi <- ggplot(
  cotwin_results_GlycA_fi %>%
    rename(b="Estimate",se="Std. Error",p="Pr(>|z|)") %>%
    mutate(twins=factor(case_when(model=="Population-level estimate (full sample)"~"Full sample",
                                  model=="Population-level estimate (DZ twins)"~"DZ twins",
                                  model=="Population-level estimate (MZ twins)"~"MZ twins",
                                  model=="Within-twin-pair estimate (DZ twins)"~"DZ twins",
                                  model=="Within-twin-pair estimate (MZ twins)"~"MZ twins"),
                        levels = c("Full sample","DZ twins","MZ twins")),
           level=factor(case_when(model=="Population-level estimate (full sample)"~"Population-level estimate",
                                  model=="Population-level estimate (DZ twins)"~"Population-level estimate",
                                  model=="Population-level estimate (MZ twins)"~"Population-level estimate",
                                  model=="Within-twin-pair estimate (DZ twins)"~"Within-twin-pair estimate",
                                  model=="Within-twin-pair estimate (MZ twins)"~"Within-twin-pair estimate"),
                        levels = c("Population-level estimate","Within-twin-pair estimate")),
           model=factor(model, levels = c("Population-level estimate (full sample)",
                                          "Population-level estimate (DZ twins)",
                                          "Population-level estimate (MZ twins)",
                                          "Within-twin-pair estimate (DZ twins)",
                                          "Within-twin-pair estimate (MZ twins)"))), 
  aes(x=level, y=b, fill=twins)) + 
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(label = round(b,2)), vjust = -0.2, hjust=-.1, position = position_dodge2(width = 0.9, preserve = "single"), size=3) +
  geom_pointrange(aes(ymin=b-1.96*se, ymax=b+1.96*se), position = position_dodge2(width = 0.9, preserve = "single"), size=.15, linewidth=.3) +
  geom_hline(yintercept=0, color="black", size=.3) +
  scale_y_continuous(limits = c(-.335,1)) +
  scale_fill_manual(values = c("#7E6148E5","#3C5488E5","#E64B35CC"), 
                    guide = guide_legend(override.aes = list(shape=NA, color="transparent"))) +
  labs(title="Co-twin control analysis for GlycA and FI in TwinGene", y="Difference in FI (%) per SD increase in GlycA") +
  theme_bw() +
  theme(plot.title = element_text(size = 9.5, face = "bold"),
        legend.position = c(.87,.91),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.45, "cm"),
        legend.text = element_text(size = 7),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.15, color="grey"),
        panel.grid.minor.y = element_line(size=.15, color="grey"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8, color="black"),
        axis.text=element_text(size=8, color="black"),
        axis.ticks.x=element_blank())
cairo_pdf("Output/Observational_results/Cotwin_TwinGene_GlycA_FI.pdf", width = 5, height = 4, pointsize = 10)
plot_cotwin_GlycA_fi
dev.off()

# Combine the two graphs
cairo_pdf("Output/Observational_results/GlycA_FI_combined.pdf", width = 8, height = 4, pointsize = 10)
plot_grid(plot_ukb_lm_GlycA_fi_subgroup, plot_cotwin_GlycA_fi, labels = c("a","b"), rel_widths = c(1, 1.05))
dev.off()



##### Forest plots for additional analysis for Creatinine-FI association ##### 

# Subgroup analysis by CKD
plot_ukb_lm_Creatinine_fi_subgroup <- ggforestplot::forestplot(
  df = ukb_lm_creatinine_fi_ckd %>% 
    mutate(subgroup="Chronic kidney disease",
           model=case_when(
             model== "No CKD" ~ "No\n(n = 84,599)",
             model== "CKD" ~ "Yes\n(n = 5,974)"
           )) %>% 
    rename(b="Estimate",se="Std. Error",p="Pr(>|t|)"),
  name = model,
  estimate = b,
  se = se,
  pvalue = p,
  xlab = "Difference in FI (%) per SD increase in creatinine",
  psignif = .05/200) +
  ggforce::facet_col(facets = ~ subgroup, scales = "free_y", space = "free" ) +
  labs(title="Subgroup analysis for creatinine and FI in UK Biobank")  +
  theme(plot.title = element_text(size = 9.5, face = "bold", hjust = .6),
        plot.title.position = "plot",
        axis.text=element_text(size=7, color="black"),
        axis.title=element_text(size=7, color="black"),
        strip.text.x = element_text(size = 7.5, face = "plain"))
cairo_pdf("Output/Observational_results/Subgroup_ukb_Creatinine_FI.pdf", width = 4.5, height = 4, pointsize = 10)
plot_ukb_lm_Creatinine_fi_subgroup
dev.off()

# Co-twin control
plot_cotwin_Creatinine_fi <- ggplot(
  cotwin_results_Creatinine_fi %>%
    rename(b="Estimate",se="Std. Error",p="Pr(>|z|)") %>%
    mutate(twins=factor(case_when(model=="Population-level estimate (full sample)"~"Full sample",
                                  model=="Population-level estimate (DZ twins)"~"DZ twins",
                                  model=="Population-level estimate (MZ twins)"~"MZ twins",
                                  model=="Within-twin-pair estimate (DZ twins)"~"DZ twins",
                                  model=="Within-twin-pair estimate (MZ twins)"~"MZ twins"),
                        levels = c("Full sample","DZ twins","MZ twins")),
           level=factor(case_when(model=="Population-level estimate (full sample)"~"Population-level estimate",
                                  model=="Population-level estimate (DZ twins)"~"Population-level estimate",
                                  model=="Population-level estimate (MZ twins)"~"Population-level estimate",
                                  model=="Within-twin-pair estimate (DZ twins)"~"Within-twin-pair estimate",
                                  model=="Within-twin-pair estimate (MZ twins)"~"Within-twin-pair estimate"),
                        levels = c("Population-level estimate","Within-twin-pair estimate")),
           model=factor(model, levels = c("Population-level estimate (full sample)",
                                          "Population-level estimate (DZ twins)",
                                          "Population-level estimate (MZ twins)",
                                          "Within-twin-pair estimate (DZ twins)",
                                          "Within-twin-pair estimate (MZ twins)"))), 
  aes(x=level, y=b, fill=twins)) + 
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(label = round(b,2)), vjust = -0.2, hjust=-.1, position = position_dodge2(width = 0.9, preserve = "single"), size=3) +
  geom_pointrange(aes(ymin=b-1.96*se, ymax=b+1.96*se), position = position_dodge2(width = 0.9, preserve = "single"), size=.15, linewidth=.3) +
  geom_hline(yintercept=0, color="black", size=.3) +
  scale_y_continuous(limits = c(-.39,1)) +
  scale_fill_manual(values = c("#7E6148E5","#3C5488E5","#E64B35CC"), 
                    guide = guide_legend(override.aes = list(shape=NA, color="transparent"))) +
  labs(title="Co-twin control analysis for creatinine and FI in TwinGene", y="Difference in FI (%) per SD increase in creatinine") +
  theme_bw() +
  theme(plot.title = element_text(size = 9.5, face = "bold"),
        legend.position = c(.87,.91),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.45, "cm"),
        legend.text = element_text(size = 7),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.15, color="grey"),
        panel.grid.minor.y = element_line(size=.15, color="grey"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8, color="black"),
        axis.text=element_text(size=8, color="black"),
        axis.ticks.x=element_blank())
cairo_pdf("Output/Observational_results/Cotwin_TwinGene_Creatinine_FI.pdf", width = 4.5, height = 4, pointsize = 10)
plot_cotwin_Creatinine_fi
dev.off()

# Combine the two graphs
cairo_pdf("Output/Observational_results/Creatinine_FI_combined.pdf", width = 8, height = 4, pointsize = 10)
plot_grid(plot_ukb_lm_Creatinine_fi_subgroup, plot_cotwin_Creatinine_fi, labels = c("a","b"), rel_widths = c(1, 1.05))
dev.off()


# Combine graphs for GlycA & creatinine
cairo_pdf("Output/Observational_results/GlycA_creatinine_FI_combined.pdf", width = 8, height = 7.5, pointsize = 10)
plot_grid(plot_ukb_lm_GlycA_fi_subgroup, plot_cotwin_GlycA_fi, plot_ukb_lm_Creatinine_fi_subgroup, plot_cotwin_Creatinine_fi, labels = c("a","b","c","d"), rel_widths = c(1, 1.05))
dev.off()


#========================== PLOTS FOR MR ANALYSIS ============================

### MR results
mr_results_metabolomics <- read.delim("Output/MR_results/mr_results_metabolomics.txt", sep=" ")
mr_results_clinical_biomarkers <- read.delim("Output/MR_results/mr_results_clinical_biomarkers.txt", sep=" ")

##### Scatter plot and forest plot for MR estimates of GlycA-FI association ##### 

# Outlier SNPs identified by MR-PRESSO
GlycA_FI_outlier <- mr_results_metabolomics[mr_results_metabolomics$exposure=="GlycA" &
                                              mr_results_metabolomics$outcome=="FI" &
                                              mr_results_metabolomics$method=="MR-PRESSO","outlier_snps"] %>%
  strsplit(., split = ", ") %>% unlist()

mr_GlycA_fi_single <- rbind(
  # Association between each GlycA-associated SNP and FI
  mr_single[mr_single$exposure=="GlycA" & mr_single$outcome=="FI" & !mr_single$SNP %in% c("All - Inverse variance weighted","All - MR Egger"),] %>%
    mutate(SNP = case_when(
      SNP %in% GlycA_FI_outlier ~ paste0(SNP, "*"),
      !SNP %in% GlycA_FI_outlier ~ SNP),
      method="") %>%
    select(exposure, outcome, SNP, method, b, se, p),
  # Main MR results for GlycA-FI association
  mr_summary_metabolomics[mr_summary_metabolomics$exposure=="GlycA"&mr_summary_metabolomics$outcome=="FI",]  %>%
    rename(p=pval) %>% 
    mutate(SNP=case_when(
      method=="Inverse variance weighted" ~ "IVW",
      method=="MR Egger" ~ "MR Egger",
      method=="Weighted median" ~ "Weighted median",
      method=="Weighted mode" ~ "Weighted mode",
      method=="MR-PRESSO" ~ "MR-PRESSO"
    )) %>%
    select(exposure, outcome, SNP, method, b, se, p)
)
mr_GlycA_fi_single$up <- mr_GlycA_fi_single$b + 1.96 * mr_GlycA_fi_single$se # confidence interval
mr_GlycA_fi_single$lo <- mr_GlycA_fi_single$b - 1.96 * mr_GlycA_fi_single$se
nom <- mr_GlycA_fi_single$SNP[grepl("rs", mr_GlycA_fi_single$SNP)]
nom <- nom[order(mr_GlycA_fi_single$b)] # sort effect size
mr_GlycA_fi_single <- rbind(mr_GlycA_fi_single, data.frame(exposure="GlycA",outcome="FI",SNP="",method="",b=NA,se=NA,p=NA,up=NA,lo=NA))
mr_GlycA_fi_single$SNP <- factor(mr_GlycA_fi_single$SNP, levels = c("MR-PRESSO", "Weighted mode", "Weighted median", "MR Egger", "IVW", "", nom))

scatterplot_GlycA_FI <- ggplot(data=mr_data_metabolomics[mr_data_metabolomics$exposure=="GlycA" & mr_data_metabolomics$outcome=="FI" & mr_data_metabolomics$mr_keep,] %>%
                                 mutate(beta.exposure_flipped=case_when(beta.exposure < 0 ~ beta.exposure * -1,
                                                                        beta.exposure >= 0 ~ beta.exposure),
                                        beta.outcome_flipped=case_when(beta.exposure < 0 ~ beta.outcome * -1,
                                                                       beta.exposure >= 0 ~ beta.outcome),
                                        outlier=case_when(SNP%in%GlycA_FI_outlier ~ 1,
                                                          !SNP%in%GlycA_FI_outlier ~ 0)),
                               aes(x=beta.exposure_flipped, y=beta.outcome_flipped)) +
  geom_errorbar(aes(ymin=beta.outcome_flipped-se.outcome, ymax=beta.outcome_flipped+se.outcome), colour="grey", width=0) +
  geom_errorbarh(aes(xmin=beta.exposure_flipped-se.exposure, xmax=beta.exposure_flipped+se.exposure), colour="grey", height=0) +
  geom_point(aes(shape=factor(outlier))) +
  geom_abline(data=mr_GlycA_fi_single[!grepl("rs", mr_GlycA_fi_single$SNP) & mr_GlycA_fi_single$SNP!="",] %>% 
                mutate(a = case_when(
                  SNP!="MR Egger" ~ 0,
                  SNP=="MR Egger" ~ mr_plt[mr_plt$exposure=="GlycA"&mr_plt$outcome=="FI","egger_intercept"]),
                  SNP = factor(SNP, levels = c("IVW", "MR Egger", "Weighted median", "Weighted mode", "MR-PRESSO"))),
              aes(intercept=a, slope=b, colour=SNP, linetype=SNP), show.legend=T)  +
  scale_shape_manual(values=c(20,8), guide="none") +  
  scale_colour_manual(labels = c("IVW\n(β = 0.37, p = 0.004)",
                                 "MR Egger\n(β = -0.02, p = 0.91)",
                                 "Weighted median\n(β = 0.34, p = 0.014)",
                                 "Weighted mode\n(β = 0.13, p = 0.51)",
                                 "MR-PRESSO\n(β = 0.47, p = 2.2E-5)"),
                      values=c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")) + 
  scale_linetype_manual(labels = c("IVW\n(β = 0.37, p = 0.004)",
                                   "MR Egger\n(β = -0.02, p = 0.91)",
                                   "Weighted median\n(β = 0.34, p = 0.014)",
                                   "Weighted mode\n(β = 0.13, p = 0.51)",
                                   "MR-PRESSO\n(β = 0.47, p = 2.2E-5)"),
                        values=c(1,6,2,4,5,3)) +
  labs(x="SNP-GlycA association", y="SNP-FI association") +
  theme_bw() +
  theme(legend.position=c(.8,.23), 
        legend.direction="vertical",
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=.1),
        legend.spacing.y = unit(2, "mm")) +
  guides(colour = guide_legend(byrow = TRUE),
         fill = guide_legend(byrow = TRUE))

singleplot_GlycA_FI <- ggplot(mr_GlycA_fi_single %>% mutate(method = factor(method, levels = c("","Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode", "MR-PRESSO"))),
                              aes(y = SNP, x = b)) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_errorbarh(aes(xmin = lo, xmax = up, size = method, colour = method), height = 0) + 
  geom_point(aes(colour = method), shape=20) + 
  geom_hline(aes(yintercept = which(levels(SNP) %in% "")), color = "grey") + 
  scale_size_manual(values = c(.3,.8,.8,.8,.8,.8)) +
  scale_colour_manual(values = c("black","#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")) + 
  labs(y = "", x = "MR effect size of GlycA on FI",
       caption = "* 6 SNPs were identified as outliers by MR-PRESSO")  +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 5, colour = "black"), 
        axis.ticks.y = element_line(size = 0), 
        axis.title.x = element_text(size = 9),
        plot.caption = element_text(size = 7))

cairo_pdf("Output/MR_results/MR_GlycA_FI.pdf", width = 8, height = 5, pointsize = 12)
scatterplot_GlycA_FI + singleplot_GlycA_FI + 
  plot_layout(widths = c(5,2)) + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face="bold"))
dev.off()


##### Scatter plot and forest plot for MR estimates of Creatinine-FI association ##### 

# Outlier SNPs identified by MR-PRESSO
Creatinine_FI_outlier <- mr_data_metabolomics[row.names(mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure=="Creatinine" & attributes(mrpresso_metabolomics)$outcome=="FI")]]$`MR-PRESSO results`$`Outlier Test`[mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure=="Creatinine" & attributes(mrpresso_metabolomics)$outcome=="FI")]]$`MR-PRESSO results`$`Outlier Test`$Pvalue<.05,]), "SNP"]

mr_Creatinine_fi_single <- rbind(
  # Association between each Creatinine-associated SNP and FI
  mr_metabolomics_results_single[mr_metabolomics_results_single$exposure=="Creatinine" & mr_metabolomics_results_single$outcome=="FI" & !mr_metabolomics_results_single$SNP %in% c("All - Inverse variance weighted","All - MR Egger"),] %>%
    mutate(SNP = case_when(
      SNP %in% Creatinine_FI_outlier ~ paste0(SNP, "*"),
      !SNP %in% Creatinine_FI_outlier ~ SNP),
      method="") %>%
    select(exposure, outcome, SNP, method, b, se, p),
  # Main MR results for Creatinine-FI association
  mr_summary_metabolomics[mr_summary_metabolomics$exposure=="Creatinine"&mr_summary_metabolomics$outcome=="FI",]  %>%
    rename(p=pval) %>% 
    mutate(SNP=case_when(
      method=="Inverse variance weighted" ~ "IVW",
      method=="MR Egger" ~ "MR Egger",
      method=="Weighted median" ~ "Weighted median",
      method=="Weighted mode" ~ "Weighted mode",
      method=="MR-PRESSO" ~ "MR-PRESSO"
    )) %>%
    select(exposure, outcome, SNP, method, b, se, p)
)
mr_Creatinine_fi_single$up <- mr_Creatinine_fi_single$b + 1.96 * mr_Creatinine_fi_single$se # confidence interval
mr_Creatinine_fi_single$lo <- mr_Creatinine_fi_single$b - 1.96 * mr_Creatinine_fi_single$se
nom <- mr_Creatinine_fi_single$SNP[grepl("rs", mr_Creatinine_fi_single$SNP)]
nom <- nom[order(mr_Creatinine_fi_single$b)] # sort effect size
mr_Creatinine_fi_single <- rbind(mr_Creatinine_fi_single, data.frame(exposure="Creatinine",outcome="FI",SNP="",method="",b=NA,se=NA,p=NA,up=NA,lo=NA))
mr_Creatinine_fi_single$SNP <- factor(mr_Creatinine_fi_single$SNP, levels = c("MR-PRESSO", "Weighted mode", "Weighted median", "MR Egger", "IVW", "", nom))

scatterplot_Creatinine_FI <- ggplot(data=mr_data_metabolomics[mr_data_metabolomics$exposure=="Creatinine" & mr_data_metabolomics$outcome=="FI" & mr_data_metabolomics$mr_keep,] %>%
                                      mutate(beta.exposure_flipped=case_when(beta.exposure < 0 ~ beta.exposure * -1,
                                                                             beta.exposure >= 0 ~ beta.exposure),
                                             beta.outcome_flipped=case_when(beta.exposure < 0 ~ beta.outcome * -1,
                                                                            beta.exposure >= 0 ~ beta.outcome),
                                             outlier=case_when(SNP%in%Creatinine_FI_outlier ~ 1,
                                                               !SNP%in%Creatinine_FI_outlier ~ 0)),
                                    aes(x=beta.exposure_flipped, y=beta.outcome_flipped)) +
  geom_errorbar(aes(ymin=beta.outcome_flipped-se.outcome, ymax=beta.outcome_flipped+se.outcome), colour="grey", width=0) +
  geom_errorbarh(aes(xmin=beta.exposure_flipped-se.exposure, xmax=beta.exposure_flipped+se.exposure), colour="grey", height=0) +
  geom_point(aes(shape=factor(outlier))) +
  geom_abline(data=mr_Creatinine_fi_single[!grepl("rs", mr_Creatinine_fi_single$SNP) & mr_Creatinine_fi_single$SNP!="",] %>% 
                mutate(a = case_when(
                  SNP!="MR Egger" ~ 0,
                  SNP=="MR Egger" ~ mr_metabolomics_plt[mr_metabolomics_plt$exposure=="Creatinine"&mr_metabolomics_plt$outcome=="FI","egger_intercept"]),
                  SNP = factor(SNP, levels = c("IVW", "MR Egger", "Weighted median", "Weighted mode", "MR-PRESSO"))),
              aes(intercept=a, slope=b, colour=SNP, linetype=SNP), show.legend=T)  +
  scale_shape_manual(values=c(20,8), guide="none") +  
  scale_colour_manual(labels = c("IVW\n(β = 0.38, p = 0.008)",
                                 "MR Egger\n(β = 0.23, p = 0.62)",
                                 "Weighted median\n(β = 0.58, p = 2.7E-5)",
                                 "Weighted mode\n(β = 0.70, p = 0.006)",
                                 "MR-PRESSO\n(β = 0.51, p = 8.3E-5)"),
                      values=c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")) + 
  scale_linetype_manual(labels = c("IVW\n(β = 0.38, p = 0.008)",
                                   "MR Egger\n(β = 0.23, p = 0.62)",
                                   "Weighted median\n(β = 0.58, p = 2.7E-5)",
                                   "Weighted mode\n(β = 0.70, p = 0.006)",
                                   "MR-PRESSO\n(β = 0.51, p = 8.3E-5)"),
                        values=c(1,6,2,4,5,3)) +
  labs(x="SNP-creatinine association", y="SNP-FI association") +
  theme_bw() +
  theme(legend.position=c(.8,.23), 
        legend.direction="vertical",
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=.1),
        legend.spacing.y = unit(2, "mm")) +
  guides(colour = guide_legend(byrow = TRUE),
         fill = guide_legend(byrow = TRUE))

singleplot_Creatinine_FI <- ggplot(mr_Creatinine_fi_single %>% mutate(method = factor(method, levels = c("","Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode", "MR-PRESSO"))),
                                   aes(y = SNP, x = b)) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_errorbarh(aes(xmin = lo, xmax = up, size = method, colour = method), height = 0) + 
  geom_point(aes(colour = method), shape=20) + 
  geom_hline(aes(yintercept = which(levels(SNP) %in% "")), color = "grey") + 
  scale_size_manual(values = c(.3,.8,.8,.8,.8,.8)) +
  scale_colour_manual(values = c("black","#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")) + 
  labs(y = "", x = "MR effect size of creatinine on FI",
       caption = "* 3 SNPs were identified as outliers by MR-PRESSO")  +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 4, colour = "black"), 
        axis.ticks.y = element_line(size = 0), 
        axis.title.x = element_text(size = 9),
        plot.caption = element_text(size = 7))

cairo_pdf("Output/MR_results/MR_Creatinine_FI.pdf", width = 8, height = 5, pointsize = 12)
scatterplot_Creatinine_FI + singleplot_Creatinine_FI + 
  plot_layout(widths = c(5,2)) + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face="bold"))
dev.off()


##### Combine plots for GlycA-FI and Creatinine-FI ##### 

cairo_pdf("Output/MR_results/MR_GlycA_Creatinine_FI.pdf", width = 8, height = 10, pointsize = 12)
scatterplot_GlycA_FI + singleplot_GlycA_FI + scatterplot_Creatinine_FI + singleplot_Creatinine_FI + 
  plot_layout(widths = c(5,2)) + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face="bold"))
dev.off()

cairo_pdf("Output/MR_results/MR_GlycA_Creatinine_FI_scatter.pdf", width = 10.5, height = 5, pointsize = 12)
scatterplot_GlycA_FI + scatterplot_Creatinine_FI +  
  plot_layout(ncol = 2) + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face="bold"))
dev.off()



##### Forest plots of MR results for modified FIs ##### 

mr_summary_mod_fi_sensitivity <- read.delim("Output/MR_results/mr_summary_mod_fi_sensitivity.txt", sep=" ")
cairo_pdf("Output/MR_results/Sensitivity/MR_sensitivity_modified_FIs.pdf", width = 8, height = 6, pointsize = 10)

plot_grid(
  ggforestplot::forestplot(
  df = mr_summary_mod_fi_sensitivity %>% 
    filter(method=="Inverse variance weighted" & !exposure%in%c("Total cholesterol","LDL-C","Triglycerides","ApoB")) %>% 
    select(outcome,exposure,b,se,pval) %>%
    mutate(Outcome=factor(case_when(
      outcome=="FI, nocancer" ~ "FI, excluded cancer items",
      outcome=="FI, nocardio" ~ "FI, excluded cardiometabolic items",
      outcome=="FI, nocranial" ~ "FI, excluded cranial items",
      outcome=="FI, nogastro" ~ "FI, excluded gastrointestinal items",
      outcome=="FI, noimmune" ~ "FI, excluded immunological items",
      outcome=="FI, noinfirmity" ~ "FI, excluded infirmity items",
      outcome=="FI, nomental" ~ "FI, excluded mental wellbeing items",
      outcome=="FI, nomusculo" ~ "FI, excluded musculoskeletal items",
      outcome=="FI, nopain" ~ "FI, excluded pain items",
      outcome=="FI, norespiratory" ~ "FI, excluded respiratory items",
      outcome=="FI, nosensory" ~ "FI, excluded sensory items"),
      levels=c("FI, excluded sensory items","FI, excluded respiratory items",
               "FI, excluded pain items","FI, excluded musculoskeletal items",
               "FI, excluded mental wellbeing items","FI, excluded infirmity items",
               "FI, excluded immunological items","FI, excluded gastrointestinal items",
               "FI, excluded cranial items","FI, excluded cardiometabolic items",
               "FI, excluded cancer items")),
      Exposure=gsub("_", "-",exposure),
      Category="NMR metabolomic biomarkers"),
  name = Exposure,
  estimate = b,
  se = se,
  pvalue = pval,
  psignif = .011,
  xlab = "IVW-MR estimate",
  colour = Outcome,
  shape = Outcome,
  xlim = c(-.45,.8),
  xtickbreaks = c(-.4,-.2,0,.2,.4,.6,.8),) + 
  scale_color_brewer(palette = "Spectral") +
  scale_shape_manual(values = c(21,22,23,24,21,22,23,24,21,22,23)) +
  ggforce::facet_col(facets = ~ Category, scales = "free_y", space = "free" ) +
  theme(legend.position="none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)),

ggforestplot::forestplot(
  df = mr_summary_mod_fi_sensitivity %>% 
    filter(method=="Inverse variance weighted" & exposure%in%c("Total cholesterol","LDL-C","Triglycerides","ApoB")) %>% 
    select(outcome,exposure,b,se,pval) %>%
    mutate(Outcome=factor(case_when(
      outcome=="FI, nocancer" ~ "FI, excluded cancer items",
      outcome=="FI, nocardio" ~ "FI, excluded cardiometabolic items",
      outcome=="FI, nocranial" ~ "FI, excluded cranial items",
      outcome=="FI, nogastro" ~ "FI, excluded gastrointestinal items",
      outcome=="FI, noimmune" ~ "FI, excluded immunological items",
      outcome=="FI, noinfirmity" ~ "FI, excluded infirmity items",
      outcome=="FI, nomental" ~ "FI, excluded mental wellbeing items",
      outcome=="FI, nomusculo" ~ "FI, excluded musculoskeletal items",
      outcome=="FI, nopain" ~ "FI, excluded pain items",
      outcome=="FI, norespiratory" ~ "FI, excluded respiratory items",
      outcome=="FI, nosensory" ~ "FI, excluded sensory items"),
      levels=c("FI, excluded sensory items","FI, excluded respiratory items",
               "FI, excluded pain items","FI, excluded musculoskeletal items",
               "FI, excluded mental wellbeing items","FI, excluded infirmity items",
               "FI, excluded immunological items","FI, excluded gastrointestinal items",
               "FI, excluded cranial items","FI, excluded cardiometabolic items",
               "FI, excluded cancer items")),
      Exposure=gsub("_", "-",exposure),
      Category="Clinical biomarkers"),
  name = Exposure,
  estimate = b,
  se = se,
  pvalue = pval,
  psignif = .011,
  xlab = "IVW-MR estimate",
  colour = Outcome,
  shape = Outcome,
  xlim = c(-.2,.8),
  xtickbreaks = c(-.2,0,.2,.4,.6,.8),) + 
  scale_color_brewer(palette = "Spectral") +
  scale_shape_manual(values = c(21,22,23,24,21,22,23,24,21,22,23)) +
  ggforce::facet_col(facets = ~ Category, scales = "free_y", space = "free" ) +
  theme(legend.position="bottom", 
        legend.direction="vertical",
        legend.text=element_text(size=9),
        legend.title=element_text(size=9.5, face="bold"),
        legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=.1),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)),

rel_widths = c(1,.95)
)

dev.off()

# =============================== END OF FILE  ===============================