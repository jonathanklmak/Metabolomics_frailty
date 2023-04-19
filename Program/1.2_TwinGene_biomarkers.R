#==============================================================================
# FILENAME: 1.2_TwinGene_biomarkers.R
# PROJECT: 	Metabolomics_frailty
# PURPOSE:  To assess associations between metabolic biomarkers and FI in TwinGene
# AUTHOR:   Jonathan Mak & Laura Kananen
# CREATED:	2022-02-18
# UPDATED: 	2023-02-21
# R VERSION: 4.1.3
#==============================================================================

### Required packages
library(dplyr)
library(ggplot2); library(patchwork); library(ggforce) # Plotting graphs
library(haven)  # Import Stata file
library(readxl) # Import Excel file
library(psych)  # For descriptive statistics
library(reshape2) # To "melt" dataset
library(ggcorrplot) # Compute a correlation matrix
library(ggforestplot) # For forest plots of metabolites
library(broom)  # For nice outputs
library(glmnet) # For penalized regression
library(drgee)  # For generalized estimating equations
library(tableone) # For descriptive statistics

#========================= IMPORT TWINGENE DATA ==============================

##### Import frailty and covariate variables from TwinGene #####

twingene_frailty <- read_dta("Data/Raw_TwinGene_data/SALT_FI.dta")
colnames(twingene_frailty); dim(twingene_frailty) # Total N = 43,800 (those missing sex and birth date were already excluded)
twingene_frailty <- twingene_frailty[,c("twinnr","pairid","bestzyg","sex","interview_age_value","birthdate","first_finish_date","bmi","eversmok","edyrs","fi_44_tv")] # Keep only relevant covariate and FI variables
colnames(twingene_frailty)[5] <- "age"
colnames(twingene_frailty)[11] <- "fi"
twingene_frailty$fi_percent <- twingene_frailty$fi*100 # FI multiplied by 100, for later regression analysis
twingene_frailty$sex <- ifelse(twingene_frailty$sex==2, 0, 1) # Recode sex; 0=women, 1=men
summary(twingene_frailty$fi) # Total-varying FI, maximum 44 items; those with >20% missing frailty items were excluded (i.e., 43,641 had FI available)

# Add additional covariate data (smoking & alcohol)
salt_covariates <- read_dta("Data/Raw_TwinGene_data/SALT_covariates.dta")
salt_covariates <- haven::as_factor(salt_covariates, levels="label") %>% as.data.frame() # Change variables to factors using labels in Stata file
  
# Combine frailty and covariate data
twingene_frailty <- merge(twingene_frailty, salt_covariates, by="twinnr", all.x=T) %>%
  # Create education level variable
  mutate(education=factor(case_when(
    edyrs<9 ~ 1,
    edyrs>=9 & edyrs <=12 ~ 2,
    edyrs>12 ~ 3), 
    labels = c("Low","Intermediate","High")))
table(twingene_frailty$smoking); table(twingene_frailty$eversmok)
twingene_frailty$smoking <- ifelse(is.na(twingene_frailty$smoking)&twingene_frailty$eversmok==0,1,twingene_frailty$smoking) # Recode 8 missing values of smoking to "Never"
twingene_frailty$smoking <- factor(twingene_frailty$smoking, labels = c("Never","Previous","Current"))
colnames(twingene_frailty)
rm(salt_covariates)


##### Import 250 metabolites variables from 3 separated datasets #####

twingene_metabolites_1 <- read.delim("Data/Raw_TwinGene_data/Results_20743a_t.txt", header=T, check.names=F)
twingene_metabolites_2 <- read.delim("Data/Raw_TwinGene_data/Results_20743B_t.txt", header=T, check.names=F)
twingene_metabolites_3 <- read.delim("Data/Raw_TwinGene_data/Results_20743C_t.txt", header=T, check.names=F)
twingene_metabolites_quality_1 <- read.delim("Data/Raw_TwinGene_data/quality_20743a_t.txt", header=T, check.names=F)
twingene_metabolites_quality_2 <- read.delim("Data/Raw_TwinGene_data/quality_20743B_t.txt", header=T, check.names=F)
twingene_metabolites_quality_3 <- read.delim("Data/Raw_TwinGene_data/quality_20743C_t.txt", header=T, check.names=F)

### Remove top rows of each dataset (these rows in TwinGene contain irrelevant information)
twingene_metabolites_1 <- twingene_metabolites_1[-(1:4),-c(2:4,255)]
twingene_metabolites_2 <- twingene_metabolites_2[-(1:4),-c(2:4,255)]
twingene_metabolites_3 <- twingene_metabolites_3[-(1:4),-c(2:4,255)]
twingene_metabolites_quality_1 <- twingene_metabolites_quality_1[-(1:3),c("TWINNR","High lactate","High pyruvate","Low glucose","Low protein")]
twingene_metabolites_quality_2 <- twingene_metabolites_quality_2[-(1:3),c("TWINNR","High lactate","High pyruvate","Low glucose","Low protein")]
twingene_metabolites_quality_3 <- twingene_metabolites_quality_3[-(1:3),c("TWINNR","High lactate","High pyruvate","Low glucose","Low protein")]

### Combine the 3 datasets
twingene_metabolites <- rbind(twingene_metabolites_1, twingene_metabolites_2, twingene_metabolites_3)
twingene_metabolites_quality <- rbind(twingene_metabolites_quality_1, twingene_metabolites_quality_2, twingene_metabolites_quality_3)
twingene_metabolites <- merge(twingene_metabolites, twingene_metabolites_quality, by="TWINNR", all.x=T)
colnames(twingene_metabolites)[1] <- "twinnr" # rename ID variable
twingene_metabolites[, 2:255] <- apply(twingene_metabolites[, 2:255], 2, function(x) as.numeric(x)) # change metabolites to numeric
twingene_metabolites <- twingene_metabolites %>%
  mutate(QCflagged = ifelse(`High lactate`==1 | `High pyruvate`==1 | `Low glucose`==1 | `Low protein`==1, 1, NA))
rm(list=c("twingene_metabolites_1", "twingene_metabolites_2", "twingene_metabolites_3","twingene_metabolites_quality_1","twingene_metabolites_quality_2","twingene_metabolites_quality_3"))


##### Import blood biomarkers #####

twingene_blood_biomarkers_long <- read_dta("Data/Raw_TwinGene_data/V_LABDATA.dta")
colnames(twingene_blood_biomarkers_long); dim(twingene_blood_biomarkers_long) # Data are in long format

### Reshape dataset from long to wide
twingene_blood_biomarkers <- dcast(as.data.frame(twingene_blood_biomarkers_long[,c("twinnr","analysisdesc","value")]),
                                   twinnr~analysisdesc ) 

### Remove missing values
twingene_blood_biomarkers[, 2:14] <- apply(twingene_blood_biomarkers[, 2:14], 2, function(x) as.numeric(x)) # change biomarkers to numeric
na.omit(twingene_blood_biomarkers[,c("S/P-CRP högkänsligt","S/P-CRP högkänsligt2")]) # n=10 have both variables available, but with same value
twingene_blood_biomarkers[["S/P-CRP högkänsligt"]] <- rowMeans(twingene_blood_biomarkers[,c("S/P-CRP högkänsligt","S/P-CRP högkänsligt2")], na.rm=T) # Combine the two variables
twingene_blood_biomarkers <- twingene_blood_biomarkers[,!colnames(twingene_blood_biomarkers) %in% "S/P-CRP h?gk?nsligt2"]
twingene_blood_biomarkers$`B-HbA1c` <- ifelse( twingene_blood_biomarkers$`B-HbA1c`<0,NA,twingene_blood_biomarkers$`B-HbA1c` ) # Remove negative values
twingene_blood_biomarkers$`B-Hemoglobin` <- ifelse( twingene_blood_biomarkers$`B-Hemoglobin`<0,NA,twingene_blood_biomarkers$`B-Hemoglobin` ) # Remove negative values

describe(twingene_blood_biomarkers) # High missingness for Anti PC and PAF AH -> remove
rm(twingene_blood_biomarkers_long)


##### Label biomarkers #####

### Import labels of biomarkers from a prepared Excel file
twingene_labels <- read_excel("Documents/Biomarker_list.xlsx", sheet = "TwinGene") 
twingene_labels <- twingene_labels[twingene_labels$Biomarker!="Anti PC" & twingene_labels$Biomarker!="PAF AH",] # Remove oestradiol and rheumatoid factor (high missingness)
twingene_labels <- twingene_labels[twingene_labels$Biomarker!="Glycerol" & twingene_labels$Biomarker!="Hemoglobin",] # Remove glycerol and hemoglobin (not available in UKB)
twingene_metabolomics_list <- twingene_labels[["Variable name"]][1:168]; twingene_metabolomics_list # Name of the 168 measured metabolomics
twingene_blood_biomarkers_list <- twingene_labels[["Variable name"]][250:258]; twingene_blood_biomarkers_list # Name of the 9 blood biomarkers

### Reorder biomarkers
twingene_metabolites <- twingene_metabolites[,c("twinnr",twingene_labels[["Other name"]][1:168],"QCflagged")] # Re-order TwinGene metabolomics variables according to order in the Excel file
twingene_blood_biomarkers <- twingene_blood_biomarkers[,c("twinnr",twingene_labels[["Other name"]][250:258])] # Re-order TwinGene blood biomarker variables according to order in the Excel file

### Rename biomarkers
colnames(twingene_metabolites) <- c("twinnr",twingene_metabolomics_list,"QCflagged"); colnames(twingene_metabolites) # Rename 168 measured metabolic biomarkers
colnames(twingene_blood_biomarkers) <- c("twinnr",twingene_blood_biomarkers_list); colnames(twingene_blood_biomarkers)  # Rename 10 blood biomarkers

### Re-scale metabolites to be consistent with the units in UK Biobank
twingene_metabolites$Creatinine <- twingene_metabolites$Creatinine/1000 # Convert creatinine from mircomole/l to mmol/l
twingene_blood_biomarkers$HBA1C_RBC <- 10.929*(twingene_blood_biomarkers$HBA1C_RBC-2.15) # Convert HBA1c from % to mmol/mol


#========================== SAMPLE SELECTION =================================

##### 0. Combine TwinGene frailty & biomarker datasets #####

twingene <- merge(twingene_metabolites, twingene_blood_biomarkers, by="twinnr", all=T) # Combine metabolomics and blood biomarker data
twingene <- merge(twingene_frailty, twingene, by="twinnr", all.y=T) # Combine frailty and metabolites data
colnames(twingene)
N_original <- dim(twingene)[1]; N_original # Total number of TwinGene participants = 12,648
rm(twingene_frailty); rm(twingene_metabolites); rm(twingene_blood_biomarkers)

##### 1. Exclude those with missing frailty index #####

twingene <- twingene[ complete.cases( twingene[, c("fi")] ), ] # Select only rows with complete data on frailty measures (N = 12,579)
N_nonmissing <- dim(twingene)[1]; N_nonmissing # After excluding those with missing frailty data, N = 12,579

##### 2. Select only those with complete biomarker data #####

dim(na.omit(twingene[,twingene_metabolomics_list])) # n = 11,842 have all 168 measured metabolomics non-missing
dim(na.omit(twingene[,twingene_blood_biomarkers_list])) # n = 11,703 have all 9 blood biomarkers non-missing
dim(na.omit(twingene[,c(twingene_metabolomics_list,twingene_blood_biomarkers_list)])) # n = 11,059 have all 168 & 9 biomarkers non-missing 

# Exclude individuals with non-complete biomarker data
twingene_complete <- twingene[ complete.cases( twingene[, c(twingene_metabolomics_list,twingene_blood_biomarkers_list)] ), ] 

##### 3. Remove metabolites with QC flags #####

QCflagged_id <- twingene_complete[which(twingene_complete$QCflagged==1), "twinnr"]  # N = 34 with QC flags
twingene_complete <- twingene_complete[!(twingene_complete$twinnr %in% QCflagged_id),] # Exclude those with any QC flags from dataset, N = 11,025

##### 4. Create dataset with standardized metabolites #####
twingene_z_complete <- cbind(twingene_complete[,c("twinnr","pairid","bestzyg","age","sex","bmi","eversmok","edyrs","smoking","alcohol","education","fi_percent")],
                             scale(twingene_complete[,c(twingene_metabolomics_list,twingene_blood_biomarkers_list)], center=T, scale=T)) # Standardization of metabolites
#rm(twingene_complete)


#======================= DESCRIPTIVE STATISTICS ==============================

##### Descriptive statistics of the 168+9 biomarkers #####
twingene_biomarkers_descriptives <- describe(twingene_complete[,c(twingene_metabolomics_list,twingene_blood_biomarkers_list)], na.rm=T, interp=F, skew=F, quant=c(.25,.5,.75))
clipr::write_clip(twingene_biomarkers_descriptives)

##### Sample characteristics #####
table_var <- c("age","sex","bmi","smoking","alcohol","education","fi_percent")
table_catvar <- c("sex","smoking","alcohol","education")
print(CreateTableOne(vars = table_var, data = twingene_complete, factorVars = table_catvar),
      formatOptions = list(big.mark = ","), quote = TRUE, noSpaces = TRUE)
rm(list=c("table_catvar","table_var"))


#====== UNIVARIATE ANALYSIS: GEE BETWEEN EACH BIOMARKER AND FRAILTY ===========

##### Linear regression between frailty and each biomarker, adjusted for age and sex #####

# Transform dataset from wide to long, for running linear regression in loops
twingene_z_fi_age_sex_biomarkers <- melt(twingene_z_complete[,c("twinnr","pairid","age","sex","fi_percent",twingene_metabolomics_list, twingene_blood_biomarkers_list)], id.vars = c("twinnr","pairid", "fi_percent", "age", "sex")) 

# Generalized estimating equation with cluster robust standard error (used in TwinGene to account for relatedness within twin pairs)
twingene_gee_univariate_fi_age_sex_biomarkers <- twingene_z_fi_age_sex_biomarkers %>% 
  group_by(variable) %>% 
  do(as_tibble(summary(gee(fi_percent ~ value + age + sex, ., clusterid="pairid", link="identity", cond=F))$coef, rownames = "term")) %>% 
  filter(term == "value") %>% 
  mutate(beta.twingene = Estimate,    # Beta-coefficients from linear regression model
         se.twingene = `Std. Error`,  # Standard errors from linear regression model
         pvalue.twingene = `Pr(>|z|)` # P-values from linear regression model
  ) %>% 
  select(beta.twingene, se.twingene, pvalue.twingene) %>% 
  as.data.frame()
twingene_gee_univariate_fi_age_sex_biomarkers$p.05 <- twingene_gee_univariate_fi_age_sex_biomarkers$pvalue<.05 # An indicator showing whether the metabolite has p<.05
# Calculate R2 of linear regression models
twingene_gee_univariate_fi_age_sex_biomarkers <- twingene_z_fi_age_sex_biomarkers %>% group_by(variable) %>% 
  do(data.frame(r2.twingene = summary(lm(fi_percent ~ value + age + sex, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(twingene_gee_univariate_fi_age_sex_biomarkers,. , by=c("variable"="variable"))
twingene_gee_univariate_fi_age_sex_biomarkers <- twingene_gee_univariate_fi_age_sex_biomarkers %>% mutate(rank_r2.twingene = rank(-r2.twingene)) # An indicator showing the rank of r2, from largest to smallest

rm(list=c("twingene_z_fi_age_sex_biomarkers"))

clipr::write_clip(twingene_gee_univariate_fi_age_sex_biomarkers)
write.table(twingene_gee_univariate_fi_age_sex_biomarkers, file = "Output/Observational_results/twingene_gee_univariate_fi_age_sex_biomarkers.txt")


##### Linear regression between frailty and each biomarker, adjusted for age, sex, BMI, smoking, alcohol, years of education #####

# Transform dataset from wide to long, for running linear regression in loops
twingene_z_fi_allcov_biomarkers <- melt(twingene_z_complete[,c("twinnr","pairid","age","sex","bmi","smoking","alcohol","edyrs","fi_percent",twingene_metabolomics_list, twingene_blood_biomarkers_list)], id.vars = c("twinnr","pairid", "fi_percent", "age", "sex","bmi","smoking","alcohol","edyrs")) 

# Generalized estimating equation with cluster robust standard error (used in TwinGene to account for relatedness within twin pairs)
twingene_gee_univariate_fi_allcov_biomarkers <- twingene_z_fi_allcov_biomarkers %>% 
  group_by(variable) %>% 
  do(as_tibble(summary(gee(fi_percent ~ value + age + sex + bmi + smoking + alcohol + edyrs, ., clusterid="pairid", link="identity", cond=F))$coef, rownames = "term")) %>% 
  filter(term == "value") %>% 
  mutate(beta.allcov.twingene = Estimate,    # Beta-coefficients from linear regression model
         se.allcov.twingene = `Std. Error`,  # Standard errors from linear regression model
         pvalue.allcov.twingene = `Pr(>|z|)` # P-values from linear regression model
  ) %>% 
  select(beta.allcov.twingene, se.allcov.twingene, pvalue.allcov.twingene) %>% 
  as.data.frame()
twingene_gee_univariate_fi_allcov_biomarkers$p.05 <- twingene_gee_univariate_fi_allcov_biomarkers$pvalue<.05 # An indicator showing whether the metabolite has p<.05
# Calculate R2 of linear regression models
twingene_gee_univariate_fi_allcov_biomarkers <- twingene_z_fi_allcov_biomarkers %>% group_by(variable) %>% 
  do(data.frame(r2.twingene = summary(lm(fi_percent ~ value + age + sex + bmi + smoking + alcohol + edyrs, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(twingene_gee_univariate_fi_allcov_biomarkers,. , by=c("variable"="variable"))
twingene_gee_univariate_fi_allcov_biomarkers <- twingene_gee_univariate_fi_allcov_biomarkers %>% mutate(rank_r2.twingene = rank(-r2.twingene)) # An indicator showing the rank of r2, from largest to smallest

rm(list=c("twingene_z_fi_allcov_biomarkers"))

clipr::write_clip(twingene_gee_univariate_fi_allcov_biomarkers)
write.table(twingene_gee_univariate_fi_allcov_biomarkers, file = "Output/Observational_results/twingene_gee_univariate_fi_allcov_biomarkers.txt")


##### Additional analysis: co-twin control analysis #####
table(twingene_z_complete$bestzyg) # 1=MZ, 2=DZ same sex, 3=unknown, 4=DZ opposite sex
### Dataset for DZ complete twin pairs
twingene_z_complete_DZ <- twingene_z_complete[ave(twingene_z_complete$pairid,twingene_z_complete$pairid,FUN=length)==2 & (twingene_z_complete$bestzyg==2 | twingene_z_complete$bestzyg==4),]
dim(twingene_z_complete_DZ) # n=5524 DZ twins
### Dataset for MZ complete twin pairs
twingene_z_complete_MZ <- twingene_z_complete[ave(twingene_z_complete$pairid,twingene_z_complete$pairid,FUN=length)==2 & twingene_z_complete$bestzyg==1,]
dim(twingene_z_complete_MZ) # n=2264 MZ twins

### Co-twin control analysis for GlycA-FI association
cotwin_results_GlycA_fi <- rbind(
  # Population estimate in full cohort
  summary(gee(fi_percent ~ GlycA + age + sex + bmi + smoking + alcohol + edyrs, twingene_z_complete, clusterid="pairid", link="identity", cond=F))$coefficient[2,],
  # Population estimate in complete DZ twins
  summary(gee(fi_percent ~ GlycA + age + sex + bmi + smoking + alcohol + edyrs, twingene_z_complete_DZ, clusterid="pairid", link="identity", cond=F))$coefficient[2,],
  # Population estimate in complete MZ twins
  summary(gee(fi_percent ~ GlycA + age + sex + bmi + smoking + alcohol + edyrs, twingene_z_complete_MZ, clusterid="pairid", link="identity", cond=F))$coefficient[2,],
  # Within-twin-pair estimate in complete DZ twins
  summary(gee(fi_percent ~ GlycA + age + sex + bmi + smoking + alcohol + edyrs, twingene_z_complete_DZ, clusterid="pairid", link="identity", cond=T))$coefficient[1,],
  # Within-twin-pair estimate in complete MZ twins
  summary(gee(fi_percent ~ GlycA + age + bmi + smoking + alcohol + edyrs, twingene_z_complete_MZ, clusterid="pairid", link="identity", cond=T))$coefficient[1,]
) %>% as.data.frame()
cotwin_results_GlycA_fi$model <- c("Population-level estimate (full sample)","Population-level estimate (DZ twins)","Population-level estimate (MZ twins)","Within-twin-pair estimate (DZ twins)", "Within-twin-pair estimate (MZ twins)")
ggforestplot::forestplot( df = cotwin_results_GlycA_fi, name = model, se = `Std. Error`, estimate = `Estimate`, pvalue = `Pr(>|z|)`,
  psignif = .05, xlab = "Difference in FI (%) per SD increase in GlycA", logodds = F,
  xtickbreaks = c(-.2,0,.2,.4,.6,.8))  +
  theme(axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA))

### Co-twin control analysis for Creatinine-FI association
cotwin_results_Creatinine_fi <- rbind(
  # Population estimate in full cohort
  summary(gee(fi_percent ~ Creatinine + age + sex + bmi + smoking + alcohol + edyrs, twingene_z_complete, clusterid="pairid", link="identity", cond=F))$coefficient[2,],
  # Population estimate in complete DZ twins
  summary(gee(fi_percent ~ Creatinine + age + sex + bmi + smoking + alcohol + edyrs, twingene_z_complete_DZ, clusterid="pairid", link="identity", cond=F))$coefficient[2,],
  # Population estimate in complete MZ twins
  summary(gee(fi_percent ~ Creatinine + age + sex + bmi + smoking + alcohol + edyrs, twingene_z_complete_MZ, clusterid="pairid", link="identity", cond=F))$coefficient[2,],
  # Within-twin-pair estimate in complete DZ twins
  summary(gee(fi_percent ~ Creatinine + age + sex + bmi + smoking + alcohol + edyrs, twingene_z_complete_DZ, clusterid="pairid", link="identity", cond=T))$coefficient[1,],
  # Within-twin-pair estimate in complete MZ twins
  summary(gee(fi_percent ~ Creatinine + age + bmi + smoking + alcohol + edyrs, twingene_z_complete_MZ, clusterid="pairid", link="identity", cond=T))$coefficient[1,]
) %>% as.data.frame()
cotwin_results_Creatinine_fi$model <- c("Population-level estimate (full sample)","Population-level estimate (DZ twins)","Population-level estimate (MZ twins)","Within-twin-pair estimate (DZ twins)", "Within-twin-pair estimate (MZ twins)")
ggforestplot::forestplot( df = cotwin_results_Creatinine_fi, name = model, se = `Std. Error`, estimate = `Estimate`, pvalue = `Pr(>|z|)`,
                          psignif = .05, xlab = "Difference in FI (%) per SD increase in creatinine", logodds = F,
                          xtickbreaks = c(-.2,0,.2,.4,.6,.8))  +
  theme(axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA))


#===================== SAVE DATA FOR FURTHER ANALYSIS ========================

save.image("Data/R_data/TwinGene_biomarkers.Rdata")
#write.table(twingene_z, file = "Data/Cleaned_data/TwinGene_metabolomics_sd.txt", sep = "\t", row.names = F, na="")

# =============================== END OF FILE  ===============================