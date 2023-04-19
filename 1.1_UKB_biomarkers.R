#==============================================================================
# FILENAME: 1.1_UKB_biomarkers.R
# PROJECT: 	Metabolomics_frailty
# PURPOSE:  To identify frailty-associated metabolic biomarkers in UK Biobank
# AUTHOR:   Jonathan Mak & Laura Kananen
# CREATED:	2021-10-14
# UPDATED: 	2023-02-21
# R VERSION: 4.1.3
#==============================================================================
### Notes:
# 1. Details of the 249 metabolites (168 measured & 81 derived) can be found on UKB website (Category 220):
#    https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=220
# 2. The 168 metabolites were available for N=~118,000, measured from plasma samples at baseline recruitment
# 3. Formulas for calculating the 81 metabolite ratios were available from:
#    https://biobank.ndph.ox.ac.uk/ukb/ukb/docs/Nightingale_ratio_calculation.tsv
# 4. The blood (Category 17518) and urine biomarkers (Category 100083) can be found on UKB website:
#    https://biobank.ctsu.ox.ac.uk/crystal/label.cgi?id=17518
#    https://biobank.ctsu.ox.ac.uk/crystal/label.cgi?id=100083
# 5. Oestradiol & rheumatoid factor were excluded due to high missingness
# 6. Frailty was primarily defined using the frailty index (0-100%) and frailty phenotype (ordinal, ranges from 0 to 5)
#    Ref: https://doi.org/10.1093/gerona/gly094; https://doi.org/10.1016/S2468-2667(18)30091-4
# 7. Only White population were included in the main analysis

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
library(tableone) # For descriptive statistics


#====================== IMPORT UK BIOBANK DATA ===============================

##### Import frailty and covariate variables #####

ukb_frailty <- read_dta("Data/Raw_UKB_data/Frailty_covariate/ukb_frailty_covariate.dta")
ukb_frailty <- haven::as_factor(ukb_frailty, levels="label") # Change variables to factors using labels in Stata file
colnames(ukb_frailty); dim(ukb_frailty) # Total N = 502,631

##### Import biomarker variables #####

# Import 168 measured metabolites
ukb_metabolites <- read.delim("Data/Raw_UKB_data/Biomarker/ukb_metabolomics.txt", header=T)
colnames(ukb_metabolites)[1] <- "eid" # rename ID variable
ukb_metabolites <- ukb_metabolites[rowSums(is.na(ukb_metabolites[,2:169])) 
                                   != ncol(ukb_metabolites[,2:169]), ] # Remove rows with all metabolites missing
colnames(ukb_metabolites); dim(ukb_metabolites) # 168 metabolites in 118,021 individuals

# Import 32 blood and urine biomarkers (remove Oestradiol & Rheumatoid factor due to high missingness)
ukb_blood_biomarkers <- read.delim("Data/Raw_UKB_data/Biomarker/ukb_blood_biomarkers.txt", header=T)
ukb_urine_biomarkers <- read.delim("Data/Raw_UKB_data/Biomarker/ukb_urine_biomarkers.txt", header=T)
ukb_blood_urine_biomarkers <- merge(ukb_blood_biomarkers, ukb_urine_biomarkers, by="f.eid", all=T)
colnames(ukb_blood_urine_biomarkers)[1] <- "eid" # rename ID variable
describe(ukb_blood_urine_biomarkers) # f.30810.0.0 (Oestradiol; n=76669) & f.30820.0.0 (Rheumatoid factor; n=41314) has high missingness -> remove
ukb_blood_urine_biomarkers <- ukb_blood_urine_biomarkers[, !colnames(ukb_blood_urine_biomarkers) 
                                                         %in% c("f.30800.0.0",  # Remove Oestradiol
                                                                "f.30820.0.0")] # Remove Rheumatoid factor
ukb_blood_urine_biomarkers <- ukb_blood_urine_biomarkers[rowSums(is.na(ukb_blood_urine_biomarkers[,2:33]))
                                                         != ncol(ukb_blood_urine_biomarkers[,2:33]), ] # Remove rows with all biomarkers missing
colnames(ukb_blood_urine_biomarkers); dim(ukb_blood_urine_biomarkers) # 32 biomarkers in 498,541 individuals
rm(ukb_blood_biomarkers);rm(ukb_urine_biomarkers)

##### Label biomarkers #####

# Import labels of biomarkers from a prepared Excel file (information was based on description in UKB website)
ukb_labels <- read_excel("Documents/Biomarker_list.xlsx", sheet = "UKB")
ukb_labels <- ukb_labels[ukb_labels$Biomarker!="Oestradiol" & ukb_labels$Biomarker!="Rheumatoid factor",] # Remove oestradiol and rheumatoid factor (high missingness)
ukb_metabolomics_list <- ukb_labels[["Variable name"]][1:168]; ukb_metabolomics_list # Name of the 249 metabolomics
ukb_blood_urine_biomarkers_list <- ukb_labels[["Variable name"]][250:281]; ukb_blood_urine_biomarkers_list # Name of the 32 serum/urine biomarkers

# Derive additional 81 metabolite ratios, based on the existing 168 metabolites (*updated 20220503: ratios are NOT used in final analysis)
ukb_metabolites <- ukb_metabolites[,c("eid",ukb_labels[["ID1"]][1:168])] # Re-order UKB metabolomics variables according to order in the Excel file
#colnames(ukb_metabolites)
#for(i in 169:249) {                                      # A loop to create the 81 metabolites (starting from metabolite #169 to #249)
#  ukb_metabolites[[ukb_labels[["Variable name"]][i]]] <- # Name of the new ratio-based metabolite
#    (ukb_metabolites[[ukb_labels[["ID1"]][i]]]  /        # Numerator metabolite
#       ukb_metabolites[[ukb_labels[["ID2"]][i]]])  *     # Denominator metabolite
#    ukb_labels[["Percent"]][i]                           # Multiply by 100 if it is a percentage
#}
#rm(i)
colnames(ukb_metabolites) <- c("eid",ukb_metabolomics_list) # Rename the first 168 metabolites
colnames(ukb_metabolites)

# Rename 32 biomarker variables
ukb_blood_urine_biomarkers <- ukb_blood_urine_biomarkers[,c("eid",ukb_labels[["ID1"]][250:281])] # Re-order UKB 32 biomarker variables according to order in the Excel file
colnames(ukb_blood_urine_biomarkers) <- c("eid",ukb_blood_urine_biomarkers_list) # Rename the 32 biomarker variables
colnames(ukb_blood_urine_biomarkers) 

##### Import biomarker QC data #####

ukb_metabolomics_QC <- read.delim("Data/Raw_UKB_data/Biomarker/ukb_metabolomics_QC.txt", header=T)
colnames(ukb_metabolomics_QC) <- c("eid","Shipment.Plate","Spectrometer","High.Lactate","High.Pyruvate","Low.Glucose","Low.Protein","Sample.Measured.Date.and.Time","Sample.Prepared.Date.and.Time","Well.position.within.plate")
ukb_metabolomics_QC <- ukb_metabolomics_QC[rowSums(is.na(ukb_metabolomics_QC[,2:9]))
                                           != ncol(ukb_metabolomics_QC[,2:9]), ] # Remove rows with all data missing
# An indicator showing whether the individual has any QC flag
ukb_metabolomics_QC <- ukb_metabolomics_QC %>%
  mutate(QCflagged = ifelse(High.Lactate==1 | High.Pyruvate==1 | Low.Glucose==1 | Low.Protein== 1, 1, NA))


#=========================== SAMPLE SELECTION ================================

##### 0. Combine UKB frailty & biomarker datasets #####

ukb <- merge(ukb_metabolites, ukb_blood_urine_biomarkers, by="eid", all=T) # Combine metabolomics and serum/urine biomarker data
ukb <- merge(ukb_frailty, ukb, by="eid", all=T) # Combine frailty and metabolites data
ukb <- merge(ukb, ukb_metabolomics_QC, by="eid", all.x=T)
colnames(ukb)
N_original <- dim(ukb)[1]; N_original # Total number of UKB participants = 502,631
rm(ukb_frailty); rm(ukb_metabolites); rm(ukb_blood_urine_biomarkers)

##### 1. Exclude individuals withdrawn from UKB study #####

withdrawn_id <- read.csv("Documents/withdrawal_list_20210809.csv", header=F) # List of participants withdrawn from UKB
ukb <- ukb[!(ukb$eid %in% withdrawn_id[,1]),] # Exclude 172 participants withdrawn from UKB
N_participants <- dim(ukb)[1]; N_participants # After excluding those withdrawn from UKB, N = 502,459
rm(withdrawn_id)

##### 2. Exclude those with missing frailty measures #####

ukb <- ukb[ complete.cases( ukb[, c("fi_percent")] ), ] # Select only rows with complete data on frailty index (N = 500,166)
N_nonmissing <- dim(ukb)[1]; N_nonmissing # After excluding those missing frailty data, N = 500,166

##### 3. Exclude non-White individuals #####

ukb_nonwhite <- subset(ukb, ukb$ethnicity!="White") # Non-white participants, N = 26,567
ukb <- subset(ukb, ukb$ethnicity=="White") # Exclude non-White participants in main analysis
N_white <- dim(ukb)[1]; N_white # After excluding non-White participants, N = 471,906

##### 4. Select only those with complete biomarker data #####

dim(na.omit(ukb[,ukb_metabolomics_list])) # n = 104,378 have all 168 metabolomics non-missing
dim(na.omit(ukb[,ukb_blood_urine_biomarkers_list])) # n = 67,488 have all 32 blood/urine biomarkers non-missing
dim(na.omit(ukb[,c(ukb_metabolomics_list,ukb_blood_urine_biomarkers_list)])) # Only n = 15,277 have all 168 + 32 biomarkers non-missing -> Separate into two sets for LASSO, one for 168 metabolomics, another for 32 blood/urine biomarkers

# Subgroup 1: complete data on 168 metabolomics
ukb_168metabolomics_complete <- ukb[ complete.cases( ukb[, c(ukb_metabolomics_list)] ), ] 
N_168metabolomics <- dim(ukb_168metabolomics_complete)[1]; N_168metabolomics # Individuals with complete data on 168 metabolomics, N = 104,378

# Subgroup 2: complete data on 32 biomarkers
ukb_32biomarkers_complete <- ukb[ complete.cases( ukb[, c(ukb_blood_urine_biomarkers_list)] ), ] 
N_32biomarkers <- dim(ukb_32biomarkers_complete)[1]; N_32biomarkers # Individuals with complete data on 32 biomarkers, N = 67,488

##### 5. Remove metabolites with QC flags from subgroup 1 #####

QCflagged_168metabolomics_id <- ukb_168metabolomics_complete[which(ukb_168metabolomics_complete$QCflagged==1), "eid"]  # N = 13,805 have at least 1 metabolite QC flagged
ukb_168metabolomics_complete <- ukb_168metabolomics_complete[!(ukb_168metabolomics_complete$eid %in% QCflagged_168metabolomics_id),] # Exclude those with any QC flags from dataset, N = 90,573

# Scale biomarkers to SD units
ukb_z_168metabolomics_complete <- cbind(ukb_168metabolomics_complete[,c("eid","age","sex","fi_percent","fp_total","assessment_centre","smoking","alcohol_3group","education","di","bmi")],
                                        scale(ukb_168metabolomics_complete[,c(ukb_metabolomics_list)], center=T, scale=T)) # Standardization of metabolites
ukb_z_32biomarkers_complete <- cbind(ukb_32biomarkers_complete[,c("eid","age","sex","fi_percent","fp_total","assessment_centre","smoking","alcohol_3group","education","di","bmi")],
                                     scale(ukb_32biomarkers_complete[,c(ukb_blood_urine_biomarkers_list)], center=T, scale=T)) # Standardization of metabolites


##### 6. Additionally exclude outliers of biomarkers (for sensitivity analysis) #####

outlier_5iqr_168metabolomics_id <- ukb_168metabolomics_complete[which(rowSums( plyr::colwise(function(x){x > quantile(x)[3] + 5*IQR(x) | x < quantile(x)[3] - 5*IQR(x)})(ukb_168metabolomics_complete[,ukb_metabolomics_list])) >= 1), "eid"]  # N = 4,048 have at least 1 metabolite outside 5 IQR
outlier_10iqr_168metabolomics_id <- ukb_168metabolomics_complete[which(rowSums( plyr::colwise(function(x){x > quantile(x)[3] + 10*IQR(x) | x < quantile(x)[3] - 10*IQR(x)})(ukb_168metabolomics_complete[,ukb_metabolomics_list])) >= 1), "eid"]  # N = 709 have at least 1 metabolite outside 10 IQR

outlier_5iqr_32biomarkers_id <- ukb_32biomarkers_complete[which(rowSums( plyr::colwise(function(x){x > quantile(x)[3] + 5*IQR(x) | x < quantile(x)[3] - 5*IQR(x)})(ukb_32biomarkers_complete[,ukb_blood_urine_biomarkers_list])) >= 1), "eid"]  # N = 10,188 have at least 1 metabolite outside 5 IQR
outlier_10iqr_32biomarkers_id <- ukb_32biomarkers_complete[which(rowSums( plyr::colwise(function(x){x > quantile(x)[3] + 10*IQR(x) | x < quantile(x)[3] - 10*IQR(x)})(ukb_32biomarkers_complete[,ukb_blood_urine_biomarkers_list])) >= 1), "eid"]  # N = 3,809 have at least 1 metabolite outside 5 IQR

### Remove metabolites outside 5 IQR from median
ukb_z_168metabolomics_complete_nooutlier5iqr <- ukb_z_168metabolomics_complete[!(ukb_z_168metabolomics_complete$eid %in% outlier_5iqr_168metabolomics_id),] # Exclude outliers (outside 5 IQR) from dataset, N = 89525
ukb_z_32biomarkers_complete_nooutlier5iqr <- ukb_z_32biomarkers_complete[!(ukb_z_32biomarkers_complete$eid %in% outlier_5iqr_32biomarkers_id),] # Exclude outliers (outside 5 IQR) from dataset, N = 57300

### Remove metabolites outside 10 IQR from median
ukb_z_168metabolomics_complete_nooutlier10iqr <- ukb_z_168metabolomics_complete[!(ukb_z_168metabolomics_complete$eid %in% outlier_10iqr_168metabolomics_id),] # Exclude outliers (outside 10 IQR) from dataset, N = 89864
ukb_z_32biomarkers_complete_nooutlier10iqr <- ukb_z_32biomarkers_complete[!(ukb_z_32biomarkers_complete$eid %in% outlier_10iqr_32biomarkers_id),] # Exclude outliers (outside 10 IQR) from dataset, N = 63679


#========================= DESCRIPTIVE STATISTICS ============================

##### Sample characteristics #####

table_var <- c("age","sex","assessment_centre","smoking","alcohol_3group","education","di","bmi","fi_percent","fi_cat","fp_total","fp_cat")
table_catvar <- c("sex","assessment_centre","smoking","alcohol_3group","education","fi_cat","fp_cat")
print(CreateTableOne(vars = table_var, data = ukb, factorVars = table_catvar),
      formatOptions = list(big.mark = ","), quote = TRUE, noSpaces = TRUE)
print(CreateTableOne(vars = table_var, data = ukb_168metabolomics_complete, factorVars = table_catvar),
      formatOptions = list(big.mark = ","), quote = TRUE, noSpaces = TRUE)
print(CreateTableOne(vars = table_var, data = ukb_32biomarkers_complete, factorVars = table_catvar),
      formatOptions = list(big.mark = ","), quote = TRUE, noSpaces = TRUE)
rm(list=c("table_catvar","table_var"))

##### Descriptive statistics of the 168+32 biomarkers #####
ukb_biomarkers_descriptives <- rbind(describe(ukb_168metabolomics_complete[,ukb_metabolomics_list],
                                              na.rm=T, interp=F, skew=F, quant=c(.25,.5,.75)),
                                     describe(ukb_32biomarkers_complete[,ukb_blood_urine_biomarkers_list],
                                              na.rm=T, interp=F, skew=F, quant=c(.25,.5,.75)))
clipr::write_clip(ukb_biomarkers_descriptives)

# Histograms of NMR metabolites
pdf("Output/Observational_results/Histogram_UKB_168metabolomics_complete_1.pdf", width = 24, height = 14, pointsize = 8)
histogram_metabolites_1 <- 
  ggplot(melt(ukb_168metabolomics_complete[,ukb_metabolomics_list[1:84]] %>%
                `colnames<-`(ukb_labels[["Biomarker"]][c(1:84)])),
         aes(x = value)) +   
  facet_wrap(~ variable, scales = "free", ncol=12) + 
  geom_histogram(bins=60)  +
  labs(y="Frequency", x="Value") +
  theme(axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12))
print(histogram_metabolites_1)
dev.off()
rm(histogram_metabolites_1)

pdf("Output/Observational_results/Histogram_UKB_168metabolomics_complete_2.pdf", width = 24, height = 14, pointsize = 8)
histogram_metabolites_2 <- 
  ggplot(melt(ukb_168metabolomics_complete[,ukb_metabolomics_list[85:168]] %>%
                `colnames<-`(ukb_labels[["Biomarker"]][c(85:168)])),
         aes(x = value)) +   
  facet_wrap(~ variable, scales = "free", ncol=12) + 
  geom_histogram(bins=60)  +
  labs(y="Frequency", x="Value") +
  theme(axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12))
print(histogram_metabolites_2)
dev.off()
rm(histogram_metabolites_2)

# Histograms of clinical biomarkers
pdf("Output/Observational_results/Histogram_UKB_32biomarkers_complete.pdf", width = 18, height = 9, pointsize = 8)
histogram_clinical_biomarker <- 
  ggplot(melt(ukb_32biomarkers_complete[,ukb_blood_urine_biomarkers_list] %>%
                `colnames<-`(ukb_labels[["Biomarker"]][250:281])),
         aes(x = value)) +   
  facet_wrap(~ variable, scales = "free", ncol=8) + 
  geom_histogram(bins=60)  +
  labs(y="Frequency", x="Value") + 
  theme(axis.text = element_text(size = 8),
        strip.text.x = element_text(size = 11.5))
print(histogram_clinical_biomarker)
dev.off()
rm(histogram_clinical_biomarker)


#============== UNIVARIATE ANALYSIS: CORRELATION AND HEATMAPS ================

##### Spearman correations of metabolites #####

# Correlation matrix for metabolites, age, FI
corr_matrix_168metabolomics <- round( cor( ukb_z_168metabolomics_complete[,c(ukb_metabolomics_list)], # Correlation matrix for 168 metabolites
                                           method="spearman" ), 9) 
corr_matrix_32biomarkers <- round( cor( ukb_z_32biomarkers_complete[,c(ukb_blood_urine_biomarkers_list)], # Correlation matrix for 32 biomarkers
                                        method="spearman" ), 9) 
colnames(corr_matrix_168metabolomics) <- c(ukb_labels[["Biomarker"]][1:168])
colnames(corr_matrix_32biomarkers) <- c(ukb_labels[["Biomarker"]][250:281])
row.names(corr_matrix_168metabolomics) <- c(ukb_labels[["Biomarker"]][1:168])
row.names(corr_matrix_32biomarkers) <- c(ukb_labels[["Biomarker"]][250:281])

# Correlation p-values matrix
corr_pmatrix_168metabolomics <- cor_pmat( ukb_z_168metabolomics_complete[,c(ukb_metabolomics_list)] )
corr_pmatrix_32biomarkers <- cor_pmat( ukb_z_32biomarkers_complete[,c(ukb_blood_urine_biomarkers_list)] )

# Heatmap of metabolites
pdf("Output/Observational_results/Heatmap_UKB_168metabolomics.pdf", width = 30, height = 30, pointsize = 10)
ggcorrplot(corr_matrix_168metabolomics,
           hc.order = T, # Hierarchical clustering
           outline.col = "white",  
           p.mat = corr_pmatrix_168metabolomics ,
           insig = "blank", # Leave blank on no significant coefficient
           pch.cex = 10, tl.cex = 14, 
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46780"),
           tl.srt = 90) +
  theme(legend.key.height = unit(2, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=24))
dev.off()

pdf("Output/Observational_results/Heatmap_UKB_32biomarkers.pdf", width = 12, height = 12, pointsize = 11)
ggcorrplot(corr_matrix_32biomarkers,
           hc.order = T, # Hierarchical clustering
           outline.col = "white",  
           p.mat = corr_pmatrix_32biomarkers ,
           insig = "blank", # Leave blank on no significant coefficient
           pch.cex = 10, tl.cex = 14, 
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46780"),
           tl.srt = 45) +
  theme(legend.key.height = unit(1.5, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

rm(list=c("corr_matrix_168metabolomics","corr_matrix_32biomarkers","corr_pmatrix_168metabolomics","corr_pmatrix_32biomarkers"))

##### Scatterplots of FI and biomarkers #####

jpeg(file="Output/Observational_results/Scatterplot_UKB_FI_168metabolomics_complete.jpg", width = 600, height = 600, units="mm", res = 300, quality=100)
scatterplot_metabolites <- 
  ggplot(melt(ukb_z_168metabolomics_complete[,c("fi_percent",ukb_metabolomics_list)] %>%
                `colnames<-`(c("fi_percent",ukb_labels[["Biomarker"]][c(1:168)])), id.vars="fi_percent"),
         aes(x = value, y = fi_percent)) +   
  facet_wrap(~ variable, scales = "free_x", ncol=14) + 
  geom_point(alpha=0.1) +
  geom_smooth(method = "loess", se = F, color="darkgreen") +
  labs(y="FI (%)", x="Biomarker z-score") + 
  theme(axis.text = element_text(size = 6),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12))
print(scatterplot_metabolites)
dev.off()
rm(scatterplot_metabolites)  

jpeg(file="Output/Observational_results/Scatterplot_UKB_FI_32biomarkers_complete.jpg", width = 180, height = 250, units="mm", res = 300, quality=100)
scatterplot_metabolites <- 
  ggplot(melt(ukb_z_32biomarkers_complete[,c("fi_percent",ukb_blood_urine_biomarkers_list)] %>%
                `colnames<-`(c("fi_percent",ukb_labels[["Biomarker"]][250:281])), id.vars="fi_percent"),
         aes(x = value, y = fi_percent)) +   
  facet_wrap(~ variable, scales = "free_x", ncol=4) + 
  geom_point(alpha=0.1) +
  geom_smooth(method = "loess", se = F, color="darkred") +
  labs(y="FI (%)", x="Biomarker z-score") + 
  theme(axis.text = element_text(size = 6))
print(scatterplot_metabolites)
dev.off()
rm(scatterplot_metabolites)


#============ MAIN ANALYSIS: LINEAR REGRESSION AND LASSO FOR FI ===============

##### 1. Univariate linear regression models between FI and each biomarker, adjusted for age and sex #####

ukb_z_fi_age_sex_168metabolomics <- melt(ukb_z_168metabolomics_complete[,c("eid","age","sex","fi_percent",ukb_metabolomics_list)], id.vars = c("eid", "fi_percent", "age", "sex")) 
ukb_z_fi_age_sex_32biomarkers <- melt(ukb_z_32biomarkers_complete[,c("eid","age","sex","fi_percent",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fi_percent", "age", "sex")) 

### 168 metabolomics
ukb_lm_univariate_fi_age_sex_168metabolomics <- ukb_z_fi_age_sex_168metabolomics %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.age.sex = estimate, # Beta-coefficients from linear regression model
         se.fi.age.sex = std.error,  # Standard errors from linear regression model
         pvalue.fi.age.sex = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.age.sex, se.fi.age.sex, pvalue.fi.age.sex) %>% 
  as.data.frame()
ukb_lm_univariate_fi_age_sex_168metabolomics$bon.p.fi.age.sex <- ukb_lm_univariate_fi_age_sex_168metabolomics$pvalue.fi.age.sex<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fi_age_sex_168metabolomics <- ukb_z_fi_age_sex_168metabolomics %>% group_by(variable) %>% 
  do(data.frame(r2.fi.age.sex = summary(lm(fi_percent ~ value + age + sex, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fi_age_sex_168metabolomics,. , by=c("variable"="variable"))
ukb_lm_univariate_fi_age_sex_168metabolomics <- ukb_lm_univariate_fi_age_sex_168metabolomics %>% mutate(rank_r2.fi.age.sex = rank(-r2.fi.age.sex)) # An indicator showing the rank of r2, from largest to smallest

### 32 biomarkers
ukb_lm_univariate_fi_age_sex_32biomarkers <- ukb_z_fi_age_sex_32biomarkers %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.age.sex = estimate, # Beta-coefficients from linear regression model
         se.fi.age.sex = std.error,  # Standard errors from linear regression model
         pvalue.fi.age.sex = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.age.sex, se.fi.age.sex, pvalue.fi.age.sex) %>% 
  as.data.frame()
ukb_lm_univariate_fi_age_sex_32biomarkers$bon.p.fi.age.sex <- ukb_lm_univariate_fi_age_sex_32biomarkers$pvalue.fi.age.sex<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fi_age_sex_32biomarkers <- ukb_z_fi_age_sex_32biomarkers %>% group_by(variable) %>% 
  do(data.frame(r2.fi.age.sex = summary(lm(fi_percent ~ value + age + sex, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fi_age_sex_32biomarkers,. , by=c("variable"="variable"))
ukb_lm_univariate_fi_age_sex_32biomarkers <- ukb_lm_univariate_fi_age_sex_32biomarkers %>% mutate(rank_r2.fi.age.sex = rank(-r2.fi.age.sex)) # An indicator showing the rank of r2, from largest to smallest

### Combine results
rm(list=c("ukb_z_fi_age_sex_168metabolomics","ukb_z_fi_age_sex_32biomarkers"))
clipr::write_clip(ukb_lm_univariate_fi_age_sex_168metabolomics)
clipr::write_clip(ukb_lm_univariate_fi_age_sex_32biomarkers)
write.table(ukb_lm_univariate_fi_age_sex_168metabolomics, file = "Output/Observational_results/ukb_lm_univariate_fi_age_sex_168metabolomics.txt")
write.table(ukb_lm_univariate_fi_age_sex_32biomarkers, file = "Output/Observational_results/ukb_lm_univariate_fi_age_sex_32biomarkers.txt")


##### 2. Univariate linear regression models between FI and each biomarker, adjusted for age, sex, assessment center, smoking, education, deprivation index, BMI #####

ukb_z_fi_allcov_168metabolomics <- melt(ukb_z_168metabolomics_complete[,c("eid","age","sex","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_metabolomics_list)], id.vars = c("eid", "fi_percent", "age", "sex","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 
ukb_z_fi_allcov_32biomarkers <- melt(ukb_z_32biomarkers_complete[,c("eid","age","sex","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fi_percent", "age", "sex","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 

### 168 metabolomics
ukb_lm_univariate_fi_allcov_168metabolomics <- ukb_z_fi_allcov_168metabolomics %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov, se.fi.allcov, pvalue.fi.allcov) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_168metabolomics$bon.p.fi.allcov <- ukb_lm_univariate_fi_allcov_168metabolomics$pvalue.fi.allcov<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fi_allcov_168metabolomics <- ukb_z_fi_allcov_168metabolomics %>% group_by(variable) %>% 
  do(data.frame(r2.fi.allcov = summary(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fi_allcov_168metabolomics,. , by=c("variable"="variable"))
ukb_lm_univariate_fi_allcov_168metabolomics <- ukb_lm_univariate_fi_allcov_168metabolomics %>% mutate(rank_r2.fi.allcov = rank(-r2.fi.allcov)) # An indicator showing the rank of r2, from largest to smallest

### 32 biomarkers
ukb_lm_univariate_fi_allcov_32biomarkers <- ukb_z_fi_allcov_32biomarkers %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov, se.fi.allcov, pvalue.fi.allcov) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_32biomarkers$bon.p.fi.allcov <- ukb_lm_univariate_fi_allcov_32biomarkers$pvalue.fi.allcov<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fi_allcov_32biomarkers <- ukb_z_fi_allcov_32biomarkers %>% group_by(variable) %>% 
  do(data.frame(r2.fi.allcov = summary(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fi_allcov_32biomarkers,. , by=c("variable"="variable"))
ukb_lm_univariate_fi_allcov_32biomarkers <- ukb_lm_univariate_fi_allcov_32biomarkers %>% mutate(rank_r2.fi.allcov = rank(-r2.fi.allcov)) # An indicator showing the rank of r2, from largest to smallest

### Summary of results
rm(list=c("ukb_z_fi_allcov_168metabolomics","ukb_z_fi_allcov_32biomarkers"))
clipr::write_clip(ukb_lm_univariate_fi_allcov_168metabolomics)
clipr::write_clip(ukb_lm_univariate_fi_allcov_32biomarkers)
write.table(ukb_lm_univariate_fi_allcov_168metabolomics, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics.txt")
write.table(ukb_lm_univariate_fi_allcov_32biomarkers, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers.txt")


##### 3. LASSO for FI and biomarkers, age, sex #####
 
### Select variables to be included in LASSO

# 168 metabolomics
lasso_168metabolomics <- ukb_z_168metabolomics_complete[, c("fi_percent", ukb_metabolomics_list, "age","sex") ]
Y_168metabolomics_fi <- lasso_168metabolomics$fi_percent; length(Y_168metabolomics_fi) # N = 90,573, the outcome in FI model
X_168metabolomics_fi <- model.matrix( fi_percent ~ ., lasso_168metabolomics[,c("fi_percent",ukb_metabolomics_list,"age","sex")])[,-1]; dim(X_168metabolomics_fi) # N = 90,573, 170 variables (168 metabolites, age, sex)

# 32 blood/urine biomarkers
lasso_32biomarkers <- ukb_z_32biomarkers_complete[, c("fi_percent", ukb_blood_urine_biomarkers_list, "age","sex") ]
Y_32biomarker_fi <- lasso_32biomarkers$fi_percent; length(Y_32biomarker_fi) # N = 67,488, the outcome in FI model
X_32biomarker_fi <- model.matrix( fi_percent ~ ., lasso_32biomarkers[,c("fi_percent",ukb_blood_urine_biomarkers_list,"age","sex")])[,-1]; dim(X_32biomarker_fi) # N = 67,488, 34 variables (32 biomarkers, age, sex)

### LASSO regression for FI, 10-fold cross-validation

# 168 metabolomics
set.seed(202204)
ukb_cvlasso_168metabolomics_fi <- cv.glmnet( x = X_168metabolomics_fi,
                                             y = Y_168metabolomics_fi,
                                             alpha=1,
                                             standardize = T,
                                             nfolds = 10 ) # Find the best lambda using 10-fold cross-validation
ukb_cvlasso_168metabolomics_fi # Results of cross-validation
# Call:  cv.glmnet(x = X_168metabolomics_fi, y = Y_168metabolomics_fi, nfolds = 10, alpha = 1, standardize = T)
#
#  Measure: Mean-Squared Error 
#
#       Lambda Index Measure     SE Nonzero
# min 0.000151   100   46.44 0.2616     165
# 1se 0.005677    61   46.69 0.2636      58

# 32 blood/urine biomarkers
set.seed(202204)
ukb_cvlasso_32biomarker_fi <- cv.glmnet( x = X_32biomarker_fi,
                                         y = Y_32biomarker_fi,
                                         alpha=1,
                                         standardize = T,
                                         nfolds = 10 ) # Find the best lambda using 10-fold cross-validation
ukb_cvlasso_32biomarker_fi # Results of cross-validation
# Call:  cv.glmnet(x = X_32biomarker_fi, y = Y_32biomarker_fi, nfolds = 10, alpha = 1, standardize = T) 
#
#  Measure: Mean-Squared Error 
#
#      Lambda Index Measure     SE Nonzero
# min 0.00165    77   48.93 0.1726      34
# 1se 0.04698    41   49.09 0.1756      23

# Plots of coefficients and MSE
pdf("Output/Observational_results/LASSO_UKB_10foldcv_FI_biomarkers_age_sex.pdf", width = 10, height = 5, pointsize = 12)
par(mfrow=c(1,2))
plot(ukb_cvlasso_168metabolomics_fi)
title("(a)  NMR metabolomic biomarkers", adj = 0, line = 3)
plot(ukb_cvlasso_32biomarker_fi)
title("(b)  Clinical biomarkers", adj = 0, line = 3)
dev.off()

### Combine LASSO results

# Beta-coefficients in LASSO regression for FI, using lambda within 1SE of min MSE
ukb_coef_lasso_fi_249metabolomics <- as.data.frame(as.matrix(coef(ukb_cvlasso_168metabolomics_fi, ukb_cvlasso_168metabolomics_fi$lambda.1se)))[-c(1,170,171),,drop=F]
ukb_coef_lasso_fi_32biomarker <- as.data.frame(as.matrix(coef(ukb_cvlasso_32biomarker_fi, ukb_cvlasso_32biomarker_fi$lambda.1se)))[-c(1,34,35),,drop=F]

ukb_coef_lasso_fi_249metabolomics$variable <- row.names(ukb_coef_lasso_fi_249metabolomics)
ukb_coef_lasso_fi_32biomarker$variable <- row.names(ukb_coef_lasso_fi_32biomarker)

colnames(ukb_coef_lasso_fi_249metabolomics)[1] <- "beta.LASSO.fi"
colnames(ukb_coef_lasso_fi_32biomarker)[1] <- "beta.LASSO.fi"

row.names(ukb_coef_lasso_fi_249metabolomics) <- NULL
row.names(ukb_coef_lasso_fi_32biomarker) <- NULL

ukb_coef_lasso_fi_249metabolomics <- ukb_coef_lasso_fi_249metabolomics[,c("variable","beta.LASSO.fi")]
ukb_coef_lasso_fi_32biomarker <- ukb_coef_lasso_fi_32biomarker[,c("variable","beta.LASSO.fi")]

rm(list=c("lasso_168metabolomics","lasso_32biomarkers",
          "Y_168metabolomics_fi","X_168metabolomics_fi","Y_32biomarker_fi","X_32biomarker_fi",
          "ukb_cvlasso_168metabolomics_fi","ukb_cvlasso_32biomarker_fi"))

clipr::write_clip(ukb_coef_lasso_fi_249metabolomics)
clipr::write_clip(ukb_coef_lasso_fi_32biomarker)
write.table(ukb_coef_lasso_fi_249metabolomics, file = "Output/Observational_results/ukb_coef_lasso_fi_249metabolomics.txt")
write.table(ukb_coef_lasso_fi_32biomarker, file = "Output/Observational_results/ukb_coef_lasso_fi_32biomarker.txt")



#=============== SENSITIVITY ANALYSIS 1: EXCLUDING OUTLIERS ===================

##### Univariate linear regression models between FI and each biomarker, adjusted for age and sex #####

ukb_z_fi_age_sex_168metabolomics_nooutlier5iqr <- melt(ukb_z_168metabolomics_complete_nooutlier5iqr[,c("eid","age","sex","fi_percent",ukb_metabolomics_list)], id.vars = c("eid", "fi_percent", "age", "sex")) 
ukb_z_fi_age_sex_32biomarkers_nooutlier5iqr <- melt(ukb_z_32biomarkers_complete_nooutlier5iqr[,c("eid","age","sex","fi_percent",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fi_percent", "age", "sex")) 

### 168 metabolomics
ukb_lm_univariate_fi_age_sex_168metabolomics_nooutlier5iqr <- ukb_z_fi_age_sex_168metabolomics_nooutlier5iqr %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.age.sex.nooutlier5iqr = estimate, # Beta-coefficients from linear regression model
         se.fi.age.sex.nooutlier5iqr = std.error,  # Standard errors from linear regression model
         pvalue.fi.age.sex.nooutlier5iqr = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.age.sex.nooutlier5iqr, se.fi.age.sex.nooutlier5iqr, pvalue.fi.age.sex.nooutlier5iqr) %>% 
  as.data.frame()
ukb_lm_univariate_fi_age_sex_168metabolomics_nooutlier5iqr$bon.p.fi.age.sex.nooutlier5iqr <- ukb_lm_univariate_fi_age_sex_168metabolomics_nooutlier5iqr$pvalue.fi.age.sex.nooutlier5iqr<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fi_age_sex_168metabolomics_nooutlier5iqr <- ukb_z_fi_age_sex_168metabolomics_nooutlier5iqr %>% group_by(variable) %>% 
  do(data.frame(r2.fi.age.sex.nooutlier5iqr = summary(lm(fi_percent ~ value + age + sex, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fi_age_sex_168metabolomics_nooutlier5iqr,. , by=c("variable"="variable"))
ukb_lm_univariate_fi_age_sex_168metabolomics_nooutlier5iqr <- ukb_lm_univariate_fi_age_sex_168metabolomics_nooutlier5iqr %>% mutate(rank_r2.fi.age.sex.nooutlier5iqr = rank(-r2.fi.age.sex.nooutlier5iqr)) # An indicator showing the rank of r2, from largest to smallest

### 32 biomarkers
ukb_lm_univariate_fi_age_sex_32biomarkers_nooutlier5iqr <- ukb_z_fi_age_sex_32biomarkers_nooutlier5iqr %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.age.sex.nooutlier5iqr = estimate, # Beta-coefficients from linear regression model
         se.fi.age.sex.nooutlier5iqr = std.error,  # Standard errors from linear regression model
         pvalue.fi.age.sex.nooutlier5iqr = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.age.sex.nooutlier5iqr, se.fi.age.sex.nooutlier5iqr, pvalue.fi.age.sex.nooutlier5iqr) %>% 
  as.data.frame()
ukb_lm_univariate_fi_age_sex_32biomarkers_nooutlier5iqr$bon.p.fi.age.sex.nooutlier5iqr <- ukb_lm_univariate_fi_age_sex_32biomarkers_nooutlier5iqr$pvalue.fi.age.sex.nooutlier5iqr<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fi_age_sex_32biomarkers_nooutlier5iqr <- ukb_z_fi_age_sex_32biomarkers_nooutlier5iqr %>% group_by(variable) %>% 
  do(data.frame(r2.fi.age.sex.nooutlier5iqr = summary(lm(fi_percent ~ value + age + sex, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fi_age_sex_32biomarkers_nooutlier5iqr,. , by=c("variable"="variable"))
ukb_lm_univariate_fi_age_sex_32biomarkers_nooutlier5iqr <- ukb_lm_univariate_fi_age_sex_32biomarkers_nooutlier5iqr %>% mutate(rank_r2.fi.age.sex.nooutlier5iqr = rank(-r2.fi.age.sex.nooutlier5iqr)) # An indicator showing the rank of r2, from largest to smallest

### Summary of results
rm(list=c("ukb_z_fi_age_sex_168metabolomics_nooutlier5iqr","ukb_z_fi_age_sex_32biomarkers_nooutlier5iqr"))
clipr::write_clip(ukb_lm_univariate_fi_age_sex_168metabolomics_nooutlier5iqr)
clipr::write_clip(ukb_lm_univariate_fi_age_sex_32biomarkers_nooutlier5iqr)
write.table(ukb_lm_univariate_fi_age_sex_168metabolomics_nooutlier5iqr, file = "Output/Observational_results/ukb_lm_univariate_fi_age_sex_168metabolomics_nooutlier5iqr.txt")
write.table(ukb_lm_univariate_fi_age_sex_32biomarkers_nooutlier5iqr, file = "Output/Observational_results/ukb_lm_univariate_fi_age_sex_32biomarkers_nooutlier5iqr.txt")


##### Univariate linear regression models between FI and each biomarker, adjusted for age, sex, assessment center, smoking, education, deprivation index, BMI #####

ukb_z_fi_allcov_168metabolomics_nooutlier5iqr <- melt(ukb_z_168metabolomics_complete_nooutlier5iqr[,c("eid","age","sex","assessment_centre","smoking","alcohol_3group","education","di","bmi","fi_percent",ukb_metabolomics_list)], id.vars = c("eid", "fi_percent", "age", "sex","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 
ukb_z_fi_allcov_32biomarkers_nooutlier5iqr <- melt(ukb_z_32biomarkers_complete_nooutlier5iqr[,c("eid","age","sex","assessment_centre","smoking","alcohol_3group","education","di","bmi","fi_percent",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fi_percent", "age", "sex","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 

### 168 metabolomics
ukb_lm_univariate_fi_allcov_168metabolomics_nooutlier5iqr <- ukb_z_fi_allcov_168metabolomics_nooutlier5iqr %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.nooutlier5iqr = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.nooutlier5iqr = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.nooutlier5iqr = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.nooutlier5iqr, se.fi.allcov.nooutlier5iqr, pvalue.fi.allcov.nooutlier5iqr) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_168metabolomics_nooutlier5iqr$bon.p.fi.allcov.nooutlier5iqr <- ukb_lm_univariate_fi_allcov_168metabolomics_nooutlier5iqr$pvalue.fi.allcov.nooutlier5iqr<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fi_allcov_168metabolomics_nooutlier5iqr <- ukb_z_fi_allcov_168metabolomics_nooutlier5iqr %>% group_by(variable) %>% 
  do(data.frame(r2.fi.allcov.nooutlier5iqr = summary(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fi_allcov_168metabolomics_nooutlier5iqr,. , by=c("variable"="variable"))
ukb_lm_univariate_fi_allcov_168metabolomics_nooutlier5iqr <- ukb_lm_univariate_fi_allcov_168metabolomics_nooutlier5iqr %>% mutate(rank_r2.fi.allcov.nooutlier5iqr = rank(-r2.fi.allcov.nooutlier5iqr)) # An indicator showing the rank of r2, from largest to smallest

### 32 biomarkers
ukb_lm_univariate_fi_allcov_32biomarkers_nooutlier5iqr <- ukb_z_fi_allcov_32biomarkers_nooutlier5iqr %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.nooutlier5iqr = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.nooutlier5iqr = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.nooutlier5iqr = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.nooutlier5iqr, se.fi.allcov.nooutlier5iqr, pvalue.fi.allcov.nooutlier5iqr) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_32biomarkers_nooutlier5iqr$bon.p.fi.allcov.nooutlier5iqr <- ukb_lm_univariate_fi_allcov_32biomarkers_nooutlier5iqr$pvalue.fi.allcov.nooutlier5iqr<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fi_allcov_32biomarkers_nooutlier5iqr <- ukb_z_fi_allcov_32biomarkers_nooutlier5iqr %>% group_by(variable) %>% 
  do(data.frame(r2.fi.allcov.nooutlier5iqr = summary(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fi_allcov_32biomarkers_nooutlier5iqr,. , by=c("variable"="variable"))
ukb_lm_univariate_fi_allcov_32biomarkers_nooutlier5iqr <- ukb_lm_univariate_fi_allcov_32biomarkers_nooutlier5iqr %>% mutate(rank_r2.fi.allcov.nooutlier5iqr = rank(-r2.fi.allcov.nooutlier5iqr)) # An indicator showing the rank of r2, from largest to smallest

### Summary of results
rm(list=c("ukb_z_fi_allcov_168metabolomics_nooutlier5iqr","ukb_z_fi_allcov_32biomarkers_nooutlier5iqr"))
clipr::write_clip(ukb_lm_univariate_fi_allcov_168metabolomics_nooutlier5iqr)
clipr::write_clip(ukb_lm_univariate_fi_allcov_32biomarkers_nooutlier5iqr)
write.table(ukb_lm_univariate_fi_allcov_168metabolomics_nooutlier5iqr, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_nooutlier5iqr.txt")
write.table(ukb_lm_univariate_fi_allcov_32biomarkers_nooutlier5iqr, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_nooutlier5iqr.txt")


#=============== SENSITIVITY ANALYSIS 2: USING FP AS OUTCOME =================

##### Univariate linear regression models between FP and each biomarker, adjusted for age and sex #####

ukb_z_fp_age_sex_168metabolomics <- melt(ukb_z_168metabolomics_complete[,c("eid","age","sex","fp_total",ukb_metabolomics_list)], id.vars = c("eid", "fp_total", "age", "sex")) 
ukb_z_fp_age_sex_32biomarkers <- melt(ukb_z_32biomarkers_complete[,c("eid","age","sex","fp_total",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fp_total", "age", "sex")) 

### 168 metabolomics
ukb_lm_univariate_fp_age_sex_168metabolomics <- ukb_z_fp_age_sex_168metabolomics %>% 
  group_by(variable) %>% 
  do(tidy(lm(fp_total ~ value + age + sex, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fp.age.sex = estimate, # Beta-coefficients from linear regression model
         se.fp.age.sex = std.error,  # Standard errors from linear regression model
         pvalue.fp.age.sex = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fp.age.sex, se.fp.age.sex, pvalue.fp.age.sex) %>% 
  as.data.frame()
ukb_lm_univariate_fp_age_sex_168metabolomics$bon.p.fp.age.sex <- ukb_lm_univariate_fp_age_sex_168metabolomics$pvalue.fp.age.sex<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fp_age_sex_168metabolomics <- ukb_z_fp_age_sex_168metabolomics %>% group_by(variable) %>% 
  do(data.frame(r2.fp.age.sex = summary(lm(fp_total ~ value + age + sex, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fp_age_sex_168metabolomics,. , by=c("variable"="variable"))
ukb_lm_univariate_fp_age_sex_168metabolomics <- ukb_lm_univariate_fp_age_sex_168metabolomics %>% mutate(rank_r2.fp.age.sex = rank(-r2.fp.age.sex)) # An indicator showing the rank of r2, from largest to smallest

### 32 biomarkers
ukb_lm_univariate_fp_age_sex_32biomarkers <- ukb_z_fp_age_sex_32biomarkers %>% 
  group_by(variable) %>% 
  do(tidy(lm(fp_total ~ value + age + sex, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fp.age.sex = estimate, # Beta-coefficients from linear regression model
         se.fp.age.sex = std.error,  # Standard errors from linear regression model
         pvalue.fp.age.sex = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fp.age.sex, se.fp.age.sex, pvalue.fp.age.sex) %>% 
  as.data.frame()
ukb_lm_univariate_fp_age_sex_32biomarkers$bon.p.fp.age.sex <- ukb_lm_univariate_fp_age_sex_32biomarkers$pvalue.fp.age.sex<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fp_age_sex_32biomarkers <- ukb_z_fp_age_sex_32biomarkers %>% group_by(variable) %>% 
  do(data.frame(r2.fp.age.sex = summary(lm(fp_total ~ value + age + sex, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fp_age_sex_32biomarkers,. , by=c("variable"="variable"))
ukb_lm_univariate_fp_age_sex_32biomarkers <- ukb_lm_univariate_fp_age_sex_32biomarkers %>% mutate(rank_r2.fp.age.sex = rank(-r2.fp.age.sex)) # An indicator showing the rank of r2, from largest to smallest

### Summary of results
rm(list=c("ukb_z_fp_age_sex_168metabolomics","ukb_z_fp_age_sex_32biomarkers"))
clipr::write_clip(ukb_lm_univariate_fp_age_sex_168metabolomics)
clipr::write_clip(ukb_lm_univariate_fp_age_sex_32biomarkers)
write.table(ukb_lm_univariate_fp_age_sex_168metabolomics, file = "Output/Observational_results/ukb_lm_univariate_fp_age_sex_168metabolomics.txt")
write.table(ukb_lm_univariate_fp_age_sex_32biomarkers, file = "Output/Observational_results/ukb_lm_univariate_fp_age_sex_32biomarkers.txt")


##### Univariate linear regression models between FP and each biomarker, adjusted for age, sex, assessment center, smoking, education, deprivation index, BMI #####

ukb_z_fp_allcov_168metabolomics <- melt(ukb_z_168metabolomics_complete[,c("eid","age","sex","fp_total","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_metabolomics_list)], id.vars = c("eid", "fp_total", "age", "sex","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 
ukb_z_fp_allcov_32biomarkers <- melt(ukb_z_32biomarkers_complete[,c("eid","age","sex","fp_total","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fp_total", "age", "sex","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 

### 168 metabolomics
ukb_lm_univariate_fp_allcov_168metabolomics <- ukb_z_fp_allcov_168metabolomics %>% 
  group_by(variable) %>% 
  do(tidy(lm(fp_total ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fp.allcov = estimate, # Beta-coefficients from linear regression model
         se.fp.allcov = std.error,  # Standard errors from linear regression model
         pvalue.fp.allcov = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fp.allcov, se.fp.allcov, pvalue.fp.allcov) %>% 
  as.data.frame()
ukb_lm_univariate_fp_allcov_168metabolomics$bon.p.fp.allcov <- ukb_lm_univariate_fp_allcov_168metabolomics$pvalue.fp.allcov<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fp_allcov_168metabolomics <- ukb_z_fp_allcov_168metabolomics %>% group_by(variable) %>% 
  do(data.frame(r2.fp.allcov = summary(lm(fp_total ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fp_allcov_168metabolomics,. , by=c("variable"="variable"))
ukb_lm_univariate_fp_allcov_168metabolomics <- ukb_lm_univariate_fp_allcov_168metabolomics %>% mutate(rank_r2.fp.allcov = rank(-r2.fp.allcov)) # An indicator showing the rank of r2, from largest to smallest

### 32 biomarkers
ukb_lm_univariate_fp_allcov_32biomarkers <- ukb_z_fp_allcov_32biomarkers %>% 
  group_by(variable) %>% 
  do(tidy(lm(fp_total ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fp.allcov = estimate, # Beta-coefficients from linear regression model
         se.fp.allcov = std.error,  # Standard errors from linear regression model
         pvalue.fp.allcov = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fp.allcov, se.fp.allcov, pvalue.fp.allcov) %>% 
  as.data.frame()
ukb_lm_univariate_fp_allcov_32biomarkers$bon.p.fp.allcov <- ukb_lm_univariate_fp_allcov_32biomarkers$pvalue.fp.allcov<.05/200 # An indicator showing whether the metabolite has p<.05/200
# Calculate R2 of the models
ukb_lm_univariate_fp_allcov_32biomarkers <- ukb_z_fp_allcov_32biomarkers %>% group_by(variable) %>% 
  do(data.frame(r2.fp.allcov = summary(lm(fp_total ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))$r.squared)) %>% as.data.frame() %>% 
  left_join(ukb_lm_univariate_fp_allcov_32biomarkers,. , by=c("variable"="variable"))
ukb_lm_univariate_fp_allcov_32biomarkers <- ukb_lm_univariate_fp_allcov_32biomarkers %>% mutate(rank_r2.fp.allcov = rank(-r2.fp.allcov)) # An indicator showing the rank of r2, from largest to smallest

### Summary of results
rm(list=c("ukb_z_fp_allcov_168metabolomics","ukb_z_fp_allcov_32biomarkers"))
clipr::write_clip(ukb_lm_univariate_fp_allcov_168metabolomics)
clipr::write_clip(ukb_lm_univariate_fp_allcov_32biomarkers)
write.table(ukb_lm_univariate_fp_allcov_168metabolomics, file = "Output/Observational_results/ukb_lm_univariate_fp_allcov_168metabolomics.txt")
write.table(ukb_lm_univariate_fp_allcov_32biomarkers, file = "Output/Observational_results/ukb_lm_univariate_fp_allcov_32biomarkers.txt")


#======================= SUBGROUP ANALYSES ===================================

##### Subgroup analysis for individuals aged above and below 60 years #####

# Indicator showing whether the person is < or >=60 years
ukb_z_168metabolomics_complete$age60 <- ukb_z_168metabolomics_complete$age >= 60
ukb_z_32biomarkers_complete$age60 <- ukb_z_32biomarkers_complete$age >= 60

# Univariate linear regression models between FI and each biomarker, adjusted for age, sex, assessment center, smoking, education, deprivation index, BMI
# Age 60-
ukb_z_fi_allcov_168metabolomics_age60minus <- melt(ukb_z_168metabolomics_complete[!ukb_z_168metabolomics_complete$age60,c("eid","age","sex","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_metabolomics_list)], id.vars = c("eid", "fi_percent", "age", "sex","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 
ukb_z_fi_allcov_32biomarkers_age60minus <- melt(ukb_z_32biomarkers_complete[!ukb_z_32biomarkers_complete$age60,c("eid","age","sex","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fi_percent", "age", "sex","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 

ukb_lm_univariate_fi_allcov_168metabolomics_age60minus <- ukb_z_fi_allcov_168metabolomics_age60minus %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.age60minus = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.age60minus = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.age60minus = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.age60minus, se.fi.allcov.age60minus, pvalue.fi.allcov.age60minus) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_168metabolomics_age60minus$bon.p.fi.allcov.age60minus <- ukb_lm_univariate_fi_allcov_168metabolomics_age60minus$pvalue.fi.allcov.age60minus<.05/200 # An indicator showing whether the metabolite has p<.05/200
write.table(ukb_lm_univariate_fi_allcov_168metabolomics_age60minus, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_age60minus.txt")

ukb_lm_univariate_fi_allcov_32biomarkers_age60minus <- ukb_z_fi_allcov_32biomarkers_age60minus %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.age60minus = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.age60minus = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.age60minus = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.age60minus, se.fi.allcov.age60minus, pvalue.fi.allcov.age60minus) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_32biomarkers_age60minus$bon.p.fi.allcov.age60minus <- ukb_lm_univariate_fi_allcov_32biomarkers_age60minus$pvalue.fi.allcov.age60minus<.05/200 # An indicator showing whether the metabolite has p<.05/200
write.table(ukb_lm_univariate_fi_allcov_32biomarkers_age60minus, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_age60minus.txt")
rm(list=c("ukb_z_fi_allcov_168metabolomics_age60minus","ukb_z_fi_allcov_32biomarkers_age60minus"))

# Age 60+
ukb_z_fi_allcov_168metabolomics_age60plus <- melt(ukb_z_168metabolomics_complete[ukb_z_168metabolomics_complete$age60,c("eid","age","sex","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_metabolomics_list)], id.vars = c("eid", "fi_percent", "age", "sex","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 
ukb_z_fi_allcov_32biomarkers_age60plus <- melt(ukb_z_32biomarkers_complete[ukb_z_32biomarkers_complete$age60,c("eid","age","sex","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fi_percent", "age", "sex","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 

ukb_lm_univariate_fi_allcov_168metabolomics_age60plus <- ukb_z_fi_allcov_168metabolomics_age60plus %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.age60plus = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.age60plus = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.age60plus = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.age60plus, se.fi.allcov.age60plus, pvalue.fi.allcov.age60plus) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_168metabolomics_age60plus$bon.p.fi.allcov.age60plus <- ukb_lm_univariate_fi_allcov_168metabolomics_age60plus$pvalue.fi.allcov.age60plus<.05/200 # An indicator showing whether the metabolite has p<.05/200
write.table(ukb_lm_univariate_fi_allcov_168metabolomics_age60plus, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_age60plus.txt")

ukb_lm_univariate_fi_allcov_32biomarkers_age60plus <- ukb_z_fi_allcov_32biomarkers_age60plus %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.age60plus = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.age60plus = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.age60plus = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.age60plus, se.fi.allcov.age60plus, pvalue.fi.allcov.age60plus) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_32biomarkers_age60plus$bon.p.fi.allcov.age60plus <- ukb_lm_univariate_fi_allcov_32biomarkers_age60plus$pvalue.fi.allcov.age60plus<.05/200 # An indicator showing whether the metabolite has p<.05/200
write.table(ukb_lm_univariate_fi_allcov_32biomarkers_age60plus, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_age60plus.txt")
rm(list=c("ukb_z_fi_allcov_168metabolomics_age60plus","ukb_z_fi_allcov_32biomarkers_age60plus"))


##### Subgroup analysis for women and men #####

# Univariate linear regression models between FI and each biomarker, adjusted for age, assessment center, smoking, education, deprivation index, BMI
# Women
ukb_z_fi_allcov_168metabolomics_women <- melt(ukb_z_168metabolomics_complete[ukb_z_168metabolomics_complete$sex=="Women",c("eid","age","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_metabolomics_list)], id.vars = c("eid", "fi_percent", "age","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 
ukb_z_fi_allcov_32biomarkers_women <- melt(ukb_z_32biomarkers_complete[ukb_z_32biomarkers_complete$sex=="Women",c("eid","age","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fi_percent", "age","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 

ukb_lm_univariate_fi_allcov_168metabolomics_women <- ukb_z_fi_allcov_168metabolomics_women %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.women = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.women = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.women = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.women, se.fi.allcov.women, pvalue.fi.allcov.women) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_168metabolomics_women$bon.p.fi.allcov.women <- ukb_lm_univariate_fi_allcov_168metabolomics_women$pvalue.fi.allcov.women<.05/200 # An indicator showing whether the metabolite has p<.05/200
write.table(ukb_lm_univariate_fi_allcov_168metabolomics_women, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_women.txt")

ukb_lm_univariate_fi_allcov_32biomarkers_women <- ukb_z_fi_allcov_32biomarkers_women %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.women = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.women = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.women = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.women, se.fi.allcov.women, pvalue.fi.allcov.women) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_32biomarkers_women$bon.p.fi.allcov.women <- ukb_lm_univariate_fi_allcov_32biomarkers_women$pvalue.fi.allcov.women<.05/200 # An indicator showing whether the metabolite has p<.05/200
write.table(ukb_lm_univariate_fi_allcov_32biomarkers_women, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_women.txt")
rm(list=c("ukb_z_fi_allcov_168metabolomics_women","ukb_z_fi_allcov_32biomarkers_women"))

# Men
ukb_z_fi_allcov_168metabolomics_men <- melt(ukb_z_168metabolomics_complete[ukb_z_168metabolomics_complete$sex=="Men",c("eid","age","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_metabolomics_list)], id.vars = c("eid", "fi_percent", "age","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 
ukb_z_fi_allcov_32biomarkers_men <- melt(ukb_z_32biomarkers_complete[ukb_z_32biomarkers_complete$sex=="Men",c("eid","age","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fi_percent", "age","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 

ukb_lm_univariate_fi_allcov_168metabolomics_men <- ukb_z_fi_allcov_168metabolomics_men %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.men = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.men = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.men = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.men, se.fi.allcov.men, pvalue.fi.allcov.men) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_168metabolomics_men$bon.p.fi.allcov.men <- ukb_lm_univariate_fi_allcov_168metabolomics_men$pvalue.fi.allcov.men<.05/200 # An indicator showing whether the metabolite has p<.05/200
write.table(ukb_lm_univariate_fi_allcov_168metabolomics_men, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_men.txt")

ukb_lm_univariate_fi_allcov_32biomarkers_men <- ukb_z_fi_allcov_32biomarkers_men %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.men = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.men = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.men = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.men, se.fi.allcov.men, pvalue.fi.allcov.men) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_32biomarkers_men$bon.p.fi.allcov.men <- ukb_lm_univariate_fi_allcov_32biomarkers_men$pvalue.fi.allcov.men<.05/200 # An indicator showing whether the metabolite has p<.05/200
write.table(ukb_lm_univariate_fi_allcov_32biomarkers_men, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_men.txt")
rm(list=c("ukb_z_fi_allcov_168metabolomics_men","ukb_z_fi_allcov_32biomarkers_men"))

 
##### Subgroup analysis for non-white ethnicity #####

ukb_nonwhite_168metabolomics_complete <- ukb_nonwhite[ complete.cases( ukb_nonwhite[, c(ukb_metabolomics_list)] ), ] 
ukb_nonwhite_168metabolomics_complete <- ukb_nonwhite_168metabolomics_complete[!(ukb_nonwhite_168metabolomics_complete$eid %in% ukb_nonwhite_168metabolomics_complete[which(ukb_nonwhite_168metabolomics_complete$QCflagged==1), "eid"]),] # Exclude those with any QC flags from dataset, N = 4785
ukb_nonwhite_z_168metabolomics_complete <- cbind(ukb_nonwhite_168metabolomics_complete[,c("eid","age","sex","fi_percent","fp_total","assessment_centre","smoking","alcohol_3group","education","di","bmi")],
                                            scale(ukb_nonwhite_168metabolomics_complete[,c(ukb_metabolomics_list)], center=T, scale=T)) # Standardization of metabolites

ukb_nonwhite_32biomarkers_complete <- ukb_nonwhite[ complete.cases( ukb_nonwhite[, c(ukb_blood_urine_biomarkers_list)] ), ] 
ukb_nonwhite_z_32biomarkers_complete <- cbind(ukb_nonwhite_32biomarkers_complete[,c("eid","age","sex","fi_percent","fp_total","assessment_centre","smoking","alcohol_3group","education","di","bmi")],
                                     scale(ukb_nonwhite_32biomarkers_complete[,c(ukb_blood_urine_biomarkers_list)], center=T, scale=T)) # Standardization of metabolites

# Univariate linear regression models between FI and each biomarker, adjusted for age, assessment center, smoking, education, deprivation index, BMI
ukb_z_fi_allcov_168metabolomics_nonwhite <- melt(ukb_nonwhite_z_168metabolomics_complete[,c("eid","age","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_metabolomics_list)], id.vars = c("eid", "fi_percent", "age","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 
ukb_z_fi_allcov_32biomarkers_nonwhite <- melt(ukb_nonwhite_z_32biomarkers_complete[,c("eid","age","fi_percent","assessment_centre","smoking","alcohol_3group","education","di","bmi",ukb_blood_urine_biomarkers_list)], id.vars = c("eid", "fi_percent", "age","assessment_centre","smoking","alcohol_3group","education","di","bmi")) 

ukb_lm_univariate_fi_allcov_168metabolomics_nonwhite <- ukb_z_fi_allcov_168metabolomics_nonwhite %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.nonwhite = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.nonwhite = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.nonwhite = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.nonwhite, se.fi.allcov.nonwhite, pvalue.fi.allcov.nonwhite) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_168metabolomics_nonwhite$bon.p.fi.allcov.nonwhite <- ukb_lm_univariate_fi_allcov_168metabolomics_nonwhite$pvalue.fi.allcov.nonwhite<.05/200 # An indicator showing whether the metabolite has p<.05/200
write.table(ukb_lm_univariate_fi_allcov_168metabolomics_nonwhite, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_168metabolomics_nonwhite.txt")

ukb_lm_univariate_fi_allcov_32biomarkers_nonwhite <- ukb_z_fi_allcov_32biomarkers_nonwhite %>% 
  group_by(variable) %>% 
  do(tidy(lm(fi_percent ~ value + age + assessment_centre + smoking + alcohol_3group + education + di + bmi, .))) %>% 
  filter(term == "value") %>% 
  mutate(beta.fi.allcov.nonwhite = estimate, # Beta-coefficients from linear regression model
         se.fi.allcov.nonwhite = std.error,  # Standard errors from linear regression model
         pvalue.fi.allcov.nonwhite = p.value # P-values from linear regression model
  ) %>% 
  select(beta.fi.allcov.nonwhite, se.fi.allcov.nonwhite, pvalue.fi.allcov.nonwhite) %>% 
  as.data.frame()
ukb_lm_univariate_fi_allcov_32biomarkers_nonwhite$bon.p.fi.allcov.nonwhite <- ukb_lm_univariate_fi_allcov_32biomarkers_nonwhite$pvalue.fi.allcov.nonwhite<.05/200 # An indicator showing whether the metabolite has p<.05/200
write.table(ukb_lm_univariate_fi_allcov_32biomarkers_nonwhite, file = "Output/Observational_results/ukb_lm_univariate_fi_allcov_32biomarkers_nonwhite.txt")
rm(list=c("ukb_z_fi_allcov_168metabolomics_nonwhite","ukb_z_fi_allcov_32biomarkers_nonwhite"))

 
##### Subgroup analysis for creatinine-FI on individuals with and without chronic kidney disease #####

ukb_ckd <- merge(ukb, 
                 read.delim("Data/Raw_UKB_data/ukb_ckd.txt", sep=" "), # A previously prepared dataset showing whether the person has prevalent CKD (based on self-reported kidney disease (coded as "1192") and ICD-10 code (N17 ("acute kidney failure"), N18 ("chronic kidney disease (CKD)), N19 ("unspecified kidney failure)) from hospital inpatient data)
                 by="eid",all.x=T)
ukb_ckd_168metabolomics_complete <- ukb_ckd[ complete.cases( ukb_ckd[, c(ukb_metabolomics_list)] ), ] 
ukb_ckd_168metabolomics_complete <- ukb_ckd_168metabolomics_complete[!(ukb_ckd_168metabolomics_complete$eid %in% QCflagged_168metabolomics_id),] # Exclude those with any QC flags from dataset, N = 90,573
ukb_ckd_z_168metabolomics_complete <- cbind(ukb_ckd_168metabolomics_complete[,c("eid","age","sex","fi_percent","fp_total","assessment_centre","smoking","alcohol_3group","education","di","bmi","ckd")],
                                            scale(ukb_ckd_168metabolomics_complete[,c(ukb_metabolomics_list)], center=T, scale=T)) # Standardization of metabolites

# Association between creatinine and FI, stratified by CKD
table(ukb_ckd_z_168metabolomics_complete$ckd) # 0=85699, 1=5974
prop.table(table(ukb_ckd_z_168metabolomics_complete$ckd))
ukb_lm_creatinine_fi_ckd <- cbind(
  data.frame(model=c("No CKD","CKD"),
             n=as.vector(table(ukb_ckd_z_168metabolomics_complete$ckd))),
  rbind(summary(lm(fi_percent ~ Creatinine + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ckd_z_168metabolomics_complete[ukb_ckd_z_168metabolomics_complete$ckd==0,]))$coefficient[2,c(1,2,4)],
        summary(lm(fi_percent ~ Creatinine + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ckd_z_168metabolomics_complete[ukb_ckd_z_168metabolomics_complete$ckd==1,]))$coefficient[2,c(1,2,4)])
)
rm(list=c("ukb_ckd","ukb_ckd_168metabolomics_complete","ukb_ckd_z_168metabolomics_complete"))


##### Subgroup analysis for GlycA-FI by LDL & CRP levels #####

# Correlations between GlycA and LDL
cor(ukb[,c("GlycA","LDL_C","LDL_serum","CRP_serum")], use="complete.obs")

# Adjust LDL levels for statin use (divide LDL by 0.684 if the participant is taking statin at baseline)
statin_id <- read.delim("Data/Raw_UKB_data/ukb_ki_fi_on_statin.dat") # Individuals who is taking statin at baseline (field 20003, defined by one of these 7 drugs 1140861958, 1140888594, 1140888648, 1141146234, 1141192410, 1140861922, 1141146138)
dim(statin_id)
table(ukb$eid %in% statin_id[,1])

# CHD prevalence
#ukb_chd <- read_dta("Data/Raw_UKB_data/Frailty_covariate/ukb_fi.dta") %>% select(eid, i9_chd)
ukb_chd_id <- read.delim("Data/Raw_UKB_data/ukb_chd_icd10_prevalence_only.dat") # Consider CHD prevalence cases as those with ICD codes I20-I25
ukb_chd_bypass_id <- read.delim("Data/Raw_UKB_data/ukb_chd_icd10_bypass_prevalence_only.dat") # Consider CHD prevalence cases as those with ICD codes I20-I25, Z95

# Create LDL and CRP subgroups
ukb_ldl <- mutate(ukb, 
                  statin = eid %in% statin_id[,1],
                  chd = eid %in% ukb_chd_id[,1],
                  chd_bypass = eid %in% ukb_chd_bypass_id[,1],
                  ldl_adjusted = ifelse(statin, LDL_serum/.684, LDL_serum),
                  high_ldl = ldl_adjusted>4.9,
                  low_ldl = ldl_adjusted<2.5,
                  ldl_cat = factor(case_when(
                    ldl_adjusted<2.59 ~ 1,
                    ldl_adjusted>=2.59 & ldl_adjusted<3.36 ~ 2,
                    ldl_adjusted>=3.36 & ldl_adjusted<4.14 ~ 3,
                    ldl_adjusted>=4.14 & ldl_adjusted<4.92 ~ 4,
                    ldl_adjusted>=4.92 ~ 5), 
                    labels = c("LDL <2.59 mmol/L (<100 mg/dL)","LDL 2.60-3.36 mmol/L (100-129 mg/dL)","LDL 3.37-4.14 mmol/L (130-159 mg/dL)","LDL 4.15-4.92 mmol/L (160-189 mg/dL)","LDL >=4.92 mmol/L (>=190 mg/dL)")),
                  high_crp = CRP_serum>=10,
                  crp_cat = factor(case_when(
                    CRP_serum<1 ~ 1,
                    CRP_serum>=1 & CRP_serum<3 ~ 2,
                    CRP_serum>=3 & CRP_serum<10 ~ 3,
                    CRP_serum>=10 ~ 4),
                    labels = c("CRP <1 mg/L","CRP 1-3 mg/L","CRP 3-10 mg/L","CRP >=10 mg/L")))
ukb_ldl_168metabolomics_complete <- ukb_ldl[ complete.cases( ukb_ldl[, c(ukb_metabolomics_list)] ), ] 
ukb_ldl_168metabolomics_complete <- ukb_ldl_168metabolomics_complete[!(ukb_ldl_168metabolomics_complete$eid %in% QCflagged_168metabolomics_id),] # Exclude those with any QC flags from dataset, N = 90,573
ukb_ldl_z_168metabolomics_complete <- cbind(ukb_ldl_168metabolomics_complete[,c("eid","age","sex","fi_percent","fp_total","assessment_centre","smoking","alcohol_3group","education","di","bmi","statin","ldl_adjusted","high_ldl","low_ldl","ldl_cat","CRP_serum","high_crp","crp_cat","chd","chd_bypass")],
                                            scale(ukb_ldl_168metabolomics_complete[,c(ukb_metabolomics_list)], center=T, scale=T)) # Standardization of metabolites

# Association between GlycA and FI, stratified by LDL subgroups
table(ukb_ldl_z_168metabolomics_complete$ldl_cat)
prop.table(table(ukb_ldl_z_168metabolomics_complete$ldl_cat))
table(ukb_ldl_z_168metabolomics_complete$ldl_cat, ukb_ldl_z_168metabolomics_complete$chd)
prop.table(table(ukb_ldl_z_168metabolomics_complete$ldl_cat, ukb_ldl_z_168metabolomics_complete$chd), 1)
table(ukb_ldl_z_168metabolomics_complete$ldl_cat, ukb_ldl_z_168metabolomics_complete$chd_bypass)
prop.table(table(ukb_ldl_z_168metabolomics_complete$ldl_cat, ukb_ldl_z_168metabolomics_complete$chd_bypass), 1)
ukb_lm_GlycA_fi_ldl <- cbind(
  data.frame(model=c("<2.59 mmol/L","2.60–3.36 mmol/L","3.37–4.14 mmol/L","4.15–4.92 mmol/L","≥4.92 mmol/L"),
             n=as.vector(table(ukb_ldl_z_168metabolomics_complete$ldl_cat))),
  rbind(summary(lm(fi_percent ~ GlycA + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ldl_z_168metabolomics_complete[ukb_ldl_z_168metabolomics_complete$ldl_cat=="LDL <2.59 mmol/L (<100 mg/dL)",]))$coefficient[2,c(1,2,4)],
        summary(lm(fi_percent ~ GlycA + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ldl_z_168metabolomics_complete[ukb_ldl_z_168metabolomics_complete$ldl_cat=="LDL 2.60-3.36 mmol/L (100-129 mg/dL)",]))$coefficient[2,c(1,2,4)],
        summary(lm(fi_percent ~ GlycA + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ldl_z_168metabolomics_complete[ukb_ldl_z_168metabolomics_complete$ldl_cat=="LDL 3.37-4.14 mmol/L (130-159 mg/dL)",]))$coefficient[2,c(1,2,4)],
        summary(lm(fi_percent ~ GlycA + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ldl_z_168metabolomics_complete[ukb_ldl_z_168metabolomics_complete$ldl_cat=="LDL 4.15-4.92 mmol/L (160-189 mg/dL)",]))$coefficient[2,c(1,2,4)],
        summary(lm(fi_percent ~ GlycA + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ldl_z_168metabolomics_complete[ukb_ldl_z_168metabolomics_complete$ldl_cat=="LDL >=4.92 mmol/L (>=190 mg/dL)",]))$coefficient[2,c(1,2,4)])
  )

# Association between GlycA and FI, stratified by CRP subgroups
table(ukb_ldl_z_168metabolomics_complete$crp_cat)
prop.table(table(ukb_ldl_z_168metabolomics_complete$crp_cat))
table(ukb_ldl_z_168metabolomics_complete$crp_cat, ukb_ldl_z_168metabolomics_complete$chd)
prop.table(table(ukb_ldl_z_168metabolomics_complete$crp_cat, ukb_ldl_z_168metabolomics_complete$chd), 1)
table(ukb_ldl_z_168metabolomics_complete$crp_cat, ukb_ldl_z_168metabolomics_complete$chd_bypass)
prop.table(table(ukb_ldl_z_168metabolomics_complete$crp_cat, ukb_ldl_z_168metabolomics_complete$chd_bypass), 1)
ukb_lm_GlycA_fi_crp <- cbind(
  data.frame(model=c("<1 mg/L","1–3 mg/L","3–10 mg/L","≥10 mg/L"),
             n=as.vector(table(ukb_ldl_z_168metabolomics_complete$crp_cat))),
  rbind(summary(lm(fi_percent ~ GlycA + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ldl_z_168metabolomics_complete[ukb_ldl_z_168metabolomics_complete$crp_cat=="CRP <1 mg/L",]))$coefficient[2,c(1,2,4)],
        summary(lm(fi_percent ~ GlycA + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ldl_z_168metabolomics_complete[ukb_ldl_z_168metabolomics_complete$crp_cat=="CRP 1-3 mg/L",]))$coefficient[2,c(1,2,4)],
        summary(lm(fi_percent ~ GlycA + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ldl_z_168metabolomics_complete[ukb_ldl_z_168metabolomics_complete$crp_cat=="CRP 3-10 mg/L",]))$coefficient[2,c(1,2,4)],
        summary(lm(fi_percent ~ GlycA + age + sex + assessment_centre + smoking + alcohol_3group + education + di + bmi, ukb_ldl_z_168metabolomics_complete[ukb_ldl_z_168metabolomics_complete$crp_cat=="CRP >=10 mg/L",]))$coefficient[2,c(1,2,4)])
)

rm(list=c("statin_id","ukb_chd_bypass_id","ukb_ldl","ukb_ldl_168metabolomics_complete","ukb_ldl_z_168metabolomics_complete"))


#==================== SAVE DATA FOR FURTHER ANALYSIS =========================

save.image("Data/R_data/UKB_biomarkers.Rdata")
#write.table(ukb_z, file = "Data/Cleaned_data/UKB_metabolomics_sd.txt", sep = "\t", row.names = F, na="")

# =============================== END OF FILE  ===============================