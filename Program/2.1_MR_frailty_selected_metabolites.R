#==============================================================================
# FILENAME: 2.1_MR_frailty_metabolomics.R
# PROJECT: 	Metabolomics_frailty
# PURPOSE:  To perform main MR analysis for the selected metabolites and frailty
# AUTHOR:   Jonathan Mak
# CREATED:	2022-07-05
# UPDATED: 	2023-02-21
# R VERSION: 4.1.3
#==============================================================================
### Note:
# A total of 37 replicated biomarkers and 12 biomarkers not available in replication 
# cohorts were selected for MR analysis

### Required packages
library(dplyr)
library(ggforestplot) # For forest plots of metabolites
library(ggplot2); library(patchwork); library(ggforce) # Plotting graphs
library(TwoSampleMR) # To perform MR analysis
library(MRPRESSO) # MR PRESSO


#========= SELECTION OF INSTRUMENTAL VARIABLES FOR NMR METABOLOMICS ===========
### Note: 
# All SNP data for metabolomics were from UK Biobank (maximum n=115,078):
# https://gwas.mrcieu.ac.uk/datasets/?gwas_id__icontains=met-d

# List of metabolomics to be tested in MR
mr_metabolomics <- c("met-d-Ala",
                     "met-d-Creatinine",
                     "met-d-Glucose",
                     "met-d-GlycA",
                     "met-d-HDL_size",
                     "met-d-IDL_CE",
                     "met-d-IDL_FC",
                     "met-d-LA",
                     "met-d-IDL_P",
                     "met-d-M_HDL_TG",
                     "met-d-M_LDL_CE",
                     "met-d-M_LDL_P",
                     "met-d-MUFA",
                     "met-d-Omega_6",
                     "met-d-Phe",
                     "met-d-Phosphatidylc",
                     "met-d-S_LDL_CE",
                     "met-d-S_LDL_L",
                     "met-d-S_LDL_P",
                     "met-d-S_LDL_PL",
                     "met-d-S_VLDL_TG",
                     "met-d-Sphingomyelins",
                     "met-d-VLDL_size",
                     "met-d-XL_HDL_FC",
                     "met-d-XS_VLDL_CE",
                     "met-d-XXL_VLDL_TG")

metabolomics_iv <- extract_instruments(mr_metabolomics, clump = TRUE, p1 = 5e-08, r2 = 0.001)

# Change names of the exposures
metabolomics_iv$exposure <- sub(".*id:met-d-", "", metabolomics_iv$exposure)

# Add sample size information for each metabolite
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-Ala",]$samplesize.exposure <- 115074
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-Creatinine",]$samplesize.exposure <- 110058
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-Glucose",]$samplesize.exposure <- 114867
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-GlycA",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-HDL_size",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-IDL_CE",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-IDL_FC",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-LA",]$samplesize.exposure <- 114999
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-IDL_P",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-M_HDL_TG",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-M_LDL_CE",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-M_LDL_P",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-MUFA",]$samplesize.exposure <- 114999
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-Omega_6",]$samplesize.exposure <- 114999
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-Phe",]$samplesize.exposure <- 115025
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-Phosphatidylc",]$samplesize.exposure <- 114999
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-S_LDL_CE",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-S_LDL_L",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-S_LDL_P",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-S_LDL_PL",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-S_VLDL_TG",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-Sphingomyelins",]$samplesize.exposure <- 114999
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-VLDL_size",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-XL_HDL_FC",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-XS_VLDL_CE",]$samplesize.exposure <- 115078
metabolomics_iv[metabolomics_iv$id.exposure=="met-d-XXL_VLDL_TG",]$samplesize.exposure <- 115078

# Number of SNPs selected
print(table(metabolomics_iv$exposure))

# Variance in exposure explained by each variant
metabolomics_iv$r2.exposure <- get_r_from_pn(metabolomics_iv$pval.exposure, metabolomics_iv$samplesize.exposure)^2

# F-statistics of each variant, F = r2 * (N-2) / (1-r2)
metabolomics_iv$f_stat.exposure <- metabolomics_iv$r2.exposure * (metabolomics_iv$samplesize.exposure) / (1-metabolomics_iv$r2.exposure)
summary(metabolomics_iv$f_stat.exposure) # all >10; mean F-statistics=119.03 (range 29.76-917.42)


#======== SELECTION OF INSTRUMENTAL VARIABLES FOR CLINICAL BIOMARKERS ========

# List of clinical biomarkers to be tested in MR
mr_clinical_biomarkers <- c("CRP",
                            "Total cholesterol",
                            "HbA1c",
                            "LDL-C",
                            "Triglycerides",
                            "ApoB",
                            "ALP",
                            "Creatinine",
                            "Cystatin C",
                            "GGT",
                            "IGF-1",
                            "Phosphate",
                            "SHBG",
                            "Testosterone",
                            "Total bilirubin",
                            "Vitamin D",
                            "Microalbumin (urine)",
                            "Potassium (urine)")

########################### SNPs for HbA1c ####################################
# Note: 
# 1. Data were obtained from the MAGIC consortium: https://magicinvestigators.org/downloads/ or
#    http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004903/HbA1c_METAL_European.txt.gz
# 2. Only those of European (N=123,665) ancestry were used. HbA1c trait values are 
#    untransformed (per 1% increase) and adjusted for age, sex and 
#    study-specific covariates.
# 3. Publication:
#    Wheeler E, Leong A, Liu CT, et al. Impact of common genetic determinants of 
#    Hemoglobin A1c on type 2 diabetes risk and diagnosis in ancestrally diverse 
#    populations: A transethnic genome-wide meta-analysis. PLoS Med. 
#    2017;14(9):e1002383. Published 2017 Sep 12. doi:10.1371/journal.pmed.1002383

### Read-in summary statistics data
hba1c_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/HbA1c_METAL_European.txt.gz", header=T)

### Select genome-wide significant SNPs
hba1c_iv <- hba1c_gwas[hba1c_gwas$pvalue<5*10^(-8),] # 821 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(hba1c_iv)[c(1,9)] <- c("rsid","pval")
hba1c_iv <- ieugwasr::ld_clump(hba1c_iv, clump_r2=.001, pop="EUR") # 37 SNPs retained

### Rename variables
colnames(hba1c_iv)
hba1c_iv <- hba1c_iv[,c(1:9)]
colnames(hba1c_iv)[c(2,3,4,5,6,8,9)] <- c("chromosome","position","effect_allele",
                                          "other_allele","eaf","se","pvalue")
hba1c_iv$N <- 123665 # Sample size
hba1c_iv$metabolite <- "HbA1c"
hba1c_iv$units <- "%"
rm(hba1c_gwas)

###################### SNPs for total cholesterol #############################
# Note: 
# 1. Data were obtained from the Global Lipids Genetics Consortium: 
#    http://csg.sph.umich.edu/willer/public/lipids2013/
# 2. Only those of European (N=187,365) ancestry were used. Values are per SD (41.75071927 mg/dL) increase.
# 3. Publication:
#    Global Lipids Genetics Consortium. Discovery and refinement of loci associated 
#    with lipid levels. Nat Genet. 2013 Nov;45(11):1274-1283. doi: 10.1038/ng.2797. 
#    Epub 2013 Oct 6. PMID: 24097068; PMCID: PMC3838666.

chol_iv <- TwoSampleMR::extract_instruments(outcomes="ieu-a-301",p1 = 5e-08,r2 = 0.001) # 88 SNPs retained
chol_iv <- chol_iv[,c("SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","samplesize.exposure")]
chol_iv$metabolite <- "Total cholesterol"
colnames(chol_iv) <- c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite")
chol_iv$units <- "SD"

###################### SNPs for LDL cholesterol #############################
# Note: 
# 1. Data were obtained from the Global Lipids Genetics Consortium: 
#    http://csg.sph.umich.edu/willer/public/lipids2013/
# 2. Only those of European (N=173,082) ancestry were used. Values are per SD (38.6745912 mg/dL) increase.
# 3. Publication:
#    Global Lipids Genetics Consortium. Discovery and refinement of loci associated 
#    with lipid levels. Nat Genet. 2013 Nov;45(11):1274-1283. doi: 10.1038/ng.2797. 
#    Epub 2013 Oct 6. PMID: 24097068; PMCID: PMC3838666.

ldl_iv <- TwoSampleMR::extract_instruments(outcomes="ieu-a-300",p1 = 5e-08,r2 = 0.001) # 81 SNPs retained
ldl_iv <- ldl_iv[,c("SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","samplesize.exposure")]
ldl_iv$metabolite <- "LDL-C"
colnames(ldl_iv) <- c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite")
ldl_iv$units <- "SD"

####################### SNPs for Triglycerides ##############################
# Note: 
# 1. Data were obtained from the Global Lipids Genetics Consortium: 
#    http://csg.sph.umich.edu/willer/public/lipids2013/
# 2. Only those of European (N=177,861) ancestry were used. Values are per SD (90.72424257 mg/dL) increase.
# 3. Publication:
#    Global Lipids Genetics Consortium. Discovery and refinement of loci associated 
#    with lipid levels. Nat Genet. 2013 Nov;45(11):1274-1283. doi: 10.1038/ng.2797. 
#    Epub 2013 Oct 6. PMID: 24097068; PMCID: PMC3838666.

tg_iv <- TwoSampleMR::extract_instruments(outcomes="ieu-a-302",p1 = 5e-08,r2 = 0.001) # 55 SNPs retained
tg_iv <- tg_iv[,c("SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","samplesize.exposure")]
tg_iv$metabolite <- "Triglycerides"
colnames(tg_iv) <- c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite")
tg_iv$units <- "SD"

####################### SNPs for C-reactive protein ###########################
# Note: 
# 1. Data were obtained from the CHARGE Consortium: 
#    https://gwas.mrcieu.ac.uk/datasets/ieu-b-35/
# 2. Only those of European (N=204,402) ancestry were used. Values are per mg/L 
#    increase and natural log transformed
# 3. Publication:
#    Said S et al. Genetic analysis of over half a million people characterises 
#    C-reactive protein loci. Nat Commun. 2022 Apr 22;13(1):2198. 
#    doi: 10.1038/s41467-022-29650-5. PMID: 35459240; PMCID: PMC9033829.

crp_iv <- TwoSampleMR::extract_instruments(outcomes="ieu-b-35",p1 = 5e-08,r2 = 0.001) # 57 SNPs retained
crp_iv <- crp_iv[,c("SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","samplesize.exposure")]
crp_iv$metabolite <- "CRP"
colnames(crp_iv) <- c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite")
crp_iv$units <- "log mg/L"

####################### SNPs for ApoB ###########################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=193,604)
# 2. Values are per SD increase

### Read-in summary statistics data
ApoB_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/ApoB_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
ApoB_iv <- ApoB_gwas[ApoB_gwas$P<5*10^(-8),] # 14769 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(ApoB_iv)[c(2,10)] <- c("rsid","pval")
ApoB_iv <- ieugwasr::ld_clump(ApoB_iv, clump_r2=.001, pop="EUR") # 102 SNPs retained

### Rename variables
colnames(ApoB_iv)
ApoB_iv <- ApoB_iv[,c(1:10)]
colnames(ApoB_iv) <- c("chromosome","rsid","position","effect_allele",
                       "other_allele","N","eaf","beta","se","pvalue")
ApoB_iv$metabolite <- "ApoB"
ApoB_iv$units <- "SD"
ApoB_iv <- ApoB_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(ApoB_gwas)

############################ SNPs for ALP #####################################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=194,576)
# 2. Values are per SD increase

### Read-in summary statistics data
ALP_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/ALP_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
ALP_iv <- ALP_gwas[ALP_gwas$P<5*10^(-8),] # 23421 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(ALP_iv)[c(2,10)] <- c("rsid","pval")
ALP_iv <- ieugwasr::ld_clump(ALP_iv, clump_r2=.001, pop="EUR") # 129 SNPs retained

### Rename variables
colnames(ALP_iv)
ALP_iv <- ALP_iv[,c(1:10)]
colnames(ALP_iv) <- c("chromosome","rsid","position","effect_allele",
                      "other_allele","N","eaf","beta","se","pvalue")
ALP_iv$metabolite <- "ALP"
ALP_iv$units <- "SD"
ALP_iv <- ALP_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(ALP_gwas)

######################## SNPs for Cystatin C ##################################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=194,545)
# 2. Values are per SD increase

### Read-in summary statistics data
CysC_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/CysC_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
CysC_iv <- CysC_gwas[CysC_gwas$P<5*10^(-8),] # 19118 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(CysC_iv)[c(2,10)] <- c("rsid","pval")
CysC_iv <- ieugwasr::ld_clump(CysC_iv, clump_r2=.001, pop="EUR") # 100 SNPs retained

### Rename variables
colnames(CysC_iv)
CysC_iv <- CysC_iv[,c(1:10)]
colnames(CysC_iv) <- c("chromosome","rsid","position","effect_allele",
                       "other_allele","N","eaf","beta","se","pvalue")
CysC_iv$metabolite <- "Cystatin C"
CysC_iv$units <- "SD"
CysC_iv <- CysC_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(CysC_gwas)

############################ SNPs for GGT #####################################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=194,466)
# 2. Values are per SD increase

### Read-in summary statistics data
GGT_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/GGT_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
GGT_iv <- GGT_gwas[GGT_gwas$P<5*10^(-8),] # 5578 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(GGT_iv)[c(2,10)] <- c("rsid","pval")
GGT_iv <- ieugwasr::ld_clump(GGT_iv, clump_r2=.001, pop="EUR") # 73 SNPs retained

### Rename variables
colnames(GGT_iv)
GGT_iv <- GGT_iv[,c(1:10)]
colnames(GGT_iv) <- c("chromosome","rsid","position","effect_allele",
                      "other_allele","N","eaf","beta","se","pvalue")
GGT_iv$metabolite <- "GGT"
GGT_iv$units <- "SD"
GGT_iv <- GGT_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(GGT_gwas)

########################### SNPs for IGF-1 ####################################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=193,517)
# 2. Values are per SD increase

### Read-in summary statistics data
IGF1_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/IGF1_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
IGF1_iv <- IGF1_gwas[IGF1_gwas$P<5*10^(-8),] # 37853 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(IGF1_iv)[c(2,10)] <- c("rsid","pval")
IGF1_iv <- ieugwasr::ld_clump(IGF1_iv, clump_r2=.001, pop="EUR") # 193 SNPs retained

### Rename variables
colnames(IGF1_iv)
IGF1_iv <- IGF1_iv[,c(1:10)]
colnames(IGF1_iv) <- c("chromosome","rsid","position","effect_allele",
                       "other_allele","N","eaf","beta","se","pvalue")
IGF1_iv$metabolite <- "IGF-1"
IGF1_iv$units <- "SD"
IGF1_iv <- IGF1_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(IGF1_gwas)

######################### SNPs for Phosphate ##################################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=177,797)
# 2. Values are per SD increase

### Read-in summary statistics data
Phosphate_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/Phosphate_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
Phosphate_iv <- Phosphate_gwas[Phosphate_gwas$P<5*10^(-8),] # 10861 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(Phosphate_iv)[c(2,10)] <- c("rsid","pval")
Phosphate_iv <- ieugwasr::ld_clump(Phosphate_iv, clump_r2=.001, pop="EUR") # 69 SNPs retained

### Rename variables
colnames(Phosphate_iv)
Phosphate_iv <- Phosphate_iv[,c(1:10)]
colnames(Phosphate_iv) <- c("chromosome","rsid","position","effect_allele",
                            "other_allele","N","eaf","beta","se","pvalue")
Phosphate_iv$metabolite <- "Phosphate"
Phosphate_iv$units <- "SD"
Phosphate_iv <- Phosphate_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(Phosphate_gwas)

########################### SNPs for SHBG #####################################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=176,414)
# 2. Values are per SD increase

### Read-in summary statistics data
SHBG_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/SHBG_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
SHBG_iv <- SHBG_gwas[SHBG_gwas$P<5*10^(-8),] # 15456 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(SHBG_iv)[c(2,10)] <- c("rsid","pval")
SHBG_iv <- ieugwasr::ld_clump(SHBG_iv, clump_r2=.001, pop="EUR") # 122 SNPs retained

### Rename variables
colnames(SHBG_iv)
SHBG_iv <- SHBG_iv[,c(1:10)]
colnames(SHBG_iv) <- c("chromosome","rsid","position","effect_allele",
                       "other_allele","N","eaf","beta","se","pvalue")
SHBG_iv$metabolite <- "SHBG"
SHBG_iv$units <- "SD"
SHBG_iv <- SHBG_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(SHBG_gwas)

######################## SNPs for Testosterone ################################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=176,307)
# 2. Values are per SD increase

### Read-in summary statistics data
Testosterone_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/Testosterone_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
Testosterone_iv <- Testosterone_gwas[Testosterone_gwas$P<5*10^(-8),] # 6192 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(Testosterone_iv)[c(2,10)] <- c("rsid","pval")
Testosterone_iv <- ieugwasr::ld_clump(Testosterone_iv, clump_r2=.001, pop="EUR") # 60 SNPs retained

### Rename variables
colnames(Testosterone_iv)
Testosterone_iv <- Testosterone_iv[,c(1:10)]
colnames(Testosterone_iv) <- c("chromosome","rsid","position","effect_allele",
                               "other_allele","N","eaf","beta","se","pvalue")
Testosterone_iv$metabolite <- "Testosterone"
Testosterone_iv$units <- "SD"
Testosterone_iv <- Testosterone_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(Testosterone_gwas)

####################### SNPs for Total bilirubin ##############################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=193,741)
# 2. Values are per SD increase

### Read-in summary statistics data
Total_bil_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/Total_bil_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
Total_bil_iv <- Total_bil_gwas[Total_bil_gwas$P<5*10^(-8),] # 10859 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(Total_bil_iv)[c(2,10)] <- c("rsid","pval")
Total_bil_iv <- ieugwasr::ld_clump(Total_bil_iv, clump_r2=.001, pop="EUR") # 46 SNPs retained

### Rename variables
colnames(Total_bil_iv)
Total_bil_iv <- Total_bil_iv[,c(1:10)]
colnames(Total_bil_iv) <- c("chromosome","rsid","position","effect_allele",
                            "other_allele","N","eaf","beta","se","pvalue")
Total_bil_iv$metabolite <- "Total bilirubin"
Total_bil_iv$units <- "SD"
Total_bil_iv <- Total_bil_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(Total_bil_gwas)

######################### SNPs for Vitamin D ##################################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=185,948)
# 2. Values are per SD increase

### Read-in summary statistics data
VitD_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/VitD_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
VitD_iv <- VitD_gwas[VitD_gwas$P<5*10^(-8),] # 6304 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(VitD_iv)[c(2,10)] <- c("rsid","pval")
VitD_iv <- ieugwasr::ld_clump(VitD_iv, clump_r2=.001, pop="EUR") # 36 SNPs retained

### Rename variables
colnames(VitD_iv)
VitD_iv <- VitD_iv[,c(1:10)]
colnames(VitD_iv) <- c("chromosome","rsid","position","effect_allele",
                       "other_allele","N","eaf","beta","se","pvalue")
VitD_iv$metabolite <- "Vitamin D"
VitD_iv$units <- "SD"
VitD_iv <- VitD_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(VitD_gwas)

#################### SNPs for Microalbumin (urine) ############################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=61,675)
# 2. Values are per SD increase

### Read-in summary statistics data
Microalbumin_urine_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/Microalbumin_urine_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
Microalbumin_urine_iv <- Microalbumin_urine_gwas[Microalbumin_urine_gwas$P<5*10^(-8),] # 6 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(Microalbumin_urine_iv)[c(2,10)] <- c("rsid","pval")
Microalbumin_urine_iv <- ieugwasr::ld_clump(Microalbumin_urine_iv, clump_r2=.001, pop="EUR") # 3 SNPs retained

### Rename variables
colnames(Microalbumin_urine_iv)
Microalbumin_urine_iv <- Microalbumin_urine_iv[,c(1:10)]
colnames(Microalbumin_urine_iv) <- c("chromosome","rsid","position","effect_allele",
                                     "other_allele","N","eaf","beta","se","pvalue")
Microalbumin_urine_iv$metabolite <- "Microalbumin (urine)"
Microalbumin_urine_iv$units <- "SD"
Microalbumin_urine_iv <- Microalbumin_urine_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(Microalbumin_urine_gwas)

###################### SNPs for Potassium (urine) #############################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=197,807)
# 2. Values are per SD increase

### Read-in summary statistics data
Potassium_urine_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/Potassium_urine_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
Potassium_urine_iv <- Potassium_urine_gwas[Potassium_urine_gwas$P<5*10^(-8),] # 396 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(Potassium_urine_iv)[c(2,10)] <- c("rsid","pval")
Potassium_urine_iv <- ieugwasr::ld_clump(Potassium_urine_iv, clump_r2=.001, pop="EUR") # 3 SNPs retained

### Rename variables
colnames(Potassium_urine_iv)
Potassium_urine_iv <- Potassium_urine_iv[,c(1:10)]
colnames(Potassium_urine_iv) <- c("chromosome","rsid","position","effect_allele",
                                  "other_allele","N","eaf","beta","se","pvalue")
Potassium_urine_iv$metabolite <- "Potassium (urine)"
Potassium_urine_iv$units <- "SD"
Potassium_urine_iv <- Potassium_urine_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(Potassium_urine_gwas)

########################## SNPs for Creatinine ################################
# Note: 
# 1. GWAS was performed in a random subsample of UKB participants (~50%, n=194,465)
# 2. Values are per SD increase

### Read-in summary statistics data
Creatinine_gwas <- read.delim("Data/GWAS_summary_statistics_data/Clinical_biomarkers/Creatinine_random50.fastGWA", header=T)

### Select genome-wide significant SNPs
Creatinine_iv <- Creatinine_gwas[Creatinine_gwas$P<5*10^(-8),] # 8254 SNPs are genome-wide significant

### LD clumping, r2 <0.001 based on European population reference panel
colnames(Creatinine_iv)[c(2,10)] <- c("rsid","pval")
Creatinine_iv <- ieugwasr::ld_clump(Creatinine_iv, clump_r2=.001, pop="EUR") # 81 SNPs retained

### Rename variables
colnames(Creatinine_iv)
Creatinine_iv <- Creatinine_iv[,c(1:10)]
colnames(Creatinine_iv) <- c("chromosome","rsid","position","effect_allele",
                             "other_allele","N","eaf","beta","se","pvalue")
Creatinine_iv$metabolite <- "Creatinine"
Creatinine_iv$units <- "SD"
Creatinine_iv <- Creatinine_iv[,c("rsid","chromosome","position","effect_allele","other_allele","eaf","beta","se","pvalue","N","metabolite","units")]
rm(Creatinine_gwas)


###############################################################################
### Combine SNPs for clinical biomarkers

clinical_biomarkers_iv <- format_data(rbind(hba1c_iv, chol_iv, ldl_iv, tg_iv, crp_iv,ApoB_iv, 
                                            ALP_iv, CysC_iv, GGT_iv, Creatinine_iv, 
                                            IGF1_iv, Phosphate_iv, SHBG_iv, 
                                            Testosterone_iv, Total_bil_iv, VitD_iv, 
                                            Microalbumin_urine_iv, Potassium_urine_iv),
                                      type = "exposure",
                                      phenotype_col = "metabolite",
                                      snp_col = "rsid",
                                      beta_col = "beta",
                                      se_col = "se",
                                      eaf_col = "eaf",
                                      effect_allele_col = "effect_allele",
                                      other_allele_col = "other_allele",
                                      pval_col = "pvalue",
                                      chr_col = "chromosome",
                                      pos_col = "position",
                                      samplesize_col = "N",
                                      units_col = "units")

# Variance in exposure explained by each variant
clinical_biomarkers_iv$r2.exposure <- get_r_from_pn(clinical_biomarkers_iv$pval.exposure, clinical_biomarkers_iv$samplesize.exposure)^2

# F-statistics of each variant, F = r2 * (N-2) / (1-r2)
clinical_biomarkers_iv$f_stat.exposure <- clinical_biomarkers_iv$r2.exposure * (clinical_biomarkers_iv$samplesize.exposure) / (1-clinical_biomarkers_iv$r2.exposure)
summary(clinical_biomarkers_iv$f_stat.exposure) # all >10



#======================= IMPORT FRAILTY GWAS DATA ============================

##### GWAS summary statistics for frailty #####
# Notes:
# (1) N=403,041, excluded individuals withdrawn from UKB, non-White ancestry,
#     high heterozygosity of genotypes, sex chromosome aneuploidy, >5% missing 
#     rate of SNPs, missing FI/FP
# (2) Linear mixed models, adjusted for age, sex, array, 10 principal components

# GWAS for full sample
fi_gwas <- read.delim("Data/GWAS_summary_statistics_data/Frailty/fi.fastGWA", header=T)
fp_gwas <- read.delim("Data/GWAS_summary_statistics_data/Frailty/fp.fastGWA", header=T)

# Additional GWAS for a random subsample of UKB without metabolomics data, n=305,791
fi_nometabolomics_gwas <- read.delim("Data/GWAS_summary_statistics_data/Frailty/fi_nometabolomics.fastGWA", header=T)
fp_nometabolomics_gwas <- read.delim("Data/GWAS_summary_statistics_data/Frailty/fp_nometabolomics.fastGWA", header=T)

# Additional GWAS for random 50% subsample of UKB, n=201,574
fi_random50_gwas <- read.delim("Data/GWAS_summary_statistics_data/Frailty/fi_random50.fastGWA", header=T)
fp_random50_gwas <- read.delim("Data/GWAS_summary_statistics_data/Frailty/fp_random50.fastGWA", header=T)


##### Format outcome data for MR analysis #####
# 1. For NMR metabolomic biomarkers
mr_outcome_nometabolomics <- format_data(rbind(fi_nometabolomics_gwas %>% mutate(outcome="FI"),
                                               fp_nometabolomics_gwas %>% mutate(outcome="FP")),
                                         type = "outcome",
                                         snps = metabolomics_iv$SNP,
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

# 2. HbA1c, Total chol, LDL, HDL, Triglycerides, CRP
mr_outcome <-   format_data(rbind(fi_gwas %>% mutate(outcome="FI"),
                                  fp_gwas %>% mutate(outcome="FP")),
                            type = "outcome",
                            snps = clinical_biomarkers_iv[clinical_biomarkers_iv$exposure %in% c("CRP","Total cholesterol","HbA1c","LDL-C","Triglycerides"),]$SNP,
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

# 3. Biomarkers in which GWAS was performed in UKB: ApoB, ALP, Cystatin C, GGT, IGF-1, Phosphate, SHBG, Testosterone, Total bilirubin, VitD, Microalbumin (urine), Potassium (urine)
mr_outcome_random50 <-   format_data(rbind(fi_random50_gwas %>% mutate(outcome="FI"),
                                           fp_random50_gwas %>% mutate(outcome="FP")),
                                     type = "outcome",
                                     snps = clinical_biomarkers_iv[clinical_biomarkers_iv$exposure %in% c("ApoB","ALP","Cystatin C","Creatinine","GGT","IGF-1","Phosphate","SHBG","Testosterone","Total bilirubin","Vitamin D","Microalbumin (urine)","Potassium (urine)"),]$SNP,
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

table(!(clinical_biomarkers_iv$SNP %in% rbind(mr_outcome,mr_outcome_random50)$SNP)) # 2 SNPs not available in outcome data
table(!(metabolomics_iv$SNP %in% mr_outcome_nometabolomics$SNP)) # 102 SNPs not available in outcome data

rm(list=c("fi_nometabolomics_gwas","fp_nometabolomics_gwas","fi_gwas","fp_gwas","fi_random50_gwas","fp_random50_gwas"))


#============================ DATA HARMONIZATION =============================

##### NMR metabolomic biomarkers #####
mr_data_metabolomics <- harmonise_data(metabolomics_iv, mr_outcome_nometabolomics, action = 2)
table(mr_data_metabolomics$mr_keep) # removed 54 SNPs for being palindromic with intermediate allele frequencies

# Calculate F-statistics for SNPs, formula: ((n-k-1)/k)*(sum(r2)/(1-sum(r2)))
F_stat_metabolomics <- data.frame(metabolite=rep("",length(mr_metabolomics)), 
                                  no_snps=rep("",length(mr_metabolomics)), 
                                  var_explained=rep("",length(mr_metabolomics)),
                                  F_stat=rep("",length(mr_metabolomics)))
for (i in mr_metabolomics) {
  met <- subset(mr_data_metabolomics, mr_data_metabolomics$id.exposure==i & mr_data_metabolomics$outcome=="FI" &
                  mr_data_metabolomics$mr_keep) # SNPs that are included in the analysis
  
  # Metabolite name
  F_stat_metabolomics$metabolite[which(mr_metabolomics==i)] <- met$exposure[1]
  
  # Number of SNPs
  F_stat_metabolomics$no_snps[which(mr_metabolomics==i)] <- dim(met)[1]
  
  # Variance explained by the SNPs
  F_stat_metabolomics$var_explained[which(mr_metabolomics==i)] <- round(sum(get_r_from_pn(met$pval.exposure, met$samplesize.exposure)^2),5)
  
  # F-statistics
  F_stat_metabolomics$F_stat[which(mr_metabolomics==i)] <- 
    round( ((met$samplesize.exposure[1] - dim(met)[1] - 1) / dim(met)[1]) *
             (sum(get_r_from_pn(met$pval.exposure, met$samplesize.exposure)^2)/
                (1-sum(get_r_from_pn(met$pval.exposure, met$samplesize.exposure)^2))), 5)
  
  rm(met)
}
clipr::write_clip(F_stat_metabolomics)


##### Clinical biomarkers #####

mr_data_clinical_biomarkers <- rbind(harmonise_data(clinical_biomarkers_iv[clinical_biomarkers_iv$exposure %in% c("CRP","Total cholesterol","HbA1c","LDL-C","Triglycerides"),], mr_outcome, action = 2),
                                     harmonise_data(clinical_biomarkers_iv[clinical_biomarkers_iv$exposure %in% c("ApoB","ALP","Cystatin C","Creatinine","GGT","IGF-1","Phosphate","SHBG","Testosterone","Total bilirubin","Vitamin D","Microalbumin (urine)","Potassium (urine)"),], mr_outcome_random50, action = 2))

table(mr_data_clinical_biomarkers$mr_keep) # removed 70 SNPs for being palindromic with intermediate allele frequencies

# Calculate F-statistics for SNPs, formula: ((n-k-1)/k)*(sum(r2)/(1-sum(r2)))
F_stat_clinical_biomarkers <- data.frame(metabolite=rep("",length(mr_clinical_biomarkers)), 
                                         no_snps=rep("",length(mr_clinical_biomarkers)), 
                                         var_explained=rep("",length(mr_clinical_biomarkers)),
                                         F_stat=rep("",length(mr_clinical_biomarkers)))
for (i in mr_clinical_biomarkers) {
  met <- subset(mr_data_clinical_biomarkers, mr_data_clinical_biomarkers$exposure==i & mr_data_clinical_biomarkers$outcome=="FI" &
                  mr_data_clinical_biomarkers$mr_keep) # SNPs that are included in the analysis
  
  # Metabolite name
  F_stat_clinical_biomarkers$metabolite[which(mr_clinical_biomarkers==i)] <- met$exposure[1]
  
  # Number of SNPs
  F_stat_clinical_biomarkers$no_snps[which(mr_clinical_biomarkers==i)] <- dim(met)[1]
  
  # Variance explained by the SNPs
  F_stat_clinical_biomarkers$var_explained[which(mr_clinical_biomarkers==i)] <- round(sum(get_r_from_pn(met$pval.exposure, met$samplesize.exposure)^2),5)
  
  # F-statistics
  F_stat_clinical_biomarkers$F_stat[which(mr_clinical_biomarkers==i)] <- 
    round( ((met$samplesize.exposure[1] - dim(met)[1] - 1) / dim(met)[1]) *
             (sum(get_r_from_pn(met$pval.exposure, met$samplesize.exposure)^2)/
                (1-sum(get_r_from_pn(met$pval.exposure, met$samplesize.exposure)^2))), 5)
  
  rm(met)
}
clipr::write_clip(F_stat_clinical_biomarkers)


#===================== TWO-SAMPLE MENDELIAN RANDOMIZATION =====================

##### NMR metabolomics #####

mr_metabolomics_results <- mr(mr_data_metabolomics, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
#mr_report(mr_data_metabolomics, output_path="Output/MR_results/") # Quick generation of MR results

# Association for each SNP
mr_metabolomics_results_single <- mr_singlesnp(mr_data_metabolomics)

# Test for heterogeneity
mr_metabolomics_het <- mr_heterogeneity(mr_data_metabolomics)

# Test for pleiotropy
mr_metabolomics_plt <- mr_pleiotropy_test(mr_data_metabolomics)

# MR PRESSO
set.seed(2022)
mrpresso_metabolomics <- run_mr_presso(mr_data_metabolomics, NbDistribution = 10000, SignifThreshold = 0.05)
for (i in 1:length(mrpresso_metabolomics)) {
  print( attributes(mrpresso_metabolomics)$exposure[i] )
  print( attributes(mrpresso_metabolomics)$outcome[i] )
  print(  mrpresso_metabolomics[[i]]$`Main MR results` )
  print( mrpresso_metabolomics[[i]]$`MR-PRESSO results` )
}
rm(i)

# Summary of results
mr_summary_metabolomics <- c()
for (j in c("FI","FP")) {
  for (i in sub(".*met-d-", "", mr_metabolomics)) {
    mr_summary_metabolomics <-
      rbind(
        mr_summary_metabolomics,
        # Results for IVW, MR Egger, Weighted median, Weighted mode
        mr_metabolomics_results[mr_metabolomics_results$outcome==j & mr_metabolomics_results$exposure==i,3:9],
        # Results for MR-PRESSO
        data.frame(outcome=j,
                   exposure=i,
                   method="MR-PRESSO",
                   nsnp=mr_metabolomics_results[mr_metabolomics_results$outcome==j&mr_metabolomics_results$exposure==i,"nsnp"][1] - sum(mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`MR-PRESSO results`$`Outlier Test`$Pvalue<.05),
                   b=ifelse(!is.na(mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`Main MR results`[2,"Causal Estimate"]),
                            mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`Main MR results`[2,"Causal Estimate"],
                            mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`Main MR results`[1,"Causal Estimate"]),
                   se=ifelse(!is.na(mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`Main MR results`[2,"Sd"]),
                             mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`Main MR results`[2,"Sd"],
                             mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`Main MR results`[1,"Sd"]),
                   pval=ifelse(!is.na(mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`Main MR results`[2,"P-value"]),
                               mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`Main MR results`[2,"P-value"],
                               mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`Main MR results`[1,"P-value"])
        ))
  }
}
# P values for heterogeneity and pleiotropy
mr_p_het_plt <- c()
for (j in c("FI","FP")) {
  for (i in sub(".*met-d-", "", mr_metabolomics)) {
    mr_p_het_plt <- 
      rbind(mr_p_het_plt,
            data.frame(outcome=rep(j,5),
                       exposure=rep(i,5),
                       p_het=c(mr_metabolomics_het[mr_metabolomics_het$outcome==j&mr_metabolomics_het$exposure==i&mr_metabolomics_het$method=="Inverse variance weighted","Q_pval"],mr_metabolomics_het[mr_metabolomics_het$outcome==j&mr_metabolomics_het$exposure==i&mr_metabolomics_het$method=="MR Egger","Q_pval"],NA,NA,NA),
                       p_plt=c(NA,mr_metabolomics_plt[mr_metabolomics_plt$outcome==j&mr_metabolomics_plt$exposure==i,"pval"],NA,NA,NA))
      )
  }
}
# Outlier SNPs from MR-PRESSO
mr_presso_metabolomics_outlier_snps <- c()
for (j in c("FI","FP")) {
  for (i in sub(".*met-d-", "", mr_metabolomics)) {
    mr_presso_metabolomics_outlier_snps <- 
      rbind(mr_presso_metabolomics_outlier_snps,
            data.frame(outcome=rep(j,5),
                       exposure=rep(i,5),
                       outlier_snps=c(NA,NA,NA,NA,
                                      paste(mr_data_metabolomics[row.names(mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`MR-PRESSO results`$`Outlier Test`[mrpresso_metabolomics[[which(attributes(mrpresso_metabolomics)$exposure==i & attributes(mrpresso_metabolomics)$outcome==j)]]$`MR-PRESSO results`$`Outlier Test`$Pvalue<.05,]), "SNP"], collapse = ", ")
                       )))
  }
}

mr_summary_metabolomics <- cbind(mr_summary_metabolomics, mr_p_het_plt[,3:4], mr_presso_metabolomics_outlier_snps[,3,drop=F])
row.names(mr_summary_metabolomics) <- NULL
rm(list=c("i","j","mr_p_het_plt"))

clipr::write_clip(mr_summary_metabolomics) # Copy results
write.table(mr_summary_metabolomics, "Output/MR_results/mr_results_metabolomics.txt", row.names = F) # Export results



##### Clinical biomarkers #####

mr_clinical_biomarkers_results <- mr(mr_data_clinical_biomarkers, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
#mr_report(mr_data_clinical_biomarkers, output_path="Output/MR_results/") # Quick generation of results

# Association for each SNP
mr_clinical_biomarkers_results_single <- mr_singlesnp(mr_data_clinical_biomarkers)

# Test for heterogeneity
mr_clinical_biomarkers_het <- mr_heterogeneity(mr_data_clinical_biomarkers)

# Test for pleiotropy
mr_clinical_biomarkers_plt <- mr_pleiotropy_test(mr_data_clinical_biomarkers)

# MR PRESSO
set.seed(2022)
mrpresso_clinical_biomarkers <- run_mr_presso(mr_data_clinical_biomarkers[!mr_data_clinical_biomarkers$exposure%in%c("Microalbumin (urine)","Potassium (urine)"),], # Too few SNPs for microalbumin and potassium in urine
                                              NbDistribution = 10000, SignifThreshold = 0.05)
for (i in 1:length(mrpresso_clinical_biomarkers)) {
  print( attributes(mrpresso_clinical_biomarkers)$exposure[i] )
  print( attributes(mrpresso_clinical_biomarkers)$outcome[i] )
  print(  mrpresso_clinical_biomarkers[[i]]$`Main MR results` )
  print( mrpresso_clinical_biomarkers[[i]]$`MR-PRESSO results` )
}

# Summary of results
mr_summary_clinical_biomarkers <- c()
for (j in c("FI","FP")) {
  for (i in mr_clinical_biomarkers[!mr_clinical_biomarkers%in%c("Microalbumin (urine)","Potassium (urine)")]) {
    mr_summary_clinical_biomarkers <-
      rbind(
        mr_summary_clinical_biomarkers,
        # Results for IVW, MR Egger, Weighted median, Weighted mode
        mr_clinical_biomarkers_results[mr_clinical_biomarkers_results$outcome==j & mr_clinical_biomarkers_results$exposure==i,3:9],
        # Results for MR-PRESSO
        data.frame(outcome=j,
                   exposure=i,
                   method="MR-PRESSO",
                   nsnp=mr_clinical_biomarkers_results[mr_clinical_biomarkers_results$outcome==j&mr_clinical_biomarkers_results$exposure==i,"nsnp"][1] - sum(mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`MR-PRESSO results`$`Outlier Test`$Pvalue<.05),
                   b=ifelse(!is.na(mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`Main MR results`[2,"Causal Estimate"]),
                            mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`Main MR results`[2,"Causal Estimate"],
                            mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`Main MR results`[1,"Causal Estimate"]),
                   se=ifelse(!is.na(mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`Main MR results`[2,"Sd"]),
                             mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`Main MR results`[2,"Sd"],
                             mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`Main MR results`[1,"Sd"]),
                   pval=ifelse(!is.na(mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`Main MR results`[2,"P-value"]),
                               mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`Main MR results`[2,"P-value"],
                               mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`Main MR results`[1,"P-value"])
        ))
  }
  
  # For microalbumin and potassium in urine, there were not enough SNPs for MR-PRESSO analysis
  for (i in c("Microalbumin (urine)","Potassium (urine)")) {
    mr_summary_clinical_biomarkers <- rbind(
      mr_summary_clinical_biomarkers,
      # Results for IVW, MR Egger, Weighted median, Weighted mode
      mr_clinical_biomarkers_results[mr_clinical_biomarkers_results$outcome==j & mr_clinical_biomarkers_results$exposure==i,3:9],
      # Results for MR-PRESSO
      data.frame(outcome=j,
                 exposure=i,
                 method="MR-PRESSO",
                 nsnp=NA,
                 b=NA,
                 se=NA,
                 pval=NA)
    )
  }
}
# P values for heterogeneity and pleiotropy
mr_p_het_plt <- c()
for (j in c("FI","FP")) {
  for (i in mr_clinical_biomarkers) {
    mr_p_het_plt <- 
      rbind(mr_p_het_plt,
            data.frame(outcome=rep(j,5),
                       exposure=rep(i,5),
                       p_het=c(mr_clinical_biomarkers_het[mr_clinical_biomarkers_het$outcome==j&mr_clinical_biomarkers_het$exposure==i&mr_clinical_biomarkers_het$method=="Inverse variance weighted","Q_pval"],mr_clinical_biomarkers_het[mr_clinical_biomarkers_het$outcome==j&mr_clinical_biomarkers_het$exposure==i&mr_clinical_biomarkers_het$method=="MR Egger","Q_pval"],NA,NA,NA),
                       p_plt=c(NA,mr_clinical_biomarkers_plt[mr_clinical_biomarkers_plt$outcome==j&mr_clinical_biomarkers_plt$exposure==i,"pval"],NA,NA,NA))
      )
  }
}
# Outlier SNPs from MR-PRESSO
mr_presso_clinical_biomarkers_outlier_snps <- c()
for (j in c("FI","FP")) {
  for (i in mr_clinical_biomarkers[!mr_clinical_biomarkers%in%c("Microalbumin (urine)","Potassium (urine)")]){
    mr_presso_clinical_biomarkers_outlier_snps <- 
      rbind(mr_presso_clinical_biomarkers_outlier_snps,
            data.frame(outcome=rep(j,5),
                       exposure=rep(i,5),
                       outlier_snps=c(NA,NA,NA,NA,
                                      paste(mr_data_clinical_biomarkers[row.names(mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`MR-PRESSO results`$`Outlier Test`[mrpresso_clinical_biomarkers[[which(attributes(mrpresso_clinical_biomarkers)$exposure==i & attributes(mrpresso_clinical_biomarkers)$outcome==j)]]$`MR-PRESSO results`$`Outlier Test`$Pvalue<.05,]), "SNP"], collapse = ", ")
                       )))
  }
  for (i in mr_clinical_biomarkers[mr_clinical_biomarkers%in%c("Microalbumin (urine)","Potassium (urine)")]){
    mr_presso_clinical_biomarkers_outlier_snps <- 
    rbind(mr_presso_clinical_biomarkers_outlier_snps,
          data.frame(outcome=rep(j,5),
                     exposure=rep(i,5),
                     outlier_snps=c(NA,NA,NA,NA,NA)
                     ))
  }
}

mr_summary_clinical_biomarkers <- cbind(mr_summary_clinical_biomarkers, mr_p_het_plt[,3:4], mr_presso_clinical_biomarkers_outlier_snps[,3,drop=F])
row.names(mr_summary_clinical_biomarkers) <- NULL
rm(list=c("i","j","mr_p_het_plt"))

clipr::write_clip(mr_summary_clinical_biomarkers) # Copy results
write.table(mr_summary_clinical_biomarkers, "Output/MR_results/mr_results_clinical_biomarkers.txt", row.names = F) # Export results


#========================== SUMMARY OF MR RESULTS ============================

### Combine results from metabolomics and clinical biomarkers
mr_results <- rbind(mr_metabolomics_results,mr_clinical_biomarkers_results)
mr_single <- rbind(mr_metabolomics_results_single, mr_clinical_biomarkers_results_single)
mr_het <- rbind(mr_metabolomics_het, mr_clinical_biomarkers_het)
mr_plt <- rbind(mr_metabolomics_plt, mr_clinical_biomarkers_plt)
mr_results_summary <- rbind(mr_summary_metabolomics,mr_summary_clinical_biomarkers)

# FDR correction for IVW
p_table <- subset_on_method(mr_results)$pval
p_table <- p_table[ order(p_table) ]
BH_table <- data.frame(p_table=p_table,BH_table=.05*1:length(p_table)/length(p_table));BH_table
fdr_p_mr <- round(max(BH_table$BH_table[BH_table$BH_table>BH_table$p_table]),3);fdr_p_mr # p=0.011
rm(list=c("p_table","BH_table"))

# SNPs data
clipr::write_clip(
  rbind(mr_data_metabolomics[mr_data_metabolomics$outcome=="FI"&mr_data_metabolomics$mr_keep,c("exposure","outcome","SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","eaf.outcome","beta.outcome","se.outcome","pval.outcome","palindromic","ambiguous","mr_keep")],
        mr_data_metabolomics[mr_data_metabolomics$outcome=="FP"&mr_data_metabolomics$mr_keep,c("exposure","outcome","SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","eaf.outcome","beta.outcome","se.outcome","pval.outcome","palindromic","ambiguous","mr_keep")])
)
clipr::write_clip(
  rbind(mr_data_clinical_biomarkers[mr_data_clinical_biomarkers$outcome=="FI"&mr_data_clinical_biomarkers$mr_keep,c("exposure","outcome","SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","eaf.outcome","beta.outcome","se.outcome","pval.outcome","palindromic","ambiguous","mr_keep")],
        mr_data_clinical_biomarkers[mr_data_clinical_biomarkers$outcome=="FP"&mr_data_clinical_biomarkers$mr_keep,c("exposure","outcome","SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","eaf.outcome","beta.outcome","se.outcome","pval.outcome","palindromic","ambiguous","mr_keep")])
)


### Scatter plots
mr_scatter_plots_metabolomics <- mr_scatter_plot(mr_metabolomics_results, mr_data_metabolomics)
for (i in 1:length(mr_scatter_plots_metabolomics)) {
  exposure <- mr_scatter_plots_metabolomics[[i]]$data$exposure[1]
  outcome <- mr_scatter_plots_metabolomics[[i]]$data$outcome[1]
  filename <- paste0("Output/MR_results/Scatter_plot/MR_metabolomics_",exposure,"_",outcome,"_scatterplot.jpg")
  jpeg(file=filename,  width = 12, height = 10, units = "cm", res = 300, quality = 100)
  print(mr_scatter_plots_metabolomics[[i]] + 
          scale_color_brewer(palette="Dark2") +
          guides(color=guide_legend(nrow=1,byrow=T)) + 
          theme_bw() + 
          theme(legend.position="top",
                legend.text = element_text(size=7),
                legend.title = element_blank(),
                legend.key.size = unit(3, "mm")))
  dev.off()
  rm(list=c("i","exposure","outcome","filename"))
}
rm(mr_scatter_plots_metabolomics)

mr_scatter_plots_clinical_biomarkers <- mr_scatter_plot(mr_clinical_biomarkers_results, mr_data_clinical_biomarkers)
for (i in 1:length(mr_scatter_plots_clinical_biomarkers)) {
  exposure <- mr_scatter_plots_clinical_biomarkers[[i]]$data$exposure[1]
  outcome <- mr_scatter_plots_clinical_biomarkers[[i]]$data$outcome[1]
  filename <- paste0("Output/MR_results/Scatter_plot/MR_clinical_biomarkers_",exposure,"_",outcome,"_scatterplot.jpg")
  jpeg(file=filename,  width = 12, height = 10, units = "cm", res = 300, quality = 100)
  print(mr_scatter_plots_clinical_biomarkers[[i]] + 
          scale_color_brewer(palette="Dark2") +
          guides(color=guide_legend(nrow=1,byrow=T)) + 
          theme_bw() + 
          theme(legend.position="top",
                legend.text = element_text(size=7),
                legend.title = element_blank(),
                legend.key.size = unit(3, "mm")))
  dev.off()
  rm(list=c("i","exposure","outcome","filename"))
}
rm(mr_scatter_plots_clinical_biomarkers)



#==================== SAVE DATA FOR FURTHER ANALYSIS =========================

save.image("Data/R_data/Metabolites_frailty_MR.Rdata")

# =============================== END OF FILE  ===============================