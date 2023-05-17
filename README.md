# Metabolomics_frailty

This document provides explanations on the codes used in the paper "Unraveling the metabolic underpinnings of frailty using multi-cohort observational and Mendelian randomization analyses".

Jonathan K. L. Mak, Department of Medical Epidemiology and Biostatistics, Karolinska Institutet, Sweden

Email: jonathan.mak@ki.se 


## Citation

Mak, J.K.L., Kananen, L., Qin, C. et al. Unraveling the metabolic underpinnings of frailty using multicohort observational and Mendelian randomization analyses. Aging cell (2023). https://doi.org/10.1111/acel.13868


## Abstract

Identifying metabolic biomarkers of frailty, an age-related state of physiological decline, is important for understanding its metabolic underpinnings and developing preventive strategies. Here, we systematically examined 168 nuclear magnetic resonance-based metabolomic biomarkers and 32 clinical biomarkers for their associations with frailty. In up to 90,573 UK Biobank participants, we identified 59 biomarkers robustly and independently associated with the frailty index (FI). Of these, 34 associations were replicated in the Swedish TwinGene study (n=11,025) and the Finnish Health 2000 Survey (n=6,073). Using two-sample Mendelian randomization, we showed that the genetically predicted level of glycoprotein acetyls, an inflammatory marker, was statistically significantly associated with an increased FI (β per SD increase=0.37%, 95% confidence interval: 0.12–0.61). Creatinine and several lipoprotein lipids were also associated with increased FI, yet their effects were mostly driven by kidney and cardiometabolic diseases, respectively. Our findings provide new insights into the causal effects of metabolites on frailty and highlight the role of chronic inflammation underlying frailty development.


## Codes and documents

#### Program/1.1_UKB_biomarkers.R

Codes used for assessing the associations between metabolic biomarkers (including 168 NMR metabolomic and 32 clinical biomarkers) and frailty (including a frailty index and frailty phenotype) in the discovery cohorts of UK Biobank.

#### Program/1.2_TwinGene_biomarkers.R

Codes used for assessing the associations between metabolic biomarkers and frailty index in the validation cohort of TwinGene.

#### Program/1.3_meta_analysis.R

Codes used for meta-analysis of the observational associations in UK Biobank, TwinGene, and Health 2000 Survey.

#### Program/2.1_MR_frailty_selected_metabolites.R

Codes used for performing the main Mendelian randomization analysis for the selected metabolites and frailty.

#### Program/2.2_MR_frailty_sensitivity.R

Codes used for performing the sensitivity MR analyses.

#### Program/3_output_figures.R

Codes used for generating the main figures in the manuscript.

#### Biomarker_list.xlsx

Excel containing the list of biomarkers used in UK Biobank and TwinGene
