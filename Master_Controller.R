#** ------ MASTER CONTROLLER ------ **#
# =================================== #
#    Zach Ribau  |  July 12, 2024    

rm(list=ls())

#**____________________________________________**#
## -------------  PARAMETERS ----------------- ##

#Set to True to run snp selection script, and false to skip (will use pre-existing filtered genotype file)
run_snp_select = FALSE

#Specify columns in each phenotype category to incorporate into model formulas
diet_response_columns <- c("kcal", "carb_prop", "prot_prop", "tfat_prop") 

sleep_response_columns <- c("Mean_TST", "Mean_WASO", "Regularity")

covariates = c("age", "sex", "bmi")

#Correlation structure in general estimating equations
correlation_structure = "exchangeable"

minor_allele_threshold <- 0.05

#Assign method for adjusting pvalues for multiple testing.
pval_adj_method = "fdr"

p_value_threshold <- 0.05

adj_OR_raw <- "adj" #set to "ajd" to consider adjusted pvalues for significant or "raw" to consider raw pvalues as sig.


#**_________________________________________________**#
## ------------- ANALYSIS PIPELINE ----------------- ##

#Run CONFIG script
project_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))


#Run the analysis (scripts 1-5)
run_all_scripts(project_dir)

#For access, datasets are saved as variables here. 
full_cohort <- Database[["Complete Datasets"]][["Full Cohort"]]
Adults <- Database[["Complete Datasets"]][["Adults"]]
Children <- Database[["Complete Datasets"]][["Childrens"]]





