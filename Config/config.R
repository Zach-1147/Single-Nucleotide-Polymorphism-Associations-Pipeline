
#**________________________________________________________**#
## --------- 2. LOCKED FUNCTIONS & VARIABLES -------------- ##

run_all_scripts <- function(project_dir) {
  
  #Conditionally run snp selection based on global variable
  if (run_snp_select) {
    cat("Running SNP Selection\n")
    source(file.path(project_dir, "1-SNP Selection.R")) 
  } else {
    cat("Skipping SNP selection.\n")
  }

  cat("\n")
  cat("Running Data Integration...\n")
  source(file.path(project_dir, "2-Data Integration.R"))
  cat("\n")
  cat("Running Dataset Cleaning...\n")
  source(file.path(project_dir, "3-Dataset Cleaning.R"))
  cat("\n")
  cat("Running SNP Stats & Filtering...\n")
  cat("\n")
  source(file.path(project_dir, "4-SNP Stats & Filtering.R"))
  cat("\n\nRunning Association Testing...\n\n\n")
  source(file.path(project_dir, "5-Association Testing.R"))
  
  cat("\n")
  
  print(paste(as.character(uniquesigSNPs), "unique significant SNPs"))
  
  complete <- Database[["Complete Datasets"]][["Full Cohort"]]
  
  adult_size <- complete %>%
    filter(cohort == "Adults")
  
  children_size <- complete %>%
    filter(cohort == "Children")
  
  adult_sample <- nrow(adult_size)
  children_sample <- nrow(children_size)
  
  cat("\n\n")
  
  print(paste(as.character(adult_sample), "adults in final analysis"))
  print(paste(as.character(children_sample), "children in final analysis"))
  
  #Now clear env other than specific things we are interested in retaining
  rm(list = (setdiff(ls(.GlobalEnv), c("Database", "sig_snps", "all_tests", "sig_snps_grouped", "snp_columns", "diet_response_columns", "sleep_response_columns", "covariates"))), envir = .GlobalEnv)
  
}

#************** VARIABLES ****************#
#This will dynamically set the project_dir to the location of the script accessing this CONFIG script.
project_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

#vector with packages
packages <- c("tidyverse", "ggplot2", "readxl", "broom", "data.table", "rentrez", "writexl", "snpStats", "LDheatmap", "lubridate", "rstudioapi", "data.table", "geepack", "car")

#FILE AND SCRIPT ID MAPPING

SNP_Selection <- 1
Data_Integration <- 2
Cleaning_Response_Variables <- 3
Criteria_Based_Filtering <- 4
SNP_Stats <- 5
LD_MAF <- 6
Association_Testing <- 7
Figures <- 8

#****************** FUNCTIONS ********************#

#Function load libraries or install where needed
install_if_needed <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name)
    library(package_name, character.only = TRUE)
  }
}

#Install libraries in packages vector when this config script is sourced
for (pkg in packages) {
  install_if_needed(pkg)
}

#This function dynamically locates a scripts input directory based on the number in the name of the folder - which must macth the number in the name of the script.
set_dir <- function(script_id) {
  
  #Assign target folder based on dir name in folder for data files
  target_folder <- "Data Files"

  #set the dir to data files
  data_files_dir <- file.path(project_dir, as.character(target_folder))
  
  #list all directories in the datafiles folder
  dirs <- list.dirs(data_files_dir, full.names = TRUE, recursive = FALSE)
  
  #look for patterns matching the script_id passed to the function as arg to find the right input dir
  dir_pattern <- paste0("^", script_id)
  #Set match as target dir
  target_dir <- dirs[grepl(dir_pattern, basename(dirs))]
  
  #rename return variable to input dir and return
  input_dir <- target_dir

  return(input_dir)
}


#Function to identify snp columns as those beginning with rs and followed by a number, provided a DF. Will be useful for filtering df's and subsetting etc...
get_snp_columns <- function(df) {
  snp_columns <- grep("^rs[0-9]+", names(df), value = TRUE)
  return(snp_columns)
}

#function to encode genotypes numerically form "AA" or "GC" type character encoding. This function will applied over all participant's for a given SNP (ie. a column) in a genotype dataframe, so we pass it a vector of strings essentially that are formatted as "AT" or "AA".
encode_genotype <- function(genotype_string, rsID, snp_stats) {
  
  # Replace missing data with NA
  genotype_string <- gsub("--", NA, genotype_string)
  
  # Lookup the ref and alt alleles for the given SNP from snp_stats
  snp_info <- snp_stats[snp_stats$rsID == rsID, ]
  
  # If the SNP isn't found in snp_stats or has invalid Ref/Alt alleles, return NA for all genotypes
  if (nrow(snp_info) == 0 || is.na(snp_info$Ref_Allele) || is.na(snp_info$Alt_Allele)) {
    print(paste("Missing or invalid allele information for rsID:", rsID))
    return(rep(NA, length(genotype_string)))
  }
  
  # Extract the reference and alternate alleles, ensuring they're uppercase
  ref_allele <- toupper(snp_info$Ref_Allele)
  alt_allele <- toupper(snp_info$Alt_Allele)
  
  
  # Apply checks over each observation to encode numerically based on alleles
  return(sapply(genotype_string, function(gt) {
    
    # For NA's
    if (is.na(gt)) {
      return(NA) # Return NA if there is nothing
    } else {
      
      # Split the genotype string into individual alleles
      alleles <- toupper(strsplit(gt, "")[[1]])
      
      
      # Check for homozygous reference (encode as 0)
      if (alleles[1] == ref_allele && alleles[2] == ref_allele) {
        return(0)  # Homozygous for reference allele
        
        # Check for homozygous alternate (encode as 2)
      } else if (alleles[1] == alt_allele && alleles[2] == alt_allele) {
        return(2)  # Homozygous for alternate allele (minor allele)
        
        # Check for heterozygous (encode as 1) regardless of order
      } else if ((alleles[1] == ref_allele && alleles[2] == alt_allele) || 
                 (alleles[1] == alt_allele && alleles[2] == ref_allele)) {
        return(1)  # Heterozygous (regardless of allele order)
        
      } else {
      }
    }
    
    return(NA)  # Return NA for unexpected cases
  }))
}



