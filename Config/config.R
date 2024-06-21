#************** PARAMETERS ****************#

#Variable to use real or scrambeled (not real) data. Set to TRUE to use scrambled data.
scramble_data <- FALSE

#************** NOTHING BELOW NEEDS TO BE CHANGED ****************#

#Function load libraries or install where needed
install_if_needed <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name)
    library(package_name, character.only = TRUE)
  }
}

#Create vector with packages
packages <- c("tidyverse", "ggplot2", "readxl", "broom", "data.table", "rentrez", "writexl", "snpStats", "LDheatmap", "lubridate", "rstudioapi")

#Nest function in loop and load libraries
for (pkg in packages) {
  install_if_needed(pkg)
}

#This will dynamically set the project_dir to the location of the script accessing this CONFIG script.
project_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

#FILE AND SCRIPT ID MAPPING

SNP_Selection <- 1
Genotype_Data_Processing <- 2
Data_Integration <- 3
Dataset_Cleaning <- 4
SNP_Stats_Filtering <- 5

#GLOBAL FUNCTIONS AND VARIABLES

#This function dynamically locates a scripts input directory based on the number in the name of the folder - which must macth the number in the name of the script.
set_dir <- function(script_id) {
  #Determine target folder based on scramble_data variable
  if (scramble_data) {
    target_folder <- "Scrambled Data"
  } else {
    target_folder <- "Data Files"
  }
  
  # Path to Data Files directory
  data_files_dir <- file.path(project_dir, as.character(target_folder))
  
  # List directories in data_files_dir
  dirs <- list.dirs(data_files_dir, full.names = TRUE, recursive = FALSE)
  
  # Find the directory that starts with the same number as script_id
  dir_pattern <- paste0("^", script_id)
  target_dir <- dirs[grepl(dir_pattern, basename(dirs))]
  
  # If a matching directory is found, set Input_Dir
  if (length(target_dir) > 0) {
    input_dir <- target_dir
  } else {
    input_dir <- NULL
  }
  
  return(input_dir)
}

#Define a function to encode genotypes numerically form "AA" or "GC" type string encoding.
encode_genotype <- function(genotype_string) {
  # Replace missing data indicated by "--" with NA
  genotype_string <- gsub("--", NA, genotype_string)
  
  # Identify the unique alleles (ignoring NA)
  alleles <- unique(na.omit(unlist(strsplit(genotype_string, ""))))
  
  # Proceed with mapping, accounting for the possibility of only one allele being available
  return(sapply(genotype_string, function(gt) {
    if (is.na(gt)) {
      return(NA)  # Missing data
    } else if (length(alleles) == 1 && gt == paste0(alleles[1], alleles[1])) {
      return(0)  # Homozygous for the only allele
    } else if (length(alleles) == 2) {
      if (gt == paste0(alleles[1], alleles[1])) {
        return(0)  # Homozygous for the first allele
      } else if (gt == paste0(alleles[2], alleles[2])) {
        return(2)  # Homozygous for the second allele
      } else if (gt == paste0(alleles[1], alleles[2]) || gt == paste0(alleles[2], alleles[1])) {
        return(1)  # Heterozygous
      }
    }
    return(NA)  # Handle non-standard data or if alleles are not exactly 1 or 2
  }))
}

#Identifying snp columns as those beginning with rs and followed by a number, provided DF
get_snp_columns <- function(df) {
  snp_columns <- grep("^rs[0-9]+", names(df), value = TRUE)
  return(snp_columns)
}
