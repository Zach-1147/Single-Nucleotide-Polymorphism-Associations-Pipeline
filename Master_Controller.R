
#** ------ MASTER CONTROLLER ------ **#
# =================================== #
#    Zach Ribau  |  June 20, 2024    

#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

rm(list=ls()) #Clear environment
project_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

#Load configuration script
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

#**_____________________________________________________**#
## ---------------- 1. PREPROCESSING ------------------- ##

source(file.path(project_dir, "3-Data Integration.R"))

source(file.path(project_dir, "4-Dataset Cleaning.R"))

#**_____________________________________________________**#
## --------------- 2. QUALITY CONTROL ------------------ ##

source(file.path(project_dir, "5-SNP Stats & Filtering.R"))


#**_____________________________________________________**#
## --------------- 3.   ------------------ ##


