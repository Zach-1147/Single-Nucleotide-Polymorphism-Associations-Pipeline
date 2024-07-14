
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

#This is filtered out because it takes quite some time to run. The end file is already in the input dir for the next script, so we can skip it.
#source(file.path(project_dir, "1-SNP Selection.R")) 

source(file.path(project_dir, "2-Data Integration.R"))

source(file.path(project_dir, "3-Dataset Cleaning.R"))

#**_____________________________________________________**#
## --------------- 2. QUALITY CONTROL ------------------ ##

source(file.path(project_dir, "4-SNP Stats & Filtering.R"))

#**_____________________________________________________**#
## ---------------- 3. ASSOCIATIONS  ------------------- ##

source(file.path(project_dir, "5-Association Testing.R"))


#Final adult sample size

complete <- Database[["Complete Datasets"]][["Full Cohort"]]

complete <- complete %>%
  group_by(fid) %>%
  summarise(
    Count = n(),
    Ad = sum(case_when(cohort == "Adults" ~ 1, TRUE ~ 0)),
    Ch = sum(case_when(cohort == "Children" ~ 1, TRUE ~ 0)),
    Yes = case_when(Ch > Ad ~ TRUE, TRUE ~ FALSE)
    
    )

complete <- complete%>%
  filter(Yes == TRUE)
