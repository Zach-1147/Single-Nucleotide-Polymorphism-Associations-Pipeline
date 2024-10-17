#** ------ CRITERIA BASED FILTERING ------- **#
# =============================================== #
#    Zach Ribau  |  Oct 17, 2024    
rm(list=ls())
#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

Input_Dir <- set_dir(Criteria_Based_Filtering)
Output_Dir <- set_dir(SNP_Stats)

setwd(Input_Dir)

#**_____________________________________________________**#
## --------------------- 1. DATA ------------------------ ##

#Lets start by loading in the diet/demo and genotype df's
adult <- read.csv("final_adult.csv")
children <- read.csv("final_child.csv")

##**_____________________________________________________**##
## ------------- 2. SLEEP RESPONSE VARIABLES ------------- ##

## DERIVING MEAN SLEEP AND WAKE TIMES ##


# WORK HERE..--....---...----
# Look at distributions, demo's, etc...




## SAVING DATA ##

setwd(Output_Dir)

#write df's to csv

write.csv(adult, "adult_processed.csv", row.names = FALSE)
write.csv(children, "children_processed.csv", row.names = FALSE)