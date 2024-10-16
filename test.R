
project_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

library(tidyverse)

full_cohort <- Database[["Complete Datasets"]][["Full Cohort"]]
Adults <- Database[["Complete Datasets"]][["Adults"]]
Children <- Database[["Complete Datasets"]][["Childrens"]]



Adults_regul_gg <- Adults %>%
  filter(rs2948276 == 2)


Adults_regul_g <- Adults %>%
  filter(rs2948276 == 1)

Adults_regul <- Adults %>%
  filter(rs2948276 == 0)