#** -------- 7. Results Investigation -------- **#
# ============================================== #
#        Zach Ribau  |  Sep 26, 2024       

library(tidyverse)

project_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(project_dir)

#First we can investigate previous associations for each the unique SNPs identified in our analysis
sig_snps_associations <- sig_snps

sig_snps <- sig_snps_grouped$term

#SNP Data aggregation
sig_snps_grouped <- sig_snps_grouped %>%
  rename(SNP = "term")

adult_snp_gen_counts <- Database[["Genetic Stats"]][["Adults"]][["SNP_Stats"]][["Allelic Stats"]] %>%
  rename(SNP = "rsID")

children_snp_gen_counts <- Database[["Genetic Stats"]][["Children"]][["SNP_Stats"]][["Allelic Stats"]] %>%
  rename(SNP = "rsID")

full_cohort_snp_gen_counts <- full_join(adult_snp_gen_counts, children_snp_gen_counts, by = "SNP", suffix = c("_adult", "_child")) %>%
  # Sum the corresponding genotype counts (handling potential NAs)
  mutate(
    hzREF_total = coalesce(hzREF_adult, 0) + coalesce(hzREF_child, 0),
    htz_total = coalesce(htz_adult, 0) + coalesce(htz_child, 0),
    hzALT_total = coalesce(hzALT_adult, 0) + coalesce(hzALT_child, 0)
  ) %>%
  # Optionally, select only the relevant columns
  select(SNP, hzREF_total, htz_total, hzALT_total)

sig_snps_associations <- sig_snps_associations %>%
  rename(SNP = "term")

sig_snps_associations_adults <- sig_snps_associations %>%
  filter(cohort == "adults")

sig_snps_associations_children <- sig_snps_associations %>%
  filter(cohort == "children")

sig_snps_associations_full <- sig_snps_associations %>%
  filter(cohort == "full")

sig_snps_associations_children <- left_join(sig_snps_associations_children, children_snp_gen_counts, by = "SNP")
sig_snps_associations_adults <- left_join(sig_snps_associations_adults, adult_snp_gen_counts, by = "SNP")
sig_snps_associations_full <- left_join(sig_snps_associations_full, full_cohort_snp_gen_counts, by = "SNP")

#access the EBI database again
ebi <- read_tsv("Data Files/1. SNP Selection/EBI.tsv")

#Load in trait terms with trait group assignments table
sleep_traits <- read.csv("Data Files/1. SNP Selection/sleep_traits.csv")
diet_traits <- read.csv("Data Files/1. SNP Selection/diet_traits.csv")

ebi_sleep <- ebi %>%
  filter(`MAPPED_TRAIT` %in% sleep_traits$Trait)

#filter ebi_sleep to only contain SNPs that were significant in our analysis
ebi_sleep <- ebi_sleep %>%
  filter(SNPS %in% sig_snps) %>%
  group_by(SNPS) %>%
  summarise(Associations = paste(MAPPED_TRAIT, collapse = " | "), Studies = paste(PUBMEDID, collapse = " | "))#show associations for each SNP

ebi_sleep <- ebi_sleep %>%
  rename(Previous_Associations = "Associations")

ebi_sleep <- ebi_sleep %>%
  rename(Previous_Studies_PUBMEDID = "Studies")

#repeat the same process for diet associations in ebi
ebi_diet <- ebi %>%
  filter(`MAPPED_TRAIT` %in% diet_traits$Trait)

ebi_diet <- ebi_diet %>%
  filter(SNPS %in% sig_snps) %>%
  group_by(SNPS) %>%
  summarise(Associations = paste(MAPPED_TRAIT, collapse = " | "), Studies = paste(PUBMEDID, collapse = " | "))#show associations for each SNP

ebi_diet <- ebi_diet %>%
  rename(Previous_Associations = "Associations")

ebi_diet <- ebi_diet %>%
  rename(Previous_Studies_PUBMEDID = "Studies")

#And out of interest, we will filter the remaining associations in EBI not already in our filtered diet and sleep dataframes, and repeat there too

ebi_other <- ebi %>%
  filter(!(`MAPPED_TRAIT` %in% diet_traits$Trait)) %>%
  filter(!(`MAPPED_TRAIT` %in% sleep_traits$Trait))

ebi_other <- ebi_other %>%
  filter(SNPS %in% sig_snps) %>%
  group_by(SNPS) %>%
  summarise(Associations = paste(MAPPED_TRAIT, collapse = " | "), Studies = paste(PUBMEDID, collapse = " | "))#show associations for each SNP

ebi_other <- ebi_other %>%
  rename(Previous_Associations = "Associations")

ebi_other <- ebi_other %>%
  rename(Previous_Studies_PUBMEDID = "Studies")




