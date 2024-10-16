#** -------- 7. Results Investigation -------- **#
# ============================================== #
#        Zach Ribau  |  Sep 26, 2024       

library(tidyverse)

project_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(project_dir)
output_dir <- "figures"

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


adult_barplotDF <- Database[["Complete Datasets"]][["Adults"]] %>%
  select(pid,bmi, age, sex, all_of(
    intersect(
      names(Database[["Complete Datasets"]][["Adults"]]),
      c(sig_snps_associations %>% filter(cohort == "adults") %>% pull(trait) %>% unique(),
        sig_snps_associations %>% filter(cohort == "adults") %>% pull(SNP) %>% unique() )
    )
  ))

children_barplotDF <- Database[["Complete Datasets"]][["Children"]] %>%
  select(pid, bmi, age, sex, all_of(
    intersect(
      names(Database[["Complete Datasets"]][["Children"]]),
      c(sig_snps_associations %>% filter(cohort == "children") %>% pull(trait) %>% unique(),
        sig_snps_associations %>% filter(cohort == "children") %>% pull(SNP) %>% unique() )
    )
  ))

full_barplotDF <- Database[["Complete Datasets"]][["Full Cohort"]] %>%
  select(pid, age, sex, all_of(
    intersect(
      names(Database[["Complete Datasets"]][["Full Cohort"]]),
      c(sig_snps_associations %>% filter(cohort == "full") %>% pull(trait) %>% unique(),
        sig_snps_associations %>% filter(cohort == "full") %>% pull(SNP) %>% unique() )
    )
  ))

library(ggpubr)

create_anova_barplot <- function(data, value_col, snp_info, output_dir = "Figures", plot_suffix = "") {
  
  value_col <- rlang::sym(value_col)
  snp_col <- rlang::sym(snp_info$SNP)
  
  # Ensure the SNP column is a factor
  data <- data %>%
    mutate(!!snp_col := as.factor(!!snp_col))
  
  # Set genotype labels dynamically based on SNP info
  ref_allele <- snp_info$Ref_Allele
  alt_allele <- snp_info$Alt_Allele
  custom_labels <- c("0" = paste0(ref_allele, ref_allele),
                     "1" = paste0(ref_allele, alt_allele),
                     "2" = paste0(alt_allele, alt_allele))
  
  custom_colors <- c("0" = "#EBDDE6", "1" = "#B599B3", "2" = "#885075")
  
  # Create the initial ggplot bar plot using the original dataset
  plot <- ggplot(data, aes(x = !!snp_col, y = !!value_col)) +
    # Plot mean values as bars
    stat_summary(fun = "mean", geom = "col", aes(fill = !!snp_col), color = "black", width = 0.85) +
    # Plot error bars for mean +/- confidence interval
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.1, color = "#22292F") +
    # Add mean labels above the bars
    stat_summary(fun = "mean", geom = "text", aes(label = round(..y.., 2)), vjust = -0.5, size = 5, color = "#313745") +
    # Set custom labels and colors for the bars
    scale_x_discrete(labels = custom_labels) +
    scale_fill_manual(values = custom_colors, labels = custom_labels) + # Set labels for the legend
    theme_minimal() +
    theme(
      plot.title = element_text(size = 22, color = "#313745", face = "bold"),
      axis.text = element_text(size = 15, color = "#313745"),
      axis.title = element_text(size = 15, color = "#313745"),
      legend.title = element_blank(), # Hide redundant legend title
      legend.text = element_text(size = 14)
    ) +
    labs(
      x = "Genotype",
      y = paste("Mean", rlang::as_string(value_col)),
      title = paste("Mean", rlang::as_string(value_col), "by", snp_info$SNP, "Genotype")
    )
  
  # Save the plot
  output_file <- file.path(output_dir, paste0(snp_info$SNP, "_", rlang::as_string(value_col), "_plot", plot_suffix, ".svg"))
  ggsave(output_file, plot, width = 12, height = 8.2, bg = "transparent")
  
  return(plot)
}

snp_info <- sig_snps_associations_adults %>% filter(SNP == "rs1801260")

plot <- create_anova_barplot(
  data = adult_barplotDF,
  value_col = snp_info$trait,
  snp_info = snp_info,
  output_dir = output_dir,
  plot_suffix = paste0("_", snp_info$SNP, "_", snp_info$trait)
)

### ANOVA For Genotype Means ####

#Formula to compute stats and perform ANOVA for given significant SNP
calculate_stats <- function(data, value_col, group_col) {
  
  value_col <- rlang::sym(value_col)
  group_col <- rlang::sym(group_col)
  
  stats <- data %>%
    group_by(!!group_col) %>%
    summarise(
      mean = mean(!!value_col, na.rm = TRUE),
      se = sd(!!value_col, na.rm = TRUE) / sqrt(n()),
      ci_lower = mean - qt(1 - (0.05 / 2), n() - 1) * se,
      ci_upper = mean + qt(1 - (0.05 / 2), n() - 1) * se,
      .groups = 'drop'
    )
  
  # Perform ANOVA
  anova_formula <- as.formula(paste(value_col, "~", group_col))
  anova_results <- summary(aov(anova_formula, data = data))
  
  return(list(stats = stats, anova = anova_results))
}


sig_snps_associations <- sig_snps_associations %>%
  mutate(barplotDF = case_when(
    cohort == "full"  ~ "full_barplotDF",
    cohort == "adults" ~ "adult_barplotDF",
    cohort == "children" ~ "children_barplotDF"
  ))

#Now apply function for all given SNPs, 

#For FULL COHORT DATA

rs2948276_regularity <- calculate_stats(full_barplotDF, "Regularity", "rs2948276")
stats_data <- rs2948276_regularity$stats
anova_results <- rs2948276_regularity$anova
ancova_result <- aov(Regularity ~ factor(rs2948276) + age + factor(sex), data = full_barplotDF)
tukey_result <- TukeyHSD(ancova_result, "factor(rs2948276)")
print(tukey_result)



rs1154155_tfat_prop <- calculate_stats(full_barplotDF, "tfat_prop", "rs1154155")
stats_data <- rs1154155_tfat_prop$stats
anova_results <- rs1154155_tfat_prop$anova
ancova_result <- aov(tfat_prop ~ factor(rs1154155) + age + factor(sex), data = full_barplotDF)
tukey_result <- TukeyHSD(ancova_result, "factor(rs1154155)")
print(tukey_result)

rs1823125_kcal <- calculate_stats(full_barplotDF, "kcal", "rs1823125")
stats_data <- rs1823125_kcal$stats
anova_results <- rs1823125_kcal$anova
ancova_result <- aov(kcal ~ factor(rs1823125) + age + factor(sex), data = full_barplotDF)
tukey_result <- TukeyHSD(ancova_result, "factor(rs1823125)")
print(tukey_result)


#Adult Data



