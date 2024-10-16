#** -------- 6. Figures -------- **#
# ========================================= #
#    Zach Ribau  |  July 24, 2024       

#Note that this script is not dynamic and will not generate figures for SNPs other than those hardocded in here. Thus, modifying parameters in the master controller and running all scripts wll throw an error once it gets to this script, should the sig snp set change. This is why it is not appart of the run_all_scripts function from the config file.

#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

Output_Dir <- set_dir(Figures)
setwd(Output_Dir)

#**______________________________________________________________________**#
## -------------------- 1. GENOTYPE SUBSETTING  ------------------------- ##

#Assign unique snp-trait pair identification column in sig snps df
sig_snps <- sig_snps %>%
  mutate(snp_trait = paste(term,trait, sep = "-"))

sig_snps_list = sig_snps_grouped$term
sig_traits_list = unique(sig_snps$trait)

sig_snps_cols = intersect(names(Database[["Complete Datasets"]][["Full Cohort"]]), sig_snps_list)
sig_traits_cols = intersect(names(Database[["Complete Datasets"]][["Full Cohort"]]), sig_traits_list)

#Subset datasets for specific snp - trait pairs to do anova tests

adult_barplotDF <- Database[["Complete Datasets"]][["Adults"]] %>%
  select(pid,bmi, age, sex, all_of(
    intersect(
      names(Database[["Complete Datasets"]][["Adults"]]),
      c(sig_snps %>% filter(cohort == "adults") %>% pull(trait) %>% unique(),
        sig_snps %>% filter(cohort == "adults") %>% pull(term) %>% unique() )
    )
  ))

children_barplotDF <- Database[["Complete Datasets"]][["Children"]] %>%
  select(pid, bmi, age, sex, all_of(
    intersect(
      names(Database[["Complete Datasets"]][["Children"]]),
      c(sig_snps %>% filter(cohort == "children") %>% pull(trait) %>% unique(),
        sig_snps %>% filter(cohort == "children") %>% pull(term) %>% unique() )
    )
  ))

full_barplotDF <- Database[["Complete Datasets"]][["Full Cohort"]] %>%
  select(pid, age, sex, all_of(
    intersect(
      names(Database[["Complete Datasets"]][["Full Cohort"]]),
      c(sig_snps %>% filter(cohort == "full") %>% pull(trait) %>% unique(),
        sig_snps %>% filter(cohort == "full") %>% pull(term) %>% unique() )
    )
  ))

#Now subsetting by genotype for each significant SNP column (0,1,2), saving to lists. We will define a function to do this.

create_snp_subsets <- function(data) {

  snp_cols <- get_snp_columns(data)
  
  subsets <- list()
  
  for (snp in snp_cols) {
    subsets[[snp]] <- list(
      recessive = data %>% filter(!!sym(snp) == 0),
      heterozygous = data %>% filter(!!sym(snp) == 1),
      dominant = data %>% filter(!!sym(snp) == 2)
    )
  }
  
  return(subsets)
}

adult_snp_subsets <- create_snp_subsets(adult_barplotDF)
children_snp_subsets <- create_snp_subsets(children_barplotDF)
full_snp_subsets <- create_snp_subsets(full_barplotDF)


#**______________________________________________________________________**#
## -------------------- 2. STATS BY GENPOTYPE  -------------------------- ##

#Loop over a given df subset list (ie. adult_snp_subsets) and compute statistics for each trait (if multiple)

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


library(dplyr)
library(ggplot2)
library(rlang)
library(ggsignif)
library(ggpubr)

## rs2948276_regularity

rs2948276_regularity <- calculate_stats(full_barplotDF, "Regularity", "rs2948276")

stats_data <- rs2948276_regularity$stats
anova_results <- rs2948276_regularity$anova


ancova_result <- aov(Regularity ~ factor(rs2948276) + age + factor(sex), data = full_barplotDF)


tukey_result <- TukeyHSD(ancova_result, "factor(rs2948276)")

print(tukey_result)

p_line_1 <- tibble(
  x = c("0", "0", "2", "2"),
  y = c(0.85, 0.94, 0.94, 0.90)
)

p_line_2 <- tibble(
  x = c("1", "1", "2", "2"),
  y = c(0.80, 0.84, 0.84, 0.80)
)

custom_colors <- c("0" = "#EBDDE6",    
                   "1" = "#B599B3",    
                   "2" = "#885075")  

custom_labels <- c("0" = "AA",   
                   "1" = "AG",   
                   "2" = "GG")  

# Add a new column for the mean values (rounded or formatted)
stats_data <- rs2948276_regularity$stats
stats_data$mean_label <- paste0(round(stats_data$mean, 2))

# Create the plot
rs2948276_regularity_plot <- ggplot(stats_data, aes(x = factor(rs2948276), y = mean)) +
  geom_col(aes(fill = factor(rs2948276)), color = "black", width = 0.85) +
  geom_errorbar(aes(ymin = ci_upper,
                    ymax = ci_lower),
                color = "#22292F",
                width = .1) +
  geom_line(data = p_line_1, 
            aes(x = x, y = y, group = 1)) +
  geom_line(data = p_line_2, 
            aes(x = x, y = y, group = 1)) +
  annotate("text", x = 2, y = 0.954, 
           label = "**",
           size = 14, color = "#313745") +
  annotate("text", x = 2.5, y = 0.843, 
           label = "*",
           size = 14, color = "#313745") +
  # Add mean labels for all bars except the last one
  geom_text(data = subset(stats_data, rs2948276 != 2), 
            aes(label = mean_label, y = mean - 0.2), 
            size = 8, color = "#313745") +
  # Hard code the mean label for the last bar (GG) in white
  geom_text(data = subset(stats_data, rs2948276 == 2), 
            aes(label = mean_label, y = mean - 0.2), 
            size = 8, color = "white") +
  scale_x_discrete(labels = custom_labels) +  
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  guides(fill = FALSE) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 38, color = "#313745",
                              face = "bold",
                              margin = margin(b = 35), hjust = .85),
    plot.margin = unit(rep(1, 4), "cm"),
    axis.text = element_text(size = 35, color = "#313745"),
    axis.title = element_text(size = 35, hjust = .50, color = "#313745"),
    axis.title.x = element_text(margin = margin(t = 25)),
    axis.title.y = element_text(margin = margin(r = 35)),
    axis.text.y = element_text(margin = margin(r = 25)),
    axis.text.x = element_text(margin = margin(t = 10)),
    plot.caption = element_text(size = 16, 
                                face = "italic",
                                color = "#76829E",
                                margin = margin(t = 15))
  ) +
  labs(
    x = "Genotype",
    y = "Sleep Regularity",
    title = "Sleep Regularity by rs2948276 Genotype"
  )

# Save the plot
svg(filename = "rs2948276_regularity_plot.svg", width = 12, height = 8.2, bg = "transparent")
rs2948276_regularity_plot
dev.off()


#COMPUTE STATS FOR NEW PLOT
rs2307111_Mean_TST <- calculate_stats(children_barplotDF, "Mean_TST", "rs2307111") 

ancova_result <- aov(Mean_TST ~ factor(rs2307111) + age + bmi + factor(sex), data = children_barplotDF)
summary(ancova_result)

tukey_result <- TukeyHSD(ancova_result, "factor(rs2307111)")

print(tukey_result)

stats_data <- rs2307111_Mean_TST$stats

p_line_1 <- tibble(
  x = c("0", "0", "2", "2"),
  y = c(605, 640, 640, 605)
)

# Add a new column for the mean values (rounded or formatted)
stats_data <- rs2307111_Mean_TST$stats
stats_data$mean_label <- paste0(round(stats_data$mean, 1))

rs2307111_Mean_TST_plot <- ggplot(stats_data, aes(x = factor(rs2307111), y = mean)) +
  geom_col(aes(fill = factor(rs2307111)), color = "black", width = 0.85) +
  geom_errorbar(aes(ymin = ci_upper,
                    ymax = ci_lower),
                color = "#22292F",
                width = .1) +
  geom_line(data = p_line_1, 
            aes(x = x, y = y, group = 1)) +
  annotate("text", x = 2, y = 643.3, 
           label = "*",
           size = 14, color = "#313745") +
  # Add mean labels for all bars except the last one
  geom_text(data = subset(stats_data, rs2307111 != 2), 
            aes(label = mean_label, y = mean - 100), 
            size = 8, color = "#313745") +
  # Hard code the mean label for the last bar (GG) in white
  geom_text(data = subset(stats_data, rs2307111 == 2), 
            aes(label = mean_label, y = mean - 100), 
            size = 8, color = "white") +
  scale_x_discrete(labels = custom_labels) +  
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(limits = c(0, 700), expand = c(0, 0)) +
  guides(fill = FALSE) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 38,color = "#313745",
                              face = "bold",
                              margin = margin(b = 35), hjust = .85),
    plot.margin = unit(rep(1, 4), "cm"),
    axis.text = element_text(size = 35, color = "#313745"),
    axis.title = element_text(size = 35, hjust = .50, color = "#313745"),
    axis.title.x = element_text(margin = margin(t = 25)),
    axis.title.y = element_text(margin = margin(r = 35)),
    axis.text.y = element_text(margin = margin(r = 25)),
    axis.text.x = element_text(margin = margin(t = 10)),
    plot.caption = element_text(size = 16, 
                                face = "italic",
                                color = "#76829E",
                                margin = margin(t = 15))
  ) +
  labs(
    x = "Genotype",
    y = "Total Sleep Time (m)",
    title = "Sleep Time by rs2307111 Genotype"
  )

rs2307111_Mean_TST_plot

svg(filename = "rs2307111_Mean_TST_plot.svg", width = 12, height = 8.2, bg = "transparent")
rs2307111_Mean_TST_plot
dev.off()


#COMPUTE STATS FOR NEW PLOT

custom_colors <- c("0" = "#EBEDF1",    
                   "1" = "#BBC1CF",    
                   "2" = "#76829E")  

rs310727_carb_prop <- calculate_stats(full_barplotDF, "carb_prop", "rs310727") 

ancova_result <- aov(carb_prop ~ factor(rs310727) + age + factor(sex), data = full_barplotDF)
summary(ancova_result)

tukey_result <- TukeyHSD(ancova_result, "factor(rs310727)")

print(tukey_result)

stats_data <- rs310727_carb_prop$stats

p_line_1 <- tibble(
  x = c("0", "0", "2", "2"),
  y = c(57, 61, 61, 57)
)

p_line_2 <- tibble(
  x = c("1", "1", "2", "2"),
  y = c(53, 55, 55, 52)
)


custom_labels <- c("0" = "CC",   
                   "1" = "CT",   
                   "2" = "TT")  

# Add a new column for the mean values (rounded or formatted)
stats_data <- rs310727_carb_prop$stats
stats_data$mean_label <- paste0(round(stats_data$mean * 100, 1), "%")



rs310727_carb_prop_plot <- ggplot(stats_data, aes(x = factor(rs310727), y = (mean*100))) +
  geom_col(aes(fill = factor(rs310727)), color = "black", width = 0.85) +
  geom_errorbar(aes(ymin = ci_upper*100,
                    ymax = ci_lower*100),
                color = "#22292F",
                width = .1) +
  geom_line(data = p_line_1, 
            aes(x = x, y = y, group = 1)) +
  annotate("text", x = 2, y = 62.4, 
           label = "*",
           size = 14, color = "#313745") +
  geom_line(data = p_line_2, 
            aes(x = x, y = y, group = 1)) +
  annotate("text", x = 2.5, y = 56.4, 
           label = "*",
           size = 14, color = "#313745") +
  geom_text(aes(label = mean_label, y = mean * 100 -10),
            size = 8, color = "#313745") +
  scale_x_discrete(labels = custom_labels) +  
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(limits = c(0, 70), expand = c(0, 0)) +
  guides(fill = FALSE) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 34,color = "#313745",
                              face = "bold",
                              margin = margin(b = 35), hjust = .85),
    plot.margin = unit(rep(1, 4), "cm"),
    axis.text = element_text(size = 35, color = "#313745"),
    axis.title = element_text(size = 28, hjust = .50, color = "#313745"),
    axis.title.x = element_text(size = 35, margin = margin(t = 25)),
    axis.title.y = element_text(margin = margin(r = 35)),
    axis.text.y = element_text(margin = margin(r = 25)),
    axis.text.x = element_text(margin = margin(t = 10)),
    plot.caption = element_text(size = 16, 
                                face = "italic",
                                color = "#76829E",
                                margin = margin(t = 15))
  ) +
  labs(
    x = "Genotype",
    y = "% Energy From Carbohydrates",
    title = "% Energy from Carbohyrates by rs310727 Genotype"
  )

rs310727_carb_prop_plot

svg(filename = "rs310727_carb_prop_plot.svg", width = 12, height = 8.2, bg = "transparent")
rs310727_carb_prop_plot
dev.off()


custom_colors <- c("0" = "#EBEDF1",    
                   "1" = "#BBC1CF",    
                   "2" = "#76829E")  

custom_labels <- c("0" = "TT",   
                   "1" = "TG",   
                   "2" = "GG")  


  #COMPUTE STATS FOR NEW PLOT
  rs1154155_tfat_prop <- calculate_stats(full_barplotDF, "tfat_prop", "rs1154155") 
  
  ancova_result <- aov(tfat_prop ~ factor(rs1154155) + age + factor(sex), data = full_barplotDF)
  summary(ancova_result)
  
  tukey_result <- TukeyHSD(ancova_result, "factor(rs1154155)")
  
  print(tukey_result)
  
  p_line_1 <- tibble(
    x = c("0", "0", "2", "2"),
    y = c(38,52,52,49)
  )
  
  
  stats_data <- rs1154155_tfat_prop$stats
  
  
  custom_labels <- c("0" = "TT",   
                     "1" = "TG",  
                     "2" = "GG")  
  
  
  # Add a new column for the mean values (rounded or formatted)
  stats_data <- rs1154155_tfat_prop$stats
  stats_data$mean_label <- paste0(round(stats_data$mean * 100, 1), "%")
  
  # Create the plot with mean labels added above the error bars
  rs1154155_tfat_prop_plot <- ggplot(stats_data, aes(x = factor(rs1154155), y = mean * 100)) +
    geom_col(aes(fill = factor(rs1154155)), color = "black", width = 0.85) +
    geom_errorbar(aes(ymin = ci_upper * 100,
                      ymax = ci_lower * 100),
                  color = "#22292F",
                  width = .1) +
    geom_line(data = p_line_1, 
              aes(x = x, y = y, group = 1)) +
    annotate("text", x = 2, y = 54.3, 
             label = "*",
             size = 14, color = "#313745") +
    geom_text(aes(label = mean_label, y = mean * 100 -10),
              size = 8, color = "#313745") +  # Adjust size and color as needed
    scale_x_discrete(labels = custom_labels) +  
    scale_fill_manual(values = custom_colors) +
    scale_y_continuous(limits = c(0, 60), expand = c(0, 0)) +
    guides(fill = FALSE) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 38, color = "#313745",
                                face = "bold",
                                margin = margin(b = 35), hjust = .85),
      plot.margin = unit(rep(1, 4), "cm"),
      axis.text = element_text(size = 35, color = "#313745"),
      axis.title = element_text(size = 35, hjust = .5, color = "#313745"),
      axis.title.x = element_text(margin = margin(t = 25)),
      axis.title.y = element_text(margin = margin(r = 35)),
      axis.text.y = element_text(margin = margin(r = 25)),
      axis.text.x = element_text(margin = margin(t = 10)),
      plot.caption = element_text(size = 16, 
                                  face = "italic",
                                  color = "#76829E",
                                  margin = margin(t = 15))
    ) +
    labs(
      x = "Genotype",
      y = "% Energy From Fat",
      title = "% Energy From Fat by rs1154155 genotype"
    )
  
  # Display the plot
  rs1154155_tfat_prop_plot
  
  # Save the plot
  svg(filename = "rs1154155_tfat_prop_plot.svg", width = 12, height = 8.2, bg = "transparent")
  rs1154155_tfat_prop_plot
  dev.off()
  
  
  
#**_____________________________________________________________________**#
## -------------------- 2. CORRELATION MATRIX -------------------------- ##


full_cohort <- Database[["Complete Datasets"]][["Full Cohort"]] %>%
    select(-pid,-fid, -cohort, all_of(diet_response_columns), all_of(sleep_response_columns), all_of(covariates), -all_of(snp_columns))
  
corr_cols <- c(diet_response_columns,sleep_response_columns, covariates)

full_cohort[-which(names(full_cohort) == "sex")] <- lapply(full_cohort[-which(names(full_cohort) == "sex")], as.numeric)

full_cohort$sex <- ifelse(full_cohort$sex == "M", 1, ifelse(full_cohort$sex == "F", 0, NA))

name_replacements <- c(
  kcal = "Energy Intake",
  carb_prop = "% Energy from Carbohydrate",
  prot_prop = "% Energy from Protein",
  tfat_prop = "% Energy from Fat",
  carb = "Total Carb",
  prot = "Total Protein",
  tfat = "Total Fat",
  Mean_TST = "Total Sleep Time",
  Mean_SE = "Sleep Efficiency",
  Mean_WASO = "Awake During Night",
  Regularity = "Sleep Regularity"
)

full_cohort$bmi <- as.numeric(full_cohort$bmi)

names(full_cohort)[names(full_cohort) %in% names(name_replacements)] <- name_replacements[names(full_cohort)]

corr_matrix <- cor(full_cohort, use="complete.obs")
library(corrplot)

#Correlation plot of response variables
custom_col <- colorRampPalette(c("#B4C9B8", "#FFFFFF", "#313745"))(200)

corrsy <- corrplot(corr_matrix, method = "circle", col = custom_col, tl.col = "black", tl.srt = 35, tl.cex = 1.5, cl.cex = 1.4)
corrsy


svg(filename = "corrplot_transparent.svg", width = 14, height = 10, bg = "transparent")
corrsy <- corrplot(corr_matrix, method = "circle", col = custom_col, tl.col = "#313745", tl.srt = 35, tl.cex = 1.5, cl.cex = 1.5)
dev.off()

