
#** ------ SNP STATS & FILTERING -------- **#
# ========================================= #
#    Zach Ribau  |  June 20, 2024       

#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

Input_Dir <- set_dir(SNP_Stats_Filtering)
Output_Dir <- set_dir(Association_Testing)

#**___________________________________________________**#
## --------------------- DATA ------------------------ ##

setwd(Input_Dir)

#Read in data
adult_preprocessed <- read.csv("adult_preprocessed.csv")
child_preprocessed <- read.csv("children_preprocessed.csv")

#**_____________________________________________________**#
## --------------- 1. BASIC SNP STATS  ----------------- ##

## --------------- 1.1 Adult Dataset ------------------- ##

snp_columns <- get_snp_columns(adult_preprocessed) #Retrieve SNP columns

#Create a long dataframe for adult cohort's genotype data
Genotype_adult <- adult_preprocessed %>%
  select(pid, all_of(snp_columns)) 

Genotype_long_adult <- Genotype_adult %>% 
  pivot_longer( #Convert to long format, getting each participant-SNP combination as a row
    cols = -pid,  
    names_to = "rsID",
    values_to = "Genotype"
  ) %>%
  mutate( #We can pull out the first character present in each row, as well as the second, to get each Allele into separate columns that we can operate on to get allele level statistics
    Allele_1 = substr(Genotype, 1, 1), 
    Allele_2 = substr(Genotype, 2, 2)
  )

#Function to compute statistics for SNPs from genotype long
genotype_summary <- function(Genotype_long_adult, snp) {
  # Filter for the SNP of interest
  filtered_df <- Genotype_long_adult %>%
    filter(rsID == snp)
  
  # Calculate allele frequencies and proportions
  combined_alleles <- bind_rows(
    filtered_df %>% group_by(Allele_1) %>% summarize(count = n()) %>% rename(allele = Allele_1),
    filtered_df %>% group_by(Allele_2) %>% summarize(count = n()) %>% rename(allele = Allele_2)
  ) %>%
    group_by(allele) %>%
    summarize(
      total_count = sum(count)
    ) %>%
    mutate(proportion = total_count / (nrow(filtered_df) * 2))
  
  # Determine reference and alternate alleles
  if (length(unique(filtered_df$Genotype)) < 2) {
    Ref_allele <- combined_alleles$allele[1]  # Assign the only allele as Ref_allele
    Ref_allele_prop <- 1  # All are reference alleles
    Alt_allele <- "/"
    Alt_allele_prop <- 0
  } else {
    # If there are at least 2 unique genotypes, calculate reference and alternate alleles and their proportions
    Ref_allele <- combined_alleles$allele[which.max(combined_alleles$proportion)]
    Alt_allele <- combined_alleles$allele[which.min(combined_alleles$proportion)]
    Ref_allele_prop <- combined_alleles$proportion[which.max(combined_alleles$proportion)]
    Alt_allele_prop <- combined_alleles$proportion[which.min(combined_alleles$proportion)]
  }
  # Calculate genotype counts
  filtered_df <- filtered_df %>%
    mutate(num_genotype = case_when(
      Allele_1 == Ref_allele & Allele_2 == Ref_allele ~ 0,  # Homozygous reference
      Allele_1 == Alt_allele & Allele_2 == Alt_allele ~ 2,  # Homozygous alternate
      Allele_1 != Allele_2 ~ 1  # Heterozygous
    ))
  
  genotype_counts <- filtered_df %>%
    group_by(num_genotype) %>%
    summarize(count = n())
  
  # Extract individual counts for each genotype
  homozygous_ref <- genotype_counts$count[genotype_counts$num_genotype == 0] %>% ifelse(length(.) == 0, 0, .)
  heterozygous <- genotype_counts$count[genotype_counts$num_genotype == 1] %>% ifelse(length(.) == 0, 0, .)
  homozygous_alt <- genotype_counts$count[genotype_counts$num_genotype == 2] %>% ifelse(length(.) == 0, 0, .)
  
  genotypes <- unique(filtered_df$Genotype)
  
  # Create the output data frame
  snp_df <- data.frame(
    rsID = snp,
    Genotypes = I(list(genotypes)),  # Wrap in list to keep as one cell
    Ref_Allele = Ref_allele,
    Alt_Allele = Alt_allele,
    Ref_allele_prop = Ref_allele_prop,
    Alt_allele_prop = Alt_allele_prop,
    hzREF = homozygous_ref,
    htz = heterozygous,
    hzALT = homozygous_alt
  )
  
  return(snp_df)
}

adult_snp_stats <- do.call(rbind, lapply(snp_columns, function(rsID) genotype_summary(Genotype_long_adult, rsID)))

## --------------- 1.2 Childrens Dataset ------------------- ##

#Repeating for childrens data

#Reassign snp columns for children dataset, although will be identical
snp_columns <- get_snp_columns(child_preprocessed)

#Create a long dataframe for childrens genotype data
Genotype_children <- child_preprocessed %>%
  select(pid, all_of(snp_columns)) 

Genotype_long_children <- Genotype_children %>% 
  pivot_longer( #Convert to long format, getting each participant-SNP combination as a row
    cols = -pid,  
    names_to = "rsID",
    values_to = "Genotype"
  ) %>%
  mutate( #We can pull out the first character present in each row, as well as the second, to get each Allele into separate columns that we can operate on to get allele level statistics
    Allele_1 = substr(Genotype, 1, 1), 
    Allele_2 = substr(Genotype, 2, 2)
  )

children_snp_stats <- do.call(rbind, lapply(snp_columns, function(rsID) genotype_summary(Genotype_long_children, rsID)))


#FILTERING BASED ON MAF
exclude_snps_adult <- adult_snp_stats %>%
  filter(Alt_allele_prop < minor_allele_threshold) %>%
  rename(MAF_adult = "Alt_allele_prop") %>%
  select(rsID, MAF_adult)

exclude_snps_children <- children_snp_stats %>%
  filter(Alt_allele_prop < minor_allele_threshold) %>%
  rename(MAF_children = "Alt_allele_prop") %>%
  select(rsID, MAF_children)

#This will be used to remove snps at the end of the script
exclude_snps_df <- full_join(exclude_snps_adult, exclude_snps_children, by = "rsID") %>%
  distinct(rsID, MAF_children, MAF_adult)

exclude_snps <- as.character(exclude_snps_df$rsID)

#**_____________________________________________________**#
## ------------ 2. LINKAGE DIS-EQUILIBRIUM ------------- ##

## ----------------- 2.1 Adult Dataset ----------------- ##

#Extract Genotype information only, keeping in wide format
Genotype_matrix_adult <- adult_preprocessed %>%
  select(all_of(snp_columns))

Genotype_matrix_adult <- as.matrix(Genotype_matrix_adult)

#Set rownames from pid in original dataset
rownames(Genotype_matrix_adult) <- adult_preprocessed$pid

# Get the intersection of rsID between both datasets
common_rsIDs <- intersect(colnames(Genotype_matrix_adult), adult_snp_stats$rsID)

# Subset Genotype_matrix_adult to keep only common SNPs
Genotype_matrix_adult <- Genotype_matrix_adult[, common_rsIDs]

# Subset adult_snp_stats to keep only common SNPs

# Choose a specific SNP (e.g., rs10995245)
snp_columns_adult <- colnames(Genotype_matrix_adult)

for (rsID in snp_columns_adult) {
  
  # Extract the column for the current SNP
  snp_column <- Genotype_matrix_adult[, rsID]
  
  #Apply the encode_genotype function and update the column with the encoded values
  Genotype_matrix_adult[, rsID] <- encode_genotype(snp_column, rsID, adult_snp_stats)
  
}


#Now using snpStats package to compute LD stats
snp_matrix <- new("SnpMatrix", Genotype_matrix_adult)
colnames(snp_matrix) <- colnames(Genotype_matrix_adult)
rownames(snp_matrix) <- rownames(Genotype_matrix_adult)

#Now, calculating ld and viewing results for all SNPs pairwise
ld_results_adult <- ld(snp_matrix, depth=3, stats=c("D.prime", "R.squared"))

#Generating dataframes from results
ld_D_prime_adult <- as.data.frame(as.matrix(ld_results_adult[["D.prime"]]))
ld_R_squared_adult <- as.data.frame(as.matrix(ld_results_adult[["R.squared"]]))

Genetic_STATS <- list(genotype_matrix_adult = as.data.frame(Genotype_matrix_adult),geno_long_adult = Genotype_long_adult, geno_long_children = Genotype_long_children, 
                      
                      SNP_Stats = list(adult = adult_snp_stats,child = children_snp_stats),
                      
                      LD_Stats = list(D_prime_adult = ld_D_prime_adult,R_sqr_adult = ld_R_squared_adult)
                      
)

ld_R_squared_adult <- ld_R_squared_adult %>% 
  rownames_to_column(var = "Variant1")

#getting long format
ld_R_squared_adult_long <- ld_R_squared_adult %>% 
  pivot_longer(cols = -Variant1, names_to = "Variant2", values_to = "R^2_score")

#remove self comparisons (not neccesary but makes things clearer)
long_ld_R_squared_adult_longdf <- ld_R_squared_adult_long %>% 
  filter(Variant1 != Variant2)

long_ld_R_squared_adult_longdf$`R^2_score`<- as.numeric(long_ld_R_squared_adult_longdf$`R^2_score`)

long_ld_R_squared_adult_longdf <- ld_R_squared_adult_long %>% 
  filter(`R^2_score` > 0.2)

## ----------------- 2.2 Childrens Dataset ----------------- ##

#Extract Genotype information only, keeping in wide format
#Extract Genotype information only, keeping in wide format
Genotype_matrix_children <- child_preprocessed %>%
  select(all_of(snp_columns))

Genotype_matrix_children <- as.matrix(Genotype_matrix_children)

#Set rownames from pid in original dataset
rownames(Genotype_matrix_children) <- child_preprocessed$pid

# Get the intersection of rsID between both datasets
common_rsIDs <- intersect(colnames(Genotype_matrix_children), children_snp_stats$rsID)

# Subset Genotype_matrix_adult to keep only common SNPs
Genotype_matrix_children <- Genotype_matrix_children[, common_rsIDs]

# Subset adult_snp_stats to keep only common SNPs

# Choose a specific SNP (e.g., rs10995245)
snp_columns_children <- colnames(Genotype_matrix_children)

for (rsID in snp_columns_children) {
  
  # Extract the column for the current SNP
  snp_column <- Genotype_matrix_children[, rsID]
  
  # Apply the encode_genotype function and update the column with the encoded values
  Genotype_matrix_children[, rsID] <- encode_genotype(snp_column, rsID, children_snp_stats)
  
}

#Now using snpStats package to compute LD stats
snp_matrix_children <- new("SnpMatrix", Genotype_matrix_children)
colnames(snp_matrix_children) <- colnames(Genotype_matrix_children)
rownames(snp_matrix_children) <- rownames(Genotype_matrix_children)

#Now, calculating ld and viewing results for all SNPs pairwise
ld_results_children <- ld(snp_matrix_children, depth=10, stats=c("D.prime", "R.squared"))

#Generating dataframes from results
ld_D_prime_child <- as.data.frame(as.matrix(ld_results_children[["D.prime"]]))
ld_R_squared_child <- as.data.frame(as.matrix(ld_results_children[["R.squared"]]))

Genotype_matrix_adult <- rownames_to_column(as.data.frame(Genotype_matrix_adult), var = "pid")
Genotype_matrix_children <- rownames_to_column(as.data.frame(Genotype_matrix_children), var = "pid")

ld_R_squared_child <- ld_R_squared_child %>% 
  rownames_to_column(var = "Variant1")

#putting in long format
ld_R_squared_child_long <- ld_R_squared_child %>% 
  pivot_longer(cols = -Variant1, names_to = "Variant2", values_to = "R^2_score")

#Removing self comparisons
long_ld_R_squared_child_longdf <- ld_R_squared_child_long %>% 
  filter(Variant1 != Variant2)

long_ld_R_squared_child_longdf$`R^2_score`<- as.numeric(long_ld_R_squared_child_longdf$`R^2_score`)

long_ld_R_squared_child_longdf <- ld_R_squared_child_long %>% 
  filter(`R^2_score` > 0.2)

#**_____________________________________________________**#
## -------------- 3. CRITERIA FILTERING ---------------- ##

#Filtering out low minor allele frequency based on config variable
Genotype_matrix_adult <- Genotype_matrix_adult %>%
  select(-all_of(exclude_snps))

#Filtering out low minor allele frequency based on config variable
Genotype_matrix_children <- Genotype_matrix_children %>%
  select(-all_of(exclude_snps))

cat("\n")
cat("\n")
print(paste(as.character(nrow(exclude_snps_df)), "SNPs excluded based on minor allele frequency < ", minor_allele_threshold))

#**_____________________________________________________**#
## -------------- 4. DATA ORGANIZATION ----------------- ##

#Save everything to orderly list in environment
Genetic_STATS <- list(
  Adults = list(
    "Genotype Table" = Genotype_adult,
    "Genotype Matrix" = as.data.frame(Genotype_matrix_adult),
    "Long Genotype" = Genotype_long_adult, 
    SNP_Stats = list(
      "Allelic Stats" = adult_snp_stats,
      "LD D-Prime" = ld_D_prime_adult,
      "LD R-Sqr" = ld_R_squared_adult
    )
  ),
  
  Children = list(
    "Genotype Table" = Genotype_children,
    "Genotype Matrix" = as.data.frame(Genotype_matrix_children),
    "Long Genotype" = Genotype_long_children, 
    SNP_Stats = list(
      "Allelic Stats" = children_snp_stats,
      "LD D-Prime" = ld_D_prime_child,
      "LD R-Sqr" = ld_R_squared_child
    )
  ), 
  
  'Excluded SNPs' = exclude_snps_df
)

Database <- list("Genetic Stats" = Genetic_STATS, "Diet Data" = `Diet Data`, "Sleep Data" = `Sleep Data`, "Demographic Data" = `Demographic Data`, "Exclusions" = Exclusions)

rm(Genetic_STATS, `Diet Data`, `Sleep Data`, `Demographic Data`, Exclusions)

#Read in genotype matrix data from environment
a_genotype <- Database[["Genetic Stats"]][["Adults"]][["Genotype Matrix"]]
a_diet <- Database[["Diet Data"]][["Adult"]]
a_sleep <- Database[["Sleep Data"]][["Adult"]]
a_demos <- Database[["Demographic Data"]][["Adult"]]

#Ensure PID is identical for all subsets of the data
all.equal(a_genotype$pid, a_diet$pid, a_sleep$pid,a_demos$pid)

#Merge into complete dataframe
adult_data <- merge(a_diet, a_sleep, by = "pid", all.x = TRUE)
adult_data <- merge(adult_data, a_demos, by = "pid", all.x = TRUE)
adult_data <- merge(adult_data, a_genotype, by = "pid", all.x = TRUE)

#And now for childrens data
#Read in genotype matrix data from environment
c_genotype <- Database[["Genetic Stats"]][["Children"]][["Genotype Matrix"]]

#Read in genotype matrix data from environment
c_diet <- Database[["Diet Data"]][["Children"]]
c_sleep <- Database[["Sleep Data"]][["Children"]]
c_demos <- Database[["Demographic Data"]][["Children"]]

all.equal(c_genotype$pid, c_diet$pid, c_sleep$pid, c_demos$pid)

children_data <- merge(c_diet, c_sleep, by = "pid", all.x = TRUE)
children_data <- merge(children_data, c_demos, by = "pid", all.x = TRUE)
children_data <- merge(children_data, c_genotype, by = "pid", all.x = TRUE)

#Clean env
rm(a_genotype,a_diet,a_demos,a_sleep,c_genotype,c_diet,c_sleep,c_demos)

#reassign snp_columns with excluded snps removed
snp_columns <- get_snp_columns(adult_data)

adult_merge <- adult_data %>%
  select(pid,fid,all_of(diet_response_columns), all_of(sleep_response_columns), all_of(snp_columns), all_of(covariates)) %>%
  mutate(cohort = "Adults") #Add dataset indicator to support population clustering of siblings and parents

children_merge <- children_data %>%
  select(pid,fid,all_of(diet_response_columns), all_of(sleep_response_columns), all_of(snp_columns), all_of(covariates)) %>%
  mutate(cohort = "Children")

full_cohort_data <- rbind(adult_merge, children_merge)

#Convert BMI column to numeric for each final dataset
adult_data$bmi <- as.numeric(adult_data$bmi)
children_data$bmi <- as.numeric(children_data$bmi)
full_cohort_data$bmi <- as.numeric(full_cohort_data$bmi)

#Adding to database list object
Database <- c(Database, list('Complete Datasets' = list("Adults" = adult_data, "Children" = children_data, "Full Cohort" = full_cohort_data)))
citation("snpStats")
rm(adult_merge,children_merge, full_cohort_data, adult_data, children_data)

#AT THE END

#Clean Env

rm(adult_snp_stats, children_snp_stats, Genotype_long_adult, Genotype_long_children,ld_D_prime_adult, ld_R_squared_adult, ld_results_adult, Genotype_matrix_adult, ld_results_children, ld_D_prime_child, ld_R_squared_child, Genotype_matrix_children, Genotype_adult, Genotype_children, adult_preprocessed, child_preprocessed, exclude_snps_adult, exclude_snps_children, exclude_snps_df)




