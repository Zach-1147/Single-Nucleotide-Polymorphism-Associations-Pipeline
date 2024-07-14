
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

#Now we can generate a summary table describing each SNP in more detail for the adult cohort first.
adult_snp_stats <- Genotype_long_adult %>%
  filter(!(Allele_1 == "-" | Allele_2 == "-")) %>%
  group_by(rsID) %>%
  summarise(
    Genotypes = paste(unique(Genotype), collapse = ", "),
    
    G_prop = (sum(Allele_1 == "G") + sum(Allele_2 == "G")) / (2 * (n())),
    A_prop = (sum(Allele_1 == "A") + sum(Allele_2 == "A")) / (2 * (n())),
    T_prop = (sum(Allele_1 == "T") + sum(Allele_2 == "T")) / (2 * (n())),
    C_prop = (sum(Allele_1 == "C") + sum(Allele_2 == "C")) / (2 * (n())),
    
    Ref_Allele = case_when(
      G_prop == max(G_prop, A_prop, T_prop, C_prop) ~ "G",
      A_prop == max(G_prop, A_prop, T_prop, C_prop) ~ "A",
      T_prop == max(G_prop, A_prop, T_prop, C_prop) ~ "T",
      C_prop == max(G_prop, A_prop, T_prop, C_prop) ~ "C"
    ),
    
    Alt_Allele = case_when(
      sum(c(G_prop, A_prop, T_prop, C_prop) > 0) == 1 ~ "/",
      sum(c(G_prop, A_prop, T_prop, C_prop) > 0) > 2 ~ {
        temp <- sort(c(G = G_prop, A = A_prop, T = T_prop, C = C_prop), decreasing = TRUE)
        paste(names(temp)[2], "/", names(temp)[3], sep = "")
      },
      TRUE ~ {
        temp <- c(G = G_prop, A = A_prop, T = T_prop, C = C_prop)
        names(sort(temp, decreasing = TRUE))[2]
      }
    ),
    
    Ref_Allele_Prop = case_when(
      Ref_Allele == "A" ~ A_prop,
      Ref_Allele == "G" ~ G_prop,
      Ref_Allele == "C" ~ C_prop,
      Ref_Allele == "T" ~ T_prop
    ),
    
    Minor_Allele_Prop = case_when(
      sum(c(G_prop, A_prop, T_prop, C_prop) > 0) == 1 ~ as.character(0.00),
      sum(c(G_prop, A_prop, T_prop, C_prop) > 0) > 2 ~ {
        temp <- sort(c(G = G_prop, A = A_prop, T = T_prop, C = C_prop), decreasing = TRUE)
        paste(as.character(temp[2]), "/", as.character(temp[3]), sep = "")
      },
      TRUE ~ {
        temp <- c(G = G_prop, A = A_prop, T = T_prop, C = C_prop)
        as.character(min(temp[temp > 0]))
      }
    )
  )

## --------------- 1.2 Childrens Dataset ------------------- ##

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

#And generating the stats dataframe
children_snp_stats <- Genotype_long_children %>%
  filter(!(Allele_1 == "-" | Allele_2 == "-")) %>%
  group_by(rsID) %>%
  summarise(
    Genotypes = paste(unique(Genotype), collapse = ", "),
    
    G_prop = (sum(Allele_1 == "G") + sum(Allele_2 == "G")) / (2 * (n())),
    A_prop = (sum(Allele_1 == "A") + sum(Allele_2 == "A")) / (2 * (n())),
    T_prop = (sum(Allele_1 == "T") + sum(Allele_2 == "T")) / (2 * (n())),
    C_prop = (sum(Allele_1 == "C") + sum(Allele_2 == "C")) / (2 * (n())),
    
    Ref_Allele = case_when(
      G_prop == max(G_prop, A_prop, T_prop, C_prop) ~ "G",
      A_prop == max(G_prop, A_prop, T_prop, C_prop) ~ "A",
      T_prop == max(G_prop, A_prop, T_prop, C_prop) ~ "T",
      C_prop == max(G_prop, A_prop, T_prop, C_prop) ~ "C"
    ),
    
    Alt_Allele = case_when(
      sum(c(G_prop, A_prop, T_prop, C_prop) > 0) == 1 ~ "/",
      sum(c(G_prop, A_prop, T_prop, C_prop) > 0) > 2 ~ {
        temp <- sort(c(G = G_prop, A = A_prop, T = T_prop, C = C_prop), decreasing = TRUE)
        paste(names(temp)[2], "/", names(temp)[3], sep = "")
      },
      TRUE ~ {
        temp <- c(G = G_prop, A = A_prop, T = T_prop, C = C_prop)
        names(sort(temp, decreasing = TRUE))[2]
      }
    ),
    
    Ref_Allele_Prop = case_when(
      Ref_Allele == "A" ~ A_prop,
      Ref_Allele == "G" ~ G_prop,
      Ref_Allele == "C" ~ C_prop,
      Ref_Allele == "T" ~ T_prop
    ),
    
    Minor_Allele_Prop = case_when(
      sum(c(G_prop, A_prop, T_prop, C_prop) > 0) == 1 ~ as.character(0.00),
      sum(c(G_prop, A_prop, T_prop, C_prop) > 0) > 2 ~ {
        temp <- sort(c(G = G_prop, A = A_prop, T = T_prop, C = C_prop), decreasing = TRUE)
        paste(as.character(temp[2]), "/", as.character(temp[3]), sep = "")
      },
      TRUE ~ {
        temp <- c(G = G_prop, A = A_prop, T = T_prop, C = C_prop)
        as.character(min(temp[temp > 0]))
      }
    )
  )

#Filtering based on MAF
exclude_snps_adult <- adult_snp_stats %>%
  filter(Minor_Allele_Prop < minor_allele_threshold) %>%
  rename(MAF_adult = "Minor_Allele_Prop") %>%
  select(rsID, MAF_adult)

exclude_snps_children <- children_snp_stats %>%
  filter(Minor_Allele_Prop < minor_allele_threshold) %>%
  rename(MAF_children = "Minor_Allele_Prop") %>%
  select(rsID, MAF_children)

#This will be used to remove snps at the end of the script
exclude_snps_df <- full_join(exclude_snps_adult, exclude_snps_children, by = "rsID") %>%
  distinct(rsID, MAF_children, MAF_adult)

exclude_snps <- as.character(exclude_snps_df$rsID)

#**_____________________________________________________**#
## ------------ 2. LINKAGE DIS-EQUILIBRIUM ------------- ##

## ----------------- 2.1 Adult Dataset ----------------- ##

#Extract Genotype information only, keeping in wide format
Genotype_matrix_adult <- as.matrix(adult_preprocessed %>%
  select(all_of(snp_columns)))

#Set rownames from pid in original dataset
rownames(Genotype_matrix_adult) <- adult_preprocessed$pid

#Apply function from CONFIG file to encode genotypes numerically. Note that this matrix will have NA wherever there is missing genotype information.
Genotype_matrix_adult <- apply(Genotype_matrix_adult, 2, encode_genotype)

#Now using snpStats package to compute LD stats
snp_matrix <- new("SnpMatrix", Genotype_matrix_adult)
colnames(snp_matrix) <- colnames(Genotype_matrix_adult)
rownames(snp_matrix) <- rownames(Genotype_matrix_adult)

#Now, calculating ld and viewing results for all SNPs pairwise
ld_results_adult <- ld(snp_matrix, depth=10, stats=c("D.prime", "R.squared"))

#Generating dataframes from results
ld_D_prime_adult <- as.data.frame(as.matrix(ld_results_adult[["D.prime"]]))
ld_R_squared_adult <- as.data.frame(as.matrix(ld_results_adult[["R.squared"]]))

Genetic_STATS <- list(genotype_matrix_adult = as.data.frame(Genotype_matrix_adult),geno_long_adult = Genotype_long_adult, geno_long_children = Genotype_long_children, 
                      
SNP_Stats = list(adult = adult_snp_stats,child = children_snp_stats),
                      
LD_Stats = list(D_prime_adult = ld_D_prime_adult,R_sqr_adult = ld_R_squared_adult)

)

## ----------------- 2.2 Childrens Dataset ----------------- ##

#Extract Genotype information only, keeping in wide format
Genotype_matrix_children <- as.matrix(child_preprocessed %>%
                                     select(all_of(snp_columns)))

#Set rownames from pid in original dataset
rownames(Genotype_matrix_children) <- child_preprocessed$pid

#Apply function from CONFIG file to encode genotypes numerically. Note that this matrix will have NA wherever there is missing genotype information.
Genotype_matrix_children <- apply(Genotype_matrix_children, 2, encode_genotype)

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

#**_____________________________________________________**#
## -------------- 3. CRITERIA FILTERING ---------------- ##

#Filtering out low minor allele frequency based on config variable
Genotype_matrix_adult <- Genotype_matrix_adult %>%
  select(-all_of(exclude_snps))

#Filtering out low minor allele frequency based on config variable
Genotype_matrix_children <- Genotype_matrix_children %>%
  select(-all_of(exclude_snps))

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

#Adding to database list object
Database <- c(Database, list('Complete Datasets' = list("Adults" = adult_data, "Children" = children_data, "Full Cohort" = full_cohort_data)))

rm(adult_merge,children_merge, full_cohort_data, adult_data, children_data)

#NEED TO REVIEW DOCUMENTATION FOR NAs in resulting dataframes here.

#AT THE END

#Clean Env

rm(adult_snp_stats, children_snp_stats, Genotype_long_adult, Genotype_long_children,ld_D_prime_adult, ld_R_squared_adult, ld_results_adult, Genotype_matrix_adult, ld_results_children, ld_D_prime_child, ld_R_squared_child, Genotype_matrix_children, Genotype_adult, Genotype_children, adult_preprocessed, child_preprocessed, exclude_snps_adult, exclude_snps_children, exclude_snps_df)
