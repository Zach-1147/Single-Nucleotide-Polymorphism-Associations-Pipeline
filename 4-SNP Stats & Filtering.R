
#** ------ SNP STATS & FILTERING -------- **#
# ========================================= #
#    Zach Ribau  |  June 20, 2024       

#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

Input_Dir <- set_dir(SNP_Stats_Filtering)

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
  )
)

#Clean Env
rm(adult_snp_stats, children_snp_stats, Genotype_long_adult, Genotype_long_children,ld_D_prime_adult, ld_R_squared_adult, ld_results_adult, Genotype_matrix_adult, ld_results_children, ld_D_prime_child, ld_R_squared_child, Genotype_matrix_children, Genotype_adult, Genotype_children)

#NEED TO REVIEW DOCUMENTATION FOR NAs in resulting dataframes here.
