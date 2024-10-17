
#** ------ SNP STATS -------- **#
# ========================================= #
#    Zach Ribau  |  June 20, 2024       
rm(list=ls())
#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

Input_Dir <- set_dir(SNP_Stats)
Output_Dir <- set_dir(LD_MAF)

#**___________________________________________________**#
## --------------------- DATA ------------------------ ##

setwd(Input_Dir)

#Read in data
adult_processed <- read.csv("adult_processed.csv")
child_processed <- read.csv("children_processed.csv")


#**_____________________________________________________**#
## --------------- 1. BASIC SNP STATS  ----------------- ##

## --------------- 1.1 Adult Dataset ------------------- ##

snp_columns <- get_snp_columns(adult_processed) #Retrieve SNP columns

#Create a long dataframe for adult cohort's genotype data
Genotype_adult <- adult_processed %>%
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
  
  #Calculate allele frequencies and proportions
  combined_alleles <- bind_rows(
    filtered_df %>% group_by(Allele_1) %>% summarize(count = n()) %>% rename(allele = Allele_1),
    filtered_df %>% group_by(Allele_2) %>% summarize(count = n()) %>% rename(allele = Allele_2)
  ) %>%
    group_by(allele) %>%
    summarize(
      total_count = sum(count)
    ) %>%
    mutate(proportion = total_count / (nrow(filtered_df) * 2))
  
  #Determine reference and alternate alleles
  if (length(unique(filtered_df$Genotype)) < 2) {
    Ref_allele <- combined_alleles$allele[1]  # Assign the only allele as Ref_allele
    Ref_allele_prop <- 1  #All are reference alleles
    Alt_allele <- "/"
    Alt_allele_prop <- 0
  } else {
    # If there are at least 2 unique genotypes, calculate reference and alternate alleles and their proportions
    Ref_allele <- combined_alleles$allele[which.max(combined_alleles$proportion)]
    Alt_allele <- combined_alleles$allele[which.min(combined_alleles$proportion)]
    Ref_allele_prop <- combined_alleles$proportion[which.max(combined_alleles$proportion)]
    Alt_allele_prop <- combined_alleles$proportion[which.min(combined_alleles$proportion)]
  }
  #Calculate genotype counts
  filtered_df <- filtered_df %>%
    mutate(num_genotype = case_when(
      Allele_1 == Ref_allele & Allele_2 == Ref_allele ~ 0,  #Homozygous reference
      Allele_1 == Alt_allele & Allele_2 == Alt_allele ~ 2,  #Homozygous alternate
      Allele_1 != Allele_2 ~ 1  #Heterozygous
    ))
  
  genotype_counts <- filtered_df %>%
    group_by(num_genotype) %>%
    summarize(count = n())
  
  #Extract individual counts for each genotype
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

#Now we will ensure that the minor allele identified here is consistent with that reported in a representative population from the 1000 genomes project, along with the relative porportion.

#now navigte to SNP_Nexus_All dir under data files in project dir
setwd(file.path(project_dir, "Data Files", "SNP_Nexus_All"))

#load in population_stats.txt
population_stats <- read.table("population_stats.txt", header = TRUE, sep = "\t")

#remove row 28
population_stats <- population_stats[-28,]

#now take only the Ref.Allele and Alt.Allele, Minor.Allele, and all columnes ending with "Frequency". 

#Ensure all column names are enclosed in '' in the select statement
population_stats <- population_stats %>% 
  select('Variation.ID','REF.Allele', 'ALT.Allele', 'Minor.Allele', ends_with('Frequency'))

#We will merge this into our adult_snp_stats df to compare the minor allele frequency
adult_snp_stats <- adult_snp_stats %>% 
  left_join(population_stats, by = c("rsID" = "Variation.ID"))

wrong_minor <- adult_snp_stats %>%  
  filter(Alt_Allele != Minor.Allele) %>%
  filter(Alt_Allele != "/")


#Flag variants not different due to complemetary base pair reads, but rather different proportions and due to close to 50/50. We want to replace these and recompute some stats.

wrong_minor <- adult_snp_stats %>%  
  filter(rsID %in% c("rs2307111", "rs800165"))

#Those that are due to complementary base pair reads will  be corrected in the final dataframes with a python script, as well as in the snp_stats df's produced from this script. 

#initialize a df to store changed rows
snp_df <- data.frame(
  rsID = character(),
  hzREF = integer(),
  htz = integer(),
  hzALT = integer(),
  stringsAsFactors = FALSE
)

for (rsID in wrong_minor$rsID) {

  #get unique occuring genotypes in adult_pre-proccessed
  genotypes <- unique(adult_processed[[rsID]])

  #get corresponding minor allele from adult_snp_stats
  minor_allele <- adult_snp_stats$Alt_Allele[adult_snp_stats$rsID == rsID]

  ref_allele <- adult_snp_stats$Ref_Allele[adult_snp_stats$rsID == rsID]

  #construct our hzREF, hzALT and htz

  hzREF_id <- paste(ref_allele, ref_allele, sep = "")
  hzALT_id <- paste(minor_allele, minor_allele, sep = "")
  hzt1 <- paste(ref_allele, minor_allele, sep = "")
  hzt2 <- paste(minor_allele, ref_allele, sep = "")

  #now count the occurences of each of these in the genotype columns

  hzREF <- sum(adult_processed[[rsID]] == hzREF_id)
  hzALT <- sum(adult_processed[[rsID]] == hzALT_id)
  htz <- sum(adult_processed[[rsID]] == hzt1) + sum(adult_processed[[rsID]] == hzt2)

  #return a df containing the rsID and each count as a column

  #Create a data frame for the current rsID
  current_snp_df <- data.frame(
    rsID = rsID,
    hzREF = hzREF,
    htz = htz,
    hzALT = hzALT,
    stringsAsFactors = FALSE
  )

  #Append the current data frame to the main snp_df
  snp_df <- rbind(snp_df, current_snp_df)

  #now change the values in the adult_snp_stats
  adult_snp_stats$hzREF[adult_snp_stats$rsID == rsID] <- hzREF
  adult_snp_stats$htz[adult_snp_stats$rsID == rsID] <- htz
  adult_snp_stats$hzALT[adult_snp_stats$rsID == rsID] <- hzALT

  #and also invert the Minor Allele and Ref Allele
  adult_snp_stats$Alt_Allele[adult_snp_stats$rsID == rsID] <- ref_allele

  #and also invert the Minor Allele and Ref Allele
  adult_snp_stats$Ref_Allele[adult_snp_stats$rsID == rsID] <- minor_allele

  #finally invert the ref_allele_prop and alt_allele_prop

  adult_snp_stats$Ref_allele_prop[adult_snp_stats$rsID == rsID] <- 1 - adult_snp_stats$Ref_allele_prop[adult_snp_stats$rsID == rsID]

  adult_snp_stats$Alt_allele_prop[adult_snp_stats$rsID == rsID] <- 1 - adult_snp_stats$Alt_allele_prop[adult_snp_stats$rsID == rsID]

}

#now output a list of minor alleles due to complenatary base pair

comp_base_adult <- adult_snp_stats %>%  
  filter(Alt_Allele != Minor.Allele) %>%
  filter(Alt_Allele != "/") %>%
  #dont include this one, because it is due to a different proportion, not a complementary base pair
  filter(rsID != "rs800165")

comp_base_adult <- list(comp_base_adult$rsID)

#write list to txt in project_dir (which is a variable)
python_dir <- file.path(project_dir, "/Data Files/For Python")

write.table(comp_base_adult, file.path(python_dir, "comp_base_adult.txt"), row.names = FALSE, col.names = FALSE)

#Changing the above, along with the snp_stats df's produced from this script, will be done with a python script before saving data to the final output location.

## --------------- 1.2 Childrens Dataset ------------------- ##

#Repeating for childrens data

#Reassign snp columns for children dataset, although will be identical
snp_columns <- get_snp_columns(child_processed)

#Create a long dataframe for childrens genotype data
Genotype_children <- child_processed %>%
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

#We will merge this into our adult_snp_stats df to compare the minor allele frequency
children_snp_stats <- children_snp_stats %>% 
  left_join(population_stats, by = c("rsID" = "Variation.ID"))

#Flag variants not different due to complemetary base pair reads, but rather different proportions and due to close to 50/50. We want to replace these and recompute some stats.
wrong_minor <- children_snp_stats %>%  
  filter(Alt_Allele != Minor.Allele) %>%
  filter(Alt_Allele != "/")

#here we have rs7155227, along with the ones from before, rs2307111 and rs800165. Filter wrong_minor for only these.

wrong_minor <- wrong_minor %>% 
  filter(rsID %in% c("rs7155227", "rs2307111", "rs800165"))

view(wrong_minor)
#Those that are due to complementary base pair reads will  be corrected in the final dataframes with a python script, as well as in the snp_stats df's produced from this script. 

#initialize a df to store changed rows
snp_df <- data.frame(
  rsID = character(),
  hzREF = integer(),
  htz = integer(),
  hzALT = integer(),
  stringsAsFactors = FALSE
)

for (rsID in wrong_minor$rsID) {

  #get unique occuring genotypes in adult_pre-proccessed
  genotypes <- unique(child_processed[[rsID]])

  #get corresponding minor allele from adult_snp_stats
  minor_allele <- children_snp_stats$Alt_Allele[children_snp_stats$rsID == rsID]

  ref_allele <- children_snp_stats$Ref_Allele[children_snp_stats$rsID == rsID]

  #construct our hzREF, hzALT and htz

  hzREF_id <- paste(ref_allele, ref_allele, sep = "")
  hzALT_id <- paste(minor_allele, minor_allele, sep = "")
  hzt1 <- paste(ref_allele, minor_allele, sep = "")
  hzt2 <- paste(minor_allele, ref_allele, sep = "")

  #now count the occurences of each of these in the genotype columns

  hzREF <- sum(child_processed[[rsID]] == hzREF_id)
  hzALT <- sum(child_processed[[rsID]] == hzALT_id)
  htz <- sum(child_processed[[rsID]] == hzt1) + sum(child_processed[[rsID]] == hzt2)

  #return a df containing the rsID and each count as a column

  #Create a data frame for the current rsID
  current_snp_df <- data.frame(
    rsID = rsID,
    hzREF = hzREF,
    htz = htz,
    hzALT = hzALT,
    stringsAsFactors = FALSE
  )

  #Append the current data frame to the main snp_df
  snp_df <- rbind(snp_df, current_snp_df)

  #now change the values in the adult_snp_stats
  children_snp_stats$hzREF[children_snp_stats$rsID == rsID] <- hzREF
  children_snp_stats$htz[children_snp_stats$rsID == rsID] <- htz
  children_snp_stats$hzALT[children_snp_stats$rsID == rsID] <- hzALT

  #and also invert the Minor Allele and Ref Allele
  children_snp_stats$Alt_Allele[children_snp_stats$rsID == rsID] <- ref_allele

  #and also invert the Minor Allele and Ref Allele
  children_snp_stats$Ref_Allele[children_snp_stats$rsID == rsID] <- minor_allele

  #finally invert the ref_allele_prop and alt_allele_prop

  children_snp_stats$Ref_allele_prop[children_snp_stats$rsID == rsID] <- 1 - children_snp_stats$Ref_allele_prop[children_snp_stats$rsID == rsID]

  children_snp_stats$Alt_allele_prop[children_snp_stats$rsID == rsID] <- 1 - children_snp_stats$Alt_allele_prop[children_snp_stats$rsID == rsID]

}

#Note, no additional SNPs with flipped base pair reads need to be added, since it is the same from the adults - as long as the pythons script corrects it in both dataframes.

#write the SNP stats df's to the python_dir as well

#remerge with pop_stats

#now ensure that the Ref_Allele and Minor_Allele are consistentbetween the two dataframes

adult_check <- adult_snp_stats %>% 
  select(rsID, Ref_Allele, Alt_Allele)

children_check <- children_snp_stats %>%
  rename(Ref_Allele_c = Ref_Allele, Alt_Allele_c = Alt_Allele) %>%
  select(rsID, Ref_Allele_c, Alt_Allele_c)

#merge these twoby rsID
check_merge <- left_join(adult_check, children_check, by = "rsID")

#Now create a column that says "Y" if Ref_Allele and Ref_Allele_c are different, or IF the Alt_Allele and Alt_Allele_c are different.

check_merge <- check_merge %>% 
  mutate(Ref_Allele_diff = ifelse(Ref_Allele != Ref_Allele_c, "Y", "N"),
         Alt_Allele_diff = ifelse(Alt_Allele != Alt_Allele_c, "Y", "N"))

#count if "Y"
sum(check_merge$Ref_Allele_diff == "Y")

#if zero, we can move along.

setwd(python_dir)

#convert Genotype col in both df's to character
adult_snp_stats$Genotypes <- as.character(adult_snp_stats$Genotypes)

children_snp_stats$Genotypes <- as.character(children_snp_stats$Genotypes)

write.csv(adult_snp_stats, "adult_snp_stats.csv", row.names = FALSE)
write.csv(children_snp_stats, "children_snp_stats.csv", row.names = FALSE)

write.csv(adult_processed, "adult_processed.csv", row.names = FALSE)
write.csv(child_processed, "child_processed.csv", row.names = FALSE)
