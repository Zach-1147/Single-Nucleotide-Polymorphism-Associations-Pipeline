#** ------ SNP SELECTION -------- **#
# ================================= #
#    Zach Ribau  |  June 23, 2024       

#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

Input_Dir <- set_dir(SNP_Selection)
Output_Dir <- set_dir(Data_Integration)

#**___________________________________________________**#
## --------------------- DATA ------------------------ ##

setwd(Input_Dir)

#Read in data

#Load EBI GWAS catalogue
EBI <- read_tsv("EBI.tsv")

#Load in trait terms with trait group assignments table
sleep_traits <- read.csv("sleep_traits.csv")
diet_traits <- read.csv("diet_traits.csv")

#Read in GFHS Genotype file
genotype_file <- read.table("genotype_file.txt", skip=9)

#And rs ID identification table
rsid_conversions <- read.table("locus_rsids.txt", header = TRUE)

#**_____________________________________________________**#
## ------------- 1. INITIAL GWAS FILTERING ------------- ##

#Create sleep associations filtering EBI to include only sleep traits

EBI_sleep <- EBI %>%
  filter(`MAPPED_TRAIT` %in% sleep_traits$Trait)

#Create diet associations df by filtering EBI to include only dietary traits as outlined in provided diet traits file

EBI_diet <- EBI %>%
  filter(`MAPPED_TRAIT` %in% diet_traits$Trait)

#Generate a studies dataframe of unique contributing studies for sleep triats with a few relevant annotations for interest.

#Now specifically sleep studies
sleep_studies <- EBI_sleep %>%
  select(STUDY, LINK, PUBMEDID, DATE, `FIRST AUTHOR`, JOURNAL) %>%
  distinct(STUDY, .keep_all = TRUE)

#Generate a dataframe with intersecting SNPs between EBI_sleep and EBI_diet dataframes, first creating a vector of intersecting SNPs from diet and sleep associations dfs

common_SNPs <- intersect(EBI_sleep$SNPS, EBI_diet$SNPS)

#Filter EBI_sleep for these common SNPs
common_sleep_SNPs <- EBI_sleep[EBI_sleep$SNPS %in% common_SNPs, c("SNPS", "MAPPED_GENE", "MAPPED_TRAIT")]

#Filter EBI_diet for these common SNPs
common_diet_SNPs <- EBI_diet[EBI_diet$SNPS %in% common_SNPs, c("SNPS", "MAPPED_GENE", "MAPPED_TRAIT")]

#Combine above dataframes and reformat
sleep_diet_SNPs <- rbind(common_sleep_SNPs, common_diet_SNPs) %>%
  group_by(SNPS, MAPPED_GENE) %>%
  summarise(traits = paste(unique(MAPPED_TRAIT), collapse = ", "),.groups = 'drop')

#Generate a simple dataframe list of the final set of SNPs of interest
snp_list <- sleep_diet_SNPs$SNPS

#**_____________________________________________________**#
## ----------- 2. GENOTYPE FILE PREPROCESSING ---------- ##

#Now create df with SNP and a new blank column to look up alt names for each SNP to ensure matches where they exist in the genotype file.
SNP_List <- data.frame(SNP = snp_list$SNPS, Alt_name = NA)

#Trim whitespace to ensure no matches are missed due to this
SNP_List$SNP <- trimws(SNP_List$SNP)

#Define a number of splits to split conversion RSID columns, as many are lists
max_splits <- max(sapply(strsplit(rsid_conversions$Most_Likely_RSID.s., ","), length))

#Split into new columns for lists in Most_Likely_RSID.s.
Conv2 <- separate(rsid_conversions, col = "Most_Likely_RSID.s.", into = paste0("PossibleRSID_", 1:max_splits), sep = ",", remove = FALSE, fill = "right")

#Convert to data.tables
setDT(SNP_List)
setDT(Conv2)

#Reshape Conv2 to long for computationally efficient match strategy.
SNPs_to_Confirm <- melt(Conv2, id.vars = "Manifest_Name", measure.vars = patterns("^PossibleRSID_"), value.name = "SNP", variable.name = "variable")[
  SNP %in% SNP_List$SNP, .(SNP, Manifest_Name)] #Filter for SNPs that are in the SNP_List of interest, and return table with look up values.

as.character(SNPs_to_Confirm$SNP)
as.character(SNP_List$SNP)

#**_____________________________________________________**#
## -------- 3. EBI-GENOTYPE FILE SNP MATCHING ---------- ##

#Populate Alt_name column of SNP_List with matching manifest names
(indices <- match(SNP_List$SNP, SNPs_to_Confirm$SNP))

#Use the indices to map Manifest_Name to Alt_name, NA where no match is found
SNP_List$Alt_name <- ifelse(!is.na(indices), SNPs_to_Confirm$Manifest_Name[indices], SNP_List$Alt_name)

genotype_file <- genotype_file %>%
  rownames_to_column(var = "SNPS")

#Filter for prescence of the values from either column to check for matches in the final report.
SNPs_Confirmed <- SNP_List %>%
  filter(SNP %in% genotype_file$SNPS | Alt_name %in% genotype_file$SNPS)

#**_____________________________________________________**#
## -------------- 4. FINAL SNP SELECTION --------------- ##

#Finally, generated a fitered_genotype file with only the snps of interest

filter_genotype <- genotype_file %>%
  filter(SNPS %in% SNPs_Confirmed$SNP)

setwd(Output_Dir)

write.csv(filter_genotype, "filtered_genotype_file.csv")

colnames(genotype_file) <- colnames(filter_genotype_file)

write.table(genotype_file, "genotype_file3.txt", row.names = FALSE, quote = FALSE, sep = "\t")
