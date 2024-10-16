#** ------ DATASET CLEANING ------- **#
# =================================== #
#    Zach Ribau  |  June 18, 2024    

#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

Input_Dir <- set_dir(Dataset_Cleaning)
Output_Dir <- set_dir(SNP_Stats_Filtering)

#**_____________________________________________________**#
## --------------------- DATA -------------------------- ##

setwd(Input_Dir)

#Read in data
adult_dataset <- read.csv("final_adult.csv")
children_dataset <- read.csv("final_child.csv")

#We will now remove participants that have notable amounts of missing data across any of the essential dietary, sleep and genetic data columns

#**_____________________________________**#
## ----- 1. MISSING DIETARY DATA ------- ##

#First make a general function to check for any form of missing data in a given column
check_missing <- function(column) {
  is.na(column) | column == "NA" | column == "" | column == "N/A (pregnant)" #count this last one for bmi specifically in adult dataset
}

#Examine how many participants are missing data in any of the dietary columns.

#Make a list of dietary columns to check
diet_columns <- c("kcal", "prot", "tfat", "carb", "caff", "carb_prop", "prot_prop", "tfat_prop")

#Now apply our function over these columns in both children and adults, extracting into filtered dataframes for examination
Missing_diet_child_df <- children_dataset %>%
  filter(apply(select(., all_of(diet_columns)), 1, function(row) any(check_missing(row))))

Missing_diet_adult_df <- adult_dataset %>%
  filter(apply(select(., all_of(diet_columns)), 1, function(row) any(check_missing(row))))

adult_dataset <- adult_dataset %>%
  rename(age = "age_ha") #provide standard name for age and other cols between adult and childrens data

children_dataset <- children_dataset %>%
  rename(age = "age_years") %>%
  rename(bmi = "bmi_z")

#Also checking on our covariate columns
Missing_demos_child_df <- children_dataset %>%
  filter(apply(select(., all_of(covariates)), 1, function(row) any(check_missing(row))))

Missing_demos_adult_df <- adult_dataset %>%
  filter(apply(select(., all_of(covariates)), 1, function(row) any(check_missing(row))))

#Now filtering for missing data in diet and demos columns

children_dataset <- children_dataset %>%
  filter(!(pid %in% Missing_diet_child_df$pid)) %>%
  filter(!(pid %in% Missing_demos_child_df$pid))

adult_dataset <- adult_dataset %>%
  filter(!(pid %in% Missing_diet_adult_df$pid)) %>%
  filter(!(pid %in% Missing_demos_adult_df$pid))

#**_____________________________________**#
## ----- 2. MISSING SLEEP DATA --------- ##

#List our sleep columns and apply function to create df's for missing sleep data
sleep_columns <- c("Nights", "Mean_TST", "Mean_SE", "Mean_WASO", "Mean_TIB", "mean_IN_TIME_24", "mean_OUT_TIME_24", "SD_sleep_onset", "SD_sleep_offset", "SD_SE", "SD_TST", "Regularity")

#Looking at missing values in each sleep column
adult_sleep_nacount <- as.data.frame(colSums(is.na(adult_dataset[sleep_columns])))
colnames(adult_sleep_nacount) <- "NA_Count"
children_sleep_nacount <- as.data.frame(colSums(is.na(children_dataset[sleep_columns])))
colnames(children_sleep_nacount) <- "NA_Count"

#These columns require a closer examination before removing participants. For example, often only a single column is missing information (ie. time in bed). 

#Making dataframes to examine closer using our missing data function
Missing_sleep_child_df <- children_dataset %>%
  filter(apply(select(., all_of(sleep_columns)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(sleep_columns))

Missing_sleep_adult_df <- adult_dataset %>%
  filter(apply(select(., all_of(sleep_columns)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(sleep_columns))

#For now, we will remove only participants missing information across all of the sleep columns and examine if there are multiple columns with missing data or if it is just time in bed.
adult_dataset <- adult_dataset[!apply(adult_dataset[sleep_columns], 1, function(x) all(is.na(x))), ]

children_dataset <- children_dataset[!apply(children_dataset[sleep_columns], 1, function(x) all(is.na(x))), ]

#Checking again missing values per sleep column

adult_sleep_nacount <- as.data.frame(colSums(is.na(adult_dataset[sleep_columns])))

children_sleep_nacount <- as.data.frame(colSums(is.na(children_dataset[sleep_columns])))

#We see there are no missing values in the adult dataset, and only the time in bed column is missing values in the children dataset - which we will accept at this time.

#We will however remove participants with fewer than 3 nights logged
too_few_nights_a <- adult_dataset %>%
  filter(Nights <= 3) %>%
  select(pid, Nights)

too_few_nights_c <- children_dataset %>%
  filter(Nights <= 3) %>%
  select(pid, Nights)

too_few_nights <- rbind(too_few_nights_a,too_few_nights_c)

#Remove from datasets
adult_dataset <- adult_dataset %>%
  filter(!(pid %in% too_few_nights$pid))

children_dataset <- children_dataset  %>%
  filter(!(pid %in% too_few_nights$pid))

#Cleaning env.
rm(Missing_sleep_child_df,Missing_sleep_adult_df,children_sleep_nacount,adult_sleep_nacount)

#**________________________________________**#
## ----- 3. MISSING GENOTYPE DATA --------- ##

#Identify snp columns as those starting with "rs" and followed by numbers. We should ensure this produces 40 columns and is identical if we do it for the children and adult dataset

snp_columns <- get_snp_columns(adult_dataset)

#Incase SNPs with missing data are encoded somehow other than blanks or NAs, we will check unqiue values across all SNPs first. Create a function that can be applied across our snp column vector

list_unique_snps <- function(dataset, snp_columns) {
  unique_values <- unique(unlist(dataset[snp_columns]))
  data.frame(Unique_SNP_Values = unique_values)
}

#Now apply to both datasets
unique_snps_adult <- list_unique_snps(adult_dataset, snp_columns)

unique_snps_children <- list_unique_snps(children_dataset, snp_columns)

#Examining these, we see that there is "--" and "NA" occurring in the snp columms, outside of the expected nucleotide bases.

#We will redefine our function to check missing to include "--".
check_missing <- function(row) {
  any(is.na(row) | row == "NA" | row == "" | row == "--")
}

#Now we will examine dataframes where participants have any of the above in snp_columnns
Missing_snp_child_df <- children_dataset %>%
  filter(apply(select(., all_of(snp_columns)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(snp_columns))

Missing_snp_adult_df <- adult_dataset %>%
  filter(apply(select(., all_of(snp_columns)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(snp_columns))

#Now, we will filter out participants with missing genotype information across ALL snps before examining those with "--" or some NAs spread out.
adult_dataset <- adult_dataset[!apply(adult_dataset[snp_columns], 1, function(x) all(is.na(x))), ]

children_dataset <- children_dataset[!apply(children_dataset[snp_columns], 1, function(x) all(is.na(x))), ]

#Adjust the missing SNP dataframe for both datasets to have counts of NA and "--" by participant
Missing_snp_adult_df <- adult_dataset %>%
  select(pid, all_of(snp_columns)) %>%
  rowwise() %>%
  mutate(
    NA_Count = sum(is.na(c_across(all_of(snp_columns)))),
    `-- Count` = sum(c_across(all_of(snp_columns)) == "--")
  ) %>%
  ungroup() %>%
  filter(NA_Count + `-- Count` > 0) %>%
  select(pid, NA_Count, `-- Count`)

Missing_snp_child_df <- children_dataset %>%
  select(pid, all_of(snp_columns)) %>%
  rowwise() %>%
  mutate(
    NA_Count = sum(is.na(c_across(all_of(snp_columns)))),
    `-- Count` = sum(c_across(all_of(snp_columns)) == "--")
  ) %>%
  ungroup() %>%
  filter(NA_Count + `-- Count` > 0) %>%
  select(pid, NA_Count, `-- Count`)


children_dataset <- children_dataset %>%
  filter(!(pid %in% Missing_snp_child_df$pid))

adult_dataset <- adult_dataset %>%
  filter(!(pid %in% Missing_snp_adult_df$pid))

setwd(Output_Dir)

#Write output files into folder for next step of pipeline
write.csv(adult_dataset, "adult_preprocessed.csv", row.names = FALSE)
write.csv(children_dataset, "children_preprocessed.csv",row.names = FALSE)

setwd(Input_Dir)

excluded_adults <- read.csv("final_adult.csv") %>%
  filter(!(pid %in% adult_dataset$pid))

excluded_children <- read.csv("final_child.csv") %>%
  filter(!(pid %in% children_dataset$pid))

Exclusions <- list(Missing_Data = list(Adult = excluded_adults, Children = excluded_children), '<= 3 Nights Sleep' = list(Adult = too_few_nights_a, Children = too_few_nights_c))

#Generate some df's for keeping tidy in our environment
adult_sleep <- adult_dataset %>%
  select(pid,all_of(sleep_columns))

children_sleep <- children_dataset %>%
  select(pid,all_of(sleep_columns))

children_diet <- children_dataset %>%
  select(pid,all_of(sleep_columns))

adult_diet <- adult_dataset %>%
  select(pid,all_of(diet_columns))

children_diet <- children_dataset %>%
  select(pid,all_of(diet_columns))

#Define remaining columns not identified as demo_columns for demographic
demo_columns_adult <- setdiff(colnames(adult_dataset) , c(snp_columns, diet_columns, sleep_columns))

demo_columns_children <- setdiff(colnames(children_dataset) , c(snp_columns, diet_columns, sleep_columns))

adult_demos <- adult_dataset %>%
  select(pid, all_of(demo_columns_adult))
 

children_demos <- children_dataset %>%
  select(pid, all_of(demo_columns_children))


children_demos_na <- children_demos %>%
  filter(is.na(bmi) | bmi == " " | bmi == "NA")
nrow(children_demos_na)

adult_demos_na <- adult_demos %>%
  filter(is.na(bmi) | bmi == " " | bmi == "NA")
nrow(adult_demos_na)



'Sleep Data' <- list(Adult = adult_sleep, Children = children_sleep)
'Diet Data' <- list(Adult = adult_diet, Children = children_diet)
'Demographic Data' <- list(Adult = adult_demos, Children = children_demos)


#Clean Env
rm(adult_dataset,children_dataset,Missing_snp_adult_df,Missing_snp_child_df, unique_snps_adult, unique_snps_children, excluded_children, excluded_adults, Missing_diet_adult_df, Missing_diet_child_df, diet_columns,sleep_columns, adult_diet,adult_sleep,children_diet,children_sleep,adult_demos,children_demos, too_few_nights_a,too_few_nights_c, too_few_nights, Missing_demos_adult_df, Missing_demos_child_df, children_demos_na, adult_demos_na)
