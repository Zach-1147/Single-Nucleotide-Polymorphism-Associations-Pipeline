
#** -------- DATA INTEGRATION --------- **#
# ======================================= #
#    Zach Ribau  |  June 18, 2024       

#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

Input_Dir <- set_dir(Data_Integration)
Output_Dir <- set_dir(Dataset_Cleaning)

#**_____________________________________________________________**#
## --------------- 1. DIET/DEMOGRAPHIC DATA -------------------- ##

## --- LOADING DATA --- ##

#Here we are loading in the diet/demo spreadsheets, doing some quality control and preparing dataframes for merge with genotype and sleep data. 


setwd(Input_Dir)

adult <- read.csv("adult_cohort.csv") #adult demographic/diet data, provided from Revital

children <- read.csv("children_cohort.csv") #Since there are duplicates of PID in this, we take rows where "t1" is in timepoint, since they contain the detailed data of interest and there appears to be a t1 and an NA for each unique pid.

#**_____________________________________________________________**#
## ----------------- 2. GENOTYPE DATA PREP --------------------- ##

#Here we are loading in the filtered genotype file and preparing it as a dataframe for merging with diet/demo and sleep datasets.

genotype_file <- read.csv("filtered_genotype_file.csv") #Prepared from genotype file processing script

#Transpose final genotype file
genotype_file <- as.data.frame(t(genotype_file))
colnames(genotype_file) <- genotype_file[1,] #Set SNP ids as column names
genotype_file <- genotype_file[-1,] #Remove SNP id column now
genotype_file <- genotype_file %>% rownames_to_column(var = "pid")
colnames(genotype_file) <- genotype_file[1, ]
genotype_file <- genotype_file[-1, ] %>%
  rename(pid = "SNPS")

#finally, replace O's with 0's for compatability with diet/demo data
genotype_file$pid <- sub("O", "0", genotype_file$pid)

#confirm structure
class(genotype_file)
dim(genotype_file)
print(head(genotype_file, 3))

#**_____________________________________________________________**#
## ----------------- 3. SLEEP DATA PREP ------------------------ ##

## LOAD IN SLEEP DATA ##

#load in excel workbook sheets by their name from "raw_sleep_excel.xlsx"
preschooler_schoolage <- read_excel("raw_sleep_excel.xlsx", sheet = "preschooler_schoolage")
toddler <- read_excel("raw_sleep_excel.xlsx", sheet = "toddler")
parent <- read_excel("raw_sleep_excel.xlsx", sheet = "parent")

## PRE-PROCCESSING ##

#Our first goal is to ensure that dates are formatted correctly and that they make sense relative to sleep and wake times reported. We will go 1 by 1 over each dataset

#---------preschoolers and school age children---------#

#Ensure the IN_DATE column is numeric
if (!is.numeric(preschooler_schoolage$IN_DATE)) {
  temp_in_date <- as.numeric(preschooler_schoolage$IN_DATE)
} else {
  temp_in_date <- preschooler_schoolage$IN_DATE
}

#Convert the numeric IN_DATE values to Date format and store in new_in_date
preschooler_schoolage$new_in_date <- as.Date(temp_in_date, origin = "1899-12-30")

#REPEAT for OUT_DATE
if (!is.numeric(preschooler_schoolage$OUT_DATE)) {
  temp_out_date <- as.numeric(preschooler_schoolage$OUT_DATE)
} else {
  temp_out_date <- preschooler_schoolage$OUT_DATE
}

#Convert the numeric OUT_DATE values to Date format and store in new_out_date
preschooler_schoolage$new_out_date <- as.Date(temp_out_date, origin = "1899-12-30")

#Now for rows where new_in_date and new_out_date are empty, it is because the original IN_DATE and OUT_DATE had 12/11/2017 format representation as strings - so we can fill those rows with the correct conversion to match the rest

preschooler_schoolage$new_in_date[is.na(preschooler_schoolage$new_in_date)] <- as.Date(preschooler_schoolage$IN_DATE[is.na(preschooler_schoolage$new_in_date)], format = "%m/%d/%Y")

#repeating for OUT_DATE
preschooler_schoolage$new_out_date[is.na(preschooler_schoolage$new_out_date)] <- as.Date(preschooler_schoolage$OUT_DATE[is.na(preschooler_schoolage$new_out_date)], format = "%m/%d/%Y")

#Ensure IN_TIME and OUT_TIME are in POSIXct format
preschooler_schoolage$IN_TIME <- as.POSIXct(preschooler_schoolage$IN_TIME, format = "%Y-%m-%d %H:%M:%S")
preschooler_schoolage$OUT_TIME <- as.POSIXct(preschooler_schoolage$OUT_TIME, format = "%Y-%m-%d %H:%M:%S")

#Extract and format the time part in AM/PM format
preschooler_schoolage$IN_TIME <- format(preschooler_schoolage$IN_TIME, format = "%I:%M %p")
preschooler_schoolage$OUT_TIME <- format(preschooler_schoolage$OUT_TIME, format = "%I:%M %p")

#filter a new df for which new_in_date and new_out_date are the same
same_day <- preschooler_schoolage %>%
  filter(new_in_date == new_out_date)

#Now retain from same_day ONLY rows for which IN_TIME and OUT_TIME contain a different AM vs PM status (ie both are not PM, both are not AM)
same_day <- same_day %>%
  filter((grepl("AM", IN_TIME) & !grepl("AM", OUT_TIME)) | (!grepl("AM", IN_TIME) & grepl("AM", OUT_TIME)))

#Also filter a df "missing_date" for which any one of new_in_date or new_out_date are missing
missing_date <- preschooler_schoolage %>%
  filter(is.na(new_in_date) | is.na(new_out_date))

# --------- TODDLERS --------- #

#Create new column called new_in_date and new_out_date for toddler.

#print the data type of IN_DATE_TIME and OUT_DATE_TIME columns
class(toddler$IN_DATE_TIME)
class(toddler$OUT_DATE_TIME)

#both are character, so we will convert 11/4/2017 19:17 formatting to just date WITHOUT time
toddler$new_in_date <- as.Date(toddler$IN_DATE_TIME, format = "%m/%d/%Y")

toddler$new_out_date <- as.Date(toddler$OUT_DATE_TIME, format = "%m/%d/%Y")

#check if any are missing data or unexpected format
missing_date_toddler <- toddler %>%
  filter(is.na(new_in_date) | is.na(new_out_date))

#Now extract just the time part from IN_DATE_TIME and OUT_DATE_TIME
toddler$IN_TIME <- format(strptime(toddler$IN_DATE_TIME, format = "%m/%d/%Y %H:%M"), format = "%I:%M %p")
toddler$OUT_TIME <- format(strptime(toddler$OUT_DATE_TIME, format = "%m/%d/%Y %H:%M"), format = "%I:%M %p")

#and make a toddler_same df with new_in_date and new_out_date the same
toddler_same <- toddler %>%
  filter(new_in_date == new_out_date)

#Now retain from same_day ONLY rows for which IN_TIME and OUT_TIME contain a different AM vs PM status (ie both are not PM, both are not AM)

toddler_same <- toddler_same %>%
  filter((grepl("AM", IN_TIME) & !grepl("AM", OUT_TIME)) | (!grepl("AM", IN_TIME) & grepl("AM", OUT_TIME)))

#--------- PARENTS ---------#

#Create new column called new_in_date and new_out_date for parent DF

#check the data type of IN_DATE and OUT_DATE
class(parent$IN_DATE_TIME)
class(parent$OUT_DATE_TIME)

#These are also character and same format seemingly as the toddler df, so we will convert to date
parent$new_in_date <- as.Date(parent$IN_DATE_TIME, format = "%m/%d/%Y")

parent$new_out_date <- as.Date(parent$OUT_DATE_TIME, format = "%m/%d/%Y")

#check if any are missing data or unexpected format
missing_date_parent <- parent %>%
  filter(is.na(new_in_date) | is.na(new_out_date))

#looks good, so now we will extract the time part from IN_DATE_TIME and OUT_DATE_TIME
parent$IN_TIME <- format(strptime(parent$IN_DATE_TIME, format = "%m/%d/%Y %H:%M"), format = "%I:%M %p")
parent$OUT_TIME <- format(strptime(parent$OUT_DATE_TIME, format = "%m/%d/%Y %H:%M"), format = "%I:%M %p")

#And now filtering for new_in_date and new_out_date the same
parent_same <- parent %>%
  filter(new_in_date == new_out_date)

#Now retain from same_day ONLY rows for which IN_TIME and OUT_TIME contain a different AM vs PM status (ie both are not PM, both are not AM)
parent_same <- parent_same %>%
  filter((grepl("AM", IN_TIME) & !grepl("AM", OUT_TIME)) | (!grepl("AM", IN_TIME) & grepl("AM", OUT_TIME)))

## COMBINING DATA ##

#Now we will attempt to combnine the data from the three df's into one main one with the same columns as the sleep df. We want to make a pid col in each df, starting with toddlers.

#We need to extract the part of the Subject column that contains either A0 or B0 followed by a series of three numbers. Often this occurs at the start of the string, but it could be in the middle or at the end and flanked by underscores.

#We will use a regular expression to extract this part of the string and store it in a new column called pid.
toddler <- toddler %>%
  mutate(pid = str_extract(Subject, "A0[0-9]{3}|B0[0-9]{3}"))

#Lets check how many unique subject IDs there are
length(unique(toddler$Subject))

#there are 83, so we should have 83 pid's
length(unique(toddler$pid))

#It looks like one Subject ID is missing a pid, so we will investgiate why that is.

#Turns out there were two "C0"s like FUL_BL_C0526, so we will correct that and reassign the pid
toddler$pid[toddler$Subject == "FUL_BL_C0576"] <- "C0526"
#And correcting the other C0 as well from before
toddler$pid[toddler$Subject == "FUL_BL_C0526"] <- "C0526"

#These should both = 83 now
length(unique(toddler$pid))
length(unique(toddler$Subject))

#Now investigating the parents (preschoolers are already formatted correctly)

#using the reg expression again this time looking for "P0" followed by 3 numbers, OR "S0"
parent <- parent %>%
  mutate(pid = str_extract(Subject, "P0[0-9]{3}|S0[0-9]{3}"))

#check how many are NAs
sum(is.na(parent$pid))

#Ensure the unique id lengths match
length(unique(parent$pid))
length(unique(parent$Subject))

#Now we will rename some columns to match the original sleep df

#In toddlers, rename the pid to username
toddler <- toddler %>%
  rename(username = pid)

#do the same for parents
parent <- parent %>%
  rename(username = pid)

#rename Subject to username in preschooler_schoolage
preschooler_schoolage <- preschooler_schoolage %>%
  rename(username = Subject)

#Now, in all three df's add a dataset column that contains parents in parent, preschoolers in preschooler_schoolage, and toddlers in toddler
toddler$dataset <- "toddler"
parent$dataset <- "parent"
preschooler_schoolage$dataset <- "preschooler_schoolage"

#Now rename our converted date columns (new_in_date and new_out_date) to IN_DATE_TIME and OUT_DATE_TIME in all df's

#Fist remove IN_DATE_TIME and OUT_DATE_TIME from each df
toddler <- toddler %>% select(-IN_DATE_TIME, -OUT_DATE_TIME)
parent <- parent %>% select(-IN_DATE_TIME, -OUT_DATE_TIME)

toddler <- toddler %>%
  rename(IN_DATE_TIME = new_in_date, OUT_DATE_TIME = new_out_date)

parent <- parent %>%
  rename(IN_DATE_TIME = new_in_date, OUT_DATE_TIME = new_out_date)

preschooler_schoolage <- preschooler_schoolage %>%
  rename(IN_DATE_TIME = new_in_date, OUT_DATE_TIME = new_out_date)

#Now we want IN_TIME_24 and OUT_TIME_24 to show the IN_TIME and OUT_TIME columns in 24-hour format, which currently are in AM/PM format.

toddler <- toddler %>%
  mutate(
    IN_TIME_24 = format(strptime(IN_TIME, format = "%I:%M %p"), format = "%H:%M"),
    OUT_TIME_24 = format(strptime(OUT_TIME, format = "%I:%M %p"), format = "%H:%M")
  )

parent <- parent %>%
  mutate(
    IN_TIME_24 = format(strptime(IN_TIME, format = "%I:%M %p"), format = "%H:%M"),
    OUT_TIME_24 = format(strptime(OUT_TIME, format = "%I:%M %p"), format = "%H:%M")
  )

preschooler_schoolage <- preschooler_schoolage %>%
  mutate(
    IN_TIME_24 = format(strptime(IN_TIME, format = "%I:%M %p"), format = "%H:%M"),
    OUT_TIME_24 = format(strptime(OUT_TIME, format = "%I:%M %p"), format = "%H:%M")
  )

#Now we should be able to combine the df's taking only the columns we need
columns <- c("username", "dataset", "Subject", "Night", "IN_DATE_TIME", "OUT_DATE_TIME", "IN_TIME", "OUT_TIME", "IN_TIME_24", "OUT_TIME_24", "SE", "TST", "TIB", "WASO")

#add a subject column and set equal to the values in username
preschooler_schoolage$Subject <- preschooler_schoolage$username

#rowbind the df's, taking only columns present in columns from each
sleep <- rbind(toddler[, columns], parent[, columns], preschooler_schoolage[, columns])

#Finally, as in the other datasets, convert O's to 0's for compatability
sleep$username <- sleep$username <- sub("O", "0", sleep$username)

# QAULITY CONTROL #

#cleaning env removing missing sleepd df's
rm(missing_date, missing_date_toddler, missing_date_parent, parent_same, same_day, toddler_same)

#**____________________________________________________________**#
## --------------- 4. COVARIATE DATA CLEANING ----------------- ##

#We want to ensure we take forward records with complete information in the required covariate columns, including age, sex, and BMI.

#First standardize age and BMI column of interest to single col name
adult <- adult %>%
  rename(age = "age_ha")

#We need to calculate age in years for the children
children <- children %>%
  mutate(age_years = round(as.numeric(difftime(date_ha, dob, units = "days")) / 365.25, 2))

children <- children %>%
  rename(age = "age_years") %>%
  rename(bmi = "bmi_z") #bmi_z is different but we are not comparing statistically here, just checking the data is present

#We can use a function to comprehensively check for missing data
check_missing <- function(column) {
  is.na(column) | column == "NA" | column == "" | column == "N/A (pregnant)" #count this last one for bmi specifically in adult dataset
}

#First check for missing values in age, sex or bmi columns now for each df. Since this may come in different forms, lets first check whats present in the col, and it's type
class(children$age)
unique(children$age)
#Take a quick look at a boxplot of age dist to ensure it is making sense
boxplot(children$age)

#now check the same for bmi, which should also be numeric but lets see whats in there first
class(children$bmi)
unique(children$bmi)
boxplot(children$bmi)

#and for sex we will just look at unique vakues
unique(children$sex) #should be M or F

#Now use the check_missing function on any of age, bmi or sex to look closer
#define our covariates 

covariates <- c("age", "sex", "bmi")

Missing_demos_child_df <- children %>%
  filter(apply(select(., all_of(covariates)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(covariates), everything()) #reorder to see missing data cols

#It looks like most of these are due to missing bmi, and often due to missing health assessments where no height is recorded. 

#Remove missing demo records here from children
children <- children %>%
  filter(!(pid %in% Missing_demos_child_df$pid))

rm(Missing_demos_child_df)

#Now repeat for the adults
class(adult$age)
unique(adult$age)
#Take a quick look at a boxplot of age dist to ensure it is making sense
boxplot(adult$age)

#and for bmi
class(adult$bmi)
#since it is character, we will want to print non numeric unique values

unique(adult$bmi[!grepl("^[0-9]", adult$bmi)])

#since we have "N/A (pregnant)" in there, we will want to note to remove records with this from adult

adult <- adult %>%
  filter(bmi != "N/A (pregnant)")

#now convert bmi to numeric
adult$bmi <- as.numeric(adult$bmi)

#Apply the check missing function to adult covariates
Missing_demos_adult_df <- adult %>%
  filter(apply(select(., all_of(covariates)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(covariates), everything()) #reorder to see missing data cols

#removing these records from adult
adult <- adult %>%
  filter(!(pid %in% Missing_demos_adult_df$pid))

#and check the boxplot
boxplot(adult$bmi)

#Since we have some high BMI values, we will want to consider removing these. Lets look at the relationship between kcal and bmi out of interest
plot(adult$kcal, adult$bmi)

#clean env
rm(Missing_demos_adult_df)

##**___________________________________________**##
## ------------- 5. SLEEP CLEANING ------------- ##

## SUBSET SLEEP DFs ##

#subset the sleep df to get preschooler and toddlers into a children sleep df, and the rest to the adult

#first rename username to pid
sleep <- sleep %>%
  rename(pid = "username")

children_sleep <- sleep %>%
  filter(dataset == "toddler" | dataset == "preschooler_schoolage")

adult_sleep <- sleep %>%
  filter(!(pid %in% children_sleep$pid))

## MISSING DATA ##

#First we will prepare the sleep data for analysis, generating summary level information with one entry per pid

#redefine our check_missing function for sleep cols now
check_missing <- function(column) {
  is.na(column) | column == "NA" | column == "" | column == 0 #count this last one for bmi specifically in adult dataset
}

sleep_cols <- c("TST", "TIB", "SE", "WASO")

#all of these should be numeric, so lets enforce this first.

#lapply as.numeric to all sleep cols in adult_sleep
adult_sleep[sleep_cols] <- lapply(adult_sleep[sleep_cols], as.numeric)

#check for missing data in sleep cols
missing_sleep_adult <- adult_sleep %>%
  filter(apply(select(., all_of(sleep_cols)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(sleep_cols), everything()) #reorder to see missing data cols

#Since all of the records in the missing_sleep adult has to do with 0's appearing in the WASO column - we can disregard this because that simply means their sleep efficiency (SE) was 100.

#Lets take a look at the values of TST, TIB SE and WASO

#Make a boxplot of each variable
boxplot(adult_sleep$TST, main="Boxplot of TST", ylab="TST (minutes)", col="lightblue")
boxplot(adult_sleep$TIB, main="Boxplot of TIB", ylab="TIB (minutes)", col="lightblue")
boxplot(adult_sleep$SE, main="Boxplot of SE", ylab="SE (%)", col="lightblue")
boxplot(adult_sleep$WASO, main="Boxplot of WASO", ylab="WASO (minutes)", col="lightblue")

#We can apply the same process to children_sleep now
children_sleep[sleep_cols] <- lapply(children_sleep[sleep_cols], as.numeric)
missing_sleep_children <- children_sleep %>%
  filter(apply(select(., all_of(sleep_cols)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(sleep_cols), everything()) #reorder to see missing data cols

#We have a ton of invalid numbers stemming from missing data in the TIB columns. While we can compute this using the SE and TST columns, we don't actually need it - our variable of interest is primarily TST and actually WASO, since TIB and SE can both be derived from TST and WASO, and WASO explains SE as well as the discrepancy between TST and TIB.

#Therefore, it turns out there is no need to remove data at this step. 

#Another check is to look for entries in sleep where the in DATE and OUT date are many days apart, indicating errors occuring at the device or data input levels.

#Now, combine our adult_sleep and children_sleep dataframes
sleep <- rbind(adult_sleep, children_sleep)

#create a filtered df to check for rows where IN_DATE_TIME is > 1 day from OUT_DATE_TIME
sleep_date_diff <- sleep %>%
  filter(as.Date(OUT_DATE_TIME) - as.Date(IN_DATE_TIME) > 1)

#These also appear as mistakes. Rather than correct them at this time, we will simply remove them from the dataset.

#Remove these specific rows from the sleep data. first add an identifier col to indicate rows to remove
sleep <- sleep %>%
  mutate(remove = ifelse(as.Date(OUT_DATE_TIME) - as.Date(IN_DATE_TIME) > 1, "Y", "N"))

#now remove rows with "Y" in the remove column
sleep <- sleep %>%
  filter(remove != "Y")

#conert IN_TIME to a time object
sleep$IN_TIME2 <- as.POSIXct(sleep$IN_TIME, format = "%I:%M %p")
class(sleep$IN_TIME2)

#Filter for rows where IN_TIME is betwen 5am and 2pm
sleep_am <- sleep %>%
  filter(hour(IN_TIME2) >= 5 & hour(IN_TIME2) <= 14)

#now remove pid's that are in the sleep_am from sleep. Here we don't just remove the rows, but we remove the participants entirely, since we are assuming that shift work will infleunce the reponse variables, and that this is not representative of the population we are interested in.
sleep <- sleep %>%
  filter(!(pid %in% sleep_am$pid))

#remove our temp columns
sleep <- sleep %>% select(-IN_TIME2, -remove)

#Write to csv in input_dir
setwd(Output_Dir)
write.csv(sleep, "sleep.csv", row.names = FALSE)

#Before we merge the sleep data with the diet/demographic and genotypes, we will need to summarize the data into mean statistics, since there are currently multiple records for a given participant. This will be done in the next script.

#**____________________________________________________________**#
## ------------------- 6. MERGING DATASETS -------------------- ##

## JOINING GENOTYPE AND DIET/DEMO DATA ##

## SUBSET GENOTYPE FILE ##

#Now we should ocmbine our sleep, genetic and demo/diet data for both the children and adult df's

#First, we should subset the genotype file into adult and parent df's. We will make a reference of pid's from the diet/demo df's.

adult_pids <- unique(adult$pid)
child_pids <- unique(children$pid)

#Now filter genotype file for childrne pid's only
children_genotype <- genotype_file %>%
  filter(pid %in% child_pids)

#and same for adults
adult_genotype <- genotype_file %>%
  filter(pid %in% adult_pids)

#First join demographic/diet data with genotype data
final_adult <- left_join(adult, adult_genotype, by = "pid")
final_child <- left_join(children, children_genotype, by = "pid")

# REMOVING MISSING DATA ##

#Now we will exclude participants at this point if they are missing any essential diet or genotype information.

#Make a vector of diet cols needed for computing response variables in the next script, or which are response variables themselves

diet_cols <- c("kcal", "prot", "carb", "tfat","g_refined", "add_sugr", "fibe")

#we can apply our check missing function over these columns

Missing_diet_adult <- final_adult %>%
  filter(apply(select(., all_of(diet_cols)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(diet_cols), everything()) #reorder to see missing data cols

Missing_diet_child <- final_child %>%
  filter(apply(select(., all_of(diet_cols)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(diet_cols), everything()) #reorder to see missing data cols

#Only a single record in the children was missing diet info. We will remove this pid now.

final_child <- final_child %>%
  filter(!(pid %in% Missing_diet_child$pid))

#Now we will check snp_columns, first identifying a list of them using the get_snp cols function from config.

#Starting with adult dataset
snp_columns <- get_snp_columns(final_adult)

Missing_genotype_adult <- final_adult %>%
  filter(apply(select(., all_of(snp_columns)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(snp_columns), everything()) #reorder to see missing data cols

#We will go ahead and remove these records from the dataset
final_adult <- final_adult %>%
  filter(!(pid %in% Missing_genotype_adult$pid))

#Now lets confirm the contents of the remaining genotype columns to see if any strange entries exist. From experience, I know we will come accross "--" as missing data as well.

#define a function to check for "--" to apply to all snp cols
check_dash <- function(column) {
   column == "--"
}

genotype_dash_adult <- final_adult %>%
  filter(apply(select(., all_of(snp_columns)), 1, function(row) any(check_dash(row)))) %>%
  select(pid, all_of(snp_columns), everything()) #reorder to see missing data cols

#Now we will remove all of these participants, which do not have a read for atleast one SNP

final_adult <- final_adult %>%
  filter(!(pid %in% genotype_dash_adult$pid))

#Repeating the above process for the children.

snp_columns <- get_snp_columns(final_child) #should be same as adult

Missing_genotype_child <- final_child %>%
  filter(apply(select(., all_of(snp_columns)), 1, function(row) any(check_missing(row)))) %>%
  select(pid, all_of(snp_columns), everything()) #reorder to see missing data cols

#removing these observations
final_child <- final_child %>%
  filter(!(pid %in% Missing_diet_child$pid))

#and now checking for "--"
genotype_dash_child <- final_child %>%
  filter(apply(select(., all_of(snp_columns)), 1, function(row) any(check_dash(row)))) %>%
  select(pid, all_of(snp_columns), everything()) #reorder to see missing data cols

#and removing these observations

final_child <- final_child %>%
  filter(!(pid %in% genotype_dash_child$pid))

#examine final structure
dim(final_adult)
dim (final_child)

#Now lets write these datasets to the output dir
set_dir(Output_Dir)

setwd(Output_Dir)

#and write to csv
write.csv(final_adult, "final_adult.csv", row.names = FALSE)
write.csv(final_child, "final_child.csv", row.names = FALSE)

