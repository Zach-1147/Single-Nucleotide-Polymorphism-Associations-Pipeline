#** ------ CLEANING & RESPONSE VARIABLES ------- **#
# =============================================== #
#    Zach Ribau  |  June 18, 2024    
rm(list=ls())
#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

Input_Dir <- set_dir(Cleaning_Response_Variables)
Output_Dir <- set_dir(Criteria_Based_Filtering)

setwd(Input_Dir)

#**_____________________________________________________**#
## --------------------- 1. DATA ------------------------ ##

#Lets start by loading in the diet/demo and genotype df's
adult <- read.csv("final_adult.csv")
children <- read.csv("final_child.csv")

#and reading in the adult and child sleep csv's
sleep <- read.csv("sleep.csv")

##**_____________________________________________________**##
## ------------- 2. SLEEP RESPONSE VARIABLES ------------- ##

## DERIVING MEAN SLEEP AND WAKE TIMES ##

#Becuase you cannot perform mathematical operations on time data right out of the gate, we will need to create columns to represent sleep and wake times as numeric by computing relative times from midnight.

#Convert the time columns to representations of time since midnight
sleep <- sleep %>%
  mutate(
    IN_DATETIME = as.POSIXct(paste(IN_DATE_TIME, IN_TIME), format="%Y-%m-%d %I:%M %p"),
    OUT_DATETIME = as.POSIXct(paste(OUT_DATE_TIME, OUT_TIME), format="%Y-%m-%d %I:%M %p"),
    IN_TIME_minutes = hour(IN_DATETIME) * 60 + minute(IN_DATETIME), #calculate minutes from midnight
    OUT_TIME_minutes = hour(OUT_DATETIME) * 60 + minute(OUT_DATETIME)
  )

#Adjust for relativity around minute for later calc of standard dev
sleep <- sleep %>%
  mutate(
    adjusted_IN_TIME_minutes = if_else(as.Date(IN_DATE_TIME) == as.Date(OUT_DATE_TIME), 
                                       1440 + IN_TIME_minutes, 
                                       IN_TIME_minutes)
  )

#Compute mean minutes from midnight for each user
mean_times <- sleep %>%
  group_by(pid) %>%
  summarise(
    mean_IN_TIME_minutes = mean(adjusted_IN_TIME_minutes, na.rm = TRUE),
    mean_OUT_TIME_minutes = mean(OUT_TIME_minutes, na.rm = TRUE)
  )

#Reconstruct hour and minute time into 24-hour time format
mean_times <- mean_times %>%
  mutate(
    mean_IN_HOUR = as.integer(floor(mean_IN_TIME_minutes / 60)),
    mean_IN_MINUTE = as.integer(round(mean_IN_TIME_minutes %% 60)),
    mean_OUT_HOUR = as.integer(floor(mean_OUT_TIME_minutes / 60)),
    mean_OUT_MINUTE = as.integer(round(mean_OUT_TIME_minutes %% 60)),
    mean_IN_TIME_24 = sprintf("%02d:%02d", mean_IN_HOUR, mean_IN_MINUTE),
    mean_OUT_TIME_24 = sprintf("%02d:%02d", mean_OUT_HOUR, mean_OUT_MINUTE)
  ) %>%
  select(pid, mean_IN_TIME_24, mean_OUT_TIME_24, mean_IN_TIME_minutes, mean_OUT_TIME_minutes)

#Now summarize with  means of sleep columns, as well as total wear nights
sleep_summary <- sleep %>%
  group_by(pid) %>%
  summarise(
    Nights = n(),
    Mean_TST = mean(TST, na.rm = TRUE),
    Mean_SE = mean(SE, na.rm = TRUE),
    Mean_WASO = mean(WASO, na.rm = TRUE), 
    dataset = first(dataset)
  )

## SPLITTING INTO ADULT AND CHILDREN ##
adult_sleep_summary <- sleep_summary %>%
  filter(dataset == "parent")

unique(sleep$dataset)

#and child
child_sleep_summary <- sleep_summary %>%
  filter(dataset == "toddler" | dataset == "preschooler_schoolage")

#Finally, we will calculate sleep regularity within each subset.

## ADULT SLEEP REGULARITY ##

sleep_adult <- sleep %>%
  filter(dataset == "parent")

#Define which columns
adult_regularity <- sleep_adult %>%
  group_by(pid) %>%
  summarise(
    SD_sleep_onset = sd(adjusted_IN_TIME_minutes, na.rm = TRUE),
    SD_sleep_offset = sd(OUT_TIME_minutes, na.rm = TRUE),
    SD_TST = sd(TST, na.rm = TRUE),
    SD_SE = sd(SE, na.rm = TRUE)
  )

#Now create an overall regularity score by taking a mean of each participants min-max scaled sd across each column - subtracting from 1 to flip the relationship so that higher values correspond to more regular sleep patterns. 
adult_regularity <- adult_regularity %>%
  mutate(
    Regularity = 1 - rowMeans(cbind(
      (SD_sleep_onset - min(SD_sleep_onset, na.rm = TRUE)) /(max(SD_sleep_onset, na.rm = TRUE) - min(SD_sleep_onset, na.rm = TRUE)),
      (SD_sleep_offset - min(SD_sleep_offset, na.rm = TRUE)) / (max(SD_sleep_offset, na.rm = TRUE) - min(SD_sleep_offset, na.rm = TRUE)),
      (SD_TST - min(SD_TST, na.rm = TRUE)) / (max(SD_TST, na.rm = TRUE) - min(SD_TST, na.rm = TRUE)),
      (SD_SE - min(SD_SE, na.rm = TRUE)) / (max(SD_SE, na.rm = TRUE) - min(SD_SE, na.rm = TRUE))), na.rm = TRUE)
  )

## CHILDREN SLEEP REGULARITY ##

sleep_children <- sleep %>%
  filter(dataset == "toddler" | dataset == "preschooler_schoolage")

child_regularity <- sleep_children %>%
  group_by(pid) %>%
  summarise(
    SD_sleep_onset = sd(adjusted_IN_TIME_minutes, na.rm = TRUE),
    SD_sleep_offset = sd(OUT_TIME_minutes, na.rm = TRUE),
    SD_TST = sd(TST, na.rm = TRUE),
    SD_SE = sd(SE, na.rm = TRUE)
  )

child_regularity <- child_regularity %>%
  mutate(
    Regularity = 1 - rowMeans(cbind(
      (SD_sleep_onset - min(SD_sleep_onset, na.rm = TRUE)) /(max(SD_sleep_onset, na.rm = TRUE) - min(SD_sleep_onset, na.rm = TRUE)),
      (SD_sleep_offset - min(SD_sleep_offset, na.rm = TRUE)) / (max(SD_sleep_offset, na.rm = TRUE) - min(SD_sleep_offset, na.rm = TRUE)),
      (SD_TST - min(SD_TST, na.rm = TRUE)) / (max(SD_TST, na.rm = TRUE) - min(SD_TST, na.rm = TRUE)),
      (SD_SE - min(SD_SE, na.rm = TRUE)) / (max(SD_SE, na.rm = TRUE) - min(SD_SE, na.rm = TRUE))), na.rm = TRUE)
  )

#Now lets merge all regularity data into the sleep summary dataframes

adult_sleep_summary <- merge(adult_sleep_summary, adult_regularity, by = "pid")

child_sleep_summary <- merge(child_sleep_summary, child_regularity, by = "pid")

#finally, we can merge the mean times into the sleep summary dataframes

adult_sleep_summary <- merge(adult_sleep_summary, mean_times, by = "pid")

child_sleep_summary <- merge(child_sleep_summary, mean_times, by = "pid")

#**_____________________________________________________**#
## ----------------- 3. DIET VARIABLES ------------------ ##

#We will now compute relative energy derived from each macronutrient as our dietary response variables, along with kcal

#define macro columns to  use. For interest, we will also look at energy derived from sugar and refined grains as well. 
macros <- c("prot", "carb", "tfat", "sugr")

#specify the multipliers for each macronutrient, according to well established values often used for research
multipliers <- c(prot = 4, carb = 4, tfat = 9, sugr = 4)

#Use a loop to calculate since it's a small vector
for (macro in macros) {
  proportion_column <- paste0(macro, "_prop")

  #Here we take the original column with grams, multiply by the multiplier to get kcal, and divide by total kcal to get the proportion of energy derived from that nutrient of interest
  adult[[proportion_column]] <- adult[[macro]] * multipliers[[macro]] / adult[["kcal"]]
  
  #Reorder the columns to place the new column next to the original column
  adult <- adult %>%
    select(1:all_of("caff"), all_of(proportion_column), everything())
}

#and for childen
for (macro in macros) {
  proportion_column <- paste0(macro, "_prop")
  children[[proportion_column]] <- children[[macro]] * multipliers[[macro]] / children[["kcal"]]
  
  #Reorder the columns to place the new column next to the original column
  children <- children %>%
    select(1:all_of("caff"), all_of(proportion_column), everything())
}

#clean up env

rm(macros, macro, multipliers, proportion_column)

#**_________________________________________________________*
*#
# -------- 6. DEMOGRAPHIC DATA EXPLORATION & SUMMARY -------- #

# NATURAL HEALTH PRODUCTS PERTAINING TO SLEEP #

#Now lets make a list of all columns in both adult and child starting with "natural_product"

all_natural_products <- grep("^natural_product", names(adult), value = TRUE)

#On investgation, the first 10 are supplements themsvels, and the rest are reasons for taking them. Lets split these.

#select first 10 from adult_natural_product_columns
nat_prod_cols <- all_natural_products[1:10]

#Now get all unique values occuring within these
natural_products <- unique(unlist(adult[nat_prod_cols]))

#On investigation, we look for supplements that may affect sleep, saving to a list below
supps <- c("Magnesium 500 mg",
           "Magnesium",
           "Magnesium gylcinate",
           "melatonin",
           "Melatonin",
           "Cannabis",
           "5-HTP",
           "Calcium and magnesium",
           "Magnesium supplement",
           "Calcium/Magnesium 500mg",
           "Magnesium Citrate Powder 1 tsp daily",
           "Magnesium Citrate. 300 mg (2-4 times/week)",
           "1 pill melatonin",
           "200mg Magnesium",
           "Magnesium 2 pills daily",
           "Magnesium Bis-Glycinate",
           "Magnesium 200 mg",
           "Melatonin, 3mg")

#add columns to indicate Y or N for presence of taking any of the above
adult <- adult %>%
  rowwise() %>%
  mutate(Supplement_Affecting_Sleep = ifelse(any(c_across(all_of(nat_prod_cols)) %in% supps), "Y", "N"))

#Now access the remaining adult_natrual_product_columns NOT in adult_nat_prods
nat_prod_reasons <- setdiff(all_natural_products, nat_prod_cols)

#and now a unique list from nat_prod_reasons
natural_product_reasons <- unique(unlist(adult[nat_prod_reasons]))

#On investigation, make a list of reasons pertaining to sleep specifically
sleep_reasons <- c("More restful sleep",
                   "help with sleep",
                   "Naturopath recommendation for sleep/muscle rest",
                   "for sleep as needed",
                   "1 capsule of 3 mg daily, to help sleep",
                   "Sleep",
                   "I want to fall asleep fast, sleep is precious to me",
                   "Help sleeping",
                   "Mood and sleep",
                   "Helping with growing pains and sleep",
                   "Sleep/stress; general health",
                   "Doctor prescribed for sleep",
                   "Sleep, anxiety, migraines, sickness etc- use when required",
                   "2-5 times per week, general health for sleep and headaches",
                   "Sleep/restless legs",
                   "Taken approximately once every 1 to 2 weeks to get a good night sleep")

#and adding Y/N cols
adult <- adult %>%
  rowwise() %>%
  mutate(Reason_For_Sleep = ifelse(any(c_across(all_of(nat_prod_reasons)) %in% sleep_reasons), "Y", "N"))

#Finally we make a master column to flag down participants taking sleep supplements

adult <- adult %>%
  mutate(Sleep_Supplement = ifelse(Reason_For_Sleep == "Y" | Supplement_Affecting_Sleep == "Y", "Y", "N"))

#now combine all nat_prod_col values into a single nat_products column, concatonating values with a semi colon

adult <- adult %>%
  rowwise() %>%
  mutate(natural_products = paste(na.omit(c_across(all_of(nat_prod_cols))[c_across(all_of(nat_prod_cols)) != ""]), collapse = "; "))

#repeat this process for the nat_prod_reasons
adult <- adult %>%
  rowwise() %>%
  mutate(natural_products_reasons = paste(na.omit(c_across(all_of(nat_prod_reasons))[c_across(all_of(nat_prod_reasons)) != ""]), collapse = "; "))

#and now we can remove the original columns
adult <- adult %>%
  select(-all_of(all_natural_products))

#Clean up env

rm(all_natural_products, nat_prod_reasons, nat_prod_cols, natural_products, sleep_reasons, supps)

# ILLNESSES AND CONDITIONS #

#Create function to make columns specific for conditions of interest and put "Y" if they are reported
create_condition_columns <- function(data, column_name, conditions) {
  for (condition in conditions) {
    #Create a valid column name by replacing spaces with underscores and removing special characters
    condition_col_name <- gsub("[^[:alnum:]_]", "", gsub(" ", "_", condition))
    
    #Use mutate to create the new column
    data <- data %>%
      mutate(!!condition_col_name := ifelse(grepl(paste0("(?i)", condition), .data[[column_name]]), "Y", "N"))
  }
  return(data)
}

#Now make a list of conditions from investigation of  columns in the df
conditions <- c("Depression", "Anorexia nervosa", "Bulimia nervosa", "Binge eating disorder",
                "Irritable bowel syndrome", "Ulcerative colitis or Crohn's disease", "Anxiety Disorder")

adult <- create_condition_columns(adult, "illnesses_selected_choice", conditions) #Apply function to list of disorders of interest

#clean env

rm(conditions)

# COLUMN CLEANING #

#first clean up the env

#Investigate column names
snp_cols <- get_snp_columns(adult)
ad_col_names <- as.data.frame(setdiff(colnames(adult), snp_cols))

cols_to_rm <- c("time_point", "time_point_survey", "parent_in_study", "experimental_group" ,"age_survey", "date_survey", "location_ha", "i_pregnancy_due_date_ha", "cleaned_ha", "phase_survey", "dob", "illnesses_selected_choice", "illnesses_other_major_illness_since_2010_please_specify_text")

#remove these columns
adult <- adult %>%
  select(-all_of(cols_to_rm))

#also remove all cols starting with ht or wc
wc_ht <- grep("^ht|^wc", names(adult), value = TRUE)

adult <- adult %>%
  select(-all_of(wc_ht))

## INVESTIGATE CHILDRENS DEMO COLUMNS HERE!! ##

#**_____________________________________________________**#
## ----------------- 7. MERGING DATA ------------------- ##

#merge adult_sleep_summary with adult df
adult <- merge(adult, adult_sleep_summary, by = "pid")
adult <- adult %>%
  select(-dataset)

#make list of sleep cols and diet cols, only those we are modelling as response variables, to appear in the begining for ease of viewing
sleep_cols <- c("Mean_TST", "Mean_SE", "Mean_WASO", "Regularity", "mean_IN_TIME_24", "mean_OUT_TIME_24")

diet_cols <- c("kcal", "prot_prop", "carb_prop", "tfat_prop", "sugr_prop", "caff")

#Now order the df to begin with pid, age, sex, and then sleep and diet cols, then snp_cols
adult <- adult %>%
  select(pid, age, sex, all_of(sleep_cols), all_of(diet_cols), all_of(snp_cols), everything())

#apply the same to child, starting with merge with sleep_summary
children <- merge(children, child_sleep_summary, by = "pid")

#and re-ordering
children <- children %>%
  select(pid, age, sex, all_of(sleep_cols), all_of(diet_cols), all_of(snp_cols), everything())

setwd(Output_Dir)

#write adult and children to csv as final_adult and final_child
write.csv(adult, "final_adult.csv", row.names = FALSE)
write.csv(children, "final_child.csv", row.names = FALSE)