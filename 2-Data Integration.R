
#** -------- DATA INTEGRATION --------- **#
# ======================================= #
#    Zach Ribau  |  June 18, 2024       

#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##

source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

Input_Dir <- set_dir(Data_Integration)
Output_Dir <- set_dir(Dataset_Cleaning)

#**_____________________________________________________**#
## ----------------------- DATA ------------------------ ##

setwd(Input_Dir)

#Read in required data files
sleep <- read.csv("sleep.csv") #Provided from Lori
genotype_file <- read.csv("filtered_genotype_file.csv") #Prepared from genotype file processing script

adult <- read.csv("adult_demos_diet.csv") #adult demographic/diet data, provided from Revital
children <- read.csv("children_demos_diet.csv") %>%
  filter(time_point == "t1")#Since there are duplicates of PID in this, we take rows where "t1" is in timepoint, since they contain the detailed data of interest

#We should ensure this is a fullproof way to remove duplicates by proving that the participant values remaining are in fact identical should we choose to filter for "t1", as opposed to those rows with blanks.
children2 <- read.csv("children_demos_diet.csv")%>%
  filter(is.na(time_point))#Since there are duplicates of PID in this, we take rows where "t1" is in timepoint, since they contain the detailed data of interest

one <- unique(children$pid)
two <- unique(children2$pid)
all.equal(one,two) #Should this print TRUE, we can proceed and remove children_dataset2

rm(children2,one,two)

#Force all participant IDs to have 0's in place of O's
adult$pid <- sub("O", "0", adult$pid)
children$pid <- sub("O", "0", children$pid)
sleep$username <- sleep$username <- sub("O", "0", sleep$username)

#**________________________________________________________**#
## ----------------- 1. GENOTYPE DATA PREP ---------------- ##

#Transpose final genotype file
genotype_file <- as.data.frame(t(genotype_file))
colnames(genotype_file) <- genotype_file[1,] #Set SNP ids as column names
genotype_file <- genotype_file[-1,] #Remove SNP id column now
genotype_file <- genotype_file %>% rownames_to_column(var = "pid")
colnames(genotype_file) <- genotype_file[1, ]
genotype_file <- genotype_file[-1, ] %>%
  rename(pid = "SNPS")

genotype_file$pid <- sub("O", "0", genotype_file$pid)

#**____________________________________________________**#
## ----------------- 2. SLEEP DATA PREP  -------------- ##

#Now we will process the time columns to support operations such as taking the mean, to express tendency for later versus early wake/sleep times.

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
    adjusted_IN_TIME_minutes = if_else(hour(IN_DATETIME) < 12, 1440 - IN_TIME_minutes, IN_TIME_minutes)
  )

sleep <- sleep %>%
  mutate(
    adjusted_OUT_TIME_minutes = if_else(hour(OUT_DATETIME) < 12, 1440 - OUT_TIME_minutes, OUT_TIME_minutes)
  )

#Compute mean minutes from midnight for each user
mean_times <- sleep %>%
  group_by(username) %>%
  summarise(
    mean_IN_TIME_minutes = mean(IN_TIME_minutes, na.rm = TRUE),
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
  select(username, mean_IN_TIME_24, mean_OUT_TIME_24)

#Now summarize with  means of sleep columns, as well as total wear nights
sleep_summary <- sleep %>%
  group_by(username) %>%
  summarise(
    Nights = n(),
    Mean_TST = mean(TST, na.rm = TRUE),
    Mean_SE = mean(SE, na.rm = TRUE),
    Mean_WASO = mean(WASO, na.rm = TRUE),
    Mean_TIB = mean(TIB, na.rm = TRUE)
  )

#We can now incorporate a measure of sleep regularity via an aggregation of variability for duration, sleep onset, sleep offset, and sleep midpoint, suing standard deviation as per Fischer et al 2021: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8503839/

#Define which columns
regularity_data <- sleep %>%
  group_by(username) %>%
  summarise(
    SD_sleep_onset = sd(adjusted_IN_TIME_minutes, na.rm = TRUE),
    SD_sleep_offset = sd(adjusted_OUT_TIME_minutes, na.rm = TRUE),
    SD_TST = sd(TST, na.rm = TRUE),
    SD_SE = sd(SE, na.rm = TRUE)
  )

#Now create an overall regularity score by taking a mean of each participants min-max scaled sd across each column - subtracting from 1 to flip the relationship so that higher values correspond to more regular sleep patterns. 
regularity_data <- regularity_data %>%
  mutate(
    Regularity = 1 - rowMeans(cbind(
      (SD_sleep_onset - min(SD_sleep_onset, na.rm = TRUE)) / 
        (max(SD_sleep_onset, na.rm = TRUE) - min(SD_sleep_onset, na.rm = TRUE)),
      (SD_sleep_offset - min(SD_sleep_offset, na.rm = TRUE)) / 
        (max(SD_sleep_offset, na.rm = TRUE) - min(SD_sleep_offset, na.rm = TRUE)),
      (SD_TST - min(SD_TST, na.rm = TRUE)) / 
        (max(SD_TST, na.rm = TRUE) - min(SD_TST, na.rm = TRUE)),
      (SD_SE - min(SD_SE, na.rm = TRUE)) / 
        (max(SD_SE, na.rm = TRUE) - min(SD_SE, na.rm = TRUE))
    ), na.rm = TRUE)
  )

#Combine all sleep data
sleep_summary <- left_join(sleep_summary, mean_times, by = "username")
sleep_summary <- left_join(sleep_summary, regularity_data, by = "username")

#Rename username column for later merge
sleep_summary <- sleep_summary %>%
  rename(pid = username)

#**__________________________________________________________**#
## ---------- 3. DEMOGRAPHIC AND DIETARY DATA PREP ---------- ##

unique_values_adult <- unique(unlist(adult[, 61:70]))

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


adult <- adult %>%
  rowwise() %>%
  mutate(Supplement_Affecting_Sleep = ifelse(any(c_across(61:70) %in% supps), "Y", "N"))

Supp_Reasons <- unique(unlist(adult[, 71:81]))

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


adult <- adult %>%
  rowwise() %>%
  mutate(Reason_For_Sleep = ifelse(any(c_across(71:81) %in% sleep_reasons), "Y", "N"))

#Finally we make a master column to flag down participants taking sleep supplements

adult <- adult %>%
  mutate(Sleep_Supplement = ifelse(Reason_For_Sleep == "Y" | Supplement_Affecting_Sleep == "Y", "Y", "N"))

#Create function to make columns specific for conditions of interest and put "Y" if they are reported
create_condition_columns <- function(data, column_name, conditions) {
  for (condition in conditions) {
    # Create a valid column name by replacing spaces with underscores and removing special characters
    condition_col_name <- gsub("[^[:alnum:]_]", "", gsub(" ", "_", condition))
    
    # Use mutate to create the new column
    data <- data %>%
      mutate(!!condition_col_name := ifelse(grepl(paste0("(?i)", condition), .data[[column_name]]), "Y", "N"))
  }
  return(data)
}

conditions <- c("Depression", "Anorexia nervosa", "Bulimia nervosa", "Binge eating disorder",
                "Irritable bowel syndrome", "Ulcerative colitis/Crohn's disease", "Anxiety Disorder")

adult <- create_condition_columns(adult, "illnesses_selected_choice", conditions) #Apply function to list of disorders of interest

#Create calculation for child age
children <- children %>%
  mutate(dob = as.Date(dob, format="%Y-%m-%d"),
         date_ha = as.Date(date_ha, format="%Y-%m-%d"))

children <- children %>%
  mutate(age_years = round(as.numeric(difftime(date_ha, dob, units = "days")) / 365.25, 2))

adult <- adult %>%
  select(-Reason_For_Sleep, -Supplement_Affecting_Sleep)

#Now prepare macro nutrient proportion columns, starting with a function we can recycle.

#Now apply to macronutrient columns for both datastes
macros <- c("prot", "carb", "tfat")
multipliers <- c(prot = 4, carb = 4, tfat = 9)

#Use a loop since it's a small vector
for (macro in macros) {
  proportion_column <- paste0(macro, "_prop")
  adult[[proportion_column]] <- adult[[macro]] * multipliers[[macro]] / adult[["kcal"]]
  
  # Reorder the columns to place the new column next to the original column
  adult <- adult %>%
    select(1:all_of("caff"), all_of(proportion_column), everything())
}

#and for childen
for (macro in macros) {
  proportion_column <- paste0(macro, "_prop")
  children[[proportion_column]] <- children[[macro]] * multipliers[[macro]] / children[["kcal"]]
  
  # Reorder the columns to place the new column next to the original column
  children <- children %>%
    select(1:all_of("caff"), all_of(proportion_column), everything())
}

#**____________________________________________**#
## ------------ 4. MERGING DATA --------------- ##

#First join demographic/diet data with sleep data
final_adult <- left_join(adult, sleep_summary, by = "pid")
final_child <- left_join(children, sleep_summary, by = "pid")

#Now merge resulting dataframes with SNP data
final_adult <- left_join(final_adult, genotype_file, by = "pid")
final_child <- left_join(final_child, genotype_file, by = "pid")

#Write final datasets to CSV in folder for next step, Dataset cleaning
setwd(Output_Dir)

write.csv(final_adult, "final_adult.csv",  row.names = FALSE)
write.csv(final_child, "final_child.csv", row.names = FALSE)

#Clean Env
rm(final_adult,final_child,genotype_file,sleep,sleep_summary,adult,children,sleep_reasons,Supp_Reasons,unique_values_adult,supps, macros, multipliers, regularity_data, mean_times)
