
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

genotype_file <- read.csv("filtered_genotype_file.csv") #Prepared from genotype file processing script

adult <- read.csv("adult_demos_diet.csv") #adult demographic/diet data, provided from Revital
children <- read.csv("children_demos_diet.csv") %>%
  filter(time_point == "t1")#Since there are duplicates of PID in this, we take rows where "t1" is in timepoint, since they contain the detailed data of interest

#We should ensure this is a fullproof way to remove duplicates by proving that the participant values remaining are in fact identical should we choose to filter for "t1", as opposed to those rows with blanks.
children2 <- read.csv("children_demos_diet.csv")%>%
  filter(is.na(time_point))#Since there are duplicates of PID in this, we take rows where "t1" is in timepoint, since they contain the detailed data of interest


FIDS <- unique(adult$fid)
FIDsc <- unique(children$fid)

combined_FIDs <- c(FIDS, FIDsc)
combined_FIDs <- as.data.frame(combined_FIDs) %>%
  distinct()

length(FIDsc)

length(FIDS)

one <- unique(children$pid)
two <- unique(children2$pid)
all.equal(one,two) #Should this print TRUE, we can proceed and remove children_dataset2

rm(children2,one,two)

#Force all participant IDs to have 0's in place of O's
adult$pid <- sub("O", "0", adult$pid)
children$pid <- sub("O", "0", children$pid)

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

## LOAD IN SLEEP DATA ##

#load in excel workbook sheets by their name from "raw_sleep_excel.xlsx"
preschooler_schoolage <- read_excel("raw_sleep_excel.xlsx", sheet = "preschooler_schoolage")
toddler <- read_excel("raw_sleep_excel.xlsx", sheet = "toddler")
parent <- read_excel("raw_sleep_excel.xlsx", sheet = "parent")

## QUALITY CONTROL ##

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
    adjusted_IN_TIME_minutes = if_else(as.Date(IN_DATE_TIME) == as.Date(OUT_DATE_TIME), 
                                       1440 + IN_TIME_minutes, 
                                       IN_TIME_minutes)
  )

#Compute mean minutes from midnight for each user
mean_times <- sleep %>%
  group_by(username) %>%
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
    SD_sleep_offset = sd(OUT_TIME_minutes, na.rm = TRUE),
    SD_TST = sd(TST, na.rm = TRUE),
    SD_SE = sd(SE, na.rm = TRUE)
  )


#Now create an overall regularity score by taking a mean of each participants min-max scaled sd across each column - subtracting from 1 to flip the relationship so that higher values correspond to more regular sleep patterns. 
regularity_data <- regularity_data %>%
  mutate(
    Regularity = 1 - rowMeans(cbind(
      (SD_sleep_onset - min(SD_sleep_onset, na.rm = TRUE)) /(max(SD_sleep_onset, na.rm = TRUE) - min(SD_sleep_onset, na.rm = TRUE)),
      (SD_sleep_offset - min(SD_sleep_offset, na.rm = TRUE)) / (max(SD_sleep_offset, na.rm = TRUE) - min(SD_sleep_offset, na.rm = TRUE)),
      (SD_TST - min(SD_TST, na.rm = TRUE)) / (max(SD_TST, na.rm = TRUE) - min(SD_TST, na.rm = TRUE)),
      (SD_SE - min(SD_SE, na.rm = TRUE)) / (max(SD_SE, na.rm = TRUE) - min(SD_SE, na.rm = TRUE))), na.rm = TRUE)
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

#Lets generate some plots to check that the data is looking right

#Boxplot of kcal from final_child df, making y axis numebrs larger and the box fill teal

par(cex.axis=2)
#Larger title font as well, 2
#teal fill by hex number with lower opcaicty
boxplot(final_child$kcal, col = "#00808033", ylab = "Kcal", main = "Kcal Distribution in Children", cex.main = 2)

library(plotly)

create_boxplot <- function(df, colname) {
  # Ensure the column is numeric
  df <- df %>%
    mutate(!!sym(colname) := as.numeric(!!sym(colname)))
  
  fig <- plot_ly(
    data = df, 
    y = ~ get(colname), 
    type = "box", 
    boxpoints = "all", 
    jitter = 0.3, 
    pointpos = -1.8
  )

  fig <- fig %>% layout(title = paste("Distribution of", colname))
  #y axis label name 
  fig <- fig %>% layout(yaxis = list(title = colname))

  return(fig)
}

#Write final datasets to CSV in folder for next step, Dataset cleaning
setwd(Output_Dir)

#combine diet_repsonse_columns and sleep_response_columns lists for one main list of response_variables
response_variables <- c(diet_response_columns, sleep_response_columns)

#Now, in the final_child, loop over each response variable and create a boxplot, then save image in output dir under boxplots

#First create boxplots folder in current dir, if it doesnt exist
if (!dir.exists("boxplots")) {
  dir.create("boxplots")
}

setwd("boxplots")

loops <- '
for (response_variable in response_variables) {
  fig <- create_boxplot(final_child, response_variable)

  # Save to boxplots dir as png
  png_filename <- paste0(response_variable, "_boxplot.png")

  fig <- fig %>% config(displayModeBar = FALSE)
  fig <- fig %>% config(displaylogo = FALSE)
  
  # Save the figure as PNG
  kaleido(fig, file = png_filename)
}
'

setwd(Output_Dir)

write.csv(final_adult, "final_adult.csv",  row.names = FALSE)
write.csv(final_child, "final_child.csv", row.names = FALSE)

#Clean Env
rm(final_adult,final_child,genotype_file,sleep,sleep_summary,adult,children,sleep_reasons,Supp_Reasons,unique_values_adult,supps, macros, multipliers, regularity_data, mean_times)
