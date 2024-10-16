## Association Testing Attempt 2

#** -------- Association Testing -------- **#
# ========================================= #
#    Zach Ribau  |  June 27, 2024       

#**_____________________________________________________**#
## -------------------- CONFIG ------------------------- ##
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Config", "config.R"))

#No specific outputs for this just yet
Output_Dir <- NA

#**_______________________________________________**#
## ---------- 1. DATA AND FUNCTIONS  ------------- ##

#Access data from database object
adult_data <- Database[["Complete Datasets"]][["Adults"]]
children_data <- Database[["Complete Datasets"]][["Children"]]
full_cohort_data <- Database[["Complete Datasets"]][["Full Cohort"]]

#Re-evaluate SNP columns using function from config script
snp_columns <- get_snp_columns(full_cohort_data)

full_cohort_data[snp_columns] <- lapply(full_cohort_data[snp_columns], as.numeric)
adult_data[snp_columns] <- lapply(adult_data[snp_columns], as.numeric)
children_data[snp_columns] <- lapply(children_data[snp_columns], as.numeric)

#Function to fit General Linear Models (GLMs) in Adult subset
fit_glms <- function(response, snp_columns, covariates, data) {
  
  #Formula construction
  snp_formula <- paste(snp_columns, collapse = " + ")
  response_formula <- paste(response, "~")
  covariates_str <- paste(covariates, collapse = " + ")
  
  #Combine components into formula that will execute model fitting when used in below
  full_formula <- as.formula(paste(response_formula, snp_formula, "+", covariates_str))
  print(full_formula)
  
  #Model fitting call using Rs glm
  model <- glm(
    formula = full_formula,
    data = data,
    family = "gaussian" 
  )
  
  return(summary(model)) #return the model summary
}

#Function to fit General Estimating Equations (GEEs) in childrens and full cohort subsets
fit_gees <- function(response, snp_columns, covariates, data, id_column, correlation_structure = "exchangeable") {
  
  #Formula construction
  snp_formula <- paste(snp_columns, collapse = " + ")
  response_formula <- paste(response, "~")
  covariates_str <- paste(covariates, collapse = " + ")
  
  
  full_formula <- as.formula(paste(response_formula, snp_formula, "+", covariates_str))
  
  #Model fitting
  gee_model <- geeglm(
    formula = full_formula,
    data = data,
    id = data[[id_column]],  #assign the cluster variable, (which is fid for family id in dataset)
    family = gaussian(),  
    corstr = correlation_structure #assign from global variable in master controller script
  )
  
  return(summary(gee_model))
}


#Function to obtain data from model summaries built with functions above: we can access model variables and in particular their p_values for all of the models built with the functions defined above.

model_dataframe <- function(model_summary, trait_category, trait, cohort, p_col) {
  coefficients_df <- as.data.frame(model_summary$coefficients)
  coefficients_df <- coefficients_df %>%
    mutate(term = rownames(coefficients_df)) %>%
    select(term, everything()) %>%
    filter(term != "(Intercept)") %>%
    filter(!(term %in% c("age", "bmi", "sexM", "sexF"))) %>%     #here we will remove covariates (Sep 16 2024 addition)
    
    #Compute adjusted p value using pval_adj_method variable which is defined in master controller script
    mutate(adjusted_p = p.adjust(as.numeric(.data[[p_col]]), method = pval_adj_method)) %>%
    arrange(adjusted_p) %>%
    
    #Reference our threshold p value parameter (master controller) to decide if significant (TRUE)
    mutate(significant = case_when(adjusted_p < p_value_threshold ~ "TRUE", TRUE ~ "FALSE")) %>%
    mutate(trait_category = trait_category) %>%
    mutate(trait = trait) %>%
    mutate(cohort = cohort)
  
  rownames(coefficients_df) <- NULL
  
  return(coefficients_df)
}

#Now we proceed with model fitting on the different data subsets

#**__________________________________________________________**#
## ---------- 2. ADULT - SNP ASSOCIATION MODELS ------------- ##

## ----------------- 2.1 Diet Phenotype Models ----------------- ##

#First we will apply our function to fit a GLM for each trait of interest within the diet response columns

kcal_adult_model <- fit_glms("kcal", snp_columns, covariates, adult_data)

kcal_adult_model


#initialize list
diet_model_dfs <- list()

#Applying the function to generate dataframes of each trait model summary and save to list
for(trait in names(diet_models_adult)) {
  
  model_summary <- diet_models_adult[[trait]]
  diet_model_dfs[[trait]] <- model_dataframe(model_summary, "diet", trait, cohort = 'adults', p_col = "Pr(>|t|)")
}

#Now we can bind all observations into a single dataframe
adult_diet_coeffs <- do.call(rbind, diet_model_dfs) %>%
  arrange(adjusted_p)
rownames(adult_diet_coeffs) <- NULL

## ----------------- 2.2 Sleep Phenotype Models ----------------- ##

#Repeating above process for sleep variables

#build modles for sleep responses
sleep_models_adult <- setNames(lapply(sleep_response_columns, function(response) {
  fit_glms(response, snp_columns, covariates, adult_data)
}), sleep_response_columns) 

sleep_model_dfs <- list()

#apply model dataframe function over traits in loop
for(trait in names(sleep_models_adult)) {
  model_summary <- sleep_models_adult[[trait]]
  sleep_model_dfs[[trait]] <- model_dataframe(model_summary, "sleep", trait, cohort = 'adults', p_col = "Pr(>|t|)")
}

#Combining
adult_sleep_coeffs <- do.call(rbind, sleep_model_dfs) %>%
  arrange(adjusted_p)
rownames(adult_sleep_coeffs) <- NULL

#**_____________________________________________________________**#
## ---------- 3. CHILDREN - SNP ASSOCIATION MODELS ------------- ##

## ----------------- 3.1 Diet Phenotype Models ----------------- ##

#Repeating steps from above, but incorporating fid column
children_data <- children_data %>%
  arrange(fid) #arrange to satisfy requirements of gee package

#build modles for traits
diet_models_children <- setNames(lapply(diet_response_columns, function(response) {
  fit_gees(response, snp_columns, covariates, children_data, id_column = "fid", correlation_structure)
}), diet_response_columns)

#initiate list
diet_model_dfs <- list()

#apply model dataframe function over traits in loop
for(trait in names(diet_models_children)) {
  model_summary <- diet_models_children[[trait]]
  diet_model_dfs[[trait]] <- model_dataframe(model_summary, "diet", trait, cohort = 'children', p_col = "Pr(>|W|)")
}

#combining into df
children_diet_coeffs <- do.call(rbind, diet_model_dfs) %>%
  arrange(adjusted_p)
rownames(children_diet_coeffs) <- NULL

## ----------------- 3.2 Sleep Phenotype Models ----------------- ##

sleep_models_children <- setNames(lapply(sleep_response_columns, function(response) {
  fit_gees(response, snp_columns, covariates, children_data, id_column = "fid", correlation_structure)
}), sleep_response_columns)

sleep_model_dfs <- list()

for(trait in names(sleep_models_children)) {
  model_summary <- sleep_models_children[[trait]]
  sleep_model_dfs[[trait]] <- model_dataframe(model_summary, "sleep", trait, cohort = 'children', p_col = "Pr(>|W|)")
}

#combining
children_sleep_coeffs <- do.call(rbind, sleep_model_dfs) %>%
  arrange(adjusted_p)
rownames(children_sleep_coeffs) <- NULL

#**________________________________________________________________**#
## ---------- 4. FULL COHORT - SNP ASSOCIATION MODELS ------------- ##

## ----------------- 4.1 Diet Phenotype Models ----------------- ##

#Repeating above steps, but removing bmi column since we cannot compare in adult/children subsets

full_cohort_data <- full_cohort_data %>%
  arrange(fid) #arrange on family id

covariates = c("age", "sex") #re-assign covariates (initially assigned in master controller script)

#Build models
diet_models_full <- setNames(lapply(diet_response_columns, function(response) {
  fit_gees(response, snp_columns, covariates, full_cohort_data, id_column = "fid", correlation_structure)
}), diet_response_columns)


#save summaries to list
for(trait in names(diet_models_full)) {
  model_summary <- diet_models_full[[trait]]
  diet_model_dfs[[trait]] <- model_dataframe(model_summary, "diet", trait, cohort = 'full', p_col = "Pr(>|W|)")
}

#combine for full diet dataframe
full_diet_coeffs <- do.call(rbind, diet_model_dfs) %>%
  arrange(adjusted_p)
rownames(full_diet_coeffs) <- NULL


## ----------------- 4.2 Sleep Phenotype Models ----------------- ##

#Repeating for sleep outcomes
sleep_models_full <- setNames(lapply(sleep_response_columns, function(response) {
  fit_gees(response, snp_columns, covariates, full_cohort_data, id_column = "fid", correlation_structure)
}), sleep_response_columns)

sleep_model_dfs <- list()

for(trait in names(sleep_models_full)) {
  model_summary <- sleep_models_full[[trait]]
  sleep_model_dfs[[trait]] <- model_dataframe(model_summary, "sleep", trait, cohort = 'full', p_col = "Pr(>|W|)")
}

full_sleep_coeffs <- do.call(rbind, sleep_model_dfs) %>%
  arrange(adjusted_p)
rownames(full_sleep_coeffs) <- NULL

#Now we can combine all dataframes into a single dataframe, since all required information exists in the columns (cohort, trait, category etc...)

#First rename the pvalue column in the glm's for the adult data to match the gee's to be able to combine. We will do this over all df's created. (Will make function to optimize this part here.)

adult_diet_coeffs <- adult_diet_coeffs %>%
  rename(pval = "Pr(>|t|)") %>%
  rename('Std.err' = "Std. Error") %>%
  select(-'t value')

adult_sleep_coeffs <- adult_sleep_coeffs %>%
  rename(pval = "Pr(>|t|)") %>%
  rename('Std.err' = "Std. Error") %>%
  select(-'t value')

children_diet_coeffs <- children_diet_coeffs %>%
  rename(pval = "Pr(>|W|)") %>%
  select(-Wald)

children_sleep_coeffs <- children_sleep_coeffs %>%
  rename(pval = "Pr(>|W|)")%>%
  select(-Wald)

full_diet_coeffs <- full_diet_coeffs %>%
  rename(pval = "Pr(>|W|)") %>%
  select(-Wald)

full_sleep_coeffs <- full_sleep_coeffs %>%
  rename(pval = "Pr(>|W|)") %>%
  select(-Wald)

#Combine to make a df of all tests conducted
all_tests <- do.call(rbind, list(adult_diet_coeffs, adult_sleep_coeffs, children_diet_coeffs, children_sleep_coeffs, full_diet_coeffs, full_sleep_coeffs))

all_tests <- all_tests %>%
  mutate(significant_raw = case_when(pval <= p_value_threshold ~ "TRUE", TRUE ~ "FALSE"))

#Now we determine which are significant and which are not conditionally based on the adj_OR_raw parameter, assigned in the master controller script.
if (adj_OR_raw == "adj") {
  sig_snps <- all_tests %>%
    filter(significant == TRUE) %>%
    filter(term %in% snp_columns) %>%
    select(-significant, -significant_raw)
} else {
  sig_snps <- all_tests %>%
    filter(significant_raw == TRUE) %>%
    filter(term %in% snp_columns) %>%
    select(-significant, -adjusted_p,-significant_raw)
}

sig_snps_grouped <- sig_snps %>%
  group_by(term) %>%
  summarise(Associations = n(),trait_categories = paste(unique(trait_category), collapse = ", "))

sig_snps_grouped <- sig_snps_grouped %>%
  mutate(trait_combined = case_when(
    trait_categories == "diet" ~ "Diet",
    trait_categories == "sleep" ~ "Sleep",
    trait_categories == "diet, sleep" | trait_categories == "sleep, diet" ~ "Both",
    TRUE ~ "Unknown"
  )) %>%
  select(-trait_categories)

num_sig_associations <- nrow(sig_snps)
uniquesigSNPs <- length(unique(sig_snps$term))

#Database[["Complete Datasets"]][["Full Cohort"]] <- full_cohort_data
#Database[["Complete Datasets"]][["Adults"]] <- adult_data
