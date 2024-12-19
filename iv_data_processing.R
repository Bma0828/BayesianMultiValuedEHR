library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

# Function to read compressed CSV files
read_gz_csv <- function(filepath) {
  fread(cmd = paste("gzip -dc", filepath))
}

# Extract vasopressor treatments
get_vasopressor_data <- function(icu_path) {

  d_items <- read_gz_csv(file.path(icu_path, "d_items.csv.gz"))
  inputevents <- read_gz_csv(file.path(icu_path, "inputevents.csv.gz"))
  
  # Define vasopressor categories
  vp_categories <- list(
    norepinephrine = "norepinephrine",
    vasopressin = "vasopressin",
    phenylephrine = "phenylephrine"
  )
  
  # Get itemids for each vasopressor
  vp_items <- lapply(vp_categories, function(x) {
    d_items %>%
      filter(grepl(x, tolower(label))) %>%
      pull(itemid)
  })
  
  # Extract vasopressor administrations
  vp_data <- inputevents %>%
    filter(itemid %in% unlist(vp_items)) %>%
    # Remove invalid times
    filter(starttime < endtime) %>%
    # Remove missing amounts
    filter(!is.na(amount)) %>%
    # Calculate rate
    mutate(
      duration_hours = as.numeric(difftime(endtime, starttime, units = "hours")),
      rate = amount / duration_hours,
      vp_type = case_when(
        itemid %in% vp_items$norepinephrine ~ "norepinephrine",
        itemid %in% vp_items$vasopressin ~ "vasopressin",
        itemid %in% vp_items$phenylephrine ~ "phenylephrine"
      )) %>%
    # Remove extreme rates
    filter(rate > 0, rate < quantile(rate, 0.99, na.rm = TRUE)) %>%
    select(stay_id, itemid, starttime, endtime, amount, vp_type, caregiver_id)
  
  return(vp_data)
}

# Calculate physician preferences (Instrumental Variable)
calculate_physician_preferences <- function(vp_data) {
  # Calculate HHI for each physician
  physician_preferences <- vp_data %>%
    # Group by physician
    group_by(caregiver_id) %>%
    summarize(
      # Count total prescriptions per physician
      total_prescriptions = n(),
      
      # Calculate HHI
      hhi = {
        drug_counts <- table(vp_type)
        market_shares <- drug_counts / sum(drug_counts)
        # Calculate HHI (sum of squared market shares)
        sum(market_shares^2)
      },
      
      # Get the most prescribed drug type
      primary_vp_type = names(which.max(table(vp_type))),
      
      # Calculate concentration ratio (share of most prescribed drug)
      concentration_ratio = max(table(vp_type)) / n(),
      
      .groups = "drop"
    ) %>%
    # Filter for physicians with minimum number of prescriptions
    filter(total_prescriptions >= 10)
  
  return(physician_preferences)
}

# Extract covariates
get_covariates <- function(hosp_path, icu_path) {
  
  patients <- read_gz_csv(file.path(hosp_path, "patients.csv.gz")) %>%
    select(subject_id, gender, anchor_age)
  admissions <- read_gz_csv(file.path(hosp_path, "admissions.csv.gz")) %>%
    select(subject_id, hadm_id, admission_type, insurance, language, marital_status)
  icustays <- read_gz_csv(file.path(icu_path, "icustays.csv.gz"))
  
  #### Function to process CHARTEVENTS in chunks #### 
  extract_chartevents_columns <- function(file_path, output_path) {
    # Define vital signs we need
    vital_items <- c(
      "Heart Rate" = 220045,
      "SBP" = 220050,
      "DBP" = 220051,
      "MAP" = 220052,
      "Respiratory Rate" = 220210,
      # "Temperature" = 223761,
      "SpO2" = 220277
    )
    
    cat("Starting to read CHARTEVENTS file...\n")
    
    tryCatch({
      # Read only necessary columns and filter for vital signs
      chartevents_subset <- fread(
        cmd = paste("gzip -dc", file_path),
        select = c("subject_id", "hadm_id", "stay_id", "charttime", "itemid", "valuenum"),
        colClasses = list(
          numeric = c("subject_id", "hadm_id", "stay_id", "itemid", "valuenum"),
          POSIXct = "charttime"
        ),
        nThread = 4  # Adjust based on your CPU
      )[itemid %in% vital_items]
      
      # Remove rows with NA valuenum
      chartevents_subset <- chartevents_subset[!is.na(valuenum)]
      
      cat(sprintf("Extracted %d rows with vital signs\n", nrow(chartevents_subset)))
      
      # Save the subset as compressed RDS
      saveRDS(chartevents_subset, file = output_path, compress = TRUE)
      cat(sprintf("Saved subset to %s\n", output_path))
      
      # Return the subset
      return(chartevents_subset)
      
    }, error = function(e) {
      cat("Error in processing:", conditionMessage(e), "\n")
      return(NULL)
    })
  }
  
  extract_labevents_columns <- function(file_path, output_path) {
    # Define key lab tests we need
    key_labs <- c(
      "Lactate" = 50813,
      "Creatinine" = 50912,
      "Glucose" = 50931,
      "Potassium" = 50971,
      "Sodium" = 50983,
      "Hemoglobin" = 51222,
      "WBC" = 51301
    )
    
    cat("Starting to read LABEVENTS file...\n")
    
    tryCatch({
      # Read only necessary columns and filter for key labs
      labevents_subset <- fread(
        cmd = paste("gzip -dc", file_path),
        select = c("subject_id", "hadm_id", "itemid", "charttime", "valuenum"),
        colClasses = list(
          numeric = c("subject_id", "hadm_id", "itemid", "valuenum"),
          POSIXct = "charttime"
        ),
        nThread = 4  # Adjust based on your CPU
      )[itemid %in% key_labs]
      
      # Remove rows with NA valuenum
      labevents_subset <- labevents_subset[!is.na(valuenum)]
      
      # Remove outliers
      labevents_subset <- labevents_subset[, {
        # Calculate quantiles for this specific lab test
        q_lower <- quantile(valuenum, 0.001, na.rm = TRUE)
        q_upper <- quantile(valuenum, 0.999, na.rm = TRUE)
        # Keep all columns but filter rows based on outlier criteria
        .SD[valuenum >= q_lower & valuenum <= q_upper]
      }, by = itemid]
      
      cat(sprintf("Extracted %d rows with key lab tests\n", nrow(labevents_subset)))
      
      # Print summary of lab tests
      cat("\nNumber of measurements for each lab test:\n")
      lab_summary <- labevents_subset[, .N, by = itemid]
      lab_summary[, lab_name := names(key_labs)[match(itemid, key_labs)]]
      print(lab_summary)
      
      # Save the subset as compressed RDS
      saveRDS(labevents_subset, file = output_path, compress = TRUE)
      cat(sprintf("Saved subset to %s\n", output_path))
      
      # Return the subset
      return(labevents_subset)
      
    }, error = function(e) {
      cat("Error in processing:", conditionMessage(e), "\n")
      return(NULL)
    })
  }
  
  # output_file <- "labevents_subset.rds"
  # labevents_subset <- extract_labevents_columns(
  #   file_path = file.path(hosp_path, "labevents.csv.gz"),
  #   output_path = output_file
  # )
  # 
  # # Check the result
  # if (!is.null(labevents_subset)) {
  #   cat("\nFirst few rows of the subset:\n")
  #   print(head(labevents_subset))
  # 
  #   cat("\nMemory usage of subset:\n")
  #   print(object.size(labevents_subset), units = "Mb")
  # 
  #   # Basic statistics for each lab test
  #   cat("\nSummary statistics for each lab test:\n")
  #   summary_stats <- labevents_subset[, .(
  #     mean = mean(valuenum, na.rm = TRUE),
  #     sd = sd(valuenum, na.rm = TRUE),
  #     median = median(valuenum, na.rm = TRUE),
  #     min = min(valuenum, na.rm = TRUE),
  #     max = max(valuenum, na.rm = TRUE),
  #     n = .N
  #   ), by = .(itemid)]
  #   summary_stats[, lab_name := names(key_labs)[match(itemid, key_labs)]]
  #   print(summary_stats)
  # }
  # 
  # output_file <- "chartevents_subset.rds"
  # chartevents_subset <- extract_chartevents_columns(
  #   file_path = file.path(icu_path, "chartevents.csv.gz"),
  #   output_path = output_file
  # )
  # 
  # # Check the result
  # if (!is.null(chartevents_subset)) {
  #   cat("\nFirst few rows of the subset:\n")
  #   print(head(chartevents_subset))
  #   
  #   cat("\nMemory usage of subset:\n")
  #   print(object.size(chartevents_subset), units = "Mb")
  #   
  #   cat("\nSummary of vital signs:\n")
  #   print(table(chartevents_subset$itemid))
  # }
  
  chartevents <- readRDS("chartevents_subset.rds")
  labevents <- readRDS("labevents_subset.rds")
  
  # Demographics and admission info
  demographics <- patients %>%
    left_join(admissions, by = "subject_id") %>%
    distinct()
  
  # Vital signs in first 24h
  vital_items <- c(
    "Heart Rate" = 220045,
    "SBP" = 220050,
    "DBP" = 220051,
    "MAP" = 220052,
    "Respiratory Rate" = 220210,
    # "Temperature" = 223761,
    "SpO2" = 220277
  )
  
  vitals_24h <- chartevents %>%
    # Join with ICU stays to get admission time
    inner_join(
      icustays %>% 
        select(subject_id, hadm_id, stay_id, intime) %>%
        distinct(),
      by = c("subject_id", "hadm_id", "stay_id")) %>%
    # Filter for valid vital signs
    filter(itemid %in% vital_items) %>%
    # Add reasonable ranges for vital signs
    filter(case_when(
      # Heart Rate: 20-300 bpm
      itemid == 220045 ~ valuenum >= 20 & valuenum <= 300,
      # SBP: 40-300 mmHg
      itemid == 220050 ~ valuenum >= 40 & valuenum <= 300,
      # DBP: 20-200 mmHg
      itemid == 220051 ~ valuenum >= 20 & valuenum <= 200,
      # MAP: 30-250 mmHg
      itemid == 220052 ~ valuenum >= 30 & valuenum <= 250,
      # Respiratory Rate: 4-60 breaths/min
      itemid == 220210 ~ valuenum >= 4 & valuenum <= 60,
      # Temperature: 32-42 Celsius
      # itemid == 223761 ~ valuenum >= 32 & valuenum <= 42,
      # SpO2: 50-100%
      itemid == 220277 ~ valuenum >= 50 & valuenum <= 100,
      TRUE ~ FALSE
    )) %>%
    group_by(stay_id, itemid) %>%
    summarise(
      mean_value = mean(valuenum, na.rm = TRUE),
      min_value = min(valuenum, na.rm = TRUE),
      max_value = max(valuenum, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      vital_name = case_when(
        itemid == 220045 ~ "HeartRate",
        itemid == 220050 ~ "SBP",
        itemid == 220051 ~ "DBP",
        itemid == 220052 ~ "MAP",
        itemid == 220210 ~ "RespRate",
        # itemid == 223761 ~ "Temperature",
        itemid == 220277 ~ "SpO2"
      )
    ) %>%
    select(-itemid) %>%
    pivot_wider(
      id_cols = stay_id,
      names_from = vital_name,
      values_from = c(mean_value, min_value, max_value),
      names_glue = "{.value}_{vital_name}"
    )
  
  # Lab values in first 24h
  key_labs <- c(
    "Lactate" = 50813,
    "Creatinine" = 50912,
    "Glucose" = 50931,
    "Potassium" = 50971,
    "Sodium" = 50983,
    "Hemoglobin" = 51222
  )
  
  # Join labevents with ICU stays
  labs_24h <- labevents %>%
    # Filter to keep only the lab items we want
    filter(itemid %in% key_labs) %>%
    # Remove any obvious outliers or invalid values
    filter(!is.na(valuenum)) %>%
    # Add reasonable ranges for lab values
    filter(case_when(
      # Lactate: 0.1-30 mmol/L
      itemid == 50813 ~ valuenum >= 0.1 & valuenum <= 30,
      # Creatinine: 0.1-20 mg/dL
      itemid == 50912 ~ valuenum >= 0.1 & valuenum <= 20,
      # Glucose: 20-1000 mg/dL
      itemid == 50931 ~ valuenum >= 20 & valuenum <= 1000,
      # Potassium: 2-10 mEq/L
      itemid == 50971 ~ valuenum >= 2 & valuenum <= 10,
      # Sodium: 100-180 mEq/L
      itemid == 50983 ~ valuenum >= 100 & valuenum <= 180,
      # Hemoglobin: 2-25 g/dL
      itemid == 51222 ~ valuenum >= 2 & valuenum <= 25,
      TRUE ~ FALSE
    )) %>%
    # Join with ICU stays
    inner_join(
      icustays %>% 
        select(subject_id, hadm_id, stay_id, intime) %>%
        distinct(),
      by = c("subject_id", "hadm_id"),
      relationship = "many-to-many") %>%
    # # Filter for first 24 hours of ICU stay
    # filter(charttime >= intime & 
    #          charttime <= intime + hours(24)) %>%
    # Add time from ICU admission
    mutate(
      hours_from_admit = as.numeric(difftime(charttime, intime, units = "hours"))
    ) %>%
    # Calculate summary statistics
    group_by(stay_id, itemid) %>%
    summarise(
      first_value = first(valuenum, order_by = charttime),
      last_value = last(valuenum, order_by = charttime),
      delta_value = last_value - first_value,
      min_value = min(valuenum, na.rm = TRUE),
      max_value = max(valuenum, na.rm = TRUE),
      mean_value = mean(valuenum, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Pivot wider to have one row per stay_id
    pivot_wider(
      id_cols = stay_id,
      names_from = itemid,
      values_from = c(
        mean_value, min_value, max_value,
        first_value, last_value, delta_value
      ),
      names_prefix = "lab_",
      values_fill = list(
        mean_value = NA_real_,
        min_value = NA_real_,
        max_value = NA_real_,
        first_value = NA_real_,
        last_value = NA_real_,
        delta_value = NA_real_
      )
    )
  
  # Add lab names to column names for clarity
  labs_24h <- labs_24h %>%
    rename_with(
      .fn = function(x) {
        sapply(x, function(col_name) {
          # Define pattern matching for all types
          patterns <- c(
            "mean_value_lab_" = "_mean",
            "min_value_lab_" = "_min",
            "max_value_lab_" = "_max",
            "first_value_lab_" = "_first",
            "last_value_lab_" = "_last",
            "delta_value_lab_" = "_delta"
          )
          
          # Check if column matches any pattern
          for (pattern in names(patterns)) {
            if (grepl(pattern, col_name)) {
              lab_id <- as.numeric(str_extract(col_name, "\\d+"))
              lab_name <- names(key_labs)[match(lab_id, key_labs)]
              return(paste0(lab_name, patterns[pattern]))
            }
          }
          return(col_name)
        })
      }
    )
  
  return(list(
    demographics = demographics,
    vitals_24h = vitals_24h,
    labs_24h = labs_24h
  ))
}

# Create final analysis dataset
create_analysis_dataset <- function(vp_data, physician_preferences, covariates, 
                                    hosp_path, icu_path) {
  # Load ICU stays and admissions for outcomes
  icustays <- read_gz_csv(file.path(icu_path, "icustays.csv.gz"))
  admissions <- read_gz_csv(file.path(hosp_path, "admissions.csv.gz"))
  
  # Determine primary vasopressor for each stay
  primary_vp <- vp_data %>%
    group_by(stay_id, vp_type, caregiver_id) %>%
    summarise(
      total_amount = sum(amount),
      .groups = "drop"
    ) %>%
    group_by(stay_id) %>%
    slice_max(total_amount, n = 1) %>%
    select(stay_id, primary_vp = vp_type, caregiver_id)
  
  # Create outcomes
  outcomes <- icustays %>%
    left_join(admissions, by = c("subject_id", "hadm_id")) %>%
    mutate(
      # Primary outcome: 28-day mortality
      mortality = case_when(
        hospital_expire_flag == 1 ~ 1,  # Died in hospital
        discharge_location %in% c("DIED", "HOSPICE", "HOSPICE-HOME", "HOSPICE-MEDICAL FACILITY") ~ 1,  # Died/hospice at discharge
        TRUE ~ 0  # Survived
      ),
      # Secondary outcomes
      los_days = as.numeric(difftime(outtime, intime, units = "days")),
      
      # Hospital-free days: 28 - LOS if survived, 0 if died
      hospital_free_days = case_when(
        mortality == 1 ~ 0,  # 0 days if died
        TRUE ~ pmin(28, 28 - los_days)  # 28 - LOS, capped at 28 for survivors
      )
    )
  
  # Combine everything
  analysis_data <- primary_vp %>%
    left_join(outcomes, by = "stay_id") %>%
    left_join(covariates$demographics, by = c("subject_id", "hadm_id")) %>%
    left_join(covariates$vitals_24h, by = "stay_id") %>%
    left_join(covariates$labs_24h, by = "stay_id") %>%
    left_join(physician_preferences, by = c("caregiver_id"))
  
  return(analysis_data)
}

# Main processing function
process_vasopressor_study <- function(hosp_path, icu_path) {
  # Get vasopressor data
  vp_data <- get_vasopressor_data(icu_path)
  
  # Calculate physician preferences
  phys_prefs <- calculate_physician_preferences(vp_data)

  # Get covariates
  covs <- get_covariates(hosp_path, icu_path)
  print(paste("Demographics rows:", nrow(covs$demographics)))
  print(paste("Vitals rows:", nrow(covs$vitals_24h)))
  print(paste("Labs rows:", nrow(covs$labs_24h)))
  
  # Create final dataset
  analysis_data <- create_analysis_dataset(vp_data, phys_prefs, covs, hosp_path, icu_path)
  print(paste("Final analysis rows:", nrow(analysis_data)))
  
  return(analysis_data)
}

# hosp_path <- "/Users/anemos/Desktop/有事/Harvard/Research/mimiciv_demo/hosp/"
# icu_path <- "/Users/anemos/Desktop/有事/Harvard/Research/mimiciv_demo/icu/"

hosp_path <- "/Users/anemos/Desktop/有事/Harvard/Research/mimiciv/hosp/"
icu_path <- "/Users/anemos/Desktop/有事/Harvard/Research/mimiciv/icu/"
analysis_data <- process_vasopressor_study(hosp_path, icu_path)
saveRDS(analysis_data, "analysis_data_full.rds")

analysis_data <- readRDS("analysis_data_full.rds")
analysis_data <- readRDS("analysis_data_demo.rds")

#######################################################
# Function to clean and recode categorical variables

# Main preprocessing function
preprocess_mimic_data <- function(data) {
  processed_data <- data %>%
    # 1. Clean categorical variables
    mutate(
      # Recode discharge location to binary
      discharge_status = case_when(
        discharge_location %in% c("DIED", "HOSPICE", "HOSPICE-HOME", 
                                  "HOSPICE-MEDICAL FACILITY") ~ "Death",
        discharge_location %in% c("HOME", "HOME HEALTH CARE", 
                                  "HOME WITH HOME IV PROVIDR") ~ "Home",
        TRUE ~ "Other"
      ),
      
      # Recode insurance to three levels
      insurance_group = case_when(
        grepl("MEDICAID|MEDICARE", insurance.x, ignore.case = TRUE) ~ "MedicareMedicaid",
        grepl("PRIVATE", insurance.x, ignore.case = TRUE) ~ "Private",
        TRUE ~ "Other"
      ),
      
      # Recode race to four categories
      race_group = case_when(
        grepl("BLACK|AFRICAN", race, ignore.case = TRUE) ~ "Black",
        grepl("WHITE", race, ignore.case = TRUE) ~ "White",
        grepl("ASIAN", race, ignore.case = TRUE) ~ "Asian",
        TRUE ~ "Other"
      )
    ) %>%
    
    # 2. Remove duplicate columns and standardize names
    dplyr::select(
      # Remove .y versions of duplicated columns
      -matches("\\.y$"),
      
      # Remove ED time columns (as they're just timestamps for ED visits)
      -edregtime, -edouttime,
      
      # Keep essential columns
      stay_id, subject_id, hadm_id, primary_vp, caregiver_id,
      first_careunit, last_careunit,
      intime, outtime, admittime, dischtime, deathtime,
      
      # Keep processed categorical variables
      discharge_status, insurance_group, race_group,
      
      # Keep other important variables
      mortality, los_days, hospital_free_days,
      gender, anchor_age, hospital_expire_flag,
      
      # Keep all lab values and vital signs
      matches("_mean$|_min$|_max$|_first$|_last$|_delta$"),
      
      # Keep IV-related variables
      hhi, primary_vp_type, concentration_ratio, total_prescriptions
    ) %>%
    
    # 3. Rename remaining .x columns
    rename_with(~str_remove(., "\\.x$")) %>%
    
    # 4. Convert categorical variables to factors
    mutate(across(c(discharge_status, insurance_group, race_group,
                    gender, first_careunit, last_careunit), as.factor)) %>%
    
    # 5. Add derived variables
    mutate(
      # Flag for weekend admission
      weekend_admission = as.factor(weekdays(admittime) %in% c("Saturday", "Sunday")),
      
      # Flag for night admission (7pm - 7am)
      night_admission = as.factor(hour(admittime) < 7 | hour(admittime) >= 19),
      
      # Season of admission
      admission_season = as.factor(case_when(
        month(admittime) %in% c(12, 1, 2) ~ "Winter",
        month(admittime) %in% c(3, 4, 5) ~ "Spring",
        month(admittime) %in% c(6, 7, 8) ~ "Summer",
        month(admittime) %in% c(9, 10, 11) ~ "Fall"
      ))
    )
  
  return(processed_data)
}

# Function to check data quality after preprocessing
check_processed_data <- function(data) {
  cat("Data Quality Report:\n\n")
  
  # Check categorical variables
  cat("1. Categorical Variables Distribution:\n")
  cat("\nDischarge Status:\n")
  print(table(data$discharge_status))
  
  cat("\nInsurance Groups:\n")
  print(table(data$insurance_group))
  
  cat("\nRace Groups:\n")
  print(table(data$race_group))
  
  # Check missing values
  cat("\n2. Missing Values Summary:\n")
  missing_summary <- colSums(is.na(data))
  print(missing_summary[missing_summary > 0])
  
  # Check numeric variables
  cat("\n3. Numeric Variables Summary:\n")
  numeric_cols <- sapply(data, is.numeric)
  print(summary(data[, numeric_cols]))
  
  # Return missing value proportions
  return(missing_summary/nrow(data))
}

# Process your data
processed_data <- preprocess_mimic_data(analysis_data)

# Check the quality of processed data
quality_report <- check_processed_data(processed_data)

processed_data_ <- processed_data %>% 
  dplyr::select(-c(stay_id, hadm_id, first_careunit, last_careunit, 
            intime, outtime, los, admittime, dischtime, deathtime, admit_provider_id,
            admission_location, discharge_location, insurance, race, total_prescriptions,
            primary_vp_type, language))

# Function to encode categorical variables
encode_categorical_variables <- function(data) {
  data <- data %>%
    mutate(
      # Encode primary_vp (0: norepinephrine, 1: phenylephrine, 2: vasopressin)
      primary_vp = case_when(
        primary_vp == "norepinephrine" ~ 0,
        primary_vp == "phenylephrine" ~ 1,
        primary_vp == "vasopressin" ~ 2,
        TRUE ~ NA_real_
      ),
      
      # Encode gender
      gender = case_when(
        gender == "M" ~ 0,
        gender == "F" ~ 1,
        TRUE ~ NA_real_
      ),
      
      # Encode admission_type (simplified to 2 categories)
      admission_type = case_when(
        is.na(admission_type) ~ NA_real_,
        admission_type %in% c("EW EMER.", "DIRECT EMER.") ~ 0,  
        TRUE ~ 1  
      ),
      
      # Encode marital_status
      marital_status = case_when(
        marital_status == "SINGLE" ~ 0,
        marital_status == "MARRIED" ~ 1,
        marital_status == "WIDOWED" ~ 2,
        marital_status == "DIVORCED" ~ 3,
        marital_status == "" ~ NA_real_,
        TRUE ~ NA_real_
      ),
      
      # Encode discharge_status
      discharge_status = case_when(
        discharge_status == "Death" ~ 0,
        discharge_status == "Home" ~ 1,
        discharge_status == "Other" ~ 2,
        TRUE ~ NA_real_
      ),
      
      # Encode insurance_group
      insurance_group = case_when(
        insurance_group == "MedicareMedicaid" ~ 0,
        insurance_group == "Private" ~ 1,
        insurance_group == "Other" ~ 2,
        TRUE ~ NA_real_
      ),
      
      # Encode race_group
      race_group = case_when(
        race_group == "White" ~ 0,
        race_group == "Black" ~ 1,
        race_group == "Asian" ~ 2,
        race_group == "Other" ~ 3,
        TRUE ~ NA_real_
      ),
      
      # Convert logical variables to numeric
      mortality = as.numeric(mortality),           # TRUE -> 1, FALSE -> 0
      weekend_admission = as.numeric(weekend_admission),  # TRUE -> 1, FALSE -> 0
      night_admission = as.numeric(night_admission),      # TRUE -> 1, FALSE -> 0
      
      # Encode admission_season
      admission_season = case_when(
        admission_season == "Winter" ~ 0,
        admission_season == "Spring" ~ 1,
        admission_season == "Summer" ~ 2,
        admission_season == "Fall" ~ 3,
        TRUE ~ NA_real_
      )
    )
  
  return(data)
}

encoded_data <- encode_categorical_variables(processed_data_)

# saveRDS(encoded_data, "final_data_full.rds")
# encoded_data <- readRDS("final_data_full.rds")
# write.csv(encoded_data, file = "filtered_prescriptions.csv", row.names = FALSE)

saveRDS(encoded_data, "final_data_demo.rds")
encoded_data <- readRDS("final_data_demo.rds")
write.csv(encoded_data, file = "filtered_prescriptions_demo.csv", row.names = FALSE)








