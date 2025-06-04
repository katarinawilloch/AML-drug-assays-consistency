# --------------------------------------------------------
# Script Name: DSS_calculation.R
# Author: Katarina Willoch
# Date: 2025-06-01
#
# Description:
# This script reads in survey data, processes and cleans it,
# generates summary statistics, and produces visualizations
# for inclusion in the quarterly report.
#
# Sections:
#   1. Load libraries and data
#   2. Data wrangling and transformation
#   3. DSS calculation
#   4. Save outputs
# --------------------------------------------------------

#Install packages
install.packages('MESS')
#Load libraries
library(MESS)
library(dplyr)

#import beat aml datasets ----
beat_aml_inhibitor <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_inhibitor.csv")
beat_aml_inhibitor$...1 <- NULL
View(beat_aml_inhibitor)

common_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/common_drugnames_pubchem.csv")


# Function to calculate AUC using MESS package
calculate_auc <- function(df, resp) {
  x <- df$well_concentration
  y <- df$avg1_response
  auc_value <- auc(x, y, type = "linear")
  return(auc_value)
}

result <- beat_aml_inhibitor %>%
  group_by(dbgap_subject_id, dbgap_dnaseq_sample, dbgap_rnaseq_sample, inhibitor, inhibitor_panel, replicate) %>%
  summarise(auc = calculate_auc(cur_data(), normalized_viability)) %>%
  ungroup()

print(result)

# Compare AUC values for different 'col2' within the same 'col1'
comparison <- result %>%
  spread(key = replicate, value = auc) %>%
  mutate(difference = abs(`1` - `2`)) %>%
  filter(difference < 100 | is.na(difference)) 

# Identify the unique combinations of 'col1' and 'col2' to be removed
groups_to_remove <- comparison %>%
  distinct(dbgap_subject_id, dbgap_dnaseq_sample, dbgap_rnaseq_sample, inhibitor, inhibitor_panel)

filtered_data <- beat_aml_inhibitor %>%
  left_join(groups_to_remove, by = c('dbgap_subject_id', 'dbgap_dnaseq_sample', 'dbgap_rnaseq_sample', 'inhibitor', 'inhibitor_panel'))

print(filtered_data)

transformed_data <- beat_aml_inhibitor %>%
  group_by(dbgap_subject_id, dbgap_dnaseq_sample, dbgap_rnaseq_sample, inhibitor, inhibitor_panel, well_concentration) %>%
  summarise(avg1_response = mean(normalized_viability)) %>%
  ungroup() %>%
  distinct()  # Remove duplicates in case there are any

print(transformed_data)


transformed_data <- transformed_data %>%
  group_by(dbgap_subject_id, dbgap_dnaseq_sample, dbgap_rnaseq_sample, inhibitor, well_concentration) %>%
  summarise(avg2_response = mean(avg1_response)) %>%
  ungroup() %>%
  distinct()  # Remove duplicates in case there are any

print(transformed_data)
write_csv(transformed_data, "~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_inhibitor_cleaned.csv")
