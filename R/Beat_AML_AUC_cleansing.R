# --------------------------------------------------------
# Script Name: Beat_AML_AUC_cleansing.R
# Author: Katarina Willoch
# Date: 2025-06-01
#
# Description:
# This script reads in normalized viability data, calculates AUC scores using the MESS package for each duplicate (clearifying that MESS package was used https://www.synapse.org/Synapse:syn20940518/discussion/threadId=6617)
# and cleans it based on the similarity of the score as described here http://vizome.org/images/Inhibitor%20data%20workflow.jpg
#
# Sections:
#   1. Load libraries and data
#   2. Calculate AUC for each replicate
#   3. Check if within-panel replicate AUC difference <= 100 if yes average AUC if no remove replicates
#   4. Check if between-panel replicate AUC difference <=75 if yes average AUC if no remove replicates
#   4. Save outputs
# --------------------------------------------------------

#Install packages
install.packages('MESS')
#Load libraries
library(MESS)
library(dplyr)

#import beat aml datasets ----
beat_aml_inhibitor <- read_csv("~/Desktop/UiO/Project 1/Data/Second run/beat_aml_inhibitor.csv")
beat_aml_inhibitor$...1 <- NULL

# Function to calculate AUC using MESS package
calculate_auc <- function(df, resp) {
  x <- df$well_concentration
  y <- df$normalized_viability #changed from avg1_response
  auc_value <- auc(x, y, type = "linear")
  return(auc_value)
}

result <- beat_aml_inhibitor %>%
  group_by(dbgap_subject_id, dbgap_dnaseq_sample, dbgap_rnaseq_sample, inhibitor, inhibitor_panel, replicate) %>%
  dplyr::summarise(auc = calculate_auc(cur_data(), normalized_viability)) %>%
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

transformed_data <- filtered_data %>% #changed from beat_aml_inhibitor
  group_by(dbgap_subject_id, dbgap_dnaseq_sample, dbgap_rnaseq_sample, inhibitor, inhibitor_panel, well_concentration) %>%
  dplyr::summarise(avg1_response = mean(normalized_viability)) %>%
  ungroup() %>%
  distinct()  # Remove duplicates in case there are any

print(transformed_data)


transformed_data <- transformed_data %>%
  group_by(dbgap_subject_id, dbgap_dnaseq_sample, dbgap_rnaseq_sample, inhibitor, well_concentration) %>%
  dplyr::summarise(avg2_response = mean(avg1_response)) %>%
  ungroup() %>%
  distinct()  # Remove duplicates in case there are any

print(transformed_data)
write_csv(transformed_data, "~/Desktop/UiO/Project 1/Data/Second run/beat_aml_inhibitor_cleaned.csv")
