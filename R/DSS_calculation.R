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

# Install and import libraries ----
library(readr)
library("ggplot2")
library(ggforce)
library(dplyr)
library(stringr)
library(pheatmap)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(readxl)
library(DataExplorer)
library(Hmisc)
library(htmltools)
library(VennDiagram)
library(grid)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace(c("sva", "pcaMethods"), quietly = TRUE))
  BiocManager::install(c("sva", "pcaMethods"))

source('~/Desktop/UiO/DSS-v2.0/DSS.R')
source('~/Desktop/UiO/DSS-v2.0/HelperFunctions.R')

#Return calculated DSS scores using https://github.com/yingjchen/DSS-v2.0/tree/main?tab=readme-ov-file
#@param df dataframe of raw normalized response values
#@param drug_name column name referring to the drugs names for restructuring the dataframe
#@param conc column name referring to the column in the dataset containing the concentration for restructuring the dataframe
#@param screen_name column name referring to the column in the dataset containing unique sample ids for restructuring the dataframe
#@param response column name referring to the column in the dataset containing normalized viability/inhibition values for restructuring the dataframe
#@param filename the location to save the resulting dataframe of dss values
#@param graph (TRUE/FALSE) if TRUE generates dose response graphs in the set_dir location; this param gets used in the CALC_METRICS function from https://github.com/yingjchen/DSS-v2.0/tree/main?tab=readme-ov-file HelperFunctions.R
#@param viability (TRUE/FALSE) refers to weather the response is given as viability or inhibition; this param gets used in the DOSE_RESPONSE_PROCESS function from https://github.com/yingjchen/DSS-v2.0/tree/main?tab=readme-ov-file HelperFunctions.R
#@param set_dir refers to the location the dose response graphs should be saved if TRUE
#@return dataframe with the dss, IC50 and AUC scores 
#@examples
#calc_dss(enserink_breeze_input, 'drug', 'dose', 'Patient.ID', 'response', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink_full_drug_set_testing.csv', viability = TRUE, graph = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/Enserink")
calc_dss <- function(df, drug_name, conc, screen_name, response, filename, graph = FALSE, viability = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures"){
  for_dss <- df[,c(drug_name, conc, screen_name, response)]
  print(for_dss)
  colnames(for_dss)[colnames(for_dss) == drug_name] <- "DRUG_NAME"
  colnames(for_dss)[colnames(for_dss) == conc] <- "CONCENTRATION_nM"
  colnames(for_dss)[colnames(for_dss) == screen_name] <- "SCREEN_NAME"
  colnames(for_dss)[colnames(for_dss) == response] <- "CELL_VIABILITY"
  
  for_dss$SCREEN_NAME <- as.character(for_dss$SCREEN_NAME)
  for_dss$CONCENTRATION_nM <- as.double(for_dss$CONCENTRATION_nM)
  
  df_dose_response <- DOSE_RESPONSE_PROCESS(for_dss, viability = viability)
  print(df_dose_response)
  if(graph == TRUE){
    setwd(set_dir)
    print(df_dose_response[[1]])
    print(df_dose_response[[2]])
    df_metrics <- CALC_METRICS(df_dose_response[[1]], df_dose_response[[2]], graph = TRUE)
    setwd("/Users/katarinawilloch/")
  }else{
    df_metrics <- CALC_METRICS(df_dose_response[[1]], df_dose_response[[2]], graph = FALSE)
  }
  
  write_csv(df_metrics, filename)
  return(df_metrics)
}

#File locations 
initial_cleansing <- "~/Desktop/UiO/Project 1/Data/Initial cleansing/"

#import common drug names dataset
common_drugs <- read_csv(paste0(initial_cleansing,"common_drugnames_pubchem.csv"))

#Beat AML
#import beat aml datasets ----
beat_aml_inhibitor <- read_csv(paste0(initial_cleansing, "beat_aml_inhibitor.csv"))
beat_aml_inhibitor$...1 <- NULL
#get only beat aml dataset for the common drug names
beat_aml_for_dss <- inner_join(common_drugs, beat_aml_inhibitor, by=c("beat_aml_drug_name"="inhibitor"))
#create distinct sample ids 
beat_aml_for_dss <- beat_aml_for_dss %>%
  mutate(patient_id = paste(dbgap_subject_id, inhibitor_panel,replicate,run_index, sep = "_")) %>% as.data.frame()
# Display the dataframe
duplicates <- beat_aml_for_dss %>%
  group_by(patient_id, pubchem_drug_name) %>%
  filter(n() > 2) %>%
  ungroup()

# Show the duplicate values
print(duplicates)
print(beat_aml_for_dss)
beat_aml_for_dss <- subset(beat_aml_for_dss, select=c('pubchem_drug_name', 'well_concentration', 'patient_id', 'normalized_viability'))
colnames(beat_aml_for_dss)[colnames(beat_aml_for_dss) == "pubchem_drug_name"] <- "DRUG_NAME"
colnames(beat_aml_for_dss)[colnames(beat_aml_for_dss) == "well_concentration"] <- "CONCENTRATION_nM"
colnames(beat_aml_for_dss)[colnames(beat_aml_for_dss) == "patient_id"] <- "SCREEN_NAME"
colnames(beat_aml_for_dss)[colnames(beat_aml_for_dss) == "normalized_viability"] <- "CELL_VIABILITY"
beat_aml_for_dss$SCREEN_NAME <- as.character(beat_aml_for_dss$SCREEN_NAME)
beat_aml_for_dss

# Display the dataframe
duplicates <- beat_aml_for_dss %>%
  arrange(SCREEN_NAME, DRUG_NAME) %>%
  group_by(SCREEN_NAME, DRUG_NAME) %>%
  filter(n() > 2) %>%
  ungroup() 

# Show the duplicate values
print(duplicates)

beat_aml_for_dss$CELL_VIABILITY <- beat_aml_for_dss$CELL_VIABILITY/100
df_dose_response_beat_aml <- DOSE_RESPONSE_PROCESS(beat_aml_for_dss, viability = TRUE)
setwd("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/BeatAML/test")
df_metrics_beat_aml <- CALC_METRICS(df_dose_response_beat_aml[[1]], df_dose_response_beat_aml[[2]], graph = TRUE)
setwd("/Users/katarinawilloch/")

write_csv(df_metrics_beat_aml, '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_beat_aml.csv')

beat_aml_raw_response <- read.csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_inhibitor_cleaned.csv')
#reformatting to get percent viability
beat_aml_raw_response$avg2_response <- beat_aml_raw_response$avg2_response/100
#making unique sample ids
beat_aml_raw_response <- beat_aml_raw_response %>%
  mutate(patient_id = paste(dbgap_subject_id, dbgap_dnaseq_sample,dbgap_rnaseq_sample, sep = "_")) %>% as.data.frame()
#transforming concentration from um to nm
beat_aml_raw_response$well_concentration <- beat_aml_raw_response$well_concentration*1000
#using function calc_dss to calculate dss values
calc_dss(beat_aml_raw_response, 'inhibitor', 'well_concentration', 'patient_id', 'avg2_response', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_beat_aml_full_drug_set_v1.csv', graph = TRUE, viability = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/BeatAML")


#Oslo
#importing the oslo dataset
enserink_breeze_input <- read_csv("/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_dose_response.csv")
#transforming concentration from m to nm
enserink_breeze_input$dose <- enserink_breeze_input$dose*10^9
#transforming inhibition values to percent viability
enserink_breeze_input$response <-1 - enserink_breeze_input$y_true/100
#using function calc_dss to get dss values
calc_dss(enserink_breeze_input, 'drug', 'dose', 'Patient.ID', 'response', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink_full_drug_set_testing.csv', viability = TRUE, graph = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/Enserink")


#Helsinki
#import helsinki dataset
fimm_raw_response <- read.csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_raw_dose_responses.csv.csv')
#transform concentrations from ... to nm
fimm_raw_response$Final.Conc <- fimm_raw_response$Final.Conc*10
#transforming viability values to percent viability
fimm_raw_response$doseResponses <- fimm_raw_response$doseResponses/100
#removing drug 'TG100-115'
fimm_raw_response_df <- subset(fimm_raw_response, drug != 'TG100-115')
#using function calc_dss to get dss values
calc_dss(fimm_raw_response_df, 'drug', 'Final.Conc', 'sample', 'doseResponses', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set_v1.csv', graph = TRUE, viability = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/FIMM")


#oslo old
for_dss <- inner_join(common_drugs, enserink_breeze_input, by=c("enserink_drug_name"="drug"))
for_dss <- subset(for_dss, select=c('pubchem_drug_name', 'dose', 'Patient.ID', 'y_true'))
colnames(for_dss)[colnames(for_dss) == "Patient.ID"] <- "SCREEN_NAME"
colnames(for_dss)[colnames(for_dss) == "pubchem_drug_name"] <- "DRUG_NAME"
colnames(for_dss)[colnames(for_dss) == "dose"] <- "CONCENTRATION_nM"
colnames(for_dss)[colnames(for_dss) == "y_true"] <- "CELL_VIABILITY"
for_dss$SCREEN_NAME <- as.character(for_dss$SCREEN_NAME)
for_dss$CONCENTRATION_nM <- as.double(for_dss$CONCENTRATION_nM)
#for_dss <- for_dss %>%mutate(CONCENTRATION_nM = log(CONCENTRATION_nM))
options(scipen = 999)
for_dss$CELL_VIABILITY <-1 - for_dss$CELL_VIABILITY/100
df_dose_response <- DOSE_RESPONSE_PROCESS(for_dss, viability = TRUE)
df_dose_response[[2]]$conc
setwd("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/")
df_metrics_enserink <- CALC_METRICS(df_dose_response[[1]], df_dose_response[[2]], graph = FALSE)
setwd("/Users/katarinawilloch/")

write_csv(df_metrics_enserink, '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink.csv')
x <- for_dss
colnames(x)[colnames(x) == "CONCENTRATION_nM"] <- "CONCENTRATION"
colnames(x)[colnames(x) == "CELL_VIABILITY"] <- "PERCENT_INHIBITION"
write_csv(x, '~/Desktop/UiO/Project 1/Data/Breeze/enserink_input_r.csv')

###################################################### calc dsss


#----FIMM calculation old----

fimm_raw_response$patient_id
summary_table <- fimm_raw_response %>% summarize(count = n_distinct(patient_id))

fimm_supplemental <- read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_normalised_reponse_raw_supplemental.csv')
fimm_supplemental_common_drugs <- inner_join(common_drugs, fimm_supplemental, by=c("fimm_drug_name"="Chemical_compound"))
colnames(fimm_supplemental_common_drugs)[colnames(fimm_supplemental_common_drugs) == "pubchem_drug_name"] <- "drug"
fimm_supplemental_common_drugs$Normalised_reponse <- fimm_supplemental_common_drugs$Normalised_reponse/100
fimm_supplemental_common_drugs <- subset(fimm_supplemental_common_drugs, drug != 'TG100-115')
fimm_supplemental_common_drugs <- fimm_supplemental_common_drugs[!grepl('Healthy', fimm_supplemental_common_drugs$Sample_ID), ]
bad_combos <- fimm_supplemental_common_drugs %>%
  filter(is.na(Normalised_reponse)) %>%
  dplyr::select(Sample_ID, drug) %>%
  distinct()

fimm_supplemental_common_drugs_clean <- fimm_supplemental_common_drugs %>%
  anti_join(bad_combos, by = c("Sample_ID", "drug"))

fimm_supplemental_common_drugs_clean[!complete.cases(fimm_supplemental_common_drugs_clean), ]

fimm_supplemental_common_drugs_clean <- fimm_supplemental_common_drugs_clean %>% group_by(Sample_ID, drug) %>% dplyr::summarize(count=n(), .groups = "drop") %>% filter(count>3) %>% left_join(fimm_supplemental_common_drugs_clean, by = c("Sample_ID", "drug"))
unique(fimm_supplemental_common_drugs_clean[,"Sample_ID"])

fimm_supplemental_common_drugs_clean <- subset(fimm_supplemental_common_drugs_clean, Sample_ID !="AML_001_01")

for(s in unique(fimm_supplemental_common_drugs_clean$Sample_ID)){
  filter <- fimm_supplemental_common_drugs_clean %>%
    filter(Sample_ID %in% s)
  print(s)
  for (d in unique(filter$drug)){
    filter <- filter %>%
      filter(drug %in% d)
    print(d)
    fimm_dss <- calc_dss(filter, 'drug', 'CONCENTRATION', 'Sample_ID', 'Normalised_reponse', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set_supplemental.csv', graph = TRUE, viability = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/FIMM")
    
  }
}
filter <- fimm_supplemental_common_drugs_clean[1:2,]
unique_ids <- unique(fimm_supplemental_common_drugs_clean$Sample_ID)
length(unique(fimm_supplemental_common_drugs_clean$Sample_ID))
next_20_ids <- unique_ids[1:1]
filter <- fimm_supplemental_common_drugs_clean %>%
  filter(Sample_ID %in% next_20_ids)

unique_ids <- unique(filter$drug)
length(unique(filter$drug))
next_20_ids <- unique_ids[1:1]
filter <- filter %>%
  filter(drug %in% next_20_ids)

fimm_dss <- calc_dss(filter, 'drug', 'CONCENTRATION', 'Sample_ID', 'Normalised_reponse', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set_supplemental.csv', graph = TRUE, viability = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/FIMM")
write_csv(fimm_dss, '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set_supplemental.csv')
as.data.frame(fimm_dss)
filter %>% group_by(Sample_ID, drug) %>% dplyr::summarize(count = n())
combo <- cbind(fimm_dss[1], fimm_dss[2], fimm_dss[3], fimm_dss[4])
test$dose_responses_grouped
unique(test$dose_responses_grouped$drug) %in% unique(filter$drug) 
length(test$dose_responses_grouped)





#----Karolinska calculation ----
karolinska_raw_response = read.csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/karolinska_normalised_reponse_raw.csv')
#karolinska_raw_response <- subset(karolinska_raw_response, DRUG_NAME != 'Trifluridine' & Patient.num != 'AML_007')
#karolinska_raw_response <- subset(karolinska_raw_response, DRUG_NAME != 'Saracatinib' & Patient.num != 'AML_002' & sample == 'fresh')
length(unique(karolinska_raw_response$Patient.num))
karolinska_raw_response$Normalised_reponse = 1 - karolinska_raw_response$Normalised_reponse/100
calc_dss(subset(karolinska_raw_response, sample == 'fresh'), 'DRUG_NAME', 'CONCENTRATION', 'Patient.num', 'Normalised_reponse', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_karolinska_full_drug_set_fresh.csv', graph = TRUE, viability = FALSE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/Karolinska/fresh")
calc_dss(subset(karolinska_raw_response, sample == 'frozen'), 'DRUG_NAME', 'CONCENTRATION', 'Patient.num', 'Normalised_reponse', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_karolinska_full_drug_set_frosen.csv', graph = TRUE, viability = FALSE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/Karolinska/Frosen")

subset(karolinska_raw_response, sample == 'fresh' & Patient.num == 'AML_001' & DRUG_NAME == '4-hydroxytamoxifen')

#Calculating AUC scores
Log10.DSS <- function(conc, resp){
  resp = resp[order(conc)]
  conc = conc[order(conc)]
  conc = log10(conc)
  
  A = 0
  for(j in 2:length(conc)){
    a <- (resp[j]) *(conc[j]-conc[j-1])/(max(conc)-min(conc))
    A = A + a
    rm(a,j)
  }
  A
}

df.rAUC <- for_dss %>%
  group_by(SCREEN_NAME, DRUG_NAME) %>%
  summarise(AUC = Log10.DSS(CONCENTRATION_nM, CELL_VIABILITY)) %>%
  as.data.frame()
