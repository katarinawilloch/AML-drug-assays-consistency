# -------------------------------------------------------- Script overview
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
# -------------------------------------------------------- Script overview

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

#Run functions from files DSS.R and HelperFunctions.R
source('~/Desktop/UiO/DSS-v2.0/DSS_V2.R')
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
initial_cleansing <- "~/Desktop/UiO/Project 1/Data/Second run/"
dss_output_locatiuon <- "~/Desktop/UiO/Project 1/Data/Response scores/"

#import common drug names dataset
common_drugs <- read_csv(paste0(initial_cleansing,"common_drugnames_pubchem.csv"))
write_csv(common_drugs, paste0(dss_output_locatiuon,'test.csv'))


#----DSS calculation----
##----Beat AML calculation ----
print('BEAT AML')
#import beat aml datasets
beat_aml_raw_response <- read.csv(paste0(initial_cleansing, 'beat_aml_inhibitor_cleaned.csv'))
#reformatting to get percent viability
beat_aml_raw_response$avg2_response <- beat_aml_raw_response$avg2_response/100
#making unique sample ids
beat_aml_raw_response <- beat_aml_raw_response %>%
  mutate(patient_id = paste(dbgap_subject_id, dbgap_dnaseq_sample,dbgap_rnaseq_sample, sep = "_")) %>% as.data.frame()
#transforming concentration from um to nm
beat_aml_raw_response$well_concentration <- beat_aml_raw_response$well_concentration*1000
#getting only 47 common drugs
beat_aml_raw_response_common_drugs <- inner_join(common_drugs, beat_aml_raw_response, by=c("beat_aml_drug_name"="inhibitor"))
#using function calc_dss to calculate dss values
#beat_aml_dss <- calc_dss(beat_aml_raw_response_common_drugs, 'pubchem_drug_name', 'well_concentration', 'patient_id', 'avg2_response', paste0(dss_output_locatiuon,'dss_github_beat_aml_full_drug_set_v1.csv'), graph = TRUE, viability = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/BeatAML")
beat_aml_dss <- read_csv(paste0(dss_output_locatiuon,'dss_github_beat_aml_full_drug_set_v1.csv'))
print(beat_aml_dss)
#annotating with what study the data is from
beat_aml_dss$lab <- 'Beat AML'
#annotating with sample type for later merging
beat_aml_dss$sample <- 'fresh'

##----Helsinki calculation ----
print('HELSINKI')
#import helsinki dataset
key <- read_excel(paste0(initial_cleansing,'key.xlsx'))
fimm_raw_response <- read.csv(paste0(initial_cleansing,'fimm_raw_dose_responses.csv.csv'))
fimm_raw_response <- left_join(fimm_raw_response, key, by = "sample")
fimm_raw_response <- subset(fimm_raw_response, !is.na(sample_id))
#transform concentrations from ... to nm
#fimm_raw_response$Final.Conc <- fimm_raw_response$Final.Conc*10
#transforming viability values to percent viability
fimm_raw_response$doseResponses <- fimm_raw_response$doseResponses/100
#removing drug 'TG100-115'
fimm_raw_response_df <- subset(fimm_raw_response, drug != 'TG100-115')
#getting only 47 common drugs 
fimm_raw_response_common_drugs <- inner_join(common_drugs, fimm_raw_response_df, by=c("fimm_drug_name"="drug"))
#using function calc_dss to get dss values
#helsinki_dss <- calc_dss(fimm_raw_response_common_drugs, 'pubchem_drug_name', 'Final.Conc', 'sample_id', 'doseResponses', paste0(dss_output_locatiuon,'dss_github_fimm_full_drug_set_v1.csv'), graph = TRUE, viability = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/FIMM")
helsinki_dss <- read_csv(paste0(dss_output_locatiuon,'dss_github_fimm_full_drug_set_v1.csv'))
print(helsinki_dss)
#annotating with what study the data is from
helsinki_dss$lab <- 'Helsinki'
#annotating with sample type for later merging will be updated with actual info when merging with annotation file
helsinki_dss$sample <- NA
helsinki_dss <- left_join(helsinki_dss, key, by = c("Patient.num" = "sample"))
helsinki_dss$Patient.num <- helsinki_dss$sample_id
helsinki_dss$sample_id <- NULL
helsinki_dss <- subset(helsinki_dss, !is.na(Patient.num))

##----Oslo calculation ----
print('OSLO')
#importing the oslo dataset
enserink_breeze_input <- read_csv("/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_dose_response.csv")
#transforming concentration from m to nm
enserink_breeze_input$dose <- enserink_breeze_input$dose*10^9
#transforming inhibition values to percent viability
enserink_breeze_input$response <-1 - enserink_breeze_input$y_true/100
#getting only 47 common drugs 
oslo_raw_response_common_drugs <- inner_join(common_drugs, enserink_breeze_input, by=c("enserink_drug_name"="drug"))
#using function calc_dss to get dss values
#oslo_dss <- calc_dss(oslo_raw_response_common_drugs, 'pubchem_drug_name', 'dose', 'Patient.ID', 'response', paste0(dss_output_locatiuon,'dss_github_enserink_full_drug_set_testing.csv'), viability = TRUE, graph = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/Enserink")
oslo_dss <- read_csv(paste0(dss_output_locatiuon,'dss_github_enserink_full_drug_set_testing.csv'))
print(oslo_dss)
#annotating with what study the data is from
oslo_dss$lab <- 'Oslo'
#annotating with sample type for later merging will be updated with actual info when merging with annotation file
oslo_dss$sample <- NA

##----Karolinska calculation ----
print('KAROLINSKA')
#import karolinska dataset
karolinska_raw_response = read.csv(paste0(initial_cleansing,'karolinska_normalised_reponse_raw.csv'))
#karolinska_raw_response <- subset(karolinska_raw_response, DRUG_NAME != 'Trifluridine' & Patient.num != 'AML_007')
#karolinska_raw_response <- subset(karolinska_raw_response, DRUG_NAME != 'Saracatinib' & Patient.num != 'AML_002' & sample == 'fresh')
#transforming inhibition values to percent viability
karolinska_raw_response$Normalised_reponse = 1 - karolinska_raw_response$Normalised_reponse/100
#getting only 47 common drugs 
karolinska_raw_response_common_drugs <- inner_join(common_drugs, karolinska_raw_response, by=c("karolinska_drug_name"="DRUG_NAME"))
#using function calc_dss to get dss values for the fresh samples
#karolinska_fresh_dss <- calc_dss(subset(karolinska_raw_response_common_drugs, sample == 'fresh'), 'pubchem_drug_name', 'CONCENTRATION', 'Patient.num', 'Normalised_reponse', paste0(dss_output_locatiuon,'dss_github_karolinska_full_drug_set_fresh.csv'), graph = TRUE, viability = FALSE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/Karolinska/fresh")
karolinska_fresh_dss <- read_csv(paste0(dss_output_locatiuon,'dss_github_karolinska_full_drug_set_fresh.csv'))
print(karolinska_fresh_dss)
#using function calc_dss to get dss values for the frozen samples
#karolinska_frozen_dss <- calc_dss(subset(karolinska_raw_response_common_drugs, sample == 'frozen'), 'pubchem_drug_name', 'CONCENTRATION', 'Patient.num', 'Normalised_reponse', paste0(dss_output_locatiuon,'dss_github_karolinska_full_drug_set_frosen.csv'), graph = TRUE, viability = FALSE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/Karolinska/Frosen")
karolinska_frozen_dss <- read_csv(paste0(dss_output_locatiuon,'dss_github_karolinska_full_drug_set_frosen.csv'))
print(karolinska_frozen_dss)
#annotating dfs with fresh and frozen samples
karolinska_fresh_dss$sample <- 'fresh'
karolinska_frozen_dss$sample <- 'frozen'
#combining the fresh and froen dss data
dss_karolinska_github <- rbind(karolinska_fresh_dss, karolinska_frozen_dss)
#annotating the combined data with wich lab the data came from
dss_karolinska_github$lab <- 'Karolinska'
#annotating sample ids with k in front to distinguish the ids from the helsinki data
dss_karolinska_github <- dss_karolinska_github %>% mutate(Patient.num = paste0('k',dss_karolinska_github$Patient.num))
dss_karolinska_github <- dss_karolinska_github %>% mutate(Patient.num = str_split_i(dss_karolinska_github$Patient.num, '_f', 1))


#----Raw AUC calculation ----
Log10.DSS <- function(conc, resp){
  resp = resp[order(conc)]
  conc = conc[order(conc)]
  conc = log10(conc)
  print(resp)
  print(conc)
  
  A = 0
  for(j in 2:length(conc)){
    a <- (resp[j]) *(conc[j]-conc[j-1])/(max(conc)-min(conc))
    A = A + a
    rm(a,j)
  }
  print(A)
  A
}
##----Beat AML calculation ----
auc_beat_aml <- beat_aml_raw_response_common_drugs %>% 
  group_by(patient_id, pubchem_drug_name) %>% 
  summarise(auc_a = Log10.DSS(well_concentration, avg2_response), .groups = 'drop') 

ggplot(auc_beat_aml, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")

##----Oslo calculation ----
auc_enserink <- oslo_raw_response_common_drugs %>% 
  group_by(Patient.ID, pubchem_drug_name) %>% 
  summarise(auc_a = Log10.DSS(dose, response), .groups = 'drop') 

ggplot(auc_enserink, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")

##----Helsinki calculation ----
auc_FIMM <- fimm_raw_response_common_drugs %>% 
  group_by(sample_id, pubchem_drug_name) %>% 
  summarise(auc_a = Log10.DSS(Final.Conc, doseResponses), .groups = 'drop') 
ggplot(auc_FIMM, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")

##----Karolinska calculation ----
auc_karolinska <- karolinska_raw_response_common_drugs %>% 
  group_by(Patient.num, sample, pubchem_drug_name) %>% 
  summarise(auc_a = Log10.DSS(CONCENTRATION, Normalised_reponse), .groups = 'drop') 

ggplot(auc_karolinska, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")

#----Merging files----
#Merge 1 - dss for all labs 
dss_values <- rbind(beat_aml_dss, helsinki_dss, oslo_dss, dss_karolinska_github)
print(dss_values)
write_csv(dss_values, paste0(dss_output_locatiuon, 'dss_values_all_labs.csv'))

#Get auc dfs ready for merging
colnames(auc_beat_aml)[colnames(auc_beat_aml) == "patient_id"] <- "Patient.num"
colnames(auc_FIMM)[colnames(auc_FIMM) == "sample_id"] <- "Patient.num"
colnames(auc_enserink)[colnames(auc_enserink) == "Patient.num"] <- "Patient.ID"
colnames(auc_enserink)[colnames(auc_enserink) == "Patient.ID"] <- "Patient.num"
auc_beat_aml$sample <- 'fresh'
auc_FIMM$sample <- NA
auc_enserink$sample <- NA
#annotating sample ids with k in front to distinguish the ids from the helsinki data
auc_karolinska <- auc_karolinska %>% mutate(Patient.num = paste0('k',auc_karolinska$Patient.num))
auc_karolinska <- auc_karolinska %>% mutate(Patient.num = str_split_i(auc_karolinska$Patient.num, '_f', 1))


#Merge 2 - raw auc for all labs 
all_labs_auc <- rbind(auc_beat_aml, auc_karolinska, auc_FIMM, auc_enserink)
colnames(all_labs_auc)[colnames(all_labs_auc) == "pubchem_drug_name"] <- "drug"
write_csv(all_labs_auc, paste0(dss_output_locatiuon,'auc_calculation_all_labs.csv'))

#Merge 3 - dss and raw auc 
all_response_metrics <- inner_join(dss_values, all_labs_auc, by=c("Patient.num", "drug", "sample"))
all_response_metrics$ID <- NULL

##Save merged file for use in lmm ----
write_csv(all_response_metrics, paste0(dss_output_locatiuon,'all_response_metrics_all_labs.csv'))











# ###################################################### old code
# #oslo old
# for_dss <- inner_join(common_drugs, enserink_breeze_input, by=c("enserink_drug_name"="drug"))
# for_dss <- subset(for_dss, select=c('pubchem_drug_name', 'dose', 'Patient.ID', 'y_true'))
# colnames(for_dss)[colnames(for_dss) == "Patient.ID"] <- "SCREEN_NAME"
# colnames(for_dss)[colnames(for_dss) == "pubchem_drug_name"] <- "DRUG_NAME"
# colnames(for_dss)[colnames(for_dss) == "dose"] <- "CONCENTRATION_nM"
# colnames(for_dss)[colnames(for_dss) == "y_true"] <- "CELL_VIABILITY"
# for_dss$SCREEN_NAME <- as.character(for_dss$SCREEN_NAME)
# for_dss$CONCENTRATION_nM <- as.double(for_dss$CONCENTRATION_nM)
# #for_dss <- for_dss %>%mutate(CONCENTRATION_nM = log(CONCENTRATION_nM))
# options(scipen = 999)
# for_dss$CELL_VIABILITY <-1 - for_dss$CELL_VIABILITY/100
# df_dose_response <- DOSE_RESPONSE_PROCESS(for_dss, viability = TRUE)
# df_dose_response[[2]]$conc
# setwd("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/")
# df_metrics_enserink <- CALC_METRICS(df_dose_response[[1]], df_dose_response[[2]], graph = FALSE)
# setwd("/Users/katarinawilloch/")
# 
# write_csv(df_metrics_enserink, '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink.csv')
# x <- for_dss
# colnames(x)[colnames(x) == "CONCENTRATION_nM"] <- "CONCENTRATION"
# colnames(x)[colnames(x) == "CELL_VIABILITY"] <- "PERCENT_INHIBITION"
# write_csv(x, '~/Desktop/UiO/Project 1/Data/Breeze/enserink_input_r.csv')
# 
# #FIMM calculation old
# 
# fimm_raw_response$patient_id
# summary_table <- fimm_raw_response %>% summarize(count = n_distinct(patient_id))
# 
# fimm_supplemental <- read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_normalised_reponse_raw_supplemental.csv')
# fimm_supplemental_common_drugs <- inner_join(common_drugs, fimm_supplemental, by=c("fimm_drug_name"="Chemical_compound"))
# colnames(fimm_supplemental_common_drugs)[colnames(fimm_supplemental_common_drugs) == "pubchem_drug_name"] <- "drug"
# fimm_supplemental_common_drugs$Normalised_reponse <- fimm_supplemental_common_drugs$Normalised_reponse/100
# fimm_supplemental_common_drugs <- subset(fimm_supplemental_common_drugs, drug != 'TG100-115')
# fimm_supplemental_common_drugs <- fimm_supplemental_common_drugs[!grepl('Healthy', fimm_supplemental_common_drugs$Sample_ID), ]
# fimm_dss <- calc_dss(fimm_supplemental_common_drugs, 'drug', 'CONCENTRATION', 'Sample_ID', 'Normalised_reponse', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set_supplemental.csv', graph = TRUE, viability = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/FIMM/supp")

# bad_combos <- fimm_supplemental_common_drugs %>%
#   filter(is.na(Normalised_reponse)) %>%
#   dplyr::select(Sample_ID, drug) %>%
#   distinct()
# 
# fimm_supplemental_common_drugs_clean <- fimm_supplemental_common_drugs %>%
#   anti_join(bad_combos, by = c("Sample_ID", "drug"))
# 
# fimm_supplemental_common_drugs_clean[!complete.cases(fimm_supplemental_common_drugs_clean), ]
# 
# fimm_supplemental_common_drugs_clean <- fimm_supplemental_common_drugs_clean %>% group_by(Sample_ID, drug) %>% dplyr::summarize(count=n(), .groups = "drop") %>% filter(count>3) %>% left_join(fimm_supplemental_common_drugs_clean, by = c("Sample_ID", "drug"))
# unique(fimm_supplemental_common_drugs_clean[,"Sample_ID"])
# 
# fimm_supplemental_common_drugs_clean <- subset(fimm_supplemental_common_drugs_clean, Sample_ID !="AML_001_01")
# 
# for(s in unique(fimm_supplemental_common_drugs_clean$Sample_ID)){
#   filter <- fimm_supplemental_common_drugs_clean %>%
#     filter(Sample_ID %in% s)
#   print(s)
#   for (d in unique(filter$drug)){
#     filter <- filter %>%
#       filter(drug %in% d)
#     print(d)
#     fimm_dss <- calc_dss(filter, 'drug', 'CONCENTRATION', 'Sample_ID', 'Normalised_reponse', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set_supplemental.csv', graph = TRUE, viability = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/FIMM")
#     
#   }
# }
# filter <- fimm_supplemental_common_drugs_clean[1:2,]
# unique_ids <- unique(fimm_supplemental_common_drugs_clean$Sample_ID)
# length(unique(fimm_supplemental_common_drugs_clean$Sample_ID))
# next_20_ids <- unique_ids[1:1]
# filter <- fimm_supplemental_common_drugs_clean %>%
#   filter(Sample_ID %in% next_20_ids)
# 
# unique_ids <- unique(filter$drug)
# length(unique(filter$drug))
# next_20_ids <- unique_ids[1:1]
# filter <- filter %>%
#   filter(drug %in% next_20_ids)
# 
# fimm_dss <- calc_dss(filter, 'drug', 'CONCENTRATION', 'Sample_ID', 'Normalised_reponse', '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set_supplemental.csv', graph = TRUE, viability = TRUE, set_dir = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/FIMM")
# write_csv(fimm_dss, '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set_supplemental.csv')
# as.data.frame(fimm_dss)
# filter %>% group_by(Sample_ID, drug) %>% dplyr::summarize(count = n())
# combo <- cbind(fimm_dss[1], fimm_dss[2], fimm_dss[3], fimm_dss[4])
# test$dose_responses_grouped
# unique(test$dose_responses_grouped$drug) %in% unique(filter$drug) 
# length(test$dose_responses_grouped)
# 
# 
# 
# 
# 
# 
# 
# #Calculating AUC scores
# Log10.DSS <- function(conc, resp){
#   resp = resp[order(conc)]
#   conc = conc[order(conc)]
#   conc = log10(conc)
#   
#   A = 0
#   for(j in 2:length(conc)){
#     a <- (resp[j]) *(conc[j]-conc[j-1])/(max(conc)-min(conc))
#     A = A + a
#     rm(a,j)
#   }
#   A
# }
# 
# df.rAUC <- for_dss %>%
#   group_by(SCREEN_NAME, DRUG_NAME) %>%
#   summarise(AUC = Log10.DSS(CONCENTRATION_nM, CELL_VIABILITY)) %>%
#   as.data.frame()
# 
# 
# #beat aml old
# #import beat aml datasets
# beat_aml_inhibitor <- read_csv(paste0(initial_cleansing, "beat_aml_inhibitor.csv"))
# beat_aml_inhibitor$...1 <- NULL
# #get only beat aml dataset for the common drug names
# beat_aml_for_dss <- inner_join(common_drugs, beat_aml_inhibitor, by=c("beat_aml_drug_name"="inhibitor"))
# #create distinct sample ids 
# beat_aml_for_dss <- beat_aml_for_dss %>%
#   mutate(patient_id = paste(dbgap_subject_id, inhibitor_panel,replicate,run_index, sep = "_")) %>% as.data.frame()
# # Display the dataframe
# duplicates <- beat_aml_for_dss %>%
#   group_by(patient_id, pubchem_drug_name) %>%
#   filter(n() > 2) %>%
#   ungroup()
# # Show the duplicate values
# print(duplicates)
# print(beat_aml_for_dss)
# beat_aml_for_dss <- subset(beat_aml_for_dss, select=c('pubchem_drug_name', 'well_concentration', 'patient_id', 'normalized_viability'))
# colnames(beat_aml_for_dss)[colnames(beat_aml_for_dss) == "pubchem_drug_name"] <- "DRUG_NAME"
# colnames(beat_aml_for_dss)[colnames(beat_aml_for_dss) == "well_concentration"] <- "CONCENTRATION_nM"
# colnames(beat_aml_for_dss)[colnames(beat_aml_for_dss) == "patient_id"] <- "SCREEN_NAME"
# colnames(beat_aml_for_dss)[colnames(beat_aml_for_dss) == "normalized_viability"] <- "CELL_VIABILITY"
# beat_aml_for_dss$SCREEN_NAME <- as.character(beat_aml_for_dss$SCREEN_NAME)
# beat_aml_for_dss
# 
# # Display the dataframe
# duplicates <- beat_aml_for_dss %>%
#   arrange(SCREEN_NAME, DRUG_NAME) %>%
#   group_by(SCREEN_NAME, DRUG_NAME) %>%
#   filter(n() > 2) %>%
#   ungroup() 
# 
# # Show the duplicate values
# print(duplicates)
# 
# beat_aml_for_dss$CELL_VIABILITY <- beat_aml_for_dss$CELL_VIABILITY/100
# df_dose_response_beat_aml <- DOSE_RESPONSE_PROCESS(beat_aml_for_dss, viability = TRUE)
# setwd("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/BeatAML/test")
# df_metrics_beat_aml <- CALC_METRICS(df_dose_response_beat_aml[[1]], df_dose_response_beat_aml[[2]], graph = TRUE)
# setwd("/Users/katarinawilloch/")
# 
# write_csv(df_metrics_beat_aml, '~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_beat_aml.csv')
