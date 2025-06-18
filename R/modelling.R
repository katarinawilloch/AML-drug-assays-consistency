# -------------------------------------------------------- Script overview
# Script Name: modelling.R
# Author: Katarina Willoch
# Date: 2025-06-01
#
# Description:
# This script reads in dose response metric scores DSS1, DSS2, DSS3, IC50, AUC, rAUC,
# and reads in experimental variable data 
# and clinical sample information
# and merges these datasets
# and generates linear mixed models for each experimental variable
# 
#
# Sections:
#   1. Load libraries and data
#   2. Merge datasets
#   3. PPCA generation for DSS2
#   4. Heatmap and PPCA generation for all metrics
#   5. ggbetweenstats for DSS2
# -------------------------------------------------------- Script overview

#Install and import libraries ----
library(lme4)
library(tidyr)
library(lattice)
library(lmerTest)
library(ggridges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(viridis)
library(hrbrthemes)
library(ComplexHeatmap)
library(readr)
library(grid)
library(car)
library(ggeffects)
library(effects)
library(merTools)
library(cowplot)
library(ez)
library(performance)
library(easystats)
library(lme4)
library(broom.mixed)
library(ggplot2)
library(sjPlot)
library(xtable)
library(htmltools)
library(forcats)
library(car)
library(influence.ME)
library(visdat)
library(see)
library(sva)
library(matrixStats)
library(ComplexHeatmap)
library(FactoMineR)
library(influence.ME)
library(formattable)
library(ggplotify)
library(plotly)
library(flexplot)
library(readxl)
library(MASS)
library(forestploter)
require(flexplot)
library(patchwork)



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace(c("sva", "pcaMethods", "ChAMP"), quietly = TRUE))
  BiocManager::install(c("sva", "pcaMethods", "ChAMP"), force = TRUE)
lapply(c("matrixStats","dplyr","reshape","reshape2", "scales", "drc", "caTools","ggplot2", "data.table", "stringr","MESS", "BiocManager","svMisc", "egg", "pheatmap", "sva", "pcaMethods"), library, character.only = T)

source('~/Desktop/UiO/Project 1/Code/R/Linear_mixed_model_resport.R')
source('~/Desktop/UiO/Project 1/Code/R/model_functions.R')

#Figure output location
figure_output <- '~/Desktop/UiO/Project 1/Figures/New karolinska data/'
initial_cleansing <- "~/Desktop/UiO/Project 1/Data/Second run/"

#Import datset containing all metric responses
all_response_metrics <- read_csv('~/Desktop/UiO/Project 1/Data/Response scores/all_response_metrics_all_labs.csv')
#Removing data without sample id
all_response_metrics <- all_response_metrics[!is.na(all_response_metrics$Patient.num),]
#Renaming drugs to first letter capital and rest lower case
all_response_metrics <- all_response_metrics %>% mutate(drug = ifelse(drug == 'PONATINIB', 'Ponatinib', drug))
all_response_metrics <- all_response_metrics %>% mutate(drug = ifelse(drug == 'MOMELOTINIB', 'Momelotinib', drug))
all_response_metrics <- all_response_metrics %>% mutate(Patient.num = ifelse(lab == 'Karolinska', paste0(Patient.num, sample), Patient.num))
all_response_metrics$drug <- tools::toTitleCase(all_response_metrics$drug)
#Convert the column to a factor with the specific order
specific_order <- c("Beat AML", "Oslo", "Helsinki", "Karolinska")
all_response_metrics$lab <- factor(all_response_metrics$lab, levels = specific_order)
#Creating new variable for negative log10 IC50
all_response_metrics$log10IC50 <- -log10(all_response_metrics$IC50)


#Merging clinical sample information with calculated metrics for each lab ----
##Karolinska
karolinska_sample_info <- read_csv(paste0(initial_cleansing,'karolinska_sample_annotation.csv'))
###Rename column headers
karolinska_sample_info <- karolinska_sample_info %>%
  dplyr::rename(
    "Patient.num" = "Running_ID",
    "Disease.status" = "diagnosis",
    "Tissue" = "Sample_type", 
    "sample" = "cohort"
  )
###Rename Tissue values
karolinska_sample_info$Tissue <- gsub('BM_blank', 'Bone marrow', karolinska_sample_info$Tissue)
karolinska_sample_info$Tissue <- gsub('BM', 'Bone marrow', karolinska_sample_info$Tissue)
karolinska_sample_info$Tissue <- gsub('PB', 'Blood', karolinska_sample_info$Tissue)
###Rename Disease.status values
karolinska_sample_info$Disease.status <- gsub('de_novo', 'Diagnosis', karolinska_sample_info$Disease.status)
karolinska_sample_info$Disease.status <- gsub('relapse', 'Relapse', karolinska_sample_info$Disease.status)
###Rename sample values
karolinska_sample_info$sample <- gsub('prospective', 'fresh', karolinska_sample_info$sample)
karolinska_sample_info$sample <- gsub('retrospective', 'frozen', karolinska_sample_info$sample)
###Rename sample ids to distinguish from Helsinki sample ids
karolinska_sample_info <- karolinska_sample_info %>% mutate(Patient.num = paste0('k',karolinska_sample_info$Patient.num,sample))
#Merge Karolinska sample annotation and response metric dataset
all_response_metrics <- left_join(all_response_metrics, karolinska_sample_info[,c('Patient.num', 'Tissue', 'Disease.status')], by=c('Patient.num'))

##Oslo
enserink_sample_annotation <- read_csv(paste0(initial_cleansing,'enserink_lab_sample_annotation.csv'))
###Rename column headers
enserink_sample_annotation <- enserink_sample_annotation %>%
  dplyr::rename(
    "Patient.num" = "Patient",
    "Cells1" = "Cells/well",
    "Plate_reader1" = "Instrument", 
    "sample" = "Fresh or frozen"
  )

###Rename Tissue values
enserink_sample_annotation$Tissue <- gsub('BM', 'Bone marrow', enserink_sample_annotation$Tissue)
enserink_sample_annotation$Tissue <- gsub('PB', 'Blood', enserink_sample_annotation$Tissue)
###Get Disease status values from Patient.num
enserink_sample_annotation <- enserink_sample_annotation %>% mutate(Disease.status = str_split_i(enserink_sample_annotation$Patient.num, '_', 2))
enserink_sample_annotation[is.na(enserink_sample_annotation$Disease.status),]$Disease.status <- 'Diagnosis'
###Rename Disease.status values
enserink_sample_annotation$Disease.status <- gsub('relapse', 'Relapse', enserink_sample_annotation$Disease.status)
enserink_sample_annotation$Disease.status <- gsub('remission2', 'remission', enserink_sample_annotation$Disease.status)
###Rename sample values 
enserink_sample_annotation$sample <- gsub('Fresh', 'fresh', enserink_sample_annotation$sample)
###Categorize number of cells per well per sample
enserink_sample_annotation <- enserink_sample_annotation %>% mutate(Cells1 = case_when(Cells1 > 5000 ~"10000", 
                                                                                       Cells1 == 5000 ~ "5000", 
                                                                                       Cells1 < 5000 ~ NA,
                                                                                       TRUE ~ NA))

all_response_metrics <- left_join(all_response_metrics, enserink_sample_annotation[,c('Patient.num', 'Tissue', 'Disease.status', 'sample')], by = "Patient.num") %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    sample = coalesce(sample.x, sample.y), 
    Disease.status = coalesce(Disease.status.x, Disease.status.y)
  ) %>%
  select(-Tissue.x, -Tissue.y, -sample.x, -sample.y, -Disease.status.x, -Disease.status.y)
##Helsinki
fimm_sample_information <- read.csv(paste0(initial_cleansing,'fimm_assay_dets.csv'))
colnames(fimm_sample_information)[colnames(fimm_sample_information) == "Medium"] <- "Medium1"
fimm_sample_information$Medium1 <- gsub("MCM", 'Mononuclear cell medium', fimm_sample_information$Medium1)
fimm_sample_information$Medium1 <- gsub("CM", 'HS-5 conditioned medium', fimm_sample_information$Medium1)
colnames(fimm_sample_information)[colnames(fimm_sample_information) == "Sample_ID"] <- "Patient.num"
fimm_sample_annotation <- read.csv(paste0(initial_cleansing,'fimm_sample_annotation.csv'))
colnames(fimm_sample_annotation)[colnames(fimm_sample_annotation) == "Sample_ID"] <- "Patient.num"
fimm_sample_annotation$Patient_ID <- gsub("FPM_", "", fimm_sample_annotation$Patient_ID)
fimm_clinical_info <- read_csv(paste0(initial_cleansing,'fimm_clinical_summary.csv'))
fimm_clinical_info$sample <- ifelse(fimm_clinical_info$Site == "Helsinki", 'fresh', 'frozen')
###Merging the halsinki sample annotation files
helsinki_sample_annotation <- merge(unique(fimm_sample_information[,c('Patient_ID', 'Patient.num', 'Medium1')]), fimm_sample_annotation[,c('Patient_ID', 'Patient.num', 'Disease.status', 'Tissue')], by = c("Patient_ID", "Patient.num"), all=TRUE)
helsinki_sample_annotation <- left_join(helsinki_sample_annotation, fimm_clinical_info[,c('Patient_ID','sample')], by = "Patient_ID")
#Merge Helsinki sample information and drug response metrics
all_response_metrics <- left_join(all_response_metrics, helsinki_sample_annotation, by = "Patient.num") %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    sample = coalesce(sample.x, sample.y), 
    Disease.status = coalesce(Disease.status.x, Disease.status.y), 
    Medium1 = coalesce(Medium1.x, Medium1.y)
  ) %>%
  select(-Tissue.x, -Tissue.y, -sample.x, -sample.y, -Disease.status.x, -Disease.status.y, -Medium1.x, -Medium1.y)

subset(all_response_metrics, lab == 'Helsinki')
###Beat AML
beat_aml_clinical_info <- read_csv(paste0(initial_cleansing,'beat_aml_clinical.csv'))
beat_aml_clinical_info <- beat_aml_clinical_info %>% mutate(Patient.num = paste0(beat_aml_clinical_info$dbgap_subject_id, '_', dbgap_dnaseq_sample, '_', dbgap_rnaseq_sample))
colnames(beat_aml_clinical_info)[colnames(beat_aml_clinical_info) == "specimenType"] <- "Tissue"
colnames(beat_aml_clinical_info)[colnames(beat_aml_clinical_info) == "diseaseStageAtSpecimenCollection"] <- "Disease.status"
beat_aml_clinical_info$Tissue <- gsub('Bone Marrow Aspirate', 'Bone marrow', beat_aml_clinical_info$Tissue)
beat_aml_clinical_info$Tissue <- gsub('Peripheral Blood', 'Blood', beat_aml_clinical_info$Tissue)
beat_aml_clinical_info$Disease.status <- gsub('Initial ', '', beat_aml_clinical_info$Disease.status)
beat_aml_clinical_info$Disease.status <- gsub('Residual', 'Refractory', beat_aml_clinical_info$Disease.status)
beat_aml_clinical_info$Disease.status <- gsub('Unknown', NA, beat_aml_clinical_info$Disease.status)
beat_aml_clinical_info$sample <- 'fresh'
#beat_aml_clinical_info <- beat_aml_clinical_info %>% mutate(Patient_ID = sub("_.*", "", Patient.num))
all_response_metrics <- left_join(all_response_metrics, subset(beat_aml_clinical_info, select = c("Patient.num", "Disease.status", "Tissue", "sample")), by = "Patient.num") %>%
  mutate(
    Tissue = coalesce(Tissue.x, Tissue.y),
    sample = coalesce(sample.x, sample.y), 
    Disease.status = coalesce(Disease.status.x, Disease.status.y)
  ) %>%
  select(-Tissue.x, -Tissue.y, -sample.x, -sample.y, -Disease.status.x, -Disease.status.y)

subset(all_response_metrics, lab == 'Beat AML')

#Getting sample counts for table in article----
##Karolinska
karolinska_count_disease_status <- subset(all_response_metrics[all_response_metrics$lab == 'Karolinska',], select = c("Patient.num", "Disease.status")) %>% distinct() %>% group_by(Disease.status) %>% dplyr::summarize(count = n())
print(karolinska_count_disease_status)
karolinska_count_tissue <- subset(all_response_metrics[all_response_metrics$lab == 'Karolinska',], select = c("Patient.num", "Tissue")) %>% distinct() %>% group_by(Tissue) %>% dplyr::summarize(count = n())
print(karolinska_count_tissue)
karolinska_count_sample <- subset(all_response_metrics[all_response_metrics$lab == 'Karolinska',], select = c("Patient.num", "sample")) %>% distinct() %>% group_by(sample) %>% dplyr::summarize(count = n())
print(karolinska_count_sample)
karolinska_count <- subset(all_response_metrics[all_response_metrics$lab == 'Karolinska',], select = c("Patient.num")) %>% distinct() %>% dplyr::summarize(count = n())
print(karolinska_count)
##Oslo
oslo_count_disease_status <- subset(all_response_metrics[all_response_metrics$lab == 'Oslo',], select = c("Patient.num", "Disease.status")) %>% distinct() %>% group_by(Disease.status) %>% dplyr::summarize(count = n())
print(oslo_count_disease_status)
oslo_count_tissue <- subset(all_response_metrics[all_response_metrics$lab == 'Oslo',], select = c("Patient.num", "Tissue")) %>% distinct() %>% group_by(Tissue) %>% dplyr::summarize(count = n())
print(oslo_count_tissue)
oslo_count_sample <- subset(all_response_metrics[all_response_metrics$lab == 'Oslo',], select = c("Patient.num", "sample")) %>% distinct() %>% group_by(sample) %>% dplyr::summarize(count = n())
print(oslo_count_sample)
oslo_count <- subset(all_response_metrics[all_response_metrics$lab == 'Oslo',], select = c("Patient.num")) %>% distinct() %>% dplyr::summarize(count = n())
print(oslo_count)
##Helsinki
helsinki_count_disease_status <- subset(all_response_metrics[all_response_metrics$lab == 'Helsinki',], select = c("Patient.num", "Disease.status")) %>% distinct() %>% group_by(Disease.status) %>% dplyr::summarize(count = n())
print(helsinki_count_disease_status)
helsinki_count_tissue <- subset(all_response_metrics[all_response_metrics$lab == 'Helsinki',], select = c("Patient.num", "Tissue")) %>% distinct() %>% group_by(Tissue) %>% dplyr::summarize(count = n())
print(helsinki_count_tissue)
helsinki_count_sample <- subset(all_response_metrics[all_response_metrics$lab == 'Helsinki',], select = c("Patient.num", "sample")) %>% distinct() %>% group_by(sample) %>% dplyr::summarize(count = n())
print(helsinki_count_sample)
helsinki_count <- subset(all_response_metrics[all_response_metrics$lab == 'Helsinki',], select = c("Patient.num")) %>% distinct() %>% dplyr::summarize(count = n())
print(helsinki_count)
##Beat AML
beat_aml_count_disease_status <- subset(all_response_metrics[all_response_metrics$lab == 'Beat AML',], select = c("Patient.num", "Disease.status")) %>% distinct() %>% group_by(Disease.status) %>% dplyr::summarize(count = n())
print(beat_aml_count_disease_status)
beat_aml_count_tissue <- subset(all_response_metrics[all_response_metrics$lab == 'Beat AML',], select = c("Patient.num", "Tissue")) %>% distinct() %>% group_by(Tissue) %>% dplyr::summarize(count = n())
print(beat_aml_count_tissue)
beat_aml_count_sample <- subset(all_response_metrics[all_response_metrics$lab == 'Beat AML',], select = c("Patient.num", "sample")) %>% distinct() %>% group_by(sample) %>% dplyr::summarize(count = n())
print(beat_aml_count_sample)
beat_aml_count <- subset(all_response_metrics[all_response_metrics$lab == 'Beat AML',], select = c("Patient.num")) %>% distinct() %>% dplyr::summarize(count = n())
print(beat_aml_count)

#Experimental data----
experimental_var <- read_delim(paste0(initial_cleansing,'experimental_variables.csv'), delim=';')
experimental_var[experimental_var$lab == 'Helsinki','medium'] <- NA

#Merging experimental data and response data----
all_response_metrics <- left_join(all_response_metrics, subset(experimental_var, select = -c(sample)), by = 'lab')
#all_response_metrics <- all_response_metrics %>% mutate(nr_of_concentration_points_cat = as.character(nr_of_concentration_points)) %>% as.data.frame()
all_response_metrics <- all_response_metrics %>% mutate(cells = as.character(cells)) %>% as.data.frame()
all_response_metrics <- all_response_metrics %>%
  mutate(medium = if_else(lab == 'Helsinki', Medium1, medium))
all_response_metrics$Cells1 <- as.character(all_response_metrics$Cells1)
all_response_metrics <- all_response_metrics %>%
  mutate(cells = if_else(lab == 'Oslo', Cells1, cells))
all_response_metrics <- all_response_metrics %>%
  mutate(plate_reader = if_else(lab == 'Oslo', Plate_reader1, plate_reader))
all_response_metrics <- all_response_metrics %>%
  mutate(time_until_sample_usage = if_else(sample == 'Frozen', NA, time_until_sample_usage))


all_response_metrics$Disease.status <- gsub('remission', 'Remission', all_response_metrics$Disease.status)
all_response_metrics$Disease.status <- gsub('Unknown', NA, all_response_metrics$Disease.status)
all_response_metrics$sample <- gsub('Cryopreserved', 'frozen', all_response_metrics$sample)
all_response_metrics$Tissue <- gsub('Leukapheresis', 'Blood', all_response_metrics$Tissue)
all_response_metrics$lab <- gsub('BeatAML', 'Beat AML', all_response_metrics$lab)


all_response_metrics$Patient_ID <- all_response_metrics$Patient.num
all_response_metrics <- all_response_metrics %>% dplyr::mutate(Patient_ID = ifelse(lab == 'Oslo' | lab == "Beat AML", gsub("_.*", "", Patient_ID), Patient_ID)) %>%
  dplyr::mutate(Patient_ID = ifelse(lab == 'Helsinki', str_extract(Patient_ID, "^[^_]+_[^_]+"), Patient_ID)) 
all_response_metrics <- all_response_metrics %>% dplyr::mutate(Patient_ID = ifelse(lab == 'Karolinska', gsub('f.*', '', Patient_ID), Patient_ID))

unique(subset(all_response_metrics, lab == 'Karolinska', select = c('Patient.num', 'Patient_ID')))

write_csv(all_response_metrics, paste0(initial_cleansing,'all_metrices_all_labs_with_experimental_var.csv'))


patient_count <- subset(all_response_metrics, select = c("Patient_ID", "lab")) %>%distinct() %>% group_by(lab) %>% dplyr::summarize(count = n())
print(patient_count)


#Modeling----
##Box Cox transformation and Z scaling----
###DSS2----
all_response_metrics$DSS2_pos <- all_response_metrics$DSS2 + 0.01
MASS::boxcox(DSS2_pos ~ time_until_sample_usage, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$DSS2_boxcox <- ((all_response_metrics$DSS2)^0.5 - 1) / 0.5

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$DSS2_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$DSS2_boxcox[all_response_metrics$drug==j])
}

###AUC----
all_response_metrics$AUC_pos <- all_response_metrics$AUC + 0.01
MASS::boxcox(AUC_pos ~ time_until_sample_usage, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$AUC_boxcox <- ((all_response_metrics$AUC)^0.5 - 1) / 0.5

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$AUC_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$AUC_boxcox[all_response_metrics$drug==j])
}

ggplot(all_response_metrics, aes(x=AUC)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=AUC_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=AUC_boxcox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")

###IC50----
all_response_metrics$IC50_pos <- all_response_metrics$IC50 + 0.01
MASS::boxcox(IC50_pos ~ time_until_sample_usage, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))


all_response_metrics <- all_response_metrics %>% mutate(IC50_boxcox = ifelse(-log10(all_response_metrics$IC50) == 0, 0, -log10(all_response_metrics$IC50)))

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$IC50_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$IC50_boxcox[all_response_metrics$drug==j])
}

ggplot(all_response_metrics, aes(x=IC50)) +
  geom_density(color="darkblue", fill="lightblue") + 
  theme_minimal() + 
  ylab("Density") + 
  xlab("IC50")

ggplot(all_response_metrics, aes(x=IC50_boxcox)) +
  geom_density(color="darkblue", fill="lightblue") + theme_minimal() + 
  ylab("Density") + 
  xlab("BoxCox transformed IC50")

ggplot(all_response_metrics, aes(x=IC50_boxcox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue") +
  theme_minimal() + 
  ylab("Density") + 
  xlab("BoxCox transformed z-Scaled IC50")


###DSS1----
all_response_metrics$DSS1_pos <- all_response_metrics$DSS1 + 0.01
MASS::boxcox(DSS1_pos ~ lab, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$DSS1_boxcox <-  ((all_response_metrics$DSS1)^0.5 - 1) / 0.5

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$DSS1_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$DSS1_boxcox[all_response_metrics$drug==j])
}


ggplot(all_response_metrics, aes(x=DSS1_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=DSS1_boxcox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")

###DSS3----
all_response_metrics$DSS3_pos <- all_response_metrics$DSS3 + 0.01
MASS::boxcox(DSS3_pos ~ lab, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$DSS3_boxcox <- ((all_response_metrics$DSS3)^0.5 - 1) / 0.5

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$DSS3_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$DSS3_boxcox[all_response_metrics$drug==j])
}


ggplot(all_response_metrics, aes(x=DSS3_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=DSS3_boxcox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")

###Raw AUC----
all_response_metrics$auc_a_pos <- all_response_metrics$auc_a - min(all_response_metrics$auc_a) + 0.1
MASS::boxcox(auc_a_pos ~ lab, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$auc_a_boxcox <- ((all_response_metrics$auc_a)^0.5 - 1) / 0.5

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$auc_a_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$auc_a_boxcox[all_response_metrics$drug==j])
}

ggplot(all_response_metrics, aes(x=auc_a_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=auc_a_boxcox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")


##LMM for each experimental variable----
# List of models and their names
variables <- c("time_until_sample_usage", "medium", "cells", "positive_control", "centrifugation_procedure", "plate_reader")
metrics <- c("IC50", "DSS1", "DSS2", "DSS3", "auc_a", "AUC")
models <- list()
for (m in metrics[1:6]){
  for (v in variables){
    model_name_plot <- case_when(v == "time_until_sample_usage" ~ "time until usage", 
                                 v == "microenvironmental_stimuli" ~ "microenvironmental stimuli", 
                                 v == "positive_control" ~ "positive control-doses-readout-cell counting", 
                                 v == "centrifugation_procedure" ~ "centrifugation procedure", 
                                 v == "plate_reader" ~ "plate reader",
                                 .default = v)
    print(m)
    print(v)
    m_box <- paste0(m, "_boxcox_sclaed2")
    # Construct formula as a string
    formula_str <- paste0(m_box, " ~ ", v, " + (1|Patient.num)", " + (",v,"|drug)")
    # Convert string to formula
    f <- as.formula(formula_str)
    model_name <- paste0(m, "_by_", v)
    models[[model_name]] <- lmer(f,all_response_metrics)
    #model_diagnostics(models[[model_name]], paste0("/Users/katarinawilloch/Desktop/UiO/Project 1/Code/Diagnostic plots/", m, "/", v, "/"), model_name = model_name_plot, metric = m) 
    print(summary(models[[model_name]]))
  }
}


##Forestplot DSS2----
dss2_cols <- grep("DSS2", names(models), value = TRUE)
dss2_models <- unlist(models[dss2_cols], use.names = FALSE)
#Get confidence intervals og p values from function in model_functions.R
df_re_intercept_models <- model_results(dss2_models, data_frame_org = all_response_metrics)
#Rename random_effect_term values
df_re_intercept_models$Random_Effect_Term <- '(1 | Patient.num) + (1 + Experimental var | drug)'
#Rename df headers
colnames(df_re_intercept_models) <- gsub("_", " ", colnames(df_re_intercept_models))     
colnames(df_re_intercept_models) <- sapply(colnames(df_re_intercept_models), tools::toTitleCase)
df_re_intercept_models <- df_re_intercept_models %>% dplyr::rename(
  `Experimental Variable` = `Fixed Effect Term`,
  `Reference Group` = `Ref Group`
)
#If corrected p value is bigger than 1 return 1
df_re_intercept_models <- df_re_intercept_models %>% mutate(`Corrected p Value` = ifelse(df_re_intercept_models$`p Value`*df_re_intercept_models$Count >= 1, 1, df_re_intercept_models$`p Value`*df_re_intercept_models$Count))
#Rename Reference group values
df_re_intercept_models <- df_re_intercept_models %>% mutate(`Reference Condition` = case_when(`Reference Group` ==  "a drug combination of flavopiridol, staurosporine and velcade" ~ "Positive Control: drug combination + Doses: 7 \n+ Readout: CellTiter96     ", 
                                                                                              `Reference Group` ==  "1-2h after receiving" ~ "Time until sample usage <2h",
                                                                                              `Reference Group` ==  "HS-5 conditioned medium" ~ "HS-5 CM",
                                                                                              `Reference Group` ==  "HS-5 CM" ~ "HS-5 CM",
                                                                                              `Reference Group` ==  "Ficoll-Paque centrifugation" ~ "Ficoll-Paque",
                                                                                              `Reference Group` ==  "10000" ~ "10000",
                                                                                              `Reference Group` ==  "Countess" ~ "Countess",
                                                                                              `Reference Group` ==  "Biotek Synergy 2" ~ "Biotek Synergy 2",
                                                                                              TRUE ~ `Reference Group`))

#Forest plot theme
tm <- forest_theme(base_size =9,
                   arrow_type = "closed",
                   base_family = "Arial",
                   footnote_gp = gpar(col = "black", cex = 1.0, size = 12, family = "Arial"), 
                   line_size = 0.5, 
                   align = "center", 
                   footnote_parse = F, 
                   xlab_gp = gpar(fontsize = 14, fontfamily = "Arial", cex=1),
                   xaxis_gp = gpar(fontsize = 14, fontfamily = "Arial", cex=1), 
                   xlab_adjust = "center")
#For bars in plot
df_re_intercept_models$` ` <- paste(rep(" ", 30), collapse = " ")
#Order by coefficient
df_re_intercept_models <- df_re_intercept_models[order(df_re_intercept_models$`Fixed Effect Coefficient`),]
#Round corrected p value to 4 decimal points
df_re_intercept_models$`p-valueᵃ` <- round(df_re_intercept_models$`Corrected p Value`, 4)
df_re_intercept_models$`p-valueᵃ` <- format(df_re_intercept_models$`p-valueᵃ`, scientific = FALSE, trim = TRUE)
df_re_intercept_models <- df_re_intercept_models %>% mutate(`p-valueᵃ` = ifelse(df_re_intercept_models$`p-valueᵃ` == '0.0000', '<0.0001', as.character(df_re_intercept_models$`p-valueᵃ`)))

p <- forest(df_re_intercept_models[,c('Experimental Variable','Reference Condition', ' ', 'p-valueᵃ')],
            est = df_re_intercept_models$`Fixed Effect Coefficient`,
            lower = df_re_intercept_models$Lower, 
            upper = df_re_intercept_models$Upper,
            sizes = 1,
            ci_column = 3,
            ref_line = 0,
            grid = F,
            #arrow_lab = c("Lower DSS2 than reference group", "Higher DSS2 than refernce group"),
            xlim = c(-1.1, 1.1),
            #ticks_at = c(-1, -0.5, 0, 0.5, 1),
            #footnote = "\n\nᵃBonferroni corrected",
            theme = tm, 
            xlab = expression('Change in Scaled DSS'[2]),
            legend_gp = gpar(fontsize = 18, fontfamily = "Arial", cex = 1),
            xaxis_gp = gpar(fontsize = 28, fontfamily = "Arial", cex=1)
) 

print(p)

###Save forestplot DSS2----
# Define the physical dimensions (cm) and resolution
output_width_cm <- 16.4#14.4 #9.5
output_height_cm <- 16#14 #10
dpi <- 300  # Resolution in dots per inch

# Convert cm to inches (1 inch = 2.54 cm)
output_width_in <- output_width_cm / 2.54
output_height_in <- output_height_cm / 2.54

# Convert inches to pixels for the PNG device
output_width_px <- output_width_in * dpi
output_height_px <- output_height_in * dpi

# Calculate scaling factors for the gtable
scale_width <- output_width_cm / 10  # Base width adjustment (10 cm as reference)
scale_height <- output_height_cm / 10  # Base height adjustment (10 cm as reference)

p <- edit_plot(p, gp = gpar(cex=1.7, fontfamily="Arial")) #1.4
p <- edit_plot(p, part = "header", gp = gpar(cex=1.7, fontfamily="Arial"))

p <- edit_plot(p, col = 4, part = "header",
               which = "text",
               hjust = unit(0.5, "npc"),
               x = unit(0.5, "npc"))
p <- edit_plot(p, col = 4, part = "body",
               which = "text",
               hjust = unit(0.5, "npc"),
               x = unit(0.5, "npc"))


# Scale the gtable layout
scaled_p <- gtable::gtable_filter(p, pattern = ".*", trim = TRUE)  # Keep all grobs
scaled_p$widths <- scaled_p$widths * scale_width
scaled_p$heights <- scaled_p$heights * scale_height

png(paste0(figure_output,'Forest_plot_DSS2.png'), height = 16, width = 46, unit = "cm", res = 300)
grid.newpage()
grid.draw(scaled_p)
dev.off()

##Forestplot all metrics----
#Get confidence intervals og p values from function in model_functions.R
df_lmm_results_all_metrics <- model_results(models, data_frame_org = all_response_metrics)
#If corrected p value is bigger than 1 return 1
df_lmm_results_all_metrics <- df_lmm_results_all_metrics %>% mutate(corrected_p_value = ifelse(df_lmm_results_all_metrics$p_value * df_lmm_results_all_metrics$count >= 1, 1, df_lmm_results_all_metrics$p_value * df_lmm_results_all_metrics$count))
#Round corrected p value to 4 decimal points
df_lmm_results_all_metrics$corrected_p_value <- round(df_lmm_results_all_metrics$corrected_p_value, 4)
df_lmm_results_all_metrics$corrected_p_value <- format(df_lmm_results_all_metrics$corrected_p_value, scientific = FALSE, trim = TRUE)
#If corrected p value is smaller than 0.0001 return <0.0001
df_lmm_results_all_metrics <- df_lmm_results_all_metrics %>% mutate(corrected_p_value = ifelse(df_lmm_results_all_metrics$corrected_p_value == '0.0000', '<0.0001', as.character(df_lmm_results_all_metrics$corrected_p_value)))
#Rename reference group values
df_lmm_results_all_metrics <- df_lmm_results_all_metrics %>% mutate(ref_group = case_when(ref_group ==  "a drug combination of flavopiridol, staurosporine and velcade" ~ "Positive Control: drug combination + Doses: 7 \n+ Readout: CellTiter96     ", 
                                          ref_group ==  "1-2h after receiving" ~ "Time until sample usage < 2h",
                                          ref_group ==  "HS-5 conditioned medium" ~ "HS-5 CM",
                                          ref_group ==  "HS-5 CM" ~ "HS-5 CM",
                                          ref_group ==  "Ficoll-Paque centrifugation" ~ "Ficoll-Paque",
                                          ref_group ==  "10000" ~ "10000",
                                          ref_group ==  "Countess" ~ "Countess",
                                          ref_group ==  "Biotek Synergy 2" ~ "Biotek Synergy 2",
                                          TRUE ~ ref_group))
#Rename random effect term values
df_lmm_results_all_metrics$Random_Effect_Term <- '(1 | Patient.num) + (1 + Experimental var | drug)'
#Rename column headers
colnames(df_lmm_results_all_metrics) <- gsub("_", " ", colnames(df_lmm_results_all_metrics))     
colnames(df_lmm_results_all_metrics) <- sapply(colnames(df_lmm_results_all_metrics), tools::toTitleCase)
df_lmm_results_all_metrics <- df_lmm_results_all_metrics %>% dplyr::rename(
  'Experimental Variable' = 'Fixed Effect Term',
  'Reference Condition' = 'Ref Group'
)
# #Rename fixed effect terms
# df_lmm_results_all_metrics <- df_lmm_results_all_metrics %>% mutate(`Fixed Effect Term` = case_when(`Fixed Effect Term` ==  "Positive Control: BzCl + Doses: 5 + Readout: CellTiter-Glo" ~ "Positive Control: BzCl + Doses: 5 \n+ Readout: CellTiter-Glo",
#                                                     `Fixed Effect Term` == "Microenvironmental Stimuli: Co-culture with activation" ~ "Microenvironmental Stimuli:\nCo-culture with activation", 
#                                                     `Fixed Effect Term` == "Centrifugation Procedure: LymphoPrepTM gradient centrifugation" ~ "Centrifugation Procedure: \nLymphoPrep\u2122 gradient", 
#                                                     TRUE ~ `Fixed Effect Term`))

#df_lmm_results_all_metrics <- df_lmm_results_all_metrics[order(df_lmm_results_all_metrics$`Fixed Effect Coefficient`),]
df_lmm_results_all_metrics <- df_lmm_results_all_metrics %>%
  arrange('Fixed Effect Coefficient', 'Response Metric', 'Experimental Variable')
colnames(df_lmm_results_all_metrics)
df_lmm_results_all_metrics_wide <- df_lmm_results_all_metrics[,c('Experimental Variable', 'Reference Condition', 'Response Metric', 'Fixed Effect Coefficient', 'Lower', 'Upper', 'Corrected p Value')] %>%
  pivot_wider(names_from = `Response Metric`, values_from = c(`Fixed Effect Coefficient`, Lower, Upper, `Corrected p Value`)) 
df_lmm_results_all_metrics_wide$` ` <- paste(rep(" ", 30), collapse = " ")

# Create a color theme vector to map colors to each dose response metric 
category_colors <- c("#8dd3c7", "#fccde5", "#bebada", "#fb8072", "#fdb462", "#b3de69")

#Define theme for the forest plot
tm <- forest_theme(
  base_size = 12,
  arrow_type = "closed",
  line_size = 1.5, 
  align = "center", 
  ci_fill = category_colors,
  footnote_gp = gpar(col = "black", cex = 1.0, fontsize = 18), 
  footnote_parse = FALSE, 
  legend_name = "   ",#Response Metric
  legend_position = "bottom",
  legend_value = c(expression("IC"[50]), expression(" DSS"[3]), expression(" DSS"[2]), expression(" DSS"[1]), " rAUC", " AUC"), 
  legend_gp = gpar(fontsize = 18, fontfamily = "Arial", cex = 1),
  xaxis_gp = gpar(fontsize = 18, fontfamily = "Arial"), 
  spacing = 3,
)

names(df_lmm_results_all_metrics_wide) <- str_replace_all(names(df_lmm_results_all_metrics_wide), '_', ' ')
names(df_lmm_results_all_metrics_wide) <- str_replace_all(names(df_lmm_results_all_metrics_wide), ':', '')
df_lmm_results_all_metrics_wide <- df_lmm_results_all_metrics_wide %>% dplyr::rename(
  `DSS2` = `Corrected p Value DSS2 Box Cox Transformed and Scaled`,
  `IC50` = `Corrected p Value IC50 boxcox sclaed2`
)
df_lmm_results_all_metrics_wide <- df_lmm_results_all_metrics_wide[order(df_lmm_results_all_metrics_wide$`Fixed Effect Coefficient DSS2 Box Cox Transformed and Scaled`),]

# Create the forest plot with different colors for each category bar
p1_all_metric <- forest(
  df_lmm_results_all_metrics_wide[, c("Experimental Variable","Reference Condition", " ", "DSS2", "IC50")],              # Main label column
  est = list(df_lmm_results_all_metrics_wide$`Fixed Effect Coefficient IC50 boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Fixed Effect Coefficient DSS3 boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Fixed Effect Coefficient DSS2 Box Cox Transformed and Scaled`, df_lmm_results_all_metrics_wide$`Fixed Effect Coefficient DSS1 boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Fixed Effect Coefficient auc a boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Fixed Effect Coefficient AUC boxcox sclaed2`),
  lower = list(df_lmm_results_all_metrics_wide$`Lower IC50 boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Lower DSS3 boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Lower DSS2 Box Cox Transformed and Scaled`, df_lmm_results_all_metrics_wide$`Lower DSS1 boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Lower auc a boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Lower AUC boxcox sclaed2`),   # Lower CIs for each category
  upper = list(df_lmm_results_all_metrics_wide$`Upper IC50 boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Upper DSS3 boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Upper DSS2 Box Cox Transformed and Scaled`, df_lmm_results_all_metrics_wide$`Upper DSS1 boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Upper auc a boxcox sclaed2`, df_lmm_results_all_metrics_wide$`Upper AUC boxcox sclaed2`),   # Upper CIs for each category
  ci_column = 3,                                      # CI column placement
  ref_line = 0,   
  sizes = 1,
  grid = F,
  xlim = c(-2.0, 2.0),
  #ticks_at = c(-2, -1, 0, 1, 2),
  footnote = "\nᵃBonferroni corrected",
  font.label = list(size = 18, family = "Arial"),
  font.ticks = list(size = 18, family = "Arial"),
  nudge_y = 0.15,
  theme = tm
)


output_width_cm <- 19.4 #12
output_height_cm <- 25 #10 19
dpi <- 300  # Resolution in dots per inch

# Convert cm to inches (1 inch = 2.54 cm)
output_width_in <- output_width_cm / 2.54
output_height_in <- output_height_cm / 2.54

# Convert inches to pixels for the PNG device
output_width_px <- output_width_in * dpi
output_height_px <- output_height_in * dpi

# Calculate scaling factors for the gtable
scale_width <- output_width_cm / 10  # Base width adjustment (10 cm as reference)
scale_height <- output_height_cm / 10  # Base height adjustment (10 cm as reference)

p1_all_metric <- edit_plot(p1_all_metric, gp = gpar(cex=2, fontfamily="Arial")) #1.4
p1_all_metric <- edit_plot(p1_all_metric, part = "header", gp = gpar(cex=2, fontfamily="Arial"))

p1_all_metric <- edit_plot(p1_all_metric, col = 4:5, part = "header",
                           which = "text",
                           hjust = unit(0.5, "npc"),
                           x = unit(0.5, "npc"))
p1_all_metric <- edit_plot(p1_all_metric, col = 4:5, part = "body",
                           which = "text",
                           hjust = unit(0.5, "npc"),
                           x = unit(0.5, "npc"))

p1_all_metric <- edit_plot(p1_all_metric,
                           col = 4,           # columns of header
                           part = "header",     # part of the gtable to edit
                           which = "text",      # what to edit in the grob
                           label = "DSS₂", # New text with subscript 2 (U+2082)
                           hjust = unit(0.5, "npc"),
                           x = unit(0.5, "npc"), 
                           gp = gpar(fontface = "bold"))
p1_all_metric <- edit_plot(p1_all_metric,
                           col = 5,           # columns of header
                           part = "header",     # part of the gtable to edit
                           which = "text",      # what to edit in the grob
                           label = "IC₅₀", # New text with subscript 2 (U+2082)
                           hjust = unit(0.5, "npc"),
                           x = unit(0.5, "npc"), 
                           gp = gpar(fontface = "bold"))
p1_all_metric <- add_text(p1_all_metric, text = "p-valueᵃ",
                          part = "header", 
                          row = 0,
                          col = 4:5,
                          gp = gpar(cex=2.5, fontface = "bold", fontfamily = "Arial"))

class(p1_all_metric$grobs[[5]])
p1_all_metric$grobs[[5]]$gp <- gpar(col="black",cex=2,fontfamily="Arial",fontsize=12,lineheight=1.2,alpha=1,font=2)
p1_all_metric$grobs[[1]]$gp
# p1_all_metric <- insert_text(p1_all_metric, text = expression("DSS"[2]), part = "header", row = 1, col = 4, before = TRUE, gp = gpar(fontface = "bold"))
legend_index <- 174
p1_all_metric$layout$t[21]
legend_index <- which(p1_all_metric$layout$name == "legend")
p1_all_metric$layout$t[legend_index] <- 12 #nrow(p1_all_metric)       # or your desired row index
p1_all_metric$layout$b[legend_index] <- 15 #nrow(p1_all_metric)
p1_all_metric$layout$l[legend_index] <- 4             # column 4
p1_all_metric$layout$r[legend_index] <- 4

# Scale the gtable layout
scaled_p1_all_metric <- gtable::gtable_filter(p1_all_metric, pattern = ".*", trim = TRUE)  # Keep all grobs
scaled_p1_all_metric$widths <- scaled_p1_all_metric$widths * scale_width
scaled_p1_all_metric$heights <- scaled_p1_all_metric$heights * scale_height
scaled_p1_all_metric

png(paste0(figure_output,'Forest_plot_all_metric.png'), height = 36, width = 56, unit = "cm", res = 300) #24.8 #height = 16, width = 46 ; height = 10, width = 28
grid.newpage()
grid.draw(scaled_p1_all_metric)
dev.off()
dev.off()


#Plots for article----
##Model visualization plot ----
medium_data <- subset(all_response_metrics, !is.na(medium))
#Renaming values for medium for vi
medium_data <- medium_data %>% mutate(Medium = case_when(medium == "RPMI + fetal bovine serum (FBS) (10%)" ~ "RPMI + FBS", 
                                           medium ==  "HS-5 conditioned medium" ~ "HS-5 CM",
                                           medium ==  "Mononuclear cell medium" ~ "MCM",                                                                                                              
                                           TRUE ~ medium))   
#Simplified medium model for visualising random slope and intercept for drug
model <- lmer(DSS2_boxcox_sclaed2 ~ Medium + (Medium|drug), medium_data)

plot <- flexplot::visualize(model, plot = "model", sample = 3)#, center = "Median + quartiles") 

plot$layers[[1]]$aes_params$alpha <- 0.5  #dot size
plot$layers[[2]]$aes_params$alpha <- 0.5  #dot size
plot$layers[[3]]$aes_params$linewidth <- 1 #linesize
plot$layers[[3]]$aes_params$linetype <- 2  #line
plot$layers[[3]]$aes_params$alpha <- 1
plot$layers[[4]]$aes_params$size <- 1.5
plot$layers[[5]]$aes_params$alpha <- 1 #average dot size
plot$layers[[5]]$aes_params$linewidth <- 1 #average line size
print(plot$layers)

plot <- plot + 
  scale_color_manual(values = c("#80b1d3", "#fb8072", "#bebada","#fccde5"), guide = guide_legend(title = "Drug")) +  
  labs(
    title = "Model Visualization",
    x = "Experimental Variable", 
    y = "Scaled DSS2", 
    color = "Drug"
  ) + coord_cartesian(ylim = c(-3, 3)) +   
  guides(
    linetype = "none", 
    shape = "none"
  ) + theme(
    legend.position = "top",  # Position legend on the right
    legend.title = element_text(family = "Arial", size = 10, color = "black"),  # Legend title Arial size 12
    legend.text = element_text(family = "Arial", size = 10, color = "black"),   # Legend content Arial size 12
    axis.text = element_text(family = "Arial", size = 10, color = "black"),     # Axis text Arial size 12
    axis.title = element_text(family = "Arial", size = 10, color = "black"),    # Axis titles Arial size 12
    axis.text.y = element_text(family = "Arial", size = 8, color = "black"),                    
    axis.text.x = element_text(family = "Arial", size = 8, color = "black"),    
    panel.grid = element_blank(),
    plot.title = element_text(
      size = 10,            # Font size
      face = "plain",        # Bold text
      hjust = 0.5,        
      family = "Arial", 
      vjust = 1
    )
  )

ggsave(paste0(figure_output,'Model.png'), plot = plot, width=9, height=6.5, dpi = 300, unit = "cm")

##Diagnostic plots----
all_response_metrics$medium <- as.factor(all_response_metrics$medium)
medium_model <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1|Patient.num) + (medium|drug), subset(all_response_metrics, Patient.num != '2493_BA2563D_BA2563R'), REML = TRUE)
pp_check_plot <- check_model(medium_model, check = "pp_check", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#80b1d3", "#fccde5", "black"))
plot_pp_check <- plot(pp_check_plot)
plot_pp_check_adjusted <- plot_pp_check$PP_CHECK + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                         plot.title = element_text(hjust = 0.5, vjust = 1),
                                                         axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                         axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                         axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                         axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                         plot.subtitle = element_blank(), 
                                                         legend.position = "top",  
                                                         legend.text = element_text(family = "Arial", color = "black", size = 10) 
) + labs(x = expression("BoxCox Transformed and Scaled DSS"[2]))
plot_pp_check_adjusted
ggsave(paste0(figure_output,"pp_check.png"), plot = plot_pp_check_adjusted,width = 9, height = 6.5, units="cm") #width = 9, height = 6.5

linearity_plot <- check_model(medium_model, check = "linearity", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
plot_linearity <- plot(linearity_plot)
plot_linearity_adjusted <- plot_linearity$NCV + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                      plot.title = element_text(hjust = 0.5, vjust = 1),
                                                      axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                      axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                      axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                      axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                      plot.subtitle = element_blank(), 
                                                      legend.position = "right",  
                                                      legend.text = element_text(family = "Arial", color = "black", size = 10) 
)
plot_linearity_adjusted
ggsave(paste0(figure_output,"linearity.png"), plot = plot_linearity_adjusted, width = 9, height = 6.5, units="cm")

homogeneity_plot <- check_model(medium_model, check = "homogeneity", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
plot_homogeneity <- plot(homogeneity_plot)
plot_homogeneity_adjusted <- plot_homogeneity$HOMOGENEITY + labs(y = expression(sqrt("|Std. Residuals|"))) + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                                                                                   plot.title = element_text(hjust = 0.5, vjust = 1),
                                                                                                                   axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                                                                                   axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                                                                                   axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                                                                                   axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                                                                                   plot.subtitle = element_blank(), 
                                                                                                                   show.legend = TRUE,
                                                                                                                   legend.key = element_rect(fill = "white"),
                                                                                                                   legend.position = "top",  
                                                                                                                   legend.text = element_text(family = "Arial", color = "black", size = 10) 
)
ggsave(paste0(figure_output,"homogeneity_check.png"), plot = plot_homogeneity_adjusted, width = 9, height = 6.5, units="cm")

influential_plot <- check_model(medium_model, check = "outliers", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3", "black"))
plot_influential <- plot(influential_plot)  
plot_influential_adjusted <- plot_influential$OUTLIERS + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                               plot.title = element_text(hjust = 0.5, vjust = 1),
                                                               axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                               axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                               axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                               axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                               plot.subtitle = element_blank(), 
                                                               legend.position = "right",  
                                                               legend.text = element_text(family = "Arial", color = "black", size = 10) 
)
plot_influential_adjusted
outlier_list <- check_outliers(medium_model)
outliers_info <- as.data.frame(outlier_list)
ggsave(paste0(figure_output,"influential_check.png"), plot = plot_influential_adjusted, width = 9, height = 6.5, units="cm")

normality_plot <- check_model(medium_model, check = "qq", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
plot_normality <- plot(normality_plot)
plot_normality_adjusted <- plot_normality$QQ + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                     plot.title = element_text(hjust = 0.5, vjust = 1),
                                                     axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                     axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                     axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                     axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                     plot.subtitle = element_blank(), 
                                                     legend.position = "right",  
                                                     legend.text = element_text(family = "Arial", color = "black", size = 10) 
)
plot_normality_adjusted
ggsave(paste0(figure_output,"normality_check.png"), plot = plot_normality_adjusted, width = 7.94, height = 7.49, units="cm")

Medium_model1 <- lmer(DSS2_boxcox_sclaed2 ~ Medium + (1|Patient.num) + (Medium|drug),medium_data)
random_effects_plot <- check_model(medium_model, check = "reqq", panel = F, title_size = 10, axis_title_size = 10, base_size = 8, colors = c("#fccde5", "#80b1d3"))
random_effects_plot$REQQ$drug$facet <- factor(random_effects_plot$REQQ$drug$facet, 
                                              levels = c("(Intercept)", "mediumMononuclear cell medium", "mediumRPMI + fetal bovine serum (FBS) (10%)"), 
                                              labels = c("Intercept", "MCM", "RPMI + FBS"))

plot_random_effects <- plot(random_effects_plot)
plot_random_effects[[2]]
plot_random_effects_adjusted_1 <- plot_random_effects[[1]] + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                                   plot.title = element_text(hjust = 0.5, vjust = 1),
                                                                   axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                                   axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                                   axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                                   axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                                   plot.subtitle = element_blank(), 
                                                                   legend.position = "right",  
                                                                   legend.text = element_text(family = "Arial", color = "black", size = 10), 
                                                                   panel.grid.major = element_blank(), 
                                                                   panel.grid.minor = element_blank(),  
                                                                   panel.background = element_blank()  
) + labs(title = "Normality of Random Effects - Patient") 
plot_random_effects_adjusted_1
ggsave(paste0(figure_output,"random_effects_check_1.png"), plot = plot_random_effects_adjusted_1, width = 9, height = 6.5, units="cm")

plot_random_effects_adjusted_2 <- plot_random_effects[[2]] + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                                   plot.title = element_text(hjust = 0.5, vjust = 1,color= "black"),
                                                                   axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                                   axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                                   axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                                   axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                                   plot.subtitle = element_blank(), 
                                                                   legend.position = "right",  
                                                                   legend.text = element_text(family = "Arial", color = "black", size = 10),
                                                                   panel.grid.major = element_blank(), 
                                                                   panel.grid.minor = element_blank(),  
                                                                   panel.background = element_blank(),
                                                                   strip.text = element_text(family = "Arial", color = "black", size = 8, face = "plain"),  
                                                                   strip.background = element_blank()
) + labs(title = "Normality of Random Effects - Drug")
plot_random_effects_adjusted_2
ggsave(paste0(figure_output,"random_effects_check_2.png"), plot = plot_random_effects_adjusted_2, width = 9, height = 6.5, units="cm")

comboplot1 <- (
  (plot | plot_homogeneity_adjusted)
) +
  plot_layout(heights = c(1, 1, 1, 2))

comboplot2 <- (
  (plot_pp_check_adjusted | plot_influential_adjusted) 
) +
  plot_layout(heights = c(1, 1, 1, 2))

comboplot3 <- (
  (plot_random_effects_adjusted_1 | plot_random_effects_adjusted_2) 
) +
  plot_layout(heights = c(1, 1, 1, 2))

ggsave(paste0(figure_output,"combo_plots.png"), plot = comboplot1, width = 20.76, height = 7.09, units="cm")
ggsave(paste0(figure_output,"combo_plots2.png"), plot = comboplot2, width = 20.76, height = 7.09, units="cm")
ggsave(paste0(figure_output,"combo_plots3.png"), plot = comboplot3, width = 20.76, height = 7.09, units="cm")

##Regression plots----
###Fresh vs Frozen plots----
####Karolinska----
all_response_metrics$sample
result <- subset(all_response_metrics, lab == 'Karolinska') %>%
  group_by(Patient_ID) %>%
  summarise(unique_x_count = n_distinct(sample)) %>%
  filter(unique_x_count == 2)

# Count how many unique ids meet the criteria
count_result <- nrow(result)

# Print the result
print(result)       # Shows the ids with exactly two unique instances of x
print(count_result) # Shows the count of such ids
print(unique(all_response_metrics$sample))

df1_filtered <- all_response_metrics[all_response_metrics$Patient_ID %in% result$Patient_ID, ]

new_df <- df1_filtered %>%
  group_by(Patient_ID, drug) %>%
  reframe(
    Fresh = DSS2[sample == 'fresh'],
    Frozen = DSS2[sample == 'frozen']
  )

# View the new data frame
print(new_df)

# Calculate adjusted R^2
adjusted_r2 <- summary(lm_model)$adj.r.squared

# Extract the p-value for Fresh
pval <- summary(lm_model)$coefficients[2, 4] 

# Apply FDR adjustment (assuming multiple tests, example 10)
pval_fdr <- p.adjust(pval, method = "fdr")

pearson_cor_test <- cor.test(new_df$Fresh, new_df$Frozen, method = "pearson")
pearson_cor <- pearson_cor_test$estimate
p_value <- pearson_cor_test$p.value

# Prepare the annotation text
annotation_text <- paste0(
  "R = ", round(pearson_cor, 3), 
  "\np ", ifelse(format(p_value, scientific = TRUE, digits = 3)<0.001,paste("=",format(p_value, scientific = TRUE, digits = 3)), "<0.001")
)

custom_colors <- c(
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", 
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", 
  "#bc80bd", "#ccebc5"  # Added two complementary colors
)
display.brewer.pal(n = 10, name = "Set3")
# Create the plot
scatterplot_fresh_vs_frozen_karolinska <- ggplot(new_df, aes(x = Fresh, y = Frozen)) +
  geom_point(aes(color = factor(Patient_ID)), size = 3) +
  ylim(0, 45) + 
  xlim(0,45) +
  coord_fixed(ratio = 1)+
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  scale_color_manual(values = custom_colors) +
  annotate(
    "text", 
    x = 0, y = 45, 
    label = annotation_text, 
    size = 3.51, family = "Arial", hjust = 0, vjust = 1
  ) +
  labs(color = "Patient ID", x= expression("Fresh DSS"[2]), y = expression("Frozen DSS"[2])) +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 10),    # All text Arial, size 10
    axis.text = element_text(size = 8, family = "Arial"), # Tick labels Arial, size 8
    legend.text = element_text(size = 10, family = "Arial"), # Legend text Arial, size 10
    legend.title = element_text(size = 10, family = "Arial"), # Legend title Arial, size 10
    legend.background = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(), # Optional: Remove background panel
    plot.background = element_blank()   # Optional: Remove plot background
  )

ggsave(paste0(figure_output,"Karolinska_fresh_vs_frozen.png"), plot = scatterplot_fresh_vs_frozen_karolinska, width = 10.76, height = 8.09, units="cm") #width = 20.76, height = 7.09

####Oslo----
###Diagnosis vs Relapse----
###Blood vs Bone marrow----
tissue_same_sample <- all_response_metrics %>%
  group_by(Patient_ID) %>%
  summarise(unique_x_count = n_distinct(Tissue)) %>%
  filter(unique_x_count >= 2)

tissue_filtered <- all_response_metrics[all_response_metrics$Patient_ID %in% tissue_same_sample$Patient_ID, ]

tissue_filtered <- tissue_filtered %>% mutate(Patient_ID == ifelse(
  Patient_ID != "2085" &
  Patient_ID != "2219" &
  Patient_ID != "2292", Patient_ID, Patient.num)) 
unique(subset(tissue_filtered, Patient_ID == '2085', select = c('Patient.num', 'Patient_ID')))

tissue_df <- subset(tissue_filtered, Patient_ID != "2085" & Patient_ID != "2219" & Patient_ID != "2292") %>%
  group_by(Patient_ID, drug) %>%
  reframe(
    Blood = DSS2[Tissue == 'Blood'],
    `Bone Marrow` = DSS2[Tissue == 'Bone marrow']
  )

# View the new data frame
print(tissue_df)


pearson_cor_test <- cor.test(tissue_df$Blood, tissue_df$`Bone Marrow`, method = "pearson")
pearson_cor <- pearson_cor_test$estimate
p_value <- pearson_cor_test$p.value

# Prepare the annotation text
annotation_text <- paste0(
  "R = ", round(pearson_cor, 3), 
  "\np ", ifelse(format(p_value, scientific = TRUE, digits = 3)<0.001,paste("=",format(p_value, scientific = TRUE, digits = 3)), "<0.001")
)

custom_colors <- c(
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", 
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", 
  "#bc80bd", "#ccebc5"  # Added two complementary colors
)
display.brewer.pal(n = 10, name = "Set3")
# Create the plot
scatterplot_blood_bone_marrow <- ggplot(tissue_df, aes(x = Blood, y = `Bone Marrow`)) +
  geom_point(aes(color = factor(Patient_ID)), size = 3) +
  ylim(0, 45) + 
  xlim(0,45) +
  coord_fixed(ratio = 1)+
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  #scale_color_manual(values = custom_colors) +
  annotate(
    "text", 
    x = 0, y = 45, 
    label = annotation_text, 
    size = 3.51, family = "Arial", hjust = 0, vjust = 1
  ) +
  labs(color = "Patient ID", x= expression("Blood DSS"[2]), y = expression("Bone Marrow DSS"[2])) +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 10),    # All text Arial, size 10
    axis.text = element_text(size = 8, family = "Arial"), # Tick labels Arial, size 8
    legend.text = element_text(size = 10, family = "Arial"), # Legend text Arial, size 10
    legend.title = element_text(size = 10, family = "Arial"), # Legend title Arial, size 10
    legend.background = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(), # Optional: Remove background panel
    plot.background = element_blank()   # Optional: Remove plot background
  )

ggsave(paste0(figure_output,"Bloor_Bone_Marrow.png"), plot = scatterplot_blood_bone_marrow, width = 10.76, height = 8.09, units="cm") #width = 20.76, height = 7.09



##Model assessment biospecimen----
all_response_metrics$Disease.status <- gsub('remission', 'Remission', all_response_metrics$Disease.status)
all_response_metrics$sample <- gsub('Cryopreserved', NA, all_response_metrics$sample) #Frozen
all_response_metrics$Tissue <- gsub('Leukapheresis', NA, all_response_metrics$Tissue) #Blood

frozen_count <- all_response_metrics %>% group_by(sample, drug) %>% dplyr::summarize(count = n())
tissue_count <- all_response_metrics %>% group_by(Tissue, drug) %>% dplyr::summarize(count = n())
disease_count <- all_response_metrics %>% group_by(Disease.status, drug) %>% dplyr::summarize(count = n())


mixed_model_sample_DSS2 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (0 + lab|drug), subset(all_response_metrics))
summary(mixed_model_sample_DSS2)
check_model(mixed_model_sample_DSS2)
icc(mixed_model_sample_DSS2)

mixed_model_sample_DSS2_1 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (lab|drug), subset(all_response_metrics))
summary(mixed_model_sample_DSS2_1)
check_model(mixed_model_sample_DSS2_1)
mixed_model_sample_DSS2_2 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (1|drug) + (1|lab), subset(all_response_metrics))
summary(mixed_model_sample_DSS2_2)
check_model(mixed_model_sample_DSS2_2)
mixed_model_sample_DSS2_3 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (1|drug) + (1|sample) + (1|Tissue) + (1|Disease.status), subset(all_response_metrics))
check_model(mixed_model_sample_DSS2_3)
#mixed_model_sample_DSS2_3 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (1|drug) + (1|sample) + (1|Tissue) + (1|Disease.status) + (1|lab), subset(all_response_metrics))
mixed_model_sample_DSS2_4 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|drug) + (1|sample) + (1|Tissue) + (1|Disease.status) + (1|lab), subset(all_response_metrics))
mixed_model_sample_DSS2_5 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (sample|drug) + (Tissue|drug) + (Disease.status|drug) + (1|lab), subset(all_response_metrics))
mixed_model_sample_DSS2_6 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (sample|drug) + (Tissue|drug) + (Disease.status|drug), subset(all_response_metrics))

anova(mixed_model_sample_DSS2, mixed_model_sample_DSS2_1, mixed_model_sample_DSS2_2, mixed_model_sample_DSS2_3, mixed_model_sample_DSS2_4, mixed_model_sample_DSS2_5, mixed_model_sample_DSS2_6, test = "LRT")

dss2_models_sample <-c(mixed_model_sample_DSS2_1)

##Forestplot biospecimen----
df_dss2_models_sample <- model_results(dss2_models_sample)
df_dss2_models_sample$Random_Effect_Term <- '(1|Patient.num) + (0 + Lab | Drug)'
colnames(df_dss2_models_sample) <- gsub("_", " ", colnames(df_dss2_models_sample))     
colnames(df_dss2_models_sample) <- sapply(colnames(df_dss2_models_sample), tools::toTitleCase)
rownames(df_dss2_models_sample) <- NULL

df_dss2_models_sample <- df_dss2_models_sample %>% mutate(`Biospecimen Type`= case_when(`Fixed Effect Term` == "TissueBone marrow" ~ "Tissue: Bone Marrow", 
                                                                                        `Fixed Effect Term` ==  "Disease.statusRefractory" ~ "Disease Status: Refractory",
                                                                                        `Fixed Effect Term` ==  "samplefrozen" ~ "Sample: Frozen",
                                                                                        `Fixed Effect Term` ==  "Disease.statusRelapse" ~ "Disease Status: Relapse",
                                                                                        `Fixed Effect Term` ==  "Disease.statusRemission" ~ "Disease Status: Remission",
                                                                                        TRUE ~ `Fixed Effect Term`))


df_dss2_models_sample <- df_dss2_models_sample %>% mutate(`Reference Group` = case_when(`Fixed Effect Term` == "samplefrozen" ~ "Sample: Fresh", 
                                                                                        str_detect(`Fixed Effect Term`, "Disease.status") ~ "Disease Starus: Diagnosis",
                                                                                        `Fixed Effect Term` == "TissueBone marrow" ~ "Tissue: Blood", 
                                                                                        TRUE ~ `Fixed Effect Term`))

#df_dss2_models_sample <- df_dss2_models_sample[,order(df_dss2_models_sample$`Fixed Effect Coefficient`)]

tm <- forest_theme(base_sixe = 12,
                   arrow_type = "closed",
                   footnote_gp = gpar(col = "black", cex = 0.6), 
                   align = "center", 
                   footnote_parse = FALSE, 
                   line_size = 1.5, 
                   xaxis_gp = gpar(fontsize = 10, fontfamily = "Arial"), 
                   spacing = 2,)
df_dss2_models_sample$` ` <- paste(rep(" ", 20), collapse = " ")
df_dss2_models_sample$`p-value` <- round(df_dss2_models_sample$`p Value`, 4)
df_dss2_models_sample <- df_dss2_models_sample %>% mutate(`n` = case_when(`Biospecimen Type` == 'Sample: Frozen' ~ '55-98',
                                                                          `Biospecimen Type` == 'Disease Status: Remission' ~ '2-21',
                                                                          `Biospecimen Type` == 'Disease Status: Relapse' ~ '50-97',
                                                                          `Biospecimen Type` == 'Tissue: Bone Marrow' ~ '270-600',
                                                                          `Biospecimen Type` == 'Disease Status: Refractory' ~ '10-145'
))

p_sample <- forest(df_dss2_models_sample[,c('Biospecimen Type', 'Reference Group', ' ', 'p-value', 'n')],
                   est = df_dss2_models_sample$`Fixed Effect Coefficient`,
                   lower = df_dss2_models_sample$Lower, 
                   upper = df_dss2_models_sample$Upper,
                   sizes = 1.0,
                   ci_column = 3,
                   ref_line = 0,
                   #arrow_lab = c("Lower DSS2 than reference group", "Higher DSS2 than refernce group"),
                   xlim = c(-1, 1),
                   #ticks_at = c(-1, -0.5, 0, 0.5, 1),
                   #footnote = "*1:LymphoPrepTM gradient centrifugation \n*2: Supernatant isolation at 300g 10min, density centrifugation at 400g for 20min without brake, afterwards always 300g 5 min",
                   theme = tm, 
                   xlab = expression('Change in Scaled DSS'[2]),
                   font.label = list(size = 10, family = "Arial"),
                   font.ticks = list(size = 9),)

# Print plot
plot(p_sample)
p_wh <- get_wh(p_sample)
pdf('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Forest_plot_DSS2_biospecimens.pdf',width = p_wh[1], height = 6)
plot(p_sample)
dev.off()

output_width_cm <- 12 #12 output_width_cm <- 20 #17
output_height_cm <- 10 #10 output_height_cm <- 16 #13

dpi <- 300  # Resolution in dots per inch

# Convert cm to inches (1 inch = 2.54 cm)
output_width_in <- output_width_cm / 2.54
output_height_in <- output_height_cm / 2.54

# Convert inches to pixels for the PNG device
output_width_px <- output_width_in * dpi
output_height_px <- output_height_in * dpi

# Calculate scaling factors for the gtable
scale_width <- output_width_cm / 10  # Base width adjustment (10 cm as reference)
scale_height <- output_height_cm / 10  # Base height adjustment (10 cm as reference)

p_sample <- edit_plot(p_sample, gp = gpar(cex=0.9, fontfamily="Arial")) #1.4
p_sample <- edit_plot(p_sample, part = "header", gp = gpar(cex=0.9, fontfamily="Arial"))

p_sample <- edit_plot(p_sample, col = 4:5, part = "header",
                      which = "text",
                      hjust = unit(0.5, "npc"),
                      x = unit(0.5, "npc"))
p_sample <- edit_plot(p_sample, col = 4:5, part = "body",
                      which = "text",
                      hjust = unit(0.5, "npc"),
                      x = unit(0.5, "npc"))

# Scale the gtable layout
scaled_p_sample_all_metric <- gtable::gtable_filter(p_sample, pattern = ".*", trim = TRUE)  # Keep all grobs
scaled_p_sample_all_metric$widths <- scaled_p_sample_all_metric$widths * scale_width
scaled_p_sample_all_metric$heights <- scaled_p_sample_all_metric$heights * scale_height


png(paste0(figure_output,'Forest_plot_blood.png'), height = 16, width = 26, unit = "cm", res = 300) #24.8 height = 26, width = 56 height = 10, width = 30
#10 28
grid.newpage()
grid.draw(scaled_p_sample_all_metric)
dev.off()
































########################################################################################################################
#----
#----


#Import datasets ----
#Karolinska
karolinska_fresh <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_karolinska_full_drug_set_fresh.csv')
#updating karolindka dataset
colnames(karolinska_fresh)[colnames(karolinska_fresh) == "DRUG_NAME"] <- "drug"
colnames(karolinska_fresh)[colnames(karolinska_fresh) == "DSS"] <- "DSS2"
#karolinska_fresh <- karolinska_fresh %>% mutate(Patient.num = paste(Patient.num, '_fresh'))
karolinska_fresh$sample <- 'fresh'
karolinska_frosen <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_karolinska_full_drug_set_frosen.csv')
#updating karolindka dataset
colnames(karolinska_frosen)[colnames(karolinska_frosen) == "DRUG_NAME"] <- "drug"
colnames(karolinska_frosen)[colnames(karolinska_frosen) == "DSS"] <- "DSS2"
#karolinska_frosen <- karolinska_frosen %>% mutate(Patient.num = paste(Patient.num, '_frozen'))
karolinska_frosen$sample <- 'frozen'
dss_karolinska_github <- rbind(karolinska_fresh, karolinska_frosen)
dss_karolinska_github$lab <- 'Karolinska'
#dss_karolinska_github$Tissue <- NA
#dss_karolinska_github$Disease.status <- NA
dss_karolinska_github <- dss_karolinska_github %>% mutate(Patient.num = paste0('k',dss_karolinska_github$Patient.num))
#dss_karolinska_github <- dss_karolinska_github %>% mutate(sample = str_split_i(dss_karolinska_github$Patient.num, '_', 3))
dss_karolinska_github <- dss_karolinska_github %>% mutate(Patient.num = str_split_i(dss_karolinska_github$Patient.num, '_f', 1))

karolinska <- rbind(karolinska_fresh, karolinska_frosen)

karolinska_sample_info <- read_excel('~/Desktop/UiO/Project 1/Data/Initial cleansing/info_Katarina_clinical_DSRT_data_AML_23092024-1.xlsx')
colnames(karolinska_sample_info)[colnames(karolinska_sample_info) == "Running_ID"] <- "Patient.num"
colnames(karolinska_sample_info)[colnames(karolinska_sample_info) == "AML.type"] <- "Disease.status"
colnames(karolinska_sample_info)[colnames(karolinska_sample_info) == "Sample_type"] <- "Tissue"
karolinska_sample_info$Tissue <- gsub('BM_blank', 'Bone marrow', karolinska_sample_info$Tissue)
karolinska_sample_info$Tissue <- gsub('BM', 'Bone marrow', karolinska_sample_info$Tissue)
karolinska_sample_info$Tissue <- gsub('PB', 'Blood', karolinska_sample_info$Tissue)
karolinska_sample_info$Disease.status <- gsub('de novo', 'Diagnosis', karolinska_sample_info$Disease.status)
karolinska_sample_info$Disease.status <- gsub('AML relaps', 'Relapse', karolinska_sample_info$Disease.status)
karolinska_sample_info$Disease.status <- gsub('s-AML', 'Diagnosis', karolinska_sample_info$Disease.status)
karolinska_sample_info$Disease.status <- gsub('t-AML', 'Diagnosis', karolinska_sample_info$Disease.status)
karolinska_sample_info$Disease.status <- gsub('AML', 'Diagnosis', karolinska_sample_info$Disease.status)
karolinska_sample_info$Disease.status <- gsub('APL', 'Diagnosis', karolinska_sample_info$Disease.status)
karolinska_sample_info$Disease.status <- gsub('Diagnosis, relaps', 'Relapse', karolinska_sample_info$Disease.status)
karolinska_sample_info <- karolinska_sample_info %>% mutate(Patient.num = paste0('k',karolinska_sample_info$Patient.num))

dss_karolinska_github <- merge(karolinska_sample_info, dss_karolinska_github, by = "Patient.num")
write_csv(as.data.frame(unique(karolinska$Patient.num)), '~/Desktop/UiO/Project 1/patient_ids_for_flora.csv')

karolinska_na_values <- subset(dss_karolinska_github, select = c("Patient.num", "Disease.status", "Tissue", "sample", "DSS2")) %>% group_by(Patient.num, Disease.status, Tissue, sample) %>% dplyr::summarize(avg_dss2 = mean(DSS2))



#Enserink data 
dss_github_enserink_full_set_drugs <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink_full_drug_set_testing.csv')
dss_github_enserink_full_set_drugs$lab <- 'Oslo'
#dss_github_enserink_full_set_drugs$sample <- 'fresh'
#dss_github_enserink_full_set_drugs$Tissue <- NA
#dss_github_enserink_full_set_drugs$Disease.status <- NA

enserink_sample_annotation <- read_excel('~/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_saple_annotation.xlsx')
colnames(enserink_sample_annotation)[colnames(enserink_sample_annotation) == "Patient"] <- "Patient.num"
colnames(enserink_sample_annotation)[colnames(enserink_sample_annotation) == "Fresh or frozen"] <- "sample"
colnames(enserink_sample_annotation)[colnames(enserink_sample_annotation) == "Cells/well"] <- "Cells1"
colnames(enserink_sample_annotation)[colnames(enserink_sample_annotation) == "Instrument"] <- "Plate_reader1"
enserink_sample_annotation$Tissue <- gsub('BM', 'Bone marrow', enserink_sample_annotation$Tissue)
enserink_sample_annotation$Tissue <- gsub('PB', 'Blood', enserink_sample_annotation$Tissue)
enserink_sample_annotation <- enserink_sample_annotation %>% mutate(Disease.status = str_split_i(enserink_sample_annotation$Patient.num, '_', 2))
enserink_sample_annotation[is.na(enserink_sample_annotation$Disease.status),]$Disease.status <- 'Diagnosis'
enserink_sample_annotation$Disease.status <- gsub('relapse', 'Relapse', enserink_sample_annotation$Disease.status)
enserink_sample_annotation$Disease.status <- gsub('remission2', 'remission', enserink_sample_annotation$Disease.status)
enserink_sample_annotation$sample <- gsub('Fresh', 'fresh', enserink_sample_annotation$sample)
enserink_sample_annotation <- enserink_sample_annotation %>% mutate(Cells1 = case_when(Cells1 > 5000 ~"10000", 
                                                                                       Cells1 == 5000 ~ "5000", 
                                                                                       Cells1 < 5000 ~ NA,
                                                                                       TRUE ~ NA))
unique(enserink_sample_annotation$Cells1)
enserink_sample_annotation %>% distinct(Patient.num, Cells1) %>% group_by(Cells1) %>% dplyr:: summarize(c = n())
count(enserink_sample_annotation[enserink_sample_annotation$Disease.status == 'Relapse',]$Disease.status)



dss_github_enserink_full_set_drugs <- merge(enserink_sample_annotation, dss_github_enserink_full_set_drugs, by = "Patient.num")
group_counts_enserink <- dss_github_enserink_full_set_drugs %>%
  group_by(lab, drug) %>%
  dplyr::summarize(count = n())


#BeatAML ----
dss_github_beat_aml <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_beat_aml_full_drug_set_v1.csv')
dss_github_beat_aml$lab <- 'Beat AML'
dss_github_beat_aml$sample <- 'fresh'
beat_aml_clinical_info <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_clinical.csv')
beat_aml_clinical_info <- beat_aml_clinical_info %>% mutate(Patient.num = paste0(beat_aml_clinical_info$dbgap_subject_id, '_', dbgap_dnaseq_sample, '_', dbgap_rnaseq_sample))
colnames(beat_aml_clinical_info)[colnames(beat_aml_clinical_info) == "specimenType"] <- "Tissue"
colnames(beat_aml_clinical_info)[colnames(beat_aml_clinical_info) == "diseaseStageAtSpecimenCollection"] <- "Disease.status"
beat_aml_clinical_info$Tissue <- gsub('Bone Marrow Aspirate', 'Bone marrow', beat_aml_clinical_info$Tissue)
beat_aml_clinical_info$Tissue <- gsub('Peripheral Blood', 'Blood', beat_aml_clinical_info$Tissue)
beat_aml_clinical_info$Disease.status <- gsub('Initial ', '', beat_aml_clinical_info$Disease.status)
beat_aml_clinical_info$Disease.status <- gsub('Residual', 'Refractory', beat_aml_clinical_info$Disease.status)
dss_github_beat_aml <- merge(subset(beat_aml_clinical_info, select = c("Patient.num", "Disease.status", "Tissue")), dss_github_beat_aml, by = "Patient.num", all.y = TRUE)

unique(beat_aml_clinical_info$Disease.status)

beat_aml_clinical_info <- beat_aml_clinical_info %>% mutate(Patient_ID = sub("_.*", "", Patient.num))

beat_aml_clinical_info$consensus_sex
beat_aml_clinical_info %>% group_by(consensus_sex) %>% dplyr::summarize(n = n_distinct(Patient_ID))

#FIMM ----
#dss_github_fimm_full_drug_set_supplemental.csv
dss_github_fimm <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set.csv') #v1 is where concentration range is from1 to 10000, not v1 is concentration range 0.1 ro 1000
dss_github_fimm$lab <- 'Helsinki'


fimm_sample_information <- read.csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_assay_dets.csv')
colnames(fimm_sample_information)[colnames(fimm_sample_information) == "Medium"] <- "Medium1"
fimm_sample_information$Medium1 <- gsub("MCM", 'Mononuclear cell medium', fimm_sample_information$Medium1)
fimm_sample_information$Medium1 <- gsub("CM", 'HS-5 conditioned medium', fimm_sample_information$Medium1)
colnames(fimm_sample_information)[colnames(fimm_sample_information) == "Sample_ID"] <- "Patient.num"
unique(fimm_sample_information$Medium1)
fimm_sample_annotation <- read.csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_sample_annotation.csv')
colnames(fimm_sample_annotation)[colnames(fimm_sample_annotation) == "Sample_ID"] <- "Patient.num"
fimm_clinical_info <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_clinical_summary.csv')
fimm_clinical_info$sample <- ifelse(fimm_clinical_info$Site == "Helsinki", 'fresh', 'frozen')

#Getting the common drugs across the 4 datasets
common_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/common_drugnames_pubchem.csv")
common_drugs$`Unnamed: 0_x` <- NULL
common_drugs$`Unnamed: 0_y` <- NULL
common_drugs$`Unnamed: 0` <- NULL

#Karolinska common drugs 
dss_karolinska_common_drugs <- inner_join(common_drugs, dss_karolinska_github, by=c("karolinska_drug_name"="drug"))
colnames(dss_karolinska_common_drugs)[colnames(dss_karolinska_common_drugs) == "pubchem_drug_name"] <- "drug"
dim(dss_karolinska_common_drugs)


#Enserink common drugs 
dss_enserink_common_drugs <- inner_join(common_drugs, dss_github_enserink_full_set_drugs, by=c("enserink_drug_name"="drug"))
colnames(dss_enserink_common_drugs)[colnames(dss_enserink_common_drugs) == "pubchem_drug_name"] <- "drug"
dim(dss_enserink_common_drugs)

group_counts_enserink <- dss_enserink_common_drugs %>%
  group_by(lab, drug) %>%
  dplyr::summarize(count = n())

#BeatAML common drugs 
dss_beat_aml_common_drugs <- inner_join(common_drugs, dss_github_beat_aml, by=c("beat_aml_drug_name"="drug"))
colnames(dss_beat_aml_common_drugs)[colnames(dss_beat_aml_common_drugs) == "pubchem_drug_name"] <- "drug"
dim(dss_beat_aml_common_drugs)

length(unique(dss_github_fimm$Patient.num))

#FIMM - My calculations from the github package
dss_fimm_common_drugs <- inner_join(common_drugs, dss_github_fimm, by=c("fimm_drug_name"="drug"))
colnames(dss_fimm_common_drugs)[colnames(dss_fimm_common_drugs) == "pubchem_drug_name"] <- "drug"
keep_letters_numbers <- function(x) {
  sub("^([^_]+_[^_]+)_.+$", "\\1", x)
}
#dss_fimm_common_drugs <- dss_fimm_common_drugs %>% mutate(Patient.num = keep_letters_numbers(Patient.num))
dim(dss_fimm_common_drugs)
dss_fimm_common_drugs$Patient_ID <- gsub("(_[^_]*){1}$", "", dss_fimm_common_drugs$Patient.num)
dss_fimm_common_drugs <- merge(subset(fimm_clinical_info, select = c(Patient_ID, sample)), dss_fimm_common_drugs, by="Patient_ID", all.y = TRUE)
dss_fimm_common_drugs$Patient_ID <- NULL
dss_fimm_common_drugs <- merge(subset(fimm_sample_annotation, select = c(Patient.num, Tissue, Disease.status)), dss_fimm_common_drugs, by="Patient.num", all.y = TRUE)
dss_fimm_common_drugs$Disease.status <- gsub('Diagnosis ', 'Diagnosis', dss_fimm_common_drugs$Disease.status)
dss_fimm_common_drugs$Disease.status <- gsub('Refractory ', 'Refractory', dss_fimm_common_drugs$Disease.status)

test <- read.csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/cm125_ctg_fo4b_fo5a_FIMM_125_AML_patients_oslo_14052024.csv')
test2 <- read_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_4_DSRT_assay_details_164S_DM.xlsx')

#fimm_org_response_long
#dss_github_fimm
#dss_fimm_common_drugs




# Merge the Data Frames on col1 and col2
# merged_df <- merge(fimm_github_dss, fimm_raw_response, by = c('Patient.num', 'drug'), suffixes = c('_df1', '_df2'))
# 
# # Check if col3 values are equal
# merged_df$col3_equal <- merged_df$DSS2_df1 == merged_df$DSS2_df2
# 
# # Show the result
# print(merged_df["col3_equal" == FALSE],)
# subset(merged_df, col3_equal == FALSE)

fimm_raw_response <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_raw_dose_responses.csv.csv')
fimm_raw_response <- fimm_raw_response[, !(names(fimm_raw_response) %in% c('Solvent', 'Sediment', 'Final.Conc', 'Dilution', 'Dest.Plate', 'Dest.Well', 'DRow', 'DCol', 'doseResponses', 'X', 'FO5No', 'Plate', 'NoDiv'))] %>% distinct() %>% as.data.frame()
fimm_raw_response <- unique(subset(fimm_raw_response, select = c(sample, drug, dss)))
colnames(fimm_raw_response)[colnames(fimm_raw_response) == "sample"] <- "Patient.num"
colnames(fimm_raw_response)[colnames(fimm_raw_response) == "dss"] <- "DSS2"
fimm_raw_response_common_drugs <- inner_join(common_drugs, fimm_raw_response, by=c("fimm_drug_name"="drug"))
colnames(fimm_raw_response_common_drugs)[colnames(fimm_raw_response_common_drugs) == "pubchem_drug_name"] <- "drug"
#fimm_raw_response_common_drugs <- fimm_raw_response_common_drugs %>% mutate(Patient.num = paste(Patient.num, '_FIMM'))
dim(fimm_raw_response_common_drugs)
fimm_raw_response_common_drugs$lab <- 'Helsinki'



fimm_org_response <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_drug_response.csv')
fimm_org_response$Drug_ID <- NULL
fimm_org_response$...1 <- NULL
fimm_org_response <- as.data.frame(fimm_org_response)
fimm_org_response_long <- gather(fimm_org_response, Patient.num, DSS2, 'AML_084_04':'Healthy_17')
length(unique(fimm_org_response_long$Patient.num))
colnames(fimm_org_response_long)[colnames(fimm_org_response_long) == "Drug_name"] <- "drug"
fimm_org_response_common_drugs <- inner_join(common_drugs, fimm_org_response_long, by=c("fimm_drug_name"="drug"))
colnames(fimm_org_response_common_drugs)[colnames(fimm_org_response_common_drugs) == "pubchem_drug_name"] <- "drug"
#fimm_org_response_common_drugs <- fimm_org_response_common_drugs %>% mutate(Patient.num = paste(Patient.num, '_FIMM'))
dim(fimm_org_response_common_drugs)
fimm_org_response_common_drugs$lab <- 'Helsinki'
fimm_org_response_common_drugs$Patient_ID <- gsub("(_[^_]*){1}$", "", fimm_org_response_common_drugs$Patient.num)
fimm_org_response_common_drugs <- merge(subset(fimm_clinical_info, select = c(Patient_ID, sample)), fimm_org_response_common_drugs, by="Patient_ID", all.y = TRUE)
fimm_org_response_common_drugs$Patient_ID <- NULL
fimm_org_response_common_drugs <- merge(subset(fimm_sample_annotation, select = c(Patient.num, Tissue, Disease.status)), fimm_org_response_common_drugs, by="Patient.num", all.y = TRUE)
fimm_org_response_common_drugs$Disease.status <- gsub('Diagnosis ', 'Diagnosis', fimm_org_response_common_drugs$Disease.status)
fimm_org_response_common_drugs$Disease.status <- gsub('Refractory ', 'Refractory', fimm_org_response_common_drugs$Disease.status)


unique(fimm_sample_annotation$Disease.status)


subset(dss_karolinska_common_drugs, DSS2 < 0)
subset(dss_enserink_common_drugs, DSS2 < 0)
subset(dss_beat_aml_common_drugs, DSS2 < 0)



#Missingnesss
# all_wide <- all_datasets %>% subset(select = c(drug, Patient.num, DSS2)) %>% pivot_wider(names_from = Patient.num, values_from = DSS2) 
# all_wide <- as.data.frame(all_wide)
# rownames(all_wide) <- all_wide$drug
# all_wide$drug <- NULL
# all_wide <- as.data.frame(t(all_wide))
# vis_dat(all_wide)
# vis_miss(all_wide)
# 
# 
# beat_aml_wide_all <- all_datasets %>% subset(lab == 'BeatAML',select = c(drug, Patient.num, DSS2)) %>% pivot_wider(names_from = Patient.num, values_from = DSS2) 
# beat_aml_wide_all <- as.data.frame(beat_aml_wide_all)
# rownames(beat_aml_wide_all) <- beat_aml_wide_all$drug
# beat_aml_wide_all$drug <- NULL
# beat_aml_wide_all <- t(beat_aml_wide_all) %>% as.data.frame()
# # Check how many values are missing in each row
# missing_per_col <- apply(beat_aml_wide_all, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in row is higher than 200
# df_filtered <- beat_aml_wide_all[,missing_per_col <= 200]
# # Check how many values are missing in each col
# missing_per_row <- apply(df_filtered, 1, function(row) sum(is.na(row)))
# #Remove col if total missing in row is higher than 100
# df_filtered <- df_filtered[missing_per_row <= 5,]
# vis_dat(df_filtered)
# vis_miss(df_filtered)
# 
# df_filtered$Patient.num <- rownames(df_filtered)
# beat_aml_dss_filtered_common_drugs <- gather(df_filtered, drug, DSS2, '781661-94-7':'AT7519', factor_key=TRUE)
# dim(beat_aml_dss_filtered_common_drugs)
# beat_aml_dss_filtered_common_drugs <- inner_join(beat_aml_dss_filtered_common_drugs, dss_beat_aml_common_drugs, by = c('drug', 'Patient.num', 'DSS2'))
# dim(beat_aml_dss_filtered_common_drugs)

#Read in experimental variabilities 
experimental_var <- read_delim('~/Desktop/UiO/Project 1/Data/Initial cleansing/experimental_variables.csv', delim=';')
experimental_var[experimental_var$lab == 'FIMM',]$lab <- 'Helsinki'
experimental_var[experimental_var$lab == 'Enserink',]$lab <- 'Oslo'
experimental_var[experimental_var$lab == 'Helsinki','medium'] <- NA
unique(experimental_var$time_until_sample_usage)

#experimental_var <- as.data.frame(experimental_var)
#experimental_var[experimental_var$lab == 'FIMM','nr_positive_control'] <- (8+10)/2
#experimental_var[is.na(experimental_var$nr_positive_control),'nr_positive_control'] <- (12+9+8)/3
#experimental_var[experimental_var$lab == 'FIMM','nr_negative_control'] <- (10+12)/2
#experimental_var[is.na(experimental_var$nr_negative_control),'nr_negative_control'] <- (12+11+8)/3

#cols_with_na <- colSums(is.na(experimental_var)) > 0
#experimental_var <- experimental_var[, !cols_with_na]

#Combine all datasets underneath eachother/append/rbind
# all_datasets <- rbind(subset(dss_karolinska_common_drugs,select=-c(conc, AUC, DSS1, DSS3, IC50, believe_DSS, ID)), subset(dss_enserink_common_drugs, select=-c(conc, AUC, DSS1, DSS3, IC50, believe_DSS, ID)), subset(dss_beat_aml_common_drugs, select=-c(conc, AUC, DSS1, DSS3, IC50, believe_DSS, ID)), fimm_org_response_common_drugs) #fimm_raw_response -> fimm_org_response_common_drugs
# all_datasets <- all_datasets %>% mutate(drug_dataset_name = paste(lab, ' ', drug))
# all_datasets <- inner_join(all_datasets, subset(experimental_var, select = -c(sample)), by = 'lab')
# all_datasets <- merge(subset(fimm_sample_information, select = c("Patient.num", "Medium")), all_datasets, by = "Patient.num",  all.y = TRUE)
# all_datasets <- all_datasets %>%
#   mutate(medium = if_else(lab == 'FIMM', Medium, medium))


#With git calculated DSS2
ncol(subset(dss_karolinska_common_drugs,select=-c(conc, AUC, DSS1, DSS3, IC50, cohort, DSRT)))
ncol(subset(dss_enserink_common_drugs, select=-c(conc, AUC, DSS1, DSS3, IC50, Cells1, Plate_reader1, Column2, Column1)))
ncol(subset(dss_beat_aml_common_drugs, select=-c(conc, AUC, DSS1, DSS3, IC50)))
ncol(subset(dss_fimm_common_drugs, select=-c(conc, AUC, DSS1, DSS3, IC50)))
all_datasets_v2 <- rbind(subset(dss_karolinska_common_drugs,select=-c(conc, AUC, DSS1, DSS3, IC50, cohort, DSRT)),
                         subset(dss_enserink_common_drugs, select=-c(conc, AUC, DSS1, DSS3, IC50, Cells1, Plate_reader1, Column2, Column1)), 
                         subset(dss_beat_aml_common_drugs, select=-c(conc, AUC, DSS1, DSS3, IC50)), 
                         subset(dss_fimm_common_drugs, select=-c(conc, AUC, DSS1, DSS3, IC50))) #fimm_raw_response -> fimm_org_response_common_drugs

all_datasets_v2 <- inner_join(all_datasets_v2, subset(experimental_var, select = -c(sample)), by = 'lab')
all_datasets_v2 <- subset(all_datasets_v2, drug != 'Sapanisertib')
all_datasets_v2 <- all_datasets_v2 %>% mutate(nr_of_concentration_points_cat = as.character(nr_of_concentration_points)) %>% as.data.frame()
all_datasets_v2 <- all_datasets_v2 %>% mutate(cells = as.character(cells)) %>% as.data.frame()
fimm_sample_information[fimm_sample_information$Patient.num == '',]
all_datasets_v2 <- merge(subset(fimm_sample_information, Patient.num != '', select = c("Patient.num", "Medium1", "Blast_cell_percentage")), all_datasets_v2, by = "Patient.num", all.y=TRUE)
fimm_sample_information$Medium1
subset(all_datasets_v2, lab == 'Helsinki')
all_datasets_v2 <- merge(unique(subset(dss_enserink_common_drugs,Patient.num != '', select = c("Patient.num", "Cells1", "Plate_reader1"))), all_datasets_v2, by = "Patient.num", all.y=TRUE)
subset(all_datasets_v2, lab == 'Oslo')
group_counts_lab <- all_datasets_v2 %>%
  group_by(lab, drug) %>%
  dplyr::summarize(count = n())

all_datasets_v2 <- all_datasets_v2 %>%
  mutate(medium = if_else(lab == 'Helsinki', Medium1, medium))
unique(all_datasets_v2[all_datasets_v2$lab == 'Helsinki', "medium"])
all_datasets_v2$Cells1 <- as.character(all_datasets_v2$Cells1)
all_datasets_v2 <- all_datasets_v2 %>%
  mutate(cells = if_else(lab == 'Oslo', Cells1, cells))
all_datasets_v2 <- all_datasets_v2 %>%
  mutate(plate_reader = if_else(lab == 'Oslo', Plate_reader1, plate_reader))
all_datasets_v2 <- subset(all_datasets_v2, !is.na(Patient.num))
all_datasets_v2$lab <- gsub('BeatAML', 'Beat AML', all_datasets_v2$lab)


subset(all_datasets_v2, lab == 'Helsinki')
 write_csv(all_datasets_v2, '~/Desktop/George/Dataset_for_George.csv')
 group_counts_lab <- all_datasets_v2 %>%
   group_by(lab, drug) %>%
   dplyr::summarize(count = n())
 
all_datasets_v2$Disease.status <- gsub('remission', 'Remission', all_datasets_v2$Disease.status)
all_datasets_v2$Disease.status <- gsub('Unknown', NA, all_datasets_v2$Disease.status)
all_datasets_v2$sample <- gsub('Cryopreserved', 'frozen', all_datasets_v2$sample)
all_datasets_v2$Tissue <- gsub('Leukapheresis', 'Blood', all_datasets_v2$Tissue)
unique(all_datasets_v2$cell_counting_method)

all_datasets_v2 %>% group_by(lab, cells) %>% dplyr::summarize(c = n())

all_response_metrics <- rbind(subset(dss_karolinska_common_drugs,select=-c(cohort, DSRT)),
                              subset(dss_enserink_common_drugs, select=-c(Cells1, Plate_reader1, Column2, Column1)), 
                              dss_beat_aml_common_drugs, 
                              dss_fimm_common_drugs) #fimm_raw_response -> fimm_org_response_common_drugs
all_response_metrics <- inner_join(all_response_metrics, subset(experimental_var, select = -c(sample)), by = 'lab')
all_response_metrics <- subset(all_response_metrics, drug != 'Sapanisertib')
all_response_metrics_df_list <- list(Karolinska = dss_karolinska_common_drugs, Enserink = dss_enserink_common_drugs, BeatAML = dss_beat_aml_common_drugs, FIMM = dss_fimm_common_drugs)
all_response_metrics <- all_response_metrics %>% mutate(nr_of_concentration_points_cat = as.character(nr_of_concentration_points)) %>% as.data.frame()
all_response_metrics <- all_response_metrics %>% mutate(cells = as.character(cells)) %>% as.data.frame()

all_lab_auc <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/auc_calculation_all_labs.csv')
all_response_metrics <- inner_join(all_response_metrics, all_lab_auc, by=c("Patient.num", "drug"))
all_response_metrics <- merge(unique(subset(fimm_sample_information, Patient.num != '', select = c("Patient.num", "Medium1"))), all_response_metrics, by = "Patient.num", all.y=TRUE)
all_response_metrics <- all_response_metrics %>%
  mutate(medium = if_else(lab == 'Helsinki', Medium1, medium))
all_response_metrics <- right_join(unique(subset(dss_enserink_common_drugs, select = c("Patient.num", "Cells1", "Plate_reader1"))), all_response_metrics, by = "Patient.num")
all_response_metrics$Cells1 <- as.character(all_response_metrics$Cells1)
all_response_metrics <- all_response_metrics %>%
  mutate(cells = if_else(lab == 'Oslo', Cells1, cells))
all_response_metrics <- all_response_metrics %>%
  mutate(plate_reader = if_else(lab == 'Oslo', Plate_reader1, plate_reader))
all_response_metrics <- subset(all_response_metrics, !is.na(Patient.num))
subset(all_response_metrics, lab =='Helsinki')
write_csv(all_response_metrics, '~/Desktop/UiO/Project 1/Data/Initial cleansing/all_metrices_all_labs.csv')
group_counts_all_metric <- all_response_metrics %>%
  group_by(lab, drug) %>%
  dplyr::summarize(count = n())

all_response_metrics$Disease.status <- gsub('remission', 'Remission', all_response_metrics$Disease.status)
all_response_metrics$Disease.status <- gsub('Unknown', NA, all_response_metrics$Disease.status)
all_response_metrics$sample <- gsub('Cryopreserved', 'frozen', all_response_metrics$sample)
all_response_metrics$Tissue <- gsub('Leukapheresis', 'Blood', all_response_metrics$Tissue)
all_response_metrics$lab <- gsub('BeatAML', 'Beat AML', all_response_metrics$lab)


list_common_dfs <- list(Karolinska = dss_karolinska_common_drugs, Enserink = dss_enserink_common_drugs, BeatAML = dss_beat_aml_common_drugs, FIMM = dss_fimm_common_drugs)

drug_names_to_include <- unique(all_response_metrics$drug)#c('cytarabine', '345627-80-7', '781661-94-7', 'Bortezomib', 'PH-797804', 'Quizartinib', 'PONATINIB', 'Ruxolitinib', 'Sorafenib', 'Cyt387')

write_csv(dss_enserink_common_drugs, '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Second cleansing/dss_enserink_common_drugs.csv')


all_response_metrics$Patient_ID <- all_response_metrics$Patient.num
all_response_metrics <- all_response_metrics %>% dplyr::mutate(Patient_ID = ifelse(lab == 'Oslo' | lab == "Beat AML", gsub("_.*", "", Patient_ID), Patient_ID)) %>%
  mutate(Patient_ID = ifelse(lab == 'Helsinki', str_extract(Patient_ID, "^[^_]+_[^_]+"), Patient_ID))

#FIMM
fimm_explore <- subset(all_datasets_v2, lab == 'Helsinki')
length(unique(fimm_explore$Patient.num))
fimm_explore %>%
  group_by(lab) %>%
  summarise()


#BeatAML
beatAML_explore <- subset(all_datasets_v2, lab == 'BeatAML')
length(unique(beatAML_explore$Patient.num))


#reference https://github.com/G-Thomson/gthor
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
  
}

#Count outliers ----

n_fun <- function(x){
  return(data.frame(y = 0.95*70,
                    label = length(x)))
}

# Function to identify and count extreme outliers
count_extreme_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  lower_bound <- q1 - 3 * iqr
  upper_bound <- q3 + 3 * iqr
  sum(x < lower_bound | x > upper_bound)
}

create_density_plots <- function(df, y, x, drug_name_list, title){
  d1 <- df %>%
    filter(get(y) %in% drug_name_list)%>%
    ggplot( aes(y=get(y), x=get(x),  fill=drug)) +
    geom_density_ridges(alpha=0.6, stat="binline", bins=20) +
    theme_ridges() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    xlab("Density") +
    ylab("Drug names") + 
    labs(title=title)
  print(d1)
  
  
  d2 <- df %>%
    filter(get(y) %in% drug_name_list)%>%
    ggplot(aes(y=get(y), x=get(x),  fill=drug)) +
    geom_density_ridges_gradient(scale = 1.7, rel_min_height = 0.01) +
    theme_ridges() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    xlab("Density") +
    ylab("Drug names") + 
    labs(title=title) + 
    xlim(0, 50)
  print(d2)
  
  #Counting number of values that are considered to be outliers
  outlier_counts <- aggregate(get(x) ~ get(y), df %>%
                                filter(get(y) %in% drug_name_list), count_extreme_outliers)
  print(outlier_counts)
  
  d3 <- df %>%
    filter(get(y) %in% drug_name_list)%>%
    ggplot(aes(y=get(y), x=get(x),  fill=drug)) +
    geom_boxplot() +
    theme_minimal()+
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8), 
      axis.title.y = element_blank()
    ) +
    stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5) +
    xlab("") +
    ylab("Drug names") + 
    labs(title=title) 
  
  print(d3)
  
  d4ensity_plots_list <- list(d1=d1, d2=d2, d3=d3)
  return(d4ensity_plots_list)
}



#Density plotting ----
# Loop through the list of data frames and create a plot for each
d_plots <- lapply(names(list_common_dfs), function(df_name) {
  create_density_plots(list_common_dfs[[df_name]], 'drug', 'DSS2', drug_names_to_include, title = paste("Density plot of selected drugs for the ", df_name, " dataset") )
})


print(d_plots)

d1_grid_plot <- plot_grid(d_plots[[1]]$d1,d_plots[[2]]$d1,d_plots[[3]]$d1,d_plots[[4]]$d1, rel_widths = c(.8,.6))
ggsave(file="~/Desktop/UiO/Project 1/Figures/density_grid_plots_boxes_all_datasets.pdf", d1_grid_plot, units="cm", width = 80, height = 50)

d2_grid_plot <- plot_grid(d_plots[[1]]$d2,d_plots[[2]]$d2,d_plots[[3]]$d2,d_plots[[4]]$d2, rel_widths = c(.8,.6))
ggsave(file="~/Desktop/UiO/Project 1/Figures/density_grid_plots_all_datasets.pdf", d2_grid_plot, units="cm", width = 80, height = 50)

d3_grid_plot <- plot_grid(d_plots[[1]]$d3,d_plots[[2]]$d3,d_plots[[3]]$d3,d_plots[[4]]$d3, nrow=1,rel_widths = c(.8,.6))
ggsave(file="~/Desktop/UiO/Project 1/Figures/box_density_grid_plots_all_datasets.pdf", d3_grid_plot, units="cm", width = 50, height = 30)



# Function to create density plots
create_density_plots <- function(df_list, y, x, drug_name_list, title) {
  # Combine data frames into one with an additional column for the dataset name
  combined_df <- bind_rows(lapply(names(df_list), function(df_name) {
    df_list[[df_name]] %>% mutate(Source = df_name)
  }))
  
  # Filter the combined data frame based on the drug names
  filtered_df <- combined_df %>%
    filter(get(y) %in% drug_name_list)
  
  # Use a pastel color palette
  pastel_colors <- c("#FFB3BA", "#FFFFBA", "#BAFFC9", "#BAE1FF")

  # Plot overlapping density plots with transparent pastel colors
  d1 <- ggplot(filtered_df, aes(x = get(x), y = get(y), color=Source)) +
    geom_density_ridges(alpha = 0, stat = "binline", bins = 20) + 
    scale_color_manual(values = pastel_colors, name = 'Lab') +
    theme_ridges() +
    theme(
      panel.spacing = unit(0.5, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    xlab("Density") +
    ylab("Drug names") + 
    labs(title = title)
  print(d1)
  
  d2 <- ggplot(filtered_df, aes(x = get(x), y = get(y))) +
    geom_density_ridges_gradient(alpha = 0, scale = 1.5, rel_min_height = 0.01, aes(fill=Source)) + 
    scale_fill_manual(values = pastel_colors, name = 'Lab') +
    theme_ridges() +
    theme(
      panel.spacing = unit(5.5, "lines"),
      strip.text.x = element_text(size = 8), 
      panel.background = element_rect(fill = "transparent"), 
      plot.background = element_rect(fill = "transparent", color = NA), 
      panel.grid = element_blank()
    ) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.2))) +
    xlab("Density") +
    ylab("Drug names") + 
    labs(title = title) + 
    xlim(0, 50)
  print(d2)
  
  # Counting number of values that are considered to be outliers
  outlier_counts <- filtered_df %>%
    group_by(get(y), Source) %>%
    summarise(outlier_count = sum(get(x) > 50)) # Adjust threshold as needed
  print(outlier_counts)
  
  # Custom function to add the number of observations on the opposite side of the axis
  add_n_obs <- function(data, x, y, yvar, xpos = 50) {
    data %>%
      group_by(!!sym(y)) %>%
      summarise(count = n()) %>%
      mutate(xpos = xpos) %>%
      ggplot(aes(x = xpos, y = !!sym(y), label = count)) +
      geom_text(hjust = 1) +
      theme_void()
  }

  # Create the boxplot with observation counts on the opposite side of the axis
  d3 <- ggplot(filtered_df, aes(x = get(x), y = get(y))) +
    geom_boxplot(alpha = 0.6, scale = 3, aes(color = Source, fill=Source)) + # Transparent pastel colors
    scale_fill_manual(values = pastel_colors, name = 'Lab') +
    scale_color_manual(values = pastel_colors, name = 'Lab') +       
    theme_minimal() +
    theme(
      panel.spacing = unit(1, "lines"),
      strip.text.x = element_text(size = 8), 
      axis.title.y = element_blank()
    ) +
    scale_y_discrete(expand = expansion(mult = c(0, 0))) +  # Adjust spacing between y levels
    xlab("DSS2") +
    ylab("Drug names") + 
    labs(title = title) 
  print(d3)
  
  density_plots_list <- list(d1 = d1, d2 = d2, d3 = d3)
  return(density_plots_list)
}

# Create density plots for the combined datasets
d_plots <- create_density_plots(list_common_dfs, 'drug', 'DSS2', drug_names_to_include, title = "Density plot of common drugs across different datasets")

ggsave(file="~/Desktop/UiO/Project 1/Figures/V2_prepare_article/density_grid_plots_boxes_all_datasets_overlapping_drugs.pdf", d_plots$d1, units="cm", width = 80, height = 50)

ggsave(file="~/Desktop/UiO/Project 1/Figures/V2_prepare_article/density_grid_plots_all_datasets_overlapping_drugs.pdf", d_plots$d2, units="cm", width = 80, height = 50)

ggsave(file="~/Desktop/UiO/Project 1/Figures/V2_prepare_article/box_density_grid_plots_all_datasets_overlapping_drugs.pdf", d_plots$d3, units="cm", width = 50, height = 30)

d_plots_alldrugs <- create_density_plots(list_common_dfs, 'drug', 'DSS2', unique(all_datasets$drug), title = "Density plot of selected drugs across different datasets")


drug_set1 <- c('Bortezomib','Navitoclax', 'Venetoclax', 'Bortezomib')
drug_set2 <- c('Bortezomib','Venetoclax')

e_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_drug_information_pubchem.csv")
dss_enserink_all_drugs_pubchem <- inner_join(e_drugs, dss_github_enserink_full_set_drugs, by=c("org_drug_name"="drug"))
colnames(dss_enserink_all_drugs_pubchem)[colnames(dss_enserink_all_drugs_pubchem) == "pubchem_drug_name"] <- "drug"
k_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/karolinska_drug_library_pubchem.csv")
dss_karolinska_all_drugs_pubchem <- inner_join(k_drugs, dss_karolinska_github, by=c("org_drug_name"="drug"))
colnames(dss_karolinska_all_drugs_pubchem)[colnames(dss_karolinska_all_drugs_pubchem) == "pubchem_drug_name"] <- "drug"
f_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_drug_library_pubchem.csv")
dss_fimm_all_drugs_pubchem <- inner_join(f_drugs, dss_github_fimm, by=c("org_drug_name"="drug"))
colnames(dss_fimm_all_drugs_pubchem)[colnames(dss_fimm_all_drugs_pubchem) == "pubchem_drug_name"] <- "drug"
ba_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/beataml_drug_information_pubchem.csv")
dss_ba_all_drugs_pubchem <- inner_join(ba_drugs, dss_github_beat_aml, by=c("org_drug_name"="drug"))
colnames(dss_ba_all_drugs_pubchem)[colnames(dss_ba_all_drugs_pubchem) == "pubchem_drug_name"] <- "drug"

list_karolinska_enserink_fimm_ba <- list(Karolinska = dss_karolinska_all_drugs_pubchem, Enserink = dss_enserink_all_drugs_pubchem, FIMM = dss_fimm_all_drugs_pubchem, BeatAML = dss_ba_all_drugs_pubchem)
subset(dss_fimm_all_drugs_pubchem, drug == 'Navitoclax')
d_plots_w_navitoclax <- create_density_plots(list_karolinska_enserink_fimm_ba, 'drug', 'DSS2', drug_set1, title = "Density plot of selected drugs across different datasets")


#Density heatmap
cbind.fill <- function(...){
  nm <- list(...) 
  dfdetect <- grepl("data.frame|matrix", unlist(lapply(nm, function(cl) paste(class(cl), collapse = " ") )))
  # first cbind vectors together 
  vec <- data.frame(nm[!dfdetect])
  n <- max(sapply(nm[dfdetect], nrow)) 
  vec <- data.frame(lapply(vec, function(x) rep(x, n)))
  if (nrow(vec) > 0) nm <- c(nm[dfdetect], list(vec))
  nm <- lapply(nm, as.data.frame)
  
  do.call(cbind, lapply(nm, function (df1) 
    rbind(df1, as.data.frame(matrix(NA, ncol = ncol(df1), nrow = n-nrow(df1), dimnames = list(NULL, names(df1))))) )) 
}

create_density_heatmaps <- function(df_list, x, drug_name_list) {
  heatmaps_list <- list()

  # Loop through each drug and create a density heatmap
  for (drug_name in drug_name_list) {
    l1 <- subset(df_list$FIMM, drug == drug_name, select=c(x))
    colnames(l1)[colnames(l1) == x] <- 'FIMM'
    
    l2 <- subset(df_list$Enserink, drug == drug_name, select=c(x))
    colnames(l2)[colnames(l2) == x] <- 'Enserink'
    
    l3 <- subset(df_list$BeatAML, drug == drug_name, select=c(x))
    colnames(l3)[colnames(l3) == x] <- 'BeatAML'
    
    l4 <- subset(df_list$Karolinska, drug == drug_name, select=c(x))
    colnames(l4)[colnames(l4) == x] <- 'Karolinska'
    
    combined_df <- cbind.fill(l1,l4, l2, l3)
    print("Check this out")
    
    print(combined_df)
    # Convert the wide dataframe to a matrix for densityHeatmap
    lab_matrix <- as.matrix(combined_df)
    
    # Create a densityHeatmap object
    heatmap_obj <- densityHeatmap(
      lab_matrix, title = drug_name
    )
    print(heatmap_obj)
    # Store the heatmap object in the list
    heatmaps_list[[drug_name]] <- heatmap_obj
  }
  
  # Return the list of heatmap objects (optional)
  return(heatmaps_list)
}

# Create density heatmaps for the combined datasets
density_heatmaps <- create_density_heatmaps(list_common_dfs, 'DSS2', drug_names_to_include)



#same drug in overlapping graphs??
#Density plot same drug


#How do I need to do the batch correction if there is only one drug???
#testing for controls as well as drugs and see if I see any differences
#Also plot variance of data 
#Understand my plots and what I should present from them
#Time versus DSS2 
#Number of negative controls versus DSS2
#Number of positive controls versus DSS2
#Concentration max or min per drug versus DSS2
#White blood cell count versus DSS2
#Batch correction then check if still significant differences for DSS2 per drug 

#All samples
model_interpretation(fimm_beataml_filtered, "DSS2", fimm_beataml_mixed_mode_all_vars, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2.html", header_1 = "FIMM BeatAML All samples")
View(fimm_beataml_filtered)

#Fresh only
fimm_beataml_mixed_mode_fresh_only <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, sample.x == 'fresh'))
model_interpretation(subset(fimm_beataml_filtered, sample.x == 'fresh'), "DSS2", fimm_beataml_mixed_mode_fresh_only, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_fresh_only.html", header_1 = "FIMM BeatAML Fresh samples")

#Blood only 
fimm_beataml_mixed_mode_blood_only <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Tissue == 'Blood'))
model_interpretation(subset(fimm_beataml_filtered, Tissue == 'Blood'), "DSS2", fimm_beataml_mixed_mode_blood_only, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_blood_only.html", header_1 = "FIMM BeatAML Blood samples")

#Blood and Fresh 
fimm_beataml_mixed_mode_blood_and_fresh <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Tissue == 'Blood' & sample.x == 'fresh'))
model_interpretation(subset(fimm_beataml_filtered, Tissue == 'Blood' & sample.x == 'fresh'), "DSS2", fimm_beataml_mixed_mode_blood_and_fresh, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_blood_and_fresh.html", header_1 = "FIMM BeatAML Blood and fresh samples")

#Bone Marrow only
fimm_beataml_mixed_mode_bone_marrow_only <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Tissue == 'Bone marrow'))
model_interpretation(subset(fimm_beataml_filtered, Tissue == 'Bone marrow'), "DSS2", fimm_beataml_mixed_mode_bone_marrow_only, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_bone_marrow_only.html", header_1 = "FIMM BeatAML Bone marrow samples")

#Bone Marrow and Fresh
fimm_beataml_mixed_mode_bone_marrow_and_fresh <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Tissue == 'Bone marrow' & sample.x == 'fresh'))
model_interpretation(subset(fimm_beataml_filtered, Tissue == 'Bone marrow' & sample.x == 'fresh'), "DSS2", fimm_beataml_mixed_mode_bone_marrow_and_fresh, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_bone_marrow_and_fresh.html", header_1 = "FIMM BeatAML Bone marrow and Fresh samples")

#Diagnosis
fimm_beataml_mixed_mode_Diagnosis_only <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Disease.status == 'Diagnosis'))
model_interpretation(subset(fimm_beataml_filtered, Disease.status == 'Diagnosis'), "DSS2", fimm_beataml_mixed_mode_Diagnosis_only, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_Diagnosis_only.html", header_1 = "FIMM BeatAML Diagnosis samples")

#Diagnosis and bone marrow 
fimm_beataml_mixed_mode_diagnosis_and_bone_marrow <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Disease.status == 'Diagnosis' & Tissue == 'Bone marrow'))
model_interpretation(subset(fimm_beataml_filtered, Disease.status == 'Diagnosis' & Tissue == 'Bone marrow'), "DSS2", fimm_beataml_mixed_mode_diagnosis_and_bone_marrow, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_diagnosis_and_bone_marrow.html", header_1 = "FIMM BeatAML Diagnosis and bone marrow samples")

#Diagnosis and blood 
fimm_beataml_mixed_mode_diagnosis_and_blood <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Disease.status == 'Diagnosis' & Tissue == 'Blood'))
model_interpretation(subset(fimm_beataml_filtered, Disease.status == 'Diagnosis' & Tissue == 'Blood'), "DSS2", fimm_beataml_mixed_mode_diagnosis_and_blood, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_diagnosis_and_blood.html", header_1 = "FIMM BeatAML Diagnosis and blood samples")


#Relapse
fimm_beataml_mixed_mode_relapse_only <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Disease.status == 'Relapse'))
model_interpretation(subset(fimm_beataml_filtered, Disease.status == 'Relapse'), "DSS2", fimm_beataml_mixed_mode_relapse_only, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_relapse_only.html", header_1 = "FIMM BeatAML Relapse samples")

#Relapse and bone marrow
fimm_beataml_mixed_mode_relapse_and_bone_marrow <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Disease.status == 'Relapse' & Tissue == 'Bone marrow'))
model_interpretation(subset(fimm_beataml_filtered, Disease.status == 'Relapse' & Tissue == 'Bone marrow'), "DSS2", fimm_beataml_mixed_mode_relapse_and_bone_marrow, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_relapse_and_bone_marrow.html", header_1 = "FIMM BeatAML Relapse and bone marrow samples")

#Relapse and bone marrow and fresh
fimm_beataml_mixed_mode_relapse_and_bone_marrow_and_fresh <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Disease.status == 'Relapse' & Tissue == 'Bone marrow' & sample.x == 'fresh'))
model_interpretation(subset(fimm_beataml_filtered, Disease.status == 'Relapse' & Tissue == 'Bone marrow' & sample.x == 'fresh'), "DSS2", fimm_beataml_mixed_mode_relapse_and_bone_marrow_and_fresh, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_relapse_and_bone_marrow_and_fresh.html", header_1 = "FIMM BeatAML Relapse, bone marrow and fresh samples")

#Relapse and Blood
fimm_beataml_mixed_mode_relapse_and_blood <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Disease.status == 'Relapse' & Tissue == 'Blood'))
model_interpretation(subset(fimm_beataml_filtered, Disease.status == 'Relapse' & Tissue == 'Blood'), "DSS2", fimm_beataml_mixed_mode_relapse_and_blood, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_relapse_and_blood.html", header_1 = "FIMM BeatAML Relapse and blood samples")

#Relapse and Blood and fresh
fimm_beataml_mixed_mode_relapse_and_blood <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Disease.status == 'Relapse' & Tissue == 'Blood'))
model_interpretation(subset(fimm_beataml_filtered, Disease.status == 'Relapse' & Tissue == 'Blood'), "DSS2", fimm_beataml_mixed_mode_relapse_and_blood, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_relapse_and_blood.html", header_1 = "FIMM BeatAML Relapse and blood samples")


#Refractory
fimm_beataml_mixed_mode_Refractory_only <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), subset(fimm_beataml_filtered, Disease.status == 'Refractory'))
model_interpretation(subset(fimm_beataml_filtered, Disease.status == 'Refractory'), "DSS2", fimm_beataml_mixed_mode_Refractory_only, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/lab_vs_dss2_Refractory_only.html", header_1 = "FIMM BeatAML Refractory samples")

unique(fimm_beataml_filtered$Disease.status)





#----ANOVA analysis per drug ###################################################################################################################################################################################
library(ggstatsplot)

grouped_ggbetweenstats(
  data = all_datasets_v2,
  x = medium,
  y = DSS2_boxcox_sclaed2,
  ylab = "Box Cox Transformed DSS2 Scaled",
  grouping.var = drug,
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = TRUE, # ANOVA or Welch ANOVA
  plot.type = "box",
  p.adjust.method = "fdr",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE, 
  package = "ggsci",
  palette = "default_jco",
  ## arguments relevant for combine_plots
  annotation.args = list(title = ""),
  plotgrid.args = list(nrow = 4)
)


#----DSS2 ANALYSIS ###############################################################################################################################################################################################
filtered_df <- all_datasets_v2 %>%
  group_by(cells) %>%
  slice_head(n = 3) %>%
  ungroup()


df_transformed_dss <- all_datasets_v2 %>%
  group_by(drug) %>%
  mutate(scale_dss = scale(DSS2)) %>%
  ungroup %>%
  as.data.frame()

df_transformed_dss_mean <- all_datasets_v2 %>%
  group_by(drug, time_until_sample_usage) %>%
  mutate(mean_dss = mean(DSS2)) %>%
  ungroup %>%
  as.data.frame()

df_transformed_dss_mean <- df_transformed_dss_mean[, !(names(df_transformed_dss_mean) %in% c("DSS2", "Patient.num"))] %>%
  group_by(drug, lab) %>%
  slice(1) %>%
  as.data.frame()


for(j in unique(all_datasets_v2$drug)){
  all_datasets_v2$DSS2_sclaed[all_datasets_v2$drug==j] <-
      scale(all_datasets_v2$DSS2_boxcox[all_datasets_v2$drug==j])
}


ggplot(df_transformed_dss, aes(Patient.num, scale_dss)) +
  geom_boxplot() 


p <- ggplot(all_datasets_v2, aes(time_until_sample_usage, DSS2_sclaed, fill=time_until_sample_usage)) +
  geom_boxplot() +
  facet_wrap(~drug)
p <- ggplotly(p)
p

p <- ggplot(all_datasets_v2, aes(time_until_sample_usage, DSS2, fill=time_until_sample_usage)) +
  geom_boxplot() +
  facet_wrap(~drug + lab)
p <- ggplotly(p)
p

p <- ggplot(all_datasets_v2, aes(time_until_sample_usage, DSS2, fill=time_until_sample_usage)) +
  geom_boxplot() 
p <- ggplotly(p)
p


ggplot(all_datasets_v2, aes(x = time_until_sample_usage, y = DSS2_sclaed)) +
  geom_point() +
  facet_wrap(~drug)


#df_transformed_dss <- subset(df_transformed_dss, !is.na(log_dss) & is.finite(log_dss))
mixed_model_time_until_sample_usage_scale_dss2 <- lmer(DSS2_sclaed ~ time_until_sample_usage + (1 + time_until_sample_usage|drug) , all_datasets_v2)
check_model(mixed_model_time_until_sample_usage_scale_dss2)
summary(mixed_model_time_until_sample_usage_scale_dss2)
visualize(mixed_model_time_until_sample_usage_scale_dss2)



#Mean DSS
##----Mean DSS2----
df_transformed_dss_mean <- all_datasets_v2 %>%
  group_by(drug, time_until_sample_usage) %>%
  mutate(mean_dss = mean(DSS2)) %>%
  ungroup %>%
  as.data.frame()

df_transformed_dss_mean <- df_transformed_dss_mean[, !(names(df_transformed_dss_mean) %in% c("DSS2", "Patient.num"))] %>%
  group_by(drug, lab) %>%
  slice(1) %>%
  as.data.frame()


mixed_model_time_until_sample_usage_mean_dss2 <- lmer(mean_dss ~ time_until_sample_usage + (1|drug), df_transformed_dss_mean)
check_model(mixed_model_time_until_sample_usage_mean_dss2)
summary(mixed_model_time_until_sample_usage_mean_dss2)
plot(mixed_model_time_until_sample_usage_mean_dss2)
visualize(mixed_model_time_until_sample_usage_mean_dss2)

df_transformed_dss_mean <- all_datasets_v2 %>%
  group_by(drug, medium) %>%
  mutate(mean_dss = mean(DSS2)) %>%
  ungroup %>%
  as.data.frame()

df_transformed_dss_mean <- df_transformed_dss_mean[, !(names(df_transformed_dss_mean) %in% c("DSS2", "Patient.num"))] %>%
  group_by(drug, lab) %>%
  slice(1) %>%
  as.data.frame()

mixed_model_medium_mean_dss2 <- lmer(mean_dss ~ medium + (1|drug), data=df_transformed_dss_mean)
summary(mixed_model_medium_mean_dss2)
check_model(mixed_model_medium_mean_dss2)
filtered_data <- df_transformed_dss_mean[df_transformed_dss_mean$drug %in% c("Bortezomib", "PH-797804", "Ruxolitinib", "cytarabine"), ]
plot <- flexplot::visualize(mixed_model_medium_mean_dss2, plot="model", color= c("red", "blue", "green", "orange", "purple", 
                                                                                 "pink", "yellow", "brown", "cyan", "black"))
plot + scale_color_manual(values = c("red", "blue", "green", "orange", "purple", 
                                     "pink", "yellow", "brown", "cyan", "black")) +  
  geom_point(aes(shape = NULL), size=3) +            
  guides(shape = "none", line="none", color = guide_legend(nrow = 5))

df_transformed_dss_mean <- all_datasets_v2 %>%
  group_by(drug, cells) %>%
  mutate(mean_dss = mean(DSS2)) %>%
  ungroup %>%
  as.data.frame()

df_transformed_dss_mean <- df_transformed_dss_mean[, !(names(df_transformed_dss_mean) %in% c("DSS2", "Patient.num"))] %>%
  group_by(drug, lab) %>%
  slice(1) %>%
  as.data.frame()

mixed_model_cells_mean_dss2 <- lmer(mean_dss ~ cells + (1|drug), data=df_transformed_dss_mean)
summary(mixed_model_cells_mean_dss2)
check_model(mixed_model_cells_mean_dss2)


df_transformed_dss_mean <- all_datasets_v2 %>%
  group_by(drug, microenvironmental_stimuli) %>%
  mutate(mean_dss = mean(DSS2)) %>%
  ungroup %>%
  as.data.frame()

df_transformed_dss_mean <- df_transformed_dss_mean[, !(names(df_transformed_dss_mean) %in% c("DSS2", "Patient.num"))] %>%
  group_by(drug, lab) %>%
  slice(1) %>%
  as.data.frame()

mixed_model_microenvironmental_stimuli_mean_dss2 <- lmer(mean_dss ~ microenvironmental_stimuli + (1|drug), data=df_transformed_dss_mean)
summary(mixed_model_microenvironmental_stimuli_mean_dss2)
check_model(mixed_model_microenvironmental_stimuli_mean_dss2)


df_transformed_dss_mean <- all_datasets_v2 %>%
  group_by(drug, cell_counting_method) %>%
  mutate(mean_dss = mean(DSS2)) %>%
  ungroup %>%
  as.data.frame()

df_transformed_dss_mean <- df_transformed_dss_mean[, !(names(df_transformed_dss_mean) %in% c("DSS2", "Patient.num"))] %>%
  group_by(drug, lab) %>%
  slice(1) %>%
  as.data.frame()

mixed_model_cell_counting_method_mean_dss2 <- lmer(mean_dss ~ cell_counting_method + (1|drug), data=df_transformed_dss_mean)
summary(mixed_model_cell_counting_method_mean_dss2)
check_model(mixed_model_cell_counting_method_mean_dss2)


df_transformed_dss_mean <- all_datasets_v2 %>%
  group_by(drug, positive_control) %>%
  mutate(mean_dss = mean(DSS2)) %>%
  ungroup %>%
  as.data.frame()

df_transformed_dss_mean <- df_transformed_dss_mean[, !(names(df_transformed_dss_mean) %in% c("DSS2", "Patient.num"))] %>%
  group_by(drug, lab) %>%
  slice(1) %>%
  as.data.frame()

mixed_model_sensitivity_readout_mean_dss2 <- lmer(mean_dss ~ positive_control + (1|drug), data=df_transformed_dss_mean)
summary(mixed_model_sensitivity_readout_mean_dss2)
check_model(mixed_model_sensitivity_readout_mean_dss2)

df_transformed_dss_mean <- all_datasets_v2 %>%
  group_by(drug, centrifugation_procedure) %>%
  mutate(mean_dss = mean(DSS2)) %>%
  ungroup %>%
  as.data.frame()

df_transformed_dss_mean <- df_transformed_dss_mean[, !(names(df_transformed_dss_mean) %in% c("DSS2", "Patient.num"))] %>%
  group_by(drug, lab) %>%
  slice(1) %>%
  as.data.frame()

mixed_model_centrifugation_procedure_mean_dss2 <- lmer(mean_dss ~ centrifugation_procedure + (1|drug), data=df_transformed_dss_mean)
summary(mixed_model_centrifugation_procedure_mean_dss2)
check_model(mixed_model_centrifugation_procedure_mean_dss2)


df_transformed_dss_mean <- all_datasets_v2 %>%
  group_by(drug, plate_reader) %>%
  mutate(mean_dss = mean(DSS2)) %>%
  ungroup %>%
  as.data.frame()

df_transformed_dss_mean <- df_transformed_dss_mean[, !(names(df_transformed_dss_mean) %in% c("DSS2", "Patient.num"))] %>%
  group_by(drug, lab) %>%
  slice(1) %>%
  as.data.frame()

mixed_model_plate_reader_mean_dss2 <- lmer(mean_dss ~ plate_reader + (1|drug), data=df_transformed_dss_mean)
summary(mixed_model_plate_reader_mean_dss2)
check_model(mixed_model_plate_reader_mean_dss2)


mean_dss2_models <- c(mixed_model_time_until_sample_usage_mean_dss2, mixed_model_medium_mean_dss2, mixed_model_cells_mean_dss2, mixed_model_microenvironmental_stimuli_mean_dss2, mixed_model_cell_counting_method_mean_dss2, mixed_model_sensitivity_readout_mean_dss2, mixed_model_centrifugation_procedure_mean_dss2, mixed_model_plate_reader_mean_dss2)
df <- model_results(mean_dss2_models)
rownames(df) <- NULL
colnames(df) <- gsub("_", " ", colnames(df))     
colnames(df) <- sapply(colnames(df), tools::toTitleCase)
df <- df[order(-df$`Fixed Effect Coefficient`), ]
color_tile_custom <- formatter("span", 
                               style = function(x) {
                                 color_scale <- col_numeric(c("lightblue", "white", "lightblue"), domain = range(df$`Fixed Effect Coefficient`))
                                 formattable::style(display = "block", padding = "0 4px", `background-color` = color_scale(x))
                               }
)
ft <- formattable(df,  list(
  `Fixed Effect Term` = formatter("span", style = ~formattable::style(color = "blue")),
  `Fixed Effect Coefficient` = color_tile_custom,
  `p Value` = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  `Significance Stars` = formatter("span", style = ~formattable::style(color = "darkred"))
))
formattable(df, list(
  Fixed_Effect_Term = formatter("span", style = ~style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~style(color = "darkred"))
))



##----Box Cox transformations ----
#Search optimal response transformation based on Box-Cox Transformations for Linear Models
all_datasets_v2$DSS2_pos <- all_datasets_v2$DSS2 + 0.01
MASS::boxcox(DSS2_pos ~ time_until_sample_usage, 
             data = all_datasets_v2,
             lambda = seq(-0.25, 2, length.out = 10))

all_datasets_v2$DSS2_boxcox <- ((all_datasets_v2$DSS2)^0.5 - 1) / 0.5

#Standardize DSS scores w.r.t. labs and drugs
for(i in unique(all_datasets_v2$lab)){
  for(j in unique(all_datasets_v2$drug)){
    all_datasets_v2$DSS2_boxcox_sclaed1[all_datasets_v2$lab==i & all_datasets_v2$drug==j] <-
      scale(all_datasets_v2$DSS2_boxcox[all_datasets_v2$lab==i & all_datasets_v2$drug==j])
  }
}

for(j in unique(all_datasets_v2$drug)){
  all_datasets_v2$DSS2_boxcox_sclaed2[all_datasets_v2$drug==j] <-
    scale(all_datasets_v2$DSS2_boxcox[all_datasets_v2$drug==j])
}

all_response_metrics$DSS2_pos <- all_response_metrics$DSS2 + 0.01
MASS::boxcox(DSS2_pos ~ time_until_sample_usage, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$DSS2_boxcox <- ((all_response_metrics$DSS2)^0.5 - 1) / 0.5

#Standardize DSS scores w.r.t. labs and drugs
for(i in unique(all_response_metrics$lab)){
  for(j in unique(all_response_metrics$drug)){
    all_response_metrics$DSS2_boxcox_sclaed1[all_response_metrics$lab==i & all_datasets_v2$drug==j] <-
      scale(all_response_metrics$DSS2_boxcox[all_response_metrics$lab==i & all_datasets_v2$drug==j])
  }
}

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$DSS2_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$DSS2_boxcox[all_response_metrics$drug==j])
}


# library(ggstatsplot)
# custom_colors <- c("Beat AML" = "#8dd3c7",  "Oslo" = "#fdb462", "Helsinki" = "#fb8072", "Karolinska"= "#80b1d3")
# all_datasets_v2$lab <- factor(all_datasets_v2$lab, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
# selected_groups <- unique(all_datasets_v2$drug)[2:5] # Select first 6 groups
# filtered_data <- all_datasets_v2[all_datasets_v2$drug %in% selected_groups, ]
# filtered_data$drug <- tools::toTitleCase(filtered_data$drug)
# 
# g <- grouped_ggbetweenstats( # paired samples
#   data = all_response_metrics,
#   x = lab,
#   y = DSS2,
#   grouping.var = drug,
#   type = "nonparametric", # for wilcoxon
#   centrality.plotting = FALSE, # remove median
#   ggstatsplot.layer = TRUE,
#   xlab = "",        
#   ylab = expression("DSS"[2]),
#   results.subtitle = TRUE,  # Removes default subtitle
#   #subtitle = "{test} (p = {p})", 
#   p.adjust.method = "bonferroni",
#   pairwise.display = "significant", #"significant"
#   #pairwise.comparisons = TRUE,
#   plotgrid.args = list(ncol = 12),
#   facet_wrap.args = list(scales = "fixed", strip.position = "top"),
#   #plot.margin = margin(0, 10, 10, 10),
#   #point.args = list(alpha = 0.9),
#   violin.args = list(alpha = 0),
#   boxplot.args = list(alpha = 0),
#   ggsignif.args = list(textsize = 3, tip_length = 0.005, step_increase = 0.05),
#   p.value.label.args = list(
#     parse = TRUE
#   ),
#   ggplot.component = list(
#     coord_cartesian(ylim = c(0, 50)),
#     #scale_y_continuous(labels = c(0:5, "")),
#     scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))),
#     scale_color_manual(values = custom_colors), 
#     #ggplot2::scale_y_continuous(limits = c(0, 50), sec.axis = ggplot2::dup_axis(name = NULL)), 
#     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1))), 
#     theme(
#           panel.background = element_blank(), 
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 0)),
#           axis.text.y.right = ggplot2::element_blank(), 
#           axis.ticks.y.right = ggplot2::element_blank(),
#           axis.ticks.x = ggplot2::element_blank(),
#           axis.text.x = element_text(family = "Arial", face = "plain",size = 13, color = "black", angle = 25, hjust = 0.6),
#           axis.title.y = element_text(family = "Arial", face = "plain", size = 13, color = "black"),
#           plot.title = element_text(family = "Arial", face = "plain", size = 13, color = "black", hjust = 0.5, margin = margin(b = 0))))
#   ) 
# annot <- g[[1]]$layers[[4]]$stat_params$annotations
# annot <- gsub("list\\(italic\\(p\\)\\[\\\"FDR\\\"\\s*-\\s*adj\\.\\] == \"", "list(FDR == ", annot)
# annot <- gsub("\"\\)", ")", annot)
# 
# for (i in 1:46) {
#   print(g[[i]])
#   print(g[[1]]$layers[[4]]$stat_params$annotations)
#   
#   # Check if the current g[[i]] contains annotations
#   if (!is.null(g[[i]]$layers[[4]]$stat_params$annotations)) {
#     
#     # Extract the annotations for this g[[i]]
#     annotations <- g[[i]]$layers[[4]]$stat_params$annotations
#     
#     # Modify the annotations by removing unwanted parts
#     annotations_modified <- unlist(lapply(annotations, function(annot) {
#       # Remove LaTeX parts (italic(p) and adj.)
#       annot <- gsub("list\\(italic\\(p\\)\\[\\\"FDR\\\"\\s*-\\s*adj\\.\\] == \"", "list(FDR == ", annot)
#       annot <- gsub("\"\\)", ")", annot)  # Remove the trailing ')'
#       annot <- sapply(annot, function(x) {
#         # Extract the numeric value
#         num <- as.numeric(gsub(".*==\\s*([0-9.eE+-]+).*", "\\1", x))
#         
#         # Check and replace
#         if (num < 0.001) {
#           gsub("==.*", "< 0.001)", x)
#         } else {
#           x
#         }
#       })
#       return(annot)
#     }))
#     
#     # Update the annotations in g[[i]]
#     g[[i]]$layers[[4]]$stat_params$annotations <- annotations_modified
#   }
#   if(!is.null(g[[i]]$labels$subtitle)){
#     subtitle <- g[[i]]$labels$subtitle
#     print(subtitle)
#     subtitle_str <- deparse(subtitle)  
#     subtitle_str <- paste(subtitle_str, collapse = " ")  # Ensure it's a single string
#     
#     # Extract only "Kruskal-Wallis: italic(p) == value"
#     subtitle_mod <- str_replace(subtitle_str, 'list\\(chi\\["Kruskal-Wallis"\\]\\^2.*?italic\\(p\\) ==      "', "Kruskal-Wallis: p = ")
#     subtitle_mod <- str_replace(subtitle_mod, '", widehat\\(epsilon\\)\\["ordinal"\\]\\^2.*', "")
#     
#     print(subtitle_mod)
#     g[[i]]$labels$subtitle <- subtitle_mod
#   }
# }
# 
# gb <- ggplot_build(g)
# gb$data
# # Now print or display the modified plot
# pdf("Desktop/UiO/Project 1/Figures/draw/Difference_DSS2_per_drug_v1.pdf", width = 80, height = 30)
# ggplot_build(g)
# dev.off()
# 
# ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Difference_DSS2_per_drug_v1.png", plot = g, width = 80, height = 30, units = "cm", limitsize = FALSE)
# ggsave("Desktop/UiO/Project 1/Figures/draw/Difference_DSS2_per_drug_sample.png", plot = g, width = 8, height = 8, dpi = 300, limitsize = FALSE)
# 
# g <- grouped_ggbetweenstats( # paired samples
#   data = subset(filtered_data, drug != "001, RAD"),
#   x = lab,
#   y = DSS2,
#   grouping.var = drug,
#   type = "nonparametric", # for wilcoxon
#   centrality.plotting = FALSE, # remove median
#   xlab = "",        
#   ylab = "DSS2",
#   results.subtitle = TRUE,  # Removes default subtitle
#   #subtitle = "{test} (p = {p})", 
#   p.adjust.method = "fdr",
#   pairwise.display = "none", #"significant"
#   #pairwise.comparisons = TRUE,
#   plotgrid.args = list(nrow = 2, ncol = 2),
#   facet_wrap.args = list(scales = "fixed", strip.position = "top"),
#   plot.margin = margin(0, 10, 10, 10),
#   #point.args = list(alpha = 0.9),
#   violin.args = list(alpha = 0),
#   boxplot.args = list(alpha = 0),
#   #ggsignif.args = list(textsize = 3, tip_length = 0.005, step_increase = 0.05),
#   ggplot.component = list(
#     coord_cartesian(ylim = c(0, 50)),
#     #scale_y_continuous(labels = c(0:5, "")),
#     #scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))),
#     scale_color_manual(values = custom_colors), 
#     #ggplot2::scale_y_continuous(limits = c(0, 50), sec.axis = ggplot2::dup_axis(name = NULL)), 
#     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1))), 
#     theme(
#       panel.background = element_blank(), 
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 0)),
#       axis.text.y.right = ggplot2::element_blank(), 
#       axis.ticks.y.right = ggplot2::element_blank(),
#       axis.ticks.x = ggplot2::element_blank(),
#       axis.text.x = element_text(family = "Arial", face = "plain",size = 12, color = "black", angle = 25, hjust = 0.6), #0.6
#       axis.title.y = element_text(family = "Arial", face = "plain", size = 12, color = "black"),
#       plot.title = element_text(family = "Arial", face = "plain", size = 12, color = "black", hjust = 0.5, margin = margin(b = 0))))
# ) 
# 
# for (i in 1:4) {
#   if(!is.null(g[[i]]$labels$subtitle)){
#     subtitle <- g[[i]]$labels$subtitle
#     print(subtitle)
#     subtitle_str <- deparse(subtitle)  
#     subtitle_str <- paste(subtitle_str, collapse = " ")  # Ensure it's a single string
#     
#     # Extract only "Kruskal-Wallis: italic(p) == value"
#     subtitle_mod <- str_replace(subtitle_str, 'list\\(chi\\["Kruskal-Wallis"\\]\\^2.*?italic\\(p\\) ==      "', "p = ")
#     subtitle_mod <- str_replace(subtitle_mod, '", widehat\\(epsilon\\)\\["ordinal"\\]\\^2.*', "")
#     
#     print(subtitle_mod)
#     g[[i]]$labels$subtitle <- subtitle_mod
#   }
# }
# print(g)
# 
# 
# g <- grouped_ggbetweenstats( # paired samples
#   data = all_response_metrics,
#   x = lab,
#   y = DSS2_boxcox_sclaed2,
#   grouping.var = drug,
#   type = "nonparametric", # for wilcoxon
#   centrality.plotting = FALSE, # remove median
#   ggstatsplot.layer = TRUE,
#   xlab = "",        
#   ylab = "Box-Cox transformed Z-scaled DSS2",
#   results.subtitle = TRUE,  # Removes default subtitle
#   #subtitle = "{test} (p = {p})", 
#   p.adjust.method = "fdr",
#   pairwise.display = "significant", #"significant"
#   #pairwise.comparisons = TRUE,
#   #plotgrid.args = list(nrow = 2, ncol = 2),
#   facet_wrap.args = list(scales = "fixed", strip.position = "top"),
#   #plot.margin = margin(0, 10, 10, 10),
#   #point.args = list(alpha = 0.9),
#   violin.args = list(alpha = 0),
#   boxplot.args = list(alpha = 0),
#   ggsignif.args = list(textsize = 3, tip_length = 0.005, step_increase = 0.05),
#   ggplot.component = list(
#     coord_cartesian(ylim = c(-4, 4)),
#     #scale_y_continuous(labels = c(0:5, "")),
#     scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))),
#     scale_color_manual(values = custom_colors), 
#     #ggplot2::scale_y_continuous(limits = c(0, 50), sec.axis = ggplot2::dup_axis(name = NULL)), 
#     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1))), 
#     theme(
#       panel.background = element_blank(), 
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 0)),
#       axis.text.y.right = ggplot2::element_blank(), 
#       axis.ticks.y.right = ggplot2::element_blank(),
#       axis.ticks.x = ggplot2::element_blank(),
#       axis.text.x = element_text(family = "Arial", face = "plain",size = 13, color = "black", angle = 25, hjust = 0.6),
#       axis.title.y = element_text(family = "Arial", face = "plain", size = 13, color = "black"),
#       plot.title = element_text(family = "Arial", face = "plain", size = 13, color = "black", hjust = 0.5, margin = margin(b = 0))))
# ) 
# 
# for (i in 1:46) {
#   print(g[[i]])
#   print(g[[1]]$layers[[4]]$stat_params$annotations)
#   
#   # Check if the current g[[i]] contains annotations
#   if (!is.null(g[[i]]$layers[[4]]$stat_params$annotations)) {
#     
#     # Extract the annotations for this g[[i]]
#     annotations <- g[[i]]$layers[[4]]$stat_params$annotations
#     
#     # Modify the annotations by removing unwanted parts
#     annotations_modified <- unlist(lapply(annotations, function(annot) {
#       # Remove LaTeX parts (italic(p) and adj.)
#       annot <- gsub("list\\(italic\\(p\\)\\[\\\"FDR\\\"\\s*-\\s*adj\\.\\] == \"", "list(FDR == ", annot)
#       annot <- gsub("\"\\)", ")", annot)  # Remove the trailing ')'
#       return(annot)
#     }))
#     
#     # Update the annotations in g[[i]]
#     g[[i]]$layers[[4]]$stat_params$annotations <- annotations_modified
#   }
#   if(!is.null(g[[i]]$labels$subtitle)){
#     subtitle <- g[[i]]$labels$subtitle
#     print(subtitle)
#     subtitle_str <- deparse(subtitle)  
#     subtitle_str <- paste(subtitle_str, collapse = " ")  # Ensure it's a single string
#     
#     # Extract only "Kruskal-Wallis: italic(p) == value"
#     subtitle_mod <- str_replace(subtitle_str, 'list\\(chi\\["Kruskal-Wallis"\\]\\^2.*?italic\\(p\\) ==      "', "Kruskal-Wallis: p-value = ")
#     subtitle_mod <- str_replace(subtitle_mod, '", widehat\\(epsilon\\)\\["ordinal"\\]\\^2.*', "")
#     
#     print(subtitle_mod)
#     g[[i]]$labels$subtitle <- subtitle_mod
#   }
# }
# 
# # Now print or display the modified plot
# print(g)
# g[[1]]$labels$subtitle
# ggsave("Desktop/UiO/Project 1/Figures/draw/Difference_box_cox_sclaed_DSS2_per_drug.png", plot = g, width = 50, height = 90, units = "cm", limitsize = FALSE)
# 
#   
# s <- grouped_ggbetweenstats( # paired samples
#   data = all_datasets_v2,
#   x = lab,
#   y = DSS2_boxcox_sclaed2,
#   grouping.var = drug,
#   type = "nonparametric", # for wilcoxon
#   centrality.plotting = FALSE # remove median
# )
# extract_stats(g)
# 
# ggsave("Desktop/UiO/Project 1/Figures/V3/Difference_DSS2_per_drug_box_cox_scaled.pdf", plot = s, width = 150, height = 250, units = "cm", limitsize = FALSE)

overall_lab_difference_model <- lmer(DSS2_boxcox_sclaed2 ~ drug + (1|Patient.num) + (lab|drug), all_datasets_v2)
summary(overall_lab_difference_model)
check_model(overall_lab_difference_model)

model <- mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept
  residuals_marginal <- residuals(model, type = "pearson")
  fitted_values <- fitted(model)
  all_datasets_v2$residuals_marginal <- residuals(model, type = "pearson")
  all_datasets_v2$fitted_values <- fitted(model)

#----Visualization good!!----
require(flexplot)
data <- subset(all_datasets_v2, !is.na(medium))
print(unique(data$medium))
data <- data %>% mutate(Medium = case_when(medium == "RPMI + fetal bovine serum (FBS) (10%)" ~ "RPMI + FBS", 
                                           medium ==  "HS-5 conditioned medium" ~ "HS-5 CM",
                                           medium ==  "Mononuclear cell medium" ~ "MCM",                                                                                                              
                                          TRUE ~ medium))   
print(unique(data$Medium))
model <- lmer(DSS2_boxcox_sclaed2 ~ Medium +(Medium|drug), data)
estimates(model)
plot <- flexplot::visualize(model, plot = "model", sample = 3)#, center = "Median + quartiles") 

plot$layers[[1]]$aes_params$alpha <- 0.5  #dot size
plot$layers[[2]]$aes_params$alpha <- 0.5  #dot size
plot$layers[[3]]$aes_params$linewidth <- 1 #linesize
plot$layers[[3]]$aes_params$linetype <- 2  #line
plot$layers[[3]]$aes_params$alpha <- 1
plot$layers[[4]]$aes_params$size <- 1.5
plot$layers[[5]]$aes_params$alpha <- 1 #average dot size
plot$layers[[5]]$aes_params$linewidth <- 1 #average line size
print(plot$layers)

plot <- plot + 
  scale_color_manual(values = c("#80b1d3", "#fb8072", "#bebada","#fccde5"), guide = guide_legend(title = "Drug")) +  
  labs(
    title = "Model Visualization",
    x = "Experimental Variable", 
    y = "Scaled DSS2", 
    color = "Drug"
  ) + coord_cartesian(ylim = c(-3, 3)) +   
  guides(
    linetype = "none", 
    shape = "none"
  ) + theme(
    legend.position = "top",  # Position legend on the right
    legend.title = element_text(family = "Arial", size = 10, color = "black"),  # Legend title Arial size 12
    legend.text = element_text(family = "Arial", size = 10, color = "black"),   # Legend content Arial size 12
    axis.text = element_text(family = "Arial", size = 10, color = "black"),     # Axis text Arial size 12
    axis.title = element_text(family = "Arial", size = 10, color = "black"),    # Axis titles Arial size 12
    axis.text.y = element_text(family = "Arial", size = 8, color = "black"),                    
    axis.text.x = element_text(family = "Arial", size = 8, color = "black"),    
    panel.grid = element_blank(),
    plot.title = element_text(
      size = 10,            # Font size
      face = "plain",        # Bold text
      hjust = 0.5,        
      family = "Arial", 
      vjust = 1
    )
  )

ggsave('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/Model.png', plot = plot, width=9, height=6.5, dpi = 300, unit = "cm")

model1 <- lmer(DSS2_boxcox_sclaed2 ~ Medium + (1|Patient.num) + (Medium|drug), data)

data$y <- predict(model1, newdata = data, re.form = ~ (1|Patient.num) + (Medium|drug))

agg_data <- data %>%
  group_by(drug, Medium) %>%
  dplyr::summarize(mean_y = mean(y), .groups = "drop")

agg_data_filtered <- agg_data %>%
  filter(Medium %in% c("RPMI + FBS", "MCM"))

plot <- ggplot(subset(agg_data_filtered, drug == 'Vismodegib' | drug == 'Crizotinib' | drug == '345627-80-7'), aes(x = Medium, y = mean_y, group = drug, color = drug)) +
  geom_point(alpha = 0) +  # Add points but make them invisible (alpha = 0)
  geom_line(linewidth = 1.2) +  # Add lines for each subject
  scale_color_manual(values = c("#80b1d3", "#fb8072", "#bebada")) +
  scale_x_discrete(
    limits = c("RPMI + FBS", "MCM"),  # Show only the desired categories
    #labels = c("Group A", "Group C")  # Custom labels for x-axis
  ) +
  labs(
    #title = "Random Intercepts and Slopes (Categorical X)",
    x = "Experimental Variable",
    y = "Response"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Position legend on the right
    #legend.title = element_text(family = "Arial", size = 17),  # Legend title Arial size 12
    #legend.text = element_text(family = "Arial", size = 17),   # Legend content Arial size 12
    axis.text = element_text(family = "Arial", size = 10, color = "black"),     # Axis text Arial size 12
    axis.title = element_text(family = "Arial", size = 10),    # Axis titles Arial size 12
    axis.text.y = element_text(family = "Arial", size = 8, color = "black"),                    
    axis.text.x = element_text(family = "Arial", size = 8, color = "black"),    
    plot.margin = margin(0, 0, 0, 0),
    panel.grid = element_blank()  # Remove grid lines
  )
plot
ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/ransom_slope_and_intercept.png", plot, width = 6.39, height = 3.54, units="cm", dpi = 300)


filtered_data <- subset(all_datasets_v2, !is.na(medium) & drug %in% c("Bortezomib", "PH-797804", "Ruxolitinib", "cytarabine"))
midpoints <- filtered_data %>%
  group_by(medium, drug) %>%
  summarise(mid_value = median(DSS2_boxcox_sclaed2), .groups = "drop") %>%
  arrange(medium, drug)
#midpoints$x_position <- as.numeric(factor(midpoints$medium))
average_line <- filtered_data %>%
  group_by(medium) %>%
  dplyr::summarize(avg_value = mean(DSS2_boxcox_sclaed2, na.rm = TRUE)) %>%
  ungroup()

average_line <- average_line %>%
  mutate(x_position = as.numeric(factor(medium)) - 0.65)  # Adjust position on x-axis for category


midpoints <- midpoints %>%
  mutate(x_position = as.numeric(factor(medium)) + (as.numeric(factor(drug))) * 0.25 - 0.65,
         # Apply a slight shift for the books that need adjustment
         x_position = case_when(
           drug == "Bortezomib" ~ x_position + 0.05,  # Adjust Book3 to the right
           drug == "cytarabine" ~ x_position + 0.05,  # Adjust Book4 to the right
           TRUE ~ x_position  # No change for other books
         ))
ggplot(filtered_data, aes(x = medium, y = DSS2_boxcox_sclaed2, color=drug, fill = drug)) +
  geom_violin(trim = FALSE) +   # trim = FALSE to show the full distribution
    scale_fill_manual(values = c("Bortezomib" = "#fccde5", "PH-797804" = "#b3de69", "Ruxolitinib" = "#ffffb3", "cytarabine" = "#bebada"), 
                      labels = c("Drug 1", "Drug 2", "Drug 3", "Drug 4")) +  # Manual colors
    scale_color_manual(values = c("Bortezomib" = "#fccde5", "PH-797804" = "#b3de69", "Ruxolitinib" = "#ffffb3", "cytarabine" = "#bebada"), 
                       labels = c("Drug 1", "Drug 2", "Drug 3", "Drug 4")) + # Manual colors
  geom_segment(data = midpoints %>% 
                 group_by(drug) %>%
                 arrange(x_position) %>%
                 mutate(next_x = lead(x_position), next_y = lead(mid_value)) %>%
                 filter(!is.na(next_x)),
               aes(x = x_position, xend = next_x, y = mid_value, yend = next_y, color = drug), size = 1) +
  theme_minimal() +  
    labs(x = "Experimental variable", y = "DSS") + 
  theme(
    legend.title = element_blank(),
    axis.text.x = element_blank(),  # Remove x-axis labels (A, B, C, D)
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis numbers
    axis.ticks.y = element_blank()   # Remove y-axis ticks
  ) 

# df_avg <- filtered_data %>%
#   group_by(medium, lab) %>%
#   dplyr::summarize(DSS2_boxcox_sclaed2 = mean(DSS2_boxcox_sclaed2, na.rm = TRUE)) %>%
#   mutate(drug = "Average")  # Add a new Book called "Average"
# 
# midpoints <- df_combined %>%
#   group_by(medium, drug) %>%
#   summarise(mid_value = median(DSS2_boxcox_sclaed2), .groups = "drop") %>%
#   arrange(medium, drug)
# 
# midpoints <- midpoints %>%
#   mutate(x_position = as.numeric(factor(medium)) + (as.numeric(factor(drug))) * 0.25 - 0.65,
#          # Apply a slight shift for the books that need adjustment
#          x_position = case_when(
#            drug == "Bortezomib" ~ x_position + 0.05,  # Adjust Book3 to the right
#            drug == "cytarabine" ~ x_position + 0.05,  # Adjust Book4 to the right
#            TRUE ~ x_position  # No change for other books
#          ))
# # Step 2: Combine the original dataset with the average drug data
# df_combined <- bind_rows(filtered_data, df_avg)
# 
# ggplot(filtered_data, aes(x = medium, y = DSS2_boxcox_sclaed2, color=drug, fill = drug)) +
#   geom_violin(trim = FALSE) +   # trim = FALSE to show the full distribution
#   scale_fill_manual(values = c("Bortezomib" = "#fccde5", "PH-797804" = "#b3de69", "Ruxolitinib" = "#ffffb3", "cytarabine" = "#bebada", "Average" = "grey"), 
#                     labels = c("Drug 1", "Drug 2", "Drug 3", "Drug 4")) +  # Manual colors
#   scale_color_manual(values = c("Bortezomib" = "#fccde5", "PH-797804" = "#b3de69", "Ruxolitinib" = "#ffffb3", "cytarabine" = "#bebada", "Average" = "grey"), 
#                      labels = c("Drug 1", "Drug 2", "Drug 3", "Drug 4")) + # Manual colors
#   geom_segment(data = midpoints %>% 
#                  group_by(drug) %>%
#                  arrange(x_position) %>%
#                  mutate(next_x = lead(x_position), next_y = lead(mid_value)) %>%
#                  filter(!is.na(next_x)),
#                aes(x = x_position, xend = next_x, y = mid_value, yend = next_y, color = drug), size = 1) +
#   theme_minimal() +  
#   labs(x = "Experimental variable", y = "DSS") + 
#   theme(
#     legend.title = element_blank(),
#     axis.text.x = element_blank(),  # Remove x-axis labels (A, B, C, D)
#     axis.ticks.x = element_blank(),  # Remove x-axis ticks
#     axis.text.y = element_blank(),  # Remove y-axis numbers
#     axis.ticks.y = element_blank()   # Remove y-axis ticks
#   ) 

# 2. Homogeneity of Residual Variance: Plot residuals vs. fitted values
  p1 <- ggplot(data = data.frame(residuals = residuals_marginal, fitted = fitted_values),
               aes(x = fitted, y = residuals)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    ggtitle("Marginal Residuals vs. Fitted Values") +
    xlab("Fitted Values") +
    ylab("Marginal Residuals")
  
  qqnorm(residuals_marginal, ylab = "Marginal residuals", xlab = "Theoretical Quantiles", main = "Q-Q Plot of Marginal residuals")
  qqline(residuals_marginal, col = 2,lwd=2,lty=2)
  
  ggplot(all_datasets_v2,
         aes(x = time_until_sample_usage, y = residuals_marginal)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    ggtitle("Marginal Residuals vs. Independent Var") +
    xlab("Fitted Values") +
    ylab("Marginal Residuals")
  
  
  # 3. Linearity: Level-1 residuals versus independent variables
  level1_residuals <- residuals(model, type = "response")
  independent_var <- model.frame(model)[[2]]  # Assuming the independent variable is the second column
  all_datasets_v2$level1_residuals <- residuals(model, type = "response")
  
  p2 <- ggplot(all_datasets_v2,
               aes(x = level1_residuals, y = fitted_values)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    ggtitle("Level-1 Residuals vs. Fitted") +
    xlab("Fitted") +
    ylab("Level-1 Residuals") + 
    facet_wrap(~drug)
  
  ggplot(all_datasets_v2,
         aes(x = time_until_sample_usage, y = level1_residuals)) +
    geom_point() +
    ggtitle("Level-1 Residuals vs. Independent variable") +
    xlab("Experimental Variable") +
    ylab("Level-1 Residuals") + 
    facet_wrap(~drug)
  
  library(lattice)
  # 4. Normality of random effects (level-2 residuals)
  re <- ranef(model, condVar = TRUE)
  ggplot(as.matrix(re$drug['(Intercept)'])) +
    stat_qq() +
    stat_qq_line() +
    ggtitle("Q-Q Plot of Random Effects")
  
  qqnorm(as.matrix(re$drug['(Intercept)']), ylab = "Random Intercept - drug", xlab = "Theoretical Quantiles", main = "Q-Q Plot of Random Effects")
  qqline(as.matrix(re$drug['(Intercept)']))
  
  qqnorm(as.matrix(re$drug['time_until_sample_usageHandled within 2-72 hours']), ylab = "Random Intercept - slope per drug", xlab = "Theoretical Quantiles", main = "Q-Q Plot of Random Effects")
  qqline(as.matrix(re$drug['time_until_sample_usageHandled within 2-72 hours']))
  
  qqnorm(as.matrix(re$Patient.num['(Intercept)']), ylab = "Random Intercept - Patients", xlab = "Theoretical Quantiles", main = "Q-Q Plot of Random Effects")
  qqline(as.matrix(re$Patient.num['(Intercept)']))
  
  # 5. Normality of level-1 residuals (conditional residuals)
  p3 <- ggplot(data = data.frame(residuals = level1_residuals),
               aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line() +
    ggtitle("Q-Q Plot of Level-1 Residuals")
  
  # 6. Homogeneity: Check within-group variance homogeneity
  residuals_df <- data.frame(residuals = level1_residuals, groups = model@flist[[1]])
  
  p4 <- ggplot(residuals_df, aes(x = groups, y = residuals)) +
    geom_boxplot() +
    ggtitle("Residuals by Groups (Homogeneity of Variance)") +
    xlab("Groups") +
    ylab("Residuals")


time_linear <- lm(DSS2_boxcox_sclaed2 ~ time_until_sample_usage, all_response_metrics)
time_mixed_1 <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1|drug), all_response_metrics)
time_mixed_2 <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient.num), all_response_metrics)
time_mixed_3 <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1 + time_until_sample_usage|drug), all_response_metrics)
time_mixed_final <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient.num) + (1 + time_until_sample_usage|drug), all_response_metrics)
time_mixed_4 <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient.num) + (1|drug), all_response_metrics)
time_mixed_5 <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient.num) + (1|drug) + (0 + time_until_sample_usage|drug), all_response_metrics)
time_mixed_6 <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient.num) + (0 + time_until_sample_usage|drug), all_response_metrics)

AIC(time_linear)
BIC(time_linear)
anova(time_mixed_1, time_mixed_2, time_mixed_3, time_mixed_final, time_mixed_4, time_mixed_5, time_mixed_6)
icc(time_mixed_1)
icc(time_mixed_2)
icc(time_mixed_3)
icc(time_mixed_final)
icc(time_mixed_4)
icc(time_mixed_5)
icc(time_mixed_6)

flexplot::model.comparison(time_mixed_1, time_mixed_2, time_mixed_3, time_mixed_final, time_mixed_4, time_mixed_5, time_mixed_6)

###----(1 | Patient.num) + (E_var | drug)----
####----Time_until_sample_usage----
mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient_ID) + (1 + time_until_sample_usage|drug), all_response_metrics, control = lmerControl(optimizer = "bobyqa"))
summary(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
model_diagnostics(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/DSS2/time_until_sample_usage/") 
  
diagnostics <- check_model(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, panel = FALSE)
confint(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
icc(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
summary(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
r2(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
check_overdispersion(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
check_singularity(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
check_heteroscedasticity(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
model_performance(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
library(brms)  # or any other modeling package you're using, e.g., lme4, etc.
pp_check(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
check_linearity(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
check_heteroscedasticity(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
outlier_list <- check_outliers(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
insight::get_data(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)[outliers_list, ]

check <- check_model(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, panel = F,title_size = 10, title_color = "black", title_hjust = 0.5, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3", "black"))
plot(check)


all_datasets_v2$medium <- as.factor(all_datasets_v2$medium)
medium_model <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1|Patient.num) + (medium|drug), subset(all_datasets_v2, Patient.num != '2493_BA2563D_BA2563R'), REML = TRUE)
pp_check_plot <- check_model(medium_model, check = "pp_check", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#80b1d3", "#fccde5", "black"))
plot_pp_check <- plot(pp_check_plot)
plot_pp_check_adjusted <- plot_pp_check$PP_CHECK + theme(text = element_text(family = "Arial", color= "black", size = 10),
                               plot.title = element_text(hjust = 0.5, vjust = 1),
                               axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                               axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                               axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                               axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                               plot.subtitle = element_blank(), 
                               legend.position = "top",  
                               legend.text = element_text(family = "Arial", color = "black", size = 10) 
                               ) + labs(x = expression("BoxCox Transformed and Scaled DSS"[2]))
plot_pp_check_adjusted
ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/pp_check.png", plot = plot_pp_check_adjusted,width = 9, height = 6.5, units="cm") #width = 9, height = 6.5
performance_accuracy(medium_model)
print(plot$layers)

linearity_plot <- check_model(medium_model, check = "linearity", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
plot_linearity <- plot(linearity_plot)
plot_linearity_adjusted <- plot_linearity$NCV + theme(text = element_text(family = "Arial", color= "black", size = 10),
                               plot.title = element_text(hjust = 0.5, vjust = 1),
                               axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                               axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                               axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                               axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                               plot.subtitle = element_blank(), 
                               legend.position = "right",  
                               legend.text = element_text(family = "Arial", color = "black", size = 10) 
)
plot_linearity_adjusted
ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/linearity.png", plot = plot_linearity_adjusted, width = 9, height = 6.5, units="cm")

homogeneity_plot <- check_model(medium_model, check = "homogeneity", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
plot_homogeneity <- plot(homogeneity_plot)
plot_homogeneity_adjusted <- plot_homogeneity$HOMOGENEITY + labs(y = expression(sqrt("|Std. Residuals|"))) + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                      plot.title = element_text(hjust = 0.5, vjust = 1),
                                                      axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                      axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                      axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                      axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                      plot.subtitle = element_blank(), 
                                                      show.legend = TRUE,
                                                      legend.key = element_rect(fill = "white"),
                                                      legend.position = "top",  
                                                      legend.text = element_text(family = "Arial", color = "black", size = 10) 
)
plot_homogeneity_adjusted$layers[[1]]$

ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/homogeneity_check.png", plot = plot_homogeneity_adjusted, width = 9, height = 6.5, units="cm")

influential_plot <- check_model(medium_model, check = "outliers", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3", "black"))
plot_influential <- plot(influential_plot)  
plot_influential_adjusted <- plot_influential$OUTLIERS + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                                  plot.title = element_text(hjust = 0.5, vjust = 1),
                                                                  axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                                  axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                                  axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                                  axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                                  plot.subtitle = element_blank(), 
                                                                  legend.position = "right",  
                                                                  legend.text = element_text(family = "Arial", color = "black", size = 10) 
)
plot_influential_adjusted
outlier_list <- check_outliers(medium_model)
outliers_info <- as.data.frame(outlier_list)
ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/influential_check.png", plot = plot_influential_adjusted, width = 9, height = 6.5, units="cm")

normality_plot <- check_model(medium_model, check = "qq", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
plot_normality <- plot(normality_plot)
plot_normality_adjusted <- plot_normality$QQ + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                               plot.title = element_text(hjust = 0.5, vjust = 1),
                                                               axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                               axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                               axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                               axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                               plot.subtitle = element_blank(), 
                                                               legend.position = "right",  
                                                               legend.text = element_text(family = "Arial", color = "black", size = 10) 
)
plot_normality_adjusted
ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/normality_check.png", plot = plot_normality_adjusted, width = 7.94, height = 7.49, units="cm")

Medium_model1 <- lmer(DSS2_boxcox_sclaed2 ~ Medium + (1|Patient.num) + (Medium|drug),data)
random_effects_plot <- check_model(medium_model, check = "reqq", panel = F, title_size = 10, axis_title_size = 10, base_size = 8, colors = c("#fccde5", "#80b1d3"))
random_effects_plot$REQQ$drug$facet <- factor(random_effects_plot$REQQ$drug$facet, 
                                        levels = c("(Intercept)", "mediumMononuclear cell medium", "mediumRPMI + fetal bovine serum (FBS) (10%)"), 
                                        labels = c("Intercept", "MCM", "RPMI + FBS"))

plot_random_effects <- plot(random_effects_plot)
plot_random_effects[[2]]
plot_random_effects_adjusted_1 <- plot_random_effects[[1]] + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                     plot.title = element_text(hjust = 0.5, vjust = 1),
                                                     axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                     axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                     axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                     axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                     plot.subtitle = element_blank(), 
                                                     legend.position = "right",  
                                                     legend.text = element_text(family = "Arial", color = "black", size = 10), 
                                                     panel.grid.major = element_blank(), 
                                                     panel.grid.minor = element_blank(),  
                                                     panel.background = element_blank()  
) + labs(title = "Normality of Random Effects - Patient") 
plot_random_effects_adjusted_1
ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/random_effects_check_1.png", plot = plot_random_effects_adjusted_1, width = 9, height = 6.5, units="cm")

plot_random_effects_adjusted_2 <- plot_random_effects[[2]] + theme(text = element_text(family = "Arial", color= "black", size = 10),
                                                                   plot.title = element_text(hjust = 0.5, vjust = 1,color= "black"),
                                                                   axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                                                                   axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                                                                   axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                                                                   axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                                                                   plot.subtitle = element_blank(), 
                                                                   legend.position = "right",  
                                                                   legend.text = element_text(family = "Arial", color = "black", size = 10),
                                                                   panel.grid.major = element_blank(), 
                                                                   panel.grid.minor = element_blank(),  
                                                                   panel.background = element_blank(),
                                                                   strip.text = element_text(family = "Arial", color = "black", size = 8, face = "plain"),  
                                                                   strip.background = element_blank()
) + labs(title = "Normality of Random Effects - Drug")
plot_random_effects_adjusted_2
ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/random_effects_check_2.png", plot = plot_random_effects_adjusted_2, width = 9, height = 6.5, units="cm")

library(patchwork)
((plot | plot_homogeneity_adjusted) / (plot_pp_check_adjusted | plot_influential_adjusted) / (plot_random_effects_adjusted_1 | plot_random_effects_adjusted_2) / patchwork::wrap_elements(p)) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) 

p_wrapped <- patchwork::wrap_elements(p)

p_wrapped <- wrap_elements(p)

comboplot1 <- (
  (plot | plot_homogeneity_adjusted)
) +
  plot_layout(heights = c(1, 1, 1, 2))

comboplot2 <- (
    (plot_pp_check_adjusted | plot_influential_adjusted) 
) +
  plot_layout(heights = c(1, 1, 1, 2))

comboplot3 <- (
    (plot_random_effects_adjusted_1 | plot_random_effects_adjusted_2) 
) +
  plot_layout(heights = c(1, 1, 1, 2))

ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/combo_plots.png", plot = comboplot1, width = 20.76, height = 7.09, units="cm")
ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/combo_plots2.png", plot = comboplot2, width = 20.76, height = 7.09, units="cm")
ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/combo_plots3.png", plot = comboplot3, width = 20.76, height = 7.09, units="cm")



random_effect_plots <- (plot_random_effects_adjusted_1 | plot_random_effects_adjusted_2)
plot_with_margin <- random_effect_plots + patchwork::plot_spacer() + patchwork::wrap_elements(p) +
  plot_layout(ncol = 1, heights = c(1, 0.1, 3))  # Adjust heights for margin
plot_with_margin

cooks_distance <- cooks.distance(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)

# Calculate standardized residuals (Pearson)
std_residuals <- residuals(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, type = "pearson")

# Calculate Hat values (Leverage)
hat_values <- hatvalues(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)

# Calculate Bonferroni p-values for residuals (using an approximate approach)
# This assumes that the number of observations is available
n <- nrow(all_datasets_v2) # Number of observations
bonferroni_p <- 2 * (1 - pnorm(abs(std_residuals)))  # Two-tailed p-values
bonferroni_p <- p.adjust(bonferroni_p, method = "bonferroni")  # Bonferroni correction

# Combine all metrics into a data frame
diagnostic_data <- data.frame(
  index = 1:n,
  Patient.num = all_datasets_v2$Patient.num,
  cooks_distance = cooks_distance,
  std_residuals = std_residuals,
  hat_values = hat_values,
  bonferroni_p = bonferroni_p
)

# Create a grid of plots with ggplot2
library(gridExtra)

cooks_threshold <- 1  # Influential if Cook's distance > 1
residuals_threshold <- 3  # Influential if standardized residuals > 3 or < -3
bonferroni_threshold <- 0.05  # Influential if Bonferroni p-value < 0.05
leverage_threshold <- 2 * mean(diagnostic_data$hat_values)  # Influential if leverage > 2 * average leverage

# Create a ggplot for each diagnostic plot with threshold lines

# Cook's Distance Plot with threshold line
p1 <- ggplot(diagnostic_data, aes(x = index, y = cooks_distance)) +
  geom_point() +
  geom_hline(yintercept = cooks_threshold, linetype = "dashed", color = "red") + 
  labs(title = "Cook's Distance", x = "Index", y = "Cook's Distance") +
  theme_minimal()

# Standardized Residuals Plot with threshold line
p2 <- ggplot(diagnostic_data, aes(x = index, y = std_residuals)) +
  geom_point() +
  geom_hline(yintercept = residuals_threshold, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -residuals_threshold, linetype = "dashed", color = "red") +
  labs(title = "Standardized Residuals", x = "Index", y = "Standardized Residuals") +
  theme_minimal()

# Bonferroni p-values Plot with threshold line
p3 <- ggplot(diagnostic_data, aes(x = index, y = bonferroni_p)) +
  geom_point() +
  geom_hline(yintercept = bonferroni_threshold, linetype = "dashed", color = "red") +
  labs(title = "Bonferroni p-values", x = "Index", y = "Bonferroni p-value") +
  theme_minimal()

# Hat-values Plot with threshold line
p4 <- ggplot(diagnostic_data, aes(x = index, y = hat_values)) +
  geom_point() +
  geom_hline(yintercept = leverage_threshold, linetype = "dashed", color = "red") +
  labs(title = "Hat-values (Leverage)", x = "Index", y = "Hat-values") +
  theme_minimal()

# Combine all plots in a grid layout
grid.arrange(p1, p2, p3, p4, nrow = 2)



ggplot(all_datasets_v2, aes(sample = residuals(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept))) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "Normality of Residuals") +
  theme_minimal()

ggplot(all_datasets_v2, aes(sample = ranef(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)$Patient.num)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "Normality of Random Effects (Patient.num)") +
  theme_minimal()

model_diagnostics(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)

library(redres)

plot_ranef(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)

dotplot(ranef(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, condVar=TRUE))
plot_redres(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, type = "raw_mar")

library(DHARMa)
res <- simulateResiduals(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, plot = T)
par(mfrow = c(1,2))
plotResiduals(res, all_datasets_v2$time_until_sample_usage)


# Extract scaled residuals
scaled_residuals <- resid(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, scaled = TRUE)

grouping_variable <- model.frame(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)$drug  # E.g., 'subject', 'site', etc.

x <- model.frame(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)$time_until_sample_usage

plot_data <- data.frame(scaled_residuals, x, grouping_variable)

ggplot(plot_data, aes(x = x, y = scaled_residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ grouping_variable) +
  labs(x = "X variable", y = "Scaled Residuals", title = "Residuals by Random Effect")


grouping_variable <- model.frame(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)$Patient.num # E.g., 'subject', 'site', etc.

x <- model.frame(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)$time_until_sample_usage

plot_data <- data.frame(scaled_residuals, x, grouping_variable)

ggplot(plot_data, aes(x = x, y = scaled_residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ grouping_variable) +
  labs(x = "X variable", y = "Scaled Residuals", title = "Residuals by Random Effect")


plot(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)

tidy(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept) %>%
  filter(term == "sd__time_until_sample_usageHandled within 2-72 hours")

coef(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)$drug
visualize(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
ranova(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
r2(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
AIC(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
ci_sig = confint(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, parm="theta_", oldNames = FALSE)
ci_sig2 = ci_sig^2
row.names(ci_sig2) = sapply(row.names(ci_sig2),function(x){gsub("sd_","var_",x)},USE.NAMES=FALSE)
row.names(ci_sig2)[length(row.names(ci_sig2))] = "sigma2"
ci_sig2

all_datasets_v2$.resid <- residuals(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
all_datasets_v2$.fitted <- fitted(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)

ggplot(all_datasets_v2 , aes(x = .fitted, y = .resid, color=drug)) +
  geom_point() +
  geom_smooth() +
  labs(title = "Residuals Level-1 - Between group homoscdasticity check") +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(color = "none") 

ggplot(all_datasets_v2 , aes(x = drug, y = .resid)) +
  geom_boxplot() +
  labs(title = "Residuals Level-1 for each drug - Within-group homoscdasticity check") +
  theme(axis.text.x = element_text(angle = 90))

ggplot(all_datasets_v2 , aes(x = time_until_sample_usage, y = .resid)) +
  geom_boxplot() +
  labs(title = "Residuals Level-1 residuals vs independent variable") +
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_wrap(~drug)

# Test for homoscedasticity
testResiduals(sim_res)
r2(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)
AIC(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)



####----Medium----
mixed_model_medium_box_cox_scaled_re_intercept <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1|Patient.num) + (medium|drug),all_response_metrics)
icc(mixed_model_medium_box_cox_scaled_re_intercept)
summary(mixed_model_medium_box_cox_scaled_re_intercept)
check_model(mixed_model_medium_box_cox_scaled_re_intercept)
model_diagnostics(mixed_model_medium_box_cox_scaled_re_intercept, "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/DSS2/medium/") 

####----Cells----
mixed_model_cells_box_cox_scaled_re_intercept <- lmer(DSS2_boxcox_sclaed2 ~ cells +(1|Patient_ID) + (1 + cells|drug),all_response_metrics)
check_model(mixed_model_cells_box_cox_scaled_re_intercept)
model_diagnostics(mixed_model_cells_box_cox_scaled_re_intercept, "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/DSS2/cells/") 

####----Microenv stimuli----
mixed_model_micro_env_stimuli_box_cox_scaled_re_intercept <- lmer(DSS2_boxcox_sclaed2 ~ microenvironmental_stimuli +(1|Patient_ID) + (1 + microenvironmental_stimuli|drug),all_response_metrics)
sum_micro_s <- summary(mixed_model_micro_env_stimuli_box_cox_scaled_re_intercept)
sum_micro_s$coefficients
model_diagnostics(mixed_model_micro_env_stimuli_box_cox_scaled_re_intercept, "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/DSS2/microenvironmental_stimuli/") 


####----Cell counting method----
unique(all_datasets_v2$cell_counting_method)
mixed_model_cell_counting_method_box_cox_scaled_re_intercept <- lmer(DSS2_boxcox_sclaed2 ~ cell_counting_method + (1|Patient_ID) + (1 + cell_counting_method|drug),all_response_metrics)
summary(mixed_model_cell_counting_method_box_cox_scaled_re_intercept)
model_diagnostics(mixed_model_cell_counting_method_box_cox_scaled_re_intercept, "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/DSS2/cell_counting_method/") 

####----Positive Control + Sensitivity Readout + Nr of Concentration Points
mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_re_intercept <- lmer(DSS2_boxcox_sclaed2 ~ positive_control + (1|Patient_ID) +  (1 + positive_control|drug),all_response_metrics)
summary(mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_re_intercept)
model_diagnostics(mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_re_intercept, "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/DSS2/positive_control/") 

####----Centrifugation Procedure ----
mixed_model_centrifugation_procedure_box_cox_scaled_re_intercept <- lmer(DSS2_boxcox_sclaed2 ~ centrifugation_procedure + (1|Patient_ID) + (1 + centrifugation_procedure|drug),all_response_metrics)
summary(mixed_model_centrifugation_procedure_box_cox_scaled_re_intercept)
check_model(mixed_model_centrifugation_procedure_box_cox_scaled_re_intercept)
model_diagnostics(mixed_model_centrifugation_procedure_box_cox_scaled_re_intercept, "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/DSS2/centrifugation_procedure/") 

####----Plate Reader----
##----here----
mixed_model_plate_reader_box_cox_scaled_re_intercept <- lmer(DSS2_boxcox_sclaed2 ~ plate_reader + (1|Patient_ID) + (1 + plate_reader|drug),all_response_metrics)
summary(mixed_model_plate_reader_box_cox_scaled_re_intercept)
model_diagnostics(mixed_model_plate_reader_box_cox_scaled_re_intercept, "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/DSS2/plate_reader/") 


####----Table of crossed re random intercept----
#df <- df %>% mutate(corrected_p_value = ifelse(df$corrected_p_value == '0.0000', '<0.0001', as.character(df$corrected_p_value)))

re_intercept_models <-c(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, mixed_model_medium_box_cox_scaled_re_intercept, mixed_model_cells_box_cox_scaled_re_intercept, mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_re_intercept, mixed_model_centrifugation_procedure_box_cox_scaled_re_intercept, mixed_model_plate_reader_box_cox_scaled_re_intercept)
df_re_intercept_models <- model_results(re_intercept_models, data_frame_org = all_datasets_v2)
df_re_intercept_models_backup <- df_re_intercept_models
df_re_intercept_models <- df_re_intercept_models_backup
df_re_intercept_models$Random_Effect_Term <- '(1 | Patient.num) + (1 + Experimental var | drug)'
colnames(df_re_intercept_models) <- gsub("_", " ", colnames(df_re_intercept_models))     
colnames(df_re_intercept_models) <- sapply(colnames(df_re_intercept_models), tools::toTitleCase)
df_re_intercept_models <- df_re_intercept_models[order(df_re_intercept_models$`p Value`),]
rownames(df_re_intercept_models) <- NULL

df_re_intercept_models <- df_re_intercept_models %>% mutate(`Corrected p Value` = ifelse(df_re_intercept_models$`p Value`*df_re_intercept_models$Count >= 1, 1, df_re_intercept_models$`p Value`*df_re_intercept_models$Count))

color_tile_custom <- formatter("span", 
                               style = function(x) {
                                 color_scale <- col_numeric(c("lightblue", "white", "lightblue"), domain = range(df_re_intercept_models$`Fixed Effect Coefficient`))
                                 formattable::style(display = "block", padding = "0 4px", `background-color` = color_scale(x))
                               }
)
ft <- formattable(df_re_intercept_models,  list(
  `Fixed Effect Term` = formatter("span", style = ~formattable::style(color = "blue")),
  `Fixed Effect Coefficient` = color_tile_custom,
  `p Value` = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  `Significance Stars` = formatter("span", style = ~formattable::style(color = "darkred"))
))

df_re_intercept_models$`Fixed Effect Term`


# df_re_intercept_models <- df_re_intercept_models %>% mutate(`Experimental Variable`= case_when(`Fixed Effect Term` == "Positive Control: BzCl + Nr of Concentration Points: 5 + Drug Sensitivity Readout: CellTiter-Glo" ~ "Positive Control-Doses-Readout   ", 
#                                                   `Fixed Effect Term` ==  "Time until sample usage: Handled within 2-72 hours" ~ "Time Until Sample Usage > 2h   ",
#                                                   `Fixed Effect Term` ==  "Medium: Mononuclear cell medium" ~ "Medium: MCM    ",
#                                                   `Fixed Effect Term` ==  "Medium: RPMI + fetal bovine serum (FBS) (10%)" ~ "Medium: RPMI + FBS   ",
#                                                   `Fixed Effect Term` ==  "microenvironmental_stimuliNone" ~ "Microenvironmental Stimuli: None    ",
#                                                   `Fixed Effect Term` ==  "Microenvironmental stimuli: Transiently cultured with feeder cells, activate samples with autologous BM T helper cells in the presence of IL-2 and a T-cell expansion cocktai" ~ "Microenvironmental Stimuli: Co-culture with activation   ",
#                                                   `Fixed Effect Term` ==  "Plate readers: Pherastar, Cytation5, Insight, Tecan" ~ "Plate Reader",
#                                                   `Fixed Effect Term` ==  "Plate Reader" ~ "Plate Reader: Pherastar",
#                                                   `Fixed Effect Term` ==  "Centrifugation procedure: LymphoPrepTM gradient centrifugation" ~ "Centrifugation Procedure: 1*    ",
#                                                   `Fixed Effect Term` ==  "Centrifugation procedure: Supernatant isolation at 300g 10min, density centrifugation at 400g for 20min without brake, afterwards always 300g 5 min" ~ "Centrifugation Procedure: 2*   ",
#                                                   `Fixed Effect Term` ==  "Cells: 5000" ~ "Number of Cells per Well: 5000   ",
#                                                   `Fixed Effect Term` ==  "Cell counting method: Trypan blue (Countess™ II FL Automated Cell Counter)" ~ "Cell Counting Method    ",
#                                                   `Fixed Effect Term` ==  "plate_readerVictorX" ~ "Plate Reader: Victor X    ",
#                                                   TRUE ~ `Fixed Effect Term`))
# df_re_intercept_models <- df_re_intercept_models %>% dplyr::mutate(`Experimental Variable`= case_when(`Experimental Variable` ==  "Plate Reader" ~ "Plate Reader: Pherastar",
#                                                              `Experimental Variable` ==  "plate_readerVictorX" ~ "Plate Reader: Victor X    ",
#                                                                                                TRUE ~ `Experimental Variable`))
df_re_intercept_models$`Experimental Variable` <- df_re_intercept_models$`Fixed Effect Term`
#df_re_intercept_models <- df_re_intercept_models %>% mutate(`Corrected p Value` = ifelse(df_re_intercept_models$`Corrected p Value` == '0.0000', '<0.0001', as.character(df_re_intercept_models$`Corrected p Value`)))
df_re_intercept_models$`Reference Group` <- df_re_intercept_models$`Ref Group`
df_re_intercept_models <- df_re_intercept_models %>% mutate(`Reference Condition` = case_when(`Reference Group` ==  "a drug combination of flavopiridol, staurosporine and velcade" ~ "Positive Control: drug combination + Doses: 7 \n+ Readout: CellTiter96     ", 
                                                                                  `Reference Group` ==  "1-2h after receiving" ~ "Time until sample usage <2h",
                                                                                  `Reference Group` ==  "HS-5 conditioned medium" ~ "HS-5 CM",
                                                                                  `Reference Group` ==  "HS-5 CM" ~ "HS-5 CM",
                                                                                  `Reference Group` ==  "Ficoll-Paque centrifugation" ~ "Ficoll-Paque",
                                                                                  `Reference Group` ==  "10000" ~ "10000",
                                                                                  `Reference Group` ==  "Countess" ~ "Countess",
                                                                                  `Reference Group` ==  "Biotek Synergy 2" ~ "Biotek Synergy 2",
                                          TRUE ~ `Reference Group`))

names(df_re_intercept_models) <- make.names(names(df_re_intercept_models), unique = TRUE)
names(df_re_intercept_models) <- gsub("\\.", " ", names(df_re_intercept_models))
colnames(df_re_intercept_models) 
duplicated(names(df_re_intercept_models))
library(forestploter)
#core = list(bg_params = list(fill = c("white"))), changes background to white
tm <- forest_theme(base_size =9,
                   arrow_type = "closed",
                   base_family = "Arial",
                   footnote_gp = gpar(col = "black", cex = 1.0, size = 12, family = "Arial"), 
                   line_size = 0.5, 
                   align = "center", 
                   footnote_parse = F, 
                   xlab_gp = gpar(fontsize = 14, fontfamily = "Arial", cex=1),
                   xaxis_gp = gpar(fontsize = 14, fontfamily = "Arial", cex=1), 
                   xlab_adjust = "center")

df_re_intercept_models$` ` <- paste(rep(" ", 30), collapse = " ")
df_re_intercept_models <- df_re_intercept_models[order(df_re_intercept_models$`Fixed Effect Coefficient`),]
df_re_intercept_models$`p-valueᵃ` <- round(df_re_intercept_models$`Corrected p Value`, 4)
df_re_intercept_models$`p-valueᵃ` <- format(df_re_intercept_models$`p-valueᵃ`, scientific = FALSE, trim = TRUE)
df_re_intercept_models <- df_re_intercept_models %>% mutate(`p-valueᵃ` = ifelse(df_re_intercept_models$`p-valueᵃ` == '0.0000', '<0.0001', as.character(df_re_intercept_models$`p-valueᵃ`)))
df_re_intercept_models$`Experimental Variable` <- df_re_intercept_models$`Fixed Effect Term` 


#df_re_intercept_models$se <- (log(df_re_intercept_models$Upper) - log(df_re_intercept_models$`Fixed Effect Coefficient`))/1.96
p <- forest(df_re_intercept_models[,c('Experimental Variable','Reference Condition', ' ', 'p-valueᵃ')],
            est = df_re_intercept_models$`Fixed Effect Coefficient`,
            lower = df_re_intercept_models$Lower, 
            upper = df_re_intercept_models$Upper,
            sizes = 1,
            ci_column = 3,
            ref_line = 0,
            grid = F,
            #arrow_lab = c("Lower DSS2 than reference group", "Higher DSS2 than refernce group"),
            xlim = c(-1.1, 1.1),
            #ticks_at = c(-1, -0.5, 0, 0.5, 1),
            footnote = "\n\nᵃBonferroni corrected",
            theme = tm, 
            xlab = expression('Change in Scaled DSS'[2]),
            legend_gp = gpar(fontsize = 18, fontfamily = "Arial", cex = 1),
            xaxis_gp = gpar(fontsize = 28, fontfamily = "Arial", cex=1)
            ) 

p
# p1 <- as.ggplot(p) +  # Set your base_size here
#   theme(
#     plot.margin = margin(10, 10, 10, 10),  # Adjust margins if necessary
#     plot.title = element_text(hjust = 0.5),
#     plot.caption = element_text(hjust = 0.5),
#     panel.grid = element_blank(),           # Remove gridlines
#     panel.background = element_blank(),     # Remove panel background
#     axis.line = element_blank(),            # Remove axis lines
#     axis.text = element_blank(),            # Remove axis labels
#     axis.ticks = element_blank(),           # Remove axis ticks
#     plot.background = element_blank(),       # Remove the overall background
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank()
#   )
# footnote_text <- "1: LymphoPrep™ gradient centrifugation\n2: Supernatant isolation, density centrifugation without brake"
# footnote_grob <- textGrob(
#   footnote_text,
#   x = 0.05, hjust = 0, gp = gpar(fontsize = 8, col = "black")
# )
# 
# 
# # Arrange the plot and footnote vertically
# grid.arrange(
#   p,
#   footnote_grob,
#   heights = unit.c(unit(1, "npc") - unit(2, "lines"), unit(2, "lines")),
#   ncol = 1
# )
# 
# 
# add_border(
#   p,
#   row = NULL,
#   col = NULL,
#   part = c("body", "header"),
#   where = c("bottom", "left", "top", "right"),
#   gp = gpar(lwd = 2)
# )
# 
# # Print plot
# plot(p)
# p_wh <- get_wh(p)
# pdf('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/Forest_plot_DSS2.pdf',width = p_wh[1], height = p_wh[2])
# plot(p)
# dev.off()
# 
# 
# dynamic_forestplot <- function(p, width_cm, height_cm, dpi = 1200, filename = "dynamic_forestplot.png") {
#   # Convert cm to inches
#   width_in <- width_cm / 2.54
#   height_in <- height_cm / 2.54
#   
#   # Convert cm to inches (1 inch = 2.54 cm)
#   output_width_in <- output_width_cm / 2.54
#   output_height_in <- output_height_cm / 2.54
#   
#   # Convert inches to pixels for the PNG device
#   output_width_px <- output_width_in * dpi
#   output_height_px <- output_height_in * dpi
#   
#   # Calculate scaling factors for the gtable
#   scale_width <- output_width_cm / 10  # Base width adjustment (10 cm as reference)
#   scale_height <- output_height_cm / 10  # Base height adjustment (10 cm as reference)
#   
#   # Scale the gtable layout
#   scaled_p <- gtable::gtable_filter(p, pattern = ".*", trim = TRUE)  # Keep all grobs
#   scaled_p$widths <- scaled_p$widths * scale_width
#   scaled_p$heights <- scaled_p$heights * scale_height
#   
#   # Save plot
#   png(filename, width = width_in * dpi, height = height_in * dpi)
#   grid.newpage()
#   grid.draw(p)
#   dev.off()
# }
# 
# # Usage
# dynamic_forestplot(p, width_cm = 4, height_cm =1, filename = '/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Forest_plot_DSS2_v2.png')


# Define the physical dimensions (cm) and resolution
output_width_cm <- 16.4#14.4 #9.5
output_height_cm <- 16#14 #10
dpi <- 300  # Resolution in dots per inch

# Convert cm to inches (1 inch = 2.54 cm)
output_width_in <- output_width_cm / 2.54
output_height_in <- output_height_cm / 2.54

# Convert inches to pixels for the PNG device
output_width_px <- output_width_in * dpi
output_height_px <- output_height_in * dpi

# Calculate scaling factors for the gtable
scale_width <- output_width_cm / 10  # Base width adjustment (10 cm as reference)
scale_height <- output_height_cm / 10  # Base height adjustment (10 cm as reference)

p <- edit_plot(p, gp = gpar(cex=1.7, fontfamily="Arial")) #1.4
p <- edit_plot(p, part = "header", gp = gpar(cex=1.7, fontfamily="Arial"))

p <- edit_plot(p, col = 4, part = "header",
               which = "text",
               hjust = unit(0.5, "npc"),
               x = unit(0.5, "npc"))
p <- edit_plot(p, col = 4, part = "body",
          which = "text",
          hjust = unit(0.5, "npc"),
          x = unit(0.5, "npc"))


# Scale the gtable layout
scaled_p <- gtable::gtable_filter(p, pattern = ".*", trim = TRUE)  # Keep all grobs
scaled_p$widths <- scaled_p$widths * scale_width
scaled_p$heights <- scaled_p$heights * scale_height

png('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/v1/Forest_plot_DSS2.png', height = 16, width = 46, unit = "cm", res = 300)
grid.newpage()
grid.draw(scaled_p)
dev.off()

ggsave('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Forest_plot_DSS2_gg.png', plot = p1, height = 20, width = 60, unit = "cm")

df_re_intercept_models <- df_re_intercept_models %>%
  mutate(
    row = factor(`Experimental Variable`, levels = rev(`Experimental Variable`)), # Reverse order for plotting
    size = 1.5, # Size for each point
    p_category = case_when(
      `Corrected P Value` < 0.05 ~ "Below 0.05",
      TRUE ~ "Above 0.05"
    )
  )

g <- ggplot(df_re_intercept_models, aes(x = `Fixed Effect Coefficient`, y = reorder(`Experimental Variable`, `Fixed Effect Coefficient`))) +  # Map color to p_value
  geom_point(shape = 15, size = 3, color = "royalblue") +  # Point for the mean estimate
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, color = "royalblue") +  # Error bars for CI
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Reference line at 0
  geom_text(aes(label = paste("P Value =", `Corrected P Value`), x = 1.1), hjust = 0, size = 3) +
  labs(y = "", x="") +
  scale_x_continuous(limits = c(-1.5, 1.5)) +
  #scale_color_manual(values = c("Below 0.05" = "royalblue",  "Above 0.05" = "black"),name = "p-value") +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 14, face = "bold")
  ) 

##########testing diagnostics####################

mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_re_intercept

ranef_data <- as.data.frame(ranef(mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_re_intercept, condVar = TRUE))

# QQ plot for random effects
ggplot(ranef_data, aes(sample = `condval`)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ grpvar, scales = "free") +  # Separate by grouping variable
  labs(title = "QQ Plot of Random Effects",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

residuals_data <- data.frame(residuals = residuals(mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_re_intercept))

# QQ plot for residuals
ggplot(residuals_data, aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot of Residuals",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()


library(htmlwidgets)

# Save formattable as HTML
saveWidget(as.htmlwidget(ft), "~/Desktop/UiO/Project 1/Figures/DSS2_formattable_table.html", selfcontained = TRUE)

df_top3 <- df_re_intercept_models %>%
  arrange(desc(`Fixed Effect Coefficient`)) %>%  # Arrange by descending Value
  head(3)

ggplot(df_top3, aes(x = `Experimental Variable`, y = `Fixed Effect Coefficient`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "black") +  # Confidence intervals
  theme_minimal() +
  labs(x = "Experimental variable", y = "Difference in Response") + 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  # Remove x-axis text (category names)
    axis.text.y = element_blank(),  # Remove y-axis text (numbers)
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),   # Remove y-axis title
    plot.background = element_blank()
  ) 


mixed_model_time_until_sample_usage_box_cox <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1 |Patient.num),all_datasets_v2)
icc(mixed_model_time_until_sample_usage_box_cox)
visualize(mixed_model_time_until_sample_usage_box_cox, formula = DSS2_boxcox_sclaed2 ~ Patient.num | time_until_sample_usage, plot="model", line.col = "blue", 
          point.col = "red", 
          ci.col = "darkgreen", 
          alpha = 0.9)
summary(mixed_model_time_until_sample_usage_box_cox)
check_model(mixed_model_time_until_sample_usage_box_cox)
dotplot(ranef(mixed_model_time_until_sample_usage_box_cox))

# Extract random effects
random_effects <- ranef(mixed_model_time_until_sample_usage_box_cox_scaled)$drug

# Convert to a data frame for plotting
random_effects_df <- as.data.frame(random_effects)
random_effects_df$drug <- rownames(random_effects_df)

# Create a dotplot using ggplot2
ggplot(random_effects_df, aes(x = `(Intercept)`, y = `time_until_sample_usageHandled within 2-72 hours`, label = drug)) +
  geom_point() +
  geom_text_repel(
    verbose = TRUE,
    max.time = 1,
    max.iter = Inf,
    size = 3
  ) + 
  labs(title = "Dotplot of Random Effects", x = "Random Intercept", y = "Random Slope") +
  theme_minimal()

library(RColorBrewer)
custom_palette <- colorRampPalette(RColorBrewer::brewer.pal(3, "Dark2"))(3)
line_types <- rep("solid", 46)


mixed_model_time_until_sample_usage_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1 + time_until_sample_usage|drug),all_datasets_v2)
visualize(mixed_model_time_until_sample_usage_box_cox_scaled, plot = "model")
model_diagnostics(mixed_model_time_until_sample_usage_box_cox_scaled, output_dir = '~/Desktop/testx/')

icc_model <- icc(mixed_model_time_until_sample_usage_box_cox_scaled)
data.frame(icc = icc_model$icc, design_effect = icc_model$design.effect)

summary(mixed_model_time_until_sample_usage_box_cox_scaled)
check_model(mixed_model_time_until_sample_usage_box_cox_scaled)
dotplot(ranef(mixed_model_time_until_sample_usage_box_cox_scaled))
random_effects_variance <- as.data.frame(VarCorr(mixed_model_time_until_sample_usage_box_cox_scaled))
print(random_effects_variance)

ggplotly(visualize(mixed_model_time_until_sample_usage_box_cox_scaled, plot = "model", line.size = 0.5)
+ scale_size(range = c(3, 5))
+ scale_color_manual(values = custom_palette)
+ guides(shape = "none"))


visualize(mixed_model_time_until_sample_usage_box_cox_scaled, plot="residuals", sample=30)
model_interpretation(all_datasets_v2, "DSS2_boxcox_sclaed", mixed_model_time_until_sample_usage_box_cox_scaled, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_time_until_sample_usage_box_cox_scaled.html", header_1 = 'Time Until Sample Usage - Box Cox sclaed')

simulationOutputL <- simulateResiduals(fittedModel = mixed_model_time_until_sample_usage_box_cox_scaled)
plot(simulationOutputL)

library(sjPlot)
plot_model(mixed_model_time_until_sample_usage_box_cox_scaled, type = "pred", terms = "time_until_sample_usage")
plot_model(mixed_model_time_until_sample_usage_box_cox_scaled, type = "slope", show.ci = TRUE)
plot(allEffects(mixed_model_time_until_sample_usage_box_cox_scaled))



library(HLMdiag)
# Extract residuals and fitted values
residuals_level_1 <- hlm_resid(mixed_model_time_until_sample_usage_box_cox_scaled, level=1)
residuals_level_2 <- hlm_resid(mixed_model_time_until_sample_usage_box_cox_scaled, level="drug")


ggplot(residuals_level_1 , aes(x = .mar.fitted, y = .mar.resid)) +
  geom_point() +
  geom_smooth() + 
  labs(title = "Margial Residuals vs Fitted Values")

qqnorm(residuals_level_1$.mar.resid)
qqline(residuals_level_1$.mar.resid, col = "darkgreen")

ggplot(residuals_level_1 , aes(x = drug, y = .resid)) +
  geom_boxplot() +
  labs(title = "Residuals Level-1 for each drug - Within-group homoscdasticity check") +
  theme(axis.text.x = element_text(angle = 90))

ggplot(residuals_level_1 , aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_smooth() + 
  labs(title = "Residuals Level-1 - Between group homoscdasticity check") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~drug)

ggplot(residuals_level_1 , aes(x = .fitted, y = .resid, color=drug)) +
  geom_point() +
  geom_smooth() + 
  labs(title = "Residuals Level-1 - Between group homoscdasticity check") +
  theme(axis.text.x = element_text(angle = 90))


library(DHARMa)

# Simulate residuals from the model
sim_res <- simulateResiduals(fittedModel = mixed_model_time_until_sample_usage_box_cox_scaled)

# Plot residuals
plot(sim_res)

# Test for homoscedasticity
testResiduals(sim_res)


qqnorm(residuals_level_1$.resid)
qqline(residuals_level_1$.resid, col = "darkgreen")

qqnorm(residuals_level_2$.ls.intercept)
qqline(residuals_level_2$.ls.intercept, col = "darkgreen")

qqnorm(residuals_level_2$.ls.time_until_sample_usage_handled_within_2_72_hours)
qqline(residuals_level_2$.ls.time_until_sample_usage_handled_within_2_72_hours, col = "darkgreen")


ggplot(residuals_level_2 , aes(x = drug, y = .ls.intercept)) +
  geom_boxplot() +
  labs(title = "Residuals vs Fitted Values")



level1_residuals <- residuals(mixed_model_time_until_sample_usage_box_cox_scaled, type = "pearson")
ggplot(all_datasets_v2, aes(x = time_until_sample_usage, y = level1_residuals)) +
  geom_point() +
  facet_wrap(~ drug) +  # Separate plot for each group
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(x = "Independent Variable", y = "Residuals",
       title = "Residuals vs Independent Variable by Group") +
  theme_minimal()


influence_data <- influence(mixed_model_time_until_sample_usage_box_cox_scaled, group = "drug")
plot(influence_data, which = "cook", main = "Influence Plot: Cook's Distance")


library(nlme)
time_nlm <- lme(DSS2_boxcox_sclaed2 ~ time_until_sample_usage, random = ~ 1 + time_until_sample_usage| drug, data=all_datasets_v2)
summary(time_nlm)
tut <- summary(time_nlm)
tabl = tut$tTable 
tabl 
plot(time_nlm)




##----time_until_sample_usage----
mixed_model_time_until_sample_usage <- lmer(DSS2 ~ time_until_sample_usage  + (1 + time_until_sample_usage|drug) , all_datasets_v2)
acf(resid(mixed_model_time_until_sample_usage))
plot(ranef(mixed_model_time_until_sample_usage), main = "Random Effects")
sim_res <- simulateResiduals(fittedModel = mixed_model_time_until_sample_usage)
plot(sim_res)
ggplot(all_datasets_v2, aes(x = drug, y = resid(mixed_model_time_until_sample_usage))) +
  geom_boxplot() +
  labs(title = "Residuals by Group", x = "Group", y = "Residuals") +
  theme_minimal()

hd <- ggplot(all_datasets_v2, aes(x = fitted(mixed_model_time_until_sample_usage), y = resid(mixed_model_time_until_sample_usage), color=drug)) +
  geom_point() +
  labs(title = "Fitted vs Residual", x = "Fitted", y = "Residuals") +
  theme_minimal()
hd + facet_grid( ~ drug)
library(DHARMa)
simulationOutputL <- simulateResiduals(fittedModel = mixed_model_time_until_sample_usage_1)
plot(simulationOutputL)
plotResiduals(simulationOutputL, all_datasets_v2$drug)
testDispersion(simulationOutputL)

visualize(mixed_model_time_until_sample_usage)
summary(mixed_model_time_until_sample_usage)
model_interpretation(all_datasets_v2, "DSS2", mixed_model_time_until_sample_usage, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_time_until_sample_usage.html", header_1 = 'Time Until Sample Usage')
check_model(mixed_model_time_until_sample_usage)
influencePlot(mixed_model_time_until_sample_usage, main = "Residuals vs Leverage")


mixed_model_time_until_sample_usage_1 <- lmer(DSS2 ~ time_until_sample_usage +(1|drug) , all_datasets_v2)
mixed_model_time_until_sample_usage_2 <- lmer(DSS2 ~ time_until_sample_usage  + (1|lab) + (1|drug), all_datasets_v2)
mixed_model_time_until_sample_usage_3 <- lmer(DSS2 ~ time_until_sample_usage  + (1|drug/lab), all_datasets_v2)
mixed_model_time_until_sample_usage_4 <- lmer(DSS2 ~ time_until_sample_usage + (1+lab|drug) , all_datasets_v2)
anova(mixed_model_time_until_sample_usage_1, mixed_model_time_until_sample_usage_2,mixed_model_time_until_sample_usage, mixed_model_time_until_sample_usage_3, mixed_model_time_until_sample_usage_4)
summary(mixed_model_time_until_sample_usage_1)
visualize(mixed_model_time_until_sample_usage_1)
check_model(mixed_model_time_until_sample_usage_1)
AIC(lm(DSS2 ~ time_until_sample_usage + drug + lab, data = all_datasets_v2), mixed_model_time_until_sample_usage_1, mixed_model_time_until_sample_usage_2, mixed_model_time_until_sample_usage_3, mixed_model_time_until_sample_usage)

ggplot(all_datasets_v2, aes(time_until_sample_usage, DSS2, color=drug)) +
  geom_point() +
  geom_line(aes(y=predict(mixed_model_time_until_sample_usage), group=drug, color=drug)) +
  facet_wrap(~drug)


ggplot(all_datasets_v2, aes(x = fitted(mixed_model_time_until_sample_usage), y = resid(mixed_model_time_until_sample_usage), color=drug)) +
  geom_point() +
  labs(title = "Fitted vs Residual", x = "Fitted", y = "Residuals") +
  theme_minimal() +
  facet_wrap(~drug)


mixed_model_time_until_sample_usage_box_cox <- lmer(DSS2_boxcox ~ time_until_sample_usage + (1 + time_until_sample_usage|drug),all_datasets_v2)


##----medium----
mixed_model_medium <- lmer(DSS2 ~ medium + (1 + medium|drug) , all_datasets_v2)
visualize(mixed_model_medium)
summary(mixed_model_medium)
model_interpretation(all_datasets_v2, "DSS2", mixed_model_medium, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_medium.html", header_1 = 'Medium')
check_model(mixed_model_medium)
influencePlot(mixed_model_medium)

mixed_model_medium_1 <- lmer(DSS2 ~ medium +(1|drug) , all_datasets_v2)
mixed_model_medium_2 <- lmer(DSS2 ~ medium  + (1|lab)+(1 |drug) , all_datasets_v2)
mixed_model_medium_3 <- lmer(DSS2 ~ medium  +(1 + lab|drug) , all_datasets_v2)
anova(mixed_model_medium_1, mixed_model_medium_2,mixed_model_medium_3)

###---Transformation box cox ----
mixed_model_medium_box_cox <- lmer(DSS2_boxcox ~ medium + (1 + medium|drug),all_datasets_v2)
summary(mixed_model_medium_box_cox)
check_model(mixed_model_medium_box_cox)

mixed_model_medium_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1 + medium|drug),all_datasets_v2)
summary(mixed_model_medium_box_cox_scaled)
check_model(mixed_model_medium_box_cox_scaled)

model_interpretation(all_datasets_v2, "DSS2_boxcox_sclaed", mixed_model_medium_box_cox_scaled, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_medium_box_cox_scaled.html", header_1 = 'Medium - Box Cox Scaled')


##----cells----
mixed_model_cells <- lmer(DSS2 ~ cells +(1 + cells|drug) , all_datasets_v2)
visualize(mixed_model_cells)
summary(mixed_model_cells)
model_interpretation(all_datasets_v2, "DSS2", mixed_model_cells, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_cells.html", header_1 = 'Nr Cells Used')
check_model(mixed_model_cells)

mixed_model_cells_1 <- lmer(DSS2 ~ cells +(1|drug) , all_datasets_v2)
mixed_model_cells_2 <- lmer(DSS2 ~ cells  + (1|lab)+(1 |drug) , all_datasets_v2)
mixed_model_cells_3 <- lmer(DSS2 ~ cells  +(1 + lab|drug) , all_datasets_v2)
anova(mixed_model_cells_1, mixed_model_cells_2,mixed_model_cells_3)


###----Transformation box cox ----
mixed_model_cells_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ cells + (1 + cells|drug),all_datasets_v2)
summary(mixed_model_cells_box_cox_scaled)
check_model(mixed_model_cells_box_cox_scaled)
visualize(mixed_model_cells_box_cox_scaled, label=NULL, sample=3, plot = 'model')
model_interpretation(all_datasets_v2, "DSS2_boxcox_sclaed", mixed_model_cells_box_cox_scaled, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_cells_box_cox_scaled.html", header_1 = 'Nr Cells Used - Box Cox Scaled')

l <- lm(DSS2_boxcox_sclaed2 ~ cells, all_datasets_v2)
check_model(l)
summary(l)

# Extract residuals and group information
residuals <- resid(mixed_model_cells_box_cox_scaled)
groups <- all_datasets_v2$drug

# Split residuals by group
residuals_by_group <- split(residuals, groups)

shapiro_results <- lapply(names(residuals_by_group), function(df_name) {
  df <- residuals_by_group[[df_name]]
  
  print(df)
  print(df_name)
  shapiro_tests <- shapiro.test(df)
  qq_norm_plot <- qqPlot(df, main = paste("Normal Q-Q plot for ", df_name))
  
  # Combine the results with the data frame name
  list(df_name = df_name, results = shapiro_tests, qqnorm = qq_norm_plot)
})

# Print the result
print(shapiro_results)

qqnomr_plot <- function(x) print(colnames(x))
  qqnorm(x, main = paste("Normal Q-Q plot for ", colnames(x)))
lapply(residuals_by_group, qqnomr_plot)

shapiro_results <- lapply(residuals_by_group, shapiro.test)

# Display results
shapiro_results


ggsave('~/Desktop/tt.jpg', p +theme(legend.position = "none"), height = 20, width = 20)
ggplotly(ggplot(all_datasets_v2, aes(x = fitted(mixed_model_cells_box_cox_scaled), y = resid(mixed_model_cells_box_cox_scaled), color=drug)) +
  geom_point() +
  labs(title = "Fitted vs Residual", x = "Fitted", y = "Residuals") +
  theme_minimal()+ 
  facet_wrap(~drug))

ggplot(all_datasets_v2, aes(x = drug, y = resid(mixed_model_cells_box_cox_scaled), color=drug)) +
  geom_point() +
  labs(title = "Residuals per drug", x = "Drug", y = "Residuals") +
  geom_smooth()+
  theme_minimal() 

ggplotly(ggplot(all_datasets_v2, aes(cells, DSS2, color=drug)) +
  geom_point() + 
  geom_line(aes(y=predict(mixed_model_cells_box_cox_scaled)))) 



library(HLMdiag)
# Extract residuals and fitted values
residuals_level_1 <- hlm_resid(mixed_model_cells_box_cox_scaled, level=1)
residuals_level_2 <- hlm_resid(mixed_model_cells_box_cox_scaled, level="drug")


ggplot(residuals_level_1 , aes(x = .mar.fitted, y = .mar.resid)) +
  geom_point() +
  geom_smooth() + 
  labs(title = "Margial Residuals vs Fitted Values")

qqnorm(residuals_level_1$.mar.resid)
qqline(residuals_level_1$.mar.resid, col = "darkgreen")

ggplot(residuals_level_1 , aes(x = drug, y = .resid)) +
  geom_boxplot() +
  labs(title = "Residuals Level-1 for each drug - Within-group homoscdasticity check") +
  theme(axis.text.x = element_text(angle = 90))

ggplot(residuals_level_1 , aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_smooth() + 
  labs(title = "Residuals Level-1 - Between group homoscdasticity check") +
  theme(axis.text.x = element_text(angle = 90))


qqnorm(residuals_level_1$.resid)
qqline(residuals_level_1$.resid, col = "darkgreen")

qqnorm(residuals_level_2$.ls.intercept)
qqline(residuals_level_2$.ls.intercept, col = "darkgreen")

qqnorm(residuals_level_2$.ls.cells5000)
qqline(residuals_level_2$.ls.cells5000, col = "darkgreen")

ggplot(residuals_level_1 , aes(x = .ls.fitted, y = .ls.resid)) +
  geom_point() +
  geom_smooth() + 
  labs(title = "Residuals vs Fitted Values")

ggplot(residuals_level_2 , aes(x = drug, y = .ls.intercept)) +
  geom_boxplot() +
  labs(title = "Residuals vs Fitted Values")


##----nr_of_concentration_points----
require(flexplot)
mixed_model_conc_points <- lmer(DSS2 ~ nr_of_concentration_points +(1 + nr_of_concentration_points|drug) , all_datasets_v2)
icc(lmer(DSS2 ~ nr_of_concentration_points +(1 + nr_of_concentration_points|drug) , all_datasets_v2))
visualize(mixed_model_conc_points)
summary(mixed_model_conc_points)
model_interpretation(all_datasets_v2, "DSS2", mixed_model_conc_points, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_nr_conc_points.html", header_1 = 'Nr Concentration Points')
check_model(mixed_model_conc_points)

mixed_model_conc_points_1 <- lmer(DSS2 ~ nr_of_concentration_points +(1|drug) , all_datasets_v2)
mixed_model_conc_points_2 <- lmer(DSS2 ~ nr_of_concentration_points  + (1|lab)+(1 |drug) , all_datasets_v2)
mixed_model_conc_points_3 <- lmer(DSS2 ~ nr_of_concentration_points  +(1 + lab|drug) , all_datasets_v2)
anova(mixed_model_conc_points_1, mixed_model_conc_points_2,mixed_model_conc_points_3)

###----Transformation box cox ----
mixed_model_conc_box_cox <- lmer(DSS2_boxcox ~ nr_of_concentration_points + (1 + nr_of_concentration_points|drug),all_datasets_v2)
summary(mixed_model_conc_box_cox)
check_model(mixed_model_conc_box_cox)

mixed_model_conc_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ nr_of_concentration_points + (1 |Patient.num),all_datasets_v2)
summary(mixed_model_conc_box_cox_scaled)
check_model(mixed_model_conc_box_cox_scaled)

model_interpretation(all_datasets_v2, "DSS2_boxcox_sclaed", mixed_model_conc_box_cox_scaled, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_nr_conc_points_box_cox_scaled.html", header_1 = 'Nr Concentration Points - Box Cox Scaled')

ggplot(all_datasets_v2, aes(x = fitted(mixed_model_conc_box_cox_scaled), y = resid(mixed_model_conc_box_cox_scaled), color=drug)) +
  geom_point() +
  labs(title = "Fitted vs Residual", x = "Fitted", y = "Residuals") +
  theme_minimal() 

ggplot(all_datasets_v2, aes(nr_of_concentration_points, DSS2, color=cells)) +
  geom_point() +
  geom_line(aes(y=predict(mixed_model_conc_box_cox_scaled))) +
  facet_wrap(~drug)




##----micro_env_stimuli----
mixed_model_micro_env_stimuli <- lmer(DSS2 ~ microenvironmental_stimuli + (1 + microenvironmental_stimuli|drug) , all_datasets_v2)
visualize(mixed_model_micro_env_stimuli)
summary(mixed_model_micro_env_stimuli)
model_interpretation(all_datasets_v2, "DSS2", mixed_model_micro_env_stimuli, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_microenvironmental_stimuli.html", header_1 = 'Microenvironmental Stimuli')
check_model(mixed_model_micro_env_stimuli)

mixed_model_micro_env_stimuli_1 <- lmer(DSS2 ~ microenvironmental_stimuli +(1|drug) , all_datasets_v2)
mixed_model_micro_env_stimuli_2 <- lmer(DSS2 ~ microenvironmental_stimuli  + (1|lab)+(1 |drug) , all_datasets_v2)
mixed_model_micro_env_stimuli_3 <- lmer(DSS2 ~ microenvironmental_stimuli  +(1 + lab|drug) , all_datasets_v2)
anova(mixed_model_micro_env_stimuli_1, mixed_model_micro_env_stimuli_2,mixed_model_micro_env_stimuli_3)

###----Transformation box cox ----
mixed_model_micro_env_stimuli_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ microenvironmental_stimuli + (1|Patient.num),all_datasets_v2)
summary(mixed_model_micro_env_stimuli_box_cox_scaled)
check_model(mixed_model_micro_env_stimuli_box_cox_scaled)

model_interpretation(all_datasets_v2, "DSS2_boxcox_sclaed", mixed_model_micro_env_stimuli_box_cox_scaled, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_microenvironmental_stimuli_box_cox_scaled.html", header_1 = 'Microenvironmental Stimuli - Box Cox Scaled')


##----cell_counting_method----
mixed_model_cell_counting_method <- lmer(DSS2 ~ cell_counting_method + (cell_counting_method|drug), all_datasets_v2)
visualize(mixed_model_cell_counting_method)
summary(mixed_model_cell_counting_method)
model_interpretation(all_datasets_v2, "DSS2", mixed_model_cell_counting_method, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_cell_counting_method.html", header_1 = 'Cell Counting Method')
check_model(mixed_model_cell_counting_method)

mixed_model_cell_counting_method_1 <- lmer(DSS2 ~ cell_counting_method +(1|drug) , all_datasets_v2)

mixed_model_cell_counting_method_2 <- lmer(DSS2 ~ cell_counting_method  + (1|lab)+(1 |drug) , all_datasets_v2)
mixed_model_cell_counting_method_3 <- lmer(DSS2 ~ cell_counting_method  +(1 + lab|drug) , all_datasets_v2)
anova(mixed_model_cell_counting_method_1, mixed_model_cell_counting_method_2,mixed_model_cell_counting_method_3)


###----Transformation box cox ----
mixed_model_cell_counting_method_box_cox <- lmer(DSS2_boxcox ~ cell_counting_method + (1|Patient.num),all_datasets_v2)
summary(mixed_model_cell_counting_method_box_cox)
check_model(mixed_model_cell_counting_method_box_cox)

mixed_model_cell_counting_method_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ cell_counting_method + (1 + lab|Patient.num),all_datasets_v2)
summary(mixed_model_cell_counting_method_box_cox_scaled)
check_model(mixed_model_cell_counting_method_box_cox_scaled)

model_interpretation(all_datasets_v2, "DSS2_boxcox_sclaed", mixed_model_cell_counting_method_box_cox_scaled, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_cell_counting_method.html", header_1 = 'Cell Counting Method')


##----sensitivity_readout_and_positive_control----
mixed_model_sensitivity_readout_and_positive_control <- lmer(DSS2 ~ positive_control + (1 + lab|drug) , all_datasets_v2)
visualize(mixed_model_sensitivity_readout_and_positive_control)
summary(mixed_model_sensitivity_readout_and_positive_control)
model_interpretation(all_datasets_v2, "DSS2", mixed_model_sensitivity_readout_and_positive_control, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_sensitivity_readout_and_positive_control.html", header_1 = 'Sensitivity Readout and Positive Control')
check_model(mixed_model_sensitivity_readout_and_positive_control)

mixed_model_sensitivity_readout_and_positive_control_1 <- lmer(DSS2 ~ sensitivity_readout + positive_control +(1|drug) , all_datasets_v2)
mixed_model_sensitivity_readout_and_positive_control_2 <- lmer(DSS2 ~ sensitivity_readout + positive_control  + (1|lab)+(1 |drug) , all_datasets_v2)
mixed_model_sensitivity_readout_and_positive_control_3 <- lmer(DSS2 ~ sensitivity_readout + positive_control  +(1 + lab|drug) , all_datasets_v2)
anova(mixed_model_sensitivity_readout_and_positive_control_1, mixed_model_sensitivity_readout_and_positive_control_2,mixed_model_sensitivity_readout_and_positive_control_3)

###----Transformation box cox ----
mixed_model_sensitivity_readout_and_positive_control_box_cox <- lmer(DSS2_boxcox ~ positive_control + (1|Patient.num),all_datasets_v2)
summary(mixed_model_sensitivity_readout_and_positive_control_box_cox)
check_model(mixed_model_sensitivity_readout_and_positive_control_box_cox)

mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ positive_control + (1 + positive_control|drug),all_datasets_v2)
summary(mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled)
check_model(mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled)

model_interpretation(all_datasets_v2, "DSS2_boxcox_sclaed", mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled.html", header_1 = 'Sensitivity Readout and Positive Control - Box Cox Scaled')



##----centrifugation_procedure----
mixed_model_centrifugation_procedure <- lmer(DSS2 ~ centrifugation_procedure + (1 + lab|drug) , all_datasets_v2)
visualize(mixed_model_centrifugation_procedure)
summary(mixed_model_centrifugation_procedure)
model_interpretation(all_datasets_v2, "DSS2", mixed_model_centrifugation_procedure, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_centrifugation_procedure.html", header_1 = 'Centrifugation Procedure')
check_model(mixed_model_centrifugation_procedure)

mixed_model_centrifugation_procedure_1 <- lmer(DSS2 ~ centrifugation_procedure +(1|drug) , all_datasets_v2)
mixed_model_centrifugation_procedure_2 <- lmer(DSS2 ~ centrifugation_procedure  + (1|lab)+(1 |drug) , all_datasets_v2)
mixed_model_centrifugation_procedure_3 <- lmer(DSS2 ~ centrifugation_procedure  +(1 + lab|drug) , all_datasets_v2)
anova(mixed_model_centrifugation_procedure_1, mixed_model_centrifugation_procedure_2,mixed_model_centrifugation_procedure_3)

###----Transformation box cox ----
mixed_model_centrifugation_procedure_box_cox <- lmer(DSS2_boxcox ~ centrifugation_procedure + (1|Patient.num),all_datasets_v2)
summary(mixed_model_centrifugation_procedure_box_cox)
check_model(mixed_model_centrifugation_procedure_box_cox)

mixed_model_centrifugation_procedure_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ centrifugation_procedure + (1|Patient.num),all_datasets_v2)
summary(mixed_model_centrifugation_procedure_box_cox_scaled)
check_model(mixed_model_centrifugation_procedure_box_cox_scaled)

model_interpretation(all_datasets_v2, "DSS2_boxcox_sclaed", mixed_model_centrifugation_procedure_box_cox_scaled, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_centrifugation_procedure_box_cox_scaled.html", header_1 = 'Centrifugation Procedure - Box Cox Scaled')



##----plate_reader----
mixed_model_plate_reader <- lmer(DSS2 ~ plate_reader + (1 + lab|drug) , all_datasets_v2)
visualize(mixed_model_plate_reader)
summary(mixed_model_plate_reader)
model_interpretation(all_datasets_v2, "DSS2", mixed_model_plate_reader, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_plate_reader.html", header_1 = 'Plate Reader')
check_model(mixed_model_plate_reader)

mixed_model_plate_reader_1 <- lmer(DSS2 ~ plate_reader +(1|drug) , all_datasets_v2)
mixed_model_plate_reader_2 <- lmer(DSS2 ~ plate_reader  + (1|lab)+(1 |drug) , all_datasets_v2)
mixed_model_plate_reader_3 <- lmer(DSS2 ~ plate_reader  +(1 + lab|drug) , all_datasets_v2)
anova(mixed_model_plate_reader_1, mixed_model_plate_reader_2,mixed_model_plate_reader_3)

###----Transformation box cox ----
mixed_model_plate_reader_box_cox <- lmer(DSS2_boxcox ~ plate_reader + (1|Patient.num),all_datasets_v2)
summary(mixed_model_plate_reader_box_cox)
check_model(mixed_model_plate_reader_box_cox)

mixed_model_plate_reader_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ plate_reader + (1|Patient.num),all_datasets_v2)
summary(mixed_model_plate_reader_box_cox_scaled)
check_model(mixed_model_plate_reader_box_cox_scaled)

model_interpretation(all_datasets_v2, "DSS2_boxcox_sclaed", mixed_model_plate_reader_box_cox_scaled, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_plate_reader_box_cox_scaled.html", header_1 = 'Plate Reader - Box Cox Scaled')






#----COMBAT ANALYSIS ###############################################################################################################################################################################################
##----Generating combat score ----
#ComBat
subset(all_datasets_v2, lab == 'Beat AML')
specific_order <- c("Helsinki","Beat AML", "Karolinska", "Oslo")

# Convert the column to a factor with the specific order
all_datasets_v2$lab <- factor(all_datasets_v2$lab, levels = specific_order)
all_datasets_for_heatmap <- pivot_wider(subset(all_datasets_v2, select=c(drug, DSS2, lab, Patient.num)), names_from = drug, values_from = DSS2, values_fn = list(DSS2 = mean))
all_datasets_for_heatmap <- as.data.frame(all_datasets_for_heatmap)
rownames(all_datasets_for_heatmap) <- all_datasets_for_heatmap$Patient.num
all_datasets_for_heatmap$Patient.num <- NULL
#all_datasets_for_heatmap$lab <- NULL


# Check how many values are missing in each row
missing_per_row <- apply(all_datasets_for_heatmap, 1, function(row) sum(is.na(row)))
#Remove row if total missing in row is higher than 200
all_datasets_for_heatmap <- all_datasets_for_heatmap[missing_per_row <= 30, ]
# Check how many values are missing in each col
missing_per_col <- apply(all_datasets_for_heatmap, 2, function(row) sum(is.na(row)))
#Remove col if total missing in col is higher than 100
all_datasets_for_heatmap <- all_datasets_for_heatmap[,missing_per_col <= 300]


# Calculate proportion of missing values per column
missing_prop <- colMeans(is.na(all_datasets_for_heatmap))

# Keep columns with <= 30% missing
all_datasets_for_heatmap <- all_datasets_for_heatmap[, missing_prop <= 0.3]

subset(all_datasets_for_heatmap, lab== 'Helsinki')
rows_with_na <- all_datasets_for_heatmap[rowSums(is.na(all_datasets_for_heatmap)) > 0, ]
print(subset(rows_with_na, lab == 'Helsinki'))
all_datasets_for_heatmap <- na.omit(all_datasets_for_heatmap)
subset(all_datasets_for_heatmap, lab== 'Helsinki')

all_datasets_for_heatmap_1 <- as.data.frame(all_datasets_for_heatmap)
#all_datasets_for_heatmap_1$lab <- NULL
all_datasets_for_heatmap_1 <- as.data.frame(all_datasets_for_heatmap_1)

dim(all_datasets_for_heatmap_1)
dim(all_datasets_for_heatmap)


# Perform PCA excluding specific columns
pca_result <- prcomp(all_datasets_for_heatmap_1, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x[, 1:2])  # Using only the first two principal components

# Combine PC scores with color information
pca_data <- cbind(pc_scores, lab = all_datasets_for_heatmap$lab)
print(pca_data)
pca_data <- cbind(pca_data, patient_id = rownames(all_datasets_for_heatmap))

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = lab)) +
  geom_point() +
  #geom_text(nudge_x = 0.2, nudge_y = 0.2, size = 3) +
  labs(color = "Labs", title = "PCA Plot") + 
  theme(legend.text = element_text(size = 12),       # Increase legend label size
                                              legend.title = element_text(size = 15), 
                                              plot.title = element_text(size = 20)) 

ggsave(plot = p, filename = file_loc)

all_datasets_for_heatmap_1 <- na.omit(as.data.frame(all_datasets_for_heatmap_1))



ppca_result <- pca(all_datasets_for_heatmap, method = "ppca", nPcs = 2)

# Extract imputed matrix
imputed_matrix <- completeObs(ppca_result)

batch <- all_datasets_for_heatmap_1$lab
all_datasets_for_heatmap_1_combat <- all_datasets_for_heatmap_1
all_datasets_for_heatmap_1$lab <- NULL
mod <- model.matrix(~1, data = data.frame(batch = batch))  # Null model (no covariates)
all_datasets_for_heatmap_combatch <- ComBat(dat = t(all_datasets_for_heatmap_1), batch = batch, mod = mod, prior.plots = TRUE)
dim(all_datasets_for_heatmap_combatch)

prcomp_before <- prcomp(t(all_datasets_for_heatmap_1), scale. = TRUE)
prcomp_after <- prcomp(t(all_datasets_for_heatmap_combatch), scale. = TRUE)

plot(prcomp_before$x[,1:2], col = batch, main = "Before ComBat")
plot(prcomp_after$x[,1:2], col = batch, main = "After ComBat")

# Perform PCA excluding specific columns
pca_result <- pca(t(all_datasets_for_heatmap_combatch), method = "ppca", nPcs = 2)
#pca_result <- prcomp(t(all_datasets_for_heatmap_combatch), scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x[, 1:2])  # Using only the first two principal components
pc_scores <- scores(pca_result)
# Combine PC scores with color information
pca_data <- cbind(pc_scores, lab = all_datasets_for_heatmap$lab)
print(pca_data)
pca_data <- cbind(pca_data, patient_id = rownames(all_datasets_for_heatmap))

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = lab, label = as.character(patient_id))) +
  geom_point() +
  geom_text(nudge_x = 0.2, nudge_y = 0.2, size = 3) +
  labs(title = "")

ggplot(pca_data, aes(x = PC1, y = PC2, color = lab)) +
  geom_point() +
  #geom_text(nudge_x = 0.2, nudge_y = 0.2, size = 3) +
  labs(color = "Labs", title = "PCA Plot") + 
  theme(legend.text = element_text(size = 12),       # Increase legend label size
        legend.title = element_text(size = 15), 
        plot.title = element_text(size = 20))

#Combining the two datasets
all_datasets_for_heatmap_1$Patient.num <- rownames(all_datasets_for_heatmap_1)
combatch_t <- as.data.frame(t(all_datasets_for_heatmap_combatch))
combatch_t$Patient.num <- rownames(combatch_t)
all_dataset_filtered <- gather(all_datasets_for_heatmap_1, drug, DSS2, 'AT7519':'Lenvatinib', factor_key=TRUE)
combatch_long <- gather(combatch_t, drug, combatch, 'AT7519':'Lenvatinib', factor_key=TRUE)
combatch_all_data <- inner_join(all_dataset_filtered, combatch_long, by = c("drug", "Patient.num"))
combatch_all_datasets <- inner_join(combatch_all_data, all_datasets_v2, by = c("drug", "Patient.num", "DSS2"))
combatch_all_datasets <- subset(combatch_all_datasets, Patient.num != 'AML_009 _frozen')

combatch_drugs <- list(Karolinska = subset(combatch_all_datasets, lab == 'Karolinska'), BeatAML = subset(combatch_all_datasets, lab == 'Beat AML'), Enserink = subset(combatch_all_datasets, lab == 'Enserink'), FIMM = subset(combatch_all_datasets, lab == 'FIMM'))
df_list <- combatch_drugs
heatmaps_list <- list()

density_heatmaps <- create_density_heatmaps(df_list, 'combatch', drug_names_to_include)
d_plots <- create_density_plots(df_list, 'drug', 'combatch', drug_names_to_include, title = "Density plot of selected drugs across different datasets")


g_combat <- grouped_ggbetweenstats( # paired samples
  data = combatch_all_datasets,
  x = lab,
  y = combatch,
  grouping.var = drug,
  type = "nonparametric", # for wilcoxon
  centrality.plotting = FALSE # remove median
)

ggplot(combatch_all_datasets, aes(x = combatch)) +
  geom_density() +
  facet_wrap(~drug)


combatch_all_datasets_heatmap <- pivot_wider(combatch_all_datasets[,c("drug", "Patient.num", "lab", "combatch")], names_from = drug, values_from = "combatch", values_fn = ~mean(.x, na.rm = TRUE))
print(combatch_all_datasets_heatmap)

ppca_metrics <- pca(combatch_all_datasets_heatmap[,-c(1, 2)], method="ppca", nPcs=3, seed=123)
## Get the estimated complete observations
ppca_scores_metrics <- scores(ppca_metrics)
ppca_data_metrics <- merge(ppca_scores_metrics, combatch_all_datasets, by = 'row.names')
#ppca_data_metrics$lab <- gsub('BeatAML', 'Beat AML', ppca_data_metrics$lab)
ppca_data_metrics$lab <- factor(ppca_data_metrics$lab, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
ppca_data_metrics$Lab <- as.factor(ppca_data_metrics$lab)
combat_ppca_plot <- ggplot(as.data.frame(ppca_data_metrics), aes(x = PC1, y = PC2, color=Lab)) +
  geom_point() +
  #geom_text_repel(data = furthest_points, aes(label = Patient.num), size = 4, box.padding = 0.5, max.overlaps = 20) +
  #scale_color_gradient(low = "lightblue", high = "blue") + 
  scale_color_manual(values = c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")) +
  labs(color = "") + 
  guides(
    color = guide_legend(
      override.aes = aes(shape = 15, size = 4, width = 1.5, height = 1),  # Adjust legend symbol size
      label.spacing = unit(0.5, "cm")  # Reduce space between text and legend symbol
    )
  ) +  
  #theme_minimal()+
  theme(
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 12, family = "Arial"),    # Center legend text
    legend.position = "top",                                
    legend.justification = "center",        # Ensure the legend is centered
    legend.box = "horizontal",
    legend.key = element_blank(),                            # Remove key background
    legend.spacing.x = unit(0.2, "cm"),  
    legend.key.size = unit(0.5, "lines"),      # Adjust size of the colored squares
    legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
    plot.margin = margin(10, 10, 10, 10),
    axis.title = element_text(size = 12, family = "Arial"),    # Axis title size
    axis.text = element_text(size = 10, family = "Arial")
  )

ppca_metrics_org_combat_data <- pca(all_datasets_for_heatmap_1_combat[,-c(1)], method="ppca", nPcs=3, seed=123)
## Get the estimated complete observations
ppca_scores_metrics_org_combat_data <- scores(ppca_metrics_org_combat_data)
ppca_scores_metrics_org_combat_data <- merge(ppca_scores_metrics_org_combat_data, all_datasets_for_heatmap_1_combat, by = 'row.names')
#ppca_data_metrics$lab <- gsub('BeatAML', 'Beat AML', ppca_data_metrics$lab)
ppca_scores_metrics_org_combat_data$lab <- factor(ppca_scores_metrics_org_combat_data$lab, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
ppca_scores_metrics_org_combat_data$Lab <- as.factor(ppca_scores_metrics_org_combat_data$lab)

all_datasets_for_heatmap_1_combat

ppca_plot_org_combat_data <- ggplot(as.data.frame(ppca_scores_metrics_org_combat_data), aes(x = PC1, y = PC2, color=Lab)) +
  geom_point() +
  #geom_text_repel(data = furthest_points, aes(label = Patient.num), size = 4, box.padding = 0.5, max.overlaps = 20) +
  #scale_color_gradient(low = "lightblue", high = "blue") + 
  scale_color_manual(values = c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")) +
  labs(color = "") + 
  guides(
    color = guide_legend(
      override.aes = aes(shape = 15, size = 4, width = 1.5, height = 1),  # Adjust legend symbol size
      label.spacing = unit(0.5, "cm")  # Reduce space between text and legend symbol
    )
  ) +  
  #theme_minimal()+
  theme(
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 12, family = "Arial"),    # Center legend text
    legend.position = "top",                                
    legend.justification = "center",        # Ensure the legend is centered
    legend.box = "horizontal",
    legend.key = element_blank(),                            # Remove key background
    legend.spacing.x = unit(0.2, "cm"),  
    legend.key.size = unit(0.5, "lines"),      # Adjust size of the colored squares
    legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
    plot.margin = margin(10, 10, 10, 10),
    axis.title = element_text(size = 12, family = "Arial"),    # Axis title size
    axis.text = element_text(size = 10, family = "Arial")
  )
#guides(color = guide_legend(override.aes = aes(label = "", alpha = 1)))

#print(ppca_plot)
#print(subset(ppca_data_metrics, Lab == "Beat AML"))
heatmap_metrics <- as.data.frame(combatch_all_datasets_heatmap)
print(rownames(heatmap_metrics))

rownames(heatmap_metrics) <- make.names(heatmap_metrics$Patient.num, unique = TRUE)
#heatmap_metrics <- t(as.matrix(heatmap_metrics[,-c(1, 2)]))
dist_no_na <- function(x) {
  #d <- daisy(x, metric = "euclidean")
  d <- as.dist(proxy::dist(x, method = "euclidean", pairwise = TRUE))  # pairwise handles NAs
  d[is.na(d)] <- max(d, na.rm = TRUE)  # Replace NA distances with max value
  as.dist(d)
}

# Perform clustering while ignoring NAs
#row_clustering <- hclust(dist_no_na(heatmap_metrics), method = "ward.D")
#col_clustering <- hclust(dist_no_na(heatmap_metrics), method = "ward.D")

col_anno <- subset(heatmap_metrics, select="lab") %>% dplyr::rename(Study = lab) %>% as.data.frame()
col_anno <- factor(col_anno$Study, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
levels(col_anno)
a_col_h <- ComplexHeatmap::HeatmapAnnotation(Study = col_anno, 
                                             col = list(Study= c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")), show_annotation_name = FALSE, gp = gpar(fonface = "plain"))
print(a_col_h)



# if(m == "IC50"){
#   heatmap_metrics <- t(scale(heatmap_metrics[,-c(1, 2)]))
# }
# else{
#   heatmap_metrics <- t(as.matrix(heatmap_metrics[,-c(1, 2)]))
# }
heatmap_metrics <- t(as.matrix(heatmap_metrics[,-c(1, 2)]))

library(circlize)

h <- ComplexHeatmap::Heatmap(
  heatmap_metrics,
  name = "ComBat",  # Legend header
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "ward.D",
  clustering_method_columns = "ward.D2",
  row_names_gp = gpar(fontsize = 9),  # Row font size 10
  column_names_gp = gpar(fontsize = 1),  # Column font size 10
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_column_names = FALSE,
  top_annotation = a_col_h,  # Column annotations
  left_annotation = NULL,  # Row annotations
  na_col = "grey",
  row_gap = unit(1, "cm"),             # Space between rows
  column_gap = unit(1, "cm"),          # Space between columns
  heatmap_legend_param = list(
    title = "Combat",  # Add header to the legend
    title_gp = gpar(fontsize = 11, fontface = "bold"),  # Legend title font
    labels_gp = gpar(fontsize = 10), 
    grid_height = unit(8, "mm")
  ),
  #annotation_legend_param = list(title = ""), 
  width = unit(8, "cm"),  # 12 #20
  height = unit(11, "cm"),  # 9.6 #16
  column_title = " ",
  column_title_gp = gpar(fontsize = 20), 
  row_dend_width = unit(2, "cm"),
  clustering_distance_rows = "euclidean",
  column_dend_height = unit(2, "cm"),
  clustering_distance_columns = dist_no_na
)

plot(
  h,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE, 
  padding = unit(c(0, 0, 0, 0), "mm")
)
library(patchwork)
heat_g <- grid.grabExpr(plot(
  h,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE, 
  padding = unit(c(0, 0, 0, 0), "mm")
))

combatch_all_datasets_heatmap_avg <- combatch_all_datasets %>% group_by(Patient.num) %>% dplyr::summarize('Average ComBat' = mean(combatch))
ppca_data_metrics <- ppca_data_metrics %>% left_join(combatch_all_datasets_heatmap_avg, by = 'Patient.num')
ppca_data_metrics$`Average ComBat` <- as.numeric(ppca_data_metrics$`Average ComBat`)

combat_ppca_plot <- ggplot(as.data.frame(ppca_data_metrics), aes(x = PC1, y = PC2, color = `Average ComBat`)) +
  # 1st Category (A)
  geom_point(aes(color = ifelse(lab == "Beat AML", `Average ComBat`, NA)), size = 4) +  
  scale_color_gradient2(low = "#ccece6", mid = "#8dd3c7", high = "#145A32", midpoint = median(ppca_data_metrics[ppca_data_metrics$lab == "Beat AML",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Beat AML", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +  
  new_scale_color() +  # Reset color scale for other points
  
  # 2nd Category (B)
  geom_point(aes(color = ifelse(lab == "Oslo", `Average ComBat`, NA)), size = 4) +  
  scale_color_gradient2(low = "#ffe6c7", mid = "#fdb462", high = "#a34700",midpoint = median(ppca_data_metrics[ppca_data_metrics$lab == "Oslo",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Oslo", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +
  new_scale_color() +
  
  # 3rd Category (C)
  geom_point(aes(color = ifelse(lab == "Helsinki", `Average ComBat`, NA)), size = 4) +  
  scale_color_gradient2(low = "#ffd2cc", mid = "#fb8072", high = "#a12b1d", ,midpoint = median(ppca_data_metrics[ppca_data_metrics$lab == "Helsinki",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Helsinki", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +
  new_scale_color() +
  
  # 4th Category (D)
  geom_point(aes(color = ifelse(lab == "Karolinska", `Average ComBat`, NA)), size = 4) +  
  scale_color_gradient2(low = "#d6ecff", mid = "#80b1d3", high = "#23537d", ,midpoint = median(ppca_data_metrics[ppca_data_metrics$lab == "Karolinska",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Karolinska", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +
  #geom_point() +
  #geom_text_repel(data = furthest_points, aes(label = Patient.num), size = 4, box.padding = 0.5, max.overlaps = 20) +
  #scale_color_gradient(low = "blue", high = "red") + 
  #scale_shape_manual(values = c("Beat AML" = 16,"Oslo"= 17,"Helsinki" = 18, "Karolinska" = 19)) +
  labs(color = "") + 
  theme_minimal()+
  theme(
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 12, family = "Arial"),    # Center legend text
    #legend.position = "right",                                
    legend.spacing.x = unit(0.2, "cm"),  
    legend.key.size = unit(0.5, "lines"),      # Adjust size of the colored squares
    legend.margin = margin(t = 0, b = 0, l = 0, r = 0),                            
    plot.margin = margin(10, 10, 10, 10),
    axis.title = element_text(size = 12, family = "Arial"),    # Axis title size
    axis.text = element_text(size = 10, family = "Arial")
  ) 
combat_ppca_plot

grid.arrange(combat_ppca_plot, heat_g, ncol = 2, widths = c(1, 1))
ggsave(paste0('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/combat_heatmap_and_ppca.png'), plot = grid.arrange(ppca_plot, heat_g, ncol = 2, widths = c(1, 1)), width = 6.47 * 2, height = 5.94 * 1, dpi = 300)

ppca_scores_metrics_org_combat_data

heatmap_metrics_org_combat_data <- as.data.frame(ppca_scores_metrics_org_combat_data)
print(rownames(heatmap_metrics_org_combat_data))

rownames(heatmap_metrics_org_combat_data) <- make.names(heatmap_metrics_org_combat_data$Row.names, unique = TRUE)
#heatmap_metrics <- t(as.matrix(heatmap_metrics[,-c(1, 2)]))
dist_no_na <- function(x) {
  #d <- daisy(x, metric = "euclidean")
  d <- as.dist(proxy::dist(x, method = "euclidean", pairwise = TRUE))  # pairwise handles NAs
  d[is.na(d)] <- max(d, na.rm = TRUE)  # Replace NA distances with max value
  as.dist(d)
}

# Perform clustering while ignoring NAs
#row_clustering <- hclust(dist_no_na(heatmap_metrics), method = "ward.D")
#col_clustering <- hclust(dist_no_na(heatmap_metrics), method = "ward.D")

col_anno <- subset(heatmap_metrics_org_combat_data, select="lab") %>% dplyr::rename(Study = lab) %>% as.data.frame()
col_anno <- factor(col_anno$Study, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
levels(col_anno)
a_col_h <- ComplexHeatmap::HeatmapAnnotation(Study = col_anno, 
                                             col = list(Study= c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")), show_annotation_name = FALSE, gp = gpar(fonface = "plain"))
print(a_col_h)



# if(m == "IC50"){
#   heatmap_metrics <- t(scale(heatmap_metrics[,-c(1, 2)]))
# }
# else{
#   heatmap_metrics <- t(as.matrix(heatmap_metrics[,-c(1, 2)]))
# }
heatmap_metrics <- t(as.matrix(heatmap_metrics_org_combat_data[,-c(1:5,45)]))

library(circlize)

h_org_combat_data <- ComplexHeatmap::Heatmap(
  heatmap_metrics,
  name = "ComBat",  # Legend header
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "ward.D",
  clustering_method_columns = "ward.D2",
  row_names_gp = gpar(fontsize = 9),  # Row font size 10
  column_names_gp = gpar(fontsize = 1),  # Column font size 10
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_column_names = FALSE,
  top_annotation = a_col_h,  # Column annotations
  left_annotation = NULL,  # Row annotations
  na_col = "grey",
  row_gap = unit(1, "cm"),             # Space between rows
  column_gap = unit(1, "cm"),          # Space between columns
  heatmap_legend_param = list(
    title = expression("DSS"[2]),  # Add header to the legend
    title_gp = gpar(fontsize = 11, fontface = "bold"),  # Legend title font
    labels_gp = gpar(fontsize = 10), 
    grid_height = unit(8, "mm")
  ),
  #annotation_legend_param = list(title = ""), 
  width = unit(8, "cm"),  # 12 #20
  height = unit(11, "cm"),  # 9.6 #16
  column_title = " ",
  column_title_gp = gpar(fontsize = 20), 
  row_dend_width = unit(2, "cm"),
  clustering_distance_rows = "euclidean",
  column_dend_height = unit(2, "cm"),
  clustering_distance_columns = dist_no_na
)

plot(
  h_org_combat_data,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE, 
  padding = unit(c(0, 0, 0, 0), "mm")
)
library(patchwork)
heat_g_org_combat_data <- grid.grabExpr(plot(
  h_org_combat_data,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE, 
  padding = unit(c(0, 0, 0, 0), "mm")
))

grid.arrange(ppca_plot_org_combat_data, heat_g_org_combat_data, ncol = 2, widths = c(1, 1))
ggsave(paste0('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/org_data_missing_used_for_combat_heatmap_and_ppca.png'), plot = grid.arrange(ppca_plot_org_combat_data, heat_g_org_combat_data, ncol = 2, widths = c(1, 1)), width = 6.47 * 2, height = 5.94 * 1, dpi = 300)


##----All labs combat mixed model ----
combatch_mixed_model <- lmer(combatch ~ lab + (1 |drug) + (1|Patient.num), combatch_all_datasets)
visualize(combatch_mixed_model)
summary(combatch_mixed_model)
#model_interpretation(combatch_all_datasets, "combatch", combatch_mixed_model, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/combatch_all_labs.html", header_1 = 'Combat vs lab - all labs')
check_model(combatch_mixed_model)

combatch_mixed_model_DSS2 <- lmer(combatch ~ DSS2 + lab + (1 + lab |drug) , combatch_all_datasets)
visualize(combatch_mixed_model_DSS2)
summary(combatch_mixed_model_DSS2)
model_interpretation(combatch_all_datasets, "combatch", combatch_mixed_model_DSS2, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/combatch_all_labs_DSS2.html", header_1 = 'Combat vs DSS2 - all labs')
check_model(combatch_mixed_model_DSS2)

combatch_mixed_model_DSS2_1 <- lmer(combatch ~ DSS2 + lab +(1|lab)+ (1 |drug) , combatch_all_datasets)
combatch_mixed_model_DSS2_2 <- lmer(combatch ~ DSS2 + lab + (1 + lab|drug) , combatch_all_datasets)
anova(combatch_mixed_model_DSS2, combatch_mixed_model_DSS2_1, combatch_mixed_model_DSS2_2)

model_formula <- formula(combatch_mixed_model)
dependent_var <- as.character(model_formula[[2]])
fixed_effects <- fixef(combatch_mixed_model)  # Get fixed effects
fixed_effect_terms <- names(fixed_effects)
fixed_effect_names <- fixed_effect_terms
fixed_effect_coef <- fixed_effects[fixed_effect_names]
coef_model <- summary(combatch_mixed_model)$coefficients
p_values <- coef_model[, "Pr(>|t|)"]
p_values_excluding_intercept <- p_values
response_metric <- dependent_var

# Create a temporary data frame for this model's fixed effects
temp_df <- data.frame(
  Fixed_Effect_Term = fixed_effect_names,
  model = paste0(model_formula[2], model_formula[1], model_formula[3]),
  Fixed_Effect_Coefficient = fixed_effect_coef,
  p_value = p_values_excluding_intercept,
  stringsAsFactors = FALSE
)

model_formula_2 <- formula(combatch_mixed_model_DSS2)
dependent_var_2 <- as.character(model_formula[[2]])
fixed_effects_2 <- fixef(combatch_mixed_model_DSS2)  # Get fixed effects
fixed_effect_terms_2 <- names(fixed_effects_2)
fixed_effect_names_2 <- fixed_effect_terms_2
fixed_effect_coef_2 <- fixed_effects_2[fixed_effect_names_2]
coef_model_2 <- summary(combatch_mixed_model_DSS2)$coefficients
p_values_2 <- coef_model_2[, "Pr(>|t|)"]

temp_df_2 <- data.frame(
  Fixed_Effect_Term = fixed_effect_names_2,
  model = paste0(model_formula_2[2], model_formula_2[1], model_formula_2[3]),
  Fixed_Effect_Coefficient = fixed_effect_coef_2,
  p_value = p_values_2,
  stringsAsFactors = FALSE
)

# Step 4: Add the temporary data frame to the main data frame
combat_analysis <- rbind(temp_df, temp_df_2)

formattable(combat_analysis, list(
  Fixed_Effect_Term = formatter("span", style = ~style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~style(color = "darkred"))
))

##----factors----
ggplot(combatch_all_datasets, aes(x=combatch)) +
  geom_density(color="darkblue", fill="lightblue")


#Combat box-cox traansform
combatch_all_datasets$combat_pos <- combatch_all_datasets$combatch + 10.1
MASS::boxcox(combat_pos ~ time_until_sample_usage, 
             data = combatch_all_datasets,
             lambda = seq(-0.25, 2, length.out = 10))

combatch_all_datasets$combatch_boxcox <- ((combatch_all_datasets$combatch)^0.5 - 1) / 0.5

#Standardize DSS scores w.r.t. labs and drugs
for(i in unique(combatch_all_datasets$lab)){
  for(j in unique(combatch_all_datasets$drug)){
    combatch_all_datasets$combatch_boxcox_sclaed1[combatch_all_datasets$lab==i & combatch_all_datasets$drug==j] <-
      scale(combatch_all_datasets$combatch_boxcox[combatch_all_datasets$lab==i & combatch_all_datasets$drug==j])
  }
}

for(j in unique(combatch_all_datasets$drug)){
  combatch_all_datasets$combatch_boxcox_sclaed2[combatch_all_datasets$drug==j] <-
    scale(combatch_all_datasets$combatch_boxcox[combatch_all_datasets$drug==j])
}

ggplot(combatch_all_datasets, aes(x=combatch)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(combatch_all_datasets, aes(x=combatch_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(combatch_all_datasets, aes(x=combatch_boxcox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")



###----Time until sample usage----
mixed_model_time_until_sample_usage_combatch <- lmer(combatch ~ DSS2+ time_until_sample_usage +(1 + time_until_sample_usage|drug) + (1|Patient.num), combatch_all_datasets)
mixed_model_time_until_sample_usage_combatch <- lmer(combatch_boxcox_sclaed2 ~ time_until_sample_usage +(1 + time_until_sample_usage|drug) + (1|Patient.num), combatch_all_datasets)

visualize(mixed_model_time_until_sample_usage_combatch)
summary(mixed_model_time_until_sample_usage_combatch)
#model_interpretation(combatch_all_datasets, "DSS2", mixed_model_time_until_sample_usage_combatch, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/mixed_model_time_until_sample_usage_combatch.html", header_1 = 'Time Until Sample Usage - combat vs DSS2', adj_nr = 9)
check_model(mixed_model_time_until_sample_usage_combatch)

mixed_model_time_until_sample_usage_1 <- lmer(combatch ~ DSS2 + time_until_sample_usage +(1|drug) , combatch_all_datasets)
mixed_model_time_until_sample_usage_2 <- lmer(combatch ~ DSS2 + time_until_sample_usage  + (1|lab)+(1 |drug) , combatch_all_datasets)
mixed_model_time_until_sample_usage_3 <- lmer(combatch ~ DSS2 + time_until_sample_usage  +(1 + lab|drug) , combatch_all_datasets)
anova(mixed_model_time_until_sample_usage_1, mixed_model_time_until_sample_usage_2,mixed_model_time_until_sample_usage_3)

model_formula_t <- formula(mixed_model_time_until_sample_usage_combatch)
dependent_var_t <- as.character(model_formula[[2]])
fixed_effects_t <- fixef(mixed_model_time_until_sample_usage_combatch)  # Get fixed effects
fixed_effect_terms_t <- names(fixed_effects_t)
fixed_effect_names_t <- fixed_effect_terms_t
fixed_effect_coef_t <- fixed_effects_t[fixed_effect_names_t]
coef_model_t <- summary(mixed_model_time_until_sample_usage_combatch)$coefficients
p_values_t <- coef_model_t[, "Pr(>|t|)"]

t_combat_model <- data.frame(
  Fixed_Effect_Term = fixed_effect_names_t,
  model = paste0(model_formula_t[2], model_formula_t[1], model_formula_t[3]),
  Fixed_Effect_Coefficient = fixed_effect_coef_t,
  p_value = p_values_t,
  stringsAsFactors = FALSE
)

formattable(t_combat_model, list(
  Fixed_Effect_Term = formatter("span", style = ~style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~style(color = "darkred"))
))

###----Medium----
mixed_model_medium_combatch <- lmer(combatch ~ DSS2 + medium + (medium|drug), combatch_all_datasets)
mixed_model_medium_combatch <- lmer(combatch_boxcox_sclaed2 ~ medium + (medium|drug) + (1|Patient.num), combatch_all_datasets)
summary(mixed_model_medium_combatch)

###----cells----
mixed_model_cells_combatch <- lmer(combatch ~ DSS2+ cells +(cells|drug), combatch_all_datasets)
mixed_model_cells_combatch <- lmer(combatch_boxcox_sclaed2 ~ cells +(cells|drug) + (1|Patient.num), combatch_all_datasets)
summary(mixed_model_cells_combatch)

###----microenvironmental_stimuli----
mixed_model_microenvironmental_stimuli_combatch <- lmer(combatch ~ DSS2+ microenvironmental_stimuli +(microenvironmental_stimuli|drug), combatch_all_datasets)
mixed_model_microenvironmental_stimuli_combatch <- lmer(combatch_boxcox_sclaed2 ~ microenvironmental_stimuli +(microenvironmental_stimuli|drug) + (1|Patient.num), combatch_all_datasets)
summary(mixed_model_microenvironmental_stimuli_combatch)

###----cell_counting_method----
mixed_model_cell_counting_method_combatch <- lmer(combatch ~ DSS2+ cell_counting_method +(cell_counting_method|drug), combatch_all_datasets)
mixed_model_cell_counting_method_combatch <- lmer(combatch_boxcox_sclaed2 ~ cell_counting_method +(cell_counting_method|drug) + (1|Patient.num), combatch_all_datasets)
summary(mixed_model_cell_counting_method_combatch)
check_model(mixed_model_cell_counting_method_combatch)

###----sensitivity_readout----
mixed_model_sensitivity_readout_combatch <- lmer(combatch ~ DSS2+ sensitivity_readout +(sensitivity_readout|drug), combatch_all_datasets)
mixed_model_sensitivity_readout_combatch <- lmer(combatch_boxcox_sclaed2 ~ positive_control +(positive_control|drug) + (1|Patient.num), combatch_all_datasets)
summary(mixed_model_sensitivity_readout_combatch)

###----centrifugation_procedure----
mixed_model_centrifugation_procedure_combatch <- lmer(combatch ~ DSS2+ centrifugation_procedure +(centrifugation_procedure|drug), combatch_all_datasets)
mixed_model_centrifugation_procedure_combatch <- lmer(combatch_boxcox_sclaed2 ~ centrifugation_procedure +(centrifugation_procedure|drug) + (1|Patient.num), combatch_all_datasets)
summary(mixed_model_centrifugation_procedure_combatch)

###----plate_reader----
mixed_model_plate_reader_combatch <- lmer(combatch ~ DSS2*plate_reader +(plate_reader|drug) + (1|Patient.num), combatch_all_datasets)
mixed_model_plate_reader_combatch <- lmer(combatch_boxcox_sclaed2 ~ plate_reader +(plate_reader|drug) + (1|Patient.num), combatch_all_datasets)
summary(mixed_model_plate_reader_combatch)

#table 
models <- c(mixed_model_time_until_sample_usage_combatch, mixed_model_medium_combatch, mixed_model_cells_combatch, mixed_model_sensitivity_readout_combatch, 
            mixed_model_centrifugation_procedure_combatch, 
            mixed_model_plate_reader_combatch)
combat_models_df <- model_results(models, data_frame_org=combatch_all_datasets)
combat_models_df <- combat_models_df[order(c(combat_models_df$Model,combat_models_df$Fixed_Effect_Coefficient)), ]
rownames(combat_models_df) <- NULL
colnames(combat_models_df) <- gsub("_", " ", colnames(combat_models_df))     
colnames(combat_models_df) <- sapply(colnames(combat_models_df), tools::toTitleCase)
combat_models_df <- combat_models_df[order(c(combat_models_df$Model,combat_models_df$Fixed_Effect_Coefficient)), ]


color_tile_custom <- formatter("span", 
                               style = function(x) {
                                 color_scale <- col_numeric(c("lightblue", "white", "lightblue"), domain = range(df$`Fixed Effect Coefficient`))
                                 formattable::style(display = "block", padding = "0 4px", `background-color` = color_scale(x))
                               }
)

combat_models_df <- combat_models_df[order(combat_models_df$`Corrected p Value`),]
ft <- formattable(combat_models_df,  list(
  `Fixed Effect Term` = formatter("span", style = ~formattable::style(color = "blue")),
  `Fixed Effect Coefficient` = color_tile_custom,
  `p Value` = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  `Significance Stars` = formatter("span", style = ~formattable::style(color = "darkred"))
))
png("~/Desktop/UiO/Project 1/Figures/V3/formattable_output.png", width = 3000, height = 2000)  # Specify dimensions
grid.table(ft)
dev.off()

combat_models_df$`Corrected p Value`
forest_plot_combat <- ggplot(combat_models_df, aes(x = `Fixed Effect Coefficient`, y = reorder(`Fixed Effect Term`, `Fixed Effect Coefficient`), color = `Corrected p Value`)) +  # Map color to p_value
  geom_point(shape = 15, size = 3) +  # Point for the mean estimate
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +  # Error bars for CI
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Reference line at 0
  labs(x = "Coefficient", y = "Model", title = "Forest Plot of Coefficients") +
  scale_color_gradient(low = "royalblue", high = "darkred", name = "p-value") +  # Color gradient
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(plot = forest_plot_combat, filename = "~/Desktop/UiO/Project 1/Figures/V3/forest_plot_combat.pdf", width = 80, height = 40, units = "cm", limitsize = FALSE)



tm <- forest_theme(base_size =9,
                   arrow_type = "closed",
                   base_family = "Arial",
                   footnote_gp = gpar(col = "black", cex = 0.6), 
                   line_size = 0.5, 
                   align = "center", 
                   footnote_parse = F)

combat_models_df$` ` <- paste(rep(" ", 30), collapse = " ")
combat_models_df <- combat_models_df[order(combat_models_df$`Fixed Effect Coefficient`),]
combat_models_df$`Corrected p-value` <- round(combat_models_df$`Corrected p Value`, 4)
combat_models_df$`Corrected p-value` <- format(combat_models_df$`Corrected p-value`, scientific = FALSE, trim = TRUE)
combat_models_df <- combat_models_df %>% mutate(`p-valueᵃ` = ifelse(combat_models_df$`Corrected p-value` == '0.0000', '<0.0001', as.character(combat_models_df$`Corrected p-value`)))
combat_models_df$`p Value` <- round(combat_models_df$`p Value`, 4)
combat_models_df$`p Value` <- format(combat_models_df$`p Value`, scientific = FALSE, trim = TRUE)
combat_models_df <- combat_models_df %>% mutate(`p-value` = ifelse(combat_models_df$`p Value`== '0.0000', '<0.0001', as.character(combat_models_df$`p Value`)))
combat_models_df$`Experimental Variable` <- combat_models_df$`Fixed Effect Term` 
combat_models_df$`Reference Condition` <- combat_models_df$`Ref Group`


#df_re_intercept_models$se <- (log(df_re_intercept_models$Upper) - log(df_re_intercept_models$`Fixed Effect Coefficient`))/1.96
p_combat <- forest(combat_models_df[,c('Experimental Variable','Reference Condition', ' ', 'p-valueᵃ')],
            est = combat_models_df$`Fixed Effect Coefficient`,
            lower = combat_models_df$Lower, 
            upper = combat_models_df$Upper,
            sizes = 1,
            ci_column = 3,
            ref_line = 0,
            grid = F,
            #arrow_lab = c("Lower DSS2 than reference group", "Higher DSS2 than refernce group"),
            xlim = c(-1, 1),
            #ticks_at = c(-1, -0.5, 0, 0.5, 1),
            footnote = "\n\nᵃBonferroni corrected",
            theme = tm, 
            xlab = expression('Change in Combat corrected DSS'[2]),
            font.label = list(size = 10, family = "Arial"),
            font.ticks = list(size = 8),
            txt_gp = list(label = gpar(fontsize = 10),  # this controls all text
                          xlab = gpar(fontsize = 10),
                          ticks = gpar(fontsize = 10))
) 

# Define the physical dimensions (cm) and resolution
output_width_cm <- 14.4 #9.5
output_height_cm <- 14 #10
dpi <- 300  # Resolution in dots per inch

# Convert cm to inches (1 inch = 2.54 cm)
output_width_in <- output_width_cm / 2.54
output_height_in <- output_height_cm / 2.54

# Convert inches to pixels for the PNG device
output_width_px <- output_width_in * dpi
output_height_px <- output_height_in * dpi

# Calculate scaling factors for the gtable
scale_width <- output_width_cm / 10  # Base width adjustment (10 cm as reference)
scale_height <- output_height_cm / 10  # Base height adjustment (10 cm as reference)

p_combat <- edit_plot(p_combat, gp = gpar(cex=1.5, fontfamily="Arial")) #1.4
p_combat <- edit_plot(p_combat, part = "header", gp = gpar(cex=1.5, fontfamily="Arial"))

p_combat <- edit_plot(p_combat, col = 4, part = "header",
               which = "text",
               hjust = unit(0.5, "npc"),
               x = unit(0.5, "npc"))
p_combat <- edit_plot(p_combat, col = 4, part = "body",
               which = "text",
               hjust = unit(0.5, "npc"),
               x = unit(0.5, "npc"))

# Scale the gtable layout
scaled_p_combat <- gtable::gtable_filter(p_combat, pattern = ".*", trim = TRUE)  # Keep all grobs
scaled_p_combat$widths <- scaled_p_combat$widths * scale_width
scaled_p_combat$heights <- scaled_p_combat$heights * scale_height

png('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/V1/Forest_plot_combat.png', height = 16, width = 46, unit = "cm", res = 300)
grid.newpage()
grid.draw(scaled_p_combat)
dev.off()




#----AUC ANALYSIS ###############################################################################################################################################################################################
##---- density plots ----
density_heatmaps_auc <- create_density_heatmaps(all_response_metrics_df_list, 'AUC', drug_names_to_include)
d_plots_auc <- create_density_plots(all_response_metrics_df_list, 'drug', 'AUC', drug_names_to_include, title = "Density plot of selected drugs across different datasets")

##---- Box Cox Scaling ----
all_response_metrics$AUC_pos <- all_response_metrics$AUC + 0.01
MASS::boxcox(AUC_pos ~ time_until_sample_usage, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$AUC_boxcox <- ((all_response_metrics$AUC)^0.5 - 1) / 0.5

#Standardize DSS scores w.r.t. labs and drugs
for(i in unique(all_response_metrics$lab)){
  for(j in unique(all_response_metrics$drug)){
    all_response_metrics$AUC_boxcox_sclaed1[all_response_metrics$lab==i & all_response_metrics$drug==j] <-
      scale(all_response_metrics$AUC_boxcox[all_response_metrics$lab==i & all_response_metrics$drug==j])
  }
}

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$AUC_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$AUC_boxcox[all_response_metrics$drug==j])
}

ggplot(all_response_metrics, aes(x=AUC)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=AUC_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=AUC_boxcox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")

##---- Factors ----
###----Time Until Sample Usage----
mixed_model_time_until_sample_usage_auc <- lmer(AUC_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient.num) + (time_until_sample_usage|drug) , all_response_metrics)
summary(mixed_model_time_until_sample_usage_auc)
check_model(mixed_model_time_until_sample_usage_auc)
#model_interpretation(all_response_metrics, "AUC", mixed_model_time_until_sample_usage_auc, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/AUC/mixed_model_time_until_sample_usage_auc.html", header_1 = 'Time Until Sample Usage - AUC', adj_nr = 9)

###----medium----
mixed_model_medium_auc <- lmer(AUC_boxcox_sclaed2 ~ medium + (1|Patient.num) + (medium|drug) , all_response_metrics)
summary(mixed_model_medium_auc)
#model_interpretation(all_response_metrics, "AUC", mixed_model_medium_auc, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/AUC/mixed_model_medium_auc.html", header_1 = 'Medium - AUC')

###----cells----
mixed_model_cells_auc <- lmer(AUC_boxcox_sclaed2 ~ cells + (1|Patient.num) + (cells|drug) , all_response_metrics)
summary(mixed_model_cells_auc)
#model_interpretation(all_response_metrics, "AUC", mixed_model_cells_auc, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/AUC/mixed_model_cells_auc.html", header_1 = 'Nr Cells Used - AUC')

###----micro_env_stimuli----
mixed_model_micro_env_stimuli_auc <- lmer(AUC_boxcox_sclaed2 ~ microenvironmental_stimuli + (1|Patient.num)+ (1 + microenvironmental_stimuli|drug) , all_response_metrics)
summary(mixed_model_micro_env_stimuli_auc)
#model_interpretation(all_response_metrics, "AUC", mixed_model_micro_env_stimuli_auc, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/AUC/mixed_model_microenvironmental_stimuli_auc.html", header_1 = 'Microenvironmental Stimuli - AUC')

###----cell_counting_method----
mixed_model_cell_counting_method_auc <- lmer(AUC_boxcox_sclaed2 ~ cell_counting_method + (1|Patient.num) + (cell_counting_method|drug) , all_response_metrics)
summary(mixed_model_cell_counting_method_auc)
#model_interpretation(all_response_metrics, "AUC", mixed_model_cell_counting_method_auc, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/AUC/mixed_model_cell_counting_method.html", header_1 = 'Cell Counting Method - AUC')

###----sensitivity_readout_and_positive_control----
mixed_model_sensitivity_readout_and_positive_control_auc <- lmer(AUC_boxcox_sclaed2 ~ positive_control + (1|Patient.num) + (positive_control|drug) , all_response_metrics)
summary(mixed_model_sensitivity_readout_and_positive_control_auc)
#model_interpretation(all_response_metrics, "AUC", mixed_model_sensitivity_readout_and_positive_control_auc, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/AUC/mixed_model_sensitivity_readout_and_positive_control.html", header_1 = 'Sensitivity Readout and Positive Control - AUC')

###----centrifugation_procedure----
mixed_model_centrifugation_procedure_auc <- lmer(AUC_boxcox_sclaed2 ~ centrifugation_procedure + (1|Patient.num) + (centrifugation_procedure|drug) , all_response_metrics)
summary(mixed_model_centrifugation_procedure_auc)
#model_interpretation(all_response_metrics, "AUC", mixed_model_centrifugation_procedure_auc, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/AUC/mixed_model_centrifugation_procedure_auc.html", header_1 = 'Centrifugation Procedure - AUC')


###----plate_reader----
mixed_model_plate_reader_auc <- lmer(AUC_boxcox_sclaed2 ~ plate_reader + (1|Patient.num) + (plate_reader|drug) , all_response_metrics)
summary(mixed_model_plate_reader_auc)
#model_interpretation(all_response_metrics, "AUC", mixed_model_plate_reader_auc, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/AUC/mixed_model_plate_reader_auc.html", header_1 = 'Plate Reader - AUC')


###----AUC Table----
auc_models <-c(mixed_model_time_until_sample_usage_auc, mixed_model_medium_auc, mixed_model_cells_auc, mixed_model_micro_env_stimuli_auc, mixed_model_cell_counting_method_auc, mixed_model_sensitivity_readout_and_positive_control_auc, mixed_model_centrifugation_procedure_auc, mixed_model_plate_reader_auc)
df_auc_models <- model_results(auc_models)
df_auc_models$Random_Effect_Term <- '(1|Patient.num) + (1 | Experimental var | drug)'
rownames(df_auc_models) <- NULL
formattable(df_auc_models,  list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))



#----IC50 ANALYSIS ###############################################################################################################################################################################################
##---- density plots ----
density_heatmaps_ic50 <- create_density_heatmaps(all_response_metrics_df_list, 'IC50', drug_names_to_include)
d_plots_ic50 <- create_density_plots(all_response_metrics_df_list, 'drug', 'IC50', drug_names_to_include, title = "Density plot of selected drugs across different datasets")

#-log10 Transform IC50
df_transformed <- all_response_metrics %>%
  mutate(log_ic50 = log10(IC50)) %>%
  as.data.frame()

df_transformed <- subset(df_transformed, !is.na(log_ic50) & is.finite(log_ic50))

##---- Box Cox Scaling ----
min(df_transformed$log_ic50)
all_response_metrics$IC50_pos <- all_response_metrics$IC50 + 0.01
MASS::boxcox(IC50_pos ~ time_until_sample_usage, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))


all_response_metrics <- all_response_metrics %>% mutate(IC50_boxcox = ifelse(-log10(all_response_metrics$IC50) == 0, 0, -log10(all_response_metrics$IC50)))
epsilon <- 1e-6
library(MASS)

# Box-Cox transformation requires all values > 0
lambda <- boxcox(lm(IC50_pos ~ 1, data = all_response_metrics), lambda = seq(-2, 2, 0.1))$x[which.max(boxcox(lm(IC50_pos ~ 1, data = all_response_metrics), lambda = seq(-2, 2, 0.1))$y)]

#all_response_metrics <- all_response_metrics %>%
  #mutate(IC50_boxcox = -(IC50_pos^lambda - 1) / lambda)
#all_response_metrics <- all_response_metrics %>%
  #mutate(IC50_boxcox = -(IC50_pos^0.5 - 1) / 0.5)
#all_response_metrics$IC50_boxcox <- -log10(all_response_metrics$IC50_pos/(1+all_response_metrics$IC50_pos))
#all_response_metrics$IC50_boxcox <- BoxCox(all_response_metrics$IC50_pos, lambda = "auto")
#all_response_metrics$IC50_boxcox <- sqrt(all_response_metrics$IC50_pos)
#all_response_metrics$IC50_boxcox <- all_response_metrics$IC50_pos^(1/3)
#all_response_metrics$IC50_boxcox <- -((all_response_metrics$IC50)^0.5 - 1) / 0.5

#Standardize DSS scores w.r.t. labs and drugs
for(i in unique(all_response_metrics$lab)){
  for(j in unique(all_response_metrics$drug)){
    all_response_metrics$IC50_boxcox_sclaed1[all_response_metrics$lab==i & all_response_metrics$drug==j] <-
      scale(all_response_metrics$IC50_boxcox[all_response_metrics$lab==i & all_response_metrics$drug==j])
  }
}

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$IC50_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$IC50_boxcox[all_response_metrics$drug==j])
}
unique(subset(all_response_metrics, is.na(IC50_boxcox_sclaed2), select = c('Patient_ID', 'IC50_boxcox_sclaed2')))
#summary(all_response_metrics$IC50)  # Check for negative values, zero, or NAs
#class(all_response_metrics$IC50) 

hist(all_response_metrics$IC50_boxcox_sclaed2, breaks = 30, main = "Histogram of Logit-Transformed IC50")
shapiro.test(all_response_metrics$IC50_boxcox_sclaed2) 

ggplot(df_transformed, aes(x=log_ic50)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=IC50)) +
  geom_density(color="darkblue", fill="lightblue") + 
  theme_minimal() + 
  ylab("Density") + 
  xlab("IC50")

ggplot(all_response_metrics, aes(x=IC50_boxcox)) +
  geom_density(color="darkblue", fill="lightblue") + theme_minimal() + 
  ylab("Density") + 
  xlab("BoxCox transformed IC50")

ggplot(all_response_metrics, aes(x=IC50_boxcox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue") +
  theme_minimal() + 
  ylab("Density") + 
  xlab("BoxCox transformed z-Scaled IC50")

skewness(all_response_metrics$IC50_boxcox_sclaed2)
kurtosis(all_response_metrics$IC50_boxcox_sclaed2)
sd(all_response_metrics$IC50_boxcox_sclaed2)


ic50_g <- grouped_ggbetweenstats( # paired samples
  data = all_response_metrics,
  x = lab,
  y = IC50_boxcox_sclaed2,
  grouping.var = drug,
  type = "nonparametric", # for wilcoxon
  centrality.plotting = FALSE, # remove median
  ggstatsplot.layer = TRUE,
  xlab = "",        
  ylab = expression("DSS"[2]),
  results.subtitle = TRUE,  # Removes default subtitle
  #subtitle = "{test} (p = {p})", 
  p.adjust.method = "fdr",
  pairwise.display = "significant", #"significant"
  #pairwise.comparisons = TRUE,
  plotgrid.args = list(ncol = 12),
  facet_wrap.args = list(scales = "fixed", strip.position = "top"),
  #plot.margin = margin(0, 10, 10, 10),
  #point.args = list(alpha = 0.9),
  violin.args = list(alpha = 0),
  boxplot.args = list(alpha = 0),
  ggsignif.args = list(textsize = 3, tip_length = 0.005, step_increase = 0.05),
  p.value.label.args = list(
    parse = TRUE
  ),
  ggplot.component = list(
    coord_cartesian(ylim = c(0, 50)),
    #scale_y_continuous(labels = c(0:5, "")),
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))),
    scale_color_manual(values = custom_colors), 
    #ggplot2::scale_y_continuous(limits = c(0, 50), sec.axis = ggplot2::dup_axis(name = NULL)), 
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1))), 
    theme(
      panel.background = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 0)),
      axis.text.y.right = ggplot2::element_blank(), 
      axis.ticks.y.right = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = element_text(family = "Arial", face = "plain",size = 13, color = "black", angle = 25, hjust = 0.6),
      axis.title.y = element_text(family = "Arial", face = "plain", size = 13, color = "black"),
      plot.title = element_text(family = "Arial", face = "plain", size = 13, color = "black", hjust = 0.5, margin = margin(b = 0))))
) 
ggsave("Desktop/UiO/Project 1/Figures/draw/Difference_IC50.png", plot = ic50_g, width = 80, height = 30, units = "cm", limitsize = FALSE)


##---- Factors ----
###----Time Until Sample Usage----
mixed_model_time_until_sample_usage_ic50 <- lmer(IC50_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient.num) + (1 + time_until_sample_usage|drug) , data = all_response_metrics)
summary(mixed_model_time_until_sample_usage_ic50)
check_model(mixed_model_time_until_sample_usage_ic50)
#model_interpretation(all_response_metrics, "IC50", mixed_model_time_until_sample_usage_ic50, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/IC50/mixed_model_time_until_sample_usage_ic50.html", header_1 = 'Time Until Sample Usage - IC50', adj_nr = 9)

#mixed_model_time_until_sample_usage_logic50 <- lmer(log_ic50 ~ time_until_sample_usage + (1 +time_until_sample_usage|drug) , df_transformed)
#summary(mixed_model_time_until_sample_usage_logic50)
#check_model(mixed_model_time_until_sample_usage_logic50)

###----medium----
unique(all_response_metrics[all_response_metrics$lab == 'Helsinki',"medium"])
mixed_model_medium_ic50 <- lmer(IC50_boxcox_sclaed2 ~ medium + (1|Patient.num) + (medium|drug) , all_response_metrics)
summary(mixed_model_medium_ic50)
check_model(mixed_model_medium_ic50)
#model_interpretation(all_response_metrics, "IC50", mixed_model_medium_ic50, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/IC50/mixed_model_medium_ic50.html", header_1 = 'Medium - IC50')

#mixed_model_medium_logic50 <- lmer(log_ic50 ~ medium + (1 + medium|drug) , df_transformed)
#summary(mixed_model_medium_logic50)

###----cells----
mixed_model_cells_ic50 <- lmer(IC50_boxcox_sclaed2 ~ cells + (1|Patient_ID) + (cells|drug) , all_response_metrics)
summary(mixed_model_cells_ic50)
check_model(mixed_model_cells_ic50)
#model_interpretation(all_response_metrics, "IC50", mixed_model_cells_ic50, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/IC50/mixed_model_cells_ic50.html", header_1 = 'Nr Cells Used - IC50')

#mixed_model_cells_logic50 <- lmer(log_ic50 ~ cells +(1 + cells|drug) , df_transformed)
#summary(mixed_model_cells_logic50)

###----micro_env_stimuli----
mixed_model_micro_env_stimuli_ic50 <- lmer(IC50_boxcox_sclaed2 ~ microenvironmental_stimuli + (1|Patient_ID) + (microenvironmental_stimuli|drug) , all_response_metrics)
summary(mixed_model_micro_env_stimuli_ic50)
check_model(mixed_model_micro_env_stimuli_ic50)
#model_interpretation(all_response_metrics, "IC50", mixed_model_micro_env_stimuli_ic50, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/IC50/mixed_model_microenvironmental_stimuli_ic50.html", header_1 = 'Microenvironmental Stimuli - IC50')

#mixed_model_micro_env_stimuli_logic50 <- lmer(log_ic50 ~ microenvironmental_stimuli + (1 + microenvironmental_stimuli|drug) , df_transformed)
#summary(mixed_model_micro_env_stimuli_logic50)

###----cell_counting_method----
mixed_model_cell_counting_method_ic50 <- lmer(IC50_boxcox_sclaed2 ~ cell_counting_method + (1|Patient.num) + (1 + cell_counting_method|drug) , all_response_metrics)
summary(mixed_model_cell_counting_method_ic50)
check_model(mixed_model_cell_counting_method_ic50)
#model_interpretation(all_response_metrics, "IC50", mixed_model_cell_counting_method_ic50, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/IC50/mixed_model_cell_counting_method_ic50.html", header_1 = 'Cell Counting Method - IC50')

#mixed_model_cell_counting_method_logic50 <- lmer(log_ic50 ~ cell_counting_method + (1 + cell_counting_method|drug) , df_transformed)
#summary(mixed_model_cell_counting_method_logic50)

###----sensitivity_readout_and_positive_control----
mixed_model_sensitivity_readout_and_positive_control_ic50 <- lmer(IC50_boxcox_sclaed2 ~ positive_control + (1|Patient_ID) + (positive_control|drug) , all_response_metrics)
summary(mixed_model_sensitivity_readout_and_positive_control_ic50)
check_model(mixed_model_sensitivity_readout_and_positive_control_ic50)
#model_interpretation(all_response_metrics, "IC50", mixed_model_sensitivity_readout_and_positive_control_ic50, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/IC50/mixed_model_sensitivity_readout_and_positive_control_ic50.html", header_1 = 'Sensitivity Readout and Positive Control - IC50')

#mixed_model_sensitivity_readout_and_positive_control_logic50 <- lmer(log_ic50 ~ positive_control + (1 + positive_control|drug) , df_transformed)
#summary(mixed_model_sensitivity_readout_and_positive_control_logic50)

###----centrifugation_procedure----
mixed_model_centrifugation_procedure_ic50 <- lmer(IC50_boxcox_sclaed2 ~ centrifugation_procedure + (1|Patient_ID) + (centrifugation_procedure|drug) , all_response_metrics)
summary(mixed_model_centrifugation_procedure_ic50)
check_model(mixed_model_centrifugation_procedure_ic50)
#model_interpretation(all_response_metrics, "IC50", mixed_model_centrifugation_procedure_ic50, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/IC50/mixed_model_centrifugation_procedure_ic50.html", header_1 = 'Centrifugation Procedure - IC50')

#mixed_model_centrifugation_procedure_logic50 <- lmer(log_ic50 ~ centrifugation_procedure + (1 + centrifugation_procedure|drug) , df_transformed)
#summary(mixed_model_centrifugation_procedure_logic50)

###----plate_reader----
mixed_model_plate_reader_ic50 <- lmer(IC50_boxcox_sclaed2 ~ plate_reader  + (1|Patient_ID)+ (plate_reader|drug) , all_response_metrics)
summary(mixed_model_plate_reader_ic50)
check_model(mixed_model_plate_reader_ic50)
#model_interpretation(all_response_metrics, "IC50", mixed_model_plate_reader_ic50, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/IC50/mixed_model_plate_reader_ic50.html", header_1 = 'Plate Reader - IC50')

#mixed_model_plate_reader_logic50 <- lmer(log_ic50 ~ plate_reader + (1 + plate_reader|drug) , df_transformed)
#summary(mixed_model_plate_reader_logic50)

IC50_models <-c(mixed_model_time_until_sample_usage_ic50, mixed_model_medium_ic50, mixed_model_cells_ic50, mixed_model_micro_env_stimuli_ic50, mixed_model_cell_counting_method_ic50, mixed_model_sensitivity_readout_and_positive_control_ic50, mixed_model_centrifugation_procedure_ic50, mixed_model_plate_reader_ic50)
df_IC50_models <- model_results(IC50_models)
df_IC50_models$Random_Effect_Term <- '(1|Patient.num) + (1 | Experimental var | drug)'
rownames(df_IC50_models) <- NULL
formattable(df_IC50_models,  list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))

#----DSS1 ANALYSIS ###############################################################################################################################################################################################
##---- density plots ----
density_heatmaps_dss1 <- create_density_heatmaps(all_response_metrics_df_list, 'DSS1', drug_names_to_include)
d_plots_dss1 <- create_density_plots(all_response_metrics_df_list, 'drug', 'DSS1', drug_names_to_include, title = "Density plot of selected drugs across different datasets")

##---- Box Cox Scaling ----
all_response_metrics$DSS1_pos <- all_response_metrics$DSS1 + 0.01
MASS::boxcox(DSS1_pos ~ lab, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$DSS1_boxcox <-  ((all_response_metrics$DSS1)^0.5 - 1) / 0.5

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$DSS1_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$DSS1_boxcox[all_response_metrics$drug==j])
}

ggplot(df_transformed, aes(x=DSS1)) +
  geom_density(color="darkblue", fill="lightblue")


ggplot(all_response_metrics, aes(x=DSS1_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=DSS1_boxcox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")

##---- Factors ----
###----Time Until Sample Usage----
mixed_model_time_until_sample_usage_dss1 <- lmer(DSS1_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient.num) + (time_until_sample_usage|drug) , all_response_metrics)
summary(mixed_model_time_until_sample_usage_dss1)
check_model(mixed_model_time_until_sample_usage_dss1)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_time_until_sample_usage_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_time_until_sample_usage_dss1.html", header_1 = 'Time Until Sample Usage - DSS1', adj_nr = 9)

###----medium----
mixed_model_medium_dss1 <- lmer(DSS1_boxcox_sclaed2 ~ medium + (1|Patient.num) + (medium|drug) , all_response_metrics)
summary(mixed_model_medium_dss1)
check_model(mixed_model_medium_dss1)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_medium_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_medium_dss1.html", header_1 = 'Medium - DSS1')

###----cells----
mixed_model_cells_dss1 <- lmer(DSS1_boxcox_sclaed2 ~ cells + (1|Patient.num) + (cells|drug) , all_response_metrics)
summary(mixed_model_cells_dss1)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_cells_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_cells_dss1.html", header_1 = 'Nr Cells Used - DSS1')

###----nr_of_concentration_points----
mixed_model_conc_points_dss1 <- lmer(DSS1_boxcox_sclaed2 ~ nr_of_concentration_points + (1|Patient.num) +(nr_of_concentration_points|drug) , all_response_metrics)
summary(mixed_model_conc_points_dss1)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_conc_points_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_nr_conc_points_dss1.html", header_1 = 'Nr Concentration Points - DSS1')

###----micro_env_stimuli----
mixed_model_micro_env_stimuli_dss1 <- lmer(DSS1_boxcox_sclaed2 ~ microenvironmental_stimuli + (1|Patient.num) + (microenvironmental_stimuli|drug) , all_response_metrics)
summary(mixed_model_micro_env_stimuli_dss1)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_micro_env_stimuli_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_microenvironmental_stimuli_dss1.html", header_1 = 'Microenvironmental Stimuli - DSS1')

###----cell_counting_method----
mixed_model_cell_counting_method_dss1 <- lmer(DSS1_boxcox_sclaed2 ~ cell_counting_method + (1|Patient.num) + (cell_counting_method|drug) , all_response_metrics)
summary(mixed_model_cell_counting_method_dss1)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_cell_counting_method_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_cell_counting_method_dss1.html", header_1 = 'Cell Counting Method - DSS1')

###----sensitivity_readout_and_positive_control----
mixed_model_sensitivity_readout_and_positive_control_dss1 <- lmer(DSS1_boxcox_sclaed2 ~ positive_control + (1|Patient.num) + (positive_control|drug) , all_response_metrics)
summary(mixed_model_sensitivity_readout_and_positive_control_dss1)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_sensitivity_readout_and_positive_control_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_sensitivity_readout_and_positive_control_dss1.html", header_1 = 'Sensitivity Readout and Positive Control - DSS1')

###----centrifugation_procedure----
mixed_model_centrifugation_procedure_dss1 <- lmer(DSS1_boxcox_sclaed2 ~ centrifugation_procedure + (1|Patient.num) + (centrifugation_procedure|drug) , all_response_metrics)
summary(mixed_model_centrifugation_procedure_dss1)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_centrifugation_procedure_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_centrifugation_procedure_dss1.html", header_1 = 'Centrifugation Procedure - DSS1')


###----plate_reader----
mixed_model_plate_reader_dss1 <- lmer(DSS1_boxcox_sclaed2 ~ plate_reader + (1|Patient.num) + (plate_reader|drug) , all_response_metrics)
summary(mixed_model_plate_reader_dss1)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_plate_reader_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_plate_reader_dss1.html", header_1 = 'Plate Reader - DSS1')

###----DSS1 Table----
dss1_models <-c(mixed_model_time_until_sample_usage_dss1, mixed_model_medium_dss1, mixed_model_cells_dss1, mixed_model_micro_env_stimuli_dss1, mixed_model_cell_counting_method_dss1, mixed_model_sensitivity_readout_and_positive_control_dss1, mixed_model_centrifugation_procedure_dss1, mixed_model_plate_reader_dss1)
df_dss1_models <- model_results(dss1_models)
df_dss1_models$Random_Effect_Term <- '(1|Patient.num) + (1 | Experimental var | drug)'
rownames(df_dss1_models) <- NULL
formattable(df_dss1_models,  list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))

#----DSS3 ANALYSIS ###############################################################################################################################################################################################
##---- density plots ----
density_heatmaps_dss3 <- create_density_heatmaps(all_response_metrics_df_list, 'DSS3', drug_names_to_include)
d_plots_dss3 <- create_density_plots(all_response_metrics_df_list, 'drug', 'DSS3', drug_names_to_include, title = "Density plot of selected drugs across different datasets")

##---- Box Cox Scaling ----
all_response_metrics$DSS3_pos <- all_response_metrics$DSS3 + 0.01
MASS::boxcox(DSS3_pos ~ lab, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$DSS3_boxcox <- ((all_response_metrics$DSS3)^0.5 - 1) / 0.5

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$DSS3_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$DSS3_boxcox[all_response_metrics$drug==j])
}

ggplot(df_transformed, aes(x=DSS3)) +
  geom_density(color="darkblue", fill="lightblue")


ggplot(all_response_metrics, aes(x=DSS3_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=DSS3_boxcox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")


##---- Factors ----
###----Time Until Sample Usage----
mixed_model_time_until_sample_usage_dss3 <- lmer(DSS3_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient.num) + (time_until_sample_usage|drug) , all_response_metrics)
summary(mixed_model_time_until_sample_usage_dss3)
check_model(mixed_model_time_until_sample_usage_dss3)
#model_interpretation(all_response_metrics, "DSS3", mixed_model_time_until_sample_usage_dss3, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS3/mixed_model_time_until_sample_usage_dss3.html", header_1 = 'Time Until Sample Usage - DSS3', adj_nr = 9)

###----medium----
mixed_model_medium_dss3 <- lmer(DSS3_boxcox_sclaed2 ~ medium + (1|Patient.num) + (medium|drug) , all_response_metrics)
summary(mixed_model_medium_dss3)
#model_interpretation(all_response_metrics, "DSS3", mixed_model_medium_dss3, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS3/mixed_model_medium_dss3.html", header_1 = 'Medium - DSS3')

###----cells----
mixed_model_cells_dss3 <- lmer(DSS3_boxcox_sclaed2 ~ cells + (1|Patient.num) + (cells|drug) , all_response_metrics)
summary(mixed_model_cells_dss3)
#model_interpretation(all_response_metrics, "DSS3", mixed_model_cells_dss3, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS3/mixed_model_cells_dss3.html", header_1 = 'Nr Cells Used - DSS3')

###----micro_env_stimuli----
mixed_model_micro_env_stimuli_dss3 <- lmer(DSS3_boxcox_sclaed2 ~ microenvironmental_stimuli + + (1|Patient.num) + (microenvironmental_stimuli|drug) , all_response_metrics)
summary(mixed_model_micro_env_stimuli_dss3)
#model_interpretation(all_response_metrics, "DSS3", mixed_model_micro_env_stimuli_dss3, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS3/mixed_model_microenvironmental_stimuli_dss3.html", header_1 = 'Microenvironmental Stimuli - DSS3')

###----cell_counting_method----
mixed_model_cell_counting_method_dss3 <- lmer(DSS3_boxcox_sclaed2 ~ cell_counting_method + + (1|Patient.num) + (cell_counting_method|drug) , all_response_metrics)
summary(mixed_model_cell_counting_method_dss3)
#model_interpretation(all_response_metrics, "DSS3", mixed_model_cell_counting_method_dss3, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS3/mixed_model_cell_counting_method_dss3.html", header_1 = 'Cell Counting Method - DSS3')

###----sensitivity_readout_and_positive_control----
mixed_model_sensitivity_readout_and_positive_control_dss3 <- lmer(DSS3_boxcox_sclaed2 ~ positive_control + + (1|Patient.num) + (positive_control|drug) , all_response_metrics)
summary(mixed_model_sensitivity_readout_and_positive_control_dss3)
#model_interpretation(all_response_metrics, "DSS3", mixed_model_sensitivity_readout_and_positive_control_dss3, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS3/mixed_model_sensitivity_readout_and_positive_control_dss3.html", header_1 = 'Sensitivity Readout and Positive Control - DSS3')

###----centrifugation_procedure----
mixed_model_centrifugation_procedure_dss3 <- lmer(DSS3_boxcox_sclaed2 ~ centrifugation_procedure + + (1|Patient.num) + (centrifugation_procedure|drug) , all_response_metrics)
summary(mixed_model_centrifugation_procedure_dss3)
#model_interpretation(all_response_metrics, "DSS3", mixed_model_centrifugation_procedure_dss3, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS3/mixed_model_centrifugation_procedure_dss3.html", header_1 = 'Centrifugation Procedure - DSS3')

###----plate_reader----
mixed_model_plate_reader_dss3 <- lmer(DSS3_boxcox_sclaed2 ~ plate_reader + (1|Patient.num) + (plate_reader|drug) , all_response_metrics)
summary(mixed_model_plate_reader_dss3)
#model_interpretation(all_response_metrics, "DSS3", mixed_model_plate_reader_dss3, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS3/mixed_model_plate_reader_dss3.html", header_1 = 'Plate Reader - DSS3')

###----DSS3 Table----
dss3_models <-c(mixed_model_time_until_sample_usage_dss3, mixed_model_medium_dss3, mixed_model_cells_dss3, mixed_model_micro_env_stimuli_dss3, mixed_model_cell_counting_method_dss3, mixed_model_sensitivity_readout_and_positive_control_dss3, mixed_model_centrifugation_procedure_dss3, mixed_model_plate_reader_dss3)
df_dss3_models <- model_results(dss3_models)
df_dss3_models$Random_Effect_Term <- '(1|Patient.num) + (1 | Experimental var | drug)'
rownames(df_dss3_models) <- NULL
formattable(df_dss3_models,  list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))


#----auc Andrea way ANALYSIS ###############################################################################################################################################################################################

##---- Box Cox Scaling ----
all_response_metrics$auc_a[all_response_metrics$auc_a < 0] 

min(all_response_metrics$auc_a)
all_response_metrics$auc_a_pos <- all_response_metrics$auc_a - min(all_response_metrics$auc_a) + 0.1
MASS::boxcox(auc_a_pos ~ lab, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$auc_a_boxcox <- -((all_response_metrics$auc_a)^0.5 - 1) / 0.5

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$auc_a_box_cox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$auc_a_boxcox[all_response_metrics$drug==j])
}

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$auc_a_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$auc_a_boxcox[all_response_metrics$drug==j])
}

ggplot(all_response_metrics, aes(x=auc_a_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=auc_a_box_cox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(all_response_metrics, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")


##---- Factors ----
###----Time Until Sample Usage----
mixed_model_time_until_sample_usage_auc_a <- lmer(auc_a_box_cox_sclaed2 ~ time_until_sample_usage + (1|Patient.num) + (time_until_sample_usage|drug) , all_response_metrics)
summary(mixed_model_time_until_sample_usage_auc_a)
check_model(mixed_model_time_until_sample_usage_auc_a)
#model_interpretation(all_response_metrics, "auc_a", mixed_model_time_until_sample_usage_auc_a, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/auc_a/mixed_model_time_until_sample_usage_auc_a.html", header_1 = 'Time Until Sample Usage - auc a', adj_nr = 9)



###----medium----
mixed_model_medium_auc_a <- lmer(auc_a_box_cox_sclaed2 ~ medium + (1 + medium|drug) , all_response_metrics)
summary(mixed_model_medium_auc_a)
check_model(mixed_model_medium_auc_a)
#model_interpretation(all_response_metrics, "auc_a", mixed_model_medium_auc_a, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/auc_a/mixed_model_medium_dss3.html", header_1 = 'Medium - auc a')

###----cells----
mixed_model_cells_auc_a <- lmer(auc_a_box_cox_sclaed2 ~ cells +(1 + cells|drug) , all_response_metrics)
summary(mixed_model_cells_auc_a)
check_model(mixed_model_cells_auc_a)
#model_interpretation(all_response_metrics, "auc_a", mixed_model_cells_auc_a, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/auc_a/mixed_model_cells_dss3.html", header_1 = 'Nr Cells Used - auc a')

###----micro_env_stimuli----
mixed_model_micro_env_stimuli_auc_a <- lmer(auc_a_box_cox_sclaed2 ~ microenvironmental_stimuli + (1 + microenvironmental_stimuli|drug) , all_response_metrics)
summary(mixed_model_micro_env_stimuli_auc_a)
check_model(mixed_model_micro_env_stimuli_auc_a)
#model_interpretation(all_response_metrics, "auc_a", mixed_model_micro_env_stimuli_auc_a, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/auc_a/mixed_model_microenvironmental_stimuli_auc_a.html", header_1 = 'Microenvironmental Stimuli - auc a')

###----cell_counting_method----
mixed_model_cell_counting_method_auc_a <- lmer(auc_a_box_cox_sclaed2 ~ cell_counting_method + (1 + cell_counting_method|drug) , all_response_metrics)
summary(mixed_model_cell_counting_method_auc_a)
check_model(mixed_model_cell_counting_method_auc_a)
#model_interpretation(all_response_metrics, "auc_a", mixed_model_cell_counting_method_auc_a, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/auc_a/mixed_model_cell_counting_method_auc_a.html", header_1 = 'Cell Counting Method - auc a')

###----sensitivity_readout_and_positive_control----
mixed_model_sensitivity_readout_and_positive_control_auc_a <- lmer(auc_a_box_cox_sclaed2 ~ positive_control + (1 + positive_control|drug) , all_response_metrics)
summary(mixed_model_sensitivity_readout_and_positive_control_auc_a)
check_model(mixed_model_sensitivity_readout_and_positive_control_auc_a)
#model_interpretation(all_response_metrics, "auc_a", mixed_model_sensitivity_readout_and_positive_control_auc_a, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/auc_a/mixed_model_sensitivity_readout_and_positive_control_auc_a.html", header_1 = 'Sensitivity Readout and Positive Control - auc a')

###----centrifugation_procedure----
mixed_model_centrifugation_procedure_auc_a <- lmer(auc_a_box_cox_sclaed2 ~ centrifugation_procedure + (1 + centrifugation_procedure|drug) , all_response_metrics)
summary(mixed_model_centrifugation_procedure_auc_a)
check_model(mixed_model_centrifugation_procedure_auc_a)
#model_interpretation(all_response_metrics, "auc_a", mixed_model_centrifugation_procedure_auc_a, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/auc_a/mixed_model_centrifugation_procedure_auc_a.html", header_1 = 'Centrifugation Procedure - auc a')

###----plate_reader----
mixed_model_plate_reader_auc_a <- lmer(auc_a_box_cox_sclaed2 ~ plate_reader + (1 + plate_reader|drug) , all_response_metrics)
summary(mixed_model_plate_reader_auc_a)
check_model(mixed_model_plate_reader_auc_a)
#model_interpretation(all_response_metrics, "auc_a", mixed_model_plate_reader_auc_a, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/auc_a/mixed_model_plate_reader_auc_a.html", header_1 = 'Plate Reader - auc a')


library(nortest)
random_effects <- ranef(mixed_model_cells_box_cox_scaled_re_intercept)
random_effects$drug
random_resid <- random_effects$Patient.num[, 1]
ad.test(random_resid)
shapiro.test(random_resid)
ks.test(random_resid, "pnorm", mean = mean(random_resid), sd = sd(random_resid))
acf(resid(mixed_model_medium_box_cox_scaled_re_intercept))
library(performance)  
x <- check_heteroscedasticity(mixed_model_medium_box_cox_scaled_re_intercept)
performance::check_normality(mixed_model_medium_box_cox_scaled_re_intercept)  


#----Should the values be transformed----
plot_theme <- theme(text = element_text(family = "Arial", color= "black", size = 10),
                    plot.title = element_text(hjust = 0.5, vjust = 1,color= "black"),
                    axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                    axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                    axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                    axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                    plot.subtitle = element_blank(), 
                    legend.position = "top",  
                    legend.text = element_text(family = "Arial", color = "black", size = 10),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),  
                    panel.background = element_blank(),
                    strip.text = element_text(family = "Arial", color = "black", size = 8, face = "plain"),  
                    strip.background = element_blank(), 
                    plot.margin = margin(10, 10, 10, 10))
metrics <- c("IC50", "DSS1", "DSS2", "DSS3", "auc_a", "AUC")
for (m in metrics[1:6]){
  
  x_lab_m <- ifelse(m == "auc_a", "rAUC", m)
  x_lab_m <- ifelse(x_lab_m == "IC50", expression("IC"[50]), x_lab_m)
  x_lab_m <- ifelse(x_lab_m == "DSS1", expression("DSS"[1]), x_lab_m)
  x_lab_m <- ifelse(x_lab_m == "DSS2", expression("DSS"[2]), x_lab_m)
  x_lab_m <- ifelse(x_lab_m == "DSS3", expression("DSS"[3]), x_lab_m)
  
  dist_untransformed <- ggplot(all_response_metrics, aes(x=get(m))) +
    geom_density(color="darkblue", fill="lightblue") + 
    theme_minimal() + 
    xlab(x_lab_m) +
    ylab("Density")
  print(dist_untransformed)
  ggsave(paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Comparisons/", m, "/dist_untransformed.png"), plot = dist_untransformed, width = 8.1, height = 8.2, units="cm")
  
  m_boxcox <- paste0(m,"_boxcox")
  dist_boxcox <- ggplot(all_response_metrics, aes(x=get(m_boxcox))) +
    geom_density(color="darkblue", fill="lightblue") + theme_minimal() +
    ylab("Density") +
    xlab(paste("BoxCox transformed", m))
  print(dist_boxcox)
  ggsave(paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Comparisons/", m, "/ddist_boxcox.png"), plot = dist_boxcox, width = 8.1, height = 8.2, units="cm")

  m_boxcox_scaled2 <- paste0(m,"_boxcox_sclaed2")
  dist_boxcox_scaled2 <- ggplot(all_response_metrics, aes(x=get(m_boxcox_scaled2))) + #IC50_boxcox_sclaed2
    geom_density(color="darkblue", fill="lightblue") +
    theme_minimal() +
    ylab("Density") +
    xlab("BoxCox transformed z-Scaled IC50")
  print(dist_boxcox_scaled2)
  ggsave(paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Comparisons/", m, "/dist_boxcox_scaled2.png"), plot = dist_boxcox_scaled2, width = 8.1, height = 8.2, units="cm")


  print(m)
  m_box <- paste0(m, "_boxcox_sclaed2")
  # Construct formula as a string
  formula_str <- paste0(m_box, " ~ medium + (1|Patient.num)", " + (medium|drug)")
  # Convert string to formula
  f <- as.formula(formula_str)
  model <- lmer(f,all_response_metrics)
  homogeneity_plot <- check_model(model, check = "homogeneity", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
  plot_homogeneity <- plot(homogeneity_plot)
  plot_homogeneity_adjusted <- plot_homogeneity$HOMOGENEITY + labs(y = expression(sqrt("|Std. Residuals|")))+ plot_theme
  print(plot_homogeneity)
  ggsave(paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Comparisons/", m, "/BoxCox_scaled_plot_homogeneity_adjusted.png"), plot = plot_homogeneity_adjusted, width = 8.1, height = 8.2, units="cm")

  pp_check_plot <- check_model(model, check = "pp_check", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3", "black"))
  plot_pp_check <- plot(pp_check_plot)
  plot_pp_check_adjusted <- plot_pp_check$PP_CHECK + plot_theme + labs(x = paste("BoxCox Transformed and Scaled", m))
  print(plot_pp_check_adjusted)
  ggsave(paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Comparisons/", m, "/BoxCox_scaled_pp_check_plot.png"), plot = plot_pp_check_adjusted,width = 8.1, height = 8.2, units="cm")

  untransformed_formula_str <- paste0(m, " ~ medium + (1|Patient.num)", " + (medium|drug)")
  # Convert string to formula
  untransformed_f <- as.formula(untransformed_formula_str)
  untransformed_model <- lmer(untransformed_f,all_response_metrics)

  untransformed_homogeneity_plot <- check_model(untransformed_model, check = "homogeneity", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
  untransformed_plot_homogeneity <- plot(untransformed_homogeneity_plot)
  untransformed_plot_homogeneity_adjusted <- untransformed_plot_homogeneity$HOMOGENEITY + labs(y = expression(sqrt("|Std. Residuals|")))+ plot_theme
  print(untransformed_plot_homogeneity)
  ggsave(paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Comparisons/",m, "/untransformed_plot_homogeneity_adjusted.png"), plot = untransformed_plot_homogeneity_adjusted, width = 8.1, height = 8.2, units="cm")

  untransformed_pp_check_plot <- check_model(untransformed_model, check = "pp_check", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3", "black"))
  untransformed_plot_pp_check <- plot(untransformed_pp_check_plot)
  untransformed_plot_pp_check_adjusted <- untransformed_plot_pp_check$PP_CHECK + plot_theme + labs(x = m)
  print(untransformed_plot_pp_check_adjusted)
  ggsave(paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Comparisons/", m, "/untransformed_pp_check_plot.png"), plot = untransformed_plot_pp_check_adjusted,width = 8.1, height = 8.2, units="cm")

  
  
}


#----Bar plot of all response scores ---- #############################################################
# List of models and their names

variables <- c("time_until_sample_usage", "medium", "cells", "positive_control", "centrifugation_procedure", "plate_reader")
metrics <- c("IC50", "DSS1", "DSS2", "DSS3", "auc_a", "AUC")
#metrics <- c("IC50")
models <- list()
for (m in metrics[1:6]){
  for (v in variables){
    model_name_plot <- case_when(v == "time_until_sample_usage" ~ "time until usage", 
                   v == "microenvironmental_stimuli" ~ "microenvironmental stimuli", 
                   v == "positive_control" ~ "positive control-doses-readout-cell counting", 
                   v == "centrifugation_procedure" ~ "centrifugation procedure", 
                   v == "plate_reader" ~ "plate reader",
                   .default = v)
    print(m)
    print(v)
    m_box <- paste0(m, "_boxcox_sclaed2")
    # Construct formula as a string
    formula_str <- paste0(m_box, " ~ ", v, " + (1|Patient.num)", " + (",v,"|drug)")
    # Convert string to formula
    f <- as.formula(formula_str)
    model_name <- paste0(m, "_by_", v)
    models[[model_name]] <- lmer(f,all_response_metrics)
    #model_diagnostics(models[[model_name]], paste0("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/", m, "/", v, "/"), model_name = model_name_plot, metric = m) 
    print(summary(models[[model_name]]))
  }
}
#models <- list("Time Until Sample Usage" = list(DSS2 = mixed_model_time_until_sample_usage, DSS3 = mixed_model_time_until_sample_usage_dss3, DSS1 = mixed_model_time_until_sample_usage_dss1, AUC = mixed_model_time_until_sample_usage_auc, IC50 = mixed_model_time_until_sample_usage_ic50),"" = list(DSS2 = ))
models <- c(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept, mixed_model_time_until_sample_usage_dss3, mixed_model_time_until_sample_usage_dss1, mixed_model_time_until_sample_usage_auc, mixed_model_time_until_sample_usage_ic50, mixed_model_time_until_sample_usage_auc_a,
            mixed_model_medium_box_cox_scaled_re_intercept, mixed_model_medium_dss3, mixed_model_medium_dss1, mixed_model_medium_auc,mixed_model_medium_ic50, mixed_model_medium_auc_a,
            mixed_model_cells_box_cox_scaled_re_intercept, mixed_model_cells_dss3, mixed_model_cells_dss1, mixed_model_cells_auc, mixed_model_cells_ic50,mixed_model_cells_auc_a,
            #mixed_model_conc_points, mixed_model_conc_points_dss3, mixed_model_conc_points_dss1, mixed_model_conc_points_auc, mixed_model_conc_points_ic50, mixed_model_conc_points_logic50,mixed_model_conc_points_auc_a,
            mixed_model_micro_env_stimuli_box_cox_scaled_re_intercept, mixed_model_micro_env_stimuli_dss3, mixed_model_micro_env_stimuli_dss1, mixed_model_micro_env_stimuli_auc, mixed_model_micro_env_stimuli_ic50, mixed_model_micro_env_stimuli_auc_a,
            mixed_model_cell_counting_method_box_cox_scaled_re_intercept, mixed_model_cell_counting_method_dss3, mixed_model_cell_counting_method_dss1, mixed_model_cell_counting_method_auc,mixed_model_cell_counting_method_ic50, mixed_model_cell_counting_method_auc_a,
            mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_re_intercept, mixed_model_sensitivity_readout_and_positive_control_dss3, mixed_model_sensitivity_readout_and_positive_control_dss1, mixed_model_sensitivity_readout_and_positive_control_auc,mixed_model_sensitivity_readout_and_positive_control_ic50,mixed_model_sensitivity_readout_and_positive_control_auc_a,
            mixed_model_centrifugation_procedure_box_cox_scaled_re_intercept, mixed_model_centrifugation_procedure_dss3, mixed_model_centrifugation_procedure_dss1, mixed_model_centrifugation_procedure_auc,mixed_model_centrifugation_procedure_ic50, mixed_model_centrifugation_procedure_auc_a,
            mixed_model_plate_reader_box_cox_scaled_re_intercept, mixed_model_plate_reader_dss3, mixed_model_plate_reader_dss1, mixed_model_plate_reader_auc, mixed_model_plate_reader_ic50, mixed_model_plate_reader_auc_a)

used_data <- model.frame(mixed_model_cells_box_cox_scaled_re_intercept)
unique(used_data$Patient.num)
dropped_rows <- setdiff(rownames(all_response_metrics), rownames(used_data))
dropped_rows <- all_response_metrics[dropped_rows, ] 

all_response_metrics$row_id <- seq_len(nrow(all_response_metrics))
mod <- lmer(DSS2_boxcox_sclaed2 ~ cells + (1 | Patient.num) + (cells|drug), data = all_response_metrics)
used_data <- as.data.frame(model.frame(mod))

# Check dropped IDs
setdiff(all_response_metrics$row_id, used_data$row_id)
length(setdiff(all_response_metrics$row_id, used_data$row_id))
summary(comparedf(used_data[,c('Patient.num', 'drug')],all_response_metrics[,c('Patient.num', 'drug')], by = 'drug'))

df <- model_results(models, data_frame_org = all_response_metrics)

#all_metric_df <- df
#df <- all_metric_df
#df <- df[df$Response_metric != 'IC50_boxcox_sclaed2',]
#df <- rbind(df, df_IC50_models)
df <- df %>% mutate(corrected_p_value = ifelse(df$p_value * df$count >= 1, 1, df$p_value * df$count))
df$corrected_p_value <- round(df$corrected_p_value, 4)
df$corrected_p_value <- format(df$corrected_p_value, scientific = FALSE, trim = TRUE)
df <- df %>% mutate(corrected_p_value = ifelse(df$corrected_p_value == '0.0000', '<0.0001', as.character(df$corrected_p_value)))

df <- df %>% mutate(ref_group = case_when(ref_group ==  "a drug combination of flavopiridol, staurosporine and velcade" ~ "Positive Control: drug combination + Doses: 7 \n+ Readout: CellTiter96     ", 
                                          ref_group ==  "1-2h after receiving" ~ "Time until sample usage < 2h",
                                          ref_group ==  "HS-5 conditioned medium" ~ "HS-5 CM",
                                          ref_group ==  "HS-5 CM" ~ "HS-5 CM",
                                          ref_group ==  "Ficoll-Paque centrifugation" ~ "Ficoll-Paque",
                                          ref_group ==  "10000" ~ "10000",
                                          ref_group ==  "Countess" ~ "Countess",
                                          ref_group ==  "Biotek Synergy 2" ~ "Biotek Synergy 2",
                                          TRUE ~ ref_group))

df$Random_Effect_Term <- '(1 | Patient.num) + (1 + Experimental var | drug)'
colnames(df) <- gsub("_", " ", colnames(df))     
colnames(df) <- sapply(colnames(df), tools::toTitleCase)
df <- df[order(df$`p Value`),]
rownames(df) <- NULL

# df <- df %>% mutate(`Experimental Variable`= case_when(`Fixed Effect Term` == "Positive Control: BzCl + Nr of Concentration Points: 5 + Drug Sensitivity Readout: CellTiter-Glo" ~ "Positive Control-Doses-Readout   ", 
#                                                        `Fixed Effect Term` ==  "Time until sample usage: Handled within 2-72 hours" ~ "Time Until Sample Usage > 2h   ",
#                                                        `Fixed Effect Term` ==  "Medium: Mononuclear cell medium" ~ "Medium: MCM    ",
#                                                        `Fixed Effect Term` ==  "Medium: RPMI + fetal bovine serum (FBS) (10%)" ~ "Medium: RPMI + FBS   ",
#                                                        `Fixed Effect Term` ==  "microenvironmental_stimuliNone" ~ "Microenvironmental Stimuli: None    ",
#                                                        `Fixed Effect Term` ==  "Microenvironmental stimuli: Transiently cultured with feeder cells, activate samples with autologous BM T helper cells in the presence of IL-2 and a T-cell expansion cocktai" ~ "Microenvironmental Stimuli: Co-culture with activation   ",
#                                                        `Fixed Effect Term` ==  "Plate readers: Pherastar, Cytation5, Insight, Tecan" ~ "Plate Reader",
#                                                        `Fixed Effect Term` ==  "Plate Reader" ~ "Plate Reader: Pherastar",
#                                                        `Fixed Effect Term` ==  "Centrifugation procedure: LymphoPrepTM gradient centrifugation" ~ "Centrifugation Procedure: 1ᵃ    ",
#                                                        `Fixed Effect Term` ==  "Centrifugation procedure: Supernatant isolation at 300g 10min, density centrifugation at 400g for 20min without brake, afterwards always 300g 5 min" ~ "Centrifugation Procedure: 2ᵇ   ",
#                                                        `Fixed Effect Term` ==  "Cells: 5000" ~ "Number of Cells per Well: 5000   ",
#                                                        `Fixed Effect Term` ==  "Cell counting method: Trypan blue (Countess™ II FL Automated Cell Counter)" ~ "Cell Counting Method: Trypan blue    ",
#                                                        `Fixed Effect Term` == "cell_counting_methodoptical density" ~ "Cell Counting Method: Optical density    ",
#                                                        `Fixed Effect Term` ==  "plate_readerVictorX" ~ "Plate Reader: Victor X    ",
#                                                        `Fixed Effect Term` == "plate_readerEnVision" ~ "Plate Reader: EnVision    ",
#                                                        `Fixed Effect Term` == "plate_readerEnSight" ~ "Plate Reader: EnSight    ",
#                                                        TRUE ~ `Fixed Effect Term`))


df <- df %>% mutate(`Fixed Effect Term` = case_when(`Fixed Effect Term` ==  "Positive Control: BzCl + Doses: 5 + Readout: CellTiter-Glo" ~ "Positive Control: BzCl + Doses: 5 \n+ Readout: CellTiter-Glo",
                                                    `Fixed Effect Term` == "Microenvironmental Stimuli: Co-culture with activation" ~ "Microenvironmental Stimuli:\nCo-culture with activation", 
                                                    `Fixed Effect Term` == "Centrifugation Procedure: LymphoPrepTM gradient centrifugation" ~ "Centrifugation Procedure: \nLymphoPrep\u2122 gradient", 
                                                    TRUE ~ `Fixed Effect Term`))
df$`Fixed Effect Term`
library(forestploter)

df$`Experimental Variable` <- df$`Fixed Effect Term`
df$`Reference Condition` <- df$`Ref Group`
df <- df[order(df$`Fixed Effect Coefficient`),]
df <- df %>%
  arrange('Fixed Effect Coefficient', 'Response Metric', 'Experimental Variable')
df_wide <- df[,c('Experimental Variable', 'Reference Condition', 'Response Metric', 'Fixed Effect Coefficient', 'Lower', 'Upper', 'Corrected p Value')] %>%
  pivot_wider(names_from = `Response Metric`, values_from = c(`Fixed Effect Coefficient`, Lower, Upper, `Corrected p Value`)) 
df_wide$` ` <- paste(rep(" ", 30), collapse = " ")

# Create a color theme vector to map colors to each category (Cat1, Cat2, etc.)
category_colors <- c("lightgreen","skyblue", "royalblue", "orange", "salmon", "pink")
category_colors <- c("#8dd3c7", "#fccde5", "#bebada", "#fb8072", "#fdb462", "#b3de69")

unique(df$`Response Metric`)
df$Color <- category_colors[df$`Response Metric`]

# Define theme for the forest plot
tm <- forest_theme(
  base_size = 12,
  arrow_type = "closed",
  line_size = 1.5, 
  align = "center", 
  ci_fill = category_colors,
  footnote_gp = gpar(col = "black", cex = 1.0, fontsize = 18), 
  footnote_parse = FALSE, 
  legend_name = "   ",#Response Metric
  legend_position = "bottom",
  legend_value = c(expression("IC"[50]), expression(" DSS"[3]), expression(" DSS"[2]), expression(" DSS"[1]), " rAUC", " AUC"), 
  legend_gp = gpar(fontsize = 18, fontfamily = "Arial", cex = 1),
  xaxis_gp = gpar(fontsize = 18, fontfamily = "Arial"), 
  spacing = 3,
)

names(df_wide) <- str_replace_all(names(df_wide), '_', ' ')
names(df_wide) <- str_replace_all(names(df_wide), ':', '')

df_wide <- df_wide %>% dplyr::rename(`DSS2` = `Corrected p Value DSS2 Box Cox Transformed and Scaled`,
                              `IC50` = `Corrected p Value IC50 boxcox sclaed2`)

df_wide <- df_wide[order(df_wide$`Fixed Effect Coefficient DSS2 Box Cox Transformed and Scaled`),]
df_wide$`p-valueᶜ` <- paste(rep(" ", 30), collapse = " ")

# Create the forest plot with different colors for each category bar
p1_all_metric <- forest(
  df_wide[, c("Experimental Variable","Reference Condition", " ", "DSS2", "IC50")],              # Main label column
  est = list(df_wide$`Fixed Effect Coefficient IC50 boxcox sclaed2`, df_wide$`Fixed Effect Coefficient DSS3 boxcox sclaed2`, df_wide$`Fixed Effect Coefficient DSS2 Box Cox Transformed and Scaled`, df_wide$`Fixed Effect Coefficient DSS1 boxcox sclaed2`, df_wide$`Fixed Effect Coefficient auc a boxcox sclaed2`, df_wide$`Fixed Effect Coefficient AUC boxcox sclaed2`),
  lower = list(df_wide$`Lower IC50 boxcox sclaed2`, df_wide$`Lower DSS3 boxcox sclaed2`, df_wide$`Lower DSS2 Box Cox Transformed and Scaled`, df_wide$`Lower DSS1 boxcox sclaed2`, df_wide$`Lower auc a boxcox sclaed2`, df_wide$`Lower AUC boxcox sclaed2`),   # Lower CIs for each category
  upper = list(df_wide$`Upper IC50 boxcox sclaed2`, df_wide$`Upper DSS3 boxcox sclaed2`, df_wide$`Upper DSS2 Box Cox Transformed and Scaled`, df_wide$`Upper DSS1 boxcox sclaed2`, df_wide$`Upper auc a boxcox sclaed2`, df_wide$`Upper AUC boxcox sclaed2`),   # Upper CIs for each category
  ci_column = 3,                                      # CI column placement
  ref_line = 0,   
  sizes = 1,
  grid = F,
  xlim = c(-1.2, 1.3),
  #ticks_at = c(-2, -1, 0, 1, 2),
  footnote = "\nᵃBonferroni corrected",
  font.label = list(size = 18, family = "Arial"),
  font.ticks = list(size = 18, family = "Arial"),
  nudge_y = 0.15,
  theme = tm
)


output_width_cm <- 19.4 #12
output_height_cm <- 25 #10 19
dpi <- 300  # Resolution in dots per inch

# Convert cm to inches (1 inch = 2.54 cm)
output_width_in <- output_width_cm / 2.54
output_height_in <- output_height_cm / 2.54

# Convert inches to pixels for the PNG device
output_width_px <- output_width_in * dpi
output_height_px <- output_height_in * dpi

# Calculate scaling factors for the gtable
scale_width <- output_width_cm / 10  # Base width adjustment (10 cm as reference)
scale_height <- output_height_cm / 10  # Base height adjustment (10 cm as reference)

p1_all_metric <- edit_plot(p1_all_metric, gp = gpar(cex=2, fontfamily="Arial")) #1.4
p1_all_metric <- edit_plot(p1_all_metric, part = "header", gp = gpar(cex=2, fontfamily="Arial"))

p1_all_metric <- edit_plot(p1_all_metric, col = 4:5, part = "header",
               which = "text",
               hjust = unit(0.5, "npc"),
               x = unit(0.5, "npc"))
p1_all_metric <- edit_plot(p1_all_metric, col = 4:5, part = "body",
               which = "text",
               hjust = unit(0.5, "npc"),
               x = unit(0.5, "npc"))

p1_all_metric <- edit_plot(p1_all_metric,
          col = 4,           # columns of header
          part = "header",     # part of the gtable to edit
          which = "text",      # what to edit in the grob
          label = "DSS₂", # New text with subscript 2 (U+2082)
          hjust = unit(0.5, "npc"),
          x = unit(0.5, "npc"), 
          gp = gpar(fontface = "bold"))
p1_all_metric <- edit_plot(p1_all_metric,
                           col = 5,           # columns of header
                           part = "header",     # part of the gtable to edit
                           which = "text",      # what to edit in the grob
                           label = "IC₅₀", # New text with subscript 2 (U+2082)
                           hjust = unit(0.5, "npc"),
                           x = unit(0.5, "npc"), 
                           gp = gpar(fontface = "bold"))
p1_all_metric <- add_text(p1_all_metric, text = "p-valueᵃ",
          part = "header", 
          row = 0,
          col = 4:5,
          gp = gpar(cex=2.5, fontface = "bold", fontfamily = "Arial"))

class(p1_all_metric$grobs[[5]])
p1_all_metric$grobs[[5]]$gp <- gpar(col="black",cex=2,fontfamily="Arial",fontsize=12,lineheight=1.2,alpha=1,font=2)
p1_all_metric$grobs[[1]]$gp
# p1_all_metric <- insert_text(p1_all_metric, text = expression("DSS"[2]), part = "header", row = 1, col = 4, before = TRUE, gp = gpar(fontface = "bold"))
legend_index <- 174
p1_all_metric$layout$t[21]
legend_index <- which(p1_all_metric$layout$name == "legend")
p1_all_metric$layout$t[legend_index] <- 12 #nrow(p1_all_metric)       # or your desired row index
p1_all_metric$layout$b[legend_index] <- 15 #nrow(p1_all_metric)
p1_all_metric$layout$l[legend_index] <- 4             # column 4
p1_all_metric$layout$r[legend_index] <- 4

# Scale the gtable layout
scaled_p1_all_metric <- gtable::gtable_filter(p1_all_metric, pattern = ".*", trim = TRUE)  # Keep all grobs
scaled_p1_all_metric$widths <- scaled_p1_all_metric$widths * scale_width
scaled_p1_all_metric$heights <- scaled_p1_all_metric$heights * scale_height
scaled_p1_all_metric


figure_output <- 'Desktop/UiO/Project 1/Figures/New karolinska data/'

png(paste0(figure_output,'Forest_plot_all_metric.png'), height = 36, width = 56, unit = "cm", res = 300) #24.8 #height = 16, width = 46 ; height = 10, width = 28
grid.newpage()
grid.draw(scaled_p1_all_metric)
dev.off()
dev.off()

add_border(
  p1,
  row = NULL,
  col = NULL,
  part = "header",
  where = "bottom",
  gp = gpar(lwd = 2)
)

plot(p)
plot(p1)
p_wh <- get_wh(p1)
pdf('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/Forest_plot_all_metric_test.pdf',width = p_wh[1], height = p_wh[2])
plot(p1)
dev.off()


library(forestplot)
df <- df[!duplicated(df$`Experimental Variable`), ]
df_tbl <- as_tibble(df)

df_tbl[,c('Experimental Variable', 'Response Metric', ' ', 'Fixed Effect Coefficient', 'Lower', 'Upper', 'Corrected p Value')] |>
  group_by(`Response Metric`) |>
  forestplot::forestplot(labeltext = `Experimental Variable`, 
             mean = `Fixed Effect Coefficient`, 
             lower = `Lower`, 
             upper = `Upper`,
             line.margin = .1,
            clip = c(-5, 5),
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,
             boxsize = 1,
             lineheight = "lines",
             xlab = "EQ-5D index") |>
  #fp_add_lines("steelblue") |>
  forestplot::fp_add_header("Variable") |>
  forestplot::fp_set_style(box = c("blue", "darkred", "salmon", "lightgreen", "violet", "orange") |> lapply(function(x) gpar(fill = x, col = "#555555")))
data(dfHRQoL)

duplicated_labels <- anyDuplicated(df$`Experimental Variable`)
if (duplicated_labels > 0) {
  print("Duplicate labels detected!")
  print(df$`Experimental Variable`[duplicated_labels])
}

# Check alignment
print(df)
####----Bar plot ----
p <- ggplot(subset(df, `Response Metric` == "DSS2: Box Cox Transformed and Scaled"), 
            aes(x = `Fixed Effect Coefficient`, y = `Fixed Effect Term`)) + 
  geom_bar(stat = "identity", position = "dodge", fill = "#80b1d3", color = "#80b1d3", width = 0.7) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), 
                 height = 0.2, color = "black", alpha = 0.5) +  # Horizontal error bars
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(size = 14, family = "Arial", color = "black"),  # Set plot title text
    panel.background = element_rect(fill = "white"),  # White background
    panel.grid = element_blank(),  # Remove gridlines
    plot.background = element_rect(fill = "white"),  # White plot background
    axis.line = element_line(color = rgb(0, 0, 0, alpha = 0.5), size = 0.5),  # Less vibrant axis lines
    axis.ticks = element_line(color = rgb(0, 0, 0, alpha = 0.5), size = 0.5),  # Same color and vibrancy for ticks
    axis.ticks.length = unit(-0.1, "cm")  # Smaller ticks inside the plot
  )

print(p)


# Save the plot to a file
ggsave("~/Desktop/UiO/Project 1/Figures/fixed_effects_plot_all_metrics.png", plot = p, width = 10, height = 15, dpi = 300)



########################################
df <- data.frame(
  Fixed_Effect_Term = character(),
  Response_metric = character(),
  Fixed_Effect_Coefficient = numeric(),
  p_value = numeric(),
  Intercept = numeric(),
  stringsAsFactors = FALSE
)

# Step 3: Extract fixed effects and their coefficients
for (model in models) {
  model_formula <- formula(model)
  dependent_var <- as.character(model_formula[[2]])
  fixed_effects <- fixef(model)  # Get fixed effects
  fixed_effect_terms <- names(fixed_effects)
  fixed_effect_names <- fixed_effect_terms[fixed_effect_terms != "(Intercept)"]
  fixed_effect_coef <- fixed_effects[fixed_effect_names]
  coef_model <- summary(model)$coefficients
  p_values <- coef_model[, "Pr(>|t|)"]
  p_values_excluding_intercept <- p_values[fixed_effect_terms != "(Intercept)"]
  intercept <- fixed_effects[fixed_effect_terms[fixed_effect_terms == "(Intercept)"]]
  response_metric <- dependent_var
  
  # Create a temporary data frame for this model's fixed effects
  temp_df <- data.frame(
    Fixed_Effect_Term = fixed_effect_names,
    Response_metric = response_metric,
    Fixed_Effect_Coefficient = fixed_effect_coef,
    p_value = p_values_excluding_intercept,
    Intercept = intercept,
    stringsAsFactors = FALSE
  )
  
  # Step 4: Add the temporary data frame to the main data frame
  df <- rbind(df, temp_df)
}

get_significance_stars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return(" ")
  }
}

# Add significance stars to the data frame
df$significance_stars <- sapply(df$p_value, get_significance_stars)
# Count the number of fixed effects per category in col1
fixed_effect_counts <- df %>%
  group_by(Response_metric) %>%
  summarise(fixed_effect_counts = n())

# Join the counts back to the original data frame
df <- df %>%
  left_join(fixed_effect_counts, by = "Response_metric")

# Create the corrected p-value column
df$corrected_p_value <- df$p_value * df$fixed_effect_counts

print(df)

####----table pretty ----
library(formattable)
rownames(df) <- NULL
df <- df %>% mutate(Fixed_Effect_Term = case_when(Fixed_Effect_Term ==  "positive_controlBzCl" ~ "Positive Control: BzCl + Nr of Concentration Points: 5 + Drug Sensitivity Readout: CellTiter-Glo", 
                                                  TRUE ~ Fixed_Effect_Term))
formattable(subset(df[order(-df$Fixed_Effect_Coefficient), ], Response_metric == "DSS2_boxcox_sclaed2"))



###################################################################################################
#----Stratified analysis----
all_response_metrics$DSS2_pos <- all_response_metrics$DSS2 + 0.01
MASS::boxcox(DSS2_pos ~ time_until_sample_usage, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$DSS2_boxcox <- ((all_response_metrics$DSS2)^0.5 - 1) / 0.5

#Standardize DSS scores w.r.t. labs and drugs
for(i in unique(all_response_metrics$lab)){
  for(j in unique(all_response_metrics$drug)){
    all_response_metrics$DSS2_boxcox_sclaed1[all_response_metrics$lab==i & all_response_metrics$drug==j] <-
      scale(all_response_metrics$DSS2_boxcox[all_response_metrics$lab==i & all_response_metrics$drug==j])
  }
}

for(j in unique(all_response_metrics$drug)){
  all_response_metrics$DSS2_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$DSS2_boxcox[all_response_metrics$drug==j])
}

df_stratified <- subset(all_response_metrics, Disease.status == 'Diagnosis' & sample == 'fresh' & Tissue == 'Bone marrow')
unique(df_stratified$medium)
subset(df_stratified, lab == 'Oslo')
length(unique(df_stratified$Patient.num))
unique(df_stratified$lab)

group_counts <- df_stratified %>%
  group_by(lab, drug) %>%
  dplyr::summarize(count = n())

##---- Factors ----
###----Time Until Sample Usage----
mixed_model_time_until_sample_usage_dss2_stratified <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1|Patient.num) + (time_until_sample_usage|drug) , df_stratified)
summary(mixed_model_time_until_sample_usage_dss2_stratified)
check_model(mixed_model_time_until_sample_usage_dss2_stratified)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_time_until_sample_usage_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_time_until_sample_usage_dss1.html", header_1 = 'Time Until Sample Usage - DSS1', adj_nr = 9)

###----medium----
mixed_model_medium_dss2_stratified <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1|Patient.num) + (medium|drug) , df_stratified)
summary(mixed_model_medium_dss2_stratified)
check_model(mixed_model_medium_dss2_stratified)

#model_interpretation(all_response_metrics, "DSS1", mixed_model_medium_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_medium_dss1.html", header_1 = 'Medium - DSS1')

###----cells----
mixed_model_cells_dss2_stratified <- lmer(DSS2_boxcox_sclaed2 ~ cells + (1|Patient.num) + (cells|drug) , df_stratified)
summary(mixed_model_cells_dss2_stratified)
check_model(mixed_model_cells_dss2_stratified)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_cells_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_cells_dss1.html", header_1 = 'Nr Cells Used - DSS1')

###----micro_env_stimuli----
mixed_model_micro_env_stimuli_dss2_stratified <- lmer(DSS2_boxcox_sclaed2 ~ microenvironmental_stimuli + (1|Patient.num) + (microenvironmental_stimuli|drug) , df_stratified)
summary(mixed_model_micro_env_stimuli_dss2_stratified)
check_model(mixed_model_micro_env_stimuli_dss2_stratified)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_micro_env_stimuli_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_microenvironmental_stimuli_dss1.html", header_1 = 'Microenvironmental Stimuli - DSS1')

###----cell_counting_method----
mixed_model_cell_counting_method_dss2_stratified <- lmer(DSS2_boxcox_sclaed2 ~ cell_counting_method + (1|Patient.num) + (cell_counting_method|drug) , df_stratified)
summary(mixed_model_cell_counting_method_dss2_stratified)
check_model(mixed_model_cell_counting_method_dss2_stratified)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_cell_counting_method_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_cell_counting_method_dss1.html", header_1 = 'Cell Counting Method - DSS1')

###----sensitivity_readout_and_positive_control----
mixed_model_sensitivity_readout_and_positive_control_dss2_stratified <- lmer(DSS2_boxcox_sclaed2 ~ positive_control + (1|Patient.num) + (positive_control|drug) , df_stratified)
summary(mixed_model_sensitivity_readout_and_positive_control_dss2_stratified)
check_model(mixed_model_sensitivity_readout_and_positive_control_dss2_stratified)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_sensitivity_readout_and_positive_control_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_sensitivity_readout_and_positive_control_dss1.html", header_1 = 'Sensitivity Readout and Positive Control - DSS1')

###----centrifugation_procedure----
mixed_model_centrifugation_procedure_dss2_stratified <- lmer(DSS2_boxcox_sclaed2 ~ centrifugation_procedure + (1|Patient.num) + (centrifugation_procedure|drug) , df_stratified)
summary(mixed_model_centrifugation_procedure_dss2_stratified)
check_model(mixed_model_centrifugation_procedure_dss2_stratified)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_centrifugation_procedure_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_centrifugation_procedure_dss1.html", header_1 = 'Centrifugation Procedure - DSS1')


###----plate_reader----
mixed_model_plate_reader_dss2_stratified <- lmer(DSS2_boxcox_sclaed2 ~ plate_reader + (1|Patient.num) + (plate_reader|drug) , df_stratified)
summary(mixed_model_plate_reader_dss2_stratified)
check_model(mixed_model_plate_reader_dss2_stratified)
#model_interpretation(all_response_metrics, "DSS1", mixed_model_plate_reader_dss1, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/DSS1/mixed_model_plate_reader_dss1.html", header_1 = 'Plate Reader - DSS1')



###----DSS2 Straified Table----
dss2_models_stratified <-c(mixed_model_time_until_sample_usage_dss2_stratified, mixed_model_medium_dss2_stratified, mixed_model_cells_dss2_stratified, mixed_model_micro_env_stimuli_dss2_stratified, mixed_model_cell_counting_method_dss2_stratified, mixed_model_sensitivity_readout_and_positive_control_dss2_stratified, mixed_model_centrifugation_procedure_dss2_stratified, mixed_model_plate_reader_dss2_stratified)
df_dss2_models_stratified <- model_results(dss2_models_stratified)
df_dss2_models_stratified$Random_Effect_Term <- '(1|Patient.num) + (1 | Experimental var | drug)'
rownames(df_dss2_models_stratified) <- NULL
formattable(df_dss2_models_stratified,  list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))
df_dss2_models_stratified_backup <- df_dss2_models_stratified
df_dss2_models_stratified$Random_Effect_Term <- '(1 | Patient.num) + (1 + Experimental var | drug)'
colnames(df_dss2_models_stratified) <- gsub("_", " ", colnames(df_dss2_models_stratified))     
colnames(df_dss2_models_stratified) <- sapply(colnames(df_dss2_models_stratified), tools::toTitleCase)
df_dss2_models_stratified <- df_dss2_models_stratified[order(df_dss2_models_stratified$`p Value`),]
rownames(df_dss2_models_stratified) <- NULL

color_tile_custom <- formatter("span", 
                               style = function(x) {
                                 color_scale <- col_numeric(c("lightblue", "white", "lightblue"), domain = range(df_re_intercept_models$`Fixed Effect Coefficient`))
                                 formattable::style(display = "block", padding = "0 4px", `background-color` = color_scale(x))
                               }
)
ft <- formattable(df_dss2_models_stratified,  list(
  `Fixed Effect Term` = formatter("span", style = ~formattable::style(color = "blue")),
  `Fixed Effect Coefficient` = color_tile_custom,
  `p Value` = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  `Significance Stars` = formatter("span", style = ~formattable::style(color = "darkred"))
))

print(ft)
library(htmlwidgets)
htmlwidgets::saveWidget(ft, "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/Table_Stratified_results.html", selfcontained = TRUE)


df_dss2_models_stratified_forest <- df_dss2_models_stratified %>% mutate(`Experimental Variable`= case_when(`Fixed Effect Term` == "Positive Control: BzCl + Nr of Concentration Points: 5 + Drug Sensitivity Readout: CellTiter-Glo" ~ "Positive Control & Dose Points & Readout   ", 
                                                                                               `Fixed Effect Term` ==  "Time until sample usage: Handled within 2-72 hours" ~ "Time Until Sample Usage   ",
                                                                                               `Fixed Effect Term` ==  "Medium: Mononuclear cell medium" ~ "Medium: MCM    ",
                                                                                               `Fixed Effect Term` ==  "Medium: RPMI + fetal bovine serum (FBS) (10%)" ~ "Medium: RPMI & FBS   ",
                                                                                               `Fixed Effect Term` ==  "microenvironmental_stimuliNone" ~ "Microenvironmental Stimuli: None    ",
                                                                                               `Fixed Effect Term` ==  "Microenvironmental stimuli: Transiently cultured with feeder cells, activate samples with autologous BM T helper cells in the presence of IL-2 and a T-cell expansion cocktai" ~ "Microenvironmental Stimuli   ",
                                                                                               `Fixed Effect Term` ==  "Plate readers: Pherastar, Cytation5, Insight, Tecan" ~ "Plate Reader",
                                                                                               `Fixed Effect Term` ==  "Centrifugation procedure: LymphoPrepTM gradient centrifugation" ~ "Centrifugation Procedure: 1*    ",
                                                                                               `Fixed Effect Term` ==  "Centrifugation procedure: Supernatant isolation at 300g 10min, density centrifugation at 400g for 20min without brake, afterwards always 300g 5 min" ~ "Centrifugation Procedure: 2*   ",
                                                                                               `Fixed Effect Term` ==  "Cells: 5000" ~ "Number of Cells per Well   ",
                                                                                               `Fixed Effect Term` ==  "Cell counting method: Trypan blue (Countess™ II FL Automated Cell Counter)" ~ "Cell Counting Method    ",
                                                                                               TRUE ~ `Fixed Effect Term`))



tm <- forest_theme(base_size = 25,
                   refline_col = "red",
                   arrow_type = "closed",
                   footnote_gp = gpar(col = "black", cex = 0.6), 
                   line_size = 10, 
                   align = "center", 
                   footnote_parse = FALSE)
df_dss2_models_stratified_forest$` ` <- paste(rep(" ", 30), collapse = " ")
df_dss2_models_stratified_forest <- df_dss2_models_stratified_forest[order(df_dss2_models_stratified_forest$`Fixed Effect Coefficient`),]
df_dss2_models_stratified_forest$`p-value` <- df_dss2_models_stratified_forest$`Corrected P Value`
#df_re_intercept_models$se <- (log(df_re_intercept_models$Upper) - log(df_re_intercept_models$`Fixed Effect Coefficient`))/1.96
df_dss2_models_stratified_forest$`Fixed Effect Term`
p <- forest(df_dss2_models_stratified_forest[,c('Experimental Variable', ' ', 'Corrected p Value')],
            est = df_dss2_models_stratified_forest$`Fixed Effect Coefficient`,
            lower = df_dss2_models_stratified_forest$Lower, 
            upper = df_dss2_models_stratified_forest$Upper,
            sizes = 1.0,
            ci_column = 2,
            ref_line = 0,
            #arrow_lab = c("Lower DSS2 than reference group", "Higher DSS2 than refernce group"),
            xlim = c(-2, 2),
            ticks_at = c(-1, -0.5, 0, 0.5, 1),
            footnote = "*1:LymphoPrepTM gradient centrifugation \n*2: Supernatant isolation at 300g 10min, density centrifugation at 400g for 20min without brake, afterwards always 300g 5 min",
            theme = tm)

# Print plot
plot(p)
p_wh <- get_wh(p)
pdf('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/Forest_plot_DSS2_stratified.pdf',width = p_wh[1], height = p_wh[2])
plot(p)
dev.off()


group_counts_sample <- all_response_metrics %>%
  group_by(lab, sample, Tissue, Disease.status) %>%
  dplyr::summarize(count = n())

#Sample Analysis##################################################################################################
##----Fresh vs Frosen & Blood vs Bone marrow & Diagnosis vs relapse----
all_response_metrics$Disease.status <- gsub('remission', 'Remission', all_response_metrics$Disease.status)
all_response_metrics$Disease.status <- gsub('Unknown', NA, all_response_metrics$Disease.status)
all_response_metrics$sample <- gsub('Cryopreserved', 'frozen', all_response_metrics$sample) #Frozen
all_response_metrics$Tissue <- gsub('Leukapheresis', 'Blood', all_response_metrics$Tissue) #Blood

frozen_count <- all_response_metrics %>% group_by(sample, drug) %>% summarize(count = n())
tissue_count <- all_response_metrics %>% group_by(Tissue, drug) %>% summarize(count = n())
disease_count <- all_response_metrics %>% group_by(Disease.status, drug) %>% dplyr::summarize(count = n())


mixed_model_sample_DSS2 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (0 + lab|drug), subset(all_response_metrics))
summary(mixed_model_sample_DSS2)
unique(all_response_metrics$sample)
check_model(mixed_model_sample_DSS2)
icc(mixed_model_sample_DSS2)
summary(mixed_model_sample_DSS2_5)
icc(mixed_model_sample_DSS2_5)

mixed_model_sample_DSS2_1 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (lab|drug), subset(all_response_metrics))
summary(mixed_model_sample_DSS2_1)
check_model(mixed_model_sample_DSS2_1)
mixed_model_sample_DSS2_2 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (1|drug) + (1|lab), subset(all_response_metrics))
summary(mixed_model_sample_DSS2_2)
check_model(mixed_model_sample_DSS2_2)
mixed_model_sample_DSS2_3 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (1|drug) + (1|sample) + (1|Tissue) + (1|Disease.status), subset(all_response_metrics))
check_model(mixed_model_sample_DSS2_3)
#mixed_model_sample_DSS2_3 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (1|drug) + (1|sample) + (1|Tissue) + (1|Disease.status) + (1|lab), subset(all_response_metrics))
mixed_model_sample_DSS2_4 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|drug) + (1|sample) + (1|Tissue) + (1|Disease.status) + (1|lab), subset(all_response_metrics))
mixed_model_sample_DSS2_5 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (sample|drug) + (Tissue|drug) + (Disease.status|drug) + (1|lab), subset(all_response_metrics))
mixed_model_sample_DSS2_6 <- lmer(DSS2_boxcox_sclaed2 ~ sample + Tissue + Disease.status + (1|Patient.num) + (sample|drug) + (Tissue|drug) + (Disease.status|drug), subset(all_response_metrics))

aic_0 <- AIC(mixed_model_sample_DSS2)
r2_0 <- r2(mixed_model_sample_DSS2)
aic_1 <- AIC(mixed_model_sample_DSS2_1)

anova(mixed_model_sample_DSS2, mixed_model_sample_DSS2_1, mixed_model_sample_DSS2_2, mixed_model_sample_DSS2_3, mixed_model_sample_DSS2_4, mixed_model_sample_DSS2_5, mixed_model_sample_DSS2_6, test = "LRT")
model.comparison(mixed_model_sample_DSS2, mixed_model_sample_DSS2_1) #, mixed_model_sample_DSS2_2,mixed_model_sample_DSS2_3, mixed_model_sample_DSS2_4, mixed_model_sample_DSS2_5, mixed_model_sample_DSS2_6)
icc(mixed_model_sample_DSS2)

dss2_models_sample <-c(mixed_model_sample_DSS2_1)
df_dss2_models_sample <- model_results(dss2_models_sample)
df_backup <- df_dss2_models_sample
df_dss2_models_sample <- df_backup
df_dss2_models_sample$Random_Effect_Term <- '(1|Patient.num) + (0 + Lab | Drug)'
colnames(df_dss2_models_sample) <- gsub("_", " ", colnames(df_dss2_models_sample))     
colnames(df_dss2_models_sample) <- sapply(colnames(df_dss2_models_sample), tools::toTitleCase)
rownames(df_dss2_models_sample) <- NULL

df_dss2_models_sample <- df_dss2_models_sample %>% mutate(`Biospecimen Type`= case_when(`Fixed Effect Term` == "TissueBone marrow" ~ "Tissue: Bone Marrow", 
                                                       `Fixed Effect Term` ==  "Disease.statusRefractory" ~ "Disease Status: Refractory",
                                                       `Fixed Effect Term` ==  "samplefrozen" ~ "Sample: Frozen",
                                                       `Fixed Effect Term` ==  "Disease.statusRelapse" ~ "Disease Status: Relapse",
                                                       `Fixed Effect Term` ==  "Disease.statusRemission" ~ "Disease Status: Remission",
                                                       TRUE ~ `Fixed Effect Term`))



get_reference_levels <- function(model) {
  # Get model frame (data used in model)
  mf <- model.frame(model)
  
  # Get the fixed effect variable names (exclude intercept)
  fixed_terms <- attr(terms(model), "term.labels")
  
  # Initialize list for reference levels
  ref_levels <- list()
  
  for (term in fixed_terms) {
    # Handle interaction terms separately
    if (grepl(":", term)) next
    
    # Try to get variable from model frame
    var <- mf[[term]]
    
    # Only care about character or factor vars (treated as categorical)
    if (is.character(var) || is.factor(var)) {
      ref_levels[[term]] <- levels(factor(var))[1]
    }
  }
  
  return(ref_levels)
}
get_reference_levels(mixed_model_sample_DSS2)
df_dss2_models_sample <- df_dss2_models_sample %>% mutate(`Reference Group` = case_when(`Fixed Effect Term` == "samplefrozen" ~ "Sample: Fresh", 
                                                                                         str_detect(`Fixed Effect Term`, "Disease.status") ~ "Disease Starus: Diagnosis",
                                                                                         `Fixed Effect Term` == "TissueBone marrow" ~ "Tissue: Blood", 
                                                                                         TRUE ~ `Fixed Effect Term`))

#df_dss2_models_sample <- df_dss2_models_sample[,order(df_dss2_models_sample$`Fixed Effect Coefficient`)]
formattable(df_dss2_models_sample,  list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))


tm <- forest_theme(base_sixe = 12,
                   arrow_type = "closed",
                   footnote_gp = gpar(col = "black", cex = 0.6), 
                   align = "center", 
                   footnote_parse = FALSE, 
                   line_size = 1.5, 
                   xaxis_gp = gpar(fontsize = 10, fontfamily = "Arial"), 
                   spacing = 2,)
df_dss2_models_sample$` ` <- paste(rep(" ", 20), collapse = " ")
df_dss2_models_sample$`p-value` <- round(df_dss2_models_sample$`p Value`, 4)
#df_dss2_models_sample$`p-value` <- round(df_dss2_models_sample$corrected_p_value, 4)
df_dss2_models_sample <- df_dss2_models_sample %>% mutate(`n` = case_when(`Biospecimen Type` == 'Sample: Frozen' ~ '55-98',
                                                                          `Biospecimen Type` == 'Disease Status: Remission' ~ '2-21',
                                                                          `Biospecimen Type` == 'Disease Status: Relapse' ~ '50-97',
                                                                          `Biospecimen Type` == 'Tissue: Bone Marrow' ~ '270-600',
                                                                          `Biospecimen Type` == 'Disease Status: Refractory' ~ '10-145'
                                                                          ))
#df_dss2_models_sample$`Biospecimen Type` <- df_dss2_models_sample$Biospecimens
#df_re_intercept_models$se <- (log(df_re_intercept_models$Upper) - log(df_re_intercept_models$`Fixed Effect Coefficient`))/1.96
p_sample <- forest(df_dss2_models_sample[,c('Biospecimen Type', 'Reference Group', ' ', 'p-value', 'n')],
            est = df_dss2_models_sample$`Fixed Effect Coefficient`,
            lower = df_dss2_models_sample$Lower, 
            upper = df_dss2_models_sample$Upper,
            sizes = 1.0,
            ci_column = 3,
            ref_line = 0,
            #arrow_lab = c("Lower DSS2 than reference group", "Higher DSS2 than refernce group"),
            xlim = c(-1, 1),
            #ticks_at = c(-1, -0.5, 0, 0.5, 1),
            #footnote = "*1:LymphoPrepTM gradient centrifugation \n*2: Supernatant isolation at 300g 10min, density centrifugation at 400g for 20min without brake, afterwards always 300g 5 min",
            theme = tm, 
            xlab = expression('Change in Scaled DSS'[2]),
            font.label = list(size = 10, family = "Arial"),
            font.ticks = list(size = 9),)

# Print plot
plot(p_sample)
p_wh <- get_wh(p_sample)
pdf('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Forest_plot_DSS2_biospecimens.pdf',width = p_wh[1], height = 6)
plot(p_sample)
dev.off()

output_width_cm <- 12 #12 output_width_cm <- 20 #17
output_height_cm <- 10 #10 output_height_cm <- 16 #13

dpi <- 300  # Resolution in dots per inch

# Convert cm to inches (1 inch = 2.54 cm)
output_width_in <- output_width_cm / 2.54
output_height_in <- output_height_cm / 2.54

# Convert inches to pixels for the PNG device
output_width_px <- output_width_in * dpi
output_height_px <- output_height_in * dpi

# Calculate scaling factors for the gtable
scale_width <- output_width_cm / 10  # Base width adjustment (10 cm as reference)
scale_height <- output_height_cm / 10  # Base height adjustment (10 cm as reference)

p_sample <- edit_plot(p_sample, gp = gpar(cex=0.9, fontfamily="Arial")) #1.4
p_sample <- edit_plot(p_sample, part = "header", gp = gpar(cex=0.9, fontfamily="Arial"))

p_sample <- edit_plot(p_sample, col = 4:5, part = "header",
                           which = "text",
                           hjust = unit(0.5, "npc"),
                           x = unit(0.5, "npc"))
p_sample <- edit_plot(p_sample, col = 4:5, part = "body",
                           which = "text",
                           hjust = unit(0.5, "npc"),
                           x = unit(0.5, "npc"))

# Scale the gtable layout
scaled_p_sample_all_metric <- gtable::gtable_filter(p_sample, pattern = ".*", trim = TRUE)  # Keep all grobs
scaled_p_sample_all_metric$widths <- scaled_p_sample_all_metric$widths * scale_width
scaled_p_sample_all_metric$heights <- scaled_p_sample_all_metric$heights * scale_height


png('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Forest_plot_blood.png', height = 16, width = 26, unit = "cm", res = 300) #24.8 height = 26, width = 56 height = 10, width = 30
#10 28
grid.newpage()
grid.draw(scaled_p_sample_all_metric)
dev.off()








model <- mixed_model_sample_DSS2 
df <- data.frame(
  model = character(),
  Fixed_Effect_Term = character(),
  Random_Eeffect_Term = character(),
  Response_metric = character(),
  Fixed_Effect_Coefficient = numeric(),
  p_value = numeric(),
  Intercept = numeric(),
  Degree_of_Freedom = numeric(),
  stringsAsFactors = FALSE, 
  lower = numeric(),
  upper = numeric()
)

  print(model)
  model_formula <- formula(model)
  dependent_var <- as.character(model_formula[[2]])
  fixed_effects <- fixef(model)  
  fixed_effect_terms <- names(fixed_effects)
  fixed_effect_names <- fixed_effect_terms[fixed_effect_terms != "(Intercept)"]
  fixed_effect_coef <- fixed_effects[fixed_effect_names]
  coef_model <- summary(model)$coefficients
  p_values <- coef_model[, "Pr(>|t|)"]
  p_values_excluding_intercept <- p_values[fixed_effect_terms != "(Intercept)"]
  intercept <- fixed_effects[fixed_effect_terms[fixed_effect_terms == "(Intercept)"]]
  degree_freedom <- coef_model[, "df"]
  degree_freedom <- degree_freedom[fixed_effect_terms != "(Intercept)"]
  print(degree_freedom)
  response_metric <- dependent_var
  #c <- confint(model, method = "profile")
  #ci_fixed <- c[rownames(c) %in% fixed_effect_names, , drop = FALSE]  # Subset only fixed effects
  #lower <- ci_fixed[, 1]
  #upper <- ci_fixed[, 2]
  
  
  #Create a temporary data frame for this model's fixed effects
  temp_df <- data.frame(
    model = paste0(model_formula[2], model_formula[1], model_formula[3]),
    Fixed_Effect_Term = fixed_effect_names,
    Random_Effect_Term = str_remove(strsplit(as.character(model_formula[3]), split = "\\+ \\(")[[1]][2], "\\)"),
    Response_metric = response_metric,
    Fixed_Effect_Coefficient = fixed_effect_coef,
    p_value = p_values_excluding_intercept,
    Intercept = intercept,
    Degree_of_Freedom = degree_freedom,
    stringsAsFactors = FALSE 
  )
  df <- rbind(df, temp_df)


get_significance_stars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return(" ")
  }
}

# Check if `df` has the necessary columns
if (!all(c("Response_metric", "model") %in% colnames(df))) {
  stop("The columns 'Response_metric' and/or 'model' are missing in the data frame.")
}

# Check the structure of `df` for potential issues
print(str(df))

# Check for missing values in critical columns
if (any(is.na(df$Response_metric)) || any(is.na(df$model))) {
  warning("There are NA values in 'Response_metric' or 'model' columns.")
}

print("######################STARS#######################")
#Add significance stars to the data frame
df$significance_stars <- sapply(df$p_value, get_significance_stars)
#Count the number of fixed effects per category
fixed_effect_counts <- df %>%
  group_by(Response_metric, model) %>%
  filter(row_number() == 1) %>%   
  ungroup() %>%                   
  group_by(Response_metric) %>%           
  dplyr::summarize(count = n()) %>%
  as.data.frame()
print(df)
df <- df %>%
  left_join(fixed_effect_counts, by = "Response_metric")
print("####################LEFT")
df$corrected_p_value <- df$p_value * df$count


print(df$model)

rownames(df) <- NULL
df <- df %>% mutate(Fixed_Effect_Term = case_when(Fixed_Effect_Term ==  "samplefrozen" ~ "Sample: Frozen", 
                                                  Fixed_Effect_Term ==  "TissueBone marrow" ~ "Tissue: Bone marrow",
                                                  Fixed_Effect_Term ==  "Disease.statusRefractory" ~ "Disease Status: Refractory",
                                                  Fixed_Effect_Term ==  "Disease.statusRelapse" ~ "Disease Status: Relapse",
                                                  Fixed_Effect_Term ==  "Disease.statusRemission" ~ "Disease Status: Remission",
                                                  TRUE ~ Fixed_Effect_Term))

df <- df %>% mutate(Response_metric = case_when(Response_metric ==  "DSS2_boxcox_sclaed2" ~ "DSS2: Box Cox Transformed and Scaled",
                                                Response_metric == "auc_a_box_cox_scaled2" ~ "AUC A: Box cox Transformed and Scaled",
                                                TRUE ~ Response_metric))
df <- df[order(-df$Fixed_Effect_Coefficient), ]


##----ScatterPlots Fresh Vs Frozen----

all_response_metrics$sample
result <- all_response_metrics %>%
  group_by(Patient.num) %>%
  summarise(unique_x_count = n_distinct(sample)) %>%
  filter(unique_x_count == 2)

# Count how many unique ids meet the criteria
count_result <- nrow(result)

# Print the result
print(result)       # Shows the ids with exactly two unique instances of x
print(count_result) # Shows the count of such ids
print(unique(all_response_metrics$sample))

df1_filtered <- all_response_metrics[all_response_metrics$Patient.num %in% result$Patient.num, ]

df1_filtered$DSS2_pos <- df1_filtered$DSS2 - min(df1_filtered$DSS2) + 0.1
MASS::boxcox(DSS2_pos ~ sample, 
             data = df1_filtered,
             lambda = seq(-0.25, 2, length.out = 10))

df1_filtered$DSS2_boxcox <- -((df1_filtered$DSS2)^0.5 - 1) / 0.5
#df1_filtered$DSS2_boxcox <- log10(df1_filtered$DSS2)


for(j in unique(df1_filtered$drug)){
  df1_filtered$DSS2_box_cox_sclaed2[df1_filtered$drug==j] <-
    scale(df1_filtered$DSS2_boxcox[df1_filtered$drug==j])
}

ggplot(enserink_fresh_and_frozen_patients, aes(x=DSS2_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(enserink_fresh_and_frozen_patients, aes(x=DSS2_box_cox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(enserink_fresh_and_frozen_patients, aes(x=DSS2)) +
  geom_density(color="darkblue", fill="lightblue")

new_df <- df1_filtered %>%
  group_by(Patient.num, drug) %>%
  reframe(
    Fresh = DSS2[sample == 'fresh'],
    Frozen = DSS2[sample == 'frozen']
  )

new_df$Patient.num <- str_replace_all(new_df$Patient.num, '_', ' ')

# View the new data frame
print(new_df)

library(ggpubr)

# Perform a linear model
lm_model <- lm(Frozen ~ Fresh, data = new_df)
lmm_model <- lmer(Fresh ~ Frozen + (1|drug), data = new_df)
visualize(lmm_model)
summary(lmm_model)
check_model(lmm_model)

check_model(lm_model)

# Calculate adjusted R^2
adjusted_r2 <- summary(lm_model)$adj.r.squared

# Extract the p-value for Fresh
pval <- summary(lm_model)$coefficients[2, 4] 

# Apply FDR adjustment (assuming multiple tests, example 10)
pval_fdr <- p.adjust(pval, method = "fdr")

pearson_cor_test <- cor.test(new_df$Fresh, new_df$Frozen, method = "pearson")
pearson_cor <- pearson_cor_test$estimate
p_value <- pearson_cor_test$p.value

# Prepare the annotation text
annotation_text <- paste0(
  "R = ", round(pearson_cor, 3), 
  "\np ", ifelse(format(p_value, scientific = TRUE, digits = 3)<0.001,paste("=",format(p_value, scientific = TRUE, digits = 3)), "<0.001")
)

custom_colors <- c(
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", 
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", 
  "#bc80bd", "#ccebc5"  # Added two complementary colors
)
display.brewer.pal(n = 10, name = "Set3")
# Create the plot
scatterplot_fresh_vs_frozen_karolinska <- ggplot(new_df, aes(x = Fresh, y = Frozen)) +
  geom_point(aes(color = factor(Patient.num)), size = 3) +
  ylim(0, 45) + 
  xlim(0,45) +
  coord_fixed(ratio = 1)+
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  scale_color_manual(values = custom_colors) +
  annotate(
    "text", 
    x = 0, y = 45, 
    label = annotation_text, 
    size = 3.51, family = "Arial", hjust = 0, vjust = 1
  ) +
  labs(color = "Patient ID", x= expression("Fresh DSS"[2]), y = expression("Frozen DSS"[2])) +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 10),    # All text Arial, size 10
    axis.text = element_text(size = 8, family = "Arial"), # Tick labels Arial, size 8
    legend.text = element_text(size = 10, family = "Arial"), # Legend text Arial, size 10
    legend.title = element_text(size = 10, family = "Arial"), # Legend title Arial, size 10
    legend.background = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(), # Optional: Remove background panel
    plot.background = element_blank()   # Optional: Remove plot background
  )

ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/Karolinska_fresh_vs_frozen.png", plot = scatterplot_fresh_vs_frozen_karolinska, width = 10.76, height = 8.09, units="cm") #width = 20.76, height = 7.09


df1_filtered <- df1_filtered %>% mutate(sample_bin = ifelse(df1_filtered$sample == 'fresh', 1, 0))

wilco_test <- wilcox_test(df1_filtered, 'drug', 'sample_bin', 'DSS2')
create_volcano_plot(wilco_test, "", y_axis = "p_value", "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/vulcano_fresh_vs_frozen_karolinska.png")

# For each drug, calculate the Pearson correlation between cell_line1 and cell_line2
cor_results <- lapply(unique(new_df$drug), function(drug) {
  # Subset data for each drug
  drug_data <- subset(new_df, drug == drug)
  
  cor_val <- cor(new_df$Fresh, new_df$Frozen)
  print(paste("Correlation for", drug, ":", cor_val))
  print(summary(new_df$Fresh))
  print(summary(new_df$Frozen))
  
  # Calculate Pearson correlation between cell_line1 and cell_line2
  cor_test <- cor.test(new_df$Fresh, new_df$Frozen, method = "spearman")
  
  # Return the correlation coefficient and p-value for this drug
  data.frame(drug = drug, Correlation = cor_test$estimate, P_Value = cor_test$p.value)
})

# Combine the results into a single data frame
cor_results_df <- do.call(rbind, cor_results)

# Display the correlation results for each drug
print(cor_results_df)

##----Enserink scatterplot fresh vs frozen ----


library(readr)
library("ggplot2")
library(dplyr)
library(readxl)
library(openxlsx)
library(stringr) # for string operations

enserink_repeated_frozen <- read_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/DSS_score/DSS/Original_DSS2_2025-02-24.xlsx')

enserink_pubchem_drug_names <- read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_drug_information_pubchem.csv')

enserink_repeated_frozen <- merge(enserink_repeated_frozen, enserink_pubchem_drug_names, by.x = "DRUG_NAME", by.y = "org_drug_name", all.y = TRUE)
enserink_repeated_frozen$...1 <- NULL

enserink_full_dss_set <- read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Second cleansing/dss_enserink_common_drugs.csv')
enserink_full_dss_set$sample

enserink_repeated_frozen_long <- enserink_repeated_frozen %>%
  pivot_longer(cols = -c('DRUG_NAME', 'ID', 'pubchem_drug_name'),  # Select all columns except the identifier column (e.g., "col1")
               names_to = "Patient.num",  # New column to store the original column names
               values_to = "DSS2")
unique(enserink_repeated_frozen_long$Patient.num)

enserink_repeated_frozen_long$sample <- 'frozen'

enserink_repeated_frozen_long <- enserink_repeated_frozen_long %>%
  mutate(Patient.num = ifelse(Patient.num == "Pt34r_frozen", "Patient 34_relapse", Patient.num))
enserink_repeated_frozen_long <- enserink_repeated_frozen_long %>%
  mutate(Patient.num = ifelse(Patient.num == "Pt30r_repeat", "Patient 30_relapse", Patient.num))
enserink_repeated_frozen_long <- enserink_repeated_frozen_long %>%
  mutate(Patient.num = ifelse(Patient.num == "Pt33_frozen", "Patient 33", Patient.num))

enserink_repeated_frozen_long$drug <- enserink_repeated_frozen_long$pubchem_drug_name
enserink_repeated_frozen_long_subset <- subset(enserink_repeated_frozen_long, select = c("drug", "DSS2", "Patient.num", "sample"))

names(enserink_full_dss_set)
enserink_full_dss_set_subset <- subset(enserink_full_dss_set, Patient.num == 'Patient 34_relapse' | Patient.num == 'Patient 30_relapse' | Patient.num == 'Patient 33', select = c("drug", "DSS2", "Patient.num", "sample"))
enserink_fresh_and_frozen_patients <- rbind(enserink_repeated_frozen_long_subset, enserink_full_dss_set_subset)


enserink_fresh_and_frozen_patients$DSS2[enserink_fresh_and_frozen_patients$DSS2 < 0] 

min(enserink_fresh_and_frozen_patients$DSS2)
enserink_fresh_and_frozen_patients$DSS2_pos <- enserink_fresh_and_frozen_patients$DSS2 - min(enserink_fresh_and_frozen_patients$DSS2) + 0.1
MASS::boxcox(DSS2_pos ~ sample, 
             data = enserink_fresh_and_frozen_patients,
             lambda = seq(-0.25, 2, length.out = 10))

enserink_fresh_and_frozen_patients$DSS2_boxcox <- -((enserink_fresh_and_frozen_patients$DSS2)^0.5 - 1) / 0.5
enserink_fresh_and_frozen_patients$DSS2_boxcox <- log10(enserink_fresh_and_frozen_patients$DSS2)


for(j in unique(enserink_fresh_and_frozen_patients$drug)){
  enserink_fresh_and_frozen_patients$DSS2_box_cox_sclaed2[enserink_fresh_and_frozen_patients$drug==j] <-
    scale(enserink_fresh_and_frozen_patients$DSS2_boxcox[enserink_fresh_and_frozen_patients$drug==j])
}

ggplot(enserink_fresh_and_frozen_patients, aes(x=DSS2_boxcox)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(enserink_fresh_and_frozen_patients, aes(x=DSS2_box_cox_sclaed2)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(enserink_fresh_and_frozen_patients, aes(x=DSS2)) +
  geom_density(color="darkblue", fill="lightblue")



enserink_fresh_vs_frozen<- enserink_fresh_and_frozen_patients %>%
  group_by(Patient.num, drug) %>%
  dplyr::summarise(
    Fresh = DSS2_boxcox[sample == 'fresh'],
    Frozen = DSS2_boxcox[sample == 'frozen'],
    .groups = "drop"
  )


enserink_fresh_vs_frozen<- enserink_fresh_and_frozen_patients %>%
  group_by(Patient.num, drug) %>% 
  reframe(Fresh = DSS2[sample == 'fresh'],
          Frozen = DSS2[sample == 'frozen'])

enserink_fresh_vs_frozen$Patient.num <- str_replace_all(enserink_fresh_vs_frozen$Patient.num, "_", " ")
enserink_fresh_vs_frozen <- enserink_fresh_vs_frozen %>%
  mutate(Patient.num = ifelse(Patient.num == "Patient 33", "Patient 33 diagnosis", Patient.num))

library(ggpubr)
library(p.adjust)

# Example data
new_df <- data.frame(
  Fresh = rnorm(100, 10, 5),
  Frozen = rnorm(100, 20, 10),
  Patient.num = sample(1:10, 100, replace = TRUE)
)

# Perform a linear model
lm_model <- lm(Frozen ~ Fresh, data = enserink_fresh_vs_frozen)

check_model <- lmer(Frozen ~ Fresh, data = enserink_fresh_vs_frozen)
visualize(check_model)
summary(check_model)
check_model(check_model)
library(MuMIn)
r.squaredGLMM(check_model)

# Calculate adjusted R^2
adjusted_r2 <- summary(lm_model)$adj.r.squared

# Extract the p-value for Fresh
pval <- summary(lm_model)$coefficients[2, 4] 

# Apply FDR adjustment (assuming multiple tests, example 10)
pval_fdr <- p.adjust(pval, method = "fdr")

pearson_cor_test <- cor.test(enserink_fresh_vs_frozen$Fresh, enserink_fresh_vs_frozen$Frozen, method = "pearson")
pearson_cor <- pearson_cor_test$estimate
p_value <- pearson_cor_test$p.value

# Prepare the annotation text
annotation_text <- paste0(
  "R = ", round(pearson_cor, 3), 
  "\np ", ifelse(format(p_value, scientific = TRUE, digits = 3)<0.001,paste("=",format(p_value, scientific = TRUE, digits = 3)), "<0.001")
)

custom_colors <- c(
  "#8dd3c7", "#bebada", "#fb8072", 
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", 
  "#bc80bd", "#ccebc5"  # Added two complementary colors
)
display.brewer.pal(n = 10, name = "Set3")
# Create the plot
scatterplot_fresh_vs_frozen_enserink <- ggplot(enserink_fresh_vs_frozen, aes(x = Fresh, y = Frozen)) +
  geom_point(aes(color = factor(Patient.num)), size = 3) +
  ylim(0, 45) +
  xlim(0,45) +
  coord_fixed(ratio = 1)+
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  scale_color_manual(values = custom_colors) +
  annotate(
    "text", 
    x = 0, y = 45, 
    label = annotation_text, 
    size = 3.51, family = "Arial", hjust = 0, vjust = 1
  ) +
  labs(color = "Patient ID", x= expression("Fresh DSS"[2]), y = expression("Frozen DSS"[2])) +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 10),    # All text Arial, size 10
    axis.text = element_text(size = 8, family = "Arial"), # Tick labels Arial, size 8
    legend.text = element_text(size = 10, family = "Arial"), # Legend text Arial, size 10
    legend.title = element_text(size = 10, family = "Arial"), # Legend title Arial, size 10
    legend.background = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(), # Optional: Remove background panel
    plot.background = element_blank(),   # Optional: Remove plot background
    legend.position = c(0.85, 0.25)
  )
scatterplot_fresh_vs_frozen_enserink
ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/scatterplot_fresh_vs_frozen_enserink.png", plot = scatterplot_fresh_vs_frozen_enserink, width = 10.76, height = 8.09, units="cm") #width = 20.76, height = 7.09 # width = 10.76, height = 8.09 width = 10.76, height = 8.09, units="cm"


enserink_fresh_and_frozen_patients <- enserink_fresh_and_frozen_patients %>% mutate(sample_bin = ifelse(enserink_fresh_and_frozen_patients$sample == 'fresh', 1, 0))

wilco_test <- wilcox_test(enserink_fresh_and_frozen_patients, 'drug', 'sample_bin', 'DSS2')
create_volcano_plot(wilco_test, "", y_axis = "p_value", "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/vulcano_fresh_vs_frozen_enserink.png")

#result_df, title, filename, threshold_pvalue = 0.05, threshold_logFC = 1, x_lab = "Change in DSS", x_axis = "FC", y_axis = "p_value", y_lab = "-log10(p-adj-value)"


t_test_result <- t.test(enserink_fresh_vs_frozen$Fresh, enserink_fresh_vs_frozen$Frozen, paired = TRUE)

# Print results
print(t_test_result)

# Pearson correlation
pearson_cor <- cor.test(enserink_fresh_vs_frozen$Fresh, enserink_fresh_vs_frozen$Frozen, method = "pearson")

# Spearman correlation (if data is not normally distributed)
spearman_cor <- cor.test(enserink_fresh_vs_frozen$Fresh, enserink_fresh_vs_frozen$Frozen, method = "spearman")

# Print results
print(pearson_cor)
print(spearman_cor)


ggplot(enserink_fresh_vs_frozen, aes(x = Fresh, y = Frozen)) +
  geom_point(color = "#1f78b4", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "Correlation Between Cell Line 1 and Cell Line 2",
       x = "IC50 (Cell Line 1)",
       y = "IC50 (Cell Line 2)") +
  theme_minimal()


##----Diagnosis vs relapse ----
all_response_metrics[all_response_metrics$lab == "Beat AML","Patient_ID"]
result <- all_response_metrics %>%
  group_by(Patient_ID) %>%
  filter(Disease.status == 'Refractory' | Disease.status == 'Diagnosis') %>%
  summarise(unique_x_count = n_distinct(Disease.status)) %>%
  filter(unique_x_count > 1)

all_response_metrics %>%
  group_by(Patient_ID) %>%
  dplyr::summarize(c = n()) 

df1_filtered <- all_response_metrics[all_response_metrics$Patient_ID %in% result$Patient_ID, ]

df1_filtered$DSS2_pos <- df1_filtered$DSS2 - min(df1_filtered$DSS2) + 0.1
MASS::boxcox(DSS2_pos ~ Disease.status, 
             data = df1_filtered,
             lambda = seq(-0.25, 2, length.out = 10))

df1_filtered$DSS2_boxcox <- -((df1_filtered$DSS2)^0.5 - 1) / 0.5
#df1_filtered$DSS2_boxcox <- log10(df1_filtered$DSS2)


for(j in unique(df1_filtered$drug)){
  df1_filtered$DSS2_box_cox_sclaed2[df1_filtered$drug==j] <-
    scale(df1_filtered$DSS2_boxcox[df1_filtered$drug==j])
}


new_df <- df1_filtered[,c("Patient_ID", "drug", "Disease.status", "DSS2")] %>%
  group_by(Patient_ID, drug, Disease.status) %>%
  summarise(DSS2 = mean(DSS2), .groups = "drop") %>%  #
  pivot_wider(
    names_from = Disease.status,
    values_from = DSS2
  )

new_df$Patient_ID <- str_replace_all(new_df$Patient_ID, '_', ' ')

# View the new data frame
print(new_df)

library(ggpubr)


# Perform a linear model
lm_model <- lm(Diagnosis ~ Refractory, data = new_df)


check_model(lm_model)

# Calculate adjusted R^2
adjusted_r2 <- summary(lm_model)$adj.r.squared

# Extract the p-value for Fresh
pval <- summary(lm_model)$coefficients[2, 4] 

# Apply FDR adjustment (assuming multiple tests, example 10)
pval_fdr <- p.adjust(pval, method = "fdr")

pearson_cor_test <- cor.test(new_df$Diagnosis, new_df$Refractory, method = "pearson")
pearson_cor <- pearson_cor_test$estimate
p_value <- pearson_cor_test$p.value

# Prepare the annotation text
annotation_text <- paste0(
  "R = ", round(pearson_cor, 3), 
  "\np = ", ifelse(format(p_value, scientific = TRUE, digits = 3) < 0.001, format(p_value, scientific = TRUE, digits = 3), '<0.001')
)

custom_colors <- c(
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", 
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", 
  "#bc80bd", "#ccebc5"  # Added two complementary colors
)
colors <- c(
  brewer.pal(12, "Set3"),
  brewer.pal(3, "Pastel1")[1:2]
)
display.brewer.pal(n = 14, name = "Set3")
# Create the plot
scatterplot_diagnosis_vs_refactory <- ggplot(new_df, aes(x = Diagnosis, y = Refractory)) +
  geom_point(aes(color = factor(Patient_ID)), size = 3) +
  ylim(0, 45) + 
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  scale_color_manual(values = colors) +
  annotate(
    "text", 
    x = 0, y = 45, 
    label = annotation_text, 
    size = 3.51, family = "Arial", hjust = 0, vjust = 1
  ) +
  labs(color = "Patient ID") +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 10),    # All text Arial, size 10
    axis.text = element_text(size = 8, family = "Arial"), # Tick labels Arial, size 8
    legend.text = element_text(size = 10, family = "Arial"), # Legend text Arial, size 10
    legend.title = element_text(size = 10, family = "Arial"), # Legend title Arial, size 10
    legend.background = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(), # Optional: Remove background panel
    plot.background = element_blank()   # Optional: Remove plot background
  )

ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/diagnosis_vs_refactory.png", plot = scatterplot_diagnosis_vs_refactory, width = 10.76, height = 8.09, units="cm") #width = 20.76, height = 7.09
