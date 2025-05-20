# Install libraries ----
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


source('~/Desktop/UiO/Project 1/code/plot_functions_project1.R')

#Import Fimm datasets ----
fimm_drug_response <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_drug_response.csv")
fimm_drug_response$...1 <- NULL
View(fimm_drug_response)
fimm_drug_response <- as.data.frame(fimm_drug_response)
fimm_drug_reponse_long <- gather(fimm_drug_response, p_id, dss_score, 'AML_084_04':'Healthy_17')


fimm_sample_annotation <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_sample_annotation.csv")
fimm_sample_annotation$...1 <- NULL
View(fimm_sample_annotation)

fimm_clinical_summary <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_clinical_summary.csv")
fimm_clinical_summary$...1 <- NULL
View(fimm_clinical_summary)

fimm_drug_library <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_drug_library.csv")
fimm_drug_library$...1 <- NULL
View(fimm_drug_library)

fimm_assay_dets <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_assay_dets.csv")
fimm_assay_dets$...1 <- NULL
View(fimm_assay_dets)

fimm_binary_mutations <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_binary_mutations.csv")
fimm_binary_mutations$...1 <- NULL
View(fimm_binary_mutations)

fimm_mutation_variant_allele_frequency <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_mutation_variant_allele_frequency.csv")
fimm_mutation_variant_allele_frequency$...1 <- NULL
View(fimm_mutation_variant_allele_frequency)

fimm_rna_raw_reads <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_rna_raw_reads.csv")
fimm_rna_raw_reads$...1 <- NULL
View(fimm_rna_raw_reads)

fimm_rna_seq_healthy <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_rna_seq_healthy.csv")
fimm_rna_seq_healthy$...1 <- NULL
View(fimm_rna_seq_healthy)

fimm_seq_sample_annotation <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_seq_sample_annotation.csv")
fimm_seq_sample_annotation$...1 <- NULL
View(fimm_seq_sample_annotation)

fimm_drug_library_pubchem <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_drug_library_pubchem.csv")
fimm_drug_library_pubchem$...1 <- NULL
fimm_drug_library_pubchem$dataset <- 'FIMM'
View(fimm_drug_library_pubchem)

#import beat aml datasets ----
beat_aml_inhibitor <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_inhibitor.csv")
beat_aml_inhibitor$...1 <- NULL
View(beat_aml_inhibitor)

beat_aml_calc <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_calc.csv")
beat_aml_calc$...1 <- NULL
View(beat_aml_calc)

beat_aml_mutations <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_mutations.csv")
beat_aml_mutations$...1 <- NULL
View(beat_aml_mutations)

#protein
beat_aml_waves_norm <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_waves_norm.csv")
beat_aml_waves_norm$...1 <- NULL
View(beat_aml_waves_norm)

beat_aml_clinical <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_clinical.csv")
beat_aml_clinical$...1 <- NULL
View(beat_aml_clinical)

beat_aml_vg_cts <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_vg_cts.csv")
beat_aml_vg_cts$...1 <- NULL
View(beat_aml_vg_cts)

beataml_drug_information_pubchem <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/beataml_drug_information_pubchem.csv")
beataml_drug_information_pubchem$...1 <- NULL
beataml_drug_information_pubchem$dataset <- 'BeatAML'
View(beataml_drug_information_pubchem)

beataml_dss_github <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_beat_aml.csv")
beataml_dss_github$...1 <- NULL
beataml_dss_github$dataset <- 'BeatAML'
View(beataml_dss_github)
rows_with_na <- beataml_dss_github[is.na(beataml_dss_github$DSS2), ]

print(subset(beataml_dss_github, beataml_dss_github$drug == '001, RAD'))

beataml_dss_github_wide <- pivot_wider(subset(beataml_dss_github, select=c('drug', 'DSS2', 'Patient.num')),names_from = drug, values_from = DSS2)
beataml_dss_github_wide <- as.data.frame(beataml_dss_github_wide)
rownames(beataml_dss_github_wide) <- beataml_dss_github_wide$Patient.num
beataml_dss_github_wide$Patient.num <- NULL


#Enserink lab datasets ---
enserink_lab_drug_information_pubchem <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_drug_information_pubchem.csv")
enserink_lab_drug_information_pubchem$...1 <- NULL
enserink_lab_drug_information_pubchem$dataset <- 'Enserink_lab'
View(enserink_lab_drug_information_pubchem)

enserink_lab_drug_sensitivity <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_drug_sensitivity.csv")
enserink_lab_drug_sensitivity$...1 <- NULL
enserink_lab_drug_sensitivity$dataset <- 'Enserink_lab'
View(enserink_lab_drug_sensitivity)

enserink_lab_drug_sensitivity_wide <- pivot_wider(subset(enserink_lab_drug_sensitivity, select=c('drug', 'DSS2', 'Patient.ID')),names_from = drug, values_from = DSS2)
enserink_lab_drug_sensitivity_wide <- as.data.frame(enserink_lab_drug_sensitivity_wide)
rownames(enserink_lab_drug_sensitivity_wide) <- enserink_lab_drug_sensitivity_wide$Patient.ID
enserink_lab_drug_sensitivity_wide$Patient.ID <- NULL

#Karolinska
karolinska_responses <- read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/karolinska_normalised_reponse_raw.csv')


#Import common drug names according to pubchem ----
common_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/common_drugnames_pubchem.csv")
common_drugs$`Unnamed: 0_x` <- NULL
common_drugs$`Unnamed: 0_y` <- NULL
common_drugs$`Unnamed: 0` <- NULL
View(common_drugs)


#Common DSS values 3 datasets---
common_enserink_dss <- inner_join(enserink_lab_drug_sensitivity, common_drugs, by=join_by('drug'=='enserink_drug_name'),keep = TRUE)


# Beat AML data exploration ----
report_loc = '/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/R'
data_exploration(beat_aml_calc, output_loc = report_loc, file_name = 'beat_aml_calc')
data_exploration(beat_aml_inhibitor, output_loc = report_loc, file_name = 'beat_aml_inhibitor')
data_exploration(beat_aml_mutations, output_loc = report_loc, file_name = 'beat_aml_mutations')
data_exploration(beat_aml_clinical, output_loc = report_loc, file_name = 'beat_aml_clinical')

# FIMM data exploration ----
data_exploration(fimm_drug_response, output_loc = report_loc, file_name = 'fimm_drug_response')
data_exploration(fimm_sample_annotation, output_loc = report_loc, file_name = 'fimm_sample_annotation')
data_exploration(fimm_clinical_summary, output_loc = report_loc, file_name = 'fimm_clinical_summary')
data_exploration(fimm_drug_library, output_loc = report_loc, file_name = 'fimm_drug_library')
data_exploration(fimm_assay_dets, output_loc = report_loc, file_name = 'fimm_assay_dets')
data_exploration(fimm_binary_mutations, output_loc = report_loc, file_name = 'fimm_binary_mutations')

# Karolinska data exploration ----
data_exploration(karolinska_responses, output_loc = report_loc, file_name = 'karolinska_drug_response')
unique(karolinska_responses$sample)

# Create a wide dataframe to show presence of IDs in col1 values
df_wide <- unique(subset(karolinska_responses, select = c(Patient.num, sample))) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = sample, values_from = value, values_fill = 0)

# Convert to long format for ggplot
df_long <- df_wide %>%
  pivot_longer(cols = -Patient.num, names_to = "sample", values_to = "value", )

totals <- df_long %>%
  group_by(sample) %>%
  summarise(total = sum(as.numeric(as.character(value))))

# Calculate combined totals for each p_id
combined_totals <- df_wide %>%
  mutate(both = fresh & frozen, either = fresh |frozen) %>%
  summarise(
    fresh = sum(fresh),
    frozen = sum(frozen),
    both = sum(both),
    either = sum(either)
  ) %>%
  pivot_longer(cols = everything(), names_to = "sample", values_to = "value") %>% as.data.frame()


p <- ggplot(df_long, aes(x = sample, y = factor(Patient.num), fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("1" = "red", "0" = "white")) +
  labs(x = "sample", y = "Patient-num", fill = "value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
print(p)



#Venn diagram overlapping drugs
beat_aml_calc$drug_name <- sub(" \\(.*", "", beat_aml_calc[["inhibitor"]])

beat_aml_drug <- beat_aml_calc[["drug_name"]]
fimm_drug <- fimm_drug_library[["Preferred name"]]
intersection <- intersect(beat_aml_drug, fimm_drug)

venn_diagram(beat_aml_drug, fimm_drug, "drug_name", "Preferred name", output_plot1 = "/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/overlapping_drugs.png")

beat_aml_calc$dbgap_subject_id <- as.character(beat_aml_calc$dbgap_subject_id)

for (drug in intersection) {
  subset_df <- subset(beat_aml_calc, drug_name == drug)
  print(subset_df)
  box_plots(subset_df, 'dbgap_subject_id', "auc", fill = "drug_name", col3 = NaN, col4 = NaN, filename = paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/", drug, ".png"))
}

for (drug in intersection) {
  subset_df <- subset(fimm_drug_reponse_long, Drug_name == drug)
  print(subset_df)
  box_plots(subset_df, 'dbgap_subject_id', "dss_score", fill = "Drug_name", col3 = NaN, col4 = NaN, filename = paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/fimm_", drug, ".png"))
}


box_plots(beataml_dss_github, 'drug', 'DSS2', fill = "dataset", col3 = "Patient.num", col4 = NaN, filename = paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V1/beataml_box_plot_dss2_github.png"))
box_plots(beataml_dss_github, 'Patient.num', 'DSS2', fill = "dataset", col3 = "drug", col4 = NaN, filename = paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V1/beataml_box_plot_dss2_github_patients.png"))

box_plots(enserink_lab_drug_sensitivity, 'drug', 'DSS2', fill = "dataset", col3 = "Patient.ID", col4 = NaN, filename = paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V1/eserink_DSS2_box_plot.png"))
box_plots(enserink_lab_drug_sensitivity, 'Patient.ID', 'DSS2', fill = "dataset", col3 = "drug", col4 = NaN, filename = paste("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V1/eserink_DSS2_box_plot_patients.png"))



fimm_drug_response_heatmap <- as.data.frame(fimm_drug_response)
rownames(fimm_drug_response_heatmap) <- fimm_drug_response_heatmap$Drug_name
fimm_drug_response_heatmap$Drug_name <- NULL
fimm_drug_response_heatmap$Drug_ID <- NULL
heatmap_plots(fimm_drug_response_heatmap, a_col = NULL, filename = '~/Desktop/UiO/Project 1/Figures/V1/fimm_drug_response_heatmap.png')


fimm_rna_seq_healthy_t <- as.data.frame(fimm_rna_seq_healthy)
rownames(fimm_rna_seq_healthy_t) <- fimm_rna_seq_healthy_t$`Unnamed: 0` 
fimm_rna_seq_healthy_t$`Unnamed: 0` <- NULL
fimm_rna_seq_healthy_pca_t<- t(fimm_rna_seq_healthy_t)
fimm_rna_seq_healthy_pca_t$usunjid <- rownames(fimm_rna_seq_healthy_pca_t)
pca_plots(fimm_rna_seq_healthy_pca_t, file_loc = '~/Desktop/UiO/Project 1/Figures/V1/pca_fimm_rna_seq_healthy_t.png', leave_out = c('usubjid'))

fimm_rna_seq_healthy_pca_t <- as.data.frame(fimm_rna_seq_healthy_pca_t)
fimm_rna_seq_healthy_pca_t <- head(fimm_rna_seq_healthy_pca_t,10) %>% mutate(group = ifelse(grepl("_01$", usubjid), "initial", "relapse"))


beataml_dss_github_wide[is.na(beataml_dss_github_wide)] <- 0

heatmap_plots(t(beataml_dss_github_wide), a_col = NULL, filename = '~/Desktop/UiO/Project 1/Figures/V1/beat_aml_dss2_drug_response_heatmap.png')
beataml_dss_github_wide$usubjid <- as.numeric(rownames(beataml_dss_github_wide))
pca_plots(beataml_dss_github_wide, file_loc = '~/Desktop/UiO/Project 1/Figures/V1/pca_beat_aml_dss2.png', leave_out = c('usubjid'))

enserink_lab_drug_sensitivity_wide$usubjid <- rownames(enserink_lab_drug_sensitivity_wide)
pca_plots(enserink_lab_drug_sensitivity_wide, file_loc = '~/Desktop/UiO/Project 1/Figures/V1/enserink.png', leave_out = c('usubjid'))



result <- beataml_dss_github %>%
  group_by(drug) %>%
  summarise(number_patients_tested = n_distinct(Patient.num)) 


print(result)


beataml_dss_github %>%
  summarise(number_patients_tested = n_distinct(Patient.num)) 

beataml_dss_github %>%
  summarise(number_patients_tested = n_distinct(drug)) 

subset(beataml_dss_github,is.na(beataml_dss_github$Patient.num))

result %>%
  mutate(drug = fct_reorder(drug, desc(number_patients_tested))) %>%
  ggplot( aes(x=drug, y=number_patients_tested)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  geom_text(aes(label=number_patients_tested), vjust=0.5, hjust=-0.5, size=3.5)+
  coord_flip() +
  xlab("") +
  theme_bw() +
  geom_hline(yintercept=569)


common_enserink_dss %>%
  summarise(number_patients_tested = n_distinct(Patient.ID)) 

common_enserink_dss%>%
  group_by(drug) %>%
  summarise(number_patients_tested = n_distinct(Patient.ID)) %>%
  mutate(drug = fct_reorder(drug, desc(number_patients_tested))) %>%
  ggplot( aes(x=drug, y=number_patients_tested)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  geom_text(aes(label=number_patients_tested), vjust=0.5, hjust=-0.5, size=3.5)+
  coord_flip() +
  xlab("") +
  theme_bw() +
  geom_hline(yintercept=69)


fimm_sample_annotation %>% 
  group_by(Disease.status) %>%
  summarise(n_disease_status = n())
