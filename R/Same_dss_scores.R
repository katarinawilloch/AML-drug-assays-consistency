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
library(arsenal)
library(ggplot2)
install.packages("remotes")
remotes::install_github("G-Thomson/gthor")
source('~/Desktop/UiO/Project 1/code/R/plot_functions_project1.R')


#Testing same scores for ENSERINK lab ----
#Getting R script generated scores
dss_github_enserink <- read_csv('~/Desktop/UiO/Project 1/Data/Response scores/dss_github_enserink.csv') #'~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink_full_drug_set_testing.csv'

#Getting precalculated scores from enserink girhub
enserink_lab_drug_sensitivity_org <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_drug_sensitivity.csv')

common_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/common_drugnames_pubchem.csv")
common_drugs$`Unnamed: 0_x` <- NULL
common_drugs$`Unnamed: 0_y` <- NULL
common_drugs$`Unnamed: 0` <- NULL

enserink_lab_drug_sensitivity <- inner_join(common_drugs, enserink_lab_drug_sensitivity_org, by=c("enserink_drug_name"="drug"))
colnames(enserink_lab_drug_sensitivity)[colnames(enserink_lab_drug_sensitivity) == "pubchem_drug_name"] <- "drug"
colnames(enserink_lab_drug_sensitivity)[colnames(enserink_lab_drug_sensitivity) == "Patient.ID"] <- "Patient.num"


compare_dss2_github_enserink <- subset(dss_github_enserink, select = c(Patient.num, drug, DSS2))
compare_dss2_github_enserink <- compare_dss2_github_enserink %>% arrange(desc(Patient.num), .by_group = FALSE)
compare_dss2_github_enserink <- compare_dss2_github_enserink %>% arrange(desc(drug), .by_group = FALSE)
compare_dss2_original_enserink <- subset(enserink_lab_drug_sensitivity, select = c(Patient.num, drug, DSS2))
compare_dss2_original_enserink <- compare_dss2_original_enserink %>% arrange(desc(Patient.num), .by_group = FALSE)
compare_dss2_original_enserink <- compare_dss2_original_enserink %>% arrange(desc(drug), .by_group = FALSE)
summary(comparedf(compare_dss2_github_enserink,compare_dss2_original_enserink))



dss_github_enserink_full_set_drugs <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink_full_drug_set_testing.csv')
dss_github_enserink_full_set_drugs <- inner_join(common_drugs, dss_github_enserink_full_set_drugs, by = c("enserink_drug_name"="drug"))
compare_dss2_github_enserink <- subset(dss_github_enserink_full_set_drugs, select = c(Patient.num, enserink_drug_name, DSS2))
compare_dss2_github_enserink <- compare_dss2_github_enserink %>% arrange(desc(Patient.num), .by_group = FALSE)
compare_dss2_github_enserink <- compare_dss2_github_enserink %>% arrange(desc(enserink_drug_name), .by_group = FALSE)
colnames(enserink_lab_drug_sensitivity_org)[colnames(enserink_lab_drug_sensitivity_org) == "pubchem_drug_name"] <- "drug"
colnames(enserink_lab_drug_sensitivity_org)[colnames(enserink_lab_drug_sensitivity_org) == "Patient.ID"] <- "Patient.num"
compare_dss2_original_enserink <- subset(enserink_lab_drug_sensitivity_org, select = c(Patient.num, drug, DSS2))
compare_dss2_original_enserink <- compare_dss2_original_enserink %>% arrange(desc(Patient.num), .by_group = FALSE)
compare_dss2_original_enserink <- compare_dss2_original_enserink %>% arrange(desc(drug), .by_group = FALSE)
summary(comparedf(compare_dss2_github_enserink,compare_dss2_original_enserink, by = c('enserink_drug_name'='drug', 'Patient.num')))



#Testing same scores for FIMM lab ----
#Getting R script generated scores
dss_github_fimm <- read_csv('~/Desktop/UiO/Project 1/Data/response scores/dss_github_fimm_full_drug_set_v1.csv')
fmm_supplementary <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set_supplemental.csv')
#Getting precalculated scores from enserink girhub
fimm_raw_response <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_raw_dose_responses.csv.csv')
fimm_raw_response <- unique(subset(fimm_raw_response, select = c(sample, drug, dss)))
colnames(fimm_raw_response)[colnames(fimm_raw_response) == "sample"] <- "Patient.num"
colnames(fimm_raw_response)[colnames(fimm_raw_response) == "dss"] <- "DSS2"

fimmm_dss_org <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_drug_response.csv')
fimmm_dss_org <- gather(fimmm_dss_org, Patient.num, DSS2, 'AML_084_04':'Healthy_17', factor_key=TRUE)
colnames(fimmm_dss_org)[colnames(fimmm_dss_org) == "Drug_name"] <- "drug"


compare_dss2_github_fimm <- subset(dss_github_fimm, select = c(Patient.num, drug, DSS2))
compare_dss2_github_fimm <- compare_dss2_github_fimm %>% arrange(desc(Patient.num), .by_group = FALSE)
compare_dss2_github_fimm <- compare_dss2_github_fimm %>% arrange(desc(drug), .by_group = FALSE)
compare_dss2_github_fimm <- compare_dss2_github_fimm %>% mutate(names = paste(Patient.num, drug)) %>% as.data.frame()
compare_dss2_github_fimm <- subset(compare_dss2_github_fimm, !is.na(Patient.num))


compare_dss2_original_fimm <- subset(fimm_raw_response, select = c(Patient.num, drug, DSS2))
compare_dss2_original_fimm <- compare_dss2_original_fimm %>% arrange(desc(Patient.num), .by_group = FALSE)
compare_dss2_original_fimm <- compare_dss2_original_fimm %>% arrange(desc(drug), .by_group = FALSE)
compare_dss2_original_fimm <- compare_dss2_original_fimm %>% mutate(names = paste(Patient.num, drug)) %>% as.data.frame()
compare_dss2_original_fimm <- subset(compare_dss2_original_fimm, !is.na(Patient.num))
summary(comparedf(compare_dss2_github_fimm,compare_dss2_original_fimm, by = 'names'))

cosine_similarity <- function(x, y) {
  return(sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))))
}

# Initialize an empty data frame to store the results
results <- data.frame(id1 = integer(), id2 = integer(), similarity = numeric())

# Compute similarity for each pair of rows
for (i in 1:nrow(fimmm_dss_org)) {
  for (j in 1:nrow(dss_github_fimm)) {
    similarity <- cosine_similarity(fimmm_dss_org[i, 'DSS2'], dss_github_fimm[j, 'DSS2'])
    results <- rbind(results, data.frame(id1 = fimmm_dss_org[i, 'Patient.num'], id2 = dss_github_fimm[j, 'Patient.num'], similarity = similarity))
  }
}

# Find the most similar pairs
most_similar_pairs <- results %>%
  group_by(id1) %>%
  top_n(1, wt = similarity)

print(most_similar_pairs)

df <- inner_join(fimmm_dss_org, dss_github_fimm, by = c('drug'))
ggplotRegression(lm(DSS2.x ~DSS2.y, df))


#Karolinska github ----
karolinska_github_fresh <- read_csv('~/Desktop/UiO/Project 1/Data/Response scores/dss_github_karolinska_full_drug_set_fresh.csv')
karolinska_github_frosen <- read_csv('~/Desktop/UiO/Project 1/Data/Response scores/dss_github_karolinska_full_drug_set_frosen.csv')


karolinska_dss2_org <- read_csv('~/Desktop/UiO/Project 1/Data/Second run/karolinska_normalised_reponse_breeze.csv')

colnames(karolinska_dss2_org)[colnames(karolinska_dss2_org) == "DRUG_NAME"] <- "drug"
colnames(karolinska_dss2_org)[colnames(karolinska_dss2_org) == "DSS"] <- "DSS2"


compare_dss2_github_enserink <- subset(karolinska_github_frosen, select = c(Patient.num, drug, DSS2))
compare_dss2_github_enserink <- compare_dss2_github_enserink %>% arrange(desc(Patient.num), .by_group = FALSE)
compare_dss2_github_enserink <- compare_dss2_github_enserink %>% arrange(desc(drug), .by_group = FALSE)
compare_dss2_github_enserink <- compare_dss2_github_enserink %>% mutate(names = paste(Patient.num, drug)) %>% as.data.frame()
rownames(compare_dss2_github_enserink) <- compare_dss2_github_enserink$names
compare_dss2_original_enserink <- subset(karolinska_dss2_org, sample == 'frozen',select = c(Patient.num, drug, DSS2))
compare_dss2_original_enserink <- compare_dss2_original_enserink %>% arrange(desc(Patient.num), .by_group = FALSE)
compare_dss2_original_enserink <- compare_dss2_original_enserink %>% arrange(desc(drug), .by_group = FALSE)
compare_dss2_original_enserink <- compare_dss2_original_enserink %>% mutate(names = paste(Patient.num, drug)) %>% as.data.frame()
rownames(compare_dss2_original_enserink) <- compare_dss2_original_enserink$names
summary(comparedf(compare_dss2_github_enserink,compare_dss2_original_enserink, by="names", tol.num.val = 10))
View(subset(karolinska_dss2_org, drug == 'Bortezomib' & Patient.num == 'AML_247'))

df <- inner_join(karolinska_github_fresh, karolinska_dss2_org, by = c('Patient.num', 'drug'))
grid_plots(df, 'DSS2.x', 'DSS2.y', 'Patient.num', "", "", "~/Desktop/UiO/Project 1/Figures/Karolinska/karolinska_dss_github_vs_org.png")

df <- inner_join(subset(karolinska_dss2_org, sample == 'fresh'), subset(karolinska_dss2_org, sample == 'frozen'), by = c('Patient.num', 'drug'), suffix = c("_fresh", "_frozen"))
grid_plots(df, 'DSS2_fresh', 'DSS2_frozen', 'Patient.num', "", "", "~/Desktop/UiO/Project 1/Figures/Karolinska/karolinska_dss_fresh_vs_frosen_org.png")
ggplotRegression(lm(DSS2_fresh ~ DSS2_frozen, data = df))


df <- inner_join(karolinska_github_fresh, karolinska_github_frosen, by = c('Patient.num', 'drug'))
grid_plots(df, 'DSS2.x', 'DSS2.y', 'Patient.num', "", "", "~/Desktop/UiO/Project 1/Figures/Karolinska/karolinska_dss_fresh_vs_frosen_github.png")

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

#Testing same scores for Beat AML conc uM and Beat AML conc nM lab ----
#Getting R script generated scores uM using concentration
dss_github_beataml <- read_csv('~/Desktop/UiO/Project 1/Data/Response scores/dss_github_beat_aml_full_drug_set_v1.csv')

#Getting R script generated scores nM using concentration
dss_github_beataml_v1 <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_beat_aml_full_drug_set_v1.csv')
dss_github_beataml_v1 <- unique(subset(dss_github_beataml_v1, select = c(sample, drug, dss)))
colnames(dss_github_beataml_v1)[colnames(dss_github_beataml_v1) == "sample"] <- "Patient.num"
colnames(dss_github_beataml_v1)[colnames(dss_github_beataml_v1) == "dss"] <- "DSS2"


compare_dss2_github_beataml <- subset(dss_github_beataml, select = c(Patient.num, drug, DSS2))
compare_dss2_github_beataml <- compare_dss2_github_beataml %>% arrange(desc(Patient.num), .by_group = FALSE)
compare_dss2_github_beataml <- compare_dss2_github_beataml %>% arrange(desc(drug), .by_group = FALSE)
compare_dss2_github_beataml <- compare_dss2_github_beataml %>% mutate(names = paste(Patient.num, drug)) %>% as.data.frame()
compare_dss2_github_beataml <- subset(compare_dss2_github_beataml, !is.na(Patient.num))


compare_dss2_beataml_v1 <- subset(dss_github_beataml_v1, select = c(Patient.num, drug, DSS2))
compare_dss2_beataml_v1 <- compare_dss2_beataml_v1 %>% arrange(desc(Patient.num), .by_group = FALSE)
compare_dss2_beataml_v1 <- compare_dss2_beataml_v1 %>% arrange(desc(drug), .by_group = FALSE)
compare_dss2_beataml_v1 <- compare_dss2_beataml_v1 %>% mutate(names = paste(Patient.num, drug)) %>% as.data.frame()
compare_dss2_beataml_v1 <- subset(compare_dss2_beataml_v1, !is.na(Patient.num))
sum <- summary(comparedf(compare_dss2_github_beataml,compare_dss2_beataml_v1, by = 'names'))
sum$diffs.table



#Compare FIMM values 

FIMM_supp <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set_supplemental.csv')
