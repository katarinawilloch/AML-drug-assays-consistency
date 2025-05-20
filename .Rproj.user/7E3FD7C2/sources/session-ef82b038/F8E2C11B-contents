if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("GRmetrics")

library(GRmetrics)
data(inputCaseA)

inputCaseA = as.data.frame(inputCaseA)

drc_output = GRfit(inputCaseA, groupingVariables = c('cell_line','agent'))

head(GRgetMetrics(drc_output))

GRgetGroupVars(drc_output)

write.table(GRgetMetrics(drc_output), file = "filename.tsv", quote = FALSE,
            sep = "\t", row.names = FALSE)



karolinska_raw_response = read.csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/karolinska_normalised_reponse_raw.csv')
subset(karolinska_raw_response, DRUG_NAME == 'Bortezomib')
karolinska_raw_response$Normalised_reponse = 1 - karolinska_raw_response$Normalised_reponse/100
karolinska_raw_response <- karolinska_raw_response %>%
  mutate(patient_id = paste('k',Patient.num, sep='')) %>% as.data.frame()

karolinska_raw_response <- inner_join(common_drugs, karolinska_raw_response, by=c("karolinska_drug_name"="DRUG_NAME"))
colnames(karolinska_raw_response)[colnames(karolinska_raw_response) == "pubchem_drug_name"] <- "drug"
colnames(karolinska_raw_response)[colnames(karolinska_raw_response) == "CONCENTRATION"] <- "concentration"
colnames(karolinska_raw_response)[colnames(karolinska_raw_response) == "Normalised_reponse"] <- "cell_count"
colnames(karolinska_raw_response)[colnames(karolinska_raw_response) == "x"] <- "cell_count__time0"
colnames(karolinska_raw_response)[colnames(karolinska_raw_response) == "X"] <- "cell_count__ctrl"


drc_output = GRfit(karolinska_raw_response, groupingVariables = c('patient_id','drug'), case = "B")
GRdrawDRC(drc_output)


enserink_breeze_input <- read_csv("/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_dose_response.csv")
enserink_breeze_input <- inner_join(common_drugs, enserink_breeze_input, by=c("enserink_drug_name"="drug"))
colnames(enserink_breeze_input)[colnames(enserink_breeze_input) == "pubchem_drug_name"] <- "drug"
enserink_breeze_input$response <-1 - enserink_breeze_input$y_true/100
colnames(enserink_breeze_input)[colnames(enserink_breeze_input) == "Patient.ID"] <- "patient_id"

subset(enserink_breeze_input, drug == 'Bortezomib')


fimm_raw_response <- read.csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_raw_dose_responses.csv.csv')
subset(fimm_raw_response, drug == 'Bortezomib')
fimm_raw_response$doseResponses <- fimm_raw_response$doseResponses/100
fimm_raw_response <- subset(fimm_raw_response, drug != 'TG100-115')
fimm_raw_response <- inner_join(common_drugs, fimm_raw_response, by=c("fimm_drug_name"="drug"))
colnames(fimm_raw_response)[colnames(fimm_raw_response) == "pubchem_drug_name"] <- "drug"
colnames(fimm_raw_response)[colnames(fimm_raw_response) == "sample"] <- "patient_id"




beat_aml_raw_response <- read.csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_inhibitor_cleaned.csv')
beat_aml_raw_response$avg2_response <- beat_aml_raw_response$avg2_response/100
beat_aml_raw_response <- beat_aml_raw_response %>%
  mutate(patient_id = paste(dbgap_subject_id, dbgap_dnaseq_sample,dbgap_rnaseq_sample, sep = "_")) %>% as.data.frame()