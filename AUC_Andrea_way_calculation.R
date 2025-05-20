
beat_aml_raw_response <- read.csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_inhibitor_cleaned.csv')
beat_aml_raw_response$avg2_response <- beat_aml_raw_response$avg2_response/100
beat_aml_raw_response <- beat_aml_raw_response %>%
  mutate(patient_id = paste(dbgap_subject_id, dbgap_dnaseq_sample,dbgap_rnaseq_sample, sep = "_")) %>% as.data.frame()

common_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/common_drugnames_pubchem.csv")
common_drugs$`Unnamed: 0_x` <- NULL
common_drugs$`Unnamed: 0_y` <- NULL
common_drugs$`Unnamed: 0` <- NULL

beat_aml_for_dss <- inner_join(common_drugs, beat_aml_raw_response, by=c("beat_aml_drug_name"="inhibitor"))
colnames(beat_aml_for_dss)[colnames(beat_aml_for_dss) == "pubchem_drug_name"] <- "drug"

# Display the dataframe
duplicates <- beat_aml_for_dss %>%
  group_by(patient_id, drug) %>%
  filter(n() > 1) %>%
  ungroup()

# Show the duplicate values
print(duplicates)
print(beat_aml_for_dss)

duplicates <- beat_aml_for_dss %>%
  arrange(patient_id, drug) %>%
  group_by(patient_id, drug) %>%
  filter(n() > 2) %>%
  ungroup() 


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
auc_beat_aml <- beat_aml_for_dss %>% 
  group_by(patient_id, drug) %>% 
  summarise(auc_a = Log10.DSS(well_concentration, avg2_response)) 

for(j in unique(beat_aml_for_dss$drug)){
  beat_aml_for_dss$auc_a[beat_aml_for_dss$drug==j] <-
    Log10.DSS(beat_aml_for_dss$well_concentration[beat_aml_for_dss$drug==j], beat_aml_for_dss$avg2_response[beat_aml_for_dss$drug==j])
}

duplicates <- auc_beat_aml %>%
  arrange(patient_id, drug) %>%
  group_by(patient_id, drug) %>%
  filter(n() > 2) %>%
  ungroup()

ggplot(beat_aml_for_dss, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")

ggplot(auc_beat_aml, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")



karolinska_raw_response = read.csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/karolinska_normalised_reponse_raw.csv')
karolinska_raw_response$Normalised_reponse = 1 - karolinska_raw_response$Normalised_reponse/100
karolinska_raw_response <- karolinska_raw_response %>%
  mutate(patient_id = paste('k',Patient.num, sep='')) %>% as.data.frame()

karolinska_raw_response <- inner_join(common_drugs, karolinska_raw_response, by=c("karolinska_drug_name"="DRUG_NAME"))
colnames(karolinska_raw_response)[colnames(karolinska_raw_response) == "pubchem_drug_name"] <- "drug"

auc_karolinska <- karolinska_raw_response %>% 
  group_by(patient_id, drug) %>% 
  summarise(auc_a = Log10.DSS(CONCENTRATION, Normalised_reponse), .groups = 'drop') 

ggplot(auc_karolinska, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")


fimm_raw_response <- read.csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_raw_dose_responses.csv.csv')
fimm_raw_response$doseResponses <- fimm_raw_response$doseResponses/100
fimm_raw_response <- subset(fimm_raw_response, drug != 'TG100-115')
fimm_raw_response <- inner_join(common_drugs, fimm_raw_response, by=c("fimm_drug_name"="drug"))
colnames(fimm_raw_response)[colnames(fimm_raw_response) == "pubchem_drug_name"] <- "drug"
colnames(fimm_raw_response)[colnames(fimm_raw_response) == "sample"] <- "patient_id"

test_fimm <- read_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_4_DSRT_assay_details_164S_DM.xlsx')

summary_table <- test_fimm %>% summarize(count = n_distinct(Sample_ID))


auc_FIMM <- fimm_raw_response %>% 
  group_by(patient_id, drug) %>% 
  summarise(auc_a = Log10.DSS(Final.Conc, doseResponses), .groups = 'drop') 

ggplot(auc_FIMM, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")


enserink_breeze_input <- read_csv("/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_dose_response.csv")
enserink_breeze_input <- inner_join(common_drugs, enserink_breeze_input, by=c("enserink_drug_name"="drug"))
colnames(enserink_breeze_input)[colnames(enserink_breeze_input) == "pubchem_drug_name"] <- "drug"
enserink_breeze_input$dose <- enserink_breeze_input$dose*10^9
enserink_breeze_input$response <-1 - enserink_breeze_input$y_true/100
colnames(enserink_breeze_input)[colnames(enserink_breeze_input) == "Patient.ID"] <- "patient_id"


auc_enserink <- enserink_breeze_input %>% 
  group_by(patient_id, drug) %>% 
  summarise(auc_a = Log10.DSS(dose, response), .groups = 'drop') 

ggplot(auc_enserink, aes(x=auc_a)) +
  geom_density(color="darkblue", fill="lightblue")

all_labs_auc <- rbind(auc_beat_aml, auc_karolinska, auc_FIMM, auc_enserink)
colnames(all_labs_auc)[colnames(all_labs_auc) == "patient_id"] <- "Patient.num"
write_csv(all_labs_auc, '~/Desktop/UiO/Project 1/Data/Initial cleansing/auc_calculation_all_labs.csv')



enserink_sensitivity <- read_csv("/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_drug_sensitivity.csv")
colnames(enserink_sensitivity)[colnames(enserink_sensitivity) == "Patient.ID"] <- "patient_id"
df1 <- subset(enserink_sensitivity, select = c("patient_id", "drug", "rAUC"))
df2 <- subset(auc_enserink, select = c("patient_id", "drug", "auc_a"))
colnames(df2)[colnames(df2) == "auc_a"] <- "rAUC"

library(arsenal)
summary(comparedf(df1,df2, by=c("patient_id", "drug")))
