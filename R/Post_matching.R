#Read libraries
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(stats)
library(ggplot2)

#Read DSS info for each dataset
oslo_response_scores <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink_full_drug_set.csv')
oslo_response_scores$lab <- 'Enserink'

fimm_response_score <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set.csv')
fimm_response_score$lab <- 'FIMM'

#Getting the common drugs across the 4 datasets
common_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/common_drugnames_pubchem.csv")
common_drugs$`Unnamed: 0_x` <- NULL
common_drugs$`Unnamed: 0_y` <- NULL
common_drugs$`Unnamed: 0` <- NULL

#Enserink common drugs 
dss_enserink_common_drugs <- inner_join(common_drugs, oslo_response_scores, by=c("enserink_drug_name"="drug"))
colnames(dss_enserink_common_drugs)[colnames(dss_enserink_common_drugs) == "pubchem_drug_name"] <- "drug"
dss_enserink_common_drugs$Patient_ID <- gsub('_.*', '', dss_enserink_common_drugs$Patient.num)
dim(dss_enserink_common_drugs)

#dss_fimm_common_drugs <- dss_fimm_common_drugs %>% mutate(Patient.num = keep_letters_numbers(Patient.num))
dim(fimm_response_score)
fimm_response_score$Patient_ID <- gsub("(_[^_]*){1}$", "", fimm_response_score$Patient.num)
dss_FIMM_common_drugs <- inner_join(common_drugs, fimm_response_score, by=c("fimm_drug_name"="drug"))
colnames(dss_FIMM_common_drugs)[colnames(dss_FIMM_common_drugs) == "pubchem_drug_name"] <- "drug"
dim(dss_FIMM_common_drugs)

matching <- read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/matching/FIMM_Oslo_matching.csv')

fimm_matching <- merge(dss_FIMM_common_drugs, matching, by.x = 'Patient_ID', by.y = 'start')
oslo_matching <- merge(dss_enserink_common_drugs, matching, by.x = 'Patient_ID', by.y = 'end')

fimm_oslo_mathcing <- merge(fimm_matching, oslo_matching, by = c('common_id', 'drug'), suffixes = c(".fimm",".oslo"))


 
fit <- lm(DSS2.fimm ~ DSS2.oslo, data = fimm_oslo_mathcing)
ggplot(fimm_oslo_mathcing, aes(x = DSS2.fimm, y = DSS2.oslo)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "blue") +
  facet_wrap(~drug) +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "Intercept =",signif(fit$coef[[1]],5 ), " Slope =",signif(fit$coef[[2]], 5), " P =",signif(summary(fit)$coef[2,4], 5)))  

get_model_stats <- function(data) {
  fit <- lm(DSS2.fimm ~ DSS2.oslo, data = data)
  return(data.frame(
    drug = unique(data$drug),
    adj_r2 = signif(summary(fit)$adj.r.squared, 5),
    intercept = signif(fit$coef[[1]], 5),
    slope = signif(fit$coef[[2]], 5),
    p_value = signif(summary(fit)$coef[2, 4], 5)
  ))
}

# Apply the function to each subset of data by 'drug'
model_stats <- fimm_oslo_mathcing %>%
  group_by(drug) %>%
  do(get_model_stats(.))

# Merge the statistics back into the original data
df <- fimm_oslo_mathcing %>%
  left_join(model_stats, by = "drug")

# Plot with facet_wrap and adjusted R^2 in the title of each facet
ggplot(df, aes(x = DSS2.fimm, y = DSS2.oslo)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "blue") +
  facet_wrap(~drug) +
  labs(title = "Linear Model for each Drug") +
  theme(strip.text = element_text(size = 12)) +
  geom_text(aes(x = Inf, y = Inf, 
                label = paste("Adj R2 =", adj_r2)), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "red", 
            inherit.aes = FALSE, check_overlap = TRUE)


#Search optimal response transformation based on Box-Cox Transformations for Linear Models
fimm_oslo_mathcing$DSS2_pos <- fimm_oslo_mathcing$DSS2.fimm + 0.01
MASS::boxcox(DSS2_pos ~ drug, 
             data = fimm_oslo_mathcing,
             lambda = seq(-0.25, 2, length.out = 10))

fimm_oslo_mathcing$DSS2_boxcox_fimm <- ((fimm_oslo_mathcing$DSS2.fimm)^0.5 - 1) / 0.5


for(j in unique(fimm_oslo_mathcing$drug)){
  fimm_oslo_mathcing$DSS2_boxcox_sclaed_fimm[fimm_oslo_mathcing$drug==j] <-
    scale(fimm_oslo_mathcing$DSS2_boxcox_fimm[fimm_oslo_mathcing$drug==j])
}


fimm_oslo_mathcing$DSS2_pos <- fimm_oslo_mathcing$DSS2.oslo + 0.01
MASS::boxcox(DSS2_pos ~ drug, 
             data = fimm_oslo_mathcing,
             lambda = seq(-0.25, 2, length.out = 10))

fimm_oslo_mathcing$DSS2_boxcox_oslo <- ((fimm_oslo_mathcing$DSS2.oslo)^0.5 - 1) / 0.5


for(j in unique(fimm_oslo_mathcing$drug)){
  fimm_oslo_mathcing$DSS2_boxcox_sclaed_oslo[fimm_oslo_mathcing$drug==j] <-
    scale(fimm_oslo_mathcing$DSS2_boxcox_oslo[fimm_oslo_mathcing$drug==j])
}


ggplot(fimm_oslo_mathcing, aes(x = DSS2_boxcox_sclaed_fimm, y = DSS2_boxcox_sclaed_oslo)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "blue") +
  facet_wrap(~drug) +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "Intercept =",signif(fit$coef[[1]],5 ), " Slope =",signif(fit$coef[[2]], 5), " P =",signif(summary(fit)$coef[2,4], 5)))  

get_model_stats <- function(data) {
  fit <- lm(DSS2_boxcox_sclaed_fimm ~ DSS2_boxcox_sclaed_oslo, data = data)
  return(data.frame(
    drug = unique(data$drug),
    adj_r2 = signif(summary(fit)$adj.r.squared, 5),
    intercept = signif(fit$coef[[1]], 5),
    slope = signif(fit$coef[[2]], 5),
    p_value = signif(summary(fit)$coef[2, 4], 5)
  ))
}

# Apply the function to each subset of data by 'drug'
model_stats <- fimm_oslo_mathcing %>%
  group_by(drug) %>%
  do(get_model_stats(.))

# Merge the statistics back into the original data
df <- fimm_oslo_mathcing %>%
  left_join(model_stats, by = "drug")

# Plot with facet_wrap and adjusted R^2 in the title of each facet
ggplot(df, aes(x = DSS2_boxcox_sclaed_fimm, y = DSS2_boxcox_sclaed_oslo)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "blue") +
  facet_wrap(~drug) +
  labs(title = "Linear Model for each Drug") +
  theme(strip.text = element_text(size = 12)) +
  geom_text(aes(x = Inf, y = Inf, 
                label = paste("Adj R2 =", adj_r2)), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "red", 
            inherit.aes = FALSE, check_overlap = TRUE)
