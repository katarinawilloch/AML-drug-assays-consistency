
library(ComplexHeatmap)

list_drugs = ['Cytarabine', '345627-80-7', '781661-94-7', 'Borzetomib', 'PH-797804', 'Quizartinib', 'PONATINIB', 'Ruxalitinib', 'Sorafenib', 'cyt387']

densityHeatmap()


ggplot(result, aes(x = lab, y = AVG_score_per_quantile, color = drug, group=interaction(drug, Quantile))) +
  geom_point() +
  geom_line() +
  labs(title = "Line plot",
       x = "Lab",
       y = "Average DSS score per drug within each quantile",
       color = "Drugs") +
  theme_minimal()

for(d in unique(result$drug)){
  s <- subset(result, drug == d, select=AVG_score_per_quantile)
  print(s)
  hist(s$AVG_score_per_quantile)
}
hist(result$AVG_score_per_quantile)

regression_model <- lm(AVG_score_per_quantile~drug, data=result)
summary(regression_model)
ggplotRegression(regression_model)
prelim_plot <- ggplot(result, aes(x = drug, y = AVG_score_per_quantile, color=lab)) +
  geom_point() +
  geom_smooth(method = "lm")
plot(regression_model)

split_plot <- ggplot(aes(AVG_score_per_quantile, medium), data = result) + 
  geom_point() + 
  facet_wrap(~ lab) + # create a facet for each mountain range
  xlab("length") + 
  ylab("test score")

dot_model <- lmer(AVG_score_per_quantile ~ 1 + (1|lab) + (1|Quantile) + (1|drug), result)
dotplot(ranef(dot_model, condVar=TRUE))

lmer_model1 <- lmer(AVG_score_per_quantile ~ medium + (1|lab), data=result)
summary(lmer_model1)
confint(lmer_model1)




###----Patient.num as random effect----
####----Time_until_sample_usage----
mixed_model_time_until_sample_usage_box_cox_scaled_patient_re <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1 |Patient.num),all_datasets_v2)

####----Medium----
mixed_model_medium_box_cox_scaled_patient_re <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1|Patient.num),all_datasets_v2)

####----Cells----
mixed_model_cells_box_cox_scaled_patient_re <- lmer(DSS2_boxcox_sclaed2 ~ cells + (1|Patient.num),all_datasets_v2)

####----Microenv stimuli----
mixed_model_micro_env_stimuli_box_cox_scaled_patient_re <- lmer(DSS2_boxcox_sclaed2 ~ microenvironmental_stimuli + (1|Patient.num),all_datasets_v2)

####----Cell counting method----
mixed_model_cell_counting_method_box_cox_scaled_patient_re <- lmer(DSS2_boxcox_sclaed2 ~ cell_counting_method + (1|Patient.num),all_datasets_v2)

####----Positive Control + Sensitivity Readout + Nr of Concentration Points
mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_patient_re <- lmer(DSS2_boxcox_sclaed2 ~ positive_control + (1|Patient.num),all_datasets_v2)

####----Centrifugation Procedure ----
mixed_model_centrifugation_procedure_box_cox_scaled_patient_re <- lmer(DSS2_boxcox_sclaed2 ~ centrifugation_procedure + (1|Patient.num),all_datasets_v2)

####----Plate Reader----
mixed_model_plate_reader_box_cox_scaled_patient_re <- lmer(DSS2_boxcox_sclaed2 ~ plate_reader + (1|Patient.num),all_datasets_v2)


####----Table of crossed re----
patient_re_models <-c(mixed_model_time_until_sample_usage_box_cox_scaled_patient_re, 
                      mixed_model_medium_box_cox_scaled_patient_re, 
                      mixed_model_cells_box_cox_scaled_patient_re, 
                      mixed_model_micro_env_stimuli_box_cox_scaled_patient_re, 
                      mixed_model_cell_counting_method_box_cox_scaled_patient_re, 
                      mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_patient_re, 
                      mixed_model_centrifugation_procedure_box_cox_scaled_patient_re, 
                      mixed_model_plate_reader_box_cox_scaled_patient_re)
df_patient_re_models <- model_results(patient_re_models)
rownames(df_patient_re_models) <- NULL
df_patient_re_models$Random_Effect_Term <- '1 | Patient.num'
formattable(df_patient_re_models, list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))

###----Drug as random effect----
####----Time_until_sample_usage----
mixed_model_time_until_sample_usage_box_cox_scaled_drug_re <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1 |drug),all_datasets_v2)

####----Medium----
mixed_model_medium_box_cox_scaled_drug_re <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1|drug),all_datasets_v2)

####----Cells----
mixed_model_cells_box_cox_scaled_drug_re <- lmer(DSS2_boxcox_sclaed2 ~ cells + (1|drug),all_datasets_v2)

####----Microenv stimuli----
mixed_model_micro_env_stimuli_box_cox_scaled_drug_re <- lmer(DSS2_boxcox_sclaed2 ~ microenvironmental_stimuli + (1|drug),all_datasets_v2)

####----Cell counting method----
mixed_model_cell_counting_method_box_cox_scaled_drug_re <- lmer(DSS2_boxcox_sclaed2 ~ cell_counting_method + (1|drug),all_datasets_v2)

####----Positive Control + Sensitivity Readout + Nr of Concentration Points
mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_drug_re <- lmer(DSS2_boxcox_sclaed2 ~ positive_control + (1|drug),all_datasets_v2)

####----Centrifugation Procedure ----
mixed_model_centrifugation_procedure_box_cox_scaled_drug_re <- lmer(DSS2_boxcox_sclaed2 ~ centrifugation_procedure + (1|drug),all_datasets_v2)

####----Plate Reader----
mixed_model_plate_reader_box_cox_scaled_drug_re <- lmer(DSS2_boxcox_sclaed2 ~ plate_reader + (1|drug),all_datasets_v2)


####----Table of crossed re----
drug_re_models <-c(mixed_model_time_until_sample_usage_box_cox_scaled_drug_re, 
                   mixed_model_medium_box_cox_scaled_drug_re, 
                   mixed_model_cells_box_cox_scaled_drug_re, 
                   mixed_model_micro_env_stimuli_box_cox_scaled_drug_re, 
                   mixed_model_cell_counting_method_box_cox_scaled_drug_re, 
                   mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_drug_re, 
                   mixed_model_centrifugation_procedure_box_cox_scaled_drug_re, 
                   mixed_model_plate_reader_box_cox_scaled_drug_re)
df_drug_re_models <- model_results(drug_re_models)
rownames(df_drug_re_models) <- NULL
df_drug_re_models$Random_Effect_Term <- '1 | drug'
formattable(df_drug_re_models, list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))

###----Experimental var as random effect----
####----Time_until_sample_usage----
mixed_model_time_until_sample_usage_box_cox_scaled_e_var_re <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1 |time_until_sample_usage),all_datasets_v2)

####----Medium----
mixed_model_medium_box_cox_scaled_e_var_re <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1|medium),all_datasets_v2)

####----Cells----
mixed_model_cells_box_cox_scaled_e_var_re <- lmer(DSS2_boxcox_sclaed2 ~ cells + (1|cells),all_datasets_v2)

####----Microenv stimuli----
mixed_model_micro_env_stimuli_box_cox_scaled_e_var_re <- lmer(DSS2_boxcox_sclaed2 ~ microenvironmental_stimuli + (1|microenvironmental_stimuli),all_datasets_v2)

####----Cell counting method----
mixed_model_cell_counting_method_box_cox_scaled_e_var_re <- lmer(DSS2_boxcox_sclaed2 ~ cell_counting_method + (1|cell_counting_method),all_datasets_v2)

####----Positive Control + Sensitivity Readout + Nr of Concentration Points
mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_e_var_re <- lmer(DSS2_boxcox_sclaed2 ~ positive_control + (1|positive_control),all_datasets_v2)

####----Centrifugation Procedure ----
mixed_model_centrifugation_procedure_box_cox_scaled_e_var_re <- lmer(DSS2_boxcox_sclaed2 ~ centrifugation_procedure + (1|centrifugation_procedure),all_datasets_v2)

####----Plate Reader----
mixed_model_plate_reader_box_cox_scaled_e_var_re <- lmer(DSS2_boxcox_sclaed2 ~ plate_reader + (1|plate_reader),all_datasets_v2)


####----Table of crossed re----
e_var_re_models <-c(mixed_model_time_until_sample_usage_box_cox_scaled_e_var_re, 
                    mixed_model_medium_box_cox_scaled_e_var_re, 
                    mixed_model_cells_box_cox_scaled_e_var_re, 
                    mixed_model_micro_env_stimuli_box_cox_scaled_e_var_re, 
                    mixed_model_cell_counting_method_box_cox_scaled_e_var_re, 
                    mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_e_var_re, 
                    mixed_model_centrifugation_procedure_box_cox_scaled_e_var_re, 
                    mixed_model_plate_reader_box_cox_scaled_e_var_re)
df_e_var_re_models <- model_results(e_var_re_models)
rownames(df_e_var_re_models) <- NULL
df_e_var_re_models$Random_Effect_Term <- '1 | Experimental var'
formattable(df_e_var_re_models, list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))

###----1 + E_var | drug----
####----Time_until_sample_usage----
mixed_model_time_until_sample_usage_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1 + time_until_sample_usage|drug),all_datasets_v2)
icc(mixed_model_time_until_sample_usage_box_cox_scaled)
ranova(mixed_model_time_until_sample_usage_box_cox_scaled)
r2(mixed_model_time_until_sample_usage_box_cox_scaled)
AIC(mixed_model_time_until_sample_usage_box_cox_scaled)

check_model(mixed_model_time_until_sample_usage_box_cox_scaled)

####----Medium----
mixed_model_medium_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1 + medium|drug),all_datasets_v2)
icc(mixed_model_medium_box_cox_scaled)
ranova(mixed_model_medium_box_cox_scaled)

####----Cells----
mixed_model_cells_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ cells + (1 + cells|drug),all_datasets_v2)

####----Microenv stimuli----
mixed_model_micro_env_stimuli_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ microenvironmental_stimuli + (1 + microenvironmental_stimuli|drug),all_datasets_v2)

####----Cell counting method----
mixed_model_cell_counting_method_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ cell_counting_method + (1 + cell_counting_method|drug),all_datasets_v2)

####----Positive Control + Sensitivity Readout + Nr of Concentration Points
mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ positive_control + (1 + positive_control|drug),all_datasets_v2)

####----Centrifugation Procedure ----
mixed_model_centrifugation_procedure_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ centrifugation_procedure + (1 + centrifugation_procedure|drug),all_datasets_v2)

####----Plate Reader----
mixed_model_plate_reader_box_cox_scaled <- lmer(DSS2_boxcox_sclaed2 ~ plate_reader + (1 + plate_reader|drug),all_datasets_v2)


####----Table of crossed re----
crossed_re_models <-c(mixed_model_time_until_sample_usage_box_cox_scaled, mixed_model_medium_box_cox_scaled, mixed_model_cells_box_cox_scaled, mixed_model_micro_env_stimuli_box_cox_scaled, mixed_model_cell_counting_method_box_cox_scaled, mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled, mixed_model_centrifugation_procedure_box_cox_scaled, mixed_model_plate_reader_box_cox_scaled)
df_crossed_re_models <- model_results(crossed_re_models)
df_crossed_re_models$Random_Effect_Term <- '1 + Experimental var | drug'
rownames(df_crossed_re_models) <- NULL
formattable(df_crossed_re_models, list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))




###---(1|E_var) + (1|drug)----
####----Time_until_sample_usage----
mixed_model_time_until_sample_usage_box_cox_scaled_two <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1 | time_until_sample_usage) + (1|drug),all_datasets_v2)
icc(mixed_model_time_until_sample_usage_box_cox_scaled_two)

####----Medium----
mixed_model_medium_box_cox_scaled_two <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1 | medium) + (1|drug),all_datasets_v2)
icc(mixed_model_medium_box_cox_scaled_two)

####----Cells----
mixed_model_cells_box_cox_scaled_two <- lmer(DSS2_boxcox_sclaed2 ~ cells + (1 | cells) + (1 | drug),all_datasets_v2)

####----Microenv stimuli----
mixed_model_micro_env_stimuli_box_cox_scaled_two <- lmer(DSS2_boxcox_sclaed2 ~ microenvironmental_stimuli + (1 | microenvironmental_stimuli) + (1|drug),all_datasets_v2)

####----Cell counting method----
mixed_model_cell_counting_method_box_cox_scaled_two <- lmer(DSS2_boxcox_sclaed2 ~ cell_counting_method + (1 | cell_counting_method) + (1 |drug),all_datasets_v2)

####----Positive Control + Sensitivity Readout + Nr of Concentration Points
mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_two <- lmer(DSS2_boxcox_sclaed2 ~ positive_control + (1 | positive_control) + (1|drug),all_datasets_v2)

####----Centrifugation Procedure ----
mixed_model_centrifugation_procedure_box_cox_scaled_two <- lmer(DSS2_boxcox_sclaed2 ~ centrifugation_procedure + (1 | centrifugation_procedure) + (1|drug),all_datasets_v2)

####----Plate Reader----
mixed_model_plate_reader_box_cox_scaled_two <- lmer(DSS2_boxcox_sclaed2 ~ plate_reader + (1 | plate_reader) + (1 |drug),all_datasets_v2)


####----Table of crossed re----
two_re_models <-c(mixed_model_time_until_sample_usage_box_cox_scaled_two, mixed_model_medium_box_cox_scaled_two, mixed_model_cells_box_cox_scaled_two, mixed_model_micro_env_stimuli_box_cox_scaled_two, mixed_model_cell_counting_method_box_cox_scaled_two, mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_two, mixed_model_centrifugation_procedure_box_cox_scaled_two, mixed_model_plate_reader_box_cox_scaled_two)
df_two_re_models <- model_results(two_re_models)
df_two_re_models$Random_Effect_Term <- '(1 | Experimental var) + (1 | drug)'
rownames(df_two_re_models) <- NULL
formattable(df_two_re_models,  list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))

###---1|E_var/drug----
####----Time_until_sample_usage----
mixed_model_time_until_sample_usage_box_cox_scaled_nested <- lmer(DSS2_boxcox_sclaed2 ~ time_until_sample_usage + (1 | time_until_sample_usage/drug),all_datasets_v2)
icc(mixed_model_time_until_sample_usage_box_cox_scaled_nested)

####----Medium----
mixed_model_medium_box_cox_scaled_nested <- lmer(DSS2_boxcox_sclaed2 ~ medium + (1 | medium/drug),all_datasets_v2)
icc(mixed_model_medium_box_cox_scaled_nested)

####----Cells----
mixed_model_cells_box_cox_scaled_nested <- lmer(DSS2_boxcox_sclaed2 ~ cells + (1 | cells/drug),all_datasets_v2)

####----Microenv stimuli----
mixed_model_micro_env_stimuli_box_cox_scaled_nested <- lmer(DSS2_boxcox_sclaed2 ~ microenvironmental_stimuli + (1 | microenvironmental_stimuli/drug),all_datasets_v2)

####----Cell counting method----
mixed_model_cell_counting_method_box_cox_scaled_nested <- lmer(DSS2_boxcox_sclaed2 ~ cell_counting_method + (1 | cell_counting_method/drug),all_datasets_v2)

####----Positive Control + Sensitivity Readout + Nr of Concentration Points
mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_nested <- lmer(DSS2_boxcox_sclaed2 ~ positive_control + (1 | positive_control/drug),all_datasets_v2)

####----Centrifugation Procedure ----
mixed_model_centrifugation_procedure_box_cox_scaled_nested <- lmer(DSS2_boxcox_sclaed2 ~ centrifugation_procedure + (1 | centrifugation_procedure/drug),all_datasets_v2)

####----Plate Reader----
mixed_model_plate_reader_box_cox_scaled_nested <- lmer(DSS2_boxcox_sclaed2 ~ plate_reader + (1 | plate_reader/drug),all_datasets_v2)


####----Table of crossed re----
nested_re_models <-c(mixed_model_time_until_sample_usage_box_cox_scaled_nested, mixed_model_medium_box_cox_scaled_nested, mixed_model_cells_box_cox_scaled_nested, mixed_model_micro_env_stimuli_box_cox_scaled_nested, mixed_model_cell_counting_method_box_cox_scaled_nested, mixed_model_sensitivity_readout_and_positive_control_box_cox_scaled_nested, mixed_model_centrifugation_procedure_box_cox_scaled_nested, mixed_model_plate_reader_box_cox_scaled_nested)
df_nested_re_models <- model_results(nested_re_models)
df_nested_re_models$Random_Effect_Term <- '1 | Experimental var / drug'
rownames(df_nested_re_models) <- NULL
formattable(df_nested_re_models,  list(
  Fixed_Effect_Term = formatter("span", style = ~formattable::style(color = "blue")),
  Fixed_Effect_Coefficient = color_tile("white", "lightblue"),
  p_value = formatter("span", style = x ~ formattable::style(color = ifelse(x < 0.05, "red", "black"))),
  significance_stars = formatter("span", style = ~formattable::style(color = "darkred"))
))


###----bar plot ----
all_models_different_re <- rbind(df_crossed_re_models, df_patient_re_models, df_drug_re_models, df_e_var_re_models, df_nested_re_models, df_two_re_models, df_re_intercept_models)
p <- ggplot(all_models_different_re, aes(x = Fixed_Effect_Term, y = Fixed_Effect_Coefficient, fill = Random_Effect_Term)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = significance_stars),
            position = position_dodge(width = 0.9),  # Adjust for dodge position
            vjust = -0.5,                           # Adjust vertical position above bars
            size = 5,                              # Text size
            color = "black") +                     # Text color
  labs(
    title = "Fixed Effects and Their Coefficients",
    x = "Fixed Effect Term",
    y = "Coefficient"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10)  # Rotate x-axis text
  )

# Print the plot
print(p)
# Save the plot to a file
ggsave("~/Desktop/UiO/Project 1/Figures/All_re_models.png", plot = p, width = 10, height = 20, dpi = 300)



###----All_Model_re----
models_re <- c(nested_re_models, crossed_re_models, two_re_models, e_var_re_models, drug_re_models, patient_re_models, re_intercept_models)

df <- data.frame(icc_df = data.frame())
for (model in re_intercept_models){
  temp_df <- model_diagnostics(model)
  df <- rbind(df, temp_df)
  r2 <- r2(model)
  
  # View the results
  print(r2)
  print(AIC(model))
}
df <- df[order(-df$icc_df.design_effect), ]
rownames(df) <- NULL
formattable(df)




# Load necessary library
library(rlang) # for use of the `:=` operator

# Compute weights based on fitted values or residuals
weights <- 1 / (fitted(mixed_model_time_until_sample_usage)^2) # Example weight calculation, adjust as needed
weights_1 <- 1 / lm(abs(resid(mixed_model_time_until_sample_usage)) ~ fitted(mixed_model_time_until_sample_usage))$fitted.values^2

# Fit the WLS model
wls_model <- lm(DSS2 ~ time_until_sample_usage, data = all_datasets_v2, weights = weights)

summary(wls_model)

check_model(wls_model)




#r <-sample(2:length(rownames(all_datasets_v2)),200, replace=F) 
View(all_datasets_v2)
all_datasets_v2$time_until_sample_usage_1 <- ifelse(all_datasets_v2$time_until_sample_usage == "1-2h after receiving",1,0)
glmm_model <- glmer(DSS2 ~ time_until_sample_usage_1 + (1 + cells|drug), data = all_datasets_v2, family = poisson, glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
fixef(glmm_model)
check_model(glmm_model)

xx <- lmer(DSS2 ~ time_until_sample_usage_1 + (1 + lab|drug) , all_datasets_v2)
check_model(xx)
summary(xx)





# Number of observations per group
n_per_group <- tapply(all_datasets_v2$DSS2_boxcox, all_datasets_v2$drug, length)

# Group-specific means and standard deviations
means <- tapply(all_datasets_v2$DSS2_boxcox, all_datasets_v2$drug, mean)
sds <- tapply(all_datasets_v2$DSS2_boxcox, all_datasets_v2$drug, sd)

# Overall (weighted) mean and standard deviation
overall_mean <- sum(means * n_per_group) / sum(n_per_group)
overall_sd <- sqrt(sum((sds^2 + means^2) * n_per_group) / sum(n_per_group) - overall_mean^2)



# Get the fixed effect coefficient for your categorical variable
fixed_effect_coeff <- fixef(mixed_model_time_until_sample_usage_box_cox_scaled_re_intercept)["time_until_sample_usageHandled within 2-72 hours"]  # Replace with your variable name

# Reverse the scaling
unscaled_coeff <- fixed_effect_coeff * overall_sd + overall_mean

# Reverse the Box-Cox transformation
original_scale_coeff <- (0.5 * unscaled_coeff + 1)^2

# The fixed effect coefficient in the original metric scale
original_scale_coeff





#---- Pairwise Mixed Effects model with everything ----
#----FIMM and Oslo---- 
ensrink_fimm <- rbind(fimm_org_response_common_drugs, subset(dss_enserink_common_drugs,select=-c(conc, AUC, DSS1, DSS3, IC50, believe_DSS, ID)))
ensrink_fimm <- inner_join(ensrink_fimm, experimental_var, by = 'lab')
list_var <- colnames(experimental_var)
ensrink_fimm <- merge(subset(fimm_sample_information, select = c("Patient.num", "Medium")), ensrink_fimm, by = "Patient.num",  all.y = TRUE)
ensrink_fimm <- ensrink_fimm %>%
  mutate(medium = if_else(lab == 'FIMM', Medium, medium))

ensrink_fimm_filtered <- ensrink_fimm[, sapply(ensrink_fimm, function(col) length(unique(col)) > 1)]
print(ensrink_fimm_filtered)
same_value_columns_enserink_fimm <- names(ensrink_fimm)[sapply(ensrink_fimm, function(col) length(unique(col)) == 1)]
print(same_value_columns_enserink_fimm)
colnames(ensrink_fimm_filtered)

enserink_fimm_mixed_mode_all_vars <- lmer(DSS2 ~ medium + microenvironmental_stimuli + (1|drug) , ensrink_fimm_filtered)
visualize(enserink_fimm_mixed_mode_all_vars)
summary(enserink_fimm_mixed_mode_all_vars)
ranef(enserink_fimm_mixed_mode_all_vars)
icc(enserink_fimm_mixed_mode_all_vars)
dotplot(ranef(enserink_fimm_mixed_mode_all_vars, condVar=TRUE))
REsim(enserink_fimm_mixed_mode_all_vars)
plotREsim(REsim(enserink_fimm_mixed_mode_all_vars)) 
View(ensrink_fimm_filtered)

#All vars
model_interpretation(ensrink_fimm_filtered, "DSS2", enserink_fimm_mixed_mode_all_vars, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-Oslo/lab_vs_dss2.html",header_1 = "FIMM Enserink All samples")

#Only fresh samples
enserink_fimm_mixed_fresh_only <- lmer(DSS2 ~ medium + microenvironmental_stimuli + (1|drug) , subset(ensrink_fimm_filtered, sample.x == 'fresh'))
model_interpretation(subset(ensrink_fimm_filtered, sample.x == 'fresh'), "DSS2", enserink_fimm_mixed_fresh_only, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-Oslo/lab_vs_dss2_fresh_only.html", header_1 = "FIMM Enserink Fresh samples")


enserink_fimm_mixed_fresh_only_test <- lmer(DSS2 ~ drug + medium + lab + (1|lab) , subset(ensrink_fimm_filtered, sample.x == 'fresh'))
model_interpretation(subset(ensrink_fimm_filtered, sample.x == 'fresh'), "DSS2", enserink_fimm_mixed_fresh_only_test, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-Oslo/lab_vs_dss2_fresh_only_test.html", header_1 = "FIMM Enserink Fresh samples - lab random effect")

#Polygeneic risk score patients only 



#----FIMM and Karolinska----
fimm_karolinska <- rbind(fimm_org_response_common_drugs, subset(dss_karolinska_common_drugs,select=-c(conc, AUC, DSS1, DSS3, IC50, believe_DSS, ID)))
fimm_karolinska <- inner_join(fimm_karolinska, experimental_var, by = 'lab')
list_var <- colnames(experimental_var)
fimm_karolinska <- merge(subset(fimm_sample_information, select = c("Patient.num", "Medium")), fimm_karolinska, by = "Patient.num",  all.y = TRUE)
fimm_karolinska <- fimm_karolinska %>%
  mutate(medium = if_else(lab == 'FIMM', Medium, medium))

fimm_karolinska_filtered <- fimm_karolinska[, sapply(fimm_karolinska, function(col) length(unique(col)) > 1)]
print(fimm_karolinska_filtered)
same_value_columns_fimm_karolinska <- names(fimm_karolinska)[sapply(fimm_karolinska, function(col) length(unique(col)) == 1)]
print(same_value_columns_fimm_karolinska)


fimm_karolinska_mixed_mode_all_vars <- lmer(DSS2 ~ medium + time_until_sample_usage  + microenvironmental_stimuli + (1|drug), fimm_karolinska_filtered)
visualize(fimm_karolinska_mixed_mode_all_vars)
summary(fimm_karolinska_mixed_mode_all_vars)
ranef(fimm_karolinska_mixed_mode_all_vars)
icc(fimm_karolinska_mixed_mode_all_vars)
dotplot(ranef(fimm_karolinska_mixed_mode_all_vars, condVar=TRUE))
REsim(fimm_karolinska_mixed_mode_all_vars)
plotREsim(REsim(fimm_karolinska_mixed_mode_all_vars)) 


model_interpretation(fimm_karolinska_filtered, "DSS2", fimm_karolinska_mixed_mode_all_vars, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-Karolinska/lab_vs_dss2.html", header_1 = "FIMM Karolinska All samples")
View(ensrink_fimm_filtered)

#Fresh only
fimm_karolinska_mixed_mode_fresh_only <- lmer(DSS2 ~ medium + time_until_sample_usage  + microenvironmental_stimuli + (1|drug), subset(fimm_karolinska_filtered, sample.x == 'fresh'))
model_interpretation(subset(fimm_karolinska_filtered, sample.x == 'fresh'), "DSS2", fimm_karolinska_mixed_mode_fresh_only, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-Karolinska/lab_vs_dss2_fresh_only.html", header_1 = "FIMM Karolinska Fresh samples")

#Frosen only
fimm_karolinska_mixed_mode_frozen_only <- lmer(DSS2 ~ medium + time_until_sample_usage  + microenvironmental_stimuli + (1|drug), subset(fimm_karolinska_filtered, sample.x == 'frozen'))
model_interpretation(subset(fimm_karolinska_filtered, sample.x == 'frozen'), "DSS2", fimm_karolinska_mixed_mode_frozen_only, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-Karolinska/lab_vs_dss2_frozen_only.html", header_1 = "FIMM Karolinska Frozen samples")


#----FIMM and BeatAML----
fimm_beataml <- rbind(fimm_org_response_common_drugs, subset(dss_beat_aml_common_drugs,select=-c(conc, AUC, DSS1, DSS3, IC50, believe_DSS, ID)))
fimm_beataml <- inner_join(fimm_beataml, experimental_var, by = 'lab')
list_var <- colnames(experimental_var)

fimm_beataml_filtered <- fimm_beataml[, sapply(fimm_beataml, function(col) length(unique(col)) > 1)]
print(fimm_beataml_filtered)
same_value_columns_fimm_beataml <- names(fimm_beataml)[sapply(fimm_beataml, function(col) length(unique(col)) == 1)]
print(same_value_columns_fimm_beataml)


fimm_beataml_mixed_mode_all_vars <- lmer(DSS2 ~ medium + sensitivity_readout  + positive_control + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug), fimm_beataml_filtered)
visualize(fimm_beataml_mixed_mode_all_vars)
summary(fimm_beataml_mixed_mode_all_vars)
ranef(fimm_beataml_mixed_mode_all_vars)
icc(fimm_beataml_mixed_mode_all_vars)
dotplot(ranef(fimm_beataml_mixed_mode_all_vars, condVar=TRUE))
REsim(fimm_beataml_mixed_mode_all_vars)
plotREsim(REsim(fimm_beataml_mixed_mode_all_vars)) 


#----Oslo and Karolinska----
ensrink_karolinska <- rbind(dss_enserink_common_drugs, dss_karolinska_common_drugs)
ensrink_karolinska <- inner_join(ensrink_karolinska, experimental_var, by = 'lab')
list_var <- colnames(experimental_var)

ensrink_karolinska_filtered <- ensrink_karolinska[, sapply(ensrink_karolinska, function(col) length(unique(col)) > 1)]
print(ensrink_karolinska_filtered)
same_value_columns_enserink_karolinska <- names(ensrink_karolinska)[sapply(ensrink_karolinska, function(col) length(unique(col)) == 1)]
print(same_value_columns_enserink_karolinska)


enserink_karolinska_mixed_mode_all_vars <- lmer(DSS2 ~ medium + time_until_sample_usage + microenvironmental_stimuli + cells + (1|drug), ensrink_karolinska_filtered)
visualize(enserink_karolinska_mixed_mode_all_vars)
summary(enserink_karolinska_mixed_mode_all_vars)
ranef(enserink_karolinska_mixed_mode_all_vars)
icc(enserink_karolinska_mixed_mode_all_vars)
dotplot(ranef(enserink_karolinska_mixed_mode_all_vars, condVar=TRUE))
REsim(enserink_karolinska_mixed_mode_all_vars)
plotREsim(REsim(enserink_karolinska_mixed_mode_all_vars))  

#All samples
model_interpretation(ensrink_karolinska_filtered, "DSS2", enserink_karolinska_mixed_mode_all_vars, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/Karolinska-Oslo/lab_vs_dss2.html", header_1 = "Enserink - Karolinska All samples")
View(ensrink_karolinska_filtered)

#Fresh only 
enserink_karolinska_mixed_mode_fresh_only <- lmer(DSS2 ~ medium + time_until_sample_usage + microenvironmental_stimuli + cells + (1|drug), subset(ensrink_karolinska_filtered, sample.x == 'fresh'))
model_interpretation(subset(ensrink_karolinska_filtered, sample.x == 'fresh'), "DSS2", enserink_karolinska_mixed_mode_fresh_only, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/Karolinska-Oslo/lab_vs_dss2_fresh_only.html", header_1 = "Enserink - Karolinska Fresh samples")

#test lab as random effect 
enserink_karolinska_mixed_mode_all_vars_test <- lmer(DSS2 ~ drug + lab + (1|lab), ensrink_karolinska_filtered)
model_interpretation(ensrink_karolinska_filtered, "DSS2", enserink_karolinska_mixed_mode_all_vars_test, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/Karolinska-Oslo/lab_vs_dss2_test.html", header_1 = "Enserink - Karolinska All samples - lab as random effect")


#----Oslo and BeatAML ----
ensrink_beataml <- rbind(dss_enserink_common_drugs, dss_beat_aml_common_drugs)
ensrink_beataml <- inner_join(ensrink_beataml, experimental_var, by = 'lab')
list_var <- colnames(ensrink_beataml)

ensrink_beataml_filtered <- ensrink_beataml[, sapply(ensrink_beataml, function(col) length(unique(col)) > 1)]
print(ensrink_beataml_filtered)
same_value_columns_ensrink_beataml <- names(ensrink_beataml)[sapply(ensrink_beataml, function(col) length(unique(col)) == 1)]
print(same_value_columns_ensrink_beataml)


ensrink_beataml_mixed_mode_all_vars <- lmer(DSS2 ~ medium + positive_control + nr_of_concentration_points + sensitivity_readout + (1|drug), ensrink_beataml_filtered)
visualize(ensrink_beataml_mixed_mode_all_vars)
summary(ensrink_beataml_mixed_mode_all_vars)
ranef(ensrink_beataml_mixed_mode_all_vars)
icc(ensrink_beataml_mixed_mode_all_vars)
dotplot(ranef(ensrink_beataml_mixed_mode_all_vars, condVar=TRUE))
REsim(ensrink_beataml_mixed_mode_all_vars)
plotREsim(REsim(ensrink_beataml_mixed_mode_all_vars)) 

#All vars
model_interpretation(ensrink_beataml_filtered, "DSS2", ensrink_beataml_mixed_mode_all_vars, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/Oslo-BeatAML/lab_vs_dss2.html", header_1 = "Enserink - BeatAML All samples")
ensrink_beataml_mixed_mode_lab <- lmer(DSS2 ~ lab + (1|drug), ensrink_beataml_filtered)
visualize(ensrink_beataml_mixed_mode_lab)
summary(ensrink_beataml_mixed_mode_lab)
as.data.frame(coef(summary(ensrink_beataml_mixed_mode_lab)))



#----Karolinska and BeatAML---- 
karolinska_beataml <- rbind(dss_karolinska_common_drugs, dss_beat_aml_common_drugs)
karolinska_beataml <- inner_join(karolinska_beataml, experimental_var, by = 'lab')
list_var <- colnames(experimental_var)

karolinska_beataml_filtered <- karolinska_beataml[, sapply(karolinska_beataml, function(col) length(unique(col)) > 1)]
print(karolinska_beataml_filtered)
same_value_columns_karolinska_beataml <- names(karolinska_beataml)[sapply(karolinska_beataml, function(col) length(unique(col)) == 1)]
print(same_value_columns_karolinska_beataml)

karolinska_beataml_filtered$Binary_lab <- ifelse(karolinska_beataml_filtered$lab == "Karolinska", 1, 0)

karolinska_beataml_mixed_mode_all_vars <- lmer(DSS2 ~ medium + sensitivity_readout + positive_control + time_until_sample_usage + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug) , karolinska_beataml_filtered)
visualize(karolinska_beataml_mixed_mode_all_vars)
summary(karolinska_beataml_mixed_mode_all_vars)
ranef(karolinska_beataml_mixed_mode_all_vars)
icc(karolinska_beataml_mixed_mode_all_vars)
dotplot(ranef(karolinska_beataml_mixed_mode_all_vars, condVar=TRUE))
REsim(karolinska_beataml_mixed_mode_all_vars)
plotREsim(REsim(karolinska_beataml_mixed_mode_all_vars)) 

#Binary just lab vs DSS2
#Assesing ransom var 
# baseline model glm
lm_check_re = lm(Binary_lab ~ 1, data = karolinska_beataml_filtered) 
# base-line mixed-model
lmer_check_re = lmer(Binary_lab ~ (1|drug), data = karolinska_beataml_filtered) 
aic.lmer <- AIC(logLik(lmer_check_re))
aic.lm <- AIC(logLik(lm_check_re))
aic.lmer; aic.lm  #Smaller AIC so including the random effect is justified
null.id = -2 * logLik(lm_check_re) + 2 * logLik(lmer_check_re)
pchisq(as.numeric(null.id), df=1, lower.tail=F) #P-value suggests random effect is varranted

karolinska_beataml_mixed <- lmer(DSS2 ~ Binary_lab + (1|drug) , karolinska_beataml_filtered)
karolinska_beataml_mixed_random_slope <- lmer(DSS2 ~ Binary_lab + (1 + Binary_lab|drug) , karolinska_beataml_filtered)
anova(karolinska_beataml_mixed_random_slope, karolinska_beataml_mixed)

karolinska_beataml_mixed_1 <- lmer(DSS2 ~ 1 + Binary_lab + (1|drug) , karolinska_beataml_filtered)
summary(karolinska_beataml_mixed)
summary(karolinska_beataml_mixed_1)
visualize(karolinska_beataml_mixed)
qqplot(fitted(karolinska_beataml_mixed), residuals(karolinska_beataml_mixed))

model_interpretation(karolinska_beataml_filtered, "DSS2", karolinska_beataml_mixed_mode_all_vars, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/Karolinska-BeatAML/Karonlinska_beataml_linear_mixed_model.html", header_1 = 'Karolinska - BeatAML All samples')
model_interpretation(karolinska_beataml_filtered, "DSS2", karolinska_beataml_mixed_random_slope, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/Karolinska-BeatAML/Karonlinska_beataml_linear_mixed_model_random_slope.html",  header_1 = 'Karolinska - BeatAML All samples - random slope')

#Fresh only 
karolinska_beataml_mixed_mode_fresh_only <- lmer(DSS2 ~ medium + sensitivity_readout + positive_control + time_until_sample_usage + microenvironmental_stimuli + cells + nr_of_concentration_points + (1|drug) , subset(karolinska_beataml_filtered, sample.x == 'fresh'))
model_interpretation(subset(karolinska_beataml_filtered, sample.x == 'fresh'), "DSS2", karolinska_beataml_mixed_random_slope, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/Karolinska-BeatAML/Karonlinska_beataml_fresh_only.html", header_1 = 'Karolinska - BeatAML Fresh samples')


#Combat values for each lab 
all_combatch_mixed_model <- lmer(combatch ~ DSS2 + lab + (1 + lab|drug) , combatch_all_datasets)
visualize(all_combatch_mixed_model)
model_interpretation(combatch_all_datasets, "combatch", all_combatch_mixed_model, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/all_labs_combatch.html", header_1 = 'Combat vs DSS2')


#FIMM vs Oslo 
combatch_mixed_model_FIMM_Enaerink <- lmer(combatch ~ DSS2 + medium + lab + (1 + lab|drug) , subset(combatch_all_datasets, lab == 'FIMM' | lab == 'Enserink'))
visualize(combatch_mixed_model_FIMM_Enaerink)
model_interpretation(subset(combatch_all_datasets, lab == 'FIMM' | lab == 'Enserink'), "combatch", combatch_mixed_model_FIMM_Enaerink, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-Oslo/FIMM_oslo_combatch.html", header_1 = 'Combat vs DSS2 - FIMM vs Oslo')

#FIMM vs Karolinska
combatch_mixed_model_FIMM_Karolinska <- lmer(combatch ~ DSS2 + medium + lab + (1 + lab|drug) , subset(combatch_all_datasets, lab == 'FIMM' | lab == 'Karolinska'))
visualize(combatch_mixed_model_FIMM_Karolinska)
model_interpretation(subset(combatch_all_datasets, lab == 'FIMM' | lab == 'Karolinska'), "combatch", combatch_mixed_model_FIMM_Karolinska, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-Karolinska/FIMM_karolinska_combatch.html", header_1 = 'Combat vs DSS2 - FIMM vs Karolinska')

#FIMM vs BeatAML
combatch_mixed_model_FIMM_BeatAML <- lmer(combatch ~ DSS2 + lab + (1 + lab|drug) , subset(combatch_all_datasets, lab == 'FIMM' | lab == 'BeatAML'))
visualize(combatch_mixed_model_FIMM_BeatAML)
model_interpretation(subset(combatch_all_datasets, lab == 'FIMM' | lab == 'BeatAML'), "combatch", combatch_mixed_model_FIMM_BeatAML, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/FIMM-BeatAML/FIMM_BeatAML_combatch.html", header_1 = 'Combat vs DSS2 - FIMM vs BeatAML')

#Oslo vs Karolinska
combatch_mixed_model_Enaerink_Karolinska <- lmer(combatch ~ DSS2 + lab + (1 + lab|drug) , subset(combatch_all_datasets, lab == 'Karolinska' | lab == 'Enserink'))
visualize(combatch_mixed_model_Enaerink_Karolinska)
model_interpretation(subset(combatch_all_datasets, lab == 'Karolinska' | lab == 'Enserink'), "combatch", combatch_mixed_model_Enaerink_Karolinska, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/Karolinska-Oslo/Karolinska_oslo_combatch.html", header_1 = 'Combat vs DSS2 - Karolinska vs Oslo')

#Oslo vs BeatAML
combatch_mixed_model_Enaerink_BeatAML <- lmer(combatch ~ DSS2 + lab + (1 + lab|drug) , subset(combatch_all_datasets, lab == 'BeatAML' | lab == 'Enserink'))
visualize(combatch_mixed_model_Enaerink_BeatAML)
model_interpretation(subset(combatch_all_datasets, lab == 'BeatAML' | lab == 'Enserink'), "combatch", combatch_mixed_model_Enaerink_BeatAML, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/Oslo-BeatAML/BeatAML_oslo_combatch.html", header_1 = 'Combat vs DSS2 - BeatAML vs Oslo')

#Karolinska vs BeatAML
combatch_mixed_model_Karolinska_BeatAML <- lmer(combatch ~ DSS2 + lab + (1 + lab|drug) , subset(combatch_all_datasets, lab == 'BeatAML' | lab == 'Karolinska'))
visualize(combatch_mixed_model_Karolinska_BeatAML)
model_interpretation(subset(combatch_all_datasets, lab == 'BeatAML' | lab == 'Karolinska'), "combatch", combatch_mixed_model_Karolinska_BeatAML, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/Karolinska-BeatAML/BeatAML_oslo_combatch.html", header_1 = 'Combat vs DSS2 - BeatAML vs Karolinska')



#DSS vs lab models
dss_vs_lab_models <- data.frame(Model = character(),
                                Fixed_Effect_Term = character(),
                                Fixed_Effect_Coefficient = numeric(),
                                stringsAsFactors = FALSE)


# List of models and their names
models <- list(FIMM_Enserink_all_samples = enserink_fimm_mixed_mode_all_vars, FIMM_Enserink_fresh_samples = enserink_fimm_mixed_fresh_only, 
               FIMM_Karolinska_all_samples = fimm_karolinska_mixed_mode_all_vars, FIMM_Karolinska_fresh_samples = fimm_karolinska_mixed_mode_fresh_only, FIMM_Karolinska_froosen_samples = fimm_karolinska_mixed_mode_frozen_only,
               FIMM_BeatAML_all_samples = fimm_beataml_mixed_mode_all_vars, FIMM_BeatAML_fresh_samples = fimm_beataml_mixed_mode_fresh_only, FIMM_BeatAML_blood_and_fresh_samples = fimm_beataml_mixed_mode_blood_and_fresh,
               FIMM_BeatAML_bone_marrow_and_fresh = fimm_beataml_mixed_mode_bone_marrow_and_fresh, FIMM_BeatAML_relapse_and_bone_marrow_and_fresh_samples = fimm_beataml_mixed_mode_relapse_and_bone_marrow_and_fresh,
               FIMM_BeatAML_blood_samples = fimm_beataml_mixed_mode_blood_only, FIMM_BeatAML_diagnosis_and_blood_samples = fimm_beataml_mixed_mode_diagnosis_and_blood, FIMM_BeatAML_bone_marrow_samples = fimm_beataml_mixed_mode_bone_marrow_only,
               FIMM_BeatAML_diagnosis_and_bone_marrow_samples = fimm_beataml_mixed_mode_diagnosis_and_bone_marrow, FIMM_BeatAML_diagnosis_samples = fimm_beataml_mixed_mode_Diagnosis_only, FIMM_BeatAML_refractory_samples = fimm_beataml_mixed_mode_Refractory_only,
               FIMM_BeatAML_relapse_samples = fimm_beataml_mixed_mode_relapse_only, FIMM_BeatAML_relapse_and_blood_samples = fimm_beataml_mixed_mode_relapse_and_blood, 
               Karolinska_BeatAML_all_samples = karolinska_beataml_mixed_mode_all_vars, Karolinska_BeatAML_fresh_samples = karolinska_beataml_mixed_mode_fresh_only, 
               Enserink_Karolinska_all_samples = enserink_karolinska_mixed_mode_all_vars, Enserink_Karolinska_fresh_samples = enserink_karolinska_mixed_mode_fresh_only)


df <- data.frame(
  Model = character(),
  Fixed_Effect_Term = character(),
  Fixed_Effect_Coefficient = numeric(),
  stringsAsFactors = FALSE
)

# Step 3: Extract fixed effects and their coefficients
for (model_name in names(models)) {
  model <- models[[model_name]]
  fixed_effects <- fixef(model)  # Get fixed effects
  fixed_effect_names <- names(fixed_effects)  # Get names of the fixed effects
  
  # Create a temporary data frame for this model's fixed effects
  temp_df <- data.frame(
    Model = model_name,
    Fixed_Effect_Term = fixed_effect_names,
    Fixed_Effect_Coefficient = fixed_effects,
    stringsAsFactors = FALSE
  )
  
  # Step 4: Add the temporary data frame to the main data frame
  df <- rbind(df, temp_df)
}

# Display the final data frame
print(df)

p <- ggplot(df, aes(x = Model, y = Fixed_Effect_Coefficient, fill = Fixed_Effect_Term)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Fixed Effects and Their Coefficients",
    x = "Fixed Effect Term",
    y = "Coefficient"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # Rotate and size x-axis labels
    plot.margin = margin(1, 1, 2, 1, "cm")  # Increase plot margins
  )

# Save the plot to a file
ggsave("~/Desktop/fixed_effects_plot.png", plot = p, width = 20, height = 30, dpi = 300)











#Using linear mixed effects model to look at variance of parients and drugs on the dataset

test_data <- all_datasets %>% filter(drug %in% drug_names_to_include)
subset(test_data, lab == 'Enserink')

plot(test_data)

qplot(test_data$drug, test_data$DSS2, xlab='DSS2', ylab='Drug')

flexplot(lab~1, data=test_data)
flexplot(drug~1, data=test_data)
unique(test_data$nr_negative_control)
model_lab<-lm(DSS2~drug+as.factor(lab), test_data)
summary(model_lab)
ggplotRegression(model_lab)

model_negative<-lm(DSS2~drug+as.factor(nr_negative_control), test_data)
summary(model_negative)
ggplotRegression(model_negative)

aov_out = ezANOVA(data = test_data, dv = .(DSS2), wid= .(lab), within = .(drug), detailed=TRUE, type=3)
ezDesign(data = test_data, medium, drug_sensitivity_readout, lab)



#Run baseline model (will use it later for comparison, we have  no controls)
model_null <- lmer(DSS2~1+(1|lab), # note how we use 1 to suggest that that we keep the slope constant and vary intercept
                   data=test_data)
visualize(model_null, plot="model")
icc(model_null)
dotplot(ranef(model_null, condVar=TRUE))

model_null_drug <- lmer(DSS2~1+(1|drug), # note how we use 1 to suggest that that we keep the slope constant and vary intercept
                        data=test_data)
visualize(model_null_drug)
icc(model_null_drug)
dotplot(ranef(model_null_drug, condVar=TRUE))

full <- lmer(DSS2~drug + (drug|lab), test_data)
reduced <- lmer(DSS2~drug + (1|lab), test_data)
visualize(full)
dotplot(ranef(full, condVar=TRUE))
icc(full)
icc(reduced)
model.comparison(full, reduced)

#Run the model (note thay we also control for the effect of time by subject)
model_mix <- lmer(DSS2 ~ drug  + medium + (1|lab), test_data)
visualize(model_mix, plot="model")
summary(model_mix)
icc(model_mix)
dotplot(ranef(model_mix, condVar=TRUE))


model.comparison(model_null, model_mix)
compare.fits(DSS2~drug | lab, data = test_data, model_null, model_mix)

#Summarise
summary(model_mix)
display(model_mix)
fixef(model_mix) 
se.fixef(model_mix)
ranef(model_mix)
se.ranef(model_mix)
coef(model_mix)
linearity<-plot(resid(model_mix),#extract the residuals
                test_data$DSS2)

qqmath(model_mix)

#Try with model_plot (argument for type can be varied)
plot_model(model_mix, type='diag')
plot_model(model_mix)


flexplot(DSS2~drug | lab, data = test_data, ghost.line = "red")
m <- flexplot(DSS2~medium | drug, data = test_data, ghost.line = "red")
ds <- flexplot(DSS2~drug_sensitivity_readout | drug, data = test_data, ghost.line = "red")


f <- DSS2 ~ drug + as.factor(medium) + as.factor(drug_sensitivity_readout) + as.factor(nr_positive_control) + as.factor(positive_control) + as.factor(spread_of_controls_on_plate) + as.factor(time_until_sample_usage) + as.factor(centrifugation_procedure) + as.factor(microenvironmental_stimuli) + as.factor(cell_culturing_conditions) + as.factor(cells) + as.factor(nr_of_concentration_points) 
m <- model.matrix(f, test_data)
x <- lmer(DSS2 ~ drug + as.factor(medium) + as.factor(drug_sensitivity_readout) + as.factor(nr_positive_control) + as.factor(positive_control) + as.factor(spread_of_controls_on_plate) + as.factor(time_until_sample_usage) + as.factor(centrifugation_procedure) + as.factor(microenvironmental_stimuli) + as.factor(cell_culturing_conditions) + as.factor(cells) + as.factor(nr_of_concentration_points) + (1|lab), test_data
)
summary(x)
model_mix_all_var <- lmer(DSS2 ~ drug + medium + sensitivity_readout + nr_positive_control + positive_control + time_until_sample_usage + centrifugation_procedure + microenvironmental_stimuli + cell_culturing_conditions + cells + nr_of_concentration_points + (1|lab) , test_data)
visualize(model_mix_all_var)
summary(model_mix_all_var)
print(summary(model_mix_all_var), correlation=TRUE)
ranef(model_mix_all_var)
icc(model_mix_all_var)
dotplot(ranef(model_mix_all_var, condVar=TRUE))
REsim(model_mix_all_var)
qqnorm(residuals(model_mix_all_var))
plot(model_mix_all_var)
plotREsim(REsim(model_mix_all_var))  



ggplot(data.frame(x1=test_data$medium,pearson=residuals(model_mix_all_var,type="pearson")),
       aes(x=x1,y=pearson)) +
  geom_point() +
  theme_bw()


df_factor <- as.data.frame(lapply(test_data, as.factor))

ggplot(data.frame(lev=hatvalues(model_mix_all_var),pearson=residuals(model_mix_all_var,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw()

levId <- which(hatvalues(model_mix_all_var) >= .172)
test_data[levId,c("DSS2","medium","drug_sensitivity_readout","lab")]
summary(test_data[,c("DSS2","medium","drug_sensitivity_readout")])

mmLev <- lmer(DSS2 ~ medium + drug_sensitivity_readout + (1|lab), data=test_data[-c(levId),])
mmLevCD <- data.frame(effect=fixef(model_mix_all_var),
                      change=(fixef(mmLev) - fixef(model_mix_all_var)),
                      se=sqrt(diag(vcov(model_mix_all_var)))
)
rownames(mmLevCD) <- names(fixef(mmLev))
mmLevCD$multiples <- abs(mmLevCD$change / mmLevCD$se)
mmLevCD

qqmath(ranef(model_mix_all_var)$lab[,1], main = "Q-Q Plot of Random Effects")

influencePlot(model_mix_all_var, main = "Residuals vs Leverage")


model.comparison(model_mix, model_mix_all_var)
anova(model_mix, model_mix_all_var)
compare.fits(DSS2~drug | lab, data = test_data, model_mix, model_mix_all_var, re=F)


visualize(model_mix_all_var_interactions)
model.comparison(model_mix_all_var_interactions, model_mix_all_var)



plot(model_mix_all_var)
qqline(resid(model_mix_all_var))

test_data$predicted <- predict(model_mix_all_var)

prelim_plot <- ggplot(test_data, aes(x = drug, y = DSS2, color=lab)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm")

mm_plot <- ggplot(test_data, aes(x = drug, y = DSS2, colour = lab)) +
  facet_wrap(~medium, nrow=2) +   # a panel for each mountain range
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_line(aes(y = predicted), size = 1) +  # adding predicted line from mixed model 
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"))

ggplot(test_data, aes(x = predicted, y = DSS2, color = lab)) +
  geom_point() +
  labs(title = "Line plot",
       x = "Predicted",
       y = "Average DSS score per drug within each quantile",
       color = "Lab") +
  theme_minimal()

ggplotRegression(lm(DSS2~predicted,data=test_data))

m3_predict <- ggpredict(model_mix_all_var, terms = c("drug"))

ggplot(m3_predict) + 
  geom_line(aes(x = x, y = predicted)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = result,                      # adding the raw data (scaled values)
             aes(x = drug, y = AVG_score_per_quantile, colour = lab)) + 
  labs(x = "Body Length (indexed)", y = "Test Score", 
       title = "Body length does not affect intelligence in dragons") + 
  theme_minimal()

ggpredict(model_mix_all_var, terms = c("drug", "lab"), type = "re") %>% 
  plot() +
  labs(x = "", y = "", title = "") + 
  theme_minimal()



#Comparing combatch with original DSS
model_simple<-lm(DSS2~combatch+as.factor(lab), combatch_all_datasets)
summary(model_simple)
ggplotRegression(model_simple)
log(combatch_all_datasets$DSS2)

test <- lmer(combatch ~ 1 + (1|lab) , combatch_all_datasets)
visualize(test)
icc(test)
dotplot(ranef(test, condVar=TRUE))

combatch_model <- lmer(DSS2 ~ combatch + (1|drug), combatch_all_datasets)
visualize(combatch_model)
dotplot(ranef(combatch_model, condVar=T))


# Calculate VIF
vif(model_mix_all_var)


df_factor <- as.data.frame(lapply(test_data, as.factor))
pca_vars <- df_factor[, c("drug", "medium", "drug_sensitivity_readout", "nr_positive_control","positive_control", "spread_of_controls_on_plate", "time_until_sample_usage", "centrifugation_procedure", "microenvironmental_stimuli", "cell_culturing_conditions", "cells", "nr_of_concentration_points")]
# Install and load FactoMineR package

# Perform MCA
mca_result <- MCA(pca_vars, graph = TRUE) 
summary(mca_result)
plot.MCA(mca_result, invisible = "ind")
mca_dimensions <- as.data.frame(mca_result$ind$coord)  # Extract coordinates of individuals (rows)
lmer_model <- lmer(DSS2 ~ . + (1 | lab), data = cbind(mca_dimensions, test_data))
summary(lmer_model)

rf_data <- cbind(mca_dimensions, test_data)
rf_model <- randomForest(DSS2 ~ ., data = rf_data)
# Extract and plot a single decision tree from Random Forest
single_tree <- rf_model$forest[[1]]  # Extract the first tree (you can choose any index)

# Plot the decision tree
rpart.plot(single_tree, type = 2, extra = 1, fallen.leaves = TRUE, tweak = 1.2)

model_matrix <- model.matrix(~ drug + medium + drug_sensitivity_readout + nr_positive_control + positive_control + spread_of_controls_on_plate + time_until_sample_usage + centrifugation_procedure + microenvironmental_stimuli + cell_culturing_conditions + cells + nr_of_concentration_points, data = test_data)


library(pls)
plsr_model <- plsr(DSS2 ~ drug + medium + drug_sensitivity_readout + positive_control + spread_of_controls_on_plate + time_until_sample_usage + centrifugation_procedure + microenvironmental_stimuli + cell_culturing_conditions + cells + nr_of_concentration_points, data = test_data, ncomp = 2)
summary(plsr_model)


library(randomForest)
rf_model <- randomForest(DSS2 ~ drug + lab + medium + positive_control + time_until_sample_usage + centrifugation_procedure + microenvironmental_stimuli + cell_culturing_conditions + cells + nr_of_concentration_points, data = all_datasets_v2)
varImpPlot(rf_model)






test_data$DSS2Unscaled <- test_data$DSS2
test_data$DSS2Scaled <- scale(test_data$DSS2, scale = F)
model_mix_all_var <- lmer(DSS2 ~ medium + positive_control + time_until_sample_usage + centrifugation_procedure + microenvironmental_stimuli + cell_culturing_conditions + cells + nr_of_concentration_points + (1|drug) , test_data)
model_interpretation(test_data, "DSS2", model_mix_all_var, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/all_var_all_labs.html")
model_mix_lab <- lmer(DSS2 ~ lab + (1|drug) , test_data)
model_interpretation(test_data, "DSS2", model_mix_lab, "/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Mixed models/lab_vs_dss2.html")



summary(model_mix_all_var)

visualize(model_mix_all_var)

init_plots(test_data, "DSS2", "drug")
step(model_mix_all_var, direction = "both")

l_model <- lm(DSS2 ~ lab, test_data)

autoplot(l_model) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw()

infl <- influence(model_mix_all_var, group = "drug")
cutoff <- 4/((nrow(test_data)-length(model_mix_all_var)-2))
# start plotting
par(mfrow = c(1, 2))           # display plots in 3 rows/2 columns
qqPlot(model_mix_all_var, main="QQ Plot")

flexplot(model_mix_all_var, which=4, cook.levels = cutoff); par(mfrow = c(1, 1))

plot(model_mix_all_var, lab ~ resid(.), abline = 0 ) # generate diagnostic plots
plot(model_mix_all_var, resid(., type = "pearson") ~ fitted(.) | lab, id = 0.05, 
     adj = -0.3, pch = 20, col = "gray40")

leveneTest(test_data$DSS2, test_data$drug, center = mean)

plot(model_mix_all_var, DSS2 ~ fitted(.), id = 0.05, adj = -0.3, cex = .8, pch = 20, col = "blue")

acf(residuals(model_mix_all_var))


anova(model_mix_all_var, l_model)




# lab as fixed effect and drug as random effect
# want random effect to have at least 5 levels usually - lab for pairwise comparison is just two
# remember here lab represents the experimental set up differences between those two labs - this is interchangable


