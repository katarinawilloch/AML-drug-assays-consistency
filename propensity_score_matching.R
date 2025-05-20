#Read libraries
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(stats)
library(ggplot2)
library(nnet)

#Oslo clinical info 
oslo_clinical_info <- read_csv('~/Desktop/UiO/Project 1/Data/Second cleansing/Oslo_clinical_info.csv')
oslo_clinical_info$lab <- 'Enserink'

#FIMM clinical infor
fimm_clinical_info <- read_csv('~/Desktop/UiO/Project 1/Data/Second cleansing/FIMM_clinical_info.csv')
fimm_clinical_info$lab <- 'FIMM'

#BeatAML clinical info 
beatAML_clinical_info <- read_csv('~/Desktop/UiO/Project 1/Data/Second cleansing/BeatAML_clinical_info.csv')
beatAML_clinical_info <- beatAML_clinical_info %>% mutate(X = paste0(beatAML_clinical_info$dbgap_subject_id, '_', beatAML_clinical_info$dbgap_sample_id, '_', beatAML_clinical_info$dbgap_rnaseq_sample))
beatAML_clinical_info$X
beatAML_clinical_info$lab <- 'Beat AML'

#bind the two clinical datasets
combined_clinical_info <- rbind(oslo_clinical_info[,c('X', 'mutations', 'FAB_class', 'sex', 'karyotype_class', 'lab')], fimm_clinical_info[,c('X', 'mutations', 'FAB_class', 'sex', 'karyotype_class', 'lab')])

#adding beatAML
combined_clinical_info_oslo_fimm_beataml <- rbind(combined_clinical_info, beatAML_clinical_info[,c('X', 'mutations', 'FAB_class', 'sex', 'karyotype_class', 'lab')])

#preparing dataset 
combined_clinical_info_oslo_fimm_beataml$mutations <- ifelse(is.na(combined_clinical_info_oslo_fimm_beataml$mutations), "None", combined_clinical_info_oslo_fimm_beataml$mutations)
combined_clinical_info_oslo_fimm_beataml$FAB_class <- ifelse(is.na(combined_clinical_info_oslo_fimm_beataml$FAB_class), "None", combined_clinical_info_oslo_fimm_beataml$FAB_class)
combined_clinical_info_oslo_fimm_beataml$karyotype_class <- ifelse(is.na(combined_clinical_info_oslo_fimm_beataml$karyotype_class), "None", combined_clinical_info_oslo_fimm_beataml$karyotype_class)


#Read DSS info for each dataset
#oslo
oslo_response_scores <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink_full_drug_set.csv')
#fimm
fimm_response_score <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set.csv')

fimm_response_score_1 <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_drug_response.csv')
fimm_response_score_1 <- gather(fimm_response_score_1, Patient.num, DSS2, 'AML_084_04':'Healthy_17', factor_key=TRUE)
fimm_response_score_1$...1 <- NULL
fimm_response_score_1$Drug_ID <- NULL
colnames(fimm_response_score_1)[colnames(fimm_response_score_1) == "Drug_name"] <- "drug"

#beataml
BeatAML_response_scores <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_beat_aml_full_drug_set.csv')
BeatAML_response_scores$lab <- 'Beat AML'
BeatAML_response_scores$X <- BeatAML_response_scores$Patient.num

#Bind the response scores for the two datasets
combined_response_score <- rbind(oslo_response_scores[,c('Patient.num', 'drug', 'DSS2')], fimm_response_score_1)
combined_response_score$X <- gsub("Patient ", "",combined_response_score$Patient.num)
combined_response_score$X <- gsub("patient ", "",combined_response_score$X)
combined_response_score$X <- gsub("_r.*", "",combined_response_score$X)

remove_after_second_underscore <- function(x) {
  parts <- unlist(strsplit(x, "_"))
  if (length(parts) > 2) {
    return(paste(parts[1:2], collapse = "_"))
  } else {
    return(x)
  }
}

# Apply the function to the data frame
combined_response_score$X<- sapply(combined_response_score$X, remove_after_second_underscore)

#Merge clinical data and response scores 
all_combined <- merge(combined_response_score, combined_clinical_info, by="X", x.all=TRUE)


#Regress group with variables and DSS2
all_combined$Binary_lab <- ifelse(all_combined$lab == "FIMM", 1, 0)
all_combined$mutations <- ifelse(is.na(all_combined$mutations), "None", all_combined$mutations)
all_combined$FAB_class <- ifelse(is.na(all_combined$FAB_class), "None", all_combined$FAB_class)
all_combined$karyotype_class <- ifelse(is.na(all_combined$karyotype_class), "None", all_combined$karyotype_class)



oslo_fimm_df <- unique(all_combined[,c('X', 'lab', 'mutations', 'FAB_class', 'Binary_lab', 'sex')])
#df <- all_combined
# Function to expand multiple columns into binary columns
expand_columns <- function(df, column_names) {
  
  # Initialize a list to store binary data frames
  binary_dfs <- list()
  
  # Loop through each column name provided
  for (column_name in column_names) {
    
    # Split comma-separated values into a list of vectors
    values_list <- strsplit(trimws(as.character(df[[column_name]])), ",")
    
    # Get unique values from all rows
    unique_values <- unique(unlist(values_list))
    
    # Create a data frame with binary columns for each unique value
    binary_df <- as.data.frame(matrix(0, nrow = nrow(df), ncol = length(unique_values)))
    colnames(binary_df) <- unique_values
    
    # Fill the binary columns
    for (i in 1:nrow(df)) {
      present_values <- values_list[[i]]
      binary_df[i, present_values] <- 1
    }
    
    # Add binary data frame to the list
    binary_dfs[[column_name]] <- binary_df
  }
  
  # Combine all binary data frames with the original data frame
  result_df <- bind_cols(df, do.call(bind_cols, binary_dfs))
  
  return(result_df)
}

# Apply the function to the data frame
data_expanded <- expand_columns(oslo_fimm_df, c("mutations", "FAB_class"))

data_expanded_oslo_fimm_beataml <- data_expanded <- expand_columns(combined_clinical_info_oslo_fimm_beataml, c("mutations", "FAB_class"))

colnames(data_expanded)

unique(subset(all_combined, lab == 'Enserink', select = 'mutations'))

#Distribution FAB classes between labs

df_long <- combined_clinical_info_oslo_fimm_beataml %>%
  separate_rows(FAB_class, sep = ",") %>%
  group_by(lab, FAB_class) %>%
  summarise(count = n())%>%
  mutate(percent = round(count / sum(count) * 100, 1))

n_obs <- combined_clinical_info_oslo_fimm_beataml %>%
  group_by(lab) %>%
  summarise(n_obs = n(), .groups = 'drop')

# Join the number of observations with the counts
df_long <- df_long %>%
  left_join(n_obs, by = "lab") %>%
  mutate(mean_count = count / n_obs)

ggplot(df_long, aes(x = FAB_class, y = percent, fill = as.factor(lab))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "FAB classes",
       x = "FAB class",
       y = "count / number of observations",
       fill = "Lab") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Distribution mutations classes between labs

df_long <- combined_clinical_info_oslo_fimm_beataml %>%
  separate_rows(mutations, sep = ",") %>%
  group_by(lab, mutations) %>%
  summarise(count = n()) %>%
  mutate(percent = round(count / sum(count) * 100, 1))

n_obs <- combined_clinical_info %>%
  group_by(lab) %>%
  summarise(n_obs = n(), .groups = 'drop')

# Join the number of observations with the counts
df_long <- df_long %>%
  left_join(n_obs, by = "lab") %>%
  mutate(mean_count = count / n_obs)

ggplot(df_long, aes(x = mutations, y = percent, fill = as.factor(lab))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mutations",
       x = "Mutation",
       y = "count / number of observations",
       fill = "Lab") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Distribution sex 
df_long <- combined_clinical_info %>%
  separate_rows(sex, sep = ",") %>%
  group_by(lab, sex) %>%
  summarise(count = n()) %>%
  mutate(percent = round(count / sum(count) * 100, 1))


ggplot(df_long, aes(x = sex, y = percent, fill = as.factor(lab))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Gender",
       x = "Gender",
       y = "Percentage of Patients per lab",
       fill = "Lab") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Test for confounders to include in the propensity score matching ... 
#The variables does not always have to be related to both labs for them to be included 
#Always better to include as many confounders as possible

library(corrplot)
colnames(data_expanded)
cor_matrix <- cor(data_expanded[, !names(data_expanded) %in% c('DSS2','drug','Patrient.num','lab', 'X', 'FAB_class', 'mutations', 'sex', "karyotype_class", "Binary_lab")])
print(cor_matrix)

corrplot(cor_matrix, 
         col = colorRampPalette(c("blue", "white", "red"))(200), # Color gradient
         type = "upper", # Show only upper triangle
         order = "hclust", # Hierarchical clustering order
         addCoef.col = NULL, # Do not add correlation coefficients
         tl.col = "black", # Text label color
         tl.srt = 45) 
#Highly correlated mutations removed


mod_test1 <- glm(as.factor(lab) ~ mutations + FAB_class + sex , data=unique(all_combined[,c('X', 'mutations', 'FAB_class', 'Binary_lab', 'sex', 'lab')]), family="binomial")
mod_test2 <- lm(DSS2 ~ mutations + FAB_class + drug, data=all_combined)
mod_test_binary <- glm(Binary_lab ~ `M4 eos` + `M1/M2` + `M4/M5` + M3 + `None...31` + M2 + M4 + M0  + M6 + M7, 
                       data=unique(data_expanded), family="binomial")


summary(mod_test1)
summary(mod_test2)
summary(mod_test_binary)   

df <- unique(all_combined[,c('X', 'mutations', 'FAB_class', 'Binary_lab', 'sex', 'lab')])
df$propensity_score <- predict(mod_test1, type = "response")
dim(df)
dim(predict(mod_test1, type = "response"))

density_values <- density(df$propensity_score)
density_range <- range(density_values$y)

y_max <- density_range[2] * 1.1

ggplot(df, aes(x = propensity_score, color = as.factor(lab))) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Propensity Scores",
       x = "Propensity Score",
       fill = "Binary_lab") +
  ylim(0, y_max)


df <- unique(data_expanded)
df$propensity_score <- predict(mod_test_binary, type = "response")


density_values <- density(df$propensity_score)
density_range <- range(density_values$y)

y_max <- density_range[2] * 1.1

ggplot(df, aes(x = propensity_score, color = as.factor(Binary_lab))) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Propensity Scores",
       x = "Propensity Score",
       fill = "Binary_lab") +
  ylim(0, y_max)


df <- unique(data_expanded)
df$propensity_score <- predict(mod_test_binary, type = "response")
density_values <- density(df$propensity_score)
density_range <- range(density_values$y)

y_max <- density_range[2] * 1.1

ggplot(df, aes(x = propensity_score, color = as.factor(lab))) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Propensity Scores",
       x = "Propensity Score",
       fill = "Binary_lab") +
  ylim(0, y_max)


# Create histogram plot with overlapping bars
ggplot(df, aes(x = propensity_score, fill = as.factor(lab))) +
  geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.05) +
  labs(title = "Distribution of Propensity Scores",
       x = "Propensity Score",
       fill = "Binary_lab") +
  ylim(0, 50) +
  theme_minimal() +
  theme(legend.position = "top")


your_data <- unique(data_expanded) %>% 
  mutate(logit = log(predict(mod_test_binary, type = "response") / (1 - predict(mod_test_binary, type = "response"))))

ggplot(your_data, aes(x = lab, y = logit)) + 
  geom_point() + 
  geom_smooth(method = "loess")


alias(mod_test1)
alias(mode_test2)
alias(mod_test_binary)


library(car)
vif(mod_test1)
vif(mode_test2)
vif(mod_test_binary)

rownames(vif(mod_test_binary))
library(pheatmap)
gvif_matrix <- matrix(vif(mod_test_binary), nrow = 1, dimnames = list(NULL, names(vif(mod_test_binary))))
pheatmap(gvif_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = colorRampPalette(c("white", "red"))(100),
         main = "Heatmap of VIF Values",
         display_numbers = TRUE,
         fontsize_number = 8, # Adjust the size of the numbers
         fontsize_row = 12,   # Adjust the size of the row names
         fontsize_col = 12) 

influencePlot(mod_test_binary)
AIC(mod_test_binary)


#Propensity score analisis
library(MatchIt)

psa_n_mutations <- matchit(Binary_lab ~ NRAS + NPM1 + `FLT3-ITD` + EVI1 + KMT2A + `FLT3-TKD` + ASXL1 + CEBPA + del20 + RUNX1 + Not_available + WT1 + KIT + MLL_AF9 + NUP,
                 data=unique(data_expanded), distance = "glm", method = "nearest", m.order = "largest", replace = FALSE)
psa_n_mutations
summary(psa_n_mutations)

#optimal pair matching
psa_0_mutations <- matchit(Binary_lab ~ NRAS + NPM1 + `FLT3-ITD` + EVI1 + KMT2A + `FLT3-TKD` + ASXL1 + CEBPA + del20 + RUNX1 + Not_available + WT1 + KIT + MLL_AF9 + NUP,
                 data=unique(data_expanded), distance = "glm", method = "optimal", ratio = 2)
psa_0_mutations
summary(psa_0_mutations)


psa_f_mutations <- matchit(Binary_lab ~ NRAS + NPM1 + `FLT3-ITD` + EVI1 + KMT2A + `FLT3-TKD` + ASXL1 + CEBPA + del20 + RUNX1 + Not_available + WT1 + KIT + MLL_AF9 + NUP + sex,
                 data=unique(data_expanded), distance = "glm", method = "full", estimand = "ATE")
psa_f_mutations
summary(psa_f_mutations)


library(cobalt)
love.plot(bal.tab(psa_n_mutations), stats = c("m", "v"), grid = TRUE, thresholds = c(m=0.25, v=1.25))

love.plot(bal.tab(psa_0_mutations), stats = c("m", "v"), grid = TRUE, thresholds = c(m=0.25, v=1.25))

love.plot(bal.tab(psa_f_mutations), stats = c("m", "v"), grid = TRUE, thresholds = c(m=0.25, v=1.25))


psa_n_dat_mutations <- match.data(psa_f_mutations)

fimm_oslo_mathcing_mutations <- merge(psa_n_dat_mutations, subset(all_combined, select= c('drug', 'DSS2', 'X')), by='X')
fimm_post_match <- subset(fimm_oslo_mathcing_mutations, lab.fimm == 'FIMM')
oslo_post_match <- subset(fimm_oslo_mathcing_mutations, lab.fimm=='Enserink')
fimm_oslo_mathcing_mutations <- merge(fimm_post_match, oslo_post_match, by=c('subclass', 'drug'), suffixes = c(".fimm",".oslo"))

# Apply the function to each subset of data by 'drug'
model_stats <- fimm_oslo_mathcing_mutations %>%
  group_by(drug) %>%
  do(get_model_stats(.))

# Merge the statistics back into the original data
df <- fimm_oslo_mathcing_mutations %>%
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


#FAB
psa_n_fab <- matchit(Binary_lab ~ `M4 eos` + `M1/M2` + `M4/M5` + M3 + `None...7` + M2 + M4 + M0  + M6 + M7,
                           data=unique(data_expanded), distance = "glm", method = "nearest", m.order = "largest", replace = FALSE)
psa_n_fab
summary(psa_n_fab)

#optimal pair matching
psa_0_fab <- matchit(Binary_lab ~ `M4 eos` + `M1/M2` + `M4/M5` + M3 + `None...7` + M2 + M4 + M0  + M6 + M7,
                           data=unique(data_expanded), distance = "glm", method = "optimal", ratio = 2)
psa_0_fab
summary(psa_0_fab)


psa_f_fab <- matchit(Binary_lab ~ `M4 eos` + `M1/M2` + `M4/M5` + M3 + `None...7` + M2 + M4 + M0  + M6 + M7,
                           data=unique(data_expanded), distance = "glm", method = "full", estimand = "ATE")
psa_f_fab
summary(psa_f_fab)


library(cobalt)
love.plot(bal.tab(psa_n_fab), stats = c("m", "v"), grid = TRUE, thresholds = c(m=0.25, v=1.25))

love.plot(bal.tab(psa_0_fab), stats = c("m", "v"), grid = TRUE, thresholds = c(m=0.25, v=1.25))

love.plot(bal.tab(psa_f_fab), stats = c("m", "v"), grid = TRUE, thresholds = c(m=0.25, v=1.25))


psa_f_dat_fab <- match.data(psa_f_fab)

fimm_oslo_mathcing_fab <- merge(psa_f_dat_fab, subset(all_combined, select= c('drug', 'DSS2', 'X')), by='X')
fimm_post_match <- subset(fimm_oslo_mathcing_fab, lab == 'FIMM')
oslo_post_match <- subset(fimm_oslo_mathcing_fab, lab=='Enserink')
fimm_oslo_mathcing_fab <- merge(fimm_post_match, oslo_post_match, by=c('subclass', 'drug'), suffixes = c(".fimm",".oslo"))

# Apply the function to each subset of data by 'drug'
model_stats <- fimm_oslo_mathcing_fab %>%
  group_by(drug) %>%
  do(get_model_stats(.))

# Merge the statistics back into the original data
df <- fimm_oslo_mathcing_fab %>%
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

#trying with beataml as well 

multinom_fit <- multinom(lab ~ NRAS + NPM1 + `FLT3-ITD` + EVI1 + KMT2A + `FLT3-TKD` + ASXL1 + CEBPA + del20 + RUNX1 + Not_available + WT1 + KIT + MLL_AF9 + NUP + 
                           `M4 eos` + `M1/M2` + `M4/M5` + M3 + `None...7` + M2 + M4 + M0  + M6 + M7
                         + sex, 
                         data=unique(data_expanded_oslo_fimm_beataml))
summary(multinom_fit)

propensity_scores <- predict(multinom_fit, type = "probs")

matchit_out <- matchit(lab ~ NRAS + NPM1 + `FLT3-ITD` + EVI1 + KMT2A + `FLT3-TKD` + ASXL1 + CEBPA + del20 + RUNX1 + Not_available + WT1 + KIT + MLL_AF9 + NUP + 
                         `M4 eos` + `M1/M2` + `M4/M5` + M3 + `None...7` + M2 + M4 + M0  + M6 + M7
                       + sex, data = unique(data_expanded_oslo_fimm_beataml), method = "nearest")




# Check matching results
summary(matchit_out)
plot(summary(matchit_out))


# Extract matched data
matched_data <- match.data(matchit_out)






#Function 
install.packages(c("ggplot2", "dplyr", "MatchIt", "cobalt", "car", "pheatmap", "htmltools"))
library(ggplot2)
library(dplyr)
library(MatchIt)
library(cobalt)
library(car)
library(pheatmap)
library(htmltools)


create_logistic_report <- function(data, formula, grouping_var, output_file = "~/Desktop/logistic_report.html") {
  # Fit the logistic regression model
  mod <- glm(formula, data = data, family = binomial())
  
  # Calculate propensity scores
  data$propensity_score <- predict(mod, type = "response")
  
  # Calculate density range for the plot
  density_values <- density(data$propensity_score)
  density_range <- range(density_values$y)
  y_max <- density_range[2] * 1.1
  
  # Generate VIF values
  vif_values <- vif(mod)
  
  # Create VIF heatmap matrix
  gvif_matrix <- matrix(vif_values, nrow = 1, dimnames = list(NULL, names(vif_values)))
  
  # Create logit values for the plot
  data <- data %>% 
    mutate(logit = log(predict(mod, type = "response") / (1 - predict(mod, type = "response"))))
  
  # Perform matching
  matchit_models <- list(
    matchit_1 = matchit(formula, data = data, method = "nearest",  distance = "glm", replace = FALSE),
    matchit_2 = matchit(formula, data = data, method = "optimal",  distance = "glm", ratio = 2),
    matchit_3 = matchit(formula, data = data, method = "full",  distance = "glm", estimand = "ATE"),
    matchit_4 = matchit(formula, data = data, method = "genetic",  distance = "glm")
  )
  
  # Create plots
  density_plot <- ggplot(data, aes(x = propensity_score, color = as.factor(data[[grouping_var]]))) +
    geom_density(alpha = 0.5) +
    labs(title = "Distribution of Propensity Scores",
         x = "Propensity Score",
         fill = grouping_var) +
    ylim(0, y_max)
  
  histogram_plot <- ggplot(data, aes(x = propensity_score, fill = as.factor(data[[grouping_var]]))) +
    geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.05) +
    labs(title = "Distribution of Propensity Scores",
         x = "Propensity Score",
         fill = grouping_var) +
    ylim(0, 50) +
    theme_minimal() +
    theme(legend.position = "top")
  
  logit_plot <- ggplot(data, aes(x = data[[grouping_var]], y = logit)) + 
    geom_point() + 
    geom_smooth(method = "loess")
  
  # Save plots to temporary files
  density_plot_file <- tempfile(fileext = ".png")
  ggsave(density_plot_file, density_plot, width = 5, height = 3, units = "in")
  
  histogram_plot_file <- tempfile(fileext = ".png")
  ggsave(histogram_plot_file, histogram_plot, width = 3, height = 2.7, units = "in")

  logit_plot_file <- tempfile(fileext = ".png")
  ggsave(logit_plot_file, logit_plot, width = 3, height = 2.7, units = "in")

  
  # VIF heatmap plot
  vif_heatmap <- tempfile(fileext = ".png")
  png(vif_heatmap)
  pheatmap(gvif_matrix, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           color = colorRampPalette(c("white", "red"))(100),
           main = "Heatmap of VIF Values",
           display_numbers = TRUE,
           fontsize_number = 8,
           fontsize_row = 12,
           fontsize_col = 12)
  dev.off()
  
  # Influence plot
  influence_plot <- tempfile(fileext = ".png")
  png(influence_plot)
  influencePlot(mod)
  dev.off()
  
  model_summary <- capture.output(summary(matchit_models[[1]]))
  love_plot_file <- tempfile(fileext = ".png")
  
  png(love_plot_file)
  love.plot(bal.tab(matchit_models[[1]]), stats = c("m", "v"), grid = TRUE, thresholds = c(m=0.25, v=1.25))
  dev.off()
  
  # Create the HTML report
  html_report <- tagList(
    tags$html(
      tags$head(
        tags$title("Logistic Regression Report")
      ),
      tags$body(
        tags$h1("Logistic Regression Summary"),
        tags$pre(paste(capture.output(summary(mod)), collapse = "\n")),
        tags$h1("Propensity Score Distribution"),
        tags$h2("Density Plot"),
        tags$img(src = density_plot_file),
        tags$h2("Histogram Plot"),
        tags$img(src = histogram_plot_file),
        tags$h2("Logit Plot"),
        tags$img(src = logit_plot_file),
        tags$h1("VIF Values"),
        tags$pre(paste(capture.output(vif_values), collapse = "\n")),
        tags$h1("VIF Heatmap"),
        tags$img(src = vif_heatmap),
        tags$h1("Influence Plot"),
        tags$img(src = influence_plot),
        tags$h1("AIC Value"),
        tags$pre(paste(capture.output(AIC(mod)), collapse = "\n")),
        tags$h1("Propensity Score Matching"),
        tags$h2("Matching Models"),
        tags$h2(paste("Matching Model testing")),
        tags$pre(paste(model_summary, collapse = "\n")),
        tags$img(src = love_plot_file)
      )
    )
  )
  
  # Add matching model summaries and love plots
  for (i in seq_along(matchit_models)) {
    model_summary <- capture.output(summary(matchit_models[[i]]))
    love_plot_file <- tempfile(fileext = ".png")
    
    png(love_plot_file)
    love.plot(bal.tab(matchit_models[[i]]), stats = c("m", "ks.statistics"), grid = TRUE, thresholds = c(m=0.25, v=1.25))
    dev.off()
    html_report <- append(html_report, list(
      tags$h2(paste("Matching Model", i)),
      tags$pre(paste(model_summary, collapse = "\n")),
      tags$img(src = love_plot_file)
    ))
  }
  save_html(html_report, file = output_file)
}

create_logistic_report(data_expanded, Binary_lab ~ NRAS + NPM1 + `FLT3-ITD` + EVI1 + KMT2A + `FLT3-TKD` + ASXL1 + CEBPA + del20 + RUNX1 + Not_available + WT1 + KIT + MLL_AF9 + NUP + 
                         `M4 eos` + `M1/M2` + `M4/M5` + M3 + `None...31` + M2 + M4 + M0  + M6 + M7
                       + sex, "Binary_lab", output_file = "~/Desktop/UiO/Project 1/Reports/psm/FIMM_Enserink.html")
dev.off()
dev.off()

create_logistic_report(data_expanded, Binary_lab ~ NRAS + NPM1 + `FLT3-ITD` + EVI1 + KMT2A + `FLT3-TKD` + ASXL1 + CEBPA + del20 + RUNX1 + Not_available + WT1 + KIT + MLL_AF9 + NUP, "Binary_lab", output_file = "~/Desktop/UiO/Project 1/Reports/psm/FIMM_Enserink_mutations_only.html")
dev.off()
dev.off()

library(WeightIt)
W.out <- weightit(Binary_lab ~ inv16 + NRAS + NPM1 + `FLT3-ITD` + t1517 + EVI1 + KMT2A + `FLT3-TKD` + ASXL1 + t911 + CEBPA + del20 + RUNX1 + STAG2 + inv3 + Not_available + t821 + WT1 + AML1_ETO + KIT + MLL_AF9 + NUP + MLL_ELL + 
                              `M4 eos` + `M1/M2` + `M4/M5` + M3 + `None...39` + Other + M2 + M4 + M0  + M6 + M7 + General
                            + sex,
                            data = data_expanded,
                            method = "glm",
                            estimand = "ATT")

bal.plot(W.out, var.name = "FLT3-ITD", which = "both")
bal.plot(W.out, var.name = "prop.score",
         which = "both",
         type = "histogram",
         mirror = TRUE)

bal.plot(W.out, var.name = "FLT3-ITD")

love.plot(W.out, stats = c("c", "ks"),
          thresholds = c(cor = .1), 
          abs = TRUE, wrap = 20,
          limits = list(ks = c(0, .5)),
          var.order = "unadjusted", line = TRUE)





#----BeatAML and FIMM ----

#BeatAML clinical info 
beatAML_clinical_info <- read_csv('~/Desktop/UiO/Project 1/Data/Second cleansing/BeatAML_clinical_info.csv')
beatAML_clinical_info <- beatAML_clinical_info %>% mutate(X = paste0(beatAML_clinical_info$dbgap_subject_id, '_', beatAML_clinical_info$dbgap_sample_id, '_', beatAML_clinical_info$dbgap_rnaseq_sample))
beatAML_clinical_info$X
beatAML_clinical_info$lab <- 'BeatAML'

#FIMM clinical infor
fimm_clinical_info <- read_csv('~/Desktop/UiO/Project 1/Data/Second cleansing/FIMM_clinical_info.csv')
fimm_clinical_info$lab <- 'FIMM'

#bind the two clinical datasets
combined_clinical_info_FIMM_beatAML <- rbind(beatAML_clinical_info[,c('X', 'mutations', 'FAB_class', 'sex', 'age', 'lab')], fimm_clinical_info[,c('X', 'mutations', 'FAB_class', 'sex', 'age', 'lab')])

#Read DSS info for each dataset
BeatAML_response_scores <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_beat_aml_full_drug_set.csv')
BeatAML_response_scores$lab <- 'BeatAML'
BeatAML_response_scores$X <- BeatAML_response_scores$Patient.num

fimm_response_score <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set.csv')
fimm_response_score$lab <- 'FIMM'

fimm_response_score_1 <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_drug_response.csv')
fimm_response_score_1 <- gather(fimm_response_score_1, Patient.num, DSS2, 'AML_084_04':'Healthy_17', factor_key=TRUE)
fimm_response_score_1$...1 <- NULL
fimm_response_score_1$Drug_ID <- NULL
colnames(fimm_response_score_1)[colnames(fimm_response_score_1) == "Drug_name"] <- "drug"
fimm_response_score_1$lab <- 'FIMM'

remove_after_second_underscore <- function(x) {
  parts <- unlist(strsplit(x, "_"))
  if (length(parts) > 2) {
    return(paste(parts[1:2], collapse = "_"))
  } else {
    return(x)
  }
}

# Apply the function to the data frame
fimm_response_score_1$X<- sapply(as.character(fimm_response_score_1$Patient.num), remove_after_second_underscore)

#Bind the response scores for the two datasets
combined_response_score_beatAML_FIMM <- rbind(BeatAML_response_scores[,c('Patient.num', 'drug', 'DSS2', 'lab', 'X')], fimm_response_score_1)

#Merge clinical data and response scores 
all_combined_FIMM_BeatAML <- merge(combined_response_score_beatAML_FIMM, combined_clinical_info_FIMM_beatAML, by="X", x.all=TRUE)


#Regress group with variables and DSS2
all_combined_FIMM_BeatAML$Binary_lab <- ifelse(all_combined_FIMM_BeatAML$lab == "FIMM", 1, 0)
all_combined_FIMM_BeatAML$mutations <- ifelse(is.na(all_combined_FIMM_BeatAML$mutations), "None", all_combined_FIMM_BeatAML$mutations)
all_combined_FIMM_BeatAML$FAB_class <- ifelse(is.na(all_combined_FIMM_BeatAML$FAB_class), "None", all_combined_FIMM_BeatAML$FAB_class)
all_combined_FIMM_BeatAML$karyotype_class <- ifelse(is.na(all_combined_FIMM_BeatAML$age), "None", all_combined_FIMM_BeatAML$age)

length(unique(all_combined_FIMM_BeatAML$X))
df_FIMM_BeatAML <- unique(all_combined_FIMM_BeatAML[,c('X', 'lab', 'mutations', 'FAB_class', 'Binary_lab', 'sex', 'age')])

# Apply the function to the data frame
data_expanded_FIMM_BeatAML <- expand_columns(df_FIMM_BeatAML, c("mutations","FAB_class"))

#all_metrices_all_labs = read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/all_metrices_all_labs.csv')

#data_expanded_FIMM_BeatAML_all_metrics <- merge(data_expanded_FIMM_BeatAML, all_metrices_all_labs, by.x="X", by.y="Patient.num")
colnames(data_expanded_FIMM_BeatAML)

unique(subset(all_combined_FIMM_BeatAML, lab == 'BeatAML', select = 'mutations'))


library(ggfortify)
data_for_pca <- data_expanded_FIMM_BeatAML
data_for_pca$Binary_sex <- ifelse(data_for_pca$sex == "Male", 1, 0)
#data_for_pca <- unique(subset(data_for_pca, select = -c(DSS2, drug)))
data_for_pca <- na.omit(data_for_pca)

# Identify columns with only one unique value
single_value_columns <- sapply(data_for_pca, function(x) length(unique(x)) == 1)

# Print names of columns with only one unique value
columns_removed <- names(single_value_columns)[single_value_columns]
print("Columns with only one unique value:")
print(columns_removed)

# Remove columns with only one unique value
data_for_pca <- data_for_pca[, !single_value_columns]

# Identify columns with only one unique value
single_value_columns_lab <- sapply(subset(data_for_pca, lab == 'FIMM'), function(x) length(unique(x)) == 1)

# Print names of columns with only one unique value
columns_removed <- names(single_value_columns_lab)[single_value_columns_lab]
print("Columns with only one unique value:")
print(columns_removed)

# Remove columns with only one unique value
data_for_pca_filtered <- data_for_pca[, !single_value_columns_lab]


# Print the cleaned data frame
print("Cleaned data frame:")
print(data_for_pca_filtered)

# Prepare data
for_pca <- subset(data_for_pca_filtered, select = -c(X, mutations, FAB_class, sex))  # Exclude the species column for PCA

# Perform PCA
pca_result <- prcomp(for_pca, center = TRUE, scale. = TRUE)

# Convert PCA result to a data frame for plotting
pca_df <- as.data.frame(pca_result$x)

# Add species labels for coloring
pca_df$lab <- data_for_pca$lab
pca_df$id <- data_for_pca$X  # Assuming `data_for_pca$id` contains unique identifiers for each point

# Plot PCA with labeled points
ggplot(pca_df, aes(x = PC1, y = PC2, color = lab)) +
  geom_point(size = 3) +
  geom_text(aes(label = id), vjust = 1, hjust = 1, size = 3) +  # Add labels to each point
  labs(title = "PCA Plot of Iris Dataset", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = lab)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot of Iris Dataset", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()



colnames(data_for_pca_filtered)
NPM1 + WT1 + KIT + None...2228 + AML1_ETO + `FLT3-ITD`+ MLL_AF9 + NUP + MLL_ELL + `M4 eos` + `M4/M5` + `None...2236` + Other + M2 + M4 + M1 + M5 + General
create_logistic_report(na.omit(unique(data_expanded_FIMM_BeatAML)), Binary_lab ~ age + sex +`None...2236` + M4 + M1 + M5
                      , "Binary_lab", output_file = '~/Desktop/UiO/Project 1/Reports/psm/FIMM_BeatAML_fab_classes.html')
create_logistic_report(na.omit(unique(data_expanded_FIMM_BeatAML)), Binary_lab ~ age +NPM1 + WT1 + KIT + AML1_ETO + `FLT3-ITD`+ MLL_AF9 + NUP + MLL_ELL
                       , "Binary_lab", output_file = '~/Desktop/UiO/Project 1/Reports/psm/FIMM_BeatAML_mutations.html')
dev.off()
dev.off()


cor_matrix <- cor(data_expanded_FIMM_BeatAML[, !names(data_expanded_FIMM_BeatAML) %in% c('DSS2','drug','Patrient.num','lab', 'X', 'FAB_class', 'mutations', 'sex', "Binary_lab")])
print(cor_matrix)

corrplot(cor_matrix, 
         col = colorRampPalette(c("blue", "white", "red"))(200), # Color gradient
         type = "upper", # Show only upper triangle
         addCoef.col = NULL, # Do not add correlation coefficients
         tl.col = "black", # Text label color
         tl.srt = 45) 


mod_test <- glm(as.factor(lab) ~ sex + age, data=data_expanded_FIMM_BeatAML_all_metrics, family="binomial")
mod_test2 <- lm(DSS2 ~ mutations + FAB_class + drug, data=all_combined_FIMM_BeatAML)
mod_test_binary <- glm(Binary_lab ~ None...2236+M1+M5,M4+ M2+ `M0/M1`+age, 
                       data=data_expanded_FIMM_BeatAML, family="binomial")

summary(mod_test1)
summary(mod_test2)
summary(mod_test_binary) 

subset(all_combined_FIMM_BeatAML, lab=='BeatAML', select=c(sex,age,Patient.num,lab))
all_combined_FIMM_BeatAML_subset <- subset(all_combined_FIMM_BeatAML, !is.na(sex) & !is.na(age))
attach(all_combined_FIMM_BeatAML)
LAB = cbind(Binary_lab)
DSS2 = cbind(DSS2)
X = cbind(sex, age)

glm1 <- glm(LAB ~ X, family=binomial, data=all_combined_FIMM_BeatAML)

library(Matching)
rr1 <- Match(Y=DSS2, Tr=LAB, X=glm1$fitted)

MatchBalance(LAB ~ X, match.out = rr1)

alias(mod_test1)
alias(mod_test2)
alias(mod_test_binary)


with(data_expanded_FIMM_BeatAML, t.test(sex ~ Binary_lab))

vif(mod_test1)
vif(mod_test2)
vif(mod_test_binary)


df <- na.omit(unique(all_combined_FIMM_BeatAML[,c('X', 'sex', 'age', 'Binary_lab', 'lab', 'mutations')]))
df$propensity_score <- predict(mod_test, type = "response")
density_values <- density(df$propensity_score)
density_range <- range(density_values$y)

y_max <- density_range[2] * 1.1

ggplot(df, aes(x = propensity_score, color = as.factor(lab))) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Propensity Scores",
       x = "Propensity Score",
       fill = "Binary_lab") +
  ylim(0, y_max)


# Create histogram plot with overlapping bars
ggplot(df, aes(x = propensity_score, fill = as.factor(lab))) +
  geom_histogram(bins = 30, color="black") +
  labs(title = "Distribution of Propensity Scores",
       x = "Propensity Score",
       fill = "Binary_lab") +
  ylim(0, 50) +
  theme_minimal() +
  theme(legend.position = "top")











#Distribution FAB classes between labs

df_long <- combined_clinical_info_FIMM_beatAML %>%
  separate_rows(FAB_class, sep = ",") %>%
  group_by(lab, FAB_class) %>%
  summarise(count = n(), .groups = 'drop')

n_obs <- combined_clinical_info_FIMM_beatAML %>%
  group_by(lab) %>%
  summarise(n_obs = n(), .groups = 'drop')

# Join the number of observations with the counts
df_long <- df_long %>%
  left_join(n_obs, by = "lab") %>%
  mutate(mean_count = count / n_obs)

ggplot(df_long, aes(x = FAB_class, y = mean_count, fill = as.factor(lab))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "FAB classes",
       x = "FAB class",
       y = "count / number of observations",
       fill = "Lab") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Distribution mutations classes between labs

df_long <- combined_clinical_info_FIMM_beatAML %>%
  separate_rows(mutations, sep = ",") %>%
  group_by(lab, mutations) %>%
  summarise(count = n(), .groups = 'drop')

n_obs <- combined_clinical_info_FIMM_beatAML %>%
  group_by(lab) %>%
  summarise(n_obs = n(), .groups = 'drop')

# Join the number of observations with the counts
df_long <- df_long %>%
  left_join(n_obs, by = "lab") %>%
  mutate(mean_count = count / n_obs)

ggplot(df_long, aes(x = mutations, y = mean_count, fill = as.factor(lab))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mutations",
       x = "Mutation",
       y = "count / number of observations",
       fill = "Lab") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
