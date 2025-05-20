
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


enserink_fresh_vs_frozen<- enserink_fresh_and_frozen_patients %>%
  group_by(Patient.num, drug) %>%
  summarise(
    Fresh = DSS2[sample == 'fresh'],
    Frozen = DSS2[sample == 'frozen'],
    .groups = "drop"
  )


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

# Calculate adjusted R^2
adjusted_r2 <- summary(lm_model)$adj.r.squared

# Extract the p-value for Fresh
pval <- summary(lm_model)$coefficients[2, 4] 

# Apply FDR adjustment (assuming multiple tests, example 10)
pval_fdr <- p.adjust(pval, method = "fdr")

# Prepare the annotation text
annotation_text <- paste0(
  "R² = ", round(adjusted_r2, 3), 
  "\nFDR p-value = ", format(pval_fdr, scientific = TRUE, digits = 3)
)

custom_colors <- c(
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", 
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", 
  "#bc80bd", "#ccebc5"  # Added two complementary colors
)
display.brewer.pal(n = 10, name = "Set3")
# Create the plot
scatterplot_fresh_vs_frozen_enserink <- ggplot(enserink_fresh_vs_frozen, aes(x = Fresh, y = Frozen)) +
  geom_point(aes(color = factor(Patient.num)), size = 3) +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  scale_color_manual(values = custom_colors) +
  annotate(
    "text", 
    x = mean(enserink_fresh_vs_frozen$Fresh), y = max(enserink_fresh_vs_frozen$Frozen), 
    label = annotation_text, 
    size = 3.51, family = "Arial", hjust = 0, vjust = 1
  ) +
  labs(color = "Patient") +
  theme_minimal() + 
  theme(
    text = element_text(family = "Arial", size = 10),    # All text Arial, size 10
    axis.text = element_text(size = 8, family = "Arial"), # Tick labels Arial, size 8
    legend.text = element_text(size = 10, family = "Arial"), # Legend text Arial, size 10
    legend.title = element_text(size = 10, family = "Arial"), # Legend title Arial, size 10
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(), # Optional: Remove background panel
    plot.background = element_blank()   # Optional: Remove plot background
  )

ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/scatterplot_fresh_vs_frozen_enserink.png", plot = scatterplot_fresh_vs_frozen_enserink, width = 10.76, height = 8.09, units="cm") #width = 20.76, height = 7.09 # width = 10.76, height = 8.09
