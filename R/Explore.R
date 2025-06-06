# -------------------------------------------------------- Script overview
# Script Name: Explore.R
# Author: Katarina Willoch
# Date: 2025-06-01
#
# Description:
# This script reads in dose response metric scores DSS1, DSS2, DSS3, IC50, AUC, rAUC,
# and generates heatmaps and PPCA visualizations
# 
#
# Sections:
#   1. Load libraries and data
#   2. Heatmap generation for DSS2
#   3. PPCA generation for DSS2
#   4. Heatmap and PPCA generation for all metrics
#   5. ggbetweenstats for DSS2
# -------------------------------------------------------- Script overview

#Install and import libraries ----
library(dplyr)
setwd("/Users/katarinawilloch/")
#source('~/Desktop/UiO/Project 1/code/R/plot_functions_project1.R')

# Install and load necessary packages
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}
if (!requireNamespace("grid", quietly = TRUE)) {
  install.packages("grid")
}
if (!requireNamespace("ggplotify", quietly = TRUE)) {
  install.packages("ggplotify")
}
library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(gridExtra)
library(grid)
library(ggplotify)
library(ComplexHeatmap)
library(pcaMethods)
library(cluster)
library(proxy)
library(ggrepel)
library(ggnewscale)
library(circlize)
library(ggstatsplot)
library(stringr)

#Figure output location
figure_output <- 'Desktop/UiO/Project 1/Figures/New karolinska data/'

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

#BoxCox transforming DSS2
all_response_metrics$DSS2_pos <- all_response_metrics$DSS2 + 0.01
MASS::boxcox(DSS2_pos ~ lab, 
             data = all_response_metrics,
             lambda = seq(-0.25, 2, length.out = 10))

all_response_metrics$DSS2_boxcox <- ((all_response_metrics$DSS2)^0.5 - 1) / 0.5
#Z-scaling (across all patients and all labs) DSS2
for(j in unique(all_response_metrics$drug)){
  all_response_metrics$DSS2_boxcox_sclaed2[all_response_metrics$drug==j] <-
    scale(all_response_metrics$DSS2_boxcox[all_response_metrics$drug==j])
}

#Creating copy of all_response_metrics so it can be displayed in the box-violin plots
copy_all_response_metrics <- all_response_metrics
subset(all_response_metrics, drug == "Everolimus" & lab == "Beat AML")
#remove drug 001, RAD because there are only 3 patients tested for it in the Beat AML cohort
all_response_metrics <- subset(all_response_metrics, drug != "Everolimus") %>% as.data.frame()
#dataframe to wide format for the heatmap and ppca
all_response_metrics_wide <- all_response_metrics[,c('Patient.num', 'drug', 'DSS2')] %>% pivot_wider(names_from = Patient.num, values_from = DSS2) %>% as.data.frame()
rownames(all_response_metrics_wide) <- all_response_metrics_wide$drug
all_response_metrics_wide$drug <- NULL


#----Heatmap creation----
#Lab annotation for heatmap creation
a_col <- unique(subset(all_response_metrics, select = c(lab, Patient.num))) %>% as.data.frame()
rownames(a_col) <- a_col$Patient.num
a_col$Patient.num <- NULL
a_col <- dplyr::rename(a_col, "Study" = "lab")
a_col$Study <- factor(a_col$Study, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
#Using ComplexHeatmap::HeatmapAnnotation to create top_annotation for heatmap
a_col_h <- ComplexHeatmap::HeatmapAnnotation(df = a_col, annotation_legend_param = list(
  title = "Study",
  title_gp = gpar(fontsize = 11, fontface = "bold"),
  labels_gp = gpar(fontsize = 11),
  legend_direction = "vertical",
  border = FALSE
),
col = list(Study= c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")), 
show_annotation_name = FALSE
)

# Custom function for calculating distance that ignores NAs
dist_no_na <- function(x) {
  d <- as.dist(proxy::dist(x, method = "euclidean", pairwise = TRUE))  #pairwise handles NAs
  d[is.na(d)] <- max(d, na.rm = TRUE)  # Replace NA distances with max value
  as.dist(d)
}

set.seed(2025)
#Creating complex heatmap for DSS2
all_labs_dss2_heatmap <- ComplexHeatmap::Heatmap(
  as.matrix(all_response_metrics_wide),
  name = "DSS2",  #Legend header
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_names_gp = gpar(fontsize = 9),  #Row font size 10
  column_names_gp = gpar(fontsize = 1),  #Column font size 10
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_column_names = FALSE,
  top_annotation = a_col_h,  #Column annotations
  left_annotation = NULL,  #No row annotations
  na_col = "grey",
  row_gap = unit(1, "cm"),             #Space between rows
  column_gap = unit(1, "cm"),          #Space between columns
  heatmap_legend_param = list(
    title = expression("DSS"[2]),  #Add header to the legend
    title_gp = gpar(fontsize = 11, fontface = "bold"),  #Legend title font
    labels_gp = gpar(fontsize = 10), 
    grid_height = unit(8, "mm")
  ),
  width = unit(8, "cm"),  # 12 #20
  height = unit(12, "cm"),  # 9.6 #16
  column_title = "",
  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  row_dend_width = unit(2, "cm"),
  clustering_distance_rows = dist_no_na,
  column_dend_height = unit(2, "cm"),
  clustering_distance_columns = dist_no_na
)
quartz()
plot(
  all_labs_dss2_heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE, 
  padding = unit(c(0, 0, 0, 0), "mm")
)

#Save the heatmap plot as a pdf; this was not the dimensions used for the figure in the article (this was done by copying the plot from rstudio)
pdf(file = paste0(figure_output,"all_labs_dss2_heatmap.pdf"), width = 8, height = 12)  # PDF output
quartz()
plot(
  all_labs_dss2_heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE, 
  padding = unit(c(0, 0, 0, 0), "mm")
)
dev.off()  #Close the graphic device


#----PPCA plot creation----
#Perform ppca
ppca <- pca(t(all_response_metrics_wide), method="ppca", nPcs=3, seed=123)

#Extract percent variance for PC1 and PC2
percent_var <- round(ppca@R2 * 100, 1)

#Create a data frame for plotting explained variance
var_df <- data.frame(
  PC = paste0("PC", seq_along(percent_var)),
  Variance = percent_var
)
#Plotting explained variance for each PC PC1 + PC2 explains more than 50% of the variance
quartz()
ggplot(var_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(round(Variance, 1), "%")), vjust = -0.5, size = 4) +
  labs(title = "Explained Variance by Principal Component",
       y = "Percent Variance Explained",
       x = "Principal Component") +
  theme_minimal()

#Get the estimated complete observations
ppca_scores <- scores(ppca)
#Get unique information for each sample by averaging DSS2
avg_dss2_all_response_metrics <- all_response_metrics[,c('Patient.num', 'lab', 'DSS2')] %>% group_by(Patient.num, lab) %>% dplyr::summarize('Average DSS2' = mean(DSS2)) %>% as.data.frame()
rownames(avg_dss2_all_response_metrics) <- avg_dss2_all_response_metrics$Patient.num
#Merging ppca results with the unique information for each sample
ppca_data <- merge(ppca_scores, avg_dss2_all_response_metrics, by = 'row.names')

ppca_data$lab <- factor(ppca_data$lab, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
#Creating PPCA plot with the first two PCs colored by lab
ppca_plot <- ggplot(as.data.frame(ppca_data), aes(x = PC1, y = PC2, color=lab)) +
  geom_point() +
  labs(
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)"),
    color = ""
  ) +
  scale_color_manual(values = c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")) +
  labs(color = "") + 
  guides(
    color = guide_legend(
      override.aes = aes(shape = 15, size = 4, width = 1.5, height = 1),  # Adjust legend symbol size
      label.spacing = unit(0.5, "cm")  # Reduce space between text and legend symbol
    )
  ) +  
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

#Printing and saving ppca plot
print(ppca_plot)
ggsave(paste0(figure_output,"PPCA_plot_all_data_all_labs.png"), plot = ppca_plot, width = 10, height = 8, dpi = 300)

##----Average DSS2 PPCA plot ----
ppca_data$`Average DSS2` <- as.numeric(ppca_data$`Average DSS2`)
#Creating ppca plot of the first two PCs colored by average DSS2 for each patient in each lab
dss2_ppca_plot <- ggplot(as.data.frame(ppca_data), aes(x = PC1, y = PC2, color = `Average DSS2`)) +
  #Beat AML coloring
  geom_point(aes(color = ifelse(lab == "Beat AML", `Average DSS2`, NA)), size = 4) +  
  scale_color_gradient2(low = "#ccece6", mid = "#8dd3c7", high = "#145A32", midpoint = median(ppca_data[ppca_data$lab == "Beat AML",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Beat AML", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +  
  new_scale_color() +  # Reset color scale for other points
  
  #Oslo coloring
  geom_point(aes(color = ifelse(lab == "Oslo", `Average DSS2`, NA)), size = 4) +  
  scale_color_gradient2(low = "#ffe6c7", mid = "#fdb462", high = "#a34700",midpoint = median(ppca_data[ppca_data$lab == "Oslo",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Oslo", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +
  new_scale_color() +
  
  #Helsinki coloring
  geom_point(aes(color = ifelse(lab == "Helsinki", `Average DSS2`, NA)), size = 4) +  
  scale_color_gradient2(low = "#ffd2cc", mid = "#fb8072", high = "#a12b1d", midpoint = median(ppca_data[ppca_data$lab == "Helsinki",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Helsinki", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +
  new_scale_color() +
  
  #Karolinska coloring
  geom_point(aes(color = ifelse(lab == "Karolinska", `Average DSS2`, NA)), size = 4) +  
  scale_color_gradient2(low = "#d6ecff", mid = "#80b1d3", high = "#23537d", midpoint = median(ppca_data[ppca_data$lab == "Karolinska",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Karolinska", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +
  labs(
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)"),
    color = ""
  ) +
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
#Printing and saving plot
print(dss2_ppca_plot)
ggsave(paste0(figure_output,"avg_dss_score_colorsPPCA_plot_all_data_all_labs_v1.png"), plot = dss2_ppca_plot, width = 10, height = 8, dpi = 300)

#Help explain the results
centroids <- ppca_data %>%
  group_by(lab) %>%
  summarise(centroid_PC1 = mean(PC1), centroid_PC2 = mean(PC2))

pca_data <- ppca_data %>%
  left_join(centroids, by = "lab") %>%
  mutate(distance_to_centroid = sqrt((PC1 - centroid_PC1)^2 + (PC2 - centroid_PC2)^2))

furthest_points <- pca_data %>%
  group_by(lab) %>%
  top_n(3, distance_to_centroid) 
#Show the 3 most distant points from each lab, shows that average DSS2 is higher for the further away points
print(furthest_points)

#----Other metrics - PPCA and Heatmap plots----
#All metrics used to analyse
metrics <- c("DSS1", "DSS2", "DSS3", "AUC", "auc_a", "IC50", "log10IC50", "DSS2_boxcox_sclaed2") #"log10IC50"
quartz()
#Creating Heatmap and PPCA plot for each metric
for (m in metrics){
  print(m)
  #Get wide format dataframe for metric m 
  specific_all_response_metrics_wide <- pivot_wider(all_response_metrics[,c("drug", "Patient.num", "lab", m)], names_from = drug, values_from = m) %>% as.data.frame()
  rownames(specific_all_response_metrics_wide) <- specific_all_response_metrics_wide$Patient.num
  specific_all_response_metrics_wide$Patient.num <- NULL
  
  #calculate ppca reaults for metric m
  ppca_metrics <- pca(specific_all_response_metrics_wide[,-c(1)], method="ppca", nPcs=3, seed=123)
  ppca_scores_metrics <- scores(ppca_metrics)
  ppca_data_metrics <- merge(ppca_scores_metrics, specific_all_response_metrics_wide, by = 'row.names')
  
  #Extract percent variance for PC1 and PC2
  percent_var <- round(ppca_metrics@R2 * 100, 1)
  
  #Create a data frame for plotting explained variance
  var_df <- data.frame(
    PC = paste0("PC", seq_along(percent_var)),
    Variance = percent_var
  )
  #Plotting explained variance for each PC PC1 + PC2 explains more than 50% of the variance
  print(ggplot(var_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = paste0(round(Variance, 1), "%")), vjust = -0.5, size = 4) +
    labs(title = "Explained Variance by Principal Component",
         y = "Percent Variance Explained",
         x = "Principal Component") +
    theme_minimal())
  
  #Creating PPCA plot of the first two PCs for metric m colored by study
  ppca_plot_metric <- ggplot(as.data.frame(ppca_data_metrics), aes(x = PC1, y = PC2, color=lab)) +
    geom_point() +
    scale_color_manual(values = c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072", "Karolinska"="#80b1d3")) +
    labs(
      x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
      y = paste0("PC2 (", round(percent_var[2], 1), "%)"),
      color = ""
    ) +
    guides(
      color = guide_legend(
        override.aes = aes(shape = 15, size = 4, width = 1.5, height = 1),  #Adjust legend symbol size
        label.spacing = unit(0.5, "cm")  #Reduce space between text and legend symbol
      )
    ) +  
    #theme_minimal()+
    theme(
      panel.background = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 12, family = "Arial"),    #Center legend text
      legend.position = "top",                                
      legend.justification = "center",        #Ensure the legend is centered
      legend.box = "horizontal",
      legend.key = element_blank(),                            #Remove key background
      legend.spacing.x = unit(0.2, "cm"),  
      legend.key.size = unit(0.5, "lines"),      #Adjust size of the colored squares
      legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
      plot.margin = margin(10, 10, 10, 10),
      axis.title = element_text(size = 12, family = "Arial"),    #Axis title size
      axis.text = element_text(size = 10, family = "Arial")
    )
  #Printing the resulting ppca plot
  print(ppca_plot_metric)
  
  #Creating function for claulating distance with missing values
  dist_no_na <- function(x) {
    d <- as.dist(proxy::dist(x, method = "euclidean", pairwise = TRUE))  # pairwise handles NAs
    d[is.na(d)] <- max(d, na.rm = TRUE)  # Replace NA distances with max value
    as.dist(d)
  }
  
  #Creating top annotation for the heatmap
  col_anno <- subset(specific_all_response_metrics_wide, select="lab") %>% dplyr::rename(Study = lab) %>% as.data.frame()
  col_anno <- factor(col_anno$Study, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
  levels(col_anno)
  a_col_h <- ComplexHeatmap::HeatmapAnnotation(Study = col_anno, 
                                               col = list(Study= c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")), show_annotation_name = FALSE, gp = gpar(fonface = "plain"))
  
  #Printing the current metric and changing the title of the legend slightly
  print(m)
  m_header <- if(m == "DSS2"){expression("DSS"[2])}
    else if(m == "DSS1"){expression("DSS"[1])}
    else if (m == "DSS3"){expression("DSS"[3])}
    else if(m == "IC50"){expression("IC"[50])}
    else if(m == "auc_a"){"rAUC"}
    else if(m == "DSS2_boxcox_sclaed2"){expression("BoxCox Z-scaled DSS"[2])}
    else{m}
  print(m_header)
  
  set.seed(2025)
  #Creating heatmap for metric m
  h <- ComplexHeatmap::Heatmap(
    as.matrix(t(specific_all_response_metrics_wide[,-c(1)])),
    name = m,  
    col = if(m == "IC50"){colorRamp2(c(10000, 1000, 0), c("blue", "white", "red"))},
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    row_names_gp = gpar(fontsize = 9),  #Row font size 10
    column_names_gp = gpar(fontsize = 1),  #Column font size 10
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    show_column_names = FALSE,
    top_annotation = a_col_h,  #Column annotations
    left_annotation = NULL,  #Row annotations
    na_col = "grey",
    row_gap = unit(1, "cm"),             #Space between rows
    column_gap = unit(1, "cm"),          #Space between columns
    heatmap_legend_param = list(
      title = m_header,  # Add header to the legend
      title_gp = gpar(fontsize = 11, fontface = "bold"),  #Legend title font
      labels_gp = gpar(fontsize = 10), 
      grid_height = unit(8, "mm")
    ),
    #annotation_legend_param = list(title = ""), 
    width = unit(8, "cm"),  # 12 #20
    height = unit(12, "cm"),  # 9.6 #16
    column_title = " ",
    column_title_gp = gpar(fontsize = 20), 
    row_dend_width = unit(2, "cm"),
    clustering_distance_rows = dist_no_na,
    column_dend_height = unit(2, "cm"),
    clustering_distance_columns = dist_no_na
  )
  
  #Showing the resulting heatmap
  print(plot(
    h,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE, 
    padding = unit(c(0, 0, 0, 0), "mm")
  ))
  #Loading relevant library
  library(patchwork)
  heat_g <- grid.grabExpr(plot(
    h,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE, 
    padding = unit(c(0, 0, 0, 0), "mm")
  ))
  #Used to adjust positioning of heatmap
  blank_space <- textGrob(" ", gp = gpar(fontsize = 24))
  #Showing combined plot
  grid.newpage()
  print(grid.arrange(ppca_plot_metric, arrangeGrob(heat_g, bottom = blank_space), ncol = 2, widths = c(1, 1)))
  grid.newpage()
  #Saving combined plot
  ggsave(paste0(figure_output, m, '_heatmap_and_ppca.png'), plot = grid.arrange(ppca_plot_metric, arrangeGrob(heat_g, bottom = blank_space), ncol = 2, widths = c(1, 1)), width = 6.47 * 2, height = 5.94 * 1, dpi = 300)
}


#----ggbetweenstats----
custom_colors <- c("Beat AML" = "#8dd3c7",  "Oslo" = "#fdb462", "Helsinki" = "#fb8072", "Karolinska"= "#80b1d3")
selected_groups <- unique(all_response_metrics$drug)[2:15] # Select first 6 groups
selected_groups <- c("5-Azacytidine", "Bortezomib", "BI 2536", "Dasatinib")
filtered_data <- all_response_metrics[all_response_metrics$drug %in% selected_groups, ]

subset(copy_all_response_metrics, drug == "001, RAD")
#Creating box-violin plots for each drug 
##Each plot shows DSS2 on y axis and lab on x axis 
##Statistical tests are performed to show the significance of the difference in DSS2 between labs
g_all_drugs <- grouped_ggbetweenstats( # paired samples
  data = copy_all_response_metrics,
  x = lab,
  y = DSS2,
  grouping.var = drug,
  type = "nonparametric", # for wilcoxon
  centrality.plotting = FALSE, # remove median
  ggstatsplot.layer = TRUE,
  xlab = "",        
  ylab = expression("DSS"[2]),
  results.subtitle = TRUE,  # Removes default subtitle
  #subtitle = "{test} (p = {p})", 
  p.adjust.method = "bonferroni",
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

#Updating labels and subtitles of plots
#47 drugs including 001, RAD
for (i in 1:47) {
  #Check if the current g[[i]] contains annotations
  if (!is.null(g_all_drugs[[i]]$layers[[4]]$stat_params$annotations)) {
    
    #Extract the annotations for this g[[i]]
    annotations <- g_all_drugs[[i]]$layers[[4]]$stat_params$annotations

    #Modify the annotations by removing unwanted parts
    annotations_modified <- unlist(lapply(annotations, function(annot) {
      #Remove LaTeX parts (italic(p) and adj.)
      annot <- gsub("list\\(italic\\(p\\)\\[\\\"Bonferroni\\\"\\s*-\\s*adj\\.\\] == \"", "list(Bonferroni == ", annot)
      annot <- gsub("\"\\)", ")", annot)  # Remove the trailing ')'
      print(annot)
      
      return(annot)
    }))
    
    # Update the annotations in g[[i]]
    g_all_drugs[[i]]$layers[[4]]$stat_params$annotations <- annotations_modified
  }
  if(!is.null(g_all_drugs[[i]]$labels$subtitle)){
    subtitle <- g_all_drugs[[i]]$labels$subtitle
    subtitle_str <- deparse(subtitle)  
    subtitle_str <- paste(subtitle_str, collapse = " ")  # Ensure it's a single string
    
    # Extract only "Kruskal-Wallis: italic(p) == value"
    subtitle_mod <- str_replace(subtitle_str, 'list\\(chi\\["Kruskal-Wallis"\\]\\^2.*?italic\\(p\\) ==      "', "Kruskal-Wallis: p = ")
    subtitle_mod <- str_replace(subtitle_mod, '", widehat\\(epsilon\\)\\["ordinal"\\]\\^2.*', "")
    p_val <- as.numeric(sub("Kruskal-Wallis: p = ", "", subtitle_mod))
    
    # Modify subtitle_mod based on value
    if (!is.na(p_val) && p_val < 0.001) {
      subtitle_mod <- "p < 0.001"
    }
    g_all_drugs[[i]]$labels$subtitle <- subtitle_mod
  }
}

# gb <- ggplot_build(g)
# gb$data
# # Now print or display the modified plot
# pdf("Desktop/UiO/Project 1/Figures/draw/Difference_DSS2_per_drug_v1.pdf", width = 80, height = 30)
# ggplot_build(g)
# dev.off()
#Print and save box-violin plots
print(g_all_drugs)
ggsave(paste0(figure_output,"All_drugs_Difference_DSS2_per_drug.png"), plot = g_all_drugs, width = 80, height = 30, units = "cm", limitsize = FALSE)

#Create box-violin plots for a few drugs for article figure
g_filtered_data <- grouped_ggbetweenstats( # paired samples
  data = filtered_data,
  x = lab,
  y = DSS2,
  grouping.var = drug,
  type = "nonparametric", # for wilcoxon
  centrality.plotting = FALSE, # remove median
  xlab = "",        
  ylab = expression("DSS"[2]),
  results.subtitle = TRUE,  # Removes default subtitle
  p.adjust.method = "bonferroni",
  pairwise.display = "none", #"significant"
  #pairwise.comparisons = TRUE,
  plotgrid.args = list(nrow = 2, ncol = 2),
  facet_wrap.args = list(scales = "fixed", strip.position = "top"),
  plot.margin = margin(0, 10, 10, 10),
  violin.args = list(alpha = 0),
  boxplot.args = list(alpha = 0),
  ggplot.component = list(
    coord_cartesian(ylim = c(0, 50)),
    scale_color_manual(values = custom_colors), 
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1))), 
    theme(
      panel.background = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 0)),
      axis.text.y.right = ggplot2::element_blank(), 
      axis.ticks.y.right = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = element_text(family = "Arial", face = "plain",size = 12, color = "black", angle = 25, hjust = 0.6), #0.6
      axis.title.y = element_text(family = "Arial", face = "plain", size = 12, color = "black"),
      plot.title = element_text(family = "Arial", face = "plain", size = 12, color = "black", hjust = 0.5, margin = margin(b = 0))))
) 

#Updating labels and subtitles of plots
for (i in 1:4) {
  if(!is.null(g_filtered_data[[i]]$labels$subtitle)){
    subtitle <- g_filtered_data[[i]]$labels$subtitle
    
    subtitle_str <- deparse(subtitle)  
    subtitle_str <- paste(subtitle_str, collapse = " ")  # Ensure it's a single string
    
    # Extract only "Kruskal-Wallis: italic(p) == value"
    subtitle_mod <- str_replace(subtitle_str, 'list\\(chi\\["Kruskal-Wallis"\\]\\^2.*?italic\\(p\\) ==      "', "p = ")
    subtitle_mod <- str_replace(subtitle_mod, '", widehat\\(epsilon\\)\\["ordinal"\\]\\^2.*', "")
    p_val <- as.numeric(sub("p = ", "", subtitle_mod))
    
    # Modify subtitle_mod based on value
    if (!is.na(p_val) && p_val < 0.001) {
      subtitle_mod <- "p < 0.001"
    } 
    g_filtered_data[[i]]$labels$subtitle <- subtitle_mod
  }
}
#Print and save box-violin plots for the 4 drugs; the plot used in the article is copied from R studio not the saved plot
print(g_filtered_data)
ggsave(paste0(figure_output, "4_drugs_Difference_DSS2_per_drug.png"), plot = g_filtered_data, width = 8, height = 8, dpi = 300, limitsize = FALSE)

#Create box-violin plots for BoxCox transformed and z-scaled DSS2 to show there are still significant differences after the transformation
g_for_transformed_dss2 <- grouped_ggbetweenstats( # paired samples
  data = copy_all_response_metrics,
  x = lab,
  y = DSS2_boxcox_sclaed2,
  grouping.var = drug,
  type = "nonparametric", # for wilcoxon
  centrality.plotting = FALSE, # remove median
  ggstatsplot.layer = TRUE,
  xlab = "",        
  ylab = "Box-Cox transformed Z-scaled DSS2",
  results.subtitle = TRUE,  # Removes default subtitle
  #subtitle = "{test} (p = {p})", 
  p.adjust.method = "bonferroni",
  pairwise.display = "significant", #"significant"
  facet_wrap.args = list(scales = "fixed", strip.position = "top"),
  violin.args = list(alpha = 0),
  boxplot.args = list(alpha = 0),
  ggsignif.args = list(textsize = 3, tip_length = 0.005, step_increase = 0.05),
  p.value.label.args = list(
    parse = TRUE
  ),
  ggplot.component = list(
    coord_cartesian(ylim = c(-4, 4)),
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))),
    scale_color_manual(values = custom_colors), 
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

#Updating labels and subtitles of plots
#47 drugs including 001, RAD
for (i in 1:47) {
  #Check if the current g[[i]] contains annotations
  if (!is.null(g_for_transformed_dss2[[i]]$layers[[4]]$stat_params$annotations)) {
    
    #Extract the annotations for this g[[i]]
    annotations <- g_for_transformed_dss2[[i]]$layers[[4]]$stat_params$annotations
    
    #Modify the annotations by removing unwanted parts
    annotations_modified <- unlist(lapply(annotations, function(annot) {
      #Remove LaTeX parts (italic(p) and adj.)
      annot <- gsub("list\\(italic\\(p\\)\\[\\\"Bonferroni\\\"\\s*-\\s*adj\\.\\] == \"", "list(Bonferroni == ", annot)
      annot <- gsub("\"\\)", ")", annot)  # Remove the trailing ')'
      p_val_str <- sub(".*==\\s*", "", annot)
      p_val_str <- sub("\\)", "", p_val_str)
      p_val <- as.numeric(p_val_str)
      # Modify subtitle_mod based on value
      if (!is.na(p_val) && p_val < 0.001) {
        annot <- "Bonferroni < 0.001"
      } 
      
      return(annot)
    }))
    
    # Update the annotations in g[[i]]
    g_for_transformed_dss2[[i]]$layers[[4]]$stat_params$annotations <- annotations_modified
  }
  if(!is.null(g_for_transformed_dss2[[i]]$labels$subtitle)){
    subtitle <- g_for_transformed_dss2[[i]]$labels$subtitle
    subtitle_str <- deparse(subtitle)  
    subtitle_str <- paste(subtitle_str, collapse = " ")  # Ensure it's a single string
    
    # Extract only "Kruskal-Wallis: italic(p) == value"
    subtitle_mod <- str_replace(subtitle_str, 'list\\(chi\\["Kruskal-Wallis"\\]\\^2.*?italic\\(p\\) ==      "', "Kruskal-Wallis: p = ")
    subtitle_mod <- str_replace(subtitle_mod, '", widehat\\(epsilon\\)\\["ordinal"\\]\\^2.*', "")
    p_val <- as.numeric(sub("Kruskal-Wallis: p = ", "", subtitle_mod))
    
    # Modify subtitle_mod based on value
    if (!is.na(p_val) && p_val < 0.001) {
      subtitle_mod <- "Kruskal-Wallis: p < 0.001"
    }
    g_for_transformed_dss2[[i]]$labels$subtitle <- subtitle_mod
  }
}

#Print and save the plot
print(g_for_transformed_dss2)
ggsave(paste0(figure_output,"Difference_box_cox_sclaed_DSS2_per_drug.png"), plot = g_for_transformed_dss2, width = 50, height = 90, units = "cm", limitsize = FALSE)


















# beat_aml_karolinska <- rbind(subset(dss_karolinska_common_drugs, select = c(drug, Patient.num, DSS2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, DSS2, lab)))
# beat_aml_karolinska_FIMM <- rbind(beat_aml_karolinska, subset(dss_fimm_common_drugs, select = c(drug, Patient.num, DSS2, lab)))
# beat_aml_karolinska_FIMM_enserink<- rbind(beat_aml_karolinska_FIMM, subset(dss_github_enserink, select = c(drug, Patient.num, DSS2,lab)))
# 
# clinical_info <- all_response_metrics %>% mutate(Patient.num = ifelse(lab == "Karolinska", paste0(as.character(Patient.num), "_", as.character(sample)), as.character(Patient.num)))
# 
# beat_aml_karolinska_FIMM_enserink <- merge(beat_aml_karolinska_FIMM_enserink, unique(clinical_info[,c('Patient.num', 'Disease.status', 'Tissue', 'sample')]), by='Patient.num')
# beat_aml_karolinska_FIMM_enserink[beat_aml_karolinska_FIMM_enserink$lab == "Karolinska",]
# specific_order <- c("Helsinki","Beat AML", "Karolinska", "Oslo")
# 
# # Convert the column to a factor with the specific order
# beat_aml_karolinska_FIMM_enserink$lab <- factor(beat_aml_karolinska_FIMM_enserink$lab, levels = specific_order)
# 
# 
# beat_aml_karolinska_FIMM_enserink <- subset(beat_aml_karolinska_FIMM_enserink, drug != 'Tanespimycin')
# beat_aml_karolinska_FIMM_enserink <- subset(beat_aml_karolinska_FIMM_enserink, drug != 'Sapanisertib')
# beat_aml_karolinska_FIMM_enserink <- subset(beat_aml_karolinska_FIMM_enserink, drug != '001, RAD')
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- pivot_wider(beat_aml_karolinska_FIMM_enserink, names_from = drug, values_from = DSS2)
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- as.data.frame(beat_aml_karolinska_FIMM_enserink_for_heatmap)
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- beat_aml_karolinska_FIMM_enserink_for_heatmap[!is.na(beat_aml_karolinska_FIMM_enserink_for_heatmap$Patient.num),]
# 
# rownames(beat_aml_karolinska_FIMM_enserink_for_heatmap) <- beat_aml_karolinska_FIMM_enserink_for_heatmap$Patient.num
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- as.data.frame(beat_aml_karolinska_FIMM_enserink_for_heatmap)
# 
# ppca <- pca(beat_aml_karolinska_FIMM_enserink_for_heatmap[,-c(1, 2)], method="ppca", nPcs=3, seed=123)
# # Get explained variance
# ppca_var <- summary(ppca)
# 
# # Extract percent variance for PC1 and PC2
# percent_var <- round(ppca@R2 * 100, 1)
# 
# # Create a data frame for plotting
# var_df <- data.frame(
#   PC = paste0("PC", seq_along(percent_var)),
#   Variance = percent_var
# )
# 
# ggplot(var_df, aes(x = PC, y = Variance)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   geom_text(aes(label = paste0(round(Variance, 1), "%")), vjust = -0.5, size = 4) +
#   labs(title = "Explained Variance by Principal Component",
#        y = "Percent Variance Explained",
#        x = "Principal Component") +
#   theme_minimal()
# 
# ## Get the estimated complete observations
# ppca_scores <- scores(ppca)
# nrow(ppca_scores)
# patient_clinical_ppca <- merge(ppca_scores, beat_aml_karolinska_FIMM_enserink_for_heatmap, by = 'row.names')
# 
# #patient_clinical_ppca <- merge(ppca_data, unique(all_response_metrics[,c('Patient.num', 'Disease.status', 'Tissue', 'sample')]), by='Patient.num')
# #patient_clinical_ppca %>% group_by(Tissue) %>% dplyr::summarize(c = n())
# #patient_clinical_ppca[is.na(patient_clinical_ppca$Tissue),]
# #all_response_metrics[is.na(all_response_metrics$Tissue),]
# clinical_ppca_plot <- ggplot(as.data.frame(patient_clinical_ppca), aes(x = PC1, y = PC2, color=`Disease.status`, shape=Tissue)) +
#   geom_point() +
#   labs(
#     x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
#     y = paste0("PC2 (", round(percent_var[2], 1), "%)"),
#     color = ""
#   ) +
#   #geom_text_repel(data = furthest_points, aes(label = Patient.num), size = 4, box.padding = 0.5, max.overlaps = 20) +
#   #scale_color_gradient(low = "lightblue", high = "blue") + 
#   #scale_color_manual(values = c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")) +
#   labs(color = "") + 
#   guides(
#     color = guide_legend(
#       override.aes = aes(shape = 15, size = 4, width = 1.5, height = 1),  # Adjust legend symbol size
#       label.spacing = unit(0.5, "cm")  # Reduce space between text and legend symbol
#     )
#   ) +  
#   #theme_minimal()+
#   theme(
#     panel.background = element_blank(), 
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 12, family = "Arial"),    # Center legend text
#     legend.position = "right",                                
#     legend.justification = "center",        # Ensure the legend is centered
#     #legend.box = "horizontal",
#     legend.key = element_blank(),                            # Remove key background
#     legend.spacing.x = unit(0.2, "cm"),  
#     legend.key.size = unit(0.5, "lines"),      # Adjust size of the colored squares
#     legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
#     plot.margin = margin(10, 10, 10, 10),
#     axis.title = element_text(size = 12, family = "Arial"),    # Axis title size
#     axis.text = element_text(size = 10, family = "Arial")
#   )
# #guides(color = guide_legend(override.aes = aes(label = "", alpha = 1)))
# 
# print(clinical_ppca_plot)


# #Heatmap all datasets 
# names(dss_karolinska)[names(dss_karolinska) == "dss2"] <- "DSS2"
# names(dss_github_enserink)[names(dss_github_enserink) == "dss2"] <- "DSS2"
# 
# beat_aml_karolinska <- rbind(subset(dss_karolinska_common_drugs, select = c(drug, Patient.num, DSS2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, DSS2, lab)))
# beat_aml_karolinska_FIMM <- rbind(beat_aml_karolinska, subset(dss_fimm_common_drugs, select = c(drug, Patient.num, DSS2, lab)))
# beat_aml_karolinska_FIMM_enserink<- rbind(beat_aml_karolinska_FIMM, subset(dss_enserink_common_drugs, select = c(drug, Patient.num, DSS2,lab)))
# beat_aml_karolinska_FIMM_enserink <- beat_aml_karolinska_FIMM_enserink[!is.na(beat_aml_karolinska_FIMM_enserink$Patient.num),]
# beat_aml_karolinska_FIMM_enserink <- beat_aml_karolinska_FIMM_enserink %>% mutate(drug = ifelse(drug == 'PONATINIB', 'Ponatinib', drug))
# 
# a_col <- unique(subset(beat_aml_karolinska_FIMM_enserink, select = c(lab, Patient.num)))
# a_col <- as.data.frame(a_col)
# rownames(a_col) <- a_col$Patient.num
# beat_aml_karolinska_FIMM_enserink <- subset(beat_aml_karolinska_FIMM_enserink, drug != '001, RAD')
# beat_aml_karolinska_FIMM_enserink$lab <- NULL
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- pivot_wider(beat_aml_karolinska_FIMM_enserink, names_from = Patient.num, values_from = DSS2)
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- as.data.frame(beat_aml_karolinska_FIMM_enserink_for_heatmap)
# rownames(beat_aml_karolinska_FIMM_enserink_for_heatmap) <- beat_aml_karolinska_FIMM_enserink_for_heatmap$drug
# beat_aml_karolinska_FIMM_enserink_for_heatmap$drug <- NULL
# beat_aml_karolinska_FIMM_enserink_for_heatmap$lab <- NULL
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- as.data.frame(beat_aml_karolinska_FIMM_enserink_for_heatmap)
# #length(missing_per_row)
# # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_karolinska_FIMM_enserink_for_heatmap, 1, function(row) sum(is.na(row)))
# #Remove row if total missing in row is higher than 200
# beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered <- beat_aml_karolinska_FIMM_enserink_for_heatmap[missing_per_row <= 800, ]
# # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_karolinska_FIMM_enserink_for_heatmap, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered <- beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered[,missing_per_col <= 25]
# #beat_aml_karolinska_for_heatmap_filtered[is.na(beat_aml_karolinska_for_heatmap_filtered)] <- 0
# 
# rownames(beat_aml_karolinska_FIMM_enserink_for_heatmap) <- tools::toTitleCase(rownames(beat_aml_karolinska_FIMM_enserink_for_heatmap))
# dss_long <- as.data.frame(beat_aml_karolinska_FIMM_enserink_for_heatmap)
# #dss_long$drug <- rownames(dss_long)
# dss_long <- gather(as.data.frame(dss_long), patient_id, dss2, 'kAML_001_fresh':'patient 29', factor_key=TRUE) %>% as.data.frame()
# a_col <- as.data.frame(a_col)
# a = inner_join(unique(subset(dss_long, select = patient_id)), unique(a_col[, c('Patient.num', 'lab')]), by = c('patient_id' = 'Patient.num'))
# a <- as.data.frame(unique(a[,c("patient_id", "lab")]))
# #rownames(a) <- a$patient_id
# #a$patient_id <- NULL
# a <- as.data.frame(a)
# heatmap_plots(as.data.frame(beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered), a_col=a, fontsize_row = FALSE,filename = "Desktop/UiO/Project 1/Figures/V3/beat_aml_karolinska_FIMM_enserink1.pdf", title = "DSS2 Heatmap of all drugs", c_method= "ward.D2", height=40, width=460)
# complexheatmap_plots(as.matrix(beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered), a_col=a, fontsize_row = FALSE,filename = "Desktop/UiO/Project 1/Figures/draw/beat_aml_karolinska_FIMM_enserink_complex.pdf", title = "DSS2 Heatmap of all drugs", c_method= "ward.D2", height=40, width=460)
# length(rownames(a))
# ncol(beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered)
# col_anno <- subset(a, patient_id != "AML_001" & patient_id != "AML_002", select="lab") %>% dplyr::rename(Lab = lab)
# col_anno$Lab <- factor(col_anno$Lab, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
# 
# a_col_h <- ComplexHeatmap::HeatmapAnnotation(df = col_anno, annotation_legend_param = list(
#   title = "Lab",
#   title_gp = gpar(fontsize = 11, fontface = "bold"),
#   labels_gp = gpar(fontsize = 11),
#   #grid_width = unit(10, "mm"),
#   #grid_height = unit(10, "mm"),
#   legend_direction = "vertical",
#   border = FALSE
# ),
# col = list(Lab= c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")), 
# show_annotation_name = FALSE
# )
# str(a)
# matrix_colnames <- colnames(beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered[,-c(1,2)])
# annotation_rownames <- rownames(a)
# beat_aml_karolinska_FIMM_enserink_for_heatmap
# 
# library(cluster)
# 
# # Compute distance matrix with Gower’s method
# dist_mat <- daisy(as.matrix(beat_aml_karolinska_FIMM_enserink_for_heatmap), metric = "gower")
# 
# # Perform clustering
# hclust_res <- hclust(dist_mat, method = "ward.D2")
# 
# 
# library(proxy)
# 
# # Custom function for distance that ignores NAs
# dist_no_na <- function(x) {
#   #d <- daisy(x, metric = "euclidean")
#   d <- as.dist(proxy::dist(x, method = "euclidean", pairwise = TRUE))  # pairwise handles NAs
#   d[is.na(d)] <- max(d, na.rm = TRUE)  # Replace NA distances with max value
#   as.dist(d)
# }
# 
# # Perform clustering while ignoring NAs
# row_clustering <- hclust(dist_no_na(as.matrix(beat_aml_karolinska_FIMM_enserink_for_heatmap)), method = "ward.D")
# col_clustering <- hclust(dist_no_na(t(as.matrix(beat_aml_karolinska_FIMM_enserink_for_heatmap))), method = "ward.D")
# 
# ComplexHeatmap::Heatmap(
#   as.matrix(beat_aml_karolinska_FIMM_enserink_for_heatmap), 
#   cluster_rows = row_clustering, 
#   cluster_columns = col_clustering, 
#   na_col = "grey")
# 
# 
# 
# h <- ComplexHeatmap::Heatmap(
#   as.matrix(beat_aml_karolinska_FIMM_enserink_for_heatmap),
#   name = "DSS2",  # Legend header
#   cluster_rows = TRUE,
#   cluster_columns = TRUE,
#   clustering_method_rows = "ward.D",
#   clustering_method_columns = "ward.D2",
#   row_names_gp = gpar(fontsize = 9),  # Row font size 10
#   column_names_gp = gpar(fontsize = 1),  # Column font size 10
#   show_row_dend = TRUE,
#   show_column_dend = TRUE,
#   show_column_names = FALSE,
#   top_annotation = a_col_h,  # Column annotations
#   left_annotation = NULL,  # Row annotations
#   na_col = "grey",
#   row_gap = unit(1, "cm"),             # Space between rows
#   column_gap = unit(1, "cm"),          # Space between columns
#   heatmap_legend_param = list(
#     title = "DSS2",  # Add header to the legend
#     title_gp = gpar(fontsize = 11, fontface = "bold"),  # Legend title font
#     labels_gp = gpar(fontsize = 10), 
#     grid_height = unit(8, "mm")
#   ),
#   #annotation_legend_param = list(title = ""), 
#   width = unit(8, "cm"),  # 12 #20
#   height = unit(12, "cm"),  # 9.6 #16
#   column_title = "",
#   column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
#   row_dend_width = unit(2, "cm"),
#   clustering_distance_rows = "euclidean",
#   column_dend_height = unit(2, "cm"),
#   clustering_distance_columns = dist_no_na
#   )
# 
# plot(h)
# 
# pdf(file = "Desktop/UiO/Project 1/Figures/draw/beat_aml_karolinska_FIMM_enserink_complex_for_schematic.pdf", width = 8, height = 12)  # PDF output
# #pushViewport(viewport(width = 10, height = 8))
# plot(
#   h,
#   heatmap_legend_side = "right",
#   annotation_legend_side = "right",
#   merge_legends = TRUE, 
#   padding = unit(c(0, 0, 0, 0), "mm")
# )
# dev.off()  # Close the graphic device
# 
# plot(
#   h,
#   heatmap_legend_side = "right",
#   annotation_legend_side = "right",
#   merge_legends = TRUE
# )
# 
# # Clamp the matrix values to 0 - 50 for color scaling
# clamped_matrix <- pmin(pmax(beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered, 0), 50)
# 
# # Define the color function for DSS2 with min = 0, mid = 25, and max = 50
# dss2_col_fun <- circlize::colorRamp2(
#   c(0, 25, 50),    # min = 0, mid = 25, max = 50
#   c("blue", "white", "red")  # Colors: Blue for low, White for mid, Red for high
# )
# 
# # Create the DSS2 legend
# dss2_legend <- ComplexHeatmap::Legend(
#   title = "DSS2",
#   col_fun = dss2_col_fun,
#   title_gp = gpar(fontsize = 12, fontface = "bold"),
#   labels_gp = gpar(fontsize = 12),
#   grid_height = unit(10, "mm")
# )
# 
# # Combine legends as before
# combined_legends <- ComplexHeatmap::packLegend(lab_legend, dss2_legend, direction = "vertical")
# 
# # Heatmap creation using the clamped matrix
# f <- ComplexHeatmap::Heatmap(
#   as.matrix(clamped_matrix),  # Use the clamped matrix for the heatmap
#   name = "DSS2",
#   cluster_rows = TRUE,
#   cluster_columns = TRUE,
#   clustering_method_rows = "ward.D2",
#   clustering_method_columns = "ward.D2",
#   row_names_gp = gpar(fontsize = 12, fontface = "bold", family = 'Arial'),
#   column_names_gp = gpar(fontsize = 12, fontface = "bold", family = 'Arial'),
#   show_row_dend = TRUE,
#   show_column_dend = TRUE,
#   show_column_names = FALSE,
#   top_annotation = a_col_h,  # Column annotations
#   left_annotation = NULL,
#   na_col = "grey",
#   row_gap = unit(1, "cm"),
#   column_gap = unit(1, "cm"),
#   col = dss2_col_fun,  # Use the adjusted color function
#   show_heatmap_legend = FALSE,
#   width = unit(12.01, "cm"),
#   height = unit(11.21, "cm"),
#   column_title = "",
#   column_title_gp = gpar(fontsize = 12, fontface = "bold"),
#   row_dend_width = unit(2, "cm"),
#   clustering_distance_rows = "euclidean",
#   column_dend_height = unit(2, "cm"),
#   clustering_distance_columns = "euclidean"
# )
# 
# # Draw heatmap with only the combined legends
# ComplexHeatmap::draw(f, annotation_legend_list = combined_legends)
# ComplexHeatmap::draw(f)
###

# #All datasets 
# beat_aml_karolinska <- rbind(subset(dss_karolinska_common_drugs, select = c(drug, Patient.num, DSS2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, DSS2, lab)))
# beat_aml_karolinska_FIMM <- rbind(beat_aml_karolinska, subset(dss_fimm_common_drugs, select = c(drug, Patient.num, DSS2, lab)))
# beat_aml_karolinska_FIMM_enserink<- rbind(beat_aml_karolinska_FIMM, subset(dss_github_enserink, select = c(drug, Patient.num, DSS2,lab)))
# specific_order <- c("Helsinki","Beat AML", "Karolinska", "Oslo")
# 
# # Convert the column to a factor with the specific order
# beat_aml_karolinska_FIMM_enserink$lab <- factor(beat_aml_karolinska_FIMM_enserink$lab, levels = specific_order)
# 
# 
# beat_aml_karolinska_FIMM_enserink <- subset(beat_aml_karolinska_FIMM_enserink, drug != 'Tanespimycin')
# beat_aml_karolinska_FIMM_enserink <- subset(beat_aml_karolinska_FIMM_enserink, drug != 'Sapanisertib')
# beat_aml_karolinska_FIMM_enserink <- subset(beat_aml_karolinska_FIMM_enserink, drug != '001, RAD')
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- pivot_wider(beat_aml_karolinska_FIMM_enserink, names_from = drug, values_from = DSS2)
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- as.data.frame(beat_aml_karolinska_FIMM_enserink_for_heatmap)
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- beat_aml_karolinska_FIMM_enserink_for_heatmap[!is.na(beat_aml_karolinska_FIMM_enserink_for_heatmap$Patient.num),]
# 
# rownames(beat_aml_karolinska_FIMM_enserink_for_heatmap) <- beat_aml_karolinska_FIMM_enserink_for_heatmap$Patient.num
# beat_aml_karolinska_FIMM_enserink_for_heatmap <- as.data.frame(beat_aml_karolinska_FIMM_enserink_for_heatmap)
# 
# # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_karolinska_FIMM_enserink_for_heatmap, 1, function(row) sum(is.na(row)))
# #Remove row if total missing in row is higher than 200
# beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered <- beat_aml_karolinska_FIMM_enserink_for_heatmap[missing_per_row <= 30, ]
# # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered <- beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered[,missing_per_col <= 300]
# beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered_x <- na.omit(beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered)
# pca_plots(beat_aml_karolinska_FIMM_enserink_for_heatmap_filtered_x, title = "PCA Plot", "~/Desktop/UiO/Project 1/Figures/V3/Karolinska_BeatAML_PCA_Plot1.png", 'Patient.num', 'lab')
# 
# ppca(beat_aml_karolinska_FIMM_enserink_for_heatmap[,-c(1, 2)], nPcs=2)
# 
# nb = estim_ncpPCA(beat_aml_karolinska_FIMM_enserink_for_heatmap[,-c(1, 2)],ncp.max=5)
# res.comp = imputePCA(beat_aml_karolinska_FIMM_enserink_for_heatmap[,-c(1, 2)],ncp=2)
# res.pca = PCA(res.comp$completeObs) 
# 
# res.comp = MIPCA(beat_aml_karolinska_FIMM_enserink_for_heatmap[,-c(1, 2)], ncp = nb$ncp, nboot = 1000)
# plot(res.comp) 
# 
# PPCA <- ppca(res.comp$completeObs, ndim=2)
# plot(PPCA$Y,  pch=19, col=label, main="PCA")
# 
# pca_result <- prcomp(res.comp$completeObs, scale. = TRUE)
# pc_scores <- as.data.frame(pca_result$x[, 1:2])  # Using only the first two principal components
# pca_data <- cbind(pc_scores, lab = beat_aml_karolinska_FIMM_enserink_for_heatmap$lab)
# nrow(pca_data)
# nrow(beat_aml_karolinska_FIMM_enserink_for_heatmap)
# ggplot(pca_data, aes(x = PC1, y = PC2, color=lab)) +
#   geom_point() +
#   #geom_text(nudge_x = 0.2, nudge_y = 0.2, size = 3) +
#   scale_color_manual(values = c("Helsinki" = "#e41a1c", "BeatAML" = "#4daf4a", "Karolinska" = "#ff7f00", "Oslo"= "#377eb8")) +
#   labs(color = "Labs") + theme(legend.text = element_text(size = 12),       # Increase legend label size
#                                legend.title = element_text(size = 15), 
#                                plot.title = element_text(size = 20))
# 
# ppca <- pca(beat_aml_karolinska_FIMM_enserink_for_heatmap[,-c(1, 2)], method="ppca", nPcs=3, seed=123)
# # Get explained variance
# ppca_var <- summary(ppca)
# 
# # Extract percent variance for PC1 and PC2
# percent_var <- round(ppca@R2 * 100, 1)
# 
# # Create a data frame for plotting
# var_df <- data.frame(
#   PC = paste0("PC", seq_along(percent_var)),
#   Variance = percent_var
# )
# 
# ggplot(var_df, aes(x = PC, y = Variance)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   geom_text(aes(label = paste0(round(Variance, 1), "%")), vjust = -0.5, size = 4) +
#   labs(title = "Explained Variance by Principal Component",
#        y = "Percent Variance Explained",
#        x = "Principal Component") +
#   theme_minimal()
# 
# ## Get the estimated complete observations
# ppca_scores <- scores(ppca)
# nrow(ppca_scores)
# ppca_data <- merge(ppca_scores, beat_aml_karolinska_FIMM_enserink_for_heatmap, by = 'row.names')
# 
# centroids <- ppca_data %>%
#   group_by(lab) %>%
#   summarise(centroid_PC1 = mean(PC1), centroid_PC2 = mean(PC2))
# 
# pca_data <- ppca_data %>%
#   left_join(centroids, by = "lab") %>%
#   mutate(distance_to_centroid = sqrt((PC1 - centroid_PC1)^2 + (PC2 - centroid_PC2)^2))
# 
# furthest_points <- pca_data %>%
#   group_by(lab) %>%
#   top_n(3, distance_to_centroid) 
# library(ggrepel)
# library(ggnewscale)
# ppca_data$lab <- gsub('BeatAML', 'Beat AML', ppca_data$lab)
# ppca_data$lab <- factor(ppca_data$lab, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
# beat_aml_karolinska_FIMM_enserink_avg <- beat_aml_karolinska_FIMM_enserink %>% group_by(Patient.num) %>% dplyr::summarize('Average DSS2' = mean(DSS2))
# ppca_data <- ppca_data %>% left_join(beat_aml_karolinska_FIMM_enserink_avg, by = 'Patient.num')
# ppca_data$`Average DSS2` <- as.numeric(ppca_data$`Average DSS2`)
# ppca_data$Lab <- as.factor(ppca_data$lab)
# ppca_plot <- ggplot(as.data.frame(ppca_data), aes(x = PC1, y = PC2, color=Lab)) +
#   geom_point() +
#   labs(
#     x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
#     y = paste0("PC2 (", round(percent_var[2], 1), "%)"),
#     color = ""
#   ) +
#   #geom_text_repel(data = furthest_points, aes(label = Patient.num), size = 4, box.padding = 0.5, max.overlaps = 20) +
#   #scale_color_gradient(low = "lightblue", high = "blue") + 
#   scale_color_manual(values = c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")) +
#   labs(color = "") + 
#   guides(
#     color = guide_legend(
#       override.aes = aes(shape = 15, size = 4, width = 1.5, height = 1),  # Adjust legend symbol size
#       label.spacing = unit(0.5, "cm")  # Reduce space between text and legend symbol
#     )
#   ) +  
#   #theme_minimal()+
#   theme(
#     panel.background = element_blank(), 
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.text = element_text(size = 12, family = "Arial"),    # Center legend text
#     legend.position = "top",                                
#     legend.justification = "center",        # Ensure the legend is centered
#     legend.box = "horizontal",
#     legend.key = element_blank(),                            # Remove key background
#     legend.spacing.x = unit(0.2, "cm"),  
#     legend.key.size = unit(0.5, "lines"),      # Adjust size of the colored squares
#     legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
#     plot.margin = margin(10, 10, 10, 10),
#     axis.title = element_text(size = 12, family = "Arial"),    # Axis title size
#     axis.text = element_text(size = 10, family = "Arial")
#   )
# #guides(color = guide_legend(override.aes = aes(label = "", alpha = 1)))
# 
# print(ppca_plot)
# 
# ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/v1/PPCA_plot_all_data_all_labs.png", plot = ppca_plot, width = 10, height = 8, dpi = 300)
# midpoint_value <- median(ppca_data[ppca_data$lab == "Beat AML",]$`Average DSS2`, na.rm = TRUE)
# ppca_data$`Average DSS2` <- as.numeric(ppca_data$`Average DSS2`)
# dss2_ppca_plot <- ggplot(as.data.frame(ppca_data), aes(x = PC1, y = PC2, color = `Average DSS2`)) +
#   # 1st Category (A)
#   geom_point(aes(color = ifelse(lab == "Beat AML", `Average DSS2`, NA)), size = 4) +  
#   scale_color_gradient2(low = "#ccece6", mid = "#8dd3c7", high = "#145A32", midpoint = median(ppca_data[ppca_data$lab == "Beat AML",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Beat AML", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +  
#   new_scale_color() +  # Reset color scale for other points
#   
#   # 2nd Category (B)
#   geom_point(aes(color = ifelse(lab == "Oslo", `Average DSS2`, NA)), size = 4) +  
#   scale_color_gradient2(low = "#ffe6c7", mid = "#fdb462", high = "#a34700",midpoint = median(ppca_data[ppca_data$lab == "Oslo",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Oslo", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +
#   new_scale_color() +
#   
#   # 3rd Category (C)
#   geom_point(aes(color = ifelse(lab == "Helsinki", `Average DSS2`, NA)), size = 4) +  
#   scale_color_gradient2(low = "#ffd2cc", mid = "#fb8072", high = "#a12b1d", ,midpoint = median(ppca_data[ppca_data$lab == "Helsinki",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Helsinki", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +
#   new_scale_color() +
#   
#   # 4th Category (D)
#   geom_point(aes(color = ifelse(lab == "Karolinska", `Average DSS2`, NA)), size = 4) +  
#   scale_color_gradient2(low = "#d6ecff", mid = "#80b1d3", high = "#23537d", ,midpoint = median(ppca_data[ppca_data$lab == "Karolinska",]$`Average DSS2`, na.rm = TRUE), na.value = NA, name = "Karolinska", guide = guide_colorbar(title.position = "top", barwidth = 1, barheight = 4)) +
#   #geom_point() +
#   #geom_text_repel(data = furthest_points, aes(label = Patient.num), size = 4, box.padding = 0.5, max.overlaps = 20) +
#   #scale_color_gradient(low = "blue", high = "red") + 
#   #scale_shape_manual(values = c("Beat AML" = 16,"Oslo"= 17,"Helsinki" = 18, "Karolinska" = 19)) +
#   labs(color = "") + 
#   theme_minimal()+
#   theme(
#     panel.background = element_blank(), 
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.text = element_text(size = 12, family = "Arial"),    # Center legend text
#     #legend.position = "right",                                
#     legend.spacing.x = unit(0.2, "cm"),  
#     legend.key.size = unit(0.5, "lines"),      # Adjust size of the colored squares
#     legend.margin = margin(t = 0, b = 0, l = 0, r = 0),                            
#     plot.margin = margin(10, 10, 10, 10),
#     axis.title = element_text(size = 12, family = "Arial"),    # Axis title size
#     axis.text = element_text(size = 10, family = "Arial")
#   ) 
# dss2_ppca_plot
# ggsave("/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/avg_dss_score_colorsPPCA_plot_all_data_all_labs_v1.png", plot = dss2_ppca_plot, width = 10, height = 8, dpi = 300)


###
# #PPCA plot Other metrics
# beat_aml_karolinska_metics <- rbind(subset(dss_karolinska_common_drugs, select = c(drug, Patient.num, DSS1, DSS2, DSS3, AUC, IC50, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, DSS1, DSS2, DSS3, AUC, IC50, lab)))
# beat_aml_karolinska_FIMM_metrics <- rbind(beat_aml_karolinska_metics, subset(dss_fimm_common_drugs, select = c(drug, Patient.num, DSS1, DSS2, DSS3, AUC, IC50, lab)))
# beat_aml_karolinska_FIMM_enserink_metrics<- rbind(beat_aml_karolinska_FIMM_metrics, subset(dss_enserink_common_drugs, select = c(drug, Patient.num, DSS1, DSS2, DSS3, AUC, IC50,lab)))
# 
# beat_aml_karolinska_FIMM_enserink_metrics <- subset(beat_aml_karolinska_FIMM_enserink_metrics, drug != 'Tanespimycin')
# beat_aml_karolinska_FIMM_enserink_metrics <- subset(beat_aml_karolinska_FIMM_enserink_metrics, drug != 'Sapanisertib')
# beat_aml_karolinska_FIMM_enserink_metrics <- subset(beat_aml_karolinska_FIMM_enserink_metrics, drug != '001, RAD')
# 
# beat_aml_karolinska_FIMM_enserink_metrics$drug <- tools::toTitleCase(beat_aml_karolinska_FIMM_enserink_metrics$drug)
# 
# specific_order <- c("Beat AML", "Oslo", "Helsinki", "Karolinska")
# 
# # Convert the column to a factor with the specific order
# beat_aml_karolinska_FIMM_enserink_metrics$lab <- factor(beat_aml_karolinska_FIMM_enserink_metrics$lab, levels = specific_order)
# 
# beat_aml_karolinska_FIMM_enserink_metrics$log10IC50 <- -log10(beat_aml_karolinska_FIMM_enserink_metrics$IC50)
# #metrics <- c("DSS1", "DSS2", "DSS3", "AUC", "IC50") #"log10IC50"
# metrics <- c("DSS2")
# for (m in metrics){
#   beat_aml_karolinska_FIMM_enserink_for_heatmap_metrics <- pivot_wider(beat_aml_karolinska_FIMM_enserink_metrics[,c("drug", "Patient.num", "lab", m)], names_from = drug, values_from = m)
#   beat_aml_karolinska_FIMM_enserink_for_heatmap_metrics <- subset(beat_aml_karolinska_FIMM_enserink_for_heatmap_metrics, !is.na(Patient.num))
#   print(beat_aml_karolinska_FIMM_enserink_for_heatmap_metrics)
#   
#   ppca_metrics <- pca(beat_aml_karolinska_FIMM_enserink_for_heatmap_metrics[,-c(1, 2)], method="ppca", nPcs=3, seed=123)
#   ## Get the estimated complete observations
#   ppca_scores_metrics <- scores(ppca_metrics)
#   ppca_data_metrics <- merge(ppca_scores_metrics, beat_aml_karolinska_FIMM_enserink_for_heatmap_metrics, by = 'row.names')
#   ppca_data_metrics$lab <- gsub('BeatAML', 'Beat AML', ppca_data_metrics$lab)
#   ppca_data_metrics$lab <- factor(ppca_data_metrics$lab, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
#   ppca_data_metrics$Lab <- as.factor(ppca_data_metrics$lab)
#   ppca_plot <- ggplot(as.data.frame(ppca_data_metrics), aes(x = PC1, y = PC2, color=Lab)) +
#     geom_point() +
#     #geom_text_repel(data = furthest_points, aes(label = Patient.num), size = 4, box.padding = 0.5, max.overlaps = 20) +
#     #scale_color_gradient(low = "lightblue", high = "blue") + 
#     scale_color_manual(values = c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")) +
#     labs(color = "") + 
#     guides(
#       color = guide_legend(
#         override.aes = aes(shape = 15, size = 4, width = 1.5, height = 1),  # Adjust legend symbol size
#         label.spacing = unit(0.5, "cm")  # Reduce space between text and legend symbol
#       )
#     ) +  
#     #theme_minimal()+
#     theme(
#       panel.background = element_blank(), 
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.text = element_text(size = 12, family = "Arial"),    # Center legend text
#       legend.position = "top",                                
#       legend.justification = "center",        # Ensure the legend is centered
#       legend.box = "horizontal",
#       legend.key = element_blank(),                            # Remove key background
#       legend.spacing.x = unit(0.2, "cm"),  
#       legend.key.size = unit(0.5, "lines"),      # Adjust size of the colored squares
#       legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
#       plot.margin = margin(10, 10, 10, 10),
#       axis.title = element_text(size = 12, family = "Arial"),    # Axis title size
#       axis.text = element_text(size = 10, family = "Arial")
#     )
#   #guides(color = guide_legend(override.aes = aes(label = "", alpha = 1)))
#   
#   #print(ppca_plot)
#   #print(subset(ppca_data_metrics, Lab == "Beat AML"))
#   heatmap_metrics <- as.data.frame(beat_aml_karolinska_FIMM_enserink_for_heatmap_metrics)
#   print(rownames(heatmap_metrics))
#   
#   rownames(heatmap_metrics) <- make.names(heatmap_metrics$Patient.num, unique = TRUE)
#   #heatmap_metrics <- t(as.matrix(heatmap_metrics[,-c(1, 2)]))
#   dist_no_na <- function(x) {
#     #d <- daisy(x, metric = "euclidean")
#     d <- as.dist(proxy::dist(x, method = "euclidean", pairwise = TRUE))  # pairwise handles NAs
#     d[is.na(d)] <- max(d, na.rm = TRUE)  # Replace NA distances with max value
#     as.dist(d)
#   }
#   
#   # Perform clustering while ignoring NAs
#   #row_clustering <- hclust(dist_no_na(heatmap_metrics), method = "ward.D")
#   #col_clustering <- hclust(dist_no_na(heatmap_metrics), method = "ward.D")
#   
#   col_anno <- subset(heatmap_metrics, select="lab") %>% dplyr::rename(Study = lab) %>% as.data.frame()
#   col_anno <- factor(col_anno$Study, levels = c("Beat AML", "Oslo", "Helsinki", "Karolinska"))
#   levels(col_anno)
#   a_col_h <- ComplexHeatmap::HeatmapAnnotation(Study = col_anno, 
#                                                col = list(Study= c("Beat AML"="#8dd3c7", "Oslo"="#fdb462", "Helsinki"= "#fb8072","Karolinska"="#80b1d3")), show_annotation_name = FALSE, gp = gpar(fonface = "plain"))
#   print(a_col_h)
#   
#   
#   
#   # if(m == "IC50"){
#   #   heatmap_metrics <- t(scale(heatmap_metrics[,-c(1, 2)]))
#   # }
#   # else{
#   #   heatmap_metrics <- t(as.matrix(heatmap_metrics[,-c(1, 2)]))
#   # }
#   heatmap_metrics <- t(as.matrix(heatmap_metrics[,-c(1, 2)]))
#   
#   print(m)
#   m_header <- if(m == "DSS2"){expression("DSS"[2])}
#   else if(m == "DSS1"){expression("DSS"[1])}
#   else if (m == "DSS3"){expression("DSS"[3])}
#   else if(m == "IC50"){expression("IC"[50])}
#   else{m}
#   
#   library(circlize)
#   print(m_header)
#   h <- ComplexHeatmap::Heatmap(
#     heatmap_metrics,
#     name = m,  # Legend header
#     col = if(m == "IC50"){colorRamp2(c(10000, 1000, 0), c("blue", "white", "red"))},
#     cluster_rows = TRUE,
#     cluster_columns = TRUE,
#     clustering_method_rows = "ward.D",
#     clustering_method_columns = "ward.D2",
#     row_names_gp = gpar(fontsize = 9),  # Row font size 10
#     column_names_gp = gpar(fontsize = 1),  # Column font size 10
#     show_row_dend = TRUE,
#     show_column_dend = TRUE,
#     show_column_names = FALSE,
#     top_annotation = a_col_h,  # Column annotations
#     left_annotation = NULL,  # Row annotations
#     na_col = "grey",
#     row_gap = unit(1, "cm"),             # Space between rows
#     column_gap = unit(1, "cm"),          # Space between columns
#     heatmap_legend_param = list(
#       title = m_header,  # Add header to the legend
#       title_gp = gpar(fontsize = 11, fontface = "bold"),  # Legend title font
#       labels_gp = gpar(fontsize = 10), 
#       grid_height = unit(8, "mm")
#     ),
#     #annotation_legend_param = list(title = ""), 
#     width = unit(8, "cm"),  # 12 #20
#     height = unit(11, "cm"),  # 9.6 #16
#     column_title = " ",
#     column_title_gp = gpar(fontsize = 20), 
#     row_dend_width = unit(2, "cm"),
#     clustering_distance_rows = "euclidean",
#     column_dend_height = unit(2, "cm"),
#     clustering_distance_columns = dist_no_na
#   )
#   
#   plot(
#     h,
#     heatmap_legend_side = "right",
#     annotation_legend_side = "right",
#     merge_legends = TRUE, 
#     padding = unit(c(0, 0, 0, 0), "mm")
#   )
#   library(patchwork)
#   heat_g <- grid.grabExpr(plot(
#     h,
#     heatmap_legend_side = "right",
#     annotation_legend_side = "right",
#     merge_legends = TRUE, 
#     padding = unit(c(0, 0, 0, 0), "mm")
#   ))
#   grid.arrange(ppca_plot, heat_g, ncol = 2, widths = c(1, 1))
#   ggsave(paste0('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/', m, '_heatmap_and_ppca.png'), plot = grid.arrange(ppca_plot, heat_g, ncol = 2, widths = c(1, 1)), width = 6.47 * 2, height = 5.94 * 1, dpi = 300)
# }
# 
# grid.arrange(dss2_ppca_plot, heat_g, ncol = 2, widths = c(1, 1))
# ggsave(paste0('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/draw/DSS2_extra_heatmap_and_ppca.png'), plot = grid.arrange(ppca_plot, heat_g, ncol = 2, widths = c(1, 1)), width = 6.47 * 2, height = 5.94 * 1, dpi = 300)
# unique(beat_aml_karolinska_FIMM_enserink_metrics[beat_aml_karolinska_FIMM_enserink_metrics$lab == 'Beat AML', c("IC50", "Patient.num")])





###
# #Importing the datasets 
# #Enserink 
# dss_github_enserink <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink.csv')
# dss_github_enserink$lab <- 'Oslo'
# #colnames(dss_github_enserink)[colnames(dss_github_enserink) == "DSS2"] <- "dss2"
# #Enserink full drug set 
# dss_github_enserink_full_set_drugs <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_enserink_full_drug_set_testing.csv')
# 
# #BeatAML
# dss_github_beat_aml <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_beat_aml_full_drug_set_v1.csv')
# dss_github_beat_aml$lab <- 'Beat AML'
# 
# #Karolinska 
# karolinska_dss2_org <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/karolinska_institute_dss2.csv')
# karolinska_dss2_org$...1 <- NULL
# karolinska_dss2 <- gather(karolinska_dss2_org, patient_id, dss2, 'AML_001':'AML_347', factor_key=TRUE)
# karolinska_dss2 <- as.data.frame(karolinska_dss2)
# colnames(karolinska_dss2)[colnames(karolinska_dss2) == "patient_id"] <- "Patient.num"
# 
# ###
# karolinska_fresh <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_karolinska_full_drug_set_fresh.csv')
# #updating karolindka dataset
# colnames(karolinska_fresh)[colnames(karolinska_fresh) == "DRUG_NAME"] <- "drug"
# colnames(karolinska_fresh)[colnames(karolinska_fresh) == "DSS"] <- "DSS2"
# #karolinska_fresh <- karolinska_fresh %>% mutate(Patient.num = paste(Patient.num, '_fresh'))
# karolinska_fresh$sample <- 'fresh'
# karolinska_fresh <- karolinska_fresh %>% mutate(Patient.num = paste0(Patient.num, '_', sample)) %>% as.data.frame()
# karolinska_frosen <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_karolinska_full_drug_set_frosen.csv')
# #updating karolindka dataset
# colnames(karolinska_frosen)[colnames(karolinska_frosen) == "DRUG_NAME"] <- "drug"
# colnames(karolinska_frosen)[colnames(karolinska_frosen) == "DSS"] <- "DSS2"
# #karolinska_frosen <- karolinska_frosen %>% mutate(Patient.num = paste(Patient.num, '_frozen'))
# karolinska_frosen$sample <- 'frozen'
# karolinska_frosen <- karolinska_frosen %>% mutate(Patient.num = paste0(Patient.num, '_', sample)) %>% as.data.frame()
# dss_karolinska_github <- rbind(karolinska_fresh, karolinska_frosen)
# dss_karolinska_github$lab <- 'Karolinska'
# dss_karolinska_github$Tissue <- NA
# dss_karolinska_github$Disease.status <- NA
# dss_karolinska_github <- dss_karolinska_github %>% mutate(Patient.num = paste0('k',dss_karolinska_github$Patient.num))
# dss_karolinska_github$drug[dss_karolinska_github$drug=='Chidamide'] <- 'Tucidinostat'
# dss_karolinska_github$drug[dss_karolinska_github$drug=='GDC-0994'] <- 'Ravoxertinib'
# dss_karolinska_github$drug[dss_karolinska_github$drug=='GSK525762'] <- 'Molibresib'
# dss_karolinska_github$drug[dss_karolinska_github$drug=='KPT-8602'] <- 'Eltanexor'
# dss_karolinska_github$drug[dss_karolinska_github$drug=='MLN-0128'] <- 'Sapanisertib'
# dss_karolinska_github$drug[dss_karolinska_github$drug=='MLN1117'] <- 'Serabelisib'
# dss_karolinska_github$drug[dss_karolinska_github$drug=='NVP-ABL001'] <- 'Asciminib'
# dss_karolinska_github$drug[dss_karolinska_github$drug=='NVP-BGJ398'] <- 'Infigratinib'
# dss_karolinska_github$drug[dss_karolinska_github$drug=='ONO-4059'] <- 'Tirabrutinib'
# dss_karolinska_github$drug[dss_karolinska_github$drug=='VX-745'] <- 'Neflamapimod'
# dss_karolinska_github <- subset(dss_karolinska_github, Patient.num != 'kAML_009_frozen')
# ###
# 
# 
# #FIMM 
# #Their calculations
# fimm_dss2 <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_drug_response.csv')
# fimm_dss2 <- gather(fimm_dss2, patient_id, dss2, 'AML_084_04':'Healthy_17', factor_key=TRUE)
# fimm_dss2 <- as.data.frame(fimm_dss2)
# colnames(fimm_dss2)[colnames(fimm_dss2) == "patient_id"] <- "Patient.num"
# 
# #My calculations from the github package
# dss_github_fimm <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_fimm_full_drug_set.csv')
# 
# #Common drug names
# #Getting the common drugs across the 4 datasets
# common_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/common_drugnames_pubchem.csv")
# common_drugs$`Unnamed: 0_x` <- NULL
# common_drugs$`Unnamed: 0_y` <- NULL
# common_drugs$`Unnamed: 0` <- NULL
# 
# #Joining DSS calculations and common drug dataset
# #All datasets are left with the 47 common drugs after the join 
# dss_beat_aml_common_drugs <- inner_join(common_drugs, dss_github_beat_aml, by=c("beat_aml_drug_name"="drug"))
# colnames(dss_beat_aml_common_drugs)[colnames(dss_beat_aml_common_drugs) == "pubchem_drug_name"] <- "drug"
# dim(dss_beat_aml_common_drugs)
# 
# dss_karolinska_common_drugs <- inner_join(common_drugs, dss_karolinska_github, by=c("karolinska_drug_name"="drug"))
# colnames(dss_karolinska_common_drugs)[colnames(dss_karolinska_common_drugs) == "pubchem_drug_name"] <- "drug"
# dss_karolinska_common_drugs$lab <- 'Karolinska'
# 
# #Their calculations
# dss_fimm <- inner_join(common_drugs, fimm_dss2, by=c("fimm_drug_name"="Drug_name"))
# colnames(dss_fimm)[colnames(dss_fimm) == "pubchem_drug_name"] <- "drug"
# dss_fimm$lab <- 'Helsinki'
# 
# #My calculations from the github package
# dss_fimm_common_drugs <- inner_join(common_drugs, dss_github_fimm, by=c("fimm_drug_name"="drug"))
# colnames(dss_fimm_common_drugs)[colnames(dss_fimm_common_drugs) == "pubchem_drug_name"] <- "drug"
# dss_fimm_common_drugs$lab <- 'Helsinki'
# 
# #Enserink common drugs 
# dss_enserink_common_drugs <- inner_join(common_drugs, dss_github_enserink_full_set_drugs, by=c("enserink_drug_name"="drug"))
# colnames(dss_enserink_common_drugs)[colnames(dss_enserink_common_drugs) == "pubchem_drug_name"] <- "drug"
# dss_enserink_common_drugs$lab <- "Oslo"
###

# 
###
# #Heatmap single datasets 
# #FIMM 
# dss_github_fimm <- subset(dss_github_fimm, !is.na(Patient.num))
# fimm_for_heatmap <- pivot_wider(subset(dss_github_fimm, select = c (Patient.num, DSS2, drug)), names_from = Patient.num, values_from = DSS2)
# fimm_for_heatmap <- as.data.frame(fimm_for_heatmap)
# rownames(fimm_for_heatmap) <- fimm_for_heatmap$drug
# fimm_for_heatmap$drug <- NULL
# heatmap_plots(data.frame(fimm_for_heatmap),NaN, fontsize_row = FALSE, fontsize = 26, filename = "Desktop/UiO/Project 1/Figures/V2_prepare_article/fimm_full_drug_set_1.pdf", title = "Heatmap all drugs - FIMM")
# complexheatmap_plots(as.matrix(fimm_for_heatmap),NaN, fontsize_row = FALSE, fontsize = 26, filename = "Desktop/UiO/Project 1/Figures/V3/fimm_full_drug_set_complex.pdf", title = "Heatmap all drugs - FIMM", height=290, width=80, c_method= "ward.D2")
# 
# 
# fimm_common_drugs_heatmap <- pivot_wider(subset(dss_fimm_common_drugs, select = c(drug, Patient.num, dss2)), names_from = Patient.num, values_from = dss2)
# fimm_common_drugs_heatmap <- as.data.frame(fimm_common_drugs_heatmap)
# rownames(fimm_common_drugs_heatmap) <- fimm_common_drugs_heatmap$drug
# fimm_common_drugs_heatmap$drug <- NULL
# fimm_common_drugs_heatmap <- as.data.frame(fimm_common_drugs_heatmap)
# #fimm_common_drugs_heatmap <- na.omit(fimm_common_drugs_heatmap)
# heatmap_plots(as.data.frame(fimm_common_drugs_heatmap),a_col=NaN, fontsize = 26, filename = "Desktop/UiO/Project 1/Figures/V2_prepare_article/FIMM_github_calculations_common_drugs1.png", title = "Heatmap common drugs - FIMM")
# 
# #Karolinska 
# karolinska_dss2_org_heat <- pivot_wider(subset(dss_karolinska_github, select = c (Patient.num, DSS2, drug)), names_from = Patient.num, values_from = DSS2)
# karolinska_dss2_org_heat <- as.data.frame(karolinska_dss2_org_heat)
# rownames(karolinska_dss2_org_heat) <- karolinska_dss2_org_heat$drug
# karolinska_dss2_org_heat$drug <- NULL
# karolinska_dss2_org_heat <- as.data.frame(karolinska_dss2_org_heat)
# heatmap_plots(karolinska_dss2_org_heat,NaN, fontsize = 26, filename = "Desktop/UiO/Project 1/Figures/V3/karolinska_github_calculations_all_drugs.pdf", title = "DSS2 Heatmap of all drugs \n Karolinska", c_method= "ward.D2", height=295, width=140)
# complexheatmap_plots(as.matrix(karolinska_dss2_org_heat),NaN, fontsize = 26, filename = "Desktop/UiO/Project 1/Figures/V3/karolinska_github_calculations_all_drugs_complex.pdf", title = "DSS2 Heatmap of all drugs \n Karolinska", c_method= "ward.D2", height=295, width=140)
# 
# 
# #Common drugs 
# karolinska_dss2_org_heat_common <- pivot_wider(subset(dss_karolinska_common_drugs, select = c (Patient.num, DSS2, drug)), names_from = Patient.num, values_from = DSS2)
# karolinska_dss2_org_heat_common <- as.data.frame(karolinska_dss2_org_heat_common)
# rownames(karolinska_dss2_org_heat_common) <- karolinska_dss2_org_heat_common$drug
# karolinska_dss2_org_heat_common$drug <- NULL
# heatmap_plots(karolinska_dss2_org_heat_common,NaN, fontsize = 26, filename = "Desktop/UiO/Project 1/Figures/V2_prepare_article/karolinska_github_calculations_common_drugs.pdf", title = "Heatmap common drugs - Karolinska")
# 
# 
# 
# 
# #Enserink 
# #Full drug set
# dss_github_enserink_full_set_drugs_for_heatmap <- pivot_wider(subset(dss_github_enserink_full_set_drugs, select = c(drug, Patient.num, DSS2)), names_from = Patient.num, values_from = DSS2)
# dss_github_enserink_full_set_drugs_for_heatmap <- as.data.frame(dss_github_enserink_full_set_drugs_for_heatmap)
# rownames(dss_github_enserink_full_set_drugs_for_heatmap) <- dss_github_enserink_full_set_drugs_for_heatmap$drug
# dss_github_enserink_full_set_drugs_for_heatmap$drug <- NULL
# dss_github_enserink_full_set_drugs_for_heatmap <- as.data.frame(dss_github_enserink_full_set_drugs_for_heatmap)
# heatmap_plots(dss_github_enserink_full_set_drugs_for_heatmap, NaN, fontsize = 26, fontsize_row = FALSE, filename = 'Desktop/UiO/Project 1/Figures/V3/Oslo_github_calculations_full_drug_set.pdf', title = "DSS2 Heatmap of all drugs \n Oslo", c_method= "ward.D2", height=200, width=80)
# complexheatmap_plots(as.matrix(dss_github_enserink_full_set_drugs_for_heatmap), NaN, fontsize = 26, fontsize_row = FALSE, filename = 'Desktop/UiO/Project 1/Figures/V3/Oslo_github_calculations_full_drug_set_complex.pdf', title = "DSS2 Heatmap of all drugs \n Oslo", c_method= "ward.D2", height=200, width=80)
# 
# 
# #Common drugs 
# dss_github_enserink_common_set_drugs_for_heatmap <- pivot_wider(subset(dss_enserink_common_drugs, select = c(drug, Patient.num, DSS2)), names_from = Patient.num, values_from = DSS2)
# dss_github_enserink_common_set_drugs_for_heatmap <- as.data.frame(dss_github_enserink_common_set_drugs_for_heatmap)
# rownames(dss_github_enserink_common_set_drugs_for_heatmap) <- dss_github_enserink_common_set_drugs_for_heatmap$drug
# dss_github_enserink_common_set_drugs_for_heatmap$drug <- NULL
# dss_github_enserink_common_set_drugs_for_heatmap <- as.data.frame(dss_github_enserink_common_set_drugs_for_heatmap)
# heatmap_plots(dss_github_enserink_common_set_drugs_for_heatmap, NaN, fontsize = 26, fontsize_row = FALSE, filename = 'Desktop/UiO/Project 1/Figures/V2_prepare_article/Enserink_github_calculations_common_drugs.pdf', title = "Heatmap common drugs - Enserink")
# 
# 
# #BeatAML
# #full drug set
# beat_aml_for_heatmap <- pivot_wider(subset(dss_github_beat_aml, select = c(drug, Patient.num, DSS2)), names_from = Patient.num, values_from = DSS2)
# beat_aml_for_heatmap <- as.data.frame(beat_aml_for_heatmap)
# rownames(beat_aml_for_heatmap) <- beat_aml_for_heatmap$drug
# beat_aml_for_heatmap$drug <- NULL
# beat_aml_for_heatmap <- as.data.frame(beat_aml_for_heatmap)
# # # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_for_heatmap, 1, function(row) sum(is.na(row)))
# # #Remove row if total missing in row is higher than 200
# df_filtered <- beat_aml_for_heatmap[missing_per_row <= 200, ]
# # # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_for_heatmap, 2, function(row) sum(is.na(row)))
# # #Remove col if total missing in col is higher than 100
# df_filtered <- df_filtered[,missing_per_col <= 100]
# # #beat_aml_for_heatmap[is.na(beat_aml_for_heatmap)] <- 0
# heatmap_plots(data.frame(df_filtered),NaN, fontsize = 26, filename = "Desktop/UiO/Project 1/Figures/V3/BeatAML_github_calculations_full_drug_set.pdf", title = "DSS2 Heatmap of all drugs \n BeatAML", c_method= "ward.D2", height=70, width=310)
# complexheatmap_plots(as.matrix(df_filtered), fontsize = 26, filename = "Desktop/UiO/Project 1/Figures/draw/BeatAML_github_calculations_full_drug_set_complex.pdf", title = "DSS2 Heatmap of all drugs \n BeatAML", c_method= "ward.D2", height=45, width=140)
# 
# #Common drugs
# beat_aml_for_heatmap_common <- pivot_wider(subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, DSS2)), names_from = Patient.num, values_from = DSS2)
# beat_aml_for_heatmap_common <- as.data.frame(beat_aml_for_heatmap_common)
# rownames(beat_aml_for_heatmap_common) <- beat_aml_for_heatmap_common$drug
# beat_aml_for_heatmap_common$drug <- NULL
# beat_aml_for_heatmap_common <- as.data.frame(beat_aml_for_heatmap_common)
# # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_for_heatmap_common, 1, function(row) sum(is.na(row)))
# #Remove row if total missing in row is higher than 200
# df_filtered <- beat_aml_for_heatmap_common[missing_per_row <= 500, ]
# # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_for_heatmap_common, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# df_filtered <- beat_aml_for_heatmap_common[,missing_per_col <= 45]
# 
# subset(dss_beat_aml_common_drugs, is.finite(DSS2))
# heatmap_plots(data.frame(beat_aml_for_heatmap_common),NaN, fontsize = 26, filename = "Desktop/UiO/Project 1/Figures/V2_prepare_article/BeatAML_github_calculations_common_drugs.pdf", title = "Heatmap common drugs - BeatAML", c_row=FALSE, c_col=FALSE)
#### 
# 
### 
# #Heatmap two datasets
# #FIMM-Karolinska
# fimm_karolinska<- rbind(subset(dss_karolinska, select = c(drug, Patient.num, DSS2,lab)), subset(dss_fimm_common_drugs, select = c(drug, Patient.num, dss2,lab)))
# 
# a_col <- unique(subset(fimm_karolinska, select = c(lab, Patient.num)))
# a_col <- as.data.frame(a_col)
# rownames(a_col) <- a_col$Patient.num
# a_col$Patient.num <- NULL
# 
# fimm_karolinska$lab <- NULL
# all_dss_for_heatmap <- pivot_wider(fimm_karolinska, names_from = Patient.num, values_from = dss2)
# all_dss_for_heatmap <- as.data.frame(all_dss_for_heatmap)
# rownames(all_dss_for_heatmap) <- all_dss_for_heatmap$drug
# all_dss_for_heatmap$drug <- NULL
# all_dss_for_heatmap <- as.data.frame(all_dss_for_heatmap)
# heatmap_plots(all_dss_for_heatmap,a_col, fontsize_row = TRUE, filename = "Desktop/UiO/Project 1/Figures/fimm_karolinska_bæææ.pdf")
# 
# #BeatAML-Enserink
# beat_aml_enserink <- rbind(subset(dss_github_enserink, select = c(drug, Patient.num, dss2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# 
# a_col <- unique(subset(beat_aml_enserink, select = c(lab, Patient.num)))
# a_col <- as.data.frame(a_col)
# rownames(a_col) <- a_col$Patient.num
# 
# beat_aml_enserink$lab <- NULL
# beat_aml_enserink_for_heatmap <- pivot_wider(beat_aml_enserink, names_from = Patient.num, values_from = dss2)
# beat_aml_enserink_for_heatmap <- as.data.frame(beat_aml_enserink_for_heatmap)
# rownames(beat_aml_enserink_for_heatmap) <- beat_aml_enserink_for_heatmap$drug
# beat_aml_enserink_for_heatmap$drug <- NULL
# beat_aml_enserink_for_heatmap$lab <- NULL
# beat_aml_enserink_for_heatmap <- as.data.frame(beat_aml_enserink_for_heatmap)
# 
# # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_enserink_for_heatmap, 1, function(row) sum(is.na(row)))
# #Remove row if total missing in row is higher than 200
# beat_aml_enserink_for_heatmap_filtered <- beat_aml_enserink_for_heatmap[missing_per_row <= 200, ]
# # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_enserink_for_heatmap_filtered, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# beat_aml_enserink_for_heatmap_filtered <- beat_aml_enserink_for_heatmap_filtered[,missing_per_col <= 30]
# #beat_aml_for_heatmap[is.na(beat_aml_for_heatmap)] <- 0
# dss_long <- as.data.frame(beat_aml_enserink_for_heatmap_filtered)
# #dss_long$drug <- rownames(dss_long)
# dss_long <- gather(as.data.frame(dss_long), patient_id, dss2, 'Patient 100':'2752_BA3448D_BA3448R', factor_key=TRUE) %>% as.data.frame()
# a_col <- as.data.frame(a_col)
# a = inner_join(unique(subset(dss_long, select = patient_id)), unique(a_col[, c('Patient.num', 'lab')]), by = c('patient_id' = 'Patient.num'))
# a <- as.data.frame(unique(a[,c("patient_id", "lab")]))
# rownames(a) <- a$patient_id
# a$patient_id <- NULL
# a <- as.data.frame(a)
# 
# a_color = list(
#   lab = c("BeatAML" = "#F4A582", "Enserink" = "#92C5DE"))
# 
# heatmap_plots(as.data.frame(beat_aml_enserink_for_heatmap_filtered), a_col=a, a_color = a_color, fontsize_row = FALSE,filename = "Desktop/UiO/Project 1/Figures/enserink_beat_aml.pdf")
# 
# #FIMM-Enserink
# fimm_enserink <- rbind(subset(dss_github_enserink, select = c(drug, Patient.num, dss2, lab)), subset(dss_fimm_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# 
# a_col <- unique(subset(fimm_enserink, select = c(lab, Patient.num)))
# a_col <- as.data.frame(a_col)
# rownames(a_col) <- a_col$Patient.num
# a_col$Patient.num <- NULL
# 
# fimm_enserink$lab <- NULL
# fimm_enserink_for_heatmap <- pivot_wider(fimm_enserink, names_from = Patient.num, values_from = dss2)
# fimm_enserink_for_heatmap <- as.data.frame(fimm_enserink_for_heatmap)
# rownames(fimm_enserink_for_heatmap) <- fimm_enserink_for_heatmap$drug
# fimm_enserink_for_heatmap$drug <- NULL
# fimm_enserink_for_heatmap$lab <- NULL
# fimm_enserink_for_heatmap <- as.data.frame(fimm_enserink_for_heatmap)
# 
# heatmap_plots(fimm_enserink_for_heatmap,a_col, fontsize_row = TRUE, filename = "Desktop/UiO/Project 1/Figures/fimm_enserink.pdf")
# 
# #Karolinska-Enserink
# karolinska_enserink<- rbind(subset(dss_karolinska, select = c(drug, Patient.num, DSS2,lab)), subset(dss_github_enserink, select = c(drug, Patient.num, DSS2,lab)))
# 
# a_col <- unique(subset(karolinska_enserink, select = c(lab, Patient.num)))
# a_col <- as.data.frame(a_col)
# rownames(a_col) <- a_col$Patient.num
# a_col$Patient.num <- NULL
# 
# karolinska_enserink$lab <- NULL
# karolinska_enserink_for_heatmap <- pivot_wider(karolinska_enserink, names_from = Patient.num, values_from = DSS2)
# karolinska_enserink_for_heatmap <- as.data.frame(karolinska_enserink_for_heatmap)
# rownames(karolinska_enserink_for_heatmap) <- karolinska_enserink_for_heatmap$drug
# karolinska_enserink_for_heatmap$drug <- NULL
# karolinska_enserink_for_heatmap <- as.data.frame(karolinska_enserink_for_heatmap)
# heatmap_plots(karolinska_enserink_for_heatmap,a_col, fontsize_row = TRUE, filename = "Desktop/UiO/Project 1/Figures/karolinska_Enserink_test.pdf", height = 70, width = 130)
# 
# #FIMM-BeatML
# beat_aml_FIMM <- rbind(subset(dss_fimm_common_drugs, select = c(drug, Patient.num, dss2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# 
# a_col <- unique(subset(beat_aml_FIMM, select = c(lab, Patient.num)))
# a_col <- as.data.frame(a_col)
# rownames(a_col) <- a_col$Patient.num
# 
# beat_aml_FIMM$lab <- NULL
# beat_aml_FIMM_for_heatmap <- pivot_wider(beat_aml_FIMM, names_from = Patient.num, values_from = dss2)
# beat_aml_FIMM_for_heatmap <- as.data.frame(beat_aml_FIMM_for_heatmap)
# rownames(beat_aml_FIMM_for_heatmap) <- beat_aml_FIMM_for_heatmap$drug
# beat_aml_FIMM_for_heatmap$drug <- NULL
# beat_aml_FIMM_for_heatmap$lab <- NULL
# beat_aml_FIMM_for_heatmap <- as.data.frame(beat_aml_FIMM_for_heatmap)
# 
# # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_FIMM_for_heatmap, 1, function(row) sum(is.na(row)))
# #Remove row if total missing in row is higher than 200
# beat_aml_FIMM_for_heatmap_filtered <- beat_aml_FIMM_for_heatmap[missing_per_row <= 200, ]
# # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_FIMM_for_heatmap_filtered, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# beat_aml_FIMM_for_heatmap_filtered <- beat_aml_FIMM_for_heatmap_filtered[,missing_per_col <= 25]
# #beat_aml_FIMM_for_heatmap_filtered[is.na(beat_aml_FIMM_for_heatmap_filtered)] <- 0
# dss_long <- as.data.frame(beat_aml_FIMM_for_heatmap_filtered)
# #dss_long$drug <- rownames(dss_long)
# dss_long <- gather(as.data.frame(dss_long), patient_id, dss2, 'BERG_470012_03022016_9999_BM':'2554_BA2594D_BA2594R', factor_key=TRUE) %>% as.data.frame()
# a_col <- as.data.frame(a_col)
# a = inner_join(unique(subset(dss_long, select = patient_id)), unique(a_col[, c('Patient.num', 'lab')]), by = c('patient_id' = 'Patient.num'))
# a <- as.data.frame(unique(a[,c("patient_id", "lab")]))
# rownames(a) <- a$patient_id
# a$patient_id <- NULL
# a <- as.data.frame(a)
# 
# heatmap_plots(as.data.frame(beat_aml_FIMM_for_heatmap_filtered), a_col=a, fontsize_row = FALSE,filename = "Desktop/UiO/Project 1/Figures/FIMM_beat_aml.pdf")
# 
# 
# #Karolinska-BeatAML
# beat_aml_karolinska <- rbind(subset(dss_karolinska, select = c(drug, Patient.num, dss2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# 
# a_col <- unique(subset(beat_aml_karolinska, select = c(lab, Patient.num)))
# a_col <- as.data.frame(a_col)
# rownames(a_col) <- a_col$Patient.num
# 
# beat_aml_karolinska$lab <- NULL
# beat_aml_karolinska_for_heatmap <- pivot_wider(beat_aml_karolinska, names_from = Patient.num, values_from = dss2)
# beat_aml_karolinska_for_heatmap <- as.data.frame(beat_aml_karolinska_for_heatmap)
# rownames(beat_aml_karolinska_for_heatmap) <- beat_aml_karolinska_for_heatmap$drug
# beat_aml_karolinska_for_heatmap$drug <- NULL
# beat_aml_karolinska_for_heatmap$lab <- NULL
# beat_aml_karolinska_for_heatmap <- as.data.frame(beat_aml_karolinska_for_heatmap)
# 
# # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_karolinska_for_heatmap, 1, function(row) sum(is.na(row)))
# #Remove row if total missing in row is higher than 200
# beat_aml_karolinska_for_heatmap_filtered <- beat_aml_karolinska_for_heatmap[missing_per_row <= 200, ]
# # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_karolinska_for_heatmap_filtered, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# beat_aml_karolinska_for_heatmap_filtered <- beat_aml_karolinska_for_heatmap_filtered[,missing_per_col <= 25]
# #beat_aml_karolinska_for_heatmap_filtered[is.na(beat_aml_karolinska_for_heatmap_filtered)] <- 0
# dss_long <- as.data.frame(beat_aml_karolinska_for_heatmap_filtered)
# #dss_long$drug <- rownames(dss_long)
# dss_long <- gather(as.data.frame(dss_long), patient_id, dss2, 'AML_001':'2554_BA2594D_BA2594R', factor_key=TRUE) %>% as.data.frame()
# a = inner_join(unique(subset(dss_long, select = patient_id)), unique(a_col[, c('Patient.num', 'lab')]), by = c('patient_id' = 'Patient.num'))
# a <- as.data.frame(unique(a[,c("patient_id", "lab")]))
# rownames(a) <- a$patient_id
# a$patient_id <- NULL
# a <- as.data.frame(a)
# 
# heatmap_plots(as.data.frame(beat_aml_karolinska_for_heatmap_filtered), a_col=a, fontsize_row = FALSE,filename = "Desktop/UiO/Project 1/Figures/Karolinska_beat_aml.pdf")
###
#PCA plot two datasets
#FIMM-Karolinska
# fimm_karolinska<- rbind(subset(dss_karolinska, select = c(drug, Patient.num, dss2,lab)), subset(dss_fimm_common_drugs, select = c(drug, Patient.num, dss2,lab)))
# all_dss_for_heatmap <- pivot_wider(fimm_karolinska, names_from = drug, values_from = dss2)
# all_dss_for_heatmap <- as.data.frame(all_dss_for_heatmap)
# 
# na_counts <- colSums(is.na(all_dss_for_heatmap))
# df_clean <-all_dss_for_heatmap[, na_counts < 30]
# df_clean <- df_clean[complete.cases(df_clean), ]
# dim(df_clean)
# 
# all_dss_for_heatmap <- na.omit(df_clean)
# pca_plots(all_dss_for_heatmap, title = "PCA Plot", "~/Desktop/UiO/Project 1/Figures/FIMM_Karolinska.png", 'Patient.num', 'lab')
# 
# #Enserink-BeatAML
# beat_aml_enserink <- rbind(subset(dss_github_enserink, select = c(drug, Patient.num, dss2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# beat_aml_enserink_for_heatmap <- pivot_wider(beat_aml_enserink, names_from = drug, values_from = dss2)
# beat_aml_enserink_for_heatmap <- as.data.frame(beat_aml_enserink_for_heatmap)
# rownames(beat_aml_enserink_for_heatmap) <- beat_aml_enserink_for_heatmap$drug
# beat_aml_enserink_for_heatmap$drug <- NULL
# beat_aml_enserink_for_heatmap <- as.data.frame(beat_aml_enserink_for_heatmap)
# 
# # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_enserink_for_heatmap, 1, function(row) sum(is.na(row)))
# #Remove row if total missing in row is higher than 200
# beat_aml_enserink_for_heatmap_filtered <- beat_aml_enserink_for_heatmap[missing_per_row <= 30, ]
# # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_enserink_for_heatmap_filtered, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# beat_aml_enserink_for_heatmap_filtered <- beat_aml_enserink_for_heatmap_filtered[,missing_per_col <= 200]
# beat_aml_enserink_for_heatmap_filtered <- na.omit(beat_aml_enserink_for_heatmap_filtered)
# #all_dss_for_heatmap <- na.omit(beat_aml_enserink_for_heatmap_filtered)
# pca_plots(beat_aml_enserink_for_heatmap_filtered, title = "PCA Plot - Enserink and BeatAML", "~/Desktop/UiO/Project 1/Figures/Enserink_BeatAML_PCA_Plot.png", 'Patient.num', 'lab')
# 
# #FIMM-Enserink
# fimm_enserink <- rbind(subset(dss_github_enserink, select = c(drug, Patient.num, dss2, lab)), subset(dss_fimm_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# 
# fimm_enserink_for_heatmap <- pivot_wider(fimm_enserink, names_from = drug, values_from = dss2)
# fimm_enserink_for_heatmap <- as.data.frame(fimm_enserink_for_heatmap)
# rownames(fimm_enserink_for_heatmap) <- fimm_enserink_for_heatmap$drug
# fimm_enserink_for_heatmap$drug <- NULL
# fimm_enserink_for_heatmap <- as.data.frame(fimm_enserink_for_heatmap)
# 
# # Check how many values are missing in each col
# missing_per_col <- apply(fimm_enserink_for_heatmap, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# fimm_enserink_for_heatmap <- fimm_enserink_for_heatmap[,missing_per_col <= 100]
# fimm_enserink_for_heatmap <- na.omit(fimm_enserink_for_heatmap)
# pca_plots(fimm_enserink_for_heatmap, title = "PCA Plot - Enserink and FIMM", "~/Desktop/UiO/Project 1/Figures/Enserink_FIMM_PCA_Plot.png", 'Patient.num', 'lab')
# 
# 
# #Karolinska-Enserink
# karolinska_enserink<- rbind(subset(dss_karolinska, select = c(drug, Patient.num, dss2,lab)), subset(dss_github_enserink, select = c(drug, Patient.num, dss2,lab)))
# 
# karolinska_enserink_for_heatmap <- pivot_wider(karolinska_enserink, names_from = drug, values_from = dss2)
# karolinska_enserink_for_heatmap <- as.data.frame(karolinska_enserink_for_heatmap)
# rownames(karolinska_enserink_for_heatmap) <- karolinska_enserink_for_heatmap$drug
# karolinska_enserink_for_heatmap$drug <- NULL
# karolinska_enserink_for_heatmap <- as.data.frame(karolinska_enserink_for_heatmap)
# 
# # Check how many values are missing in each col
# missing_per_col <- apply(karolinska_enserink_for_heatmap, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# karolinska_enserink_for_heatmap <- karolinska_enserink_for_heatmap[,missing_per_col <= 100]
# karolinska_enserink_for_heatmap <- na.omit(karolinska_enserink_for_heatmap)
# pca_plots(karolinska_enserink_for_heatmap, title = "PCA Plot - Karolinska and Enserink", "~/Desktop/UiO/Project 1/Figures/Enserink_Karolinska_PCA_Plot.png", 'Patient.num', 'lab')
# 
# 
# #FIMM-BeatAML
# beat_aml_FIMM <- rbind(subset(dss_fimm_common_drugs, select = c(drug, Patient.num, dss2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# beat_aml_FIMM <- subset(beat_aml_FIMM, drug != 'Sapanisertib')
# beat_aml_FIMM_for_heatmap <- pivot_wider(beat_aml_FIMM, names_from = drug, values_from = dss2)
# beat_aml_FIMM_for_heatmap <- as.data.frame(beat_aml_FIMM_for_heatmap)
# rownames(beat_aml_FIMM_for_heatmap) <- beat_aml_FIMM_for_heatmap$Patient.num
# beat_aml_FIMM_for_heatmap <- as.data.frame(beat_aml_FIMM_for_heatmap)
# 
# # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_FIMM_for_heatmap, 1, function(row) sum(is.na(row)))
# #Remove row if total missing in row is higher than 200
# beat_aml_FIMM_for_heatmap_filtered <- beat_aml_FIMM_for_heatmap[missing_per_row <= 30, ]
# # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_FIMM_for_heatmap_filtered, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# beat_aml_FIMM_for_heatmap_filtered <- beat_aml_FIMM_for_heatmap_filtered[,missing_per_col <= 200]
# beat_aml_FIMM_for_heatmap_filtered <- na.omit(beat_aml_FIMM_for_heatmap_filtered)
# pca_plots(beat_aml_FIMM_for_heatmap_filtered, title = "PCA Plot - BeatAML and FIMM", "~/Desktop/UiO/Project 1/Figures/FIMM_BeatAML_PCA_Plot.png", 'Patient.num', 'lab')
# 
# 
# #Karolinska-BeatAML
# beat_aml_karolinska <- rbind(subset(dss_karolinska, select = c(drug, Patient.num, dss2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# beat_aml_karolinska <- subset(beat_aml_karolinska, drug != 'Tanespimycin')
# beat_aml_karolinska_for_heatmap <- pivot_wider(beat_aml_karolinska, names_from = drug, values_from = dss2)
# beat_aml_karolinska_for_heatmap <- as.data.frame(beat_aml_karolinska_for_heatmap)
# rownames(beat_aml_karolinska_for_heatmap) <- beat_aml_karolinska_for_heatmap$Patient.num
# beat_aml_karolinska_for_heatmap <- as.data.frame(beat_aml_karolinska_for_heatmap)
# 
# # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_karolinska_for_heatmap, 1, function(row) sum(is.na(row)))
# #Remove row if total missing in row is higher than 200
# beat_aml_karolinska_for_heatmap_filtered <- beat_aml_karolinska_for_heatmap[missing_per_row <= 30, ]
# # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_karolinska_for_heatmap_filtered, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# beat_aml_karolinska_for_heatmap_filtered <- beat_aml_karolinska_for_heatmap_filtered[,missing_per_col <= 200]
# beat_aml_karolinska_for_heatmap_filtered <- na.omit(beat_aml_karolinska_for_heatmap_filtered)
# pca_plots(beat_aml_karolinska_for_heatmap_filtered, title = "PCA Plot - BeatAML and FIMM", "~/Desktop/UiO/Project 1/Figures/Karolinska_BeatAML_PCA_Plot.png", 'Patient.num', 'lab')
###
# #Wilcoxon rank sum test - see weather there is a difference between the dss values for the two datasets
# #Enserink-BeatAML
# #Common drugs
# beat_aml_enserink_w <- rbind(subset(dss_github_enserink, select = c(drug, Patient.num, dss2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# beat_aml_enserink_w$binary_lab<- ifelse(beat_aml_enserink_w$lab == "Enserink", 1, 0)
# #ighv both single and combination drug
# beat_aml_enserink_wilcoxon_test <- wilcox_test(beat_aml_enserink_w, 'drug', "binary_lab", 'dss2')
# # Create volcano plot using the result dataframe
# volcano_beat_aml_enserink <- create_volcano_plot(beat_aml_enserink_wilcoxon_test, title ='Volcano Plot for IGHV mutations', filename='/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/enserink_beat_aml_wilcoxon.png')
# # Print the volcano plot
# print(volcano_beat_aml_enserink)
# 
# 
# fimm_enserink <- rbind(subset(dss_fimm, select = c(drug, Patient.num, dss2,lab)), subset(dss_github_enserink, select = c(drug, Patient.num, dss2,lab)))
# 
# fimm_enserink$binary_lab<- ifelse(fimm_enserink$lab == "Enserink", 1, 0)
# #ighv both single and combination drug
# ighv_df <- wilcox_test(fimm_enserink, 'drug', "binary_lab", 'dss2')
# # Create volcano plot using the result dataframe
# volcano_plot_ighv <- create_volcano_plot(ighv_df, title ='Volcano Plot for FIMM Enserink', filename='/Users/katarinawilloch/Desktop/UiO/Project 1/fimm_Enserink.png')
# # Print the volcano plot
# print(volcano_plot_ighv)
# 
# 
# #FIMM-Karolinska
# fimm_karolinska<- rbind(subset(dss_karolinska, select = c(drug, Patient.num, dss2,lab)), subset(dss_fimm_common_drugs, select = c(drug, Patient.num, dss2,lab)))
# fimm_karolinska$binary_lab<- ifelse(fimm_karolinska$lab == "FIMM", 1, 0)
# fimm_karolinska <- na.omit(fimm_karolinska)
# binary_lab_fimm_karolinska <- wilcox_test(subset(fimm_karolinska, drug != 'Sapanisertib'), 'drug', "binary_lab", 'dss2')
# # Create volcano plot using the result dataframe
# volcano_plot_fimm_karolinska <- create_volcano_plot(binary_lab_fimm_karolinska, title ='Volcano Plot Wilcoxon rank sum test - DSS2 Karolinska vs FIMM', filename='/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/FIMM_Karolinska_vulcano_plot.png')
# # Print the volcano plot
# print(volcano_plot_fimm_karolinska)
# 
# #FIMM-Enserink
# fimm_enserink <- rbind(subset(dss_github_enserink, select = c(drug, Patient.num, dss2, lab)), subset(dss_fimm_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# fimm_enserink$binary_lab<- ifelse(fimm_enserink$lab == "FIMM", 1, 0)
# fimm_enserink <- na.omit(fimm_enserink)
# binary_lab_fimm_enserink <- wilcox_test(fimm_enserink, 'drug', "binary_lab", 'dss2')
# # Create volcano plot using the result dataframe
# volcano_plot_fimm_enserink <- create_volcano_plot(binary_lab_fimm_enserink, title ='Volcano Plot Wilcoxon rank sum test - DSS2 Enserink vs FIMM', filename='/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/FIMM_Enserink_vulcano_plot.png')
# # Print the volcano plot
# print(volcano_plot_fimm_enserink)
# 
# 
# #Karolinska-Enserink
# karolinska_enserink<- rbind(subset(dss_karolinska, select = c(drug, Patient.num, dss2,lab)), subset(dss_github_enserink, select = c(drug, Patient.num, dss2,lab)))
# karolinska_enserink$binary_lab<- ifelse(karolinska_enserink$lab == "Enserink", 1, 0)
# karolinska_enserink <- na.omit(karolinska_enserink)
# 
# binary_lab_karolinska_enserink <- wilcox_test(karolinska_enserink, 'drug', "binary_lab", 'dss2')
# # Create volcano plot using the result dataframe
# volcano_plot_karolinska_enserink <- create_volcano_plot(binary_lab_karolinska_enserink, title ='Volcano Plot Wilcoxon rank sum test - DSS2 Karolinska vs Enserink', filename='/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/Karolinska_Enserink_vulcano_plot.png')
# # Print the volcano plot
# print(volcano_plot_karolinska_enserink)
# 
# #FIMM-BeatAML
# beat_aml_FIMM <- rbind(subset(dss_fimm_common_drugs, select = c(drug, Patient.num, dss2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# beat_aml_FIMM$binary_lab<- ifelse(beat_aml_FIMM$lab == "FIMM", 1, 0)
# beat_aml_FIMM <- na.omit(beat_aml_FIMM)
# 
# binary_lab_beat_aml_FIMM <- wilcox_test(beat_aml_FIMM, 'drug', "binary_lab", 'dss2')
# # Create volcano plot using the result dataframe
# volcano_plot_beat_aml_FIMM <- create_volcano_plot(binary_lab_beat_aml_FIMM, title ='Volcano Plot Wilcoxon rank sum test - DSS2 FIMM vs BeatAML', filename='/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/FIMM_BeatAML_vulcano_plot.png')
# # Print the volcano plot
# print(volcano_plot_beat_aml_FIMM)
# 
# 
# #Karolinska-BeatAML
# beat_aml_karolinska <- rbind(subset(dss_karolinska, select = c(drug, Patient.num, dss2, lab)), subset(dss_beat_aml_common_drugs, select = c(drug, Patient.num, dss2, lab)))
# beat_aml_karolinska$binary_lab<- ifelse(beat_aml_karolinska$lab == "Karolinska", 1, 0)
# beat_aml_karolinska <- na.omit(beat_aml_karolinska)
# 
# binary_lab_beat_aml_karolinska <- wilcox_test(beat_aml_karolinska, 'drug', "binary_lab", 'dss2')
# # Create volcano plot using the result dataframe
# volcano_plot_beat_aml_karolinska <- create_volcano_plot(binary_lab_beat_aml_karolinska, title ='Volcano Plot Wilcoxon rank sum test - DSS2 BeatAML vs Karolinska', filename='/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/Karolinska_BeatAML_vulcano_plot.png')
# # Print the volcano plot
# print(volcano_plot_beat_aml_karolinska)
###





# xicor <- function(X, Y, ties = TRUE){
#   n <- length(X)
#   r <- rank(Y[order(X)], ties.method = "random")
#   set.seed(42)
#   if(ties){
#     l <- rank(Y[order(X)], ties.method = "max")
#     return( 1 - n*sum( abs(r[-1] - r[-n]) ) / (2*sum(l*(n - l))) )
#   } else {
#     return( 1 - 3 * sum( abs(r[-1] - r[-n]) ) / (n^2 - 1) )    
#   }
# }
# 
# grid_plots <- function(df, col1, col2, label,title1, title2, filename){
#   library(gridExtra)
#   require(qqplotr)
#   library(ggplot2)
#   library(grid)
#   library("psychTools")
#   # Calculate R2 value
#   lm_model <- lm(get(col2) ~ get(col1), data = df)
#   r2 <- summary(lm_model)$r.squared
#   p_value <- summary(lm_model)$coefficients[2,4]
#   
#   # Calculate xicor value
#   xicor_value <- xicor(df[[col1]], df[[col2]])
#   # Line chart
#   p3 <- ggplot(df, aes(x = get(col1), y = get(col2))) +
#     geom_line(color = "gray20") +
#     xlab(paste0('', col1)) +
#     ylab(paste0('', col2))
#   
#   # Scatter plot with R2 value
#   p4 <- ggplot(df, aes(x = get(col1), y = get(col2), label=get(label))) +
#     geom_point(color = "grey20") +
#     geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
#     xlab(paste0('', col1)) +
#     ylab(paste0('', col2)) +
#     geom_text() +
#     geom_label(aes(x = Inf, y = Inf, label = paste("R² = ", round(r2, 2), "\nP-value = ", round(p_value,3),"\nXicor = ", round(xicor_value, 2))),
#                hjust = 1.1, vjust = 2, size = 15, color = "blue", fontface = "bold",
#                fill = "white", label.size = NA) 
#   
#   gg1 <- ggplot(data = df, mapping = aes(sample = get(col1))) +
#     geom_qq_band(bandType = "ks", mapping = aes(fill = "KS"), alpha = 0.5) +
#     geom_qq_band(bandType = "ts", mapping = aes(fill = "TS"), alpha = 0.5) +
#     geom_qq_band(bandType = "pointwise", mapping = aes(fill = "Normal"), alpha = 0.5) +
#     geom_qq_band(bandType = "boot", mapping = aes(fill = "Bootstrap"), alpha = 0.5) +
#     stat_qq_line() +
#     stat_qq_point() +
#     labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
#     scale_fill_discrete("Bandtype") +
#     ggtitle(title1)
#   
#   gg2 <- ggplot(data = df, mapping = aes(sample = get(col2))) +
#     geom_qq_band(bandType = "ks", mapping = aes(fill = "KS"), alpha = 0.5) +
#     geom_qq_band(bandType = "ts", mapping = aes(fill = "TS"), alpha = 0.5) +
#     geom_qq_band(bandType = "pointwise", mapping = aes(fill = "Normal"), alpha = 0.5) +
#     geom_qq_band(bandType = "boot", mapping = aes(fill = "Bootstrap"), alpha = 0.5) +
#     stat_qq_line() +
#     stat_qq_point() +
#     labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
#     scale_fill_discrete("Bandtype") +
#     ggtitle(title2)
#   
#   g <- arrangeGrob(p3, p4, gg1, gg2) 
#   ggsave(file=filename, g, units="cm", width = 100, height = 100)
#   return(g)
# }
# 
# 
# 
# beat_aml_calc <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_calc.csv')
# beat_aml_calc <- beat_aml_calc %>% mutate(patient_id = paste(dbgap_subject_id, dbgap_dnaseq_sample,dbgap_rnaseq_sample, sep = "_")) %>% as.data.frame()
# beat_aml_calc_wide <- pivot_wider(subset(beat_aml_calc, select = c(inhibitor, patient_id, auc)), names_from = patient_id, values_from = auc)
# dim(beat_aml_calc)
# common_drugs <- read_csv("~/Desktop/UiO/Project 1/Data/Initial cleansing/common_drugnames_pubchem.csv")
# common_drugs$`Unnamed: 0_x` <- NULL
# common_drugs$`Unnamed: 0_y` <- NULL
# common_drugs$`Unnamed: 0` <- NULL
# 
# beat_aml_scores <- inner_join(common_drugs, beat_aml_calc, by=c("beat_aml_drug_name"="inhibitor"))
# colnames(beat_aml_scores)[colnames(beat_aml_scores) == "pubchem_drug_name"] <- "drug"
# beat_aml_scores$dbgap_subject_id <- as.character(beat_aml_scores$dbgap_subject_id)
# 
# dss_github_beat_aml_full_set <- read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/dss_github_beat_aml_full_drug_set.csv')
# dss_github_beat_aml_full_set <- inner_join(common_drugs, dss_github_beat_aml_full_set, by=c("beat_aml_drug_name"="inhibitor"))
# colnames(dss_github_beat_aml_full_set)[colnames(dss_github_beat_aml_full_set) == "pubchem_drug_name"] <- "drug"
# colnames(dss_github_beat_aml_full_set)[colnames(dss_github_beat_aml_full_set) == "Patient.num"] <- "patient_id"
# #dss_github_beat_aml_full_set$dbgap_subject_id <- as.character(dss_github_beat_aml_full_set$dbgap_subject_id)
# 
# beat_aml_scores_and_dss <- inner_join(beat_aml_scores, dss_github_beat_aml_full_set, by=c("drug"="drug", "patient_id"="patient_id"), relationship = "many-to-many")
# 
# for (i in unique(beat_aml_scores_and_dss$drug)) {
#   drug_col <- i
#   print(drug_col)
#   df <- subset(beat_aml_scores_and_dss, drug == drug_col)
#   print(df)
#   grid_plots(df, 'auc', 'dss2', 'patient_id',"", "", paste0("~/Desktop/UiO/Project 1/Figures/BeatAML/", drug_col, "_Beataml_dss_correlation_auc_test.png"))
# }
# grid_plots(beat_aml_scores_and_dss, 'auc', 'DSS2', 'patient_id',"", "", paste0("~/Desktop/UiO/Project 1/Figures/BeatAML/Beataml_dss_correlation_aic_test.png"))
# 
# 
# beat_aml_clinical <- read_csv("Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_clinical.csv")
# beat_aml_calc <- read_csv("Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_calc.csv")
# beat_aml_clinical <- beat_aml_clinical %>%
#   mutate(patient_id = paste(dbgap_subject_id, dbgap_dnaseq_sample,dbgap_rnaseq_sample, sep = "_")) %>% as.data.frame()
# df <- inner_join(dss_github_beat_aml_full_set, beat_aml_clinical, by=c('patient_id'))
# 
# multiple_combinations <- df %>%
#   arrange(., drug) %>%
#   group_by(patient_id, drug) %>%
#   filter(n() > 1) %>%
#   ungroup() %>%
#   as.data.frame()
# 
# bm_p_blood <- subset(multiple_combinations, select = c(patient_id, specimenType, drug, DSS2)) %>%
#   pivot_wider(names_from = specimenType, values_from = DSS2, values_fn = mean) %>%
#   as.data.frame()
# bm_p_blood <- subset(bm_p_blood, select = -Leukapheresis)
# bm_p_blood <- na.omit(bm_p_blood)
# grid_plots(bm_p_blood, 'Bone Marrow Aspirate', 'Peripheral Blood', 'patient_id',"", "", paste0("~/Desktop/UiO/Project 1/Figures/BeatAML/p_blood_bm/Beataml_p_bllod_bm_v1.png"))
# 
# 
# d1 <- subset(dss_github_fimm, select = c(Patient.num, drug, DSS2))
# d1$dataset <- "github"
# colnames(d1)[colnames(d1) == "DSS2"] <- "dss2"
# colnames(fimm_dss2)[colnames(fimm_dss2) == "Drug_name"] <- "drug"
# d2 <- subset(fimm_dss2, select = c(drug, Patient.num, dss2))
# d2$dataset <- "org"
# df <- rbind(d1, d2)
# bm_p_blood <- subset(df, select = c(dataset, drug, dss2)) %>%
#   pivot_wider(names_from = dataset, values_from = dss2, values_fn = mean) %>%
#   as.data.frame()
# grid_plots(bm_p_blood, 'github', 'org', 'drug',"", "", paste0("~/Desktop/UiO/Project 1/Figures/BeatAML/p_blood_bm/V1/fimm_dss_github_vs_org.png"))
# 
# fimm_vs_enserink <- subset(fimm_enserink, select = c(lab, drug, dss2)) %>%
#   pivot_wider(names_from = lab, values_from = dss2, values_fn = mean) %>%
#   as.data.frame()
# grid_plots(fimm_vs_enserink, 'Enserink', 'FIMM', 'drug', "", "", "~/Desktop/UiO/Project 1/Figures/fimm_vs_Enserink_correlation_plot.png")
# 
###
# #Density explorations
# #Density heatmap
# 
# #BeatAML
# beat_aml_for_heatmap <- pivot_wider(subset(dss_github_beat_aml, select = c(drug, Patient.num, dss2)), names_from = Patient.num, values_from = dss2)
# beat_aml_for_heatmap <- as.data.frame(beat_aml_for_heatmap)
# rownames(beat_aml_for_heatmap) <- beat_aml_for_heatmap$drug
# beat_aml_for_heatmap$drug <- NULL
# beat_aml_for_heatmap <- as.data.frame(beat_aml_for_heatmap)
# # Check how many values are missing in each row
# missing_per_row <- apply(beat_aml_for_heatmap, 1, function(row) sum(is.na(row)))
# #Remove row if total missing in row is higher than 200
# df_filtered <- beat_aml_for_heatmap[missing_per_row <= 200, ]
# # Check how many values are missing in each col
# missing_per_col <- apply(beat_aml_for_heatmap, 2, function(row) sum(is.na(row)))
# #Remove col if total missing in col is higher than 100
# df_filtered <- df_filtered[,missing_per_col <= 100]
# beat_aml_for_heatmap[is.na(beat_aml_for_heatmap)] <- -1
# densityHeatmap(t(beat_aml_for_heatmap))
# 
# 
# m = matrix(rnorm(1000), nc = 10)
# lt = apply(m, 2, function(x) data.frame(density(x)[c("x", "y")]))
# ha = rowAnnotation(foo = anno_joyplot(lt, width = unit(4, "cm"), 
#                                       gp = gpar(fill = 1:10), transparency = 0.75))
