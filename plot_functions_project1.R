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

data_exploration <- function(df, output_loc, file_name){
  summary_stats <- describe(df)
  plot(summary_stats)
  print(summary_stats)
  df_long <- melt(is.na(df))
  p <- ggplot(df_long, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "white")) +
    labs(x = "Columns", y = "Rows", fill = "Missing") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  report = create_report(df, output_dir = output_loc, output_file = file_name)
}

venn_diagram <- function(data1, data2, column1, column2, output_plot1, output_plot2){
  # Calculate the intersection of values in the corresponding columns
  intersection <- intersect(data1, data2)
  # Create Venn diagram
  diagram <- venn.diagram(
    x = list(column1 = data1, column2 = data2),
    category.names = c(column1, column2),
    main = "Venn Diagram of Overlapping Rows",
    sub = paste("Number of overlapping rows:", length(intersection)),
    filename = output_plot1,  # Set filename to NULL to prevent saving to a file
    output = FALSE
  )
  
  
  # Find overlapping and non-overlapping values
  overlap <- intersect(data1, data2)
  only_in_data1 <- setdiff(data1, data2)
  only_in_data2 <- setdiff(data2, data1)
  
  # Create a table of overlapping and non-overlapping values
  overlap_df <- data.frame(Value = overlap, In_both = rep("Yes", length(overlap)))
  only_in_data1_df <- data.frame(Value = only_in_data1, In_both = rep("No (only in Data 1)", length(only_in_data1)))
  only_in_data2_df <- data.frame(Value = only_in_data2, In_both = rep("No (only in Data 2)", length(only_in_data2)))
  combined_df <- rbind(overlap_df, only_in_data1_df, only_in_data2_df)
  combined_df <- combined_df[order(combined_df$Value), ]
  print(combined_df)
}

# Creating box plot function ----
box_plots <- function(df, col1, col2, fill = NULL, col3 = NaN, col4 = NaN, filename = NaN) {
  df <- df[order(df[[fill]]), ] 
  df[[col1]] <- factor(df[[col1]], levels = unique(df[[col1]]))
  
  # Create a box plot with colors based on combinations of categories
  if (is.null(fill)) {
    p <- ggplot(df, aes_string(x = col1, y = col2)) +
      geom_boxplot() +
      labs(title = paste("Box Plot of ", col2, " by ", col1), x = col1, y = col2) +
      scale_x_discrete(breaks = unique(df[[col1]])[seq(1, length(unique(df[[col1]])), length.out = 5)]) +
      expand_limits(y = c(min(df[[col2]]), max(df[[col2]]))) +
      theme(axis.text = element_text(size = 14))
  } else {
    p <- ggplot(df, aes_string(x = col1, y = col2, fill = fill)) +
      geom_boxplot() +
      labs(title = paste("Box Plot of ", col2, " by ", col1), x = col1, y = col2) +
      scale_x_discrete(breaks = unique(df[[col1]])[seq(1, length(unique(df[[col1]])), length.out = 5)]) +
      expand_limits(y = c(min(df[[col2]]), max(df[[col2]]))) +
      theme(axis.text = element_text(size = 14))
  }
  
  df_max <- df[which(df[[col2]] == tapply(df[[col2]], df[[col1]], max)[df[[col1]]]), ]
  
  if (!is.nan(col3) & !is.nan(col4)) {
    labels <- ifelse(!is.nan(df_max[[col4]]), paste0(df_max[[col3]], " ", df_max[[col4]]), as.character(df_max[[col3]]))
    p <- p + geom_text(data = df_max, aes(label = labels), vjust = -1)
  } else if (!is.nan(col3)) {
    p <- p + geom_text(data = df_max, aes_string(label = col3), vjust = -1)
  } else if (!is.nan(col4)) {
    p <- p + geom_text(data = df_max, aes_string(label = col4), vjust = -1)
  }
  
  if (!is.nan(filename)) {
    ggsave(filename, plot = p + theme(axis.text = element_text(size = 24),   # Adjust axis text size
                                      axis.title = element_text(size = 24),  # Adjust axis title size
                                      plot.title = element_text(size = 30)), width = 60, height = 50, units = "cm", dpi = 300)
  }
  
  return(p)
}


violin_plots <- function(df, col1, col2, fill = NULL, col3 = NaN, col4 = NaN, filename = NaN){
  df <- df[order(df[[fill]]), ]
  df[[col1]] <- factor(df[[col1]], levels = unique(df[[col1]]))
  
  # Create a violin plot with colors based on combinations of categories
  if (is.null(fill)) {
    p <- ggplot(df, aes_string(x = col1, y = col2)) +
      geom_violin() +
      geom_jitter(width=0.15, alpha=0.5) +
      labs(title = paste("Box Plot of ", col2, " by ", col1), x = col1, y = col2) +
      scale_x_discrete(breaks = unique(df[[col1]])[seq(1, length(unique(df[[col1]])), length.out = 5)]) +
      expand_limits(y = c(min(df[[col2]]), max(df[[col2]]))) 
  } else {
    p <- ggplot(df, aes_string(x = col1, y = col2, fill = fill)) +
      geom_violin() +
      geom_jitter(width=0.15, alpha=0.5) +
      labs(title = paste("Box Plot of ", col2, " by ", col1), x = col1, y = col2) +
      scale_x_discrete(breaks = unique(df[[col1]])[seq(1, length(unique(df[[col1]])), length.out = 5)]) +
      expand_limits(y = c(min(df[[col2]]), max(df[[col2]]))) 
  }
  
  df_max <- df[which(df[[col2]] == tapply(df[[col2]], df[[col1]], max)[df[[col1]]]), ]
  
  if (!is.nan(col3) & !is.nan(col4)) {
    labels <- ifelse(!is.nan(df_max[[col4]]), paste0(df_max[[col3]], " ", df_max[[col4]]), as.character(df_max[[col3]]))
    p <- p + geom_text(data = df_max, aes(label = labels), vjust = -1)
  } else if (!is.nan(col3)) {
    p <- p + geom_text(data = df_max, aes_string(label = col3), vjust = -1)
  } else if (!is.nan(col4)) {
    p <- p + geom_text(data = df_max, aes_string(label = col4), vjust = -1)
  }
  
  if (!is.nan(filename)) {
    ggsave(filename, plot = p+ theme(axis.text = element_text(size = 24),   # Adjust axis text size
                                     axis.title = element_text(size = 24),  # Adjust axis title size
                                     plot.title = element_text(size = 30)), width = 60, height = 60, units = "cm")
  }
  
  return(p)
}


dot_plots <- function(df, col1, col2, groups = NULL,thresholds=NULL, filename = NULL) {
  # Create a box plot with colors based on combinations of categories
  if (is.null(fill)) {
    p <- ggplot(df, aes_string(x = col1, y = col2, shape = groups)) +
      geom_point() +
      labs(title = paste("Dot Plot of ", col2, " by ", col1), x = col1, y = col2) +
      expand_limits(y = c(min(df[[col2]]), max(df[[col2]]))) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  } else {
    p <- ggplot(df, aes_string(x = col1, y = col2, shape = groups)) +
      geom_point() +
      labs(title = paste("Dot Plot of ", col2, " by ", col1), x = col1, y = col2) +
      expand_limits(y = c(min(df[[col2]]), max(df[[col2]]))) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  
  # Add threshold lines and assign colors
  if (!is.null(thresholds)) {
    num_thresholds <- length(thresholds)
    color_palette <- colorRampPalette(c("red", "blue"))(num_thresholds)
    
    for (i in 1:num_thresholds) {
      p <- p + geom_hline(yintercept = thresholds[i], linetype = "dashed", color = color_palette[i])
    }
  }
  
  if (!is.nan(filename)) {
    ggsave(filename, plot = p, width = 20, height = 20, units = "cm")
  }
  
  return(p)
}

#Pie plot ----
pie_plots <- function(df, group_var, filename) {
  category_counts <- df %>%
    count({{ group_var }}) %>%
    mutate(percent = freq / sum(freq) * 100)
  
  print(category_counts)
  
  pie_chart <- ggplot(category_counts, aes(x="", y = freq, fill = !!as.name(names(category_counts)[1]))) +
    geom_col() +
    coord_polar(theta = "y", start=0) + 
    theme_void() +
    labs(title = paste("Category Counts for ", group_var)) +
    geom_text(aes(label = paste0(round(percent, 1), '%')), position = position_stack(vjust = 0.5))+
    scale_fill_brewer(palette = 'Spectral')
  
  print(pie_chart)
  ggsave(filename, plot=pie_chart)
  
}

make_bold_names <- function(mat, rc_fun, row=TRUE) {
  bold_names <- rc_fun(mat)
  ids <- rownames(mat) %>% match(rc_fun(mat))
  if(row==FALSE){
    ids <- colnames(mat) %>% match(rc_fun(mat))
  }
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

#Heatmap plot function ----
heatmap_plots <- function(df, a_col = NA, a_row = NA, filename, fontsize = 18, fontsize_row = TRUE, a_color = NULL, title = 'Heatmap', c_row = TRUE, c_col = TRUE, height=150, width = 100) {
  if (fontsize_row == FALSE){
    fontsize_r = 7
  }else if (fontsize_row == TRUE){
    fontsize_r = nrow(df) / 3
  }
    # Create heatmap
    h <- pheatmap(
      df,
      annotation_row = a_row,
      annotation_col = a_col,
      #distfun = 'euclidean',
      #color = hcl.colors(50, "BluYl"),
      #color = colorRampPalette(brewer.pal(3, "RdBu"))(256), # Define color range
      na_color = "grey",  # Set color for NA values
      annotation_colors = a_color,
      labels_row = make_bold_names(df, rownames),
      labels_col = make_bold_names(df, colnames, FALSE),
      main = title,
      fontsize_row=18,
      fontsize_col=8,
      fontsize = fontsize,
      cutree_rows = 2,
      cutree_cols = 2, 
      cluster_rows=c_row, cluster_cols=c_col
      #clustering_distance_rows = "euclidean",
      #clustering_distance_cols = "euclidean"
    )
    if (!is.nan(filename)) {
      ggsave(filename, plot = h, width = width, height = nrow(df), units = "cm", limitsize = FALSE)
    }
    return(h)
}

heatmap_plots <- function(df, 
                          a_col = NA, 
                          a_row = NA, 
                          filename, 
                          fontsize = 24,          # Increased base font size
                          fontsize_row = TRUE, 
                          a_color = NULL, 
                          c_method = 'average', 
                          title = "Heatmap", height=150, width = 100) {
  
  # Dynamically calculate row font size
  if (fontsize_row == FALSE){
    fontsize_r = 10                # Set a minimum font size for rows if not scaling
  } else {
    fontsize_r = min(nrow(df) / 2, 12) # Dynamically adjust font size but limit max size
  }
  
  # Create the heatmap
  h <- pheatmap(
    df,
    annotation_row = a_row,
    annotation_col = a_col,
    distfun = function(x) dist(x, method = "euclidean"),
    na_color = "grey",  # Set color for NA values
    annotation_colors = a_color,
    labels_row = make_bold_names(df, rownames),
    labels_col = make_bold_names(df, colnames, FALSE),
    fontsize_row = 16,      # Increased row label font size
    fontsize_col = 16,              # Increased column label font size
    fontsize = fontsize,            # Base font size for other elements, including legend
    cellwidth = 15,                 # Adjust cell width for visibility
    cellheight = 15,                # Adjust cell height for visibility
    cutree_rows = 3,
    cutree_cols = 3,
    clustering_method = c_method, 
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    main = title,                   # Set title for the heatmap
    treeheight_row = 100,           # Increase row dendrogram height
    legend = TRUE,                  # Ensure legend is shown
    legend_labels = "Custom Legend Title"  # Correct way to set the legend title
    #legend_labels_gp = gpar(fontsize = 16), # Legend labels font size
    #legend_width = unit(4, "cm"),   # Increase the width of the legend (color bar)
    #legend_height = unit(8, "cm")   # Increase the height of the color bar
  )
  
  # Now, increase the overall plot size for a larger output
  if (!is.nan(filename)) {
    ggsave(filename, plot = h, width = width, height = height, units = "cm", limitsize = FALSE)  # Increased output size
  }
  
  return(h)
}

# Function to create a ComplexHeatmap
complexheatmap_plots <- function(df, 
                          a_col = NA, 
                          a_row = NA, 
                          filename = NA, 
                          fontsize = 24,  # Base font size for heatmap
                          fontsize_row = TRUE, 
                          a_color = NULL, 
                          c_method = 'average', 
                          title = "Heatmap", 
                          height = 150, 
                          width = 100) {
  library(ComplexHeatmap)
  library(circlize)  # Needed for colorRamp2
  # Dynamically calculate row font size
  if (fontsize_row == FALSE) {
    fontsize_r = 15  # Set minimum font size
  } else {
    fontsize_r = min(nrow(df) / 2, 12)  # Adjust font size with a max limit
  }
  
  print(a_col)
  ha_col <- HeatmapAnnotation(df = a_col, col = a_color)
  
  ha_row <- if (!is.null(a_row) | !is.na(a_row)) {
    rowAnnotation(df = a_row, col = a_color)
  } else {
    NULL
  }
  
  # Create the heatmap
  h <- ComplexHeatmap::Heatmap(
    df,
    name = "DSS2",  # Legend header
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_method_rows = c_method,
    clustering_method_columns = c_method,
    row_names_gp = gpar(fontsize = fontsize_r, fontface = "bold"),  # Row font size
    column_names_gp = gpar(fontsize = fontsize_r, fontface = "bold"),  # Column font size
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    top_annotation = NULL,  # Column annotations
    left_annotation = NULL,  # Row annotations
    na_col = "grey",
    row_gap = unit(1, "cm"),             # Space between rows
    column_gap = unit(1, "cm"),          # Space between columns
    heatmap_legend_param = list(
      title = "DSS2 \n ",  # Add header to the legend
      title_gp = gpar(fontsize = 12, fontface = "bold"),  # Legend title font
      labels_gp = gpar(fontsize = fontsize - 6)  # Legend label font size
    ),
    width = unit(width, "cm"),  # Adjust heatmap width
    height = unit(height, "cm"),  # Adjust heatmap height
    column_title = title,
    column_title_gp = gpar(fontsize = fontsize, fontface = "bold"), 
    row_dend_width = unit(8, "cm"),
    clustering_distance_rows = "euclidean",
    column_dend_height = unit(8, "cm"),
    clustering_distance_columns = "euclidean"
    
  )
  
  # Save heatmap if filename is specified
  if (!is.na(filename)) {
    pdf(file = filename, width = width, height = height)  # PDF output
    plot(h)  # Use ComplexHeatmap's draw method
    dev.off()  # Close the graphic device
  } else {
    # Display the heatmap interactively
    plot(h)
  }
  
  return(h)
}

#pca plots function ----
pca_plots <- function(df, title = "PCA Plot", file_loc, label, color){
  print(df[[label]])
  # Perform PCA excluding specific columns
  pca_result <- prcomp(df[, !names(df) %in% c(label, color)], scale. = TRUE)
  
  # Extract PC scores
  pc_scores <- as.data.frame(pca_result$x[, 1:2])  # Using only the first two principal components
  
  # Combine PC scores with color information
  pca_data <- cbind(pc_scores, lab = df[label])
  print(pca_data)
  pca_data <- cbind(pca_data, patient_id = df[color])
  
  # Plot PCA
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = get(color), label = get(label))) +
    geom_point() +
    geom_text(nudge_x = 0.2, nudge_y = 0.2, size = 3) +
    labs(color = "Labs", title = title) + theme(legend.text = element_text(size = 12),       # Increase legend label size
            legend.title = element_text(size = 15), 
            plot.title = element_text(size = 20)) 
  ggsave(plot = p, filename = file_loc)
  return(p)
}

xicor <- function(X, Y, ties = TRUE){
  n <- length(X)
  r <- rank(Y[order(X)], ties.method = "random")
  set.seed(42)
  if(ties){
    l <- rank(Y[order(X)], ties.method = "max")
    return( 1 - n*sum( abs(r[-1] - r[-n]) ) / (2*sum(l*(n - l))) )
  } else {
    return( 1 - 3 * sum( abs(r[-1] - r[-n]) ) / (n^2 - 1) )    
  }
}

# Function to create correlation matrix using xicor function
create_correlation_matrix <- function(data, ties = TRUE) {
  n_vars <- ncol(data)
  cor_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars, dimnames = list(colnames(data), colnames(data)))
  
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      cor_matrix[i, j] <- xicor(data[, i], data[, j], ties = ties)
    }
  }
  
  return(cor_matrix)
}

correlation_plot <- function(df, method = 'pearson', title = 'Correlation matrix plot', filesave = "correlation_matrix.pdf", width = 80, height = 80, label_size = 0.15, bar_size = 5){
  library(corrplot)
  if (method == 'xicor'){
    cor_df <- create_correlation_matrix(df)
  }
  else{
    cor_df <- cor(df, method=method)
  }
  
  # Open a PDF device for saving with specified width and height
  pdf(filesave, width = width, height = height)
  par(cex.lab = label_size, cex = bar_size)
  # Create the correlation plot
  p <- corrplot(cor_df, method = 'color', order = 'alphabet', sig.level = 0.05, mar=c(0,0,1,0),
                tl.col='black', tl.srt = 60, bg = "grey", title = title)
  
  # Close the PDF device
  dev.off()
  
  return(p)
}

grid_plots <- function(df, col1, col2, label,title1, title2, filename){
  library(gridExtra)
  require(qqplotr)
  library(ggplot2)
  library(grid)
  library("psychTools")
  # Calculate R2 value
  lm_model <- lm(get(col2) ~ get(col1), data = df)
  r2 <- summary(lm_model)$r.squared
  p_value <- summary(lm_model)$coefficients[2,4]
  
  # Calculate xicor value
  xicor_value <- cor(df[[col1]], df[[col2]], method = 'spearman')
  # Line chart
  p3 <- ggplot(df, aes(x = get(col1), y = get(col2))) +
    geom_line(color = "gray20") +
    xlab(paste0('', col1)) +
    ylab(paste0('', col2))
  
  # Scatter plot with R2 value
  p4 <- ggplot(df, aes(x = get(col1), y = get(col2), label=get(label))) +
    geom_point(color = "grey20") +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    xlab(paste0('', col1)) +
    ylab(paste0('', col2)) +
    geom_text() +
    geom_label(aes(x = Inf, y = Inf, label = paste("R² = ", round(r2, 2), "\nP-value = ", round(p_value,3),"\nXicor = ", round(xicor_value, 2))),
               hjust = 1.1, vjust = 2, size = 15, color = "blue", fontface = "bold",
               fill = "white", label.size = NA) 
  
  gg1 <- ggplot(data = df, mapping = aes(sample = get(col1))) +
    geom_qq_band(bandType = "ks", mapping = aes(fill = "KS"), alpha = 0.5) +
    geom_qq_band(bandType = "ts", mapping = aes(fill = "TS"), alpha = 0.5) +
    geom_qq_band(bandType = "pointwise", mapping = aes(fill = "Normal"), alpha = 0.5) +
    geom_qq_band(bandType = "boot", mapping = aes(fill = "Bootstrap"), alpha = 0.5) +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    scale_fill_discrete("Bandtype") +
    ggtitle(title1)
  
  gg2 <- ggplot(data = df, mapping = aes(sample = get(col2))) +
    geom_qq_band(bandType = "ks", mapping = aes(fill = "KS"), alpha = 0.5) +
    geom_qq_band(bandType = "ts", mapping = aes(fill = "TS"), alpha = 0.5) +
    geom_qq_band(bandType = "pointwise", mapping = aes(fill = "Normal"), alpha = 0.5) +
    geom_qq_band(bandType = "boot", mapping = aes(fill = "Bootstrap"), alpha = 0.5) +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    scale_fill_discrete("Bandtype") +
    ggtitle(title2)
  
  g <- arrangeGrob(p3, p4, gg1, gg2) 
  ggsave(file=filename, g, units="cm", width = 100, height = 100)
  return(g)
}

wilcox_test <- function(df, comparison_group, testing_group, value){
  # Create an empty list to store results
  results <- list()
  result_df <- data.frame(
    comparison_group = character(),
    p_value = numeric(),
    p_adj_value = numeric(),
    logFC = numeric(), 
    FC = numeric
  )
  skipped_comp_group = c()
  # Loop through comparison groups
  for (comp_group in unique(df[[comparison_group]])) {
    # Subset dataframe for each comparison group
    subset_df <- subset(df, df[[comparison_group]] == comp_group)
    print(comp_group)
    # Subset data for testing_group = 1 and 0
    l1 <- subset(subset_df, subset_df[[testing_group]] == 1, select = value)
    l1 <- l1[order(l1[[value]]),]
    print('l1')
    print(l1)
    l0 <- subset(subset_df, subset_df[[testing_group]] == 0, select = value)
    l0 <- l0[order(l0[[value]]),]
    print('l0')
    print(l0)
    if (dim(l1)[1] == 0 | dim(l0)[1] == 0){
      skipped_comp_group <- comp_group
      next
    }
    
    # Perform Wilcoxon rank sum test
    test_result <- wilcox.test(l1[[value]], l0[[value]])
    logFC <- log(median(l1[[value]]) / median(l0[[value]]))
    print(median(l1[[value]]))
    FC <- median(l1[[value]]) - median(l0[[value]])
    # Store the result
    results[[as.character(comp_group)]] <- test_result
    p_value <- results[[comp_group]]$p.value
    p_adj_value <- NaN
    result_df <- rbind(result_df, data.frame(comparison_group = comp_group, p_value = p_value, p_adj_value = p_adj_value, logFC = logFC, FC=FC))
  }
  
  result_df$p_adj_value <- p.adjust(result_df$p_value, "fdr")
  print(paste0("Groups skipped because of not enough observations ", skipped_comp_group))
  return(result_df)
}

#Create Volcano plots
create_volcano_plot <- function(result_df, title, filename, threshold_pvalue = 0.05, threshold_logFC = 1, x_lab = "Change in DSS", x_axis = "FC", y_axis = "p_adj_value", y_lab = "-log10(p-adj-value)") {
  library(ggplot2)
  library(ggrepel)
  
  # Define a new column for color based on significance
  result_df$point_color <- ifelse(result_df[[y_axis]] < threshold_pvalue, "significant", "non-significant")
  
  # Define a new column for label colors
  result_df$label_color <- ifelse(grepl('clax', result_df$comparison_group), "blue", 
                                  ifelse(grepl('ib', result_df$comparison_group), "red", "black"))
  
  # Create the plot
  volcano_plot <- ggplot(result_df, aes(x = get(x_axis), y = -log10(get(y_axis)), color = point_color)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("significant" = "#FFD966", "non-significant" = "black")) +
    geom_hline(yintercept = -log10(threshold_pvalue), linetype = "dashed", color = "darkgrey", size = 1) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 1) +
    labs(x = x_lab, y = y_lab, title = title) +
    theme(legend.position = "none")
  
  # Add labels for significant points with specified colors
  volcano_plot <- volcano_plot +
    geom_text_repel(data = subset(result_df, point_color == "significant"),
                    aes(label = comparison_group), color = result_df$label_color[result_df$point_color == "significant"], 
                    hjust = -0.1, vjust = -0.5, size = 3, max.overlaps = Inf)
  
  ggsave(filename, plot = volcano_plot, width = 20, height = 10)
  return(volcano_plot)
}
