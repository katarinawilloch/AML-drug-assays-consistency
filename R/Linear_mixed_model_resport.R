library(lme4)
library(broom.mixed)
library(ggplot2)
library(sjPlot)
library(xtable)
library(htmltools)
library(forcats)


init_plots <- function(df, pred, group){
  p1 <- ggplot(df, aes(x = get(group), y = get(pred))) +
    geom_point() +
    geom_smooth(method = "lm", se = F, color = "red", linetype = "dashed") +
    theme_bw() +
    labs(y = "Frequency\n(Prepositions)")
  p2 <- df %>%
    dplyr::group_by(get(group)) %>%
    dplyr::mutate(Mean = round(mean(get(pred)), 1), SD = round(sd(get(pred)), 1)) %>%
    ggplot(aes(get(group), get(pred))) + 
    geom_boxplot() +
    geom_text(aes(label = paste("M = ", Mean, sep = ""), y = 1)) +
    geom_text(aes(label = paste("SD = ", SD, sep = ""), y = 0)) +
    theme_bw(base_size = 15) +
    labs(x = group) +                      
    labs(y = pred, cex = .75) 
  p3 <- ggplot(df, aes(get(pred))) +
    geom_histogram() +
    theme_bw() + 
    labs(y = "Count", x = "Frequency (Prepositions)")
  grid.arrange(grobs = list(p1, p2, p3), widths = c(1, 1), layout_matrix = rbind(c(1, 1), c(2, 3)))
}

fit_model <- function(df, formula){
  mixed_model <- lmer(formula, data = df)
  return(mixed_model)
}

model_interpretation <- function(df, pred, mem_model, output_file, header_1 = "Linear Mixed Model Interpretation", adj_nr = 1){
  # Summary of the model
  model_summary <- summary(mem_model)
  tidy_model <- tidy(mem_model, effects = "fixed")
  tidy_model$corrected_p_value <- tidy_model$p.value*adj_nr
  tidy_model$nr_models_adj_for <- adj_nr
  
  # Convert the model summary to an HTML table
  model_summary_df <- as.data.frame(coef(summary(mem_model)))
  print(model_summary_df)
  p_value_list <- model_summary$coefficients[, "Pr(>|t|)"]
  model_summary_df$corrected_p_value <- p_value_list*adj_nr
  model_summary_df$nr_models_adj_for <- adj_nr
  #model_summary_df$corrected_p_value 
  model_summary_html <- print(xtable(model_summary_df), type = "html")
  
  # Convert the fixed effects to an HTML table
  fixed_effects_html <- print(xtable(tidy_model), type = "html")
  
  # Extract random effects and convert to data frame
  random_effects_df <- as.data.frame(ranef(mem_model))
  
  # Convert the random effects to an HTML table
  random_effects_html <- print(xtable(random_effects_df), type = "html")
  
  # Function to save a check_model plot
  save_check_model_plot <- function(check_model_output, file_path) {
    tryCatch({
      # Render the check_model object
      plot_output <- plot(check_model_output)
      
      # Save the plot using the appropriate method
      if (inherits(plot_output, "ggplot")) {
        ggsave(file_path, plot = plot_output, width = 20, height = 15)
      } else if (inherits(plot_output, "html")) {
        writeLines(plot_output, file_path)
      } else {
        stop("Unsupported plot type")
      }
      
      paste0('<img src="', file_path, '" style="width: 700px; height: 500px;">')
    }, error = function(e) {
      paste0('<p>Error saving plot: ', e$message, '</p>')
    })
  }
  
  # Example usage
  check_model_output <- check_model(mem_model)
  check_model_plot_path <- tempfile(fileext = ".png")
  html_check_model_plot <- save_check_model_plot(check_model_output, check_model_plot_path)
  
  
  
  # Function to create diagnostic plots
  create_diagnostic_plots <- function(model) {
    # Extract residuals and fitted values
    residuals <- residuals(model)
    fitted <- fitted(model)
    
    p3 <- ggplot(data = data.frame(residuals = residuals), aes(x = residuals)) +
      geom_histogram(bins = 30, fill = "blue", color = "black") +
      ggtitle("Histogram of Residuals") +
      xlab("Residuals") +
      ylab("Frequency")
    
    list(p3)
  }
  
  # Function to create additional plots and analyses
  create_additional_plots <- function(df, pred, model) {
    # ICC Calculation
    icc_value <- performance::icc(model)
    icc_html <- paste0('<p>Intraclass Correlation Coefficient (ICC): ', round(icc_value, 4), '</p>')
    
    # REsim
    sim_results <- REsim(model)
    
    # Plot REsim
    re_sim_plot <- plotREsim(sim_results)
    
    # Predicted vs Actual Plot
    predicted_values <- fitted(model)
    actual_values <- df[[pred]]
    predicted_vs_actual <- tryCatch({ggplot(data = data.frame(predicted = predicted_values, actual = actual_values), aes(x = predicted, y = actual)) +
        geom_point() +
        geom_smooth(method = "loess") +
        ggtitle("Predicted vs Actual Values") +
        xlab("Predicted values") +
        ylab("Actual values")}, error = function(e){
          ggplot() + ggtitle(paste("Error in dotplot:", e$message))
        })
    
    # Dotplot of Random Effects
    dotplot_re <- tryCatch({
      plot <- dotplot(ranef(model, condVar = TRUE), scales = list(x = list(relation = "free")))
      cat("Class of dotplot_re: ", class(plot), "\n")
      plot
    }, error = function(e) {
      ggplot() + ggtitle(paste("Error in dotplot:", e$message))
    })
    
    # Q-Q plot of Random Effects
    qq_plot_re <- tryCatch({
      # Extract random effects
      random_effects <- ranef(model)
      
      # Flatten the random effects to a data frame
      random_effects_df <- as.data.frame(random_effects)
      
      # Ensure random effects column is present
      if (ncol(random_effects_df) == 0) {
        stop("No random effects data to plot.")
      }
      
      # Create Q-Q plot
      qq_plot <- ggplot(random_effects_df, aes(sample = .[[1]])) +  # Adjust column index as necessary
        stat_qq() +
        stat_qq_line() +
        ggtitle("Q-Q Plot of Random Effects") +
        xlab("Theoretical Quantiles") +
        ylab("Sample Quantiles")
      
      qq_plot
    }, error = function(e) {
      ggplot() + ggtitle(paste("Error in Q-Q plot:", e$message))
    })
    
    
    model_plot <-  tryCatch({
      plot <- visualize(model)
      if (inherits(plot, "gg")) {
        plot <- plot + ggtitle("Model Visualization") +
          xlab("Variables") + ylab("Values")
      }
      plot
    }, error = function(e) {
      ggplot() + ggtitle(paste("Error in visualize:", e$message))
    })
    
    list(
      icc = icc_html,
      re_sim_plot = re_sim_plot,
      predicted_vs_actual = predicted_vs_actual,
      dotplot_re = dotplot_re,
      qq_plot_re = qq_plot_re, 
      model_plot = model_plot
    )
  }
  
  # Create diagnostic plots
  diagnostic_plots <- create_diagnostic_plots(mem_model)
  
  # Create additional plots and analyses
  additional_plots <- create_additional_plots(df, pred, mem_model)
  
  save_html_check_model <- function(plot, file_path) {
    tryCatch({
      html <- as.character(plot)  # Convert to HTML if possible
      writeLines(html, file_path)
      paste0('<iframe src="', file_path, '" style="width: 700px; height: 500px;"></iframe>')
    }, error = function(e) {
      paste0('<p>Error saving HTML plot: ', e$message, '</p>')
    })
  }
  
  save_plot_to_html <- function(plot) {
    tryCatch({
      cat("Class of plot: ", class(plot), "\n")  # Debugging line
      if (inherits(plot, "ggplot")) {
        plot_path <- tempfile(fileext = ".png")
        ggsave(plot_path, plot = plot, width = 7, height = 5)
        paste0('<img src="', plot_path, '" style="width: 700px; height: 500px;">')
      } else if (inherits(plot, "trellis")) {
        save_lattice_plot_to_html(plot)
      } else if (inherits(plot, "check_model")) {
        plot_path <- tempfile(fileext = ".png")
        ggsave(plot_path, plot = plot, width = 7, height = 5)
        paste0('<img src="', plot_path, '" style="width: 700px; height: 500px;">')
      } else if (is.list(plot)) {
        plot_paths <- lapply(plot, function(p) {
          if (inherits(p, "trellis")) {
            save_lattice_plot_to_html(p)
          } else if (inherits(p, "ggplot")) {
            plot_path <- tempfile(fileext = ".png")
            ggsave(plot_path, plot = p, width = 7, height = 5)
            paste0('<img src="', plot_path, '" style="width: 700px; height: 500px;">')
          } else {
            paste0('<p>Unable to handle plot of class: ', class(p), '</p>')
          }
        })
        paste(plot_paths, collapse = "\n")
      } else if (inherits(plot, "icc") || inherits(plot, "data.frame")) {
        paste0('<p>', as.character(plot), '</p>')
      } else if (is.character(plot)) {
        plot
      } else {
        paste0('<p>', as.character(plot), '</p>')
      }
    }, error = function(e) {
      paste0('<p>Error saving plot: ', e$message, '</p>')
    })
  }
  
  # Convert all plots to HTML
  diagnostic_plots_html <- lapply(diagnostic_plots, save_plot_to_html)
  additional_plots_html <- lapply(additional_plots, save_plot_to_html)
  
  
  # Combine all plots into a single HTML string
  diagnostic_plots_html <- paste(diagnostic_plots_html, collapse = "\n")
  additional_plots_html <- paste(unlist(additional_plots_html), collapse = "\n")
  
  # Create the HTML report
  report_html <- htmltools::tagList(
    htmltools::tags$html(
      htmltools::tags$head(
        htmltools::tags$title("Linear Mixed Model Interpretation")
      ),
      htmltools::tags$body(
        htmltools::tags$h1(header_1),
        htmltools::tags$h2("Model Summary"),
        htmltools::HTML(model_summary_html),
        htmltools::tags$h2("Fixed Effects"),
        htmltools::HTML(fixed_effects_html),
        htmltools::tags$h2("Random Effects"),
        htmltools::HTML(random_effects_html),
        htmltools::tags$h2("Additional Plots and Analyses"),
        htmltools::HTML(additional_plots_html),
        htmltools::tags$h2("Check Model Output"),
        htmltools::HTML(html_check_model_plot),
        htmltools::tags$h2("Model Diagnostics"),
        htmltools::HTML(diagnostic_plots_html)
      )
    )
  )
  
  # Save the HTML report to a file
  report_path <- output_file
  htmltools::save_html(report_html, report_path)
  
  # Print the location of the report
  cat("HTML report saved to:", report_path, "\n")
}
