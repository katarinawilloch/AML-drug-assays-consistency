library(lme4)
library(lmerTest)
library(flexplot)
library(RColorBrewer)
library(plotly)
library(performance)


#Functions
#----Model----
create_model <- function(formula, df, model_spec="lmer", random_effects){
  
  if(model_spec == "lmer"){
    model <- lmer(formula, data = df)
  }
  else if(model_spec == "nlme"){
    model <- lme(formula, random = random_effects, data=df)
  }
  else if(model_spec == "lm"){
    model <- lm(formula, data=df)
  }
  else if(model_spec == "binomial"){
    model <- glm(formula, data=df, family = "binomial")
  }
  else{
    "No model specification given"
  }
  
  return(c(model, model_spec))
}

#----Model results----
model_results <- function(models, group = "Response_metric", data_frame_org, scaling_group = "drug"){
  df <- data.frame(
    model = character(),
    Fixed_Effect_Term = character(),
    Random_Eeffect_Term = character(),
    Response_metric = character(),
    Fixed_Effect_Coefficient = numeric(),
    p_value = numeric(),
    Intercept = numeric(),
    ref_group = character(),
    Degree_of_Freedom = numeric(),
    stringsAsFactors = FALSE, 
    lower = numeric(),
    upper = numeric()
  )
  # confints <- pbmclapply(models, function(mod) {
  #   print(model)
  #   model_formula <- formula(model)
  #   dependent_var <- as.character(model_formula[[2]])
  #   fixed_effects <- fixef(model)  
  #   fixed_effect_terms <- names(fixed_effects)
  #   fixed_effect_names <- fixed_effect_terms[fixed_effect_terms != "(Intercept)"]
  #   fixed_effect_coef <- fixed_effects[fixed_effect_names]
  #   coef_model <- summary(model)$coefficients
  #   p_values <- coef_model[, "Pr(>|t|)"]
  #   p_values_excluding_intercept <- p_values[fixed_effect_terms != "(Intercept)"]
  #   intercept <- fixed_effects[fixed_effect_terms[fixed_effect_terms == "(Intercept)"]]
  #   degree_freedom <- coef_model[, "df"]
  #   degree_freedom <- degree_freedom[fixed_effect_terms != "(Intercept)"]
  #   print(degree_freedom)
  #   response_metric <- dependent_var
  #   c <- confint(model, method = "Wald")
  #   ci_fixed <- c[rownames(c) %in% fixed_effect_names, , drop = FALSE]  # Subset only fixed effects
  #   lower <- ci_fixed[, 1]
  #   upper <- ci_fixed[, 2]
  #   
  #   
  #   #Create a temporary data frame for this model's fixed effects
  #   temp_df <- data.frame(
  #     model = paste0(model_formula[2], model_formula[1], model_formula[3]),
  #     Fixed_Effect_Term = fixed_effect_names,
  #     Random_Effect_Term = str_remove(strsplit(as.character(model_formula[3]), split = "\\+ \\(")[[1]][2], "\\)"),
  #     Response_metric = response_metric,
  #     Fixed_Effect_Coefficient = fixed_effect_coef,
  #     p_value = p_values_excluding_intercept,
  #     Intercept = intercept,
  #     Degree_of_Freedom = degree_freedom,
  #     stringsAsFactors = FALSE, 
  #     lower = lower, 
  #     upper = upper
  #   )
  #   
  #   # n_per_group <- tapply(data_frame_org[,response_metric], data_frame_org[,scaling_group], length)
  #   # print(n_per_group)
  #   # 
  #   # means <- tapply(data_frame_org[,response_metric], data_frame_org[,scaling_group], mean)
  #   # sds <- tapply(data_frame_org[,response_metric], data_frame_org[,scaling_group], sd)
  #   # 
  #   # overall_mean <- sum(means * n_per_group) / sum(n_per_group)
  #   # overall_sd <- sqrt(sum((sds^2 + means^2) * n_per_group) / sum(n_per_group) - overall_mean^2)
  #   # 
  #   # print(fixed_effect_names)
  #   # fixed_effect_coeff <- fixef(model)[fixed_effect_names]  
  #   # print(fixed_effect_coeff[fixed_effect_names])
  #   # print('fixed_effect_coef')
  #   # unscaled_coeff <- fixed_effect_coeff * overall_sd + overall_mean
  #   # print(overall_mean)
  #   # print(unscaled_coeff)
  #   # untransformed_coef <-  (0.5 * unscaled_coeff + 1)^2
  #   # print(untransformed_coef)
  #   df <- rbind(df, temp_df)
  # }, mc.cores = parallel::detectCores() - 1)
  #Extract fixed effects and their coefficients
  get_reference_levels <- function(model) {
    # Get model frame (data used in model)
    mf <- model.frame(model)
    
    # Get the fixed effect variable names (exclude intercept)
    fixed_terms <- attr(terms(model), "term.labels")
    
    # Initialize list for reference levels
    ref_levels <- list()
    
    for (term in fixed_terms) {
      # Handle interaction terms separately
      if (grepl(":", term)) next
      
      # Try to get variable from model frame
      var <- mf[[term]]
      
      # Only care about character or factor vars (treated as categorical)
      if (is.character(var) || is.factor(var)) {
        ref_levels[[term]] <- levels(factor(var))[1]
      }
    }
    
    return(ref_levels[[1]])
  }
  
  
  for (model in models) {
    model_formula <- formula(model)
    dependent_var <- as.character(model_formula[[2]])
    fixed_effects <- fixef(model)
    fixed_effect_terms <- names(fixed_effects)
    fixed_effect_names <- fixed_effect_terms[fixed_effect_terms != "(Intercept)"]
    fixed_effect_coef <- fixed_effects[fixed_effect_names]
    coef_model <- summary(model)$coefficients
    p_values <- coef_model[, "Pr(>|t|)"]
    p_values_excluding_intercept <- p_values[fixed_effect_terms != "(Intercept)"]
    intercept <- fixed_effects[fixed_effect_terms[fixed_effect_terms == "(Intercept)"]]
    degree_freedom <- coef_model[, "df"]
    degree_freedom <- degree_freedom[fixed_effect_terms != "(Intercept)"]
    response_metric <- dependent_var
    c <- confint(model, method = "Wald")
    ci_fixed <- c[rownames(c) %in% fixed_effect_names, , drop = FALSE]  # Subset only fixed effects
    lower <- ci_fixed[, 1]
    upper <- ci_fixed[, 2]
    ref_group <- get_reference_levels(model)
    print(ref_group)

    #Create a temporary data frame for this model's fixed effects
    temp_df <- data.frame(
      model = paste0(model_formula[2], model_formula[1], model_formula[3]),
      Fixed_Effect_Term = fixed_effect_names,
      Random_Effect_Term = str_remove(strsplit(as.character(model_formula[3]), split = "\\+ \\(")[[1]][2], "\\)"),
      Response_metric = response_metric,
      Fixed_Effect_Coefficient = fixed_effect_coef,
      p_value = p_values_excluding_intercept,
      Intercept = intercept,
      ref_group = ref_group,
      Degree_of_Freedom = degree_freedom,
      stringsAsFactors = FALSE,
      lower = lower,
      upper = upper
    )

    # n_per_group <- tapply(data_frame_org[,response_metric], data_frame_org[,scaling_group], length)
    # print(n_per_group)
    #
    # means <- tapply(data_frame_org[,response_metric], data_frame_org[,scaling_group], mean)
    # sds <- tapply(data_frame_org[,response_metric], data_frame_org[,scaling_group], sd)
    #
    # overall_mean <- sum(means * n_per_group) / sum(n_per_group)
    # overall_sd <- sqrt(sum((sds^2 + means^2) * n_per_group) / sum(n_per_group) - overall_mean^2)
    #
    # print(fixed_effect_names)
    # fixed_effect_coeff <- fixef(model)[fixed_effect_names]
    # print(fixed_effect_coeff[fixed_effect_names])
    # print('fixed_effect_coef')
    # unscaled_coeff <- fixed_effect_coeff * overall_sd + overall_mean
    # print(overall_mean)
    # print(unscaled_coeff)
    # untransformed_coef <-  (0.5 * unscaled_coeff + 1)^2
    # print(untransformed_coef)
    df <- rbind(df, temp_df)
  }

  get_significance_stars <- function(p) {
    if (p < 0.001) {
      return("***")
    } else if (p < 0.01) {
      return("**")
    } else if (p < 0.05) {
      return("*")
    } else {
      return(" ")
    }
  }
  
  # Check if `df` has the necessary columns
  if (!all(c("Response_metric", "model") %in% colnames(df))) {
    stop("The columns 'Response_metric' and/or 'model' are missing in the data frame.")
  }
  
  
  # Check for missing values in critical columns
  if (any(is.na(df$Response_metric)) || any(is.na(df$model))) {
    warning("There are NA values in 'Response_metric' or 'model' columns.")
  }
  
  print("######################STARS#######################")
  #Add significance stars to the data frame
  df$significance_stars <- sapply(df$p_value, get_significance_stars)
  #Count the number of fixed effects per category
  fixed_effect_counts <- df %>%
    group_by(Response_metric, model) %>%
    filter(row_number() == 1) %>%   
    ungroup() %>%                   
    group_by(Response_metric) %>%           
    dplyr::summarize(count = n()) %>%
    as.data.frame()
  df <- df %>%
    left_join(fixed_effect_counts, by = "Response_metric")
  print("####################LEFT")
  df <- df %>% mutate(corrected_p_value = ifelse(df$p_value * df$count >= 1, 1, df$p_value * df$count))
  
  
  
  rownames(df) <- NULL
  df <- df %>% mutate(Fixed_Effect_Term = case_when(Fixed_Effect_Term ==  "positive_controlBzCl" ~ "Positive Control: BzCl + Doses: 5 \n+ Readout: CellTiter-Glo\n+ Cell Counting: Trypan blue", 
                                                    Fixed_Effect_Term ==  "time_until_sample_usageHandled within 2-72 hours" ~ "Time until sample usage >2h",
                                                    Fixed_Effect_Term ==  "time_until_sample_usageHandled within 24 hours" ~ "Time until sample usage < 24h",
                                                    Fixed_Effect_Term ==  "mediumMononuclear cell medium" ~ "Medium: MCM",
                                                    Fixed_Effect_Term ==  "mediumRPMI + fetal bovine serum (FBS) (10%)" ~ "Medium: RPMI + FBS",
                                                    Fixed_Effect_Term ==  "microenvironmental_stimuliNone" ~ "Microenvironmental stimuli: None",
                                                    Fixed_Effect_Term ==  "microenvironmental_stimulitransiently cultured with feeder cells, activate samples with autologous BM T helper cells in the presence of IL-2 and a T-cell expansion cocktail" ~ "Microenvironmental Stimuli: \nCo-culture with activation",
                                                    Fixed_Effect_Term ==  "plate_readerPherastar, Cytation5, Insight, Tecan" ~ "Plate Reader: Pherastar\n or Cytation5",
                                                    Fixed_Effect_Term ==  "centrifugation_procedureLymphoPrepTM gradient centrifugation " ~ "Centrifugation Procedure: \nLymphoPrep\u2122 gradient",
                                                    Fixed_Effect_Term ==  "centrifugation_procedureSupernatant isolation at 300g 10min, density centrifugation at 400g for 20min without brake, afterwards always 300g 5 min" ~ "Centrifugation Procedure: 2ᵇ",
                                                    Fixed_Effect_Term ==  "cells5000" ~ "Number of Cells per Well: 5000",
                                                    Fixed_Effect_Term ==  "cell_counting_methodTrypan blue (Countess™ II FL Automated Cell Counter)" ~ "Cell Counting Method: Trypan blue",
                                                    Fixed_Effect_Term ==  "cell_counting_methodoptical density" ~ "Cell Counting Method: Optical density",
                                                    Fixed_Effect_Term ==  "plate_readerVictorX" ~ "Plate Reader: Victor X",
                                                    Fixed_Effect_Term ==  "plate_readerEnVision" ~ "Plate Reader: EnVision",
                                                    Fixed_Effect_Term ==  "plate_readerEnSight" ~ "Plate Reader: EnSight",
                                                    TRUE ~ Fixed_Effect_Term))

  df <- df %>% mutate(ref_group = case_when(ref_group ==  "a drug combination of flavopiridol, staurosporine and velcade" ~ "Positive Control: drug combination + Doses: 7 \n+ Readout: CellTiter96\n+ Cell Counting: Optical density     ", 
                                            ref_group ==  "1-2h after receiving" ~ "Time until sample usage <2h",
                                            ref_group ==  "HS-5 conditioned medium" ~ "HS-5 CM",
                                            ref_group ==  "HS-5 CM" ~ "HS-5 CM",
                                            ref_group ==  "Ficoll-Paque centrifugation" ~ "Ficoll-Paque",
                                            ref_group ==  "10000" ~ "10000",
                                            ref_group ==  "Countess" ~ "Countess",
                                            ref_group ==  "Biotek Synergy 2" ~ "Biotek Synergy 2",
                                                    TRUE ~ ref_group))
  
  df <- df %>% mutate(Response_metric = case_when(Response_metric ==  "DSS2_boxcox_sclaed2" ~ "DSS2: Box Cox Transformed and Scaled",
                                                  Response_metric == "auc_a_box_cox_scaled2" ~ "AUC A: Box cox Transformed and Scaled",
                                                  TRUE ~ Response_metric))
  df <- df[order(-df$Fixed_Effect_Coefficient), ]
  
  #returns df with summarised results
  return(df)
}

#----Model visualisation----
model_visualisation <- function(model){
  custom_palette <- colorRampPalette(RColorBrewer::brewer.pal(3, "Dark2"))(3)
  line_types <- rep("solid", 46)
  ggplotly(visualize(mixed_model_time_until_sample_usage_box_cox_scaled, plot = "model", line.size = 0.5)
           + scale_size(range = c(3, 5))
           + scale_color_manual(values = custom_palette)
           + guides(shape = "none"))
}

#----Model Diagnostics----
model_diagnostics <- function(model, output_dir = NULL, check_model = TRUE, model_name = "medium", metric = ""){
  print(output_dir)
  if (!is.null(output_dir)){
    if (!dir.exists(output_dir)){
      dir.create(output_dir, recursive = TRUE)
    } else {
      print("Dir already exists!")
    }
  }
  
  icc_model <- icc(model)
  icc_df <- data.frame(model = paste0(formula(model)[2], formula(model)[1], formula(model)[3]), icc = icc_model$icc, design_effect = icc_model$design.effect)
  print(formattable(icc_df))
  
  overall_df <- data.frame()
  
  if (check_model == TRUE){
    pp_check_plot <- check_model(model, check = "pp_check", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3", "black"))
    plot_pp_check <- plot(pp_check_plot)
    plot_theme <- theme(text = element_text(family = "Arial", color= "black", size = 10),
                        plot.title = element_text(hjust = 0.5, vjust = 1,color= "black"),
                        axis.title.x = element_text(family = "Arial", color = "black", size = 10),  
                        axis.title.y = element_text(family = "Arial", color = "black", size = 10),
                        axis.text.x = element_text(family = "Arial", color = "black", size = 8),  
                        axis.text.y = element_text(family = "Arial", color = "black", size = 8), 
                        plot.subtitle = element_blank(), 
                        legend.position = "top",  
                        legend.text = element_text(family = "Arial", color = "black", size = 10),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),  
                        panel.background = element_blank(),
                        strip.text = element_text(family = "Arial", color = "black", size = 8, face = "plain"),  
                        strip.background = element_blank(), 
                        plot.margin = margin(10, 10, 10, 10))
    plot_pp_check_adjusted <- plot_pp_check$PP_CHECK + plot_theme + labs(x = expression("BoxCox Transformed and Scaled DSS"[2]))
    print(plot_pp_check_adjusted)
    ggsave(paste0(output_dir, "pp_check_plot.png"), plot = plot_pp_check_adjusted,width = 9, height = 6.5, units="cm")
    
    linearity_plot <- check_model(model, check = "linearity", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
    plot_linearity <- plot(linearity_plot)
    plot_linearity_adjusted <- plot_linearity$NCV + plot_theme
    print(plot_linearity_adjusted)
    ggsave(paste0(output_dir, "plot_linearity_adjusted.png"), plot = plot_linearity_adjusted, width = 9, height = 6.5, units="cm")
    
    homogeneity_plot <- check_model(model, check = "homogeneity", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
    plot_homogeneity <- plot(homogeneity_plot)
    plot_homogeneity_adjusted <- plot_homogeneity$HOMOGENEITY + labs(y = expression(sqrt("|Std. Residuals|")))+ plot_theme
    print(plot_homogeneity)
    ggsave(paste0(output_dir, "plot_homogeneity_adjusted.png"), plot = plot_homogeneity_adjusted, width = 9, height = 6.5, units="cm")
    
    influential_plot <- check_model(model, check = "outliers", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3", "black"))
    plot_influential <- plot(influential_plot)  
    plot_influential_adjusted <- plot_influential$OUTLIERS + plot_theme
    print(plot_influential_adjusted)
    ggsave(paste0(output_dir, "plot_influential_adjusted.png"), plot = plot_influential_adjusted, width = 9, height = 6.5, units="cm")
    
    normality_plot <- check_model(model, check = "qq", panel = FALSE, title_size = 10, axis_title_size = 10, base_size = 0, colors = c("#fccde5", "#80b1d3"))
    plot_normality <- plot(normality_plot)
    plot_normality_adjusted <- plot_normality$QQ + plot_theme
    print(plot_normality_adjusted)
    ggsave(paste0(output_dir, "plot_normality_adjusted.png"), plot = plot_normality_adjusted, width = 7.94, height = 7.49, units="cm")
    
    random_effects_plot <- check_model(model, check = "reqq", panel = F, title_size = 10, axis_title_size = 10, base_size = 8, colors = c("#fccde5", "#80b1d3"))
    
    
    random_effects_plot$REQQ$drug$facet <- case_when("time_until_sample_usageHandled within 2-72 hours" %in% random_effects_plot$REQQ$drug$facet ~ factor(random_effects_plot$REQQ$drug$facet, 
                                                                    levels = c("(Intercept)", "time_until_sample_usageHandled within 2-72 hours"), 
                                                                    labels = c("Intercept", "Time until sample usage >2h")), 
                                                     "mediumMononuclear cell medium" %in% random_effects_plot$REQQ$drug$facet ~ factor(random_effects_plot$REQQ$drug$facet, 
                                                                                                                                       levels = c("(Intercept)", "mediumMononuclear cell medium", "mediumRPMI + fetal bovine serum (FBS) (10%)"), 
                                                                                                                                       labels = c("Intercept", "MCM", "RPMI + FBS")),
                                                     "cells5000" %in% random_effects_plot$REQQ$drug$facet ~ factor(random_effects_plot$REQQ$drug$facet, 
                                                                                                                                                          levels = c("(Intercept)", "cells5000"), 
                                                                                                                                                          labels = c("Intercept", "Cells: 5000")), 
                                                     "positive_controlBzCl" %in% random_effects_plot$REQQ$drug$facet ~ factor(random_effects_plot$REQQ$drug$facet, 
                                                                                                                   levels = c("(Intercept)", "positive_controlBzCl"), 
                                                                                                                   labels = c("Intercept", "Positive Control: BzCl + Doses: 5 \n+ Readout: CellTiter-Glo\n+ Cell Counting: Trypan blue")),
                                                     "microenvironmental_stimulitransiently cultured with feeder cells, activate samples with autologous BM T helper cells in the presence of IL-2 and a T-cell expansion cocktail" %in% random_effects_plot$REQQ$drug$facet ~ factor(random_effects_plot$REQQ$drug$facet, 
                                                                                                                              levels = c("(Intercept)", "microenvironmental_stimulitransiently cultured with feeder cells, activate samples with autologous BM T helper cells in the presence of IL-2 and a T-cell expansion cocktail"), 
                                                                                                                              labels = c("Intercept", "Microenvironmental Stimuli: \nCo-culture with activation")),
                                                     "centrifugation_procedureLymphoPrepTM gradient centrifugation" %in% random_effects_plot$REQQ$drug$facet ~ factor(random_effects_plot$REQQ$drug$facet, 
                                                                                                                   levels = c("(Intercept)", "centrifugation_procedureLymphoPrepTM gradient centrifugation"), 
                                                                                                                   labels = c("Intercept", "Centrifugation Procedure: \nLymphoPrep\u2122 gradient")),
                                                     "cells5000" %in% random_effects_plot$REQQ$drug$facet ~ factor(random_effects_plot$REQQ$drug$facet, 
                                                                                                                   levels = c("(Intercept)", "cells5000"), 
                                                                                                                   labels = c("Intercept", "Cells: 5000")),
                                                     .default = random_effects_plot$REQQ$drug$facet)
    plot_random_effects <- plot(random_effects_plot)
    plot_random_effects_adjusted_1 <- plot_random_effects[[1]] + plot_theme + labs(title = "Normality of Random Effects - Patient") 
    print(plot_random_effects_adjusted_1)
    ggsave(paste0(output_dir, "random_effects_check_1.png"), plot = plot_random_effects_adjusted_1, width = 9, height = 6.5, units="cm")
    
    plot_random_effects_adjusted_2 <- plot_random_effects[[2]] + plot_theme + labs(title = "Normality of Random Effects - Drug")
    print(plot_random_effects_adjusted_2)
    ggsave(paste0(output_dir, "random_effects_check_2.png"), plot = plot_random_effects_adjusted_2, width = 9, height = 6.5, units="cm")
    
    
    marginal_residuals <- model.response(model.frame(model)) - predict(model, re.form = NA)
    
    # Create the QQ plot styled like check_model
    plot_marginal_residuals <- ggplot(data = data.frame(resid = marginal_residuals), aes(sample = resid)) +
      stat_qq(size = 0.8, color = "#80b1d3") +
      stat_qq_line(color = "#fccde5") +
      labs(
        title = "Normality of Marginal Residuals",
        x = "Theoretical Quantiles",
        y = "Sample Quantiles"
      ) +
      theme_minimal(base_size = 8, base_family = "Arial") +
      plot_theme
    
    # Show the plot
    print(plot_marginal_residuals)
    

    label = paste0(
                 "Diagnostic plots for modelling differences in BoxCox transformed DSS2 when using different ", model_name, ".\na: Posterior predictive check for comparing distributions of observed scaled DSS2 data with that simulated by the\nposterior distribution of the ", model_name, "-specific LMM.\nb: Assessment of influential variables in the ",model_name,"-specific LMM.\nThe black labelled points indicate marginally influential points, and the pink dotted curves show Cook’s distance threshold of 0.5.\nc: Evaluation of the model assumption of homogeneity of variance in the ",model_name,"-specific LMM.\nd: Linearity plot checks the assumption of linear relationship, shows weather the predictor has a linear relationship\nwith the outcome, a straight horizontal line indicates that the model specification is good.\ne-f: Q-Q plots for confirming the normality of marginal and conditional residuals.\ng-h: Q-Q plots for confirming the normality of random effects for patients (g)\nand drugs (h)  in the ",model_name, "-specific LMM."
               )
    
    title_plot <- paste(metric, model_name) 
    
    library(patchwork)
    comboplot2 <- (
      (plot_pp_check_adjusted | plot_influential_adjusted) / 
      (plot_homogeneity_adjusted | plot_linearity_adjusted)/
      (plot_normality_adjusted | plot_marginal_residuals)/
      (plot_random_effects_adjusted_1 | plot_random_effects_adjusted_2)
    )+ plot_annotation(tag_levels = "a", title = title_plot, caption = label, theme = theme(
      plot.title = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0, size = 10)  # hjust = 1 = right aligned
    ))
    #plot_layout(heights = c(1, 1, 1, 1,1), widths = c(1,1))
    print(comboplot2)
    ggsave(paste0(output_dir, "final_plot.png"), plot = comboplot2, width = 20.76, height = 28.36, units="cm")
    print(random_effects_plot$REQQ$drug$facet)
  }
  # # Extract residuals
  # residuals <- residuals(model)
  # 
  # # Extract random effects
  # random_effects <- ranef(model)
  # 
  # # Shapiro-Wilk test for random intercepts
  # random_intercepts_patient <- random_effects[[1]][, "(Intercept)"]
  # shapiro_intercept_patient <- shapiro.test(random_intercepts_patient)
  # print(shapiro_intercept_patient)
  # 
  # random_intercepts_drug <- random_effects[[2]][, "(Intercept)"]
  # shapiro_intercept_drug <- shapiro.test(random_intercepts_drug)
  # print(shapiro_intercept_drug)
  # 
  # # Shapiro-Wilk test for random slopes (if they exist)
  # random_slopes <- random_effects[[2]]
  # shapiro_slope <- NA  # Default in case no slope
  # for(n in names(random_slopes)){
  #   print(n)
  #   if(n != "(Intercept)"){
  #     shapiro_slope <- shapiro.test(random_slopes[,n])
  #     print(shapiro_slope)
  #   }
  # }
  # 
  # 
  # # Homogeneity of variance test (Levene's test)
  # fitted_values <- fitted(model)
  # levene_test <- leveneTest(residuals ~ cut(fitted_values, breaks = 4))  # Groups fitted values into 4 bins
  # 
  # # Create and return results as a dataframe
  # diagnostics_df <- data.frame(
  #   Test = c("Shapiro-Wilk Intercepts Patient", "Shapiro-Wilk Intercepts Drug", "Shapiro-Wilk Slopes", "Levene's Test"),
  #   Statistic = c(shapiro_intercept_patient$statistic, shapiro_intercept_drug$statistic, shapiro_slope$statistic, levene_test$`F value`[1]),
  #   P_Value = c(shapiro_intercept_patient$p.value, shapiro_intercept_drug$p.value, shapiro_slope$p.value, levene_test$`Pr(>F)`[1])
  # )
  # 
  # return(diagnostics_df)
}

# 
# library(HLMdiag)
# # Extract residuals and fitted values
# residuals_level_1 <- hlm_resid(model, level=1)
# residuals_level_2 <- hlm_resid(model, level="drug")
# 
# 
# ggplot(residuals_level_1 , aes(x = .mar.fitted, y = .mar.resid)) +
#   geom_point() +
#   geom_smooth() + 
#   labs(title = "Margial Residuals vs Fitted Values")
# 
# qqnorm(residuals_level_1$.mar.resid)
# qqline(residuals_level_1$.mar.resid, col = "darkgreen")
# 
# ggplot(residuals_level_1 , aes(x = drug, y = .resid)) +
#   geom_boxplot() +
#   labs(title = "Residuals Level-1 for each drug - Within-group homoscdasticity check") +
#   theme(axis.text.x = element_text(angle = 90))
# 
# ggplot(residuals_level_1 , aes(x = .fitted, y = .resid)) +
#   geom_point() +
#   geom_smooth() + 
#   labs(title = "Residuals Level-1 - Between group homoscdasticity check") +
#   theme(axis.text.x = element_text(angle = 90)) +
#   facet_wrap(~drug)
# 
# ggplot(residuals_level_1 , aes(x = .fitted, y = .resid, color=drug)) +
#   geom_point() +
#   geom_smooth() + 
#   labs(title = "Residuals Level-1 - Between group homoscdasticity check") +
#   theme(axis.text.x = element_text(angle = 90))
# 
# 
# library(DHARMa)
# 
# # Simulate residuals from the model
# sim_res <- simulateResiduals(fittedModel = model)
# 
# # Plot residuals
# plot(sim_res)
# 
# # Test for homoscedasticity
# testResiduals(sim_res)
# 
# 
# qqnorm(residuals_level_1$.resid)
# qqline(residuals_level_1$.resid, col = "darkgreen")
# 
# qqnorm(residuals_level_2$.ls.intercept)
# qqline(residuals_level_2$.ls.intercept, col = "darkgreen")
# 
# qqnorm(residuals_level_2$.ls.time_until_sample_usage_handled_within_2_72_hours)
# qqline(residuals_level_2$.ls.time_until_sample_usage_handled_within_2_72_hours, col = "darkgreen")
# 
# 
# ggplot(residuals_level_2 , aes(x = drug, y = .ls.intercept)) +
#   geom_boxplot() +
#   labs(title = "Residuals vs Fitted Values")
# 
# 
# influence_data <- influence(model, group = "drug")
# plot(influence_data, which = "cook", main = "Influence Plot: Cook's Distance")


