# server_gmm.R
# This module contains the logic for the "Subpopulation Detection (GMM)" tab.
# It handles data upload, running the GMM analysis, and rendering the results.

# =========================================================================
# UTILITY FUNCTIONS FOR GMM ANALYSIS
# =========================================================================

# Z-transform a numeric vector (standardization)
# @param x: A numeric vector.
# @return: The standardized numeric vector.
z_transform <- function(x) {
  if (sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

# Function to apply Yeo-Johnson Transformation conditionally
# It checks the skewness of a data vector and applies the transformation if the absolute skewness is above a threshold.
# @param data_vector: A numeric vector to transform.
# @param skewness_threshold: The threshold for applying the transformation.
# @return: A list containing the transformed data and a boolean indicating if a transformation was applied.
apply_conditional_yeo_johnson <- function(data_vector, skewness_threshold = 0.5) {
  transformed_data <- data_vector
  transformation_applied <- FALSE
  skew <- moments::skewness(data_vector, na.rm = TRUE)

  if (abs(skew) > skewness_threshold) {
    tryCatch({
      pt_result <- powerTransform(data_vector)
      lambda <- pt_result$lambda
      transformed_data <- car::yjPower(data_vector, lambda)
      transformation_applied <- TRUE
      message(paste("Yeo-Johnson transformation applied (skewness=", round(skew, 2), ")"))
    }, error = function(e) {
      warning(paste("Could not apply Yeo-Johnson transformation:", e$message))
    })
  } else {
    message(paste("Yeo-Johnson transformation not needed (skewness=", round(skew, 2), ")"))
  }

  return(list(transformed_data = transformed_data, transformation_applied = transformation_applied))
}

# Function to run GMM analysis using mclust with a selectable criterion (BIC or ICL)
# @param data_mat: A numeric matrix or data frame for clustering.
# @param G_range: A range of component numbers to test (e.g., 2:5).
# @param criterion: The model selection criterion ("BIC" or "ICL").
# @return: An Mclust object representing the best-fit model.
run_gmm_with_criterion <- function(data_mat, G_range = 2:5, criterion = "BIC") {
  if (!is.matrix(data_mat) && !is.data.frame(data_mat)) {
    stop("Input data_mat must be a matrix or data frame for GMM analysis.")
  }
  if (!all(sapply(data_mat, is.numeric))) {
    stop("All columns in data_mat must be numeric.")
  }
  if (any(is.na(data_mat))) {
    stop("Input data_mat contains NA values. Please remove or impute before clustering.")
  }
  
  multivariate_model_names <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EVE", "VEE", "VVV")
  
  tryCatch({
    if (criterion == "BIC") {
      gmm_model <- Mclust(data_mat, G = G_range, modelNames = multivariate_model_names)
      return(gmm_model)
    } else if (criterion == "ICL") {
      mclust_all <- Mclust(data_mat, G = G_range, modelNames = multivariate_model_names, warn = FALSE)
      if (is.null(mclust_all)) {
        stop("Mclust could not be fitted to the data.")
      }
      icl_values <- mclustICL(mclust_all)
      
      best_model_index <- which.max(icl_values)
      
      best_G <- icl_values[best_model_index, "G"]
      best_model_name <- rownames(icl_values)[best_model_index]

      gmm_model <- Mclust(data_mat, G = best_G, modelNames = best_model_name)
      return(gmm_model)

    } else {
      stop("Invalid criterion selected. Please choose 'BIC' or 'ICL'.")
    }
  }, error = function(e) {
    warning(paste("GMM run failed:", e$message))
    return(NULL)
  })
}

# Function to assign clusters back to the original data frame
# It adds a 'cluster' column to the data frame based on the GMM model's classification.
# @param df: The original data frame.
# @param gmm_model: The Mclust model object with a 'classification' property.
# @return: The data frame with an added 'cluster' column.
assign_clusters <- function(df, gmm_model) {
  if (is.null(gmm_model) || is.null(gmm_model$classification)) {
    warning("GMM model or classification is NULL. Cannot assign clusters.")
    return(df)
  }
  df$cluster <- gmm_model$classification
  return(df)
}

# Function to plot age vs HGB colored by cluster
# It generates a ggplot scatter plot with confidence ellipses and cluster means.
# @param df: The data frame with 'Age', 'HGB', 'Gender', and 'cluster' columns.
# @param male_hgb_transformed: A boolean indicating if male HGB data was transformed.
# @param female_hgb_transformed: A boolean indicating if female HGB data was transformed.
# @return: A ggplot object.
plot_age_hgb <- function(df, male_hgb_transformed, female_hgb_transformed) {
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No GMM data available for plotting.", size = 6, color = "grey50"))
  }

  plot_title <- "HGB vs Age by Subpopulation Cluster"
  plot_subtitle <- ""
  if (male_hgb_transformed) {
    plot_subtitle <- paste(plot_subtitle, "Male HGB transformed for GMM.", sep="\n")
  }
  if (female_hgb_transformed) {
    plot_subtitle <- paste(plot_subtitle, "Female HGB transformed for GMM.", sep="\n")
  }
  
  cluster_means <- df %>%
    group_by(Gender, cluster) %>%
    summarise(mean_Age = mean(Age, na.rm = TRUE),
              mean_HGB = mean(HGB, na.rm = TRUE),
              .groups = 'drop')

  ggplot(df, aes(x = Age, y = HGB, color = factor(cluster))) +
    geom_point(position = position_jitter(width = 0.2, height = 0.2), alpha = 0.6) +
    stat_ellipse(geom = "polygon", aes(fill = factor(cluster)), alpha = 0.2, show.legend = FALSE, level = 0.95) +
    geom_point(data = cluster_means, aes(x = mean_Age, y = mean_HGB), shape = 4, size = 5, color = "red", stroke = 2) +
    facet_wrap(~Gender) +
    theme_minimal() +
    labs(title = plot_title,
         subtitle = plot_subtitle,
         x = "Age", y = "HGB", color = "Cluster") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "bottom")
}

# =========================================================================
# MAIN SERVER LOGIC
# =========================================================================

gmmServer <- function(input, output, session, gmm_uploaded_data_rv, gmm_processed_data_rv, gmm_transformation_details_rv, gmm_models_rv, message_rv, analysis_running_rv) {

  # Reactive values to hold models for each criterion
  gmm_models_bic_rv <- reactiveVal(list(male = NULL, female = NULL))
  gmm_models_icl_rv <- reactiveVal(list(male = NULL, female = NULL))

  # Helper function to guess column names (could be moved to a shared utils file)
  guess_column <- function(cols_available, common_names) {
    for (name in common_names) {
      match_idx <- grep(paste0("^", name, "$"), cols_available, ignore.case = TRUE)
      if (length(match_idx) > 0) {
        return(cols_available[match_idx[1]])
      }
    }
    return("")
  }

  # Observer for GMM file upload
  # Reads the uploaded Excel file and updates column selectors based on likely names
  observeEvent(input$gmm_file_upload, {
    req(input$gmm_file_upload)
    tryCatch({
      data <- readxl::read_excel(input$gmm_file_upload$datapath)
      gmm_uploaded_data_rv(data)
      message_rv(list(text = "GMM data uploaded successfully.", type = "success"))

      col_names <- colnames(data)
      all_col_choices_with_none <- c("None" = "", col_names)
      
      updateSelectInput(session, "gmm_hgb_col", choices = all_col_choices_with_none, selected = guess_column(col_names, c("HGB", "hgb", "HB", "hb")))
      updateSelectInput(session, "gmm_age_col", choices = all_col_choices_with_none, selected = guess_column(col_names, c("Age", "age", "leeftijd")))
      updateSelectInput(session, "gmm_gender_col", choices = all_col_choices_with_none, selected = guess_column(col_names, c("Gender", "gender", "Sex", "sex", "geslacht")))
      
    }, error = function(e) {
      message_rv(list(text = paste("Error reading GMM file:", e$message), type = "error"))
      gmm_uploaded_data_rv(NULL)
    })
  })

  # Observer for GMM analysis button
  # This is the core logic for the GMM tab, running the analysis with progress updates
  observeEvent(input$run_gmm_analysis_btn, {
    req(gmm_uploaded_data_rv(), input$gmm_hgb_col, input$gmm_age_col, input$gmm_gender_col)

    if (analysis_running_rv()) {
      message_rv(list(text = "An analysis is already running. Please wait.", type = "warning"))
      return(NULL)
    }

    analysis_running_rv(TRUE)
    shinyjs::disable("tabs")

    withProgress(message = 'Running GMM Analysis', value = 0, {
      incProgress(0.1, detail = "Loading data...")

      data <- gmm_uploaded_data_rv()
      hgb_col <- input$gmm_hgb_col
      age_col <- input$gmm_age_col
      gender_col <- input$gmm_gender_col

      if (!all(c(hgb_col, age_col, gender_col) %in% names(data))) {
        message_rv(list(text = "Selected columns not found in data. Please check selections.", type = "error"))
        analysis_running_rv(FALSE)
        shinyjs::enable("tabs")
        return(NULL)
      }

      gmm_data <- data %>%
        dplyr::select(HGB = !!sym(hgb_col), Age = !!sym(age_col), Gender_orig = !!sym(gender_col)) %>%
        na.omit()

      if (nrow(gmm_data) == 0) {
        message_rv(list(text = "No complete rows for GMM after NA removal. Check data or selections.", type = "error"))
        analysis_running_rv(FALSE)
        shinyjs::enable("tabs")
        return(NULL)
      }

      incProgress(0.2, detail = "Splitting data by gender and transforming...")

      gmm_data <- gmm_data %>%
        mutate(Gender = case_when(
          grepl("male|m", Gender_orig, ignore.case = TRUE) ~ "Male",
          grepl("female|f", Gender_orig, ignore.case = TRUE) ~ "Female",
          TRUE ~ "Other"
        )) %>%
        filter(Gender %in% c("Male", "Female"))

      male_data <- gmm_data %>% filter(Gender == "Male")
      female_data <- gmm_data %>% filter(Gender == "Female")

      combined_clustered_data_bic <- tibble()
      combined_clustered_data_icl <- tibble()
      male_hgb_transformed_flag <- FALSE
      female_hgb_transformed_flag <- FALSE
      male_gmm_model_bic <- NULL
      female_gmm_model_bic <- NULL
      male_gmm_model_icl <- NULL
      female_gmm_model_icl <- NULL
      
      if (nrow(male_data) > 0) {
        yj_result_male <- apply_conditional_yeo_johnson(male_data$HGB)
        male_data$HGB_transformed <- yj_result_male$transformed_data
        male_hgb_transformed_flag <- yj_result_male$transformation_applied
        male_data$HGB_z <- z_transform(male_data$HGB_transformed)
        male_data$Age_z <- z_transform(male_data$Age)
        
        incProgress(0.2, detail = "Running GMM for Male data (BIC)...")
        tryCatch({
          male_gmm_model_bic <- run_gmm_with_criterion(male_data %>% dplyr::select(HGB = HGB_z, Age = Age_z), criterion = "BIC")
          if (!is.null(male_gmm_model_bic)) {
            male_data_bic <- assign_clusters(male_data, male_gmm_model_bic)
            male_data_bic$cluster <- as.factor(male_data_bic$cluster)
            combined_clustered_data_bic <- bind_rows(combined_clustered_data_bic, male_data_bic %>% dplyr::select(HGB, Age, Gender, cluster))
          }
        }, error = function(e) {
          message_rv(list(text = paste("Error running BIC GMM for male data:", e$message), type = "error"))
        })

        incProgress(0.2, detail = "Running GMM for Male data (ICL)...")
        tryCatch({
          male_gmm_model_icl <- run_gmm_with_criterion(male_data %>% dplyr::select(HGB = HGB_z, Age = Age_z), criterion = "ICL")
          if (!is.null(male_gmm_model_icl)) {
            male_data_icl <- assign_clusters(male_data, male_gmm_model_icl)
            male_data_icl$cluster <- as.factor(male_data_icl$cluster)
            combined_clustered_data_icl <- bind_rows(combined_clustered_data_icl, male_data_icl %>% dplyr::select(HGB, Age, Gender, cluster))
          }
        }, error = function(e) {
          message_rv(list(text = paste("Error running ICL GMM for male data:", e$message), type = "error"))
        })
      }

      if (nrow(female_data) > 0) {
        yj_result_female <- apply_conditional_yeo_johnson(female_data$HGB)
        female_data$HGB_transformed <- yj_result_female$transformed_data
        female_hgb_transformed_flag <- yj_result_female$transformation_applied
        female_data$HGB_z <- z_transform(female_data$HGB_transformed)
        female_data$Age_z <- z_transform(female_data$Age)
        
        incProgress(0.2, detail = "Running GMM for Female data (BIC)...")
        tryCatch({
          female_gmm_model_bic <- run_gmm_with_criterion(female_data %>% dplyr::select(HGB = HGB_z, Age = Age_z), criterion = "BIC")
          if (!is.null(female_gmm_model_bic)) {
            female_data_bic <- assign_clusters(female_data, female_gmm_model_bic)
            female_data_bic$cluster <- as.factor(female_data_bic$cluster)
            combined_clustered_data_bic <- bind_rows(combined_clustered_data_bic, female_data_bic %>% dplyr::select(HGB, Age, Gender, cluster))
          }
        }, error = function(e) {
          message_rv(list(text = paste("Error running BIC GMM for female data:", e$message), type = "error"))
        })

        incProgress(0.2, detail = "Running GMM for Female data (ICL)...")
        tryCatch({
          female_gmm_model_icl <- run_gmm_with_criterion(female_data %>% dplyr::select(HGB = HGB_z, Age = Age_z), criterion = "ICL")
          if (!is.null(female_gmm_model_icl)) {
            female_data_icl <- assign_clusters(female_data, female_gmm_model_icl)
            female_data_icl$cluster <- as.factor(female_data_icl$cluster)
            combined_clustered_data_icl <- bind_rows(combined_clustered_data_icl, female_data_icl %>% dplyr::select(HGB, Age, Gender, cluster))
          }
        }, error = function(e) {
          message_rv(list(text = paste("Error running ICL GMM for female data:", e$message), type = "error"))
        })
      }
      
      gmm_models_bic_rv(list(male = male_gmm_model_bic, female = female_gmm_model_bic))
      gmm_models_icl_rv(list(male = male_gmm_model_icl, female = female_gmm_model_icl))

      gmm_transformation_details_rv(list(male_hgb_transformed = male_hgb_transformed_flag, female_hgb_transformed = female_hgb_transformed_flag))

      if (nrow(combined_clustered_data_bic) > 0 || nrow(combined_clustered_data_icl) > 0) {
        gmm_processed_data_rv(list(bic = combined_clustered_data_bic, icl = combined_clustered_data_icl))
        message_rv(list(text = "GMM analysis complete!", type = "success"))
      } else {
        message_rv(list(text = "No data available after GMM processing for plotting/summary.", type = "error"))
        gmm_processed_data_rv(NULL)
      }

      incProgress(0.1, detail = "Generating plots and summaries...")
    })

    analysis_running_rv(FALSE)
    shinyjs::enable("tabs")
  })

  # Observer for the Reset button on the GMM tab
  observeEvent(input$reset_gmm_analysis_btn, {
    gmm_uploaded_data_rv(NULL)
    gmm_processed_data_rv(NULL)
    gmm_transformation_details_rv(list(male_hgb_transformed = FALSE, female_hgb_transformed = FALSE))
    gmm_models_bic_rv(list(male = NULL, female = NULL))
    gmm_models_icl_rv(list(male = NULL, female = NULL))
    shinyjs::reset("gmm_file_upload")
    output$gmm_results_ui <- renderUI(NULL)
    message_rv(list(text = "GMM data and results reset.", type = "info"))
    
    updateSelectInput(session, "gmm_hgb_col", choices = c("None" = ""), selected = "")
    updateSelectInput(session, "gmm_age_col", choices = c("None" = ""), selected = "")
    updateSelectInput(session, "gmm_gender_col", choices = c("None" = ""), selected = "")
  })

  output$gmm_results_ui <- renderUI({
    plot_data_bic <- gmm_processed_data_rv()$bic
    plot_data_icl <- gmm_processed_data_rv()$icl

    if (is.null(plot_data_bic) && is.null(plot_data_icl)) {
      return(NULL)
    }

    tagList(
      fluidRow(
        column(6,
          div(class = "output-box",
              h4("BIC Criterion Results"),
              plotOutput("gmm_model_selection_plot_bic", height = "400px"),
              plotOutput("plot_output_gmm_bic", height = "600px"),
              verbatimTextOutput("gmm_summary_output_bic"),
              tableOutput("gmm_age_group_summary_output_bic")
          )
        ),
        column(6,
          div(class = "output-box",
              h4("ICL Criterion Results"),
              plotOutput("gmm_model_selection_plot_icl", height = "400px"),
              plotOutput("plot_output_gmm_icl", height = "600px"),
              verbatimTextOutput("gmm_summary_output_icl"),
              tableOutput("gmm_age_group_summary_output_icl")
          )
        )
      )
    )
  })

  output$gmm_model_selection_plot_bic <- renderPlot({
    models <- gmm_models_bic_rv()
    if (is.null(models$male) && is.null(models$female)) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No GMM models available for plotting.", size = 6, color = "grey50"))
    }
    
    if (!is.null(models$male) && !is.null(models$female)) {
      par(mfrow = c(1, 2))
      plot(models$male, what = "BIC", main = "Male - BIC Plot")
      plot(models$female, what = "BIC", main = "Female - BIC Plot")
      par(mfrow = c(1, 1))
    } else if (!is.null(models$male)) {
      plot(models$male, what = "BIC", main = "Male - BIC Plot")
    } else if (!is.null(models$female)) {
      plot(models$female, what = "BIC", main = "Female - BIC Plot")
    }
  })

  output$gmm_model_selection_plot_icl <- renderPlot({
    models <- gmm_models_icl_rv()
    if (is.null(models$male) && is.null(models$female)) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No GMM models available for plotting.", size = 6, color = "grey50"))
    }
    
    if (!is.null(models$male) && !is.null(models$female)) {
      par(mfrow = c(1, 2))
      plot(mclustICL(models$male), main = "Male - ICL Plot")
      plot(mclustICL(models$female), main = "Female - ICL Plot")
      par(mfrow = c(1, 1))
    } else if (!is.null(models$male)) {
      plot(mclustICL(models$male), main = "Male - ICL Plot")
    } else if (!is.null(models$female)) {
      plot(mclustICL(models$female), main = "Female - ICL Plot")
    }
  })

  output$plot_output_gmm_bic <- renderPlot({
    plot_data <- gmm_processed_data_rv()$bic
    if (is.null(plot_data) || nrow(plot_data) == 0) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No GMM data available for plotting.", size = 6, color = "grey50"))
    }
    plot_age_hgb(plot_data,
                 male_hgb_transformed = gmm_transformation_details_rv()$male_hgb_transformed,
                 female_hgb_transformed = gmm_transformation_details_rv()$female_hgb_transformed)
  })

  output$plot_output_gmm_icl <- renderPlot({
    plot_data <- gmm_processed_data_rv()$icl
    if (is.null(plot_data) || nrow(plot_data) == 0) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No GMM data available for plotting.", size = 6, color = "grey50"))
    }
    plot_age_hgb(plot_data,
                 male_hgb_transformed = gmm_transformation_details_rv()$male_hgb_transformed,
                 female_hgb_transformed = gmm_transformation_details_rv()$female_hgb_transformed)
  })


  output$gmm_summary_output_bic <- renderPrint({
    plot_data <- gmm_processed_data_rv()$bic
    models <- gmm_models_bic_rv()
    
    if (is.null(plot_data) || nrow(plot_data) == 0) {
      return("No GMM analysis results to display.")
    }

    cat("--- GMM Analysis Summary (BIC Criterion) ---\n")
    
    if (!is.null(models$male) && !inherits(models$male, "try-error")) {
        cat("\n--- Male Subpopulations ---\n")
        
        print(summary(models$male))

        num_clusters <- models$male$G
        for (i in 1:num_clusters) {
            cat(paste0("Cluster ", i, ":\n"))
            cat(paste0("  Proportion: ", round(models$male$parameters$pro[i], 3), "\n"))
            
            male_cluster_data <- plot_data %>% filter(Gender == "Male", cluster == i)
            mean_hgb <- mean(male_cluster_data$HGB, na.rm = TRUE)
            mean_age <- mean(male_cluster_data$Age, na.rm = TRUE)
            sd_hgb <- sd(male_cluster_data$HGB, na.rm = TRUE)
            sd_age <- sd(male_cluster_data$Age, na.rm = TRUE)
            
            cat(paste0("  Mean HGB: ", round(mean_hgb, 3), "\n"))
            cat(paste0("  Mean Age: ", round(mean_age, 3), "\n"))
            cat(paste0("  Std Dev HGB: ", round(sd_hgb, 3), "\n"))
            cat(paste0("  Std Dev Age: ", round(sd_age, 3), "\n"))
            
            if (!is.na(sd_age)) {
              lower_age <- round(mean_age - 2 * sd_age, 1)
              upper_age <- round(mean_age + 2 * sd_age, 1)
              cat(paste0("  Estimated Age Range (Mean +/- 2SD): [", max(0, lower_age), " to ", upper_age, "] years\n"))
            } else {
              cat("  Estimated Age Range: N/A (Std Dev Age problematic)\n")
            }
            cat("\n")
        }
    } else {
        cat("No male subpopulations detected.\n")
    }
    
    if (!is.null(models$female) && !inherits(models$female, "try-error")) {
        cat("\n--- Female Subpopulations ---\n")
        
        print(summary(models$female))

        num_clusters <- models$female$G
        for (i in 1:num_clusters) {
            cat(paste0("Cluster ", i, ":\n"))
            cat(paste0("  Proportion: ", round(models$female$parameters$pro[i], 3), "\n"))
            
            female_cluster_data <- plot_data %>% filter(Gender == "Female", cluster == i)
            mean_hgb <- mean(female_cluster_data$HGB, na.rm = TRUE)
            mean_age <- mean(female_cluster_data$Age, na.rm = TRUE)
            sd_hgb <- sd(female_cluster_data$HGB, na.rm = TRUE)
            sd_age <- sd(female_cluster_data$Age, na.rm = TRUE)
            
            cat(paste0("  Mean HGB: ", round(mean_hgb, 3), "\n"))
            cat(paste0("  Mean Age: ", round(mean_age, 3), "\n"))
            cat(paste0("  Std Dev HGB: ", round(sd_hgb, 3), "\n"))
            cat(paste0("  Std Dev Age: ", round(sd_age, 3), "\n"))
            
            if (!is.na(sd_age)) {
              lower_age <- round(mean_age - 2 * sd_age, 1)
              upper_age <- round(mean_age + 2 * sd_age, 1)
              cat(paste0("  Estimated Age Range (Mean +/- 2SD): [", max(0, lower_age), " to ", upper_age, "] years\n"))
            } else {
              cat("  Estimated Age Range: N/A (Std Dev Age problematic)\n")
            }
            cat("\n")
        }
    } else {
        cat("No female subpopulations detected.\n")
    }

    if (gmm_transformation_details_rv()$male_hgb_transformed || gmm_transformation_details_rv()$female_hgb_transformed) {
      cat("\nNote: HGB values were transformed (Yeo-Johnson) for GMM input due to skewness. Reported HGB values are original.\n")
    }
  })

  output$gmm_summary_output_icl <- renderPrint({
    plot_data <- gmm_processed_data_rv()$icl
    models <- gmm_models_icl_rv()
    
    if (is.null(plot_data) || nrow(plot_data) == 0) {
      return("No GMM analysis results to display.")
    }

    cat("--- GMM Analysis Summary (ICL Criterion) ---\n")
    
    if (!is.null(models$male) && !inherits(models$male, "try-error")) {
        cat("\n--- Male Subpopulations ---\n")
        
        print(summary(models$male))

        num_clusters <- models$male$G
        for (i in 1:num_clusters) {
            cat(paste0("Cluster ", i, ":\n"))
            cat(paste0("  Proportion: ", round(models$male$parameters$pro[i], 3), "\n"))
            
            male_cluster_data <- plot_data %>% filter(Gender == "Male", cluster == i)
            mean_hgb <- mean(male_cluster_data$HGB, na.rm = TRUE)
            mean_age <- mean(male_cluster_data$Age, na.rm = TRUE)
            sd_hgb <- sd(male_cluster_data$HGB, na.rm = TRUE)
            sd_age <- sd(male_cluster_data$Age, na.rm = TRUE)
            
            cat(paste0("  Mean HGB: ", round(mean_hgb, 3), "\n"))
            cat(paste0("  Mean Age: ", round(mean_age, 3), "\n"))
            cat(paste0("  Std Dev HGB: ", round(sd_hgb, 3), "\n"))
            cat(paste0("  Std Dev Age: ", round(sd_age, 3), "\n"))
            
            if (!is.na(sd_age)) {
              lower_age <- round(mean_age - 2 * sd_age, 1)
              upper_age <- round(mean_age + 2 * sd_age, 1)
              cat(paste0("  Estimated Age Range (Mean +/- 2SD): [", max(0, lower_age), " to ", upper_age, "] years\n"))
            } else {
              cat("  Estimated Age Range: N/A (Std Dev Age problematic)\n")
            }
            cat("\n")
        }
    } else {
        cat("No male subpopulations detected.\n")
    }
    
    if (!is.null(models$female) && !inherits(models$female, "try-error")) {
        cat("\n--- Female Subpopulations ---\n")
        
        print(summary(models$female))

        num_clusters <- models$female$G
        for (i in 1:num_clusters) {
            cat(paste0("Cluster ", i, ":\n"))
            cat(paste0("  Proportion: ", round(models$female$parameters$pro[i], 3), "\n"))
            
            female_cluster_data <- plot_data %>% filter(Gender == "Female", cluster == i)
            mean_hgb <- mean(female_cluster_data$HGB, na.rm = TRUE)
            mean_age <- mean(female_cluster_data$Age, na.rm = TRUE)
            sd_hgb <- sd(female_cluster_data$HGB, na.rm = TRUE)
            sd_age <- sd(female_cluster_data$Age, na.rm = TRUE)
            
            cat(paste0("  Mean HGB: ", round(mean_hgb, 3), "\n"))
            cat(paste0("  Mean Age: ", round(mean_age, 3), "\n"))
            cat(paste0("  Std Dev HGB: ", round(sd_hgb, 3), "\n"))
            cat(paste0("  Std Dev Age: ", round(sd_age, 3), "\n"))
            
            if (!is.na(sd_age)) {
              lower_age <- round(mean_age - 2 * sd_age, 1)
              upper_age <- round(mean_age + 2 * sd_age, 1)
              cat(paste0("  Estimated Age Range (Mean +/- 2SD): [", max(0, lower_age), " to ", upper_age, "] years\n"))
            } else {
              cat("  Estimated Age Range: N/A (Std Dev Age problematic)\n")
            }
            cat("\n")
        }
    } else {
        cat("No female subpopulations detected.\n")
    }

    if (gmm_transformation_details_rv()$male_hgb_transformed || gmm_transformation_details_rv()$female_hgb_transformed) {
      cat("\nNote: HGB values were transformed (Yeo-Johnson) for GMM input due to skewness. Reported HGB values are original.\n")
    }
  })
}