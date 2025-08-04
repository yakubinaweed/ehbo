# server_main.R
# This module contains the logic for the "Main Analysis" tab.
# It handles data upload, filtering, running the refineR model, and rendering the results.

# =========================================================================
# UTILITY FUNCTIONS FOR MAIN ANALYSIS
# =========================================================================

# Function to filter data based on gender and age
# It takes a data frame, gender choice, age range, and column names for gender and age.
# It returns a filtered data frame with a standardized 'Gender_Standardized' column.
filter_data <- function(data, gender_choice, age_min, age_max, col_gender, col_age) {
  if (!col_gender %in% names(data) || !col_age %in% names(data)) {
    stop("Gender or age column not found in data.")
  }

  filtered_data <- data %>%
    filter(!!sym(col_age) >= age_min & !!sym(col_age) <= age_max)
  
  filtered_data <- filtered_data %>%
    mutate(Gender_Standardized = case_when(
      grepl("male|m|man|jongen(s)?|heren|mannelijk(e)?", !!sym(col_gender), ignore.case = TRUE) ~ "Male",
      grepl("female|f|vrouw(en)?|v|meisje(s)?|dame|mevr|vrouwelijke", !!sym(col_gender), ignore.case = TRUE) ~ "Female",
      TRUE ~ "Other"
    ))
  
  if (gender_choice != "Both") {
    filtered_data <- filtered_data %>%
      filter(Gender_Standardized == case_when(
        gender_choice == "M" ~ "Male",
        gender_choice == "F" ~ "Female"
      ))
  }

  return(filtered_data)
}

# Function to generate a safe filename for plots
# It creates a unique filename using a title, datestamp, and timestamp.
generate_safe_filename <- function(plot_title, base_path, extension = "png") {
  safe_title <- gsub("[^a-zA-Z0-9_-]", "_", plot_title)
  datestamp <- format(Sys.Date(), "%Y%m%d")
  timestamp <- format(Sys.time(), "%H%M%S")
  file.path(base_path, paste0(safe_title, "_", datestamp, "-", timestamp, ".", extension))
}

# =========================================================================
# MAIN SERVER LOGIC
# =========================================================================

mainServer <- function(input, output, session, data_reactive, selected_dir_reactive, message_rv, analysis_running_rv) {

  # Reactive value to hold the refineR model result
  refiner_model_rv <- reactiveVal(NULL)
  # Reactive value to hold the plot title
  plot_title_rv <- reactiveVal("")

  # Helper function to guess column names based on common keywords
  guess_column <- function(cols_available, common_names) {
    for (name in common_names) {
      match_idx <- grep(paste0("^", name, "$"), cols_available, ignore.case = TRUE)
      if (length(match_idx) > 0) {
        return(cols_available[match_idx[1]])
      }
    }
    return("")
  }

  # Observer for file upload: reads the uploaded Excel file and updates column selectors
  observeEvent(input$data_file, {
    req(input$data_file)
    tryCatch({
      data <- read_excel(input$data_file$datapath)
      data_reactive(data)
      message_rv(list(type = "success", text = "Data file uploaded and loaded successfully."))

      col_names <- colnames(data)
      all_col_choices_with_none <- c("None" = "", col_names)

      # Automatically select likely columns
      updateSelectInput(session, "col_value", choices = all_col_choices_with_none, selected = guess_column(col_names, c("HB_value", "Value", "Result", "Measurement", "Waarde")))
      updateSelectInput(session, "col_age", choices = all_col_choices_with_none, selected = guess_column(col_names, c("leeftijd", "age", "AgeInYears", "Years")))
      updateSelectInput(session, "col_gender", choices = all_col_choices_with_none, selected = guess_column(col_names, c("geslacht", "gender", "sex", "Gender", "Sex")))
    }, error = function(e) {
      message_rv(list(type = "error", text = paste("Error loading file:", e$message)))
      data_reactive(NULL)
    })
  })

  # Observer for the Reset button on the Main Analysis tab
  observeEvent(input$reset_btn, {
    # Reset all inputs and outputs to their initial state
    shinyjs::reset("data_file")
    data_reactive(NULL)
    refiner_model_rv(NULL) # Reset the model
    plot_title_rv("") # Reset the title
    message_rv(list(type = "", text = ""))
    output$result_text <- renderPrint({ cat("") })
    output$result_plot <- renderPlot(plot.new())

    updateSelectInput(session, "col_value", choices = c("None" = ""), selected = "")
    updateSelectInput(session, "col_age", choices = c("None" = ""), selected = "")
    updateSelectInput(session, "col_gender", choices = c("None" = ""), selected = "")
  })

  # Observer for directory selection using shinyFiles
  shinyFiles::shinyDirChoose(
    input, id = 'select_dir_btn',
    roots = c(home = '~', wd = '.'), session = session
  )

  # Updates the reactive value with the selected directory path
  observeEvent(input$select_dir_btn, {
    if (!is.integer(input$select_dir_btn)) {
      path <- shinyFiles::parseDirPath(c(home = '~', wd = '.'), input$select_dir_btn)
      if (length(path) > 0) {
        selected_dir_reactive(path)
        message_rv(list(type = "success", text = paste("Output directory selected:", path)))
      } else {
        selected_dir_reactive(NULL)
        message_rv(list(type = "warning", text = "Directory selection cancelled."))
      }
    }
  })
  
  # Reactive expression for filtered data
  # This filters the raw data based on user selections (gender and age range)
  filtered_data_reactive <- reactive({
    req(data_reactive(), input$col_value, input$col_age, input$col_gender)
    if (input$col_value == "" || input$col_age == "" || input$col_gender == "") {
      return(NULL)
    }
    filter_data(data_reactive(), input$gender_choice, input$age_range[1], input$age_range[2], input$col_gender, input$col_age)
  })

  # Observer for the Analyze button
  # This is the core logic for the main analysis, running the refineR model
  observeEvent(input$analyze_btn, {
    if (analysis_running_rv()) {
      message_rv(list(text = "Analysis is already running. Please wait or reset.", type = "warning"))
      return()
    }

    filtered_data <- filtered_data_reactive()
    req(filtered_data)
    
    if (nrow(filtered_data) == 0) {
      message_rv(list(text = "Filtered dataset is empty. Please adjust your filtering criteria.", type = "danger"))
      return()
    }
    
    # Disable button and change text when analysis starts
    shinyjs::disable("analyze_btn")
    shinyjs::runjs("$('#analyze_btn').text('Analyzing...');")
    
    analysis_running_rv(TRUE)
    message_rv(list(text = "Analysis started...", type = "info"))
    session$sendCustomMessage('analysisStatus', TRUE)

    # Isolate inputs to prevent re-running the analysis on every change
    isolated_inputs <- isolate({
      list(
        gender_choice = input$gender_choice,
        age_range = input$age_range,
        col_value = input$col_value,
        col_age = input$col_age,
        col_gender = input$col_gender,
        nbootstrap_speed = input$nbootstrap_speed,
        unit_input = input$unit_input,
        ref_low = input$ref_low,
        ref_high = input$ref_high,
        enable_directory = input$enable_directory
      )
    })
    
    refiner_model <- NULL
    
    tryCatch({
      nbootstrap_value <- switch(isolated_inputs$nbootstrap_speed, "Fast" = 1, "Medium" = 50, "Slow" = 200, 1)

      # Run the main RefineR function
      refiner_model <- refineR::findRI(Data = filtered_data[[isolated_inputs$col_value]], NBootstrap = nbootstrap_value)
      
      if (is.null(refiner_model) || inherits(refiner_model, "try-error")) {
        stop("RefineR model could not be generated. Check your input data and parameters.")
      }
      
      refiner_model_rv(refiner_model)
      plot_title_rv(paste0("Estimated Reference Intervals for ", isolated_inputs$col_value, 
                           " (Gender: ", isolated_inputs$gender_choice, 
                           ", Age: ", isolated_inputs$age_range[1], "-", isolated_inputs$age_range[2], ")"))
      
      # Render the text summary of the model
      output$result_text <- renderPrint({
        print(refiner_model)
      })

      # If auto-save is enabled, save the plot to the selected directory
      if (isolated_inputs$enable_directory && !is.null(selected_dir_reactive())) {
        filename <- generate_safe_filename("RefineR_Plot", selected_dir_reactive(), "png")
        png(filename, width = 800, height = 600)
        
        plot(refiner_model, showCI = TRUE, RIperc = c(0.025, 0.975), showPathol = FALSE,
             title = plot_title_rv(),
             xlab = sprintf("%s [%s]", isolated_inputs$col_value, isolated_inputs$unit_input))

        usr <- par("usr")
        y_max <- usr[4]
        y_label_pos <- y_max * 0.95

        # Add manual reference limit lines
        if (!is.na(isolated_inputs$ref_low) && is.numeric(isolated_inputs$ref_low)) {
          abline(v = isolated_inputs$ref_low, col = "red", lty = 2, lwd = 2)
          text(x = isolated_inputs$ref_low, y = y_label_pos,
               labels = round(isolated_inputs$ref_low, 2),
               col = "red", cex = 1.1, pos = 4)
        }

        if (!is.na(isolated_inputs$ref_high) && is.numeric(isolated_inputs$ref_high)) {
          abline(v = isolated_inputs$ref_high, col = "blue", lty = 2, lwd = 2)
          text(x = isolated_inputs$ref_high, y = y_label_pos,
               labels = round(isolated_inputs$ref_high, 2),
               col = "blue", cex = 1.1, pos = 2)
        }
        
        dev.off()
        message_rv(list(text = paste0("Plot saved to ", selected_dir_reactive()), type = "success"))
      }

      message_rv(list(text = "Analysis complete!", type = "success"))

    }, error = function(e) {
      error_message <- paste("Analysis Error:", e$message)
      message_rv(list(text = error_message, type = "danger"))
      output$result_text <- renderPrint({ cat(error_message) })
      output$result_plot <- renderPlot(plot.new())
      print(error_message)
    }, finally = {
      # Re-enable button and restore text when analysis finishes
      analysis_running_rv(FALSE)
      session$sendCustomMessage('analysisStatus', FALSE)
      shinyjs::enable("analyze_btn")
      shinyjs::runjs("$('#analyze_btn').text('Analyze');")
    })
  })

  # Renders the live-updating plot output that depends on reactive inputs
  output$result_plot <- renderPlot({
    refiner_model <- refiner_model_rv()
    req(refiner_model) # Requires the model to be present before plotting
    plot_title <- plot_title_rv()
    
    plot(refiner_model, showCI = TRUE, RIperc = c(0.025, 0.975), showPathol = FALSE,
         title = plot_title,
         xlab = sprintf("%s [%s]", input$col_value, input$unit_input))

    usr <- par("usr")
    y_max <- usr[4]
    y_label_pos <- y_max * 0.95

    # Live update manual lower limit line and text
    if (!is.na(input$ref_low) && is.numeric(input$ref_low)) {
      abline(v = input$ref_low, col = "red", lty = 2, lwd = 2)
      text(x = input$ref_low, y = y_label_pos,
           labels = round(input$ref_low, 2),
           col = "red", cex = 1.1, pos = 4)
    }

    # Live update manual upper limit line and text
    if (!is.na(input$ref_high) && is.numeric(input$ref_high)) {
      abline(v = input$ref_high, col = "blue", lty = 2, lwd = 2)
      text(x = input$ref_high, y = y_label_pos,
           labels = round(input$ref_high, 2),
           col = "blue", cex = 1.1, pos = 2)
    }
  })
}