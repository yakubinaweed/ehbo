library(shiny)
library(bslib)
library(refineR)
library(readxl)
library(moments)
library(shinyjs)
library(shinyWidgets)
library(shinyFiles)
library(ggplot2)

ui <- navbarPage(
  title = "RefineR Reference Interval Estimation",
  id = "tabs",
  # Sets the visual theme and fonts for the Shiny app
  theme = bs_theme(version = 4, base_font = font_google("Inter"), heading_font = font_google("Rethink Sans"), font_scale = 1.1, bootswatch = "default"),

  # First tab for the main RefineR analysis
  tabPanel(
    title = "Main Analysis",
    useShinyjs(),
    tags$head(
      # Includes the custom CSS from the 'www' directory
      includeCSS("www/styles.css")
    ),
    # Custom JavaScript to handle disabling the other tab during analysis
    tags$script(HTML("
      var analysisRunning = false;
      Shiny.addCustomMessageHandler('analysisStatus', function(status) {
        analysisRunning = status;
        if (status) {
          $('a[data-value=\"Subpopulation Detection (GMM)\"]').attr('data-toggle', 'disabled').addClass('disabled-tab-link');
        } else {
          $('a[data-value=\"Subpopulation Detection (GMM)\"]').attr('data-toggle', 'tab').removeClass('disabled-tab-link');
        }
      });
      $(document).on('click', 'a[data-toggle=\"tab\"]', function(event) {
        if (analysisRunning) {
          event.preventDefault();
          Shiny.setInputValue('tab_switch_blocked', new Date().getTime());
          return false;
        }
      });
    ")),
    sidebarLayout(
      sidebarPanel(
        style = "padding-right: 15px;",
        # User inputs for data filtering and analysis parameters
        selectInput(inputId = "gender_choice", label = "Select Gender:", choices = c("Male" = "M", "Female" = "F", "Both" = "Both"), selected = "Both"),
        sliderInput(inputId = "age_range", label = "Age Range:", min = 0, max = 100, value = c(0, 100), step = 1),
        fileInput(inputId = "data_file", label = "Upload Data (Excel File)", accept = c(".xlsx")),
        # Dynamic inputs for selecting data columns
        selectInput(inputId = "col_value", label = "Select Column for Values:", choices = c("None" = ""), selected = ""),
        selectInput(inputId = "col_age", label = "Select Column for Age:", choices = c("None" = ""), selected = ""),
        selectInput(inputId = "col_gender", label = "Select Column for Gender:", choices = c("None" = ""), selected = ""),
        radioButtons(inputId = "nbootstrap_speed", label = "Select Computation Speed:", choices = c("Fast", "Medium", "Slow"), selected = "Fast", inline = TRUE),
        # Action buttons for the analysis
        actionButton("analyze_btn", "Analyze", class = "btn-primary"),
        actionButton("reset_btn", "Reset File", class = "btn-secondary"),
        shinyFiles::shinyDirButton(id = "select_dir_btn", label = "Select Output Directory", title = "Select a directory to save plots", style = "margin-top: 5px;"),
        div(style = "margin-top: 5px; display: flex; align-items: center; justify-content: flex-start; width: 100%;",
            prettySwitch(inputId = "enable_directory", label = "Auto-Save Graph", status = "success", fill = TRUE, inline = TRUE)
        ),
        uiOutput("app_message"), # Placeholder for displaying app messages
        hr(),
        # Inputs for manual reference limits and units for the plot
        numericInput("ref_low", "Reference Lower Limit:", value = NA),
        numericInput("ref_high", "Reference Upper Limit:", value = NA),
        textInput(inputId = "unit_input", label = "Unit of Measurement", value = "mmol/L", placeholder = "ex. g/L")
      ),
      mainPanel(
        # Outputs for the main analysis results
        plotOutput("result_plot"),
        verbatimTextOutput("result_text")
      )
    )
  ),

  # Second tab for Gaussian Mixture Model (GMM) subpopulation detection
  tabPanel(
    title = "Subpopulation Detection (GMM)",
    useShinyjs(),
    h4("Detect Subpopulations using HGB and Age"),
    p("Gaussian Mixture Models aim to detect hidden subpopulations within your data based on HGB and Age. The system will automatically select the best model and number of components based on the BIC criterion. For each detected subpopulation, estimated age ranges will be provided directly from the model's characteristics, avoiding predefined bins. While increasing the number of components can improve model fit, it also increases the risk of overfitting, where the model learns noise rather than true underlying patterns."),
    sidebarLayout(
      sidebarPanel(
        fileInput(inputId = "gmm_file_upload", label = "Upload Data (Excel File)", accept = c(".xlsx")),
        hr(),
        # New radio buttons for gender selection
        radioButtons(inputId = "gmm_gender_choice", label = "Select Gender:", choices = c("Male" = "Male", "Female" = "Female", "Both" = "Both"), selected = "Both", inline = TRUE),
        # Dynamic inputs for selecting HGB, Age, and Gender columns for GMM
        selectInput(inputId = "gmm_hgb_col", label = "Select HGB Column:", choices = c("None" = ""), selected = ""),
        selectInput(inputId = "gmm_age_col", label = "Select Age Column:", choices = c("None" = ""), selected = ""),
        selectInput(inputId = "gmm_gender_col", label = "Select Gender Column:", choices = c("None" = ""), selected = ""),
        hr(),
        # Action buttons for the GMM analysis
        actionButton("run_gmm_analysis_btn", "Run Subpopulation Detection", class = "btn-primary"),
        actionButton("reset_gmm_analysis_btn", "Reset GMM Data", class = "btn-secondary")
      ),
      mainPanel(
        # Renders the UI for GMM results dynamically
        uiOutput("gmm_results_ui")
      )
    )
  ),

  # Footer of the application with copyright and a link to the author's GitHub
  footer = tags$footer(
    HTML('© 2025 <a href="https://github.com/yakubinaweed/refineR-reference-interval" target="_blank">Naweed Yakubi</a> • All rights reserved.'),
    style = "
      position: relative;
      bottom: 0;
      width: 100%;
      text-align: center;
      padding: 10px 0;
      color: #777;
      font-size: 0.8em;"
  )
)