# CompLi: Computational Lipidomics
# Shiny app
# By Qassim Alkassir


library(shiny)
library(shinyjs)
library(shinydashboard)
library(zip)
library(reticulate)
library(here)
library(shinybusy)
library(dplyr)
library(DT)
library(httr)

# Set Python environment
#Sys.setenv(RETICULATE_PYTHON = "C:/Users/Miroslava/anaconda3/envs/py311_test/python.exe")
#use_condaenv("C:/Users/Miroslava/anaconda3/envs/py311_test/python.exe")
use_python("/home/almalinux/condaenv/compli/bin/python", required  = TRUE)

# Source all required functions and modules
source("r_shiny/modules/ui_components.R")
source("r_shiny/functions/organism_data.R")
source("r_shiny/functions/rhea_links.R")
source("r_shiny/functions/table_renderers.R")
source("r_shiny/functions/input.R", local = TRUE)
source("r_shiny/functions/processing.R", local = TRUE)
source("r_shiny/functions/download.R", local = TRUE)
source("r_shiny/functions/enzyme_table.R", local = TRUE)
source("r_shiny/functions/sendNetworkToAPI.R", local = TRUE)
source("r_shiny/functions/single_search.R", local = TRUE)


valid_users <- data.frame(
  username = c("nrluser"),
  password = c("E+fv29Xy"),
  stringsAsFactors = FALSE
)

`%||%` <- function(a, b) if (is.null(a)) b else a



# Cache organisms at app startup - BEFORE UI definition
ALL_ORGANISMS <- get_all_organisms()
ORGANISM_CHOICES <- create_organism_choices(ALL_ORGANISMS)

# UI Definition
ui <- dashboardPage(
  dashboardHeader(
    title = tags$a(
      href = "https://nsilicobiology.ca",
      tags$span("In Silico Biology", style = "color: white; text-decoration: none;"),
      target = "_blank",
      style = "text-decoration: none;"
    )
  ),
  sidebar(),
  body(ORGANISM_CHOICES),
  skin = "black",
  tags$head(
    includeCSS(here("www", "complimetGUI.css"))
  )
)

# Server Logic
server <- function(input, output, session) {





user_authenticated <- reactiveVal(FALSE)

  show_login_modal <- function(error_msg = NULL) {
    showModal(modalDialog(
      title = "MWMVisual v1.0 - Sign In",
      textInput("login_username", "Username"),
      passwordInput("login_password", "Password"),
      if (!is.null(error_msg)) {
        tags$div(style = "color:#d9534f; font-weight:bold;", error_msg)
      },
      footer = actionButton("login_btn", "Sign In", class = "btn-primary"),
      easyClose = FALSE
    ))
  }

  observe({
    if (!isTRUE(user_authenticated())) show_login_modal()
  })

  observeEvent(input$login_btn, {
    user <- trimws(input$login_username %||% "")
    pass <- input$login_password %||% ""
    if (length(which(valid_users$username == user & valid_users$password == pass)) == 1) {
      user_authenticated(TRUE)
      removeModal()
    } else {
      removeModal()
      show_login_modal("Invalid username or password.")
    }
  })


  # Initialize reactive values
  rv <- reactiveValues(
    download_message = "Progress indicator - This will be greyed out when analysis is running.",
    processed_file = NULL,
    processing_complete = FALSE,
    enzyme_table = NULL,
    enzyme_table_loaded = FALSE,
    arena3d_link = NULL,
    single_lipid_results_df = NULL,
    single_lipid_timestamp = NULL,
    single_enzyme_results_df = NULL,
    single_enzyme_timestamp = NULL
  )

  # ============================================================================
  # EVENT OBSERVERS
  # ============================================================================

  # File upload handler
  observeEvent(input$input_lipid_data, {
    handleFileUpload(input, rv)
  }, ignoreInit = TRUE)

  # Matrix creation handler
  observeEvent(input$createAdjMatrix, {
    print("createAdjMatrix button clicked")
    print(paste("Before handleMatrixCreation, rv$processed_file:", rv$processed_file$filepath))
    handleMatrixCreation(input, rv)
  }, ignoreInit = TRUE)

  # Single lipid search handler
  observeEvent(input$runSingleLipidSearch, {
    handleSingleLipidSearch(input, rv)
  })

  # Single enzyme search handler
  observeEvent(input$runSingleEnzymeSearch, {
    handleSingleEnzymeSearch(input, rv)
  })

  # Load enzyme table after processing is complete
  observe({
    req(rv$processing_complete)
    loadEnzymeTable(rv)
  })

  # ============================================================================
  # RENDER OUTPUTS - TABLES
  # ============================================================================

  # Main enzyme table
  output$enzymeTable <- DT::renderDataTable({
    renderEnzymeTable(rv$enzyme_table, rv$enzyme_table_loaded)
  })

  # Single lipid results table
  output$singleLipidResultsTable <- DT::renderDataTable({
    req(rv$single_lipid_results_df)
    render_single_search_table(
      rv$single_lipid_results_df,
      "Single lipid search results. Rhea IDs link to detailed reaction information."
    )
  })

  # Single enzyme results table
  output$singleEnzymeResultsTable <- DT::renderDataTable({
    req(rv$single_enzyme_results_df)
    render_single_search_table(
      rv$single_enzyme_results_df,
      "Single enzyme search results. Rhea IDs link to detailed reaction information."
    )
  })

  # Supported classes table (static)
  output$supportedClassesTable <- DT::renderDataTable({
    create_supported_classes_table()
  }, server = FALSE)

  # ============================================================================
  # RENDER OUTPUTS - UI ELEMENTS
  # ============================================================================

  # Enzyme table message
  output$enzymeTableMessage <- renderUI({
    if (!rv$enzyme_table_loaded) {
      return(HTML("<p>Enzyme table will appear here after processing.</p>"))
    } else if (is.null(rv$enzyme_table)) {
      return(HTML("<p>No enzyme data available.</p>"))
    } else {
      return(NULL)
    }
  })

  # Arena3D link
  output$arena3dLink <- renderUI({
    renderArena3dLink(rv$arena3d_link)
  })

  # Download message
  output$downloadMessage <- renderText({
    rv$download_message
  })

  # Download button UIs for single searches
  output$downloadSingleLipidUI <- renderUI({
    req(rv$single_lipid_timestamp)
    downloadButton("downloadSingleLipidResults", "Download Lipid Results", class = "btn-orange")
  })

  output$downloadSingleEnzymeUI <- renderUI({
    req(rv$single_enzyme_timestamp)
    downloadButton("downloadSingleEnzymeResults", "Download Enzyme Results", class = "btn-orange")
  })

  # ============================================================================
  # DOWNLOAD HANDLERS
  # ============================================================================

  # Main matrix download
  output$downloadMatrix <- downloadHandler(
    filename = "LipidREAD_files.zip",
    content = function(file) {
      handleDownload(file, input, rv)
    }
  )

  # Enzyme table download
  output$downloadEnzymeTable <- downloadHandler(
    filename = function() {
      paste0("enzyme_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      downloadEnzymeTable(file, rv)
    }
  )

  # Single lipid results download
  output$downloadSingleLipidResults <- downloadHandler(
    filename = function() {
      req(rv$single_lipid_timestamp)
      paste0(rv$single_lipid_timestamp, "_single_lipid_analysis_single_lipid_search_report.tsv")
    },
    content = function(file) {
      req(rv$single_lipid_timestamp)
      report_filename <- paste0(rv$single_lipid_timestamp, "_single_lipid_analysis_single_lipid_search_report.tsv")
      report_path <- file.path(tempdir(), report_filename)
      handleSingleDownload(file, report_path)
    }
  )

  # Single enzyme results download
  output$downloadSingleEnzymeResults <- downloadHandler(
    filename = function() {
      req(rv$single_enzyme_timestamp)
      paste0(rv$single_enzyme_timestamp, "_single_enzyme_analysis_single_enzyme_search_report.tsv")
    },
    content = function(file) {
      req(rv$single_enzyme_timestamp)
      report_filename <- paste0(rv$single_enzyme_timestamp, "_single_enzyme_analysis_single_enzyme_search_report.tsv")
      report_path <- file.path(tempdir(), report_filename)
      handleSingleDownload(file, report_path)
    }
  )

  # Sample data download
  output$downloadBatchSample <- downloadHandler(
    filename = function() {
      "LipidREAD_Sample_Data.xlsx"
    },
    content = function(file) {
      file.copy("www/LipidREAD_Sample_Data.xlsx", file)
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )

  # ============================================================================
  # REACTIVE UI UPDATES
  # ============================================================================

  # Show/hide download buttons based on completion status
  observe({
    shinyjs::toggle("downloadMatrix", condition = rv$processing_complete)
    shinyjs::toggle("downloadEnzymeTable", condition = rv$enzyme_table_loaded)
  })



}



# Run the application
shinyApp(ui = ui, server = server)
