handleMatrixCreation <- function(input, rv) {
  req(rv$processed_file)

  # Check if organism ID is specified
  organism_id <- get_organism_id(input$select_organism, input$text_select_organism)
  if (is.null(organism_id)) {
    showNotification("Please select or enter an organism ID.", type = "error")
    return()
  }

  shinyjs::show("spinner")
  rv$download_message <- "Please wait while the matrix is computed..."
  rv$processing_complete <- FALSE
  rv$arena3d_link <- NULL
  rv$enzyme_table <- NULL  # Reset the enzyme table
  rv$enzyme_table_loaded <- FALSE  # Reset the enzyme table loaded state

  src_path <- here("src", "core")
  print(src_path)
  py <- reticulate::import_from_path("lipid_analysis_workflow", src_path)
  workflow <- py$LipidAnalysisWorkflow(here("config", "lipid_config.yaml"), tempdir())

  tryCatch({
    workflow$process_files(
      list(rv$processed_file$filepath),
      organism_id,
      TRUE,
      input$processNetwork
    )
    rv$download_message <- "Analysis complete. You can now download the results."
    rv$processing_complete <- TRUE

    # Check if network JSON file exists
    network_file <- file.path(tempdir(), paste0(rv$processed_file$timestamp, "_network.json"))
    if (input$processNetwork && file.exists(network_file)) {
      arena3d_url <- sendNetworkToAPI(network_file)
      if (!is.null(arena3d_url)) {
        rv$arena3d_link <- arena3d_url
      } else {
        showNotification("Failed to get Arena3D URL from API", type = "warning")
      }
    }
  }, error = function(e) {
    rv$download_message <- paste("Error:", e$message)
    rv$processing_complete <- FALSE
  })
  shinyjs::hide("spinner")
}


get_organism_id <- function(ui_organism_id, text_organism_id) {
  # Handle separator option
  if (!is.null(ui_organism_id) && ui_organism_id == "separator") {
    ui_organism_id <- " -- "
  }

  if (ui_organism_id == " -- " || is.null(ui_organism_id)) {
    if (is.null(text_organism_id) || text_organism_id == "") {
      return(NULL)
    }
    return(as.numeric(text_organism_id))
  } else {
    return(as.numeric(ui_organism_id))
  }
}
