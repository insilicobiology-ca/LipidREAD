handleSingleLipidSearch <- function(input, rv) {
  lipids_raw <- input$lipid_input_single
  if (is.null(lipids_raw) || lipids_raw == "") {
    showNotification("Please enter at least one lipid name.", type = "error")
    return()
  }
  lipids_list <- trimws(strsplit(lipids_raw, "\n")[[1]])
  lipids_list <- lipids_list[lipids_list != ""]

  organism_id <- get_organism_id(input$select_organism_lipid, input$text_select_organism_lipid)
  if (is.null(organism_id)) {
    showNotification("Please select or enter an organism ID for the lipid search.", type = "error")
    return()
  }

  shinyjs::show("spinner")
  rv$single_lipid_results_df <- NULL
  rv$single_lipid_timestamp <- NULL

  tryCatch({
    py <- reticulate::import_from_path("lipid_analysis_workflow", here("src", "core"))
    workflow <- py$LipidAnalysisWorkflow(here("config", "lipid_config.yaml"), tempdir())

    timestamp <- format(Sys.time(), "%m%d%Y%H%M%S")
    base_name <- paste0(timestamp, "_single_lipid_analysis")

    workflow$analyze_single_lipid(lipids_list, organism_id, base_name)

    report_filename <- paste0(base_name, "_single_lipid_search_report.tsv")
    report_path <- file.path(tempdir(), report_filename)

    if (file.exists(report_path)) {
      rv$single_lipid_results_df <- read.delim(report_path, sep = "\t", check.names = FALSE)
      rv$single_lipid_timestamp <- timestamp
    } else {
      showNotification("Analysis completed, but the output report file was not found.", type = "warning")
    }

  }, error = function(e) {
    showNotification(paste("An error occurred during analysis:", e$message), type = "error")
  })

  shinyjs::hide("spinner")
}


handleSingleEnzymeSearch <- function(input, rv) {
  enzyme_id <- input$enzyme_identifier_single
  if (is.null(enzyme_id) || enzyme_id == "") {
    showNotification("Please enter an enzyme identifier.", type = "error")
    return()
  }

  organism_id <- get_organism_id(input$select_organism_enzyme, input$text_select_organism_enzyme)
  if (is.null(organism_id)) {
    showNotification("Please select or enter an organism ID for the enzyme search.", type = "error")
    return()
  }

  shinyjs::show("spinner")
  rv$single_enzyme_results_df <- NULL
  rv$single_enzyme_timestamp <- NULL

  tryCatch({
    py <- reticulate::import_from_path("lipid_analysis_workflow", here("src", "core"))
    workflow <- py$LipidAnalysisWorkflow(here("config", "lipid_config.yaml"), tempdir())

    timestamp <- format(Sys.time(), "%m%d%Y%H%M%S")
    base_name <- paste0(timestamp, "_single_enzyme_analysis")

    workflow$search_reactions_by_enzyme(enzyme_id, organism_id, base_name)

    report_filename <- paste0(base_name, "_single_enzyme_search_report.tsv")
    report_path <- file.path(tempdir(), report_filename)

    if (file.exists(report_path)) {
      rv$single_enzyme_results_df <- read.delim(report_path, sep = "\t", check.names = FALSE)
      rv$single_enzyme_timestamp <- timestamp
    } else {
      # Handle case where Python finds no results and doesn't create a file
      py_results <- workflow$search_reactions_by_enzyme(enzyme_id, organism_id, base_name)
      if (!py_results$found) {
        showNotification("Enzyme identifier not found.", type = "warning")
      } else {
        showNotification("Analysis completed, but the output report file was not found.", type = "warning")
      }
    }

  }, error = function(e) {
    showNotification(paste("An error occurred during analysis:", e$message), type = "error")
  })

  shinyjs::hide("spinner")
}