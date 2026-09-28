handleFileUpload <- function(input, rv) {
  req(input$input_lipid_data)
  processed_file <- process_input_file(input$input_lipid_data, tempdir())
  rv$processed_file <- processed_file
  rv$processing_complete <- FALSE
  rv$download_message <- "Progress indicator - This will be greyed out when analysis is running."
  print(paste("Processed file:", processed_file$filepath))  # Add this for debugging
  print(paste("rv$processed_file:", rv$processed_file$filepath))  # Additional debug line
}

process_input_file <- function(input_file, temp_dir) {
  timestamp <- format(Sys.time(), "%m%d%Y%H%M%S")
  new_filename <- paste0(timestamp, ".", tools::file_ext(input_file$name))
  new_filepath <- file.path(temp_dir, new_filename)
  file.copy(input_file$datapath, new_filepath, overwrite = TRUE)
  print(new_filepath)
  list(timestamp = timestamp, filepath = new_filepath)

}