handleDownload <- function(file, input, rv) {
  req(rv$processed_file)

  files_to_zip <- c(
    paste0(rv$processed_file$timestamp, "_enzyme_adjacency_matrix.csv"),
    paste0(rv$processed_file$timestamp, "_binary_adjacency_matrix.csv"),
    paste0(rv$processed_file$timestamp, "_full_reaction_list.tsv"),
    paste0(rv$processed_file$timestamp, "_translation_tracker.tsv"),
    paste0(rv$processed_file$timestamp, "_pair_list.tsv"),
    paste0(rv$processed_file$timestamp, "_enzyme_table.csv")
  )

  # if (input$translate) {
  #   files_to_zip <- c(files_to_zip,
  #                     paste0(rv$processed_file$timestamp, "_translation_tracker.tsv"))
  # }

  zip::zip(file, files = files_to_zip, root = tempdir(), mode = "cherry-pick")
}

handleSingleDownload <- function(file, file_path) {
  # Check if the file actually exists before trying to copy it.
  if (file.exists(file_path)) {
    file.copy(file_path, file)
  } else {
    # This is a fallback in case the file path is set but the file is missing.
    write("Error: The report file could not be found on the server.", file)
  }
}