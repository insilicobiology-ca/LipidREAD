loadEnzymeTable <- function(rv) {
  req(rv$processed_file)

  enzyme_file <- paste0(rv$processed_file$timestamp, "_enzyme_table.csv")
  enzyme_path <- file.path(tempdir(), enzyme_file)

  if (file.exists(enzyme_path)) {
    enzyme_data <- read.csv(enzyme_path, stringsAsFactors = FALSE, check.names = FALSE)
    rv$enzyme_table <- enzyme_data
    rv$enzyme_table_loaded <- TRUE
  } else {
    rv$enzyme_table <- NULL
    rv$enzyme_table_loaded <- FALSE
  }
}

# Enhanced version with better styling and error handling
renderEnzymeTable <- function(enzyme_table, enzyme_table_loaded) {
  if (!enzyme_table_loaded || is.null(enzyme_table)) {
    return(NULL)
  } else {

    # Create a copy of the table to modify
    display_table <- enzyme_table

    # Find UniProt ID column (flexible column name matching)
    uniprot_columns <- grep("uniprot|UniProt|UNIPROT", colnames(enzyme_table),
                           value = TRUE, ignore.case = TRUE)

    # Process each UniProt column found
    for (col_name in uniprot_columns) {
      display_table[[col_name]] <- sapply(enzyme_table[[col_name]], function(uniprot_id) {
        # Handle various empty/null cases
        if (is.na(uniprot_id) || uniprot_id == "" ||
            uniprot_id == "NULL" || uniprot_id == "NA") {
          return("—")  # Em dash for missing data
        }

        # Handle multiple UniProt IDs separated by various delimiters
        ids <- strsplit(as.character(uniprot_id), "[;,|\\s]+")[[1]]
        ids <- trimws(ids)  # Remove whitespace
        ids <- ids[nchar(ids) > 0]  # Remove empty strings

        if (length(ids) == 0) {
          return("—")
        }

        # Create links for each ID
        links <- sapply(ids, function(id) {
          # Basic validation - UniProt IDs are typically 6-10 characters
          if (nchar(id) >= 4 && nchar(id) <= 15) {
            paste0('<a href="https://www.uniprot.org/uniprotkb/',
                   id,
                   '" target="_blank" ',
                   'title="View ', id, ' on UniProt" ',
                   'style="color: #007bff; text-decoration: none; margin-right: 8px;">',
                   id,
                   '</a>')
          } else {
            # If it doesn't look like a UniProt ID, just return as text
            paste0('<span style="margin-right: 8px;">', id, '</span>')
          }
        })

        # Join multiple links
        paste(links, collapse = "")
      })
    }

    # Create the datatable
    dt <- DT::datatable(display_table,
                       options = list(
                         pageLength = 10,
                         scrollX = TRUE,
                         scrollY = "300px",
                         columnDefs = list(
                           list(targets = "_all", className = "dt-center")
                         ),
                         # Add some styling
                         initComplete = DT::JS(
                           "function(settings, json) {",
                           "$(this.api().table().header()).css({'background-color': '#f8f9fa', 'color': '#495057'});",
                           "}"
                         )
                       ),
                       rownames = FALSE,
                       # Enable HTML rendering
                       escape = FALSE,
                       # Add caption
                       caption = paste("Enzymes found in your dataset. UniProt IDs link to uniprot.org"))

    return(dt)
  }
}

downloadEnzymeTable <- function(file, rv) {
  req(rv$enzyme_table_loaded)

  if (!is.null(rv$enzyme_table)) {
    write.csv(rv$enzyme_table, file, row.names = FALSE)
  }
}
