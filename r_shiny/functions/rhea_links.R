# Rhea database link generation functions
# This file handles the creation of clickable links to Rhea database entries

#' Create clickable Rhea database links
#' @param rhea_id Rhea ID value (can be NA, empty, or numeric)
#' @return HTML string with clickable link or dash for missing data
create_rhea_links <- function(rhea_id) {
  # Handle various empty/null cases
  if (is.na(rhea_id) || rhea_id == "" || rhea_id == "NULL" || rhea_id == "NA") {
    return("-")  # Em dash for missing data
  }
  
  # Clean the Rhea ID (remove any extra whitespace)
  rhea_id <- trimws(as.character(rhea_id))
  
  # Basic validation - Rhea IDs are typically numeric
  if (nchar(rhea_id) > 0 && grepl("^[0-9]+$", rhea_id)) {
    paste0('<a href="https://www.rhea-db.org/rhea/',
           rhea_id,
           '" target="_blank" ',
           'title="View reaction ', rhea_id, ' on Rhea database" ',
           'style="color: #007bff; text-decoration: none;">',
           rhea_id,
           '</a>')
  } else {
    # If it doesn't look like a Rhea ID, just return as text
    as.character(rhea_id)
  }
}

#' Process dataframe columns that should contain Rhea ID links
#' @param df Data frame to process
#' @return Data frame with Rhea ID columns converted to HTML links
process_specific_rhea_columns <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  
  # Target columns that should have Rhea ID links
  # Add any additional column name patterns as needed
  target_columns <- c(
    "Rhea ID",
    "Sourced from Rhea Example",
    "rhea_id",
    "Rhea_ID",
    "RheaID",
    "Rhea.ID"  # In case dots are used instead of spaces
  )
  
  # Apply Rhea link creation to matching columns
  for (col_name in target_columns) {
    if (col_name %in% colnames(df)) {
      df[[col_name]] <- sapply(df[[col_name]], create_rhea_links, USE.NAMES = FALSE)
    }
  }
  
  return(df)
}

#' Process all Rhea-related columns (more flexible version)
#' @param df Data frame to process
#' @return Data frame with all potential Rhea columns converted to HTML links
process_all_rhea_columns <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  
  # Find columns that likely contain Rhea IDs based on name patterns
  rhea_column_pattern <- "rhea|Rhea|RHEA"
  rhea_columns <- grep(rhea_column_pattern, colnames(df), value = TRUE, ignore.case = TRUE)
  
  # Apply Rhea link creation to all identified columns
  for (col_name in rhea_columns) {
    df[[col_name]] <- sapply(df[[col_name]], create_rhea_links, USE.NAMES = FALSE)
  }
  
  return(df)
}