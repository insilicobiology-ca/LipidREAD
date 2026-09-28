get_all_organisms <- function() {
  tryCatch({
    py <- reticulate::import_from_path("lipid_analysis_workflow", here("src", "core"))
    workflow <- py$LipidAnalysisWorkflow(here("config", "lipid_config.yaml"), tempdir())

    # FIXED: Use correct column name 'taxon_scientific_name'
    query <- paste(
      "SELECT organism_id, taxon_scientific_name",
      "FROM lipograph.organism",
      "WHERE organism_id IS NOT NULL",
      "AND taxon_scientific_name IS NOT NULL",
      "AND taxon_scientific_name != ''",
      "ORDER BY taxon_scientific_name"
    )

    workflow$data_access$cursor$execute(query)
    results <- workflow$data_access$cursor$fetchall()

    if (length(results) == 0) {
      cat("No organisms found in database\n")
      return(NULL)
    }

    # Build choices vector safely
    organism_choices <- c()

    for (result in results) {
      if (!is.null(result$organism_id) && !is.null(result$taxon_scientific_name)) {
        name <- as.character(result$taxon_scientific_name)
        id <- as.numeric(result$organism_id)

        if (!is.na(name) && !is.na(id) && nchar(name) > 0) {
          organism_choices[name] <- id
        }
      }
    }

    if (length(organism_choices) > 0) {
      cat("Successfully loaded", length(organism_choices), "organisms\n")
      return(organism_choices)
    } else {
      cat("No valid organisms after filtering\n")
      return(NULL)
    }

  }, error = function(e) {
    cat("Error loading organisms:", e$message, "\n")
    return(NULL)
  })
}

# Create the dropdown choices with priority organisms at top
create_organism_choices <- function(all_organisms) {
  if (is.null(all_organisms)) {
    # Fallback to hardcoded if database query fails
    return(c(
      " -- " = " -- ",
      "Homo sapiens" = 9606,
      "Mus musculus" = 10090,
      "Rattus norvegicus" = 10116
    ))
  }

  # Priority organisms
  priority_organisms <- c(
    "Homo sapiens" = 9606,
    "Mus musculus" = 10090,
    "Rattus norvegicus" = 10116
  )

  # Remove priority organisms from full list to avoid duplicates
  remaining_organisms <- all_organisms[!names(all_organisms) %in% names(priority_organisms)]

  # Combine: placeholder, priority organisms, separator, then all others
  choices <- c(
    " -- " = " -- ",
    priority_organisms,
    "-----------------" = "separator",
    remaining_organisms
  )

  return(choices)
}

# # FIXED: Match your original working structure
# load_organism_choices <- function() {
#   # Step 1: Get all organisms from database (just like your ALL_ORGANISMS)
#   all_organisms <- get_all_organisms()
#
#   # Step 2: Create choices with priority and separator (just like your ORGANISM_CHOICES)
#   organism_choices <- create_organism_choices(all_organisms)
#
#   return(organism_choices)
# }