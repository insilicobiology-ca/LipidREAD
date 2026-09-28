render_single_search_table <- function(df, caption_text) {
  req(df)

  # Process the dataframe to add Rhea ID links
  display_df <- process_specific_rhea_columns(df)

  DT::datatable(
    display_df,
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      columnDefs = list(
        list(targets = "_all", className = "dt-center")
      ),
      # Add styling to header
      initComplete = DT::JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#f8f9fa', 'color': '#495057'});",
        "}"
      )
    ),
    rownames = FALSE,
    escape = FALSE,  # Critical for HTML links to render
    caption = caption_text
  )
}

#' Create the supported lipid classes table from TSV file
#' @param tsv_path Path to the supported lipid classes TSV file
#' @return DT datatable object with lipid class information
create_supported_classes_table <- function(tsv_path = "www/supported_lipid_classes.tsv") {

  # Try to read the TSV file
  tryCatch({
    lipid_classes_data <- read.delim(tsv_path, sep = "\t", stringsAsFactors = FALSE, encoding = "UTF-8")

    # Validate required columns
    required_cols <- c("Category", "LipidMAPS.Sub.class", "Examples")
    if (!all(required_cols %in% colnames(lipid_classes_data))) {
      warning("TSV file missing required columns. Using fallback data.")
      return(create_fallback_supported_classes_table())
    }

    # Clean up column names for display
    colnames(lipid_classes_data) <- c("Category", "LipidMAPS Sub-class", "Examples")

  }, error = function(e) {
    warning(paste("Could not read TSV file:", e$message, ". Using fallback data."))
    return(create_fallback_supported_classes_table())
  })

  # Create and style the datatable
  dt <- DT::datatable(
    lipid_classes_data,
    options = list(
      pageLength = 25,
      scrollX = TRUE,
      scrollY = "400px",
      dom = 'ftp',  # Include filter, table, and pagination
      columnDefs = list(
        list(width = '150px', targets = 0),
        list(width = '300px', targets = 1),
        list(width = '300px', targets = 2)
      ),
      # Add search functionality
      search = list(regex = TRUE, caseInsensitive = TRUE)
    ),
    rownames = FALSE,
    escape = FALSE,  # Allow HTML rendering for links
    caption = "Comprehensive list of lipid classes supported by LipidCRED with LipidMAPS classification",
    filter = "top"  # Add column filters at the top
  ) %>%
    DT::formatStyle(
      'Category',
      target = 'row',
      backgroundColor = DT::styleEqual(
        unique(lipid_classes_data$Category),
        generate_category_colors(unique(lipid_classes_data$Category))
      )
    )

  return(dt)
}

#' Generate colors for different categories
#' @param categories Vector of category names
#' @return Vector of color codes
generate_category_colors <- function(categories) {
  # Define a color palette
  color_palette <- c(
    '#f8f9fa',  # Light gray
    '#e3f2fd',  # Light blue
    '#f3e5f5',  # Light purple
    '#fff3e0',  # Light orange
    '#e8f5e8',  # Light green
    '#fff8e1',  # Light yellow
    '#fce4ec',  # Light pink
    '#e0f2f1'   # Light teal
  )

  # Assign colors cyclically
  colors <- color_palette[1:length(categories)]
  if (length(categories) > length(color_palette)) {
    colors <- rep(color_palette, length.out = length(categories))
  }

  return(colors)
}

#' Fallback function that creates the table with hardcoded data
#' @return DT datatable object with fallback lipid class information
create_fallback_supported_classes_table <- function() {

  # Fallback data in case TSV file is not available
  lipid_classes_data <- data.frame(
    Category = c(
      rep("Sphingolipids", 12),
      rep("Glycerophospholipids", 15),
      rep("Glycerolipids", 6),
      rep("Sterol Lipids", 4),
      rep("Prenol Lipids", 3)
    ),
    `LipidMAPS Sub-class` = c(
      # Sphingolipids
      "Ceramides (Cer)", "Sphingomyelins (SM)", "Glucosylceramides (beta-GlcCer)",
      "Galactosylceramides (beta-GalCer)", "Lactosylceramides (Lac-Cer)",
      "Sphingosines (Sph)", "Sphingosine-1-phosphates (S1P)",
      "Ceramide-1-phosphates (C1P)", "Gangliosides (GM1, GM2, GM3)",
      "Sulfatides (ST)", "Psychosines (Psy)", "Dihydroceramides (dhCer)",

      # Glycerophospholipids
      "Phosphatidylcholines (PC)", "Phosphatidylethanolamines (PE)",
      "Phosphatidylserines (PS)", "Phosphatidylglycerols (PG)",
      "Phosphatidylinositols (PI)", "Phosphatidic acids (PA)",
      "Lysophosphatidylcholines (LPC)", "Lysophosphatidylethanolamines (LPE)",
      "Lysophosphatidylserines (LPS)", "Lysophosphatidylglycerols (LPG)",
      "Lysophosphatidylinositols (LPI)", "Lysophosphatidic acids (LPA)",
      "Cardiolipins (CL)", "Phosphatidylinositol phosphates (PIP, PIP2)",
      "Ether-linked phospholipids (PE-O, PC-O)",

      # Glycerolipids
      "Triacylglycerols (TAG)", "Diacylglycerols (DAG)",
      "Monoacylglycerols (MAG)", "Galactosyldiacylglycerols (MGDG)",
      "Digalactosyldiacylglycerols (DGDG)", "Sulfoquinovosyldiacylglycerols (SQDG)",

      # Sterol Lipids
      "Cholesterol (Chol)", "Cholesteryl esters (CE)",
      "Bile acids (BA)", "Steroid hormones",

      # Prenol Lipids
      "Ubiquinones (CoQ)", "Dolichols", "Prenyl phosphates"
    ),
    Examples = c(
      # Sphingolipids examples
      "Cer(d18:1/16:0), Cer(d18:1/24:0)", "SM(d18:1/16:0), SM(d18:1/24:1)",
      "beta-GlcCer(d18:1/16:0)", "beta-GalCer(d18:1/18:0)",
      "LacCer(d18:1/16:0)", "Sph(d18:1), Sph(d18:0)",
      "S1P(d18:1)", "C1P(d18:1/16:0)", "GM1(d18:1/18:0)",
      "ST(d18:1/24:0)", "Psy(d18:1)", "dhCer(d18:0/16:0)",

      # Glycerophospholipids examples
      "PC(16:0/18:1), PC(32:1)", "PE(16:0/18:1), PE(34:1)",
      "PS(18:0/18:1), PS(36:1)", "PG(16:0/18:1), PG(34:1)",
      "PI(16:0/18:1), PI(34:1)", "PA(16:0/18:1), PA(34:1)",
      "LPC(16:0), LPC(18:1)", "LPE(16:0), LPE(18:1)",
      "LPS(18:0), LPS(18:1)", "LPG(16:0), LPG(18:1)",
      "LPI(18:0), LPI(18:1)", "LPA(16:0), LPA(18:1)",
      "CL(18:1)4, CL(72:4)", "PIP(34:1), PIP2(34:1)",
      "PE(O-16:0/18:1), PC(P-16:0/18:1)",

      # Glycerolipids examples
      "TAG(16:0/18:1/18:2)", "DAG(16:0/18:1)",
      "MAG(18:1)", "MGDG(16:0/18:3)",
      "DGDG(16:0/18:3)", "SQDG(16:0/18:3)",

      # Sterol examples
      "Cholesterol", "CE(18:1), CE(20:4)",
      "Cholic acid, CDCA", "Cortisol, Testosterone",

      # Prenol examples
      "CoQ10, CoQ9", "Dolichol-P", "Geranyl-PP"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Create the datatable with fallback data
  dt <- DT::datatable(
    lipid_classes_data,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      scrollY = "400px",
      dom = 'ft',
      columnDefs = list(
        list(width = '150px', targets = 0),
        list(width = '250px', targets = 1),
        list(width = '300px', targets = 2)
      )
    ),
    rownames = FALSE,
    caption = "Comprehensive list of lipid classes supported by LipidCRED (fallback data)"
  ) %>%
    DT::formatStyle(
      'Category',
      target = 'row',
      backgroundColor = DT::styleEqual(
        c('Sphingolipids', 'Glycerophospholipids', 'Glycerolipids', 'Sterol Lipids', 'Prenol Lipids'),
        c('#f8f9fa', '#e3f2fd', '#f3e5f5', '#fff3e0', '#e8f5e8')
      )
    )

  return(dt)
}