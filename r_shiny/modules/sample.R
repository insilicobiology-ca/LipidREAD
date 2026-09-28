samplePageUI <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "sample",
    box(width = 12,
        column(
          width = 12,
          title = "",
          id = ns("sampleBox"),
          div(
            h1("Sample Datasets"),
            p("Download example datasets to test LipidREAD functionality and learn proper formatting."),

            fluidRow(
              box(width = 6, title = "Batch Processing Sample", status = "warning",
                  h4("Comprehensive Dataset Example"),
                  p("Sample Excel file with multiple lipid classes for batch analysis demonstration."),
                  tags$ul(
                    tags$li("50 diverse lipids across major classes"),
                    tags$li("Proper formatting for batch upload"),
                    tags$li("Expected processing time: ~2 minutes"),
                    tags$li("Generates all output file types")
                  ),
                  downloadButton("downloadBatchSample", "Download Batch Sample (.xlsx)",
                                class = "btn-orange")
              ),

              box(width = 6, title = "Single Search Examples", status = "warning",
                  h4("Individual Entity Examples"),
                  p("Example inputs for single lipid and enzyme searches."),

                  h5("Single Lipid Examples:"),
                  tags$ul(
                    tags$li("PC(16:0/18:1)"),
                    tags$li("Cer(d18:1/16:0)"),
                    tags$li("beta-GlcCer(d18:1/16:0)")
                  ),

                  h5("Single Enzyme Examples:"),
                  tags$ul(
                    tags$li("Q86VZ5 (SGMS1)"),
                    tags$li("SPHK1"),
                    tags$li("CERK")
                  ),

                  p("Copy these examples directly into the search fields.")
              )
            ),

            h3("Dataset Descriptions"),

            h4("Batch Processing Sample Details:"),
            p("The batch processing sample contains representative lipids from:"),
            tags$ul(
              tags$li("Sphingolipids: Ceramides, sphingomyelins, glycosphingolipids"),
              tags$li("Glycerophospholipids: PC, PE, PS, PG, PI and their lyso forms"),
              tags$li("Various chain lengths and saturation levels"),
              tags$li("Both systematic and common nomenclature examples")
            ),

            h4("Expected Results:"),
            p("Processing the sample dataset will generate:"),
            tags$ul(
              tags$li("~200-300 reactions in the full reaction list"),
              tags$li("~50-100 unique enzymes"),
              tags$li("Network connections between related lipids"),
              tags$li("Translation mappings for standardized names")
            ),

            h3("Using Sample Data"),
            tags$ol(
              tags$li("Download the appropriate sample file"),
              tags$li("Navigate to the Analyze tab"),
              tags$li("Upload the sample file using the file input"),
              tags$li("Select 'Homo sapiens' as the organism"),
              tags$li("Click 'Run Analysis' and wait for completion"),
              tags$li("Download and examine the results")
            )
          )
        )
    )
  )
}