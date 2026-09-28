troubleshootPageUI <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "troubleshoot",
    box(width = 12,
        column(
          width = 12,
          title = "",
          id = ns("troubleshootBox"),
          div(
            h1("Troubleshooting Guide"),
            p("Common issues and solutions for LipidREAD users. If your problem persists, contact ",
              a("ldomic@uottawa.ca", href = "mailto:ldomic@uottawa.ca"),
              " with your dataset and problem description."),

            h3("Data Input Issues"),

            h4("Problem: Lipids not recognized or 'not found' errors"),
            p(strong("Possible causes and solutions:")),
            tags$ul(
              tags$li(strong("Nomenclature mismatch:"), " Check the FAQ page for supported lipid classes and naming conventions"),
              tags$li(strong("Formatting errors:"), " Ensure proper spacing, parentheses, and separators (e.g., '/' or '_')"),
              tags$li(strong("Unsupported class:"), " Verify your lipids belong to supported categories (see FAQ)"),
              tags$li(strong("Spelling mistakes:"), " Double-check lipid names against standard nomenclature"),
              tags$li(strong("Special characters:"), " Avoid unusual characters or encoding issues")
            ),

            box(width = 12, status = "warning",
                h5("Quick Test:"),
                p("Use the Single Lipid Search to test problematic lipids individually before batch processing.")
            ),

            h4("Problem: File upload fails"),
            tags$ul(
              tags$li("Ensure file is in .xlsx format (not .xls or .csv)"),
              tags$li("Check that lipids are in the first row starting from column A"),
              tags$li("Verify file size is reasonable (<10MB recommended)"),
              tags$li("Try removing empty columns/rows from your spreadsheet")
            ),

            h3("Performance Issues"),

            h4("Problem: Analysis running very slowly (>10 minutes)"),
            p(strong("Common causes:")),
            tags$ul(
              tags$li(strong("Large dataset:"), " >500 lipids significantly increase processing time"),
              tags$li(strong("Complex lipids:"), " Species-level lipids (e.g., PC(16:0/18:1)) require less database queries than sum compositions (e.g., PC(34:1))"),
              tags$li(strong("Translation overhead:"), " Non-standard names require additional processing"),
              tags$li(strong("Server load:"), " Shared server resources during peak usage")
            ),

            p(strong("Solutions:")),
            tags$ul(
              tags$li("Try analysis during off-peak hours"),
              tags$li("Consider using Single Lipid Search for exploratory analysis")
            ),

            h3("Output Issues"),

            h4("Problem: Empty adjacency matrices"),
            p(HTML("Adjacency matrix generation is limited to datasets with &le;2000.
              The complete reaction list will still contain all identified reactions.")),
            p(strong("Solutions:")),
            tags$ul(
              tags$li("Use the reaction lists for network analysis"),
              tags$li("Filter your dataset to <2000 lipids for matrix generation"),
              tags$li("Focus on specific lipid classes of interest")
            ),

            h4("Problem: No reactions found"),
            tags$ul(
              tags$li("Check organism selection - ensure appropriate species is selected"),
              tags$li("Verify lipid nomenclature (see FAQ)"),
              tags$li("Consider that some lipids may not have known enzymatic reactions"),
            ),

            h3("Browser and Technical Issues"),

            h4("Problem: Interface not loading or freezing"),
            tags$ul(
              tags$li("Refresh the browser page"),
              tags$li("Clear browser cache and cookies"),
              tags$li("Try a different browser (Chrome, Firefox, Safari)"),
              tags$li("Disable browser extensions that might interfere"),
              tags$li("Check internet connection stability")
            ),

            h4("Problem: Download fails"),
            tags$ul(
              tags$li("Ensure popup blockers are disabled"),
              tags$li("Try right-clicking download button and 'Save link as...'"),
              tags$li("Check available disk space"),
              tags$li("Verify file permissions in download directory")
            ),

            h3("Getting Help"),

            box(width = 12, status = "warning",
                h4("When contacting support, please include:"),
                tags$ul(
                  tags$li("Your input dataset (or representative sample)"),
                  tags$li("Detailed description of the problem"),
                  tags$li("Organism selected"),
                  tags$li("Browser and operating system information"),
                  tags$li("Any error messages received"),
                  tags$li("Screenshots if applicable")
                )
            ),
          )
        )
    )
  )
}