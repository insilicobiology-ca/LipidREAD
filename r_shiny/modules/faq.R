faqPageUI <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "faq",
    box(width = 12,
        column(
          width = 12,
          title = "",
          id = ns("faqBox"),
          div(
            h1("Frequently Asked Questions"),

            h3("Supported Lipid Classes"),
            p("LipidREAD supports the following major lipid categories:"),

            fluidRow(
              box(width = 12, title = "Supported Lipid Classes", status = "warning", solidHeader = TRUE,
                  DT::dataTableOutput("supportedClassesTable")
              )
            ),



            h3("Nomenclature Guidelines"),

            h4("Q: What lipid name formats are accepted?"),
            p("LipidREAD accepts multiple nomenclature standards:"),
            tags$ul(
              tags$li(strong("Shorthand notation:"), " PC(16:0/18:1), Cer(d18:1/16:0)"),
              tags$li(strong("Sum composition:"), " PC(34:1), Cer(34:1)"),
              tags$li(strong("Systematic names:"), " Various IUPAC and common names"),
              tags$li(strong("Hex/Hex2 conversions:"), " Automatically converts Hex to beta-Glc/beta-Gal and Hex2 to Lac/Gala")
            ),

            h4("Q: My lipid wasn't found. What should I do?"),
            tags$ul(
              tags$li("Check spelling and formatting against examples above"),
              tags$li("Ensure the lipid class is supported (see list above)"),
              tags$li("Use the Single Lipid Search to test individual lipids"),
              tags$li("Contact us if you believe a supported lipid should be recognized")
            ),

            h3("Performance and Limitations"),

            h4("Q: Why is my analysis running slowly?"),
            p("Several factors affect processing time:"),
            tags$ul(
              tags$li(strong("Dataset size:"), " Larger datasets (>500 lipids) require more time"),
              tags$li(strong("Lipid complexity:"), " Species-level lipids require more database queries"),
              tags$li(strong("Translation needs:"), " Non-standard names require additional processing"),
              tags$li(strong("Server load:"), " Shared resources may cause delays")
            ),

            h4("Q: Why are my adjacency matrices empty?"),
            p(HTML("Adjacency matrix generation is limited to datasets with &le;2000 unique lipids.
              For larger datasets, matrices will be empty, but the complete reaction list will still be generated with all identified reactions.")),

            h3("Technical Questions"),

            h4("Q: What organisms are supported?"),
            p("LipidREAD includes enzymatic data for multiple organisms. Use the dropdown menu or enter any NCBI Taxon ID.
              Common organisms include:"),
            tags$ul(
              tags$li("Homo sapiens (9606)"),
              tags$li("Mus musculus (10090)"),
              tags$li("Rattus norvegicus (10116)")
            ),

            h4("Q: How current is the database?"),
            p("LipidREAD database is updated monthly from SwissLipids and Rhea databases to ensure current biochemical knowledge."),
          )
        )
    )
  )
}