startPageUI <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "start",
    box(width = 12,
        column(
          12,
          h1("Overview of LipidREAD"),
          p(
            "LipidREAD (Lipids Reactions and Enzymes Annotation Database) is a lipidomics analysis software and
            relational database that empowers researchers
              to uncover potential metabolic pathways within their datasets. It integrates and generalizes
              lipid reaction data from public sources, includes orthologous enzymatic drivers, and handles
              various lipid nomenclatures. Critically, LipidREAD employs intelligent carbon-chain matching algorithms
              to predict specific molecular products, enabling the connection of general biochemical knowledge to
              experimentally identified lipids."
          ),

          fluidRow(
            box(width = 6, title = "Database Statistics", status = "warning", solidHeader = TRUE,
                tags$ul(
                  tags$li("2,796 lipid reactions from SwissLipids"),
                  tags$li("4,281 enzymes across multiple organisms"),
                  tags$li("700 supported organisms (including Homo sapiens, Mus musculus, Rattus norvegicus)"),
                  tags$li("XX lipid classes and subclasses supported")
                )
            ),
            box(width = 6, title = "Analysis Modes", status = "warning", solidHeader = TRUE,
                tags$ul(
                  tags$li(strong("Batch Processing:"), " Upload datasets for comprehensive network analysis"),
                  tags$li(strong("Single Lipid Search:"), " Find all reactions involving specific lipids"),
                  tags$li(strong("Single Enzyme Search:"), " Discover all reactions catalyzed by specific enzymes")
                )
            )
          ),
          br(),

          h2("Analysis Workflow"),

          # Processing mode tabs
          tabBox(
            title = "",
            width = 12,

            tabPanel("Batch Processing Mode",
                     h3("For comprehensive dataset analysis:"),

                     h4("Data Preparation:"),
                     tags$ol(
                       tags$li("Prepare your lipid data in Excel format (.xlsx)"),
                       tags$li("List lipids in the first row, starting from column A"),
                       tags$li("Ensure lipid names follow standard nomenclature (see FAQ for supported formats)")
                     ),

                     # Example table
                     h4("Example Input Format:"),
                     fluidRow(
                       box(width = 8,
                           tags$table(class = "table table-striped",
                                      tags$thead(
                                        tags$tr(
                                          tags$th(""),
                                          tags$th("A"),
                                          tags$th("B"),
                                          tags$th("C"),
                                          tags$th("D"),
                                          tags$th("E"),
                                          tags$th("...")
                                        )
                                      ),
                                      tags$tbody(
                                        tags$tr(
                                          tags$td("1"),
                                          tags$td("Cer(d18:1/16:0)"),
                                          tags$td("Cer(d18:1/20:0)"),
                                          tags$td("PC(16:0/18:0)"),
                                          tags$td("LPC(16:0/0:0)"),
                                          tags$td("PS(18:0/20:0)"),
                                          tags$td("...")
                                        )
                                      )
                           )
                       )
                     ),

                     h4("Analysis Steps:"),
                     tags$ol(
                       tags$li("Upload your formatted dataset"),
                       tags$li("Select organism or enter taxon ID"),
                       tags$li("Click 'Run Analysis'"),
                       tags$li("Download comprehensive results package")
                     ),

                     h4("Output Files:"),
                     tags$ul(
                       tags$li(strong("Adjacency Matrices:"), " Network connectivity (binary and enzyme-annotated)"),
                       tags$li(strong("Reaction Lists:"), " Complete reactions and simplified pairs"),
                       tags$li(strong("Enzyme Table:"), " All enzymes involved in your dataset"),
                       tags$li(strong("Translation Tracker:"), " Maps original to standardized lipid names")
                     )
            ),

            tabPanel("Single Entity Search",
                     h3("For targeted analysis:"),

                     h4("Single Lipid Search:"),
                     tags$ul(
                       tags$li("Enter one or more lipid names (one per line)"),
                       tags$li("Finds all reactions where lipids participate as reactants or products"),
                       tags$li("No constraints on product availability"),
                       tags$li("Results include reaction details and enzymatic information")
                     ),

                     h4("Single Enzyme Search:"),
                     tags$ul(
                       tags$li("Enter enzyme identifier (UniProt ID or gene name)"),
                       tags$li("Discovers all reactions catalyzed by the enzyme"),
                       tags$li("Results aggregated across multiple reaction databases"),
                       tags$li("Reactions shown at lipid class level for clarity")
                     ),

                     h4("Output Format:"),
                     p("Results displayed in interactive tables with download options as TSV files.")
            ),
          )
    ))
  )
}
