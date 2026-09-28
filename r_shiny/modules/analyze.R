analyzePageUI <- function(id, organism_choices = NULL) {
  ns <- NS(id)
  # Use fallback if organism_choices is NULL
  if (is.null(organism_choices)) {
    organism_choices <- c(
      " -- " = " -- ",
      "Homo sapiens" = 9606,
      "Mus musculus" = 10090,
      "Rattus norvegicus" = 10116
    )
  }
  tabItem(
    tabName = "analyze",
    fluidRow(
      tabBox(
        title = "",
        id = ns("thisBox"),
        width = 12,
        tabPanel("Batch Processing",
                 h2("Comprehensive Dataset Analysis"),
                 p("Process entire lipid datasets to discover metabolic networks and enzymatic relationships."),

                 h3("Step 1: Data Upload"),
                 p("Upload your lipid dataset in Excel format (.xlsx). See the Getting Started page for formatting requirements."),
                 fileInput(inputId = "input_lipid_data",
                           multiple = TRUE,
                           label = NULL,
                           width = "50%",
                           buttonLabel = "Browse...",
                           placeholder = "No file selected",
                           accept = ".xlsx"
                 ),

                 h3("Step 2: Configure Analysis"),
                 selectizeInput(
                   inputId = "select_organism",
                   label = "Select organism (type to search)",
                   choices = organism_choices,
                   selected = " -- ",
                   multiple = FALSE,
                   options = list(
                     placeholder = "Type organism name to search...",
                     create = FALSE
                   )
                 ),
                 p("Otherwise please enter the taxon ID of the desired organism below:"),

                 textInput(inputId = "text_select_organism",
                           label = NULL,
                           value = "",
                           width = '100px',
                           placeholder = "Taxon ID"
                 ),

                 # p(
                 #   "This will take between 1 and 5 minutes depending on input size and if translation is necessary."
                 # ),
                 # checkboxInput("translate",
                 #               "Translate lipids?",
                 #               value = FALSE),
                 #checkboxInput("processNetwork", "Create 3D network visualization?", value = FALSE), TODO: uncomment when ready
                 tags$div(style = "display: none;",
                          checkboxInput("processNetwork", "Create 3D network visualization?",
                                        value = FALSE)),
                 h3("Step 3: Run Analysis"),
                 p("Processing time varies from 1-5 minutes depending on dataset size and complexity."),
                 actionButton("createAdjMatrix",
                              label = "Run",
                              icon = icon("play"),
                              class = "btn-warning",
                              style = "color: #fff;"),
                 textOutput("downloadMessage"),

                 h3("Step 4: Download Results"),

                 downloadButton("downloadMatrix",
                                label = "Download matrices + reaction list + enzyme list",
                                icon = icon("download"),
                                class = "btn-warning",
                                style = "color: #fff;"
                 ),
        ),
        tabPanel("Enzyme list",
                 fluidRow(
                   column(12,
                          h2("Enzyme Table"),
                          uiOutput("enzymeTableMessage"),
                          DT::dataTableOutput("enzymeTable"),
                          br(),
                          downloadButton("downloadEnzymeTable", "Download Enzyme Table",
                                         class = "btn-warning", style = "color: #fff;" )
                   )
                 )

        ),
        # --- NEW TAB PANEL FOR SINGLE SEARCH ---
        tabPanel("Single Entity Search",
                 # Use fluidRow and column structure instead of nested tabBox
                 fluidRow(
                   column(12,
                          h3("Single Entity Search"),
                          p("Use the sections below to search for a single lipid or a single enzyme."),
                          br()
                   )
                 ),

                 # --- Single Lipid Search Section ---
                 fluidRow(
                   column(12,
                          box(
                            title = "Single Lipid Search",
                            status = "warning",
                            solidHeader = TRUE,
                            width = 12,
                            collapsible = TRUE,
                            collapsed = FALSE,

                            # --- Inputs Section ---
                            h4("Search Inputs"),
                            p("Enter one or more lipid names (one per line) to find all reactions they participate in (both as a reactant and as a product)."),

                            fluidRow(
                              column(6,
                                     textAreaInput(
                                       inputId = "lipid_input_single",
                                       label = "Lipid Name(s):",
                                       placeholder = "e.g., PC(16:0/18:1)\nCer(d18:1/16:0)",
                                       rows = 4,
                                       width = "100%"
                                     )
                              ),
                              column(6,
                                     selectizeInput(
                                       inputId = "select_organism_lipid",
                                       label = "Select organism",
                                       choices = organism_choices,
                                       selected = " -- ",
                                       multiple = FALSE,
                                       options = list(
                                         placeholder = "Type organism name...",
                                         create = FALSE
                                       )
                                     ),
                                     textInput(
                                       inputId = "text_select_organism_lipid",
                                       label = "Or enter Taxon ID",
                                       placeholder = "Taxon ID",
                                       width = '150px'
                                     )
                              )
                            ),

                            actionButton(
                              inputId = "runSingleLipidSearch",
                              label = "Search Lipid(s)",
                              icon = icon("search"),
                              class = "btn-warning",
                              style = "color: #fff;"
                            ),

                            # --- Results Section ---
                            br(),
                            hr(),
                            h4("Lipid Search Results"),
                            DT::dataTableOutput("singleLipidResultsTable"),
                            br(),
                            uiOutput("downloadSingleLipidUI")
                          )
                   )
                 ),

                 # --- Single Enzyme Search Section ---
                 fluidRow(
                   column(12,
                          box(
                            title = "Single Enzyme Search",
                            status = "warning",
                            solidHeader = TRUE,
                            width = 12,
                            collapsible = TRUE,
                            collapsed = TRUE,

                            # --- Inputs Section ---
                            h4("Search Inputs"),
                            p("Enter an enzyme identifier (UniProt ID or gene name) to find all associated reactions."),

                            fluidRow(
                              column(6,
                                     textInput(
                                       inputId = "enzyme_identifier_single",
                                       label = "Enzyme Identifier:",
                                       placeholder = "e.g., P04075 or SMS1",
                                       width = "100%"
                                     )
                              ),
                              column(6,
                                     selectizeInput(
                                       inputId = "select_organism_enzyme",
                                       label = "Select organism",
                                       choices = organism_choices,
                                       selected = " -- ",
                                       multiple = FALSE,
                                       options = list(
                                         placeholder = "Type organism name...",
                                         create = FALSE
                                       )
                                     ),
                                     textInput(
                                       inputId = "text_select_organism_enzyme",
                                       label = "Or enter Taxon ID",
                                       placeholder = "Taxon ID",
                                       width = '150px'
                                     )
                              )
                            ),

                            actionButton(
                              inputId = "runSingleEnzymeSearch",
                              label = "Search Enzyme",
                              icon = icon("search"),
                              class = "btn-warning",
                              style = "color: #fff;"
                            ),

                            # --- Results Section ---
                            br(),
                            hr(),
                            h4("Enzyme Search Results"),
                            DT::dataTableOutput("singleEnzymeResultsTable"),
                            br(),
                            uiOutput("downloadSingleEnzymeUI")
                          )
                   )
                 )
        )
        # ----------------------------------------
        # tabPanel("2. Arena 3D Visualization",
        #          h2("3D Network Visualization"),
        #          uiOutput("arena3dLink")

        # tabPanel("2. Arena 3D Visualization", TODO: uncomment when ready
        #          h2("3D Network Visualization"),
        #          uiOutput("arena3dLink")
                 # Everything below is commented out and unnecessary
                 # fileInput(inputId = "reactionList",
                 #           multiple = FALSE,
                 #           label = NULL,
                 #           width = "50%",
                 #           buttonLabel = "Browse...",
                 #           placeholder = "No file selected",
                 #           accept = ".xlsx"
                 # ),
                 # actionButton("sendRequest",
                 #              label = "Send JSON to Arena3D",
                 #              icon = icon("person-running"),
                 #              class = "btn-warning",
                 #              style = "color: #fff;"
                 #
                 # ),
                 # textOutput("response_text")
        # ),

      )
    )
  )
}
