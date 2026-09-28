source("r_shiny/modules/home.R")
source("r_shiny/modules/sample.R")
source("r_shiny/modules/analyze.R")
source("r_shiny/modules/troubleshoot.R")
source("r_shiny/modules/faq.R")
source("r_shiny/modules/cite.R")


sidebar <- function() {
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      div(
        align = "center",
        br(),
        tags$figure(
          tags$img(src = "LipidREAD_just_logo.svg", width = 100, alt = "CompLi logo")
        ),
        h3("LipidREAD v1.0"),
        br()
      ),
      menuItem("Getting started", tabName = "start", icon = icon("gauge")),
      menuItem("Download sample data", tabName = "sample", icon = icon("download")),
      menuItem("Analyze", tabName = "analyze", icon = icon("magnifying-glass")),
      menuItem("FAQ", tabName = "faq", icon = icon("question")),
      menuItem("Troubleshoot", tabName = "troubleshoot", icon = icon("screwdriver-wrench")),
      menuItem("Authors and citing", tabName = "cite", icon = icon("pencil")),
      menuItem("Return to In Silico Biology", icon = icon("house"), href = "https://insilicobiology.ca")
    )
  )
}


body <- function(organism_choices = NULL) {
  # Set default if not provided
  if (is.null(organism_choices)) {
    organism_choices <- c(
      " -- " = " -- ",
      "Homo sapiens" = 9606,
      "Mus musculus" = 10090,
      "Rattus norvegicus" = 10116
    )
  }
  dashboardBody(
    useShinyjs(),

    use_busy_spinner(
      spin = "fading-circle",
      position = "bottom-right",
      spin_id = "spinner" # Use the same ID
    ),

    tabItems(
      tabItem(tabName = "start", startPageUI("start_page")),
      tabItem(tabName = "sample", samplePageUI("sample_page")),
      tabItem(tabName = "analyze", analyzePageUI("analyze_page", organism_choices)),
      tabItem(tabName = "faq", faqPageUI("faq_page")),
      tabItem(tabName = "troubleshoot", troubleshootPageUI("troubleshoot_page")),
      tabItem(tabName = "cite", citePageUI("cite_page")),
      tabItem(tabName = "house")
    ),
    tags$head(tags$script(HTML(
      '$(document).ready(function() {
        $("header").find("nav").append(\'<span class="myClass">Lipids Reactions and Enzymes Annotation Database</span>\');
      })'
    )))
  )
}
