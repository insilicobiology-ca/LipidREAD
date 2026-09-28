citePageUI <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "cite",
    box(width = 12,
        column(
          width = 12,
          title = "",
          id = ns("citeBox"),
          div(
            h3("Contact and Support"),
            p("For questions, bug reports, or collaboration inquiries:"),
            tags$ul(
              tags$li("Email: ", a("ldomic@uottawa.ca", href = "mailto:ldomic@uottawa.ca")),
              tags$li("Include your dataset and detailed problem description for technical support")
            ),
            br(),
            h3("Cite the use of LipidREAD in a publication"),
            box(width = 12, status = "warning",
                p(em("Alkassir, Q., [Author List], Cuperlovic-Culf, M. & Bennett, S.A.L. (2025).
                   LipidREAD: A Comprehensive Lipid Reaction Enzyme Annotation Database
                   for Lipidomics Analysis. [Journal Name], [Volume], [Pages].")),

                p("Manuscript in preparation - please check our website for the most current citation information.")
            ),
            br(),
            h3("Public Server"),
            p("LipidREAD: ", a("https://insilicobiology.ca/LipidREAD/", href = "https://insilicobiology.ca/shiny/LipidREAD/")),
            br(),
            h3("Content of LipidREAD"),
            p("LipidREAD is provided 'AS IS', without warranty of any kind, express or implied")),
            br(),
            h3("Software License"),
            p("All LipidREAD® databases are licensed under a ",a("Creative Commons Attribution 4.0 International License.",
            href="http://creativecommons.org/licenses/by/4.0/"),

          )
        )
    )
  )
}
