sendNetworkToAPI <- function(network_file_path) {
  req(network_file_path)

  network_json <- readLines(network_file_path, warn = FALSE)
  network_json <- paste(network_json, collapse = "")

  response <- httr::POST(
    url = "https://bib.fleming.gr/bib/api/arena3dweb",
    body = network_json,
    httr::add_headers("Content-Type" = "application/json")
  )

  if (httr::status_code(response) == 200) {
    response_content <- httr::content(response, "parsed")
    if (!is.null(response_content$url)) {
      return(response_content$url)
    } else {
      warning("API response does not contain a URL")
      return(NULL)
    }
  } else {
    warning(paste("API request failed with status code:", httr::status_code(response)))
    return(NULL)
  }
}


renderArena3dLink <- function(arena3d_link) {
  if (is.null(arena3d_link)) {
    return(HTML("<p>Arena3D link will appear here after processing if you checked the option.</p>"))
  } else {
    return(tags$a(href = arena3d_link, "Open in Arena3D", target = "_blank", class = "btn btn-primary"))
  }
}
