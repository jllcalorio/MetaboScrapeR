# ================================================================
# Load all required libraries
# ================================================================
library(shiny)
library(shinyjs)
library(DT)
library(writexl)

# Core processing libraries
library(rvest)
library(stringr)
library(dplyr)
library(purrr)
library(httr)
library(jsonlite)

# Parallel processing libraries for performance
library(future)
library(furrr)

# ================================================================
# Helper Function: Fetch Metadata via API
# ================================================================

#' Fetches HMDB compound metadata via the public API in parallel
#'
#' @param compound_ids A character vector of HMDB compound IDs (e.g., "HMDB0000001").
#' @return A data frame of compound metadata.
#' @noRd
fetch_metadata_via_api <- function(compound_ids) {

  if (is.null(compound_ids) || length(compound_ids) == 0) {
    return(data.frame())
  }

  # Define the API endpoint
  api_url <- "https://hmdb.ca/metabolites/"

  # Safely extract a value from a list or return NA
  safe_extract <- function(x, name) {
    purrr::pluck(x, name, .default = NA_character_)
  }

  # Use furrr to map over IDs in parallel
  metadata_df <- furrr::future_map_dfr(compound_ids, ~ {
    hmdb_id <- .x
    full_url <- paste0(api_url, hmdb_id, ".json")

    # Add timestamped log message
    message(sprintf("[%s] Fetching API metadata for %s...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), hmdb_id))

    tryCatch({
      # Make the API request
      resp <- httr::GET(full_url, httr::accept_json())

      # Check for successful request
      if (httr::status_code(resp) != 200) {
        warning(sprintf("[%s] API request failed for HMDB ID %s with status %d", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), hmdb_id, httr::status_code(resp)), call. = FALSE)
        return(NULL) # Return NULL, map_dfr will skip it
      }

      # Parse the JSON content
      content_json <- httr::content(resp, as = "parsed", type = "application/json")

      # Extract required fields using the safe helper
      tibble::tibble(
        Compound = hmdb_id,
        URL = full_url,
        Description = safe_extract(content_json, "description"),
        Kingdom = safe_extract(content_json, "taxonomy", "kingdom"),
        SuperClass = safe_extract(content_json, "taxonomy", "super_class"),
        Class = safe_extract(content_json, "taxonomy", "class"),
        Subclass = safe_extract(content_json, "taxonomy", "sub_class"),
        Direct_Parent = safe_extract(content_json, "taxonomy", "direct_parent"),
        Alternative_Parent = paste(safe_extract(content_json, "taxonomy", "alternative_parents"), collapse = "; "),
        Substituents = paste(safe_extract(content_json, "taxonomy", "substituents"), collapse = "; "),
        Molecular_Framework = safe_extract(content_json, "taxonomy", "molecular_framework"),
        External_Descriptors = safe_extract(content_json, "taxonomy", "external_descriptors"),
        Pathway_Name = paste(purrr::map_chr(content_json$pathways, "name"), collapse = "; "),
        Disposition = safe_extract(content_json, "ontology", "disposition"),
        Process = safe_extract(content_json, "ontology", "physiological_relevance")
      )

    }, error = function(e) {
      warning(sprintf("[%s] Could not fetch API data for HMDB ID %s: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), hmdb_id, e$message), call. = FALSE)
      return(NULL) # Return NULL on error
    })
  }, .options = furrr::furrr_options(seed = NULL))

  return(metadata_df)
}

# ================================================================
# Core Function: Query HMDB
# ================================================================

#' Query HMDB MS Spectra and Retrieve Compound Metadata (Optimized)
#'
#' @param mz_rt Character string or vector. Feature m/z and RT in the format
#'   `"m/z@RT"`. Defaults to `NULL`. Can accept multiple values.
#' @param ion_mode Character string. MS ionization mode. One of
#'   `"Positive"`, `"Negative"`, or `"Neutral"`. Default = `"Positive"`.
#' @param filter_adduct_type Character vector. Specific adducts to retain.
#'   Defaults to `NULL`. If provided, results are filtered accordingly.
#' @param tolerance Numeric. Mass tolerance value. Default = `5`.
#' @param tolerance_units Character string. Mass tolerance unit.
#'   One of `"Da"` or `"ppm"`. Default = `"ppm"`.
#' @param ccs_predictors Character string or vector. CCS prediction model(s).
#'   Examples: `"AllCCS"`, `"DarkChem"`, `"DeepCCS"`. Default = `NULL`.
#' @param ccs_tolerance Numeric. CCS tolerance value. Choose between c(1, 3, 5, 10). Default = `NULL`.
#'
#' @return A list with two data frames:
#'   \itemize{
#'     \item \code{compound_metadata}: Compound details retrieved via the HMDB API,
#'       merged with the initial search results.
#'     \item \code{results_table_excluded}: Table of results where adducts
#'       do not match the provided filter list.
#'   }
#'
query_hmdb_ms <- function(mz_rt = NULL,
                          ion_mode = "Positive",
                          filter_adduct_type = NULL,
                          tolerance = 5,
                          tolerance_units = "ppm",
                          ccs_predictors = NULL,
                          ccs_tolerance = NULL) {

  # --- Input Validation ---
  if (is.null(mz_rt) || length(mz_rt) == 0) {
    stop("Argument 'mz_rt' must be provided and non-empty.", call. = FALSE)
  }
  invalid_mz_rt <- !grepl("^[0-9.]+@[0-9.]+$", mz_rt)
  if (any(invalid_mz_rt)) {
    stop("Invalid mz_rt format. Expected 'mz@rt' format. Invalid entries: ", paste(mz_rt[invalid_mz_rt], collapse = ", "), call. = FALSE)
  }
  valid_ion_modes <- c("Positive", "Negative", "Neutral")
  if (!ion_mode %in% valid_ion_modes) {
    stop("ion_mode must be one of: ", paste(valid_ion_modes, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(tolerance) || tolerance <= 0) stop("tolerance must be a positive numeric value.", call. = FALSE)
  valid_tolerance_units <- c("Da", "ppm")
  if (!tolerance_units %in% valid_tolerance_units) {
    stop("tolerance_units must be one of: ", paste(valid_tolerance_units, collapse = ", "), call. = FALSE)
  }

  # Helper to process a single mz_rt entry
  process_mz_rt <- function(mz_entry) {

    # Add timestamped log message
    message(sprintf("[%s] === Processing feature: %s ===", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), mz_entry))

    tryCatch({
      # 1. Submit the search form
      form_values <- list(
        query_masses = as.numeric(stringr::str_extract(mz_entry, "^[^@]+")),
        ms_search_ion_mode = tolower(ion_mode),
        tolerance = tolerance,
        tolerance_units = tolerance_units,
        ccs_predictors = if (is.null(ccs_predictors)) "" else paste(ccs_predictors, collapse = ","),
        ccs_tolerance = if (is.null(ccs_tolerance)) "" else as.character(ccs_tolerance),
        "commit" = "Search"
      )

      message(sprintf("[%s] Submitting search to HMDB...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

      resp <- httr::POST("https://www.hmdb.ca/spectra/ms/search", body = form_values, encode = "form", httr::timeout(60))
      httr::stop_for_status(resp)

      results_page_raw <- httr::content(resp, as = "text")
      results_page_parsed <- rvest::read_html(results_page_raw)

      message(sprintf("[%s] Search submitted. Parsing results page...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

      # 2. Find the table and check if it exists
      results_table_node <- rvest::html_element(results_page_parsed, ".table")

      if (inherits(results_table_node, "xml_missing")) {
        page_text <- rvest::html_text(results_page_parsed, trim = TRUE)
        if (grepl("No results were found", page_text, ignore.case = TRUE)) {
          message(sprintf("[%s] No results found for feature '%s'.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), mz_entry))
        } else {
          warning(sprintf("[%s] Could not find results table for feature '%s'. The page layout may have changed.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), mz_entry), call. = FALSE)
        }
        return(list(compound_metadata = data.frame(), results_table_excluded = data.frame()))
      }

      # 3. Parse the table
      results_table <- rvest::html_table(results_table_node, fill = TRUE) %>%
        setNames(c("Compound", "Name", "Formula", "Monoisotopic Mass", "Adduct", "Adduct M/Z", "Delta (ppm)", "CCS"))

      message(sprintf("[%s] Found %d initial matches for feature '%s'.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), nrow(results_table), mz_entry))

      # 4. Clean and filter the table
      results_table$`Adduct M/Z` <- trimws(stringr::str_replace(results_table$`Adduct M/Z`, "m/z calculator", ""))

      # Filter for excluded results
      excluded <- if (!is.null(filter_adduct_type)) {
        dplyr::filter(results_table, !Adduct %in% filter_adduct_type)
      } else {
        data.frame() # If no filter, nothing is excluded by adduct
      }

      # Filter for included results
      included <- if (!is.null(filter_adduct_type)) {
        dplyr::filter(results_table, Adduct %in% filter_adduct_type)
      } else {
        results_table
      }

      if (!is.null(filter_adduct_type)) {
        message(sprintf("[%s] Filtered %d matches by adduct type. %d matches remaining.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), nrow(excluded), nrow(included)))
      }

      # 5. Fetch detailed metadata using the API
      compound_ids <- included$Compound

      if(length(compound_ids) > 0) {
        message(sprintf("[%s] Fetching detailed metadata for %d compound(s) via API...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(compound_ids)))
        compound_metadata <- fetch_metadata_via_api(compound_ids)
      } else {
        message(sprintf("[%s] No compounds remaining after filtering. Skipping API fetch.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
        compound_metadata <- data.frame()
      }

      # 6. Merge and add query metadata
      final_metadata <- if (nrow(compound_metadata) > 0) {
        dplyr::left_join(included, compound_metadata, by = "Compound")
      } else if (nrow(included) > 0) {
        # If API fetch failed but we have included results, return them
        included
      } else {
        data.frame()
      }

      query_meta <- data.frame(
        Ion = ion_mode,
        MZ = form_values$query_masses,
        RT = as.numeric(stringr::str_extract(mz_entry, "[^@]+$")),
        Tolerance = tolerance,
        Unit = tolerance_units,
        CCS_Predictors = paste(ccs_predictors, collapse = ","),
        CCS_Tolerance = paste(ccs_tolerance, collapse = ","),
        stringsAsFactors = FALSE
      )

      if (nrow(final_metadata) > 0) {
        final_metadata <- cbind(query_meta, final_metadata)
      }

      # Add MZ and RT to excluded results
      if (nrow(excluded) > 0) {
        excluded <- cbind(
          data.frame(
            MZ = form_values$query_masses,
            RT = as.numeric(stringr::str_extract(mz_entry, "[^@]+$")),
            stringsAsFactors = FALSE
          ),
          excluded
        )
      }

      message(sprintf("[%s] === Finished processing feature: %s. Found %d valid compounds. ===\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), mz_entry, nrow(final_metadata)))

      return(list(compound_metadata = final_metadata, results_table_excluded = excluded))

    }, error = function(e) {
      warning(sprintf("[%s] Failed to process feature '%s': %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), mz_entry, e$message), call. = FALSE)
      return(list(compound_metadata = data.frame(), results_table_excluded = data.frame()))
    })
  }

  # Process all mz_rt inputs and combine results
  results_list <- lapply(mz_rt, process_mz_rt)

  compound_metadata_all <- do.call(rbind, lapply(results_list, `[[`, "compound_metadata"))
  results_excluded_all <- do.call(rbind, lapply(results_list, `[[`, "results_table_excluded"))

  # Report
  total_found <- nrow(compound_metadata_all) + nrow(results_excluded_all)
  message(sprintf("[%s] === Total compounds found before adduct filtering: %d ===", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), total_found))

  n_excluded <- nrow(results_excluded_all)
  if (n_excluded > 0) {
    message(sprintf("[%s] Excluded %d compound(s) based on adduct filter. %d remaining.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n_excluded, nrow(compound_metadata_all)))
  } else {
    message(sprintf("[%s] No compounds were excluded by the adduct filter.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  }

  return(list(compound_metadata = compound_metadata_all,
              results_table_excluded = results_excluded_all))
}

# ================================================================
# UI
# ================================================================
ui <- fluidPage(
  useShinyjs(),
  titlePanel("HMDB MS Spectra Query App"),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      helpText("Enter m/z and RT values to query the HMDB."),

      textAreaInput("mz_rt_input",
                    "Features (m/z@RT)",
                    value = "",
                    rows = 4,
                    placeholder = "Example: 626.3547@95.2\n515.2597@10.7733"),
      helpText("Separate multiple features with newlines, commas, or both."),

      selectInput("ion_mode_input",
                  "Ion Mode",
                  choices = c("Positive", "Negative", "Neutral"),
                  selected = "Positive"),

      radioButtons(
        "adduct_input_mode",
        "Adducts to Keep:",
        choices = c("Select from list" = "select", "Paste comma-separated" = "paste"),
        selected = "select",
        inline = TRUE
      ),

      conditionalPanel(
        condition = "input.adduct_input_mode == 'select'",
        selectizeInput(
          "adduct_types_select",
          NULL, # No label, radio button acts as label
          # New comprehensive list of choices
          choices = c(
            "2M+2H+3H2O", "2M+ACN+H", "2M+ACN+Na", "2M+FA-H", "2M+H", "2M+H-H2O",
            "2M+K", "2M+Na", "2M-H", "3M-H", "M+2ACN+2H", "M+2ACN+H", "M+2H",
            "M+2H+Na", "M+2K-H", "M+2Na", "M+2Na-H", "M+3ACN+2H", "M+3H",
            "M+3Na", "M+ACN+2H", "M+ACN+H", "M+ACN+Na", "M+Br", "M+CH3OH+H",
            "M+Cl", "M+DMSO+H", "M+F", "M+FA-H", "M+H", "M+H+2K", "M+H+2Na",
            "M+H+HCOONa", "M+H+K", "M+H+Na", "M+H-2H2O", "M+H-H2O", "M+Hac-H",
            "M+IsoProp+H", "M+K", "M+K-2H", "M+Li", "M+NH4", "M+NH4-H2O", "M+Na",
            "M+Na-2H", "M+TFA-H", "M-2H", "M-3H", "M-H", "M-H+HCOONa", "M-H20-H",
            "Unknown"
          ),
          selected = NULL, # Server will set this dynamically
          multiple = TRUE,
          options = list(placeholder = 'Select adducts...')
        )
      ),

      conditionalPanel(
        condition = "input.adduct_input_mode == 'paste'",
        textAreaInput(
          "adduct_types_paste",
          NULL, # No label
          value = "", # Server will set this dynamically
          rows = 3,
          placeholder = "Example: M+H, M+Na, M+K"
        )
      ),

      helpText("See Help section for the full list of adducts."),

      numericInput("tolerance_input",
                   "Tolerance",
                   value = 5,
                   min = 1,
                   step = 1),

      selectInput("tolerance_units_input",
                  "Tolerance Units",
                  choices = c("ppm", "Da"),
                  selected = "ppm"),

      checkboxGroupInput("ccs_predictors_input",
                         "CCS Predictors (Optional)",
                         choices = c("AllCCS", "DarkChem", "DeepCCS")),

      checkboxGroupInput("ccs_tolerance_input",
                         "CCS Tolerance (Optional)",
                         choices = c("1", "3", "5", "10")),

      actionButton("run_query", "Run Query", class = "btn-primary", icon = icon("play")),

      hr(),

      downloadButton("download_excel", "Download Results (.xlsx)", class = "btn-success", style = "width: 100%;"),

      hr(),

      div(
        class = "well",
        HTML('<p>Source code on GitHub: <a href="https://github.com/jllcalorio/MetaboScrapeR" target="_blank">MetaboScrapeR</a></p>')
      )
    ),

    mainPanel(
      width = 8,
      shiny::fluidRow(
        shiny::column(width = 12,
                      align = "center",
                      div(id = "loading-message",
                          HTML("<div class='spinner'></div>"),
                          tags$p("Querying HMDB... This may take several minutes per feature."),
                          style = "display: none; color: #555; font-style: italic; margin-top: 15px;")
        )
      ),
      tags$head(
        tags$style(
          HTML("
            .spinner {
              border: 5px solid #f3f3f3;
              border-radius: 50%;
              border-top: 5px solid #3498db;
              width: 50px;
              height: 50px;
              -webkit-animation: spin 2s linear infinite;
              animation: spin 2s linear infinite;
              margin: 0 auto 10px auto; /* Centered */
            }
            @-webkit-keyframes spin {
              0% { -webkit-transform: rotate(0deg); }
              100% { -webkit-transform: rotate(360deg); }
            }
            @keyframes spin {
              0% { transform: rotate(0deg); }
              100% { transform: rotate(360deg); }
            }
            #compound_metadata_table table, #excluded_results_table table {
              white-space: nowrap !important;
            }
            /* Help tab styling */
            .help-content { padding: 10px; }
            .help-content h3 { border-bottom: 2px solid #eee; padding-bottom: 5px; margin-top: 20px;}
            .help-content h4 { margin-top: 15px; }
            .help-content code { background-color: #f5f5f5; border-radius: 3px; padding: 2px 4px; }
          ")
        )
      ),

      tabsetPanel(
        id = "main_tabs",

        tabPanel("Compound Metadata",
                 h4("Filtered Compound Metadata"),
                 DT::DTOutput("compound_metadata_table")),

        tabPanel("Excluded Results",
                 h4("Excluded Results Table (Adducts Not Filtered)"),
                 DT::DTOutput("excluded_results_table")),

        tabPanel("Help",
                 div(class = "help-content",
                     h3("Overview"),
                     p("This app queries the HMDB (Human Metabolome Database) MS search interface (",
                       tags$a(href = "https://www.hmdb.ca/spectra/ms/search", target = "_blank", "https://www.hmdb.ca/spectra/ms/search"),
                       ") to retrieve compound information based on mass spectrometry features. It uses the official HMDB API for fast and reliable data retrieval."),

                     h3("Input Format"),
                     h4("Features (m/z@RT)"),
                     p("Enter your features in the format: ", code("MZ@RT"), " where:"),
                     tags$ul(
                       tags$li(code("MZ"), " is the mass-to-charge ratio (e.g., 626.3547)"),
                       tags$li(code("RT"), " is the retention time (e.g., 95.2)"),
                       tags$li("Use ", code("@"), " (single at sign) as the separator"),
                       tags$li("Separate multiple features with newlines, commas, or both")
                     ),
                     p("Examples:"),
                     tags$pre("626.3547@95.2\n615.3753@120.5\n450.2315@80.0"),
                     p("OR:"),
                     tags$pre("626.3547@95.2, 615.3753@120.5, 450.2315@80.0"),
                     p("OR:"),
                     tags$pre("626.3547@95.2, 615.3753@120.5\n450.2315@80.0"),

                     h4("Ion Mode"),
                     p("Select the ionization mode used in your experiment. This will automatically update the default 'Adducts to Keep' list."),

                     h4("Adducts to Keep"),
                     p("You can either:"),
                     tags$ul(
                       tags$li(strong("Select from list:"), " Choose from predefined adduct types (recommended)."),
                       tags$li(strong("Paste comma-separated:"), " Paste your own list of adducts separated by commas.")
                     ),
                     p("If left blank, all results will be returned in the 'Compound Metadata' tab."),

                     p("Here are the adducts available in HMDB:"),
                     tags$ul(
                       tags$li(strong("Negative Mode:"),
                               tags$ul(
                                 tags$li("M-H, M-H20-H, M+Na-2H, M+Cl, M+K-2H, M+FA-H, M-H+HCOONa, 2M-H, 2M+FA-H, 3M-H, M-2H, M-3H")
                               )),
                       tags$li(strong("Positive Mode:"),
                               tags$ul(
                                 tags$li("M+H, M+H-2H2O, M+H-H2O, M+Na, M+CH3OH+H, M+K, M+ACN+H, M+2Na-H, M+ACN+Na, M+2K-H, M+2ACN+H, M+H+HCOONa, 2M+H, 2M+Na, 2M+2H+3H2O, 2M+K, 2M+ACN+H, 2M+ACN+Na, 2M+H-H2O, M+2H, M+H+Na, M+H+K, M+ACN+2H, M+2Na, M+2ACN+2H, M+3ACN+2H, M+3H, M+2H+Na, M+H+2Na, M+3Na, M+H+2K")
                               ))
                     ),
                     p("The default adducts selected by this app when you change Ion Mode are a curated list of the most common adducts from the lists above."),

                     h4("Tolerance"),
                     p("Set the mass tolerance for matching."),

                     h3("Packages Used"),
                     p("This app is built using the following R packages:"),
                     tags$ul(
                       tags$li(code("shiny"), " - Web application framework"),
                       tags$li(code("shinyjs"), " - JavaScript operations"),
                       tags$li(code("DT"), " - Interactive data tables"),
                       tags$li(code("writexl"), " - Excel file export"),
                       tags$li(code("rvest"), " - Web scraping"),
                       tags$li(code("stringr"), " - String manipulation"),
                       tags$li(code("dplyr"), " - Data manipulation"),
                       tags$li(code("purrr"), " - Functional programming"),
                       tags$li(code("httr"), " - HTTP requests"),
                       tags$li(code("jsonlite"), " - JSON parsing"),
                       tags$li(code("future"), ", ", code("furrr"), " - Parallel processing")
                     ),
                     p("Note: When accessing this app via shinyapps.io, all packages are pre-installed on the server. No installation is required on your end."),

                     h3("Output Tabs"),
                     tags$ul(
                       tags$li(strong("Compound Metadata:"), " Shows compounds that matched your search criteria AND your selected adducts. Includes detailed metadata from the HMDB API."),
                       tags$li(strong("Excluded Results:"), " Shows compounds that matched your search criteria but were filtered out because their adducts did not match your 'Adducts to Keep' list."),
                       tags$li(strong("Logs:"), " Shows a detailed, timestamped log of the entire query process.")
                     )
                 )
        ),

        tabPanel("Logs",
                 h4("Query Log"),
                 verbatimTextOutput("query_log"))
      )
    )
  )
)

# ================================================================
# Server
# ================================================================
server <- function(input, output, session) {

  # --- Set up parallel processing plan ---
  workers_available <- future::availableCores() - 1
  if (workers_available < 1) workers_available <- 1
  log_init_msg <- sprintf("[%s] Parallel backend initializing with %d workers...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), workers_available)
  future::plan(future::multisession, workers = workers_available)

  # --- Define Default Adduct Lists ---
  positive_adducts <- c("M+H", "M+H-2H2O", "M+H-H2O", "M+Na", "M+CH3OH+H", "M+K", "M+ACN+H", "M+2Na-H", "M+ACN+Na", "M+2K-H", "M+2ACN+H", "M+H+HCOONa", "2M+H", "2M+Na", "2M+2H+3H2O", "2M+K", "2M+ACN+H", "2M+ACN+Na", "2M+H-H2O", "M+2H", "M+H+Na", "M+H+K", "M+ACN+2H", "M+2Na", "M+2ACN+2H", "M+3ACN+2H", "M+3H", "M+2H+Na", "M+H+2Na", "M+3Na", "M+H+2K")
  negative_adducts <- c("M-H", "M-H20-H", "M+Na-2H", "M+Cl", "M+K-2H", "M+FA-H", "M-H+HCOONa", "2M-H", "2M+FA-H", "3M-H", "M-2H", "M-3H")

  # --- Reactive value to store log messages ---
  log_data <- reactiveVal(character(0))

  # --- Helper function to append log messages ---
  append_log <- function(msg) {
    current_log <- log_data()
    log_data(c(current_log, msg))
  }

  # --- Dynamic Adduct Updater ---
  observeEvent(input$ion_mode_input, {
    ion_mode <- input$ion_mode_input

    if (ion_mode == "Positive") {
      defaults <- positive_adducts
    } else if (ion_mode == "Negative") {
      defaults <- negative_adducts
    } else {
      defaults <- "" # c() # Neutral - clear selection
    }

    # Update the selectize input
    updateSelectizeInput(session, "adduct_types_select",
                         selected = defaults)

    # Update the paste input
    updateTextAreaInput(session, "adduct_types_paste",
                        value = paste(defaults, collapse = ", "))

    # Reset to "select" mode
    updateRadioButtons(session, "adduct_input_mode",
                       selected = "select")

  }, ignoreInit = FALSE) # Run on app start to set positive defaults

  # --- Main Query Logic ---
  results_data <- eventReactive(input$run_query, {
    req(input$mz_rt_input)
    shinyjs::show("loading-message")

    # Clear log and add init message
    log_data(c(log_init_msg, sprintf("[%s] Starting query...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))))

    results <- withCallingHandlers({

      # Parse features (handles comma or newline or both)
      features_text <- input$mz_rt_input
      features_split <- strsplit(features_text, "[\n,]")[[1]]
      mz_rt_vector <- trimws(features_split)
      mz_rt_vector <- mz_rt_vector[mz_rt_vector != ""]

      if(length(mz_rt_vector) == 0) {
        stop("No valid features provided. Please check input.", call. = FALSE)
      }

      append_log(sprintf("[%s] Parsed %d features.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(mz_rt_vector)))

      # Get adduct types
      if (input$adduct_input_mode == "select") {
        adducts_vector <- input$adduct_types_select
      } else {
        adduct_text <- input$adduct_types_paste
        adducts_vector <- trimws(strsplit(adduct_text, ",")[[1]])
        adducts_vector <- adducts_vector[adducts_vector != ""]
      }

      # Handle NULL/empty adducts vector
      if (length(adducts_vector) == 0) {
        adducts_vector <- NULL
        append_log(sprintf("[%s] No adduct filter specified. Will return all results.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      } else {
        append_log(sprintf("[%s] Using %d adduct type(s) for filtering: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(adducts_vector), paste(adducts_vector, collapse = ", ")))
      }

      # Get CCS
      ccs_predictors <- if(length(input$ccs_predictors_input) > 0) input$ccs_predictors_input else NULL
      ccs_tolerance <- if(length(input$ccs_tolerance_input) > 0) as.numeric(input$ccs_tolerance_input) else NULL

      if(!is.null(ccs_predictors)) {
        append_log(sprintf("[%s] Using CCS Predictors: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(ccs_predictors, collapse = ", ")))
      }
      if(!is.null(ccs_tolerance)) {
        append_log(sprintf("[%s] Using CCS Tolerance: %s %%", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(ccs_tolerance, collapse = ", ")))
      }

      # Run the query function
      query_hmdb_ms(
        mz_rt = mz_rt_vector,
        ion_mode = input$ion_mode_input,
        filter_adduct_type = adducts_vector,
        tolerance = input$tolerance_input,
        tolerance_units = input$tolerance_units_input,
        ccs_predictors = ccs_predictors,
        ccs_tolerance = ccs_tolerance
      )
    },
    # Capture messages (for logging)
    message = function(m) {
      append_log(conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    # Capture warnings (for logging)
    warning = function(w) {
      append_log(paste("WARNING:", conditionMessage(w)))
      invokeRestart("muffleWarning")
    },
    # Capture errors (for logging)
    error = function(e) {
      append_log(paste("ERROR:", conditionMessage(e)))
      NULL # Return NULL on error
    })

    shinyjs::hide("loading-message")

    if(is.null(results)) {
      append_log(sprintf("[%s] Query failed. See error message above.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    } else {
      append_log(sprintf("[%s] Query finished successfully.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    }

    results
  })

  # --- Render Outputs ---

  output$compound_metadata_table <- DT::renderDT({
    req(results_data()$compound_metadata)

    if(nrow(results_data()$compound_metadata) == 0) {
      return(data.frame(Message = "No compounds found matching your criteria."))
    }

    DT::datatable(results_data()$compound_metadata,
                  options = list(scrollX = TRUE, pageLength = 25, autoWidth = TRUE),
                  rownames = FALSE,
                  selection = 'none'
    )
  })

  output$excluded_results_table <- DT::renderDT({
    req(results_data()$results_table_excluded)

    if(nrow(results_data()$results_table_excluded) == 0) {
      return(data.frame(Message = "No results were excluded by the adduct filter."))
    }

    DT::datatable(results_data()$results_table_excluded,
                  options = list(scrollX = TRUE, pageLength = 25, autoWidth = TRUE),
                  rownames = FALSE,
                  selection = 'none'
    )
  })

  output$query_log <- renderPrint({
    # This will automatically update as log_data() changes
    cat(paste(log_data(), collapse = "\n"))
  })

  # --- Render Downloads ---

  output$download_excel <- downloadHandler(
    filename = function() {
      paste("hmdb_results-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      res_data <- results_data()
      current_log <- log_data()

      if (is.null(res_data) || (nrow(res_data$compound_metadata) == 0 && nrow(res_data$results_table_excluded) == 0)) {
        showNotification("No data available to download.", type = "warning", duration = 5)
        return(NULL)
      }

      # Prepare log data as a data frame
      log_df <- if (length(current_log) > 0) {
        data.frame(Log = current_log, stringsAsFactors = FALSE)
      } else {
        data.frame(Log = "No log data available.", stringsAsFactors = FALSE)
      }

      data_list <- list(
        "Compound_Metadata" = if(nrow(res_data$compound_metadata) == 0) data.frame(Message="No compounds found matching filter") else res_data$compound_metadata,
        "Excluded_Results" = if(nrow(res_data$results_table_excluded) == 0) data.frame(Message="No results were excluded") else res_data$results_table_excluded,
        "Logs" = log_df
      )

      writexl::write_xlsx(data_list, path = file)
    }
  )

}

# ================================================================
# Run the app
# ================================================================
shinyApp(ui, server)
