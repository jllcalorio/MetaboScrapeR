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

# Parallel processing libraries for performance and non-blocking operation
library(future)
library(furrr)
library(later)

# ================================================================
# Helper Function: Fetch Metadata via Web Scraping
# ================================================================

fetch_metadata_via_scraping <- function(compound_ids) {
  if (is.null(compound_ids) || length(compound_ids) == 0) {
    return(data.frame())
  }
  safe_text <- function(page, sel) {
    res <- page %>% rvest::html_element(sel)
    if (length(res) == 0 || is.na(res)) return(NA)
    rvest::html_text(res, trim = TRUE)
  }
  
  # Use an anonymous function inside future_map_dfr to ensure cleanup
  metadata_df <- furrr::future_map_dfr(compound_ids, ~ {
    compound_id <- .x
    url <- paste0("https://www.hmdb.ca/metabolites/", compound_id)
    
    con <- base::url(url, "rb") # Open the connection in binary mode
    on.exit(close(con)) # Use on.exit to ensure the connection is closed even if an error occurs
    
    tryCatch({
      # Use read_html on the connection instead of the URL string
      page <- rvest::read_html(con) 

      data.frame(
        Compound = compound_id,
        URL = url,
        Description = safe_text(page, ".met-desc"),
        Kingdom = safe_text(page, "#taxonomy~ tr:nth-child(24) td"),
        SuperClass = safe_text(page, "tbody:nth-child(1) tr:nth-child(25) td"),
        Class = safe_text(page, "tbody:nth-child(1) tr:nth-child(26) td"),
        Subclass = safe_text(page, "tbody:nth-child(1) tr:nth-child(27) td"),
        Direct_Parent = safe_text(page, "tbody:nth-child(1) tr:nth-child(28) td"),
        Alternative_Parent = safe_text(page, "tr:nth-child(29) ul"),
        Substituents = safe_text(page, "tr:nth-child(30) ul"),
        Molecular_Framework = safe_text(page, "tbody:nth-child(1) tr:nth-child(31) td"),
        External_Descriptors = safe_text(page, "tbody:nth-child(1) tr:nth-child(32) td"),
        Pathway_Name = safe_text(page, "#metabolite-pathway-links th:nth-child(1)"),
        Pathway_SMPDB_PathBank = safe_text(page, "#metabolite-pathway-links th:nth-child(2)"),
        Pathway_KEGG = safe_text(page, "#metabolite-pathway-links th~ th+ th"),
        Disposition = safe_text(page, "#ontology~ tr:nth-child(35) td"),
        Process = safe_text(page, "tr:nth-child(36) td"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning(sprintf("Could not fetch data for HMDB ID %s: %s", compound_id, e$message), call. = FALSE)
      return(NULL)
    })
  }, .options = furrr::furrr_options(seed = NULL))
  return(metadata_df)
}

# ================================================================
# Core Function: Query HMDB (same as original, kept intact)
# ================================================================

query_hmdb_ms <- function(mz_rt,
                          ion_mode,
                          filter_adduct_type,
                          tolerance,
                          tolerance_units,
                          ccs_predictors,
                          ccs_tolerance,
                          log_function) {
  if (is.null(mz_rt) || length(mz_rt) == 0) {
    log_function("ERROR: Argument 'mz_rt' must be provided and non-empty.")
    stop("Argument 'mz_rt' must be provided and non-empty.", call. = FALSE)
  }
  if (length(mz_rt) > 700) {
    log_function("ERROR: The number of features exceeds the HMDB limit of 700. Please reduce your input.")
    stop("The number of features exceeds the HMDB limit of 700.", call. = FALSE)
  }
  process_mz_rt <- function(mz_entry) {
    log_function(sprintf("=== Processing feature: %s ===", mz_entry))
    tryCatch({
      form_values <- list(
        query_masses = as.numeric(stringr::str_extract(mz_entry, "^[^@]+")),
        ms_search_ion_mode = tolower(ion_mode),
        tolerance = tolerance,
        tolerance_units = tolerance_units,
        ccs_predictors = if (is.null(ccs_predictors)) "" else ccs_predictors,
        ccs_tolerance = if (is.null(ccs_tolerance)) "" else as.character(ccs_tolerance),
        "commit" = "Search"
      )
      log_function("Submitting search to HMDB...")
      Sys.sleep(0.5)
      resp <- httr::POST("https://www.hmdb.ca/spectra/ms/search", body = form_values, encode = "form", httr::timeout(60))
      httr::stop_for_status(resp)
      results_page_raw <- httr::content(resp, as = "text")
      results_page_parsed <- rvest::read_html(results_page_raw)
      log_function("Search submitted. Parsing results page...")
      results_table_node <- rvest::html_element(results_page_parsed, ".table")
      if (inherits(results_table_node, "xml_missing")) {
        page_text <- rvest::html_text(results_page_parsed, trim = TRUE)
        if (grepl("No results were found", page_text, ignore.case = TRUE)) {
          log_function(sprintf("No results found for feature '%s'.", mz_entry))
        } else {
          log_function(sprintf("WARNING: Could not find results table for feature '%s'. The page layout may have changed.", mz_entry))
        }
        return(list(compound_metadata = data.frame(), results_table_excluded = data.frame()))
      }
      results_table <- rvest::html_table(results_table_node, fill = TRUE) %>%
        setNames(c("Compound", "Name", "Formula", "Monoisotopic Mass", "Adduct", "Adduct M/Z", "Delta (ppm)", "CCS"))
      log_function(sprintf("Found %d initial matches for feature '%s'.", nrow(results_table), mz_entry))
      results_table$`Adduct M/Z` <- trimws(stringr::str_replace(results_table$`Adduct M/Z`, "m/z calculator", ""))
      excluded <- if (!is.null(filter_adduct_type)) {
        dplyr::filter(results_table, !Adduct %in% filter_adduct_type)
      } else {
        data.frame()
      }
      included <- if (!is.null(filter_adduct_type)) {
        dplyr::filter(results_table, Adduct %in% filter_adduct_type)
      } else {
        results_table
      }
      if (!is.null(filter_adduct_type)) {
        log_function(sprintf("Filtered %d matches by adduct type. %d matches remaining.", nrow(excluded), nrow(included)))
      }
      compound_ids <- included$Compound
      if(length(compound_ids) > 0) {
        log_function(sprintf("Fetching detailed metadata for %d compound(s) via web scraping...", length(compound_ids)))
        Sys.sleep(0.5)
        compound_metadata <- fetch_metadata_via_scraping(compound_ids)
      } else {
        log_function("No compounds remaining after filtering. Skipping metadata fetch.")
        compound_metadata <- data.frame()
      }
      final_metadata <- if (nrow(compound_metadata) > 0) {
        dplyr::left_join(included, compound_metadata, by = "Compound")
      } else if (nrow(included) > 0) {
        included
      } else {
        data.frame()
      }
      query_meta <- data.frame(
        Ion = ion_mode,
        MZ = form_values$query_masses,
        RT = as.numeric(stringr::str_extract(mz_entry, "[^@]+$")),
        stringsAsFactors = FALSE
      )
      if (nrow(final_metadata) > 0) {
        final_metadata <- cbind(query_meta, final_metadata)
      }
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
      log_function(sprintf("=== Finished processing feature: %s. Found %d valid compounds. ===\n", mz_entry, nrow(final_metadata)))
      return(list(compound_metadata = final_metadata, results_table_excluded = excluded))
    }, error = function(e) {
      log_function(sprintf("ERROR: Failed to process feature '%s': %s", mz_entry, e$message))
      return(list(compound_metadata = data.frame(), results_table_excluded = data.frame()))
    })
  }
  results_list <- lapply(mz_rt, process_mz_rt)
  compound_metadata_all <- do.call(rbind, lapply(results_list, `[[`, "compound_metadata"))
  results_excluded_all <- do.call(rbind, lapply(results_list, `[[`, "results_table_excluded"))
  total_found <- nrow(compound_metadata_all) + nrow(results_excluded_all)
  log_function(sprintf("=== Total compounds found before adduct filtering: %d ===", total_found))
  n_excluded <- nrow(results_excluded_all)
  if (n_excluded > 0) {
    log_function(sprintf("Excluded %d compound(s) based on adduct filter. %d remaining.", n_excluded, nrow(compound_metadata_all)))
  } else {
    log_function("No compounds were excluded by the adduct filter.")
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
  tags$head(
    tags$style(
      HTML(
        "
        /* Make sidebar scrollable independently */
        .sidebar-container {
          position: fixed;
          top: 60px;
          left: 0;
          width: 33.33%;
          height: calc(100vh - 60px);
          overflow-y: auto;
          padding: 15px;
          background-color: #f5f5f5;
          border-right: 1px solid #ddd;
        }
        /* Make main panel scrollable independently */
        .main-container {
          margin-left: 33.33%;
          padding: 15px;
          height: calc(100vh - 60px);
          overflow-y: auto;
        }
        /* Spinner styling */
        .spinner {
          border: 5px solid #f3f3f3;
          border-radius: 50%;
          border-top: 5px solid #3498db;
          width: 50px;
          height: 50px;
          -webkit-animation: spin 2s linear infinite;
          animation: spin 2s linear infinite;
          margin: 0 auto 10px auto;
        }
        @-webkit-keyframes spin {
          0% { -webkit-transform: rotate(0deg); }
          100% { -webkit-transform: rotate(360deg); }
        }
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
        /* DataTables styles for column width/scroll */
        .dataTables_scrollBody {
          overflow-x: auto !important;
        }
        .dataTables_scrollHeadInner, .dataTable {
          width: 100% !important;
        }
        /* Help tab styling */
        .help-content { padding: 10px; }
        .help-content h3 { border-bottom: 2px solid #eee; padding-bottom: 5px; margin-top: 20px;}
        .help-content h4 { margin-top: 15px; }
        .help-content code { background-color: #f5f5f5; border-radius: 3px; padding: 2px 4px; }
        /* Force inputs to 100% width */
        .shiny-input-container {
            width: 100% !important;
        }
        /* Author info styling */
        .author-info {
          font-size: 0.9em;
          line-height: 1.6;
        }
        .author-info a {
          color: #337ab7;
          text-decoration: none;
        }
        .author-info a:hover {
          text-decoration: underline;
        }
        /* Prevent horizontal resize on specific textareas, allow vertical */
        #mz_rt_input, #adduct_types_paste {
          resize: vertical; /* Allow vertical resize, but not horizontal */
          width: 100% !important; /* Force to 100% of parent */
        }
        /* Note box styling for Help */
        .note-box {
          border-left: 4px solid #2b7cff;
          background: #f0f8ff;
          padding: 10px;
          margin-top: 12px;
          border-radius: 4px;
        }
        "
      )
    )
  ),
  
  # Sidebar Panel (Scrollable)
  div(class = "sidebar-container",
      helpText("Enter m/z and RT values to query the HMDB."),
      textAreaInput("mz_rt_input",
                    "Features (m/z@RT)",
                    value = "",
                    rows = 4,
                    placeholder = "Example: 626.3547@95.2\n515.2597@10.7733",
                    width = '100%'),
      helpText("Separate multiple features with newlines, commas, or both. Maximum 700 features allowed."),
      radioButtons("ion_mode_input",
                   "Ion Mode",
                   choices = c("Positive", "Negative", "Neutral"),
                   selected = "Positive",
                   inline = TRUE,
                   width = '100%'),
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
          NULL,
          choices = c(
            "2M+2H+3H2O", "2M+ACN+H", "2M+ACN+Na", "2M+FA-H", "2M+H", "2M+H-H2O",
            "2M+K", "2M+Na", "2M-H", "3M-H", "M", "M+2ACN+2H", "M+2ACN+H", "M+2H",
            "M+2H+Na", "M+2K-H", "M+2Na", "M+2Na-H", "M+3ACN+2H", "M+3H",
            "M+3Na", "M+ACN+2H", "M+ACN+H", "M+ACN+Na", "M+Br", "M+CH3OH+H",
            "M+Cl", "M+DMSO+H", "M+F", "M+FA-H", "M+H", "M+H+2K", "M+H+2Na",
            "M+H+HCOONa", "M+H+K", "M+H+Na", "M+H-2H2O", "M+H-H2O", "M+Hac-H",
            "M+IsoProp+H", "M+K", "M+K-2H", "M+Li", "M+NH4", "M+NH4-H2O", "M+Na",
            "M+Na-2H", "M+TFA-H", "M-2H", "M-3H", "M-H", "M-H+HCOONa", "M-H20-H",
            "Unknown"
          ),
          selected = NULL,
          multiple = TRUE,
          options = list(placeholder = 'Select adducts...'),
          width = '100%'
        )
      ),
      conditionalPanel(
        condition = "input.adduct_input_mode == 'paste'",
        textAreaInput(
          "adduct_types_paste",
          NULL,
          value = "",
          rows = 3,
          placeholder = "Example: M+H, M+Na, M+K",
          width = '100%'
        )
      ),
      helpText("See Help section for the full list of adducts."),
      radioButtons("tolerance_units_input",
                   "Tolerance Units",
                   choices = c("ppm", "Da"),
                   selected = "ppm",
                   inline = TRUE,
                   width = '100%'),
      numericInput("tolerance_input",
                   "Mass Tolerance Value",
                   value = 5,
                   min = 0,
                   step = 1,
                   width = '100%'),
      radioButtons("ccs_predictors_input",
                   "CCS Predictor (Optional)",
                   choices = c("None" = "", "AllCCS", "DarkChem", "DeepCCS"),
                   selected = "",
                   inline = TRUE,
                   width = '100%'),
      radioButtons("ccs_tolerance_input",
                   "CCS Tolerance (Optional)",
                   choices = c("None" = "", "1", "3", "5", "10"),
                   selected = "",
                   inline = TRUE,
                   width = '100%'),
      div(style = "display: flex; justify-content: space-between; gap: 10px;",
          actionButton("run_query", "Run Query", class = "btn-primary", icon = icon("play"), style = "width: 100%;")
      ),
      hr(),
      downloadButton("download_excel", "Download Results (.xlsx)", class = "btn-success", style = "width: 100%;"),
      hr(),
      div(
        class = "well author-info",
        HTML('<p><strong>Source code on GitHub:</strong> <a href="https://github.com/jllcalorio/MetaboScrapeR" target="_blank">MetaboScrapeR</a></p>
              <p><strong>Author:</strong> John Lennon L. Calorio</p>
              <p><strong>Email:</strong> <a href="mailto:jllcalorio@gmail.com">jllcalorio@gmail.com</a></p>
              <p><strong>GitHub:</strong> <a href="https://github.com/jllcalorio/" target="_blank">github.com/jllcalorio/</a></p>
              <p><strong>LinkedIn:</strong> <a href="https://linkedin.com/in/caloriojohnlennon/" target="_blank">linkedin.com/in/caloriojohnlennon/</a></p>
              <p><strong>Instagram:</strong> <a href="https://instagram.com/johnlennoncalorio/" target="_blank">instagram.com/johnlennoncalorio/</a></p>
              <p><strong>Google Scholar:</strong> <a href="https://scholar.google.com/citations?user=5rBOxykAAAAJ&hl=en" target="_blank">scholar.google.com/citations?user=5rBOxykAAAAJ&hl=en</a></p>')
      )
  ),
  
  # Main Panel (Scrollable independently)
  div(class = "main-container",
      shiny::fluidRow(
        shiny::column(width = 12,
                      align = "center",
                      div(id = "loading-message",
                          HTML("<div class='spinner'></div>"),
                          tags$p("Querying HMDB... This may take several minutes per feature. Do not close the tab."),
                          style = "display: none; color: #555; font-style: italic; margin-top: 15px;")
        )
      ),
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Compound Metadata",
                 h4("Filtered Compound Metadata (Draggable columns, sortable)"),
                 DT::DTOutput("compound_metadata_table")),
        tabPanel("Excluded Results",
                 h4("Excluded Results Table (Adducts Not Filtered)"),
                 DT::DTOutput("excluded_results_table")),
        tabPanel("Help",
                 div(class = "help-content",
                     h3("Overview"),
                     p("This app queries the HMDB (Human Metabolome Database) MS search interface (",
                       tags$a(href = "https://www.hmdb.ca/spectra/ms/search", target = "_blank", "https://www.hmdb.ca/spectra/ms/search"),
                       ") to retrieve compound information based on mass spectrometry features. It uses web scraping for reliable data retrieval."),
                     
                     h3("Input Format"),
                     h4("Features (m/z@RT)"),
                     p("Enter your features in the format: ", code("MZ@RT"), " where:"),
                     tags$ul(
                       tags$li(code("MZ"), " is the mass-to-charge ratio (e.g., 626.3547)"),
                       tags$li(code("RT"), " is the retention time (e.g., 95.2)"),
                       tags$li("Use ", code("@"), " (single at sign) as the separator"),
                       tags$li("Separate multiple features with newlines, commas, or both"),
                       tags$li(strong("Maximum 700 features allowed per query"))
                     ),
                     p("Examples:"),
                     tags$pre("626.3547@95.2\n615.3753@120.5\n450.2315@80.0"),
                     p("OR:"),
                     tags$pre("626.3547@95.2, 615.3753@120.5, 450.2315@80.0"),
                     p("OR:"),
                     tags$pre("626.3547@95.2, 615.3753@120.5\n450.2315@80.0"),
                     
                     # Ion Mode description and Adducts to Keep
                     h4("Ion Mode"),
                     p("Select the ionization mode used in your experiment. This will automatically update the default 'Adducts to Keep' list."),
                     
                     h4("Adducts to Keep"),
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
                               )),
                       tags$li(strong("Neutral Mode:"),
                               tags$ul(
                                 tags$li("Unknown, M")
                               ))
                     ),
                     p("The default adducts selected by this app when you change Ion Mode are a curated list of the most common adducts from the lists above."),
                     
                     h3("R Packages Used"),
                     p("This application uses the following R packages:"),
                     tags$ul(
                       tags$li(code("shiny"), " - Web application framework"),
                       tags$li(code("shinyjs"), " - JavaScript operations in Shiny"),
                       tags$li(code("DT"), " - Interactive DataTables"),
                       tags$li(code("writexl"), " - Excel file export"),
                       tags$li(code("rvest"), " - Web scraping"),
                       tags$li(code("stringr"), " - String manipulation"),
                       tags$li(code("dplyr"), " - Data manipulation"),
                       tags$li(code("purrr"), " - Functional programming tools"),
                       tags$li(code("httr"), " - HTTP requests"),
                       tags$li(code("jsonlite"), " - JSON parsing"),
                       tags$li(code("future"), " - Parallel processing framework"),
                       tags$li(code("furrr"), " - Parallel mapping with purrr"),
                       tags$li(code("later"), " - Deferred execution")
                     ),
                     h3("Output Tabs"),
                     tags$ul(
                       tags$li(strong("Compound Metadata:"), " Shows compounds that matched your search criteria AND your selected adducts. Includes detailed metadata from HMDB."),
                       tags$li(strong("Excluded Results:"), " Shows compounds that matched your search criteria but were filtered out because their adducts did not match your 'Adducts to Keep' list."),
                       tags$li(strong("Logs:"), " Shows a detailed, timestamped log of the entire query process, updated in real-time.")
                     ),
                     
                     # Note-box
                     div(class = "note-box",
                         tags$p(style = "margin: 0;",
                                icon("info-circle"), " ", strong("Need Additional Features?"),
                                tags$br(),
                                p("If you need additional columns extracted from HMDB, have suggestions for improvements, or encounter any issues, please feel free to reach out! You can contact me using the details provided in the sidebar (GitHub, Email, LinkedIn, etc.). If you want to suggest and make revisions yourself, make full use of", tags$a(href = "https://github.com/jllcalorio/MetaboScrapeR/issues", target = "_blank", "GitHub's 'Issues' feature here"), ".")
                         )
                     )
                 )
        ),
        tabPanel("Logs",
                 h4("Query Log (Real-time update)"),
                 verbatimTextOutput("query_log"))
      )
  )
)

# ================================================================
# Server
# ================================================================
server <- function(input, output, session) {
  workers_available <- future::availableCores() - 1
  if (workers_available < 1) workers_available <- 1
  log_init_msg <- sprintf("[%s] Parallel backend initializing with %d workers...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), workers_available)
  future::plan(future::multisession, workers = workers_available)
  
  rv <- reactiveValues(
    is_running = FALSE,
    results_data = NULL,
    log_data = character(0),
    query_successful = FALSE,
    query_error = FALSE
  )
  
  append_log <- function(msg) {
    current_log <- rv$log_data
    rv$log_data <- c(current_log, msg)
  }
  
  positive_adducts <- c("M+H", "M+H-2H2O", "M+H-H2O", "M+Na", "M+CH3OH+H", "M+K", "M+ACN+H", "M+2Na-H", "M+ACN+Na", "M+2K-H", "M+2ACN+H", "M+H+HCOONa", "2M+H", "2M+Na", "2M+2H+3H2O", "2M+K", "2M+ACN+H", "2M+ACN+Na", "2M+H-H2O", "M+2H", "M+H+Na", "M+H+K", "M+ACN+2H", "M+2Na", "M+2ACN+2H", "M+3ACN+2H", "M+3H", "M+2H+Na", "M+H+2Na", "M+3Na", "M+H+2K")
  negative_adducts <- c("M-H", "M-H20-H", "M+Na-2H", "M+Cl", "M+K-2H", "M+FA-H", "M-H+HCOONa", "2M-H", "2M+FA-H", "3M-H", "M-2H", "M-3H")
  neutral_adducts <- c("Unknown", "M")
  
  observeEvent(input$ion_mode_input, {
    ion_mode <- input$ion_mode_input
    if (ion_mode == "Positive") {
      defaults <- positive_adducts
    } else if (ion_mode == "Negative") {
      defaults <- negative_adducts
    } else if (ion_mode == "Neutral") {
      defaults <- neutral_adducts
    } else {
      defaults <- ""
    }
    updateSelectizeInput(session, "adduct_types_select",
                         selected = defaults)
    updateTextAreaInput(session, "adduct_types_paste",
                        value = paste(defaults, collapse = ", "))
    updateRadioButtons(session, "adduct_input_mode",
                       selected = "select")
  }, ignoreInit = FALSE)
  
  observeEvent(input$adduct_input_mode, {
    if (input$adduct_input_mode == "paste") {
      selected_adducts <- input$adduct_types_select
      if (length(selected_adducts) > 0) {
        updateTextAreaInput(session, "adduct_types_paste",
                            value = paste(selected_adducts, collapse = ", "))
      }
    } else {
      paste_text <- input$adduct_types_paste
      if (!is.null(paste_text) && nchar(trimws(paste_text)) > 0) {
        adducts_vector <- trimws(strsplit(paste_text, ",")[[1]])
        adducts_vector <- adducts_vector[adducts_vector != ""]
        updateSelectizeInput(session, "adduct_types_select",
                             selected = adducts_vector)
      }
    }
  })
  
  # --- Dynamic tolerance value and step based on units ---
  observeEvent(input$tolerance_units_input, {
    if (input$tolerance_units_input == "ppm") {
      updateNumericInput(session, "tolerance_input", value = 5, step = 1)
    } else if (input$tolerance_units_input == "Da") {
      updateNumericInput(session, "tolerance_input", value = 0.05, step = 0.01)
    }
  })
  
  # --- Main Query Logic (Blocking) ---
  results_data <- eventReactive(input$run_query, {
    
    # 1. Setup State and UI
    req(input$mz_rt_input)
    shinyjs::show("loading-message")
    rv$is_running <- TRUE # Set immediately
    rv$query_successful <- FALSE
    rv$query_error <- FALSE
    rv$results_data <- NULL # Clear previous results
    rv$log_data <- c(log_init_msg, sprintf("[%s] Starting query...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    # 2. Prepare arguments
    features_text <- input$mz_rt_input
    features_split <- strsplit(features_text, "[\n,]")[[1]]
    mz_rt_vector <- trimws(features_split)
    mz_rt_vector <- mz_rt_vector[mz_rt_vector != ""]
    
    if(length(mz_rt_vector) == 0) {
      append_log(sprintf("[%s] ERROR: No valid features provided.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      rv$is_running <- FALSE
      shinyjs::hide("loading-message")
      return(NULL)
    }
    if (length(mz_rt_vector) > 700) {
      append_log(sprintf("[%s] ERROR: The number of features exceeds the HMDB limit of 700. Please reduce your input.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      rv$is_running <- FALSE
      shinyjs::hide("loading-message")
      return(NULL)
    }
    append_log(sprintf("[%s] Parsed %d features.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(mz_rt_vector)))
    
    adducts_vector <- if (input$adduct_input_mode == "select") {
      input$adduct_types_select
    } else {
      adduct_text <- input$adduct_types_paste
      trimws(strsplit(adduct_text, ",")[[1]])
    }
    adducts_vector <- adducts_vector[adducts_vector != ""]
    if (length(adducts_vector) == 0) {
      adducts_vector <- NULL
      append_log(sprintf("[%s] No adduct filter specified. Will return all results.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    } else {
      append_log(sprintf("[%s] Using %d adduct type(s) for filtering: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(adducts_vector), paste(adducts_vector, collapse = ", ")))
    }
    
    ccs_predictors <- if(input$ccs_predictors_input == "") NULL else input$ccs_predictors_input
    ccs_tolerance <- if(input$ccs_tolerance_input == "") NULL else as.numeric(input$ccs_tolerance_input)
    if(!is.null(ccs_predictors)) {
      append_log(sprintf("[%s] Using CCS Predictor: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ccs_predictors))
    }
    if(!is.null(ccs_tolerance)) {
      append_log(sprintf("[%s] Using CCS Tolerance: %s %%", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ccs_tolerance))
    }
    
    # 3. Blocking Query Execution
    res <- tryCatch({
      query_hmdb_ms(
        mz_rt = mz_rt_vector,
        ion_mode = input$ion_mode_input,
        filter_adduct_type = adducts_vector,
        tolerance = input$tolerance_input,
        tolerance_units = input$tolerance_units_input,
        ccs_predictors = ccs_predictors,
        ccs_tolerance = ccs_tolerance,
        # Use real-time log function for blocking query
        log_function = append_log 
      )
    }, error = function(e) {
      append_log(sprintf("[%s] ERROR during query execution: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), conditionMessage(e)))
      NULL
    })
    
    # 4. Final UI cleanup and state update
    rv$is_running <- FALSE # Set when query is finished
    shinyjs::hide("loading-message")
    
    # 5. Process results and return
    if (is.null(res)) {
      rv$query_error <- TRUE
      rv$query_successful <- FALSE
      return(NULL)
    } else {
      append_log(sprintf("[%s] Query finished successfully. Found %d compound metadata rows and %d excluded result rows.", 
                         format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
                         nrow(res$compound_metadata), 
                         nrow(res$results_table_excluded)))
      rv$query_successful <- TRUE
      rv$results_data <- res # Set the final results
      rv$query_error <- FALSE
      return(res) # return results to results_data()
    }
  })
  
  observeEvent(input$stop_query, {
    if (rv$is_running) {
      rv$stop_requested <- TRUE
      append_log("User requested query interruption. Resetting UI state.")
      shinyjs::hide("loading-message")
      rv$is_running <- FALSE
    }
  })
  
  # --- UI State Observer (Disable/Enable buttons) ---
  observe({
    has_features <- !is.null(input$mz_rt_input) && nchar(trimws(input$mz_rt_input)) > 0
    
    # Only toggle the run button based on running state and feature input
    shinyjs::toggleState("run_query", !rv$is_running && has_features)
  })
  
  # --- Render Outputs ---
  output$query_log <- renderPrint({
    cat(paste(rv$log_data, collapse = "\n"))
  })
  
  output$compound_metadata_table <- DT::renderDT({
    req(results_data())
    if(is.null(results_data()$compound_metadata) || nrow(results_data()$compound_metadata) == 0) {
      return(data.frame(Message = "No compounds found matching your criteria."))
    }
    DT::datatable(results_data()$compound_metadata,
                  options = list(
                    scrollX = TRUE,
                    pageLength = 25,
                    autoWidth = TRUE,
                    columnDefs = list(list(className = 'dt-head-center', targets = '_all'))
                  ),
                  rownames = FALSE,
                  width = 150,
                  height = 50,
                  selection = 'none'
    )
  })
  
  output$excluded_results_table <- DT::renderDT({
    req(results_data())
    if(is.null(results_data()$results_table_excluded) || nrow(results_data()$results_table_excluded) == 0) {
      return(data.frame(Message = "No results were excluded by the adduct filter."))
    }
    DT::datatable(results_data()$results_table_excluded,
                  options = list(
                    scrollX = TRUE,
                    pageLength = 25,
                    autoWidth = TRUE,
                    columnDefs = list(list(className = 'dt-head-center', targets = '_all'))
                  ),
                  rownames = FALSE,
                  width = 150,
                  height = 50,
                  selection = 'none'
    )
  })
  
  # --- Render Downloads ---
  output$download_excel <- downloadHandler(
    filename = function() {
      paste("hmdb_results-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      res_data <- results_data()
      current_log <- rv$log_data
      if (is.null(res_data) || (nrow(res_data$compound_metadata) == 0 && nrow(res_data$results_table_excluded) == 0)) {
        showNotification("No data available to download.", type = "warning", duration = 5)
        return(NULL)
      }
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
