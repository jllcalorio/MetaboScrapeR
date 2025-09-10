#' Query HMDB MS Spectra and Retrieve Compound Metadata
#'
#' This function queries the Human Metabolome Database (HMDB) Mass Spectrometry
#' (MS) search interface using provided m/z and retention time (RT), ion mode,
#' adduct filters, and optional collision cross section (CCS) parameters.
#' It returns both detailed compound metadata from individual HMDB metabolite
#' pages and the excluded results table (those adducts not in the filter list).
#'
#' @param mz_rt Character string or vector. Feature m/z and RT in the format
#'   `"m/z@RT"`. Defaults to `NULL`. Can accept multiple values.
#' @param ion_mode Character string. MS ionization mode. One of
#'   `"Positive"`, `"Negative"`, or `"Neutral"`. Default = `"Positive"`.
#' @param filter_adduct_type Character vector. Specific adducts to retain.
#'   Defaults to `NULL`. If provided, results are filtered accordingly.
#' @param tolerance Numeric. Mass tolerance value.
#'   Default = `5`.
#' @param tolerance_units Character string. Mass tolerance unit.
#'   One of `"Da"` or `"ppm"`. Default = `"ppm"`.
#' @param ccs_predictors Character string or vector. CCS prediction model(s).
#'   Examples: `"AllCCS"`, `"DarkChem"`, `"DeepCCS"`. Default = `NULL`.
#' @param ccs_tolerance Numeric. CCS tolerance value. Default = `NULL`.
#'
#' @return A list with two data frames:
#'   \itemize{
#'     \item \code{compound_metadata}: Compound details scraped from individual
#'       HMDB metabolite pages (name, taxonomy, framework, etc.).
#'     \item \code{results_table_excluded}: Table of results where adducts
#'       do not match the provided filter list.
#'   }
#'
#' @details
#' The function submits form data to the HMDB MS Spectra Search interface,
#' retrieves the results table, filters adducts (if specified), and then
#' scrapes additional information for each compound from its dedicated page.
#'
#' @importFrom rvest read_html html_form html_form_set html_form_submit html_element html_table html_text
#' @importFrom stringr str_extract
#' @importFrom dplyr filter
#' @importFrom purrr map_df
#' @importFrom httr content
#'
#' @examples
#' \dontrun{
#' results <- query_hmdb_ms(
#'   mz_rt = "515.2597@10.7733",
#'   ion_mode = "Positive",
#'   filter_adduct_type = c("M+H", "M+Na"),
#'   tolerance = 5,
#'   tolerance_units = "ppm"
#' )
#'
#' head(results$compound_metadata)
#' head(results$results_table_excluded)
#' }
#'
#' @seealso \url{https://www.hmdb.ca/spectra/ms/search}
#' @author Your Name
#' @export
query_hmdb_ms <- function(mz_rt = NULL,
                          ion_mode = "Positive",
                          filter_adduct_type = NULL,
                          tolerance = 5,
                          tolerance_units = "ppm",
                          ccs_predictors = NULL,
                          ccs_tolerance = NULL) {

  # Force messages to output in console
  utils::flush.console()

  # Load packages safely
  if (!requireNamespace("rvest", quietly = TRUE)) stop("Package 'rvest' required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' required.")
  if (!requireNamespace("httr", quietly = TRUE)) stop("Package 'httr' required.")
  if (!requireNamespace("purrr", quietly = TRUE)) stop("Package 'purrr' required.")

  # Early exit if no mz_rt is given
  if (is.null(mz_rt)) stop("Argument 'mz_rt' must be provided.")

  # Check if a container for selenium/standalone-chrome is already running
  # We use system() to execute a shell command and capture its output
  container_id <- system("docker run -d -p 4444:4444 -p 5900:5900 selenium/standalone-chrome", intern = TRUE)

  # The intern = TRUE argument captures the output as a character vector.
  # If the vector is empty, no running container was found.
  if (length(container_id) == 0) {
    # No running container found, so run the command
    system("docker run -d -p 4444:4444 -p 5900:5900 selenium/standalone-chrome")
    print("Started a new selenium/standalone-chrome container.")
  } else {
    # A container is already running
    print("A selenium/standalone-chrome container is already running.")
  }


  # Read website and extract form
  static <- rvest::read_html("https://www.hmdb.ca/spectra/ms/search")
  form <- rvest::html_form(static)[[2]]

  # Helper: extract metadata for one compound
  get_compound_info <- function(compound_id) {
    url <- paste0("https://www.hmdb.ca/metabolites/", compound_id)
    page <- tryCatch(rvest::read_html(url), error = function(e) return(NULL))
    if (is.null(page)) return(NULL)

    safe_text <- function(sel) {
      res <- page %>% rvest::html_element(sel)
      if (length(res) == 0 || is.na(res)) return(NA)
      rvest::html_text(res, trim = TRUE)
    }

    # The ones in double quotation marks are CSS selectors found using the SelectorGadget tool
    data.frame(
      Compound = compound_id,
      URL = url,
      Description = safe_text(".met-desc"),
      Kingdom = safe_text("#taxonomy~ tr:nth-child(24) td"),
      SuperClass = safe_text("tbody:nth-child(1) tr:nth-child(25) td"),
      Class = safe_text("tbody:nth-child(1) tr:nth-child(26) td"),
      Subclass = safe_text("tbody:nth-child(1) tr:nth-child(27) td"),
      Direct_Parent = safe_text("tbody:nth-child(1) tr:nth-child(28) td"),
      Alternative_Parent = safe_text("tr:nth-child(29) ul"),
      Substituents = safe_text("tr:nth-child(30) ul"),
      Molecular_Framework = safe_text("tbody:nth-child(1) tr:nth-child(31) td"),
      External_Descriptors = safe_text("tbody:nth-child(1) tr:nth-child(32) td"),
      Disposition = safe_text("#ontology~ tr:nth-child(35) td"),
      stringsAsFactors = FALSE
    )
  }

  # Process each mz_rt input
  results_list <- lapply(mz_rt, function(mz_entry) {
    mz <- stringr::str_extract(mz_entry, "^[^@]+") %>% as.numeric()
    rt <- stringr::str_extract(mz_entry, "[^@]+$") %>% as.numeric()

    form_set <- rvest::html_form_set(
      form,
      query_masses = mz,
      ms_search_ion_mode = tolower(ion_mode),
      tolerance = tolerance,
      tolerance_units = tolerance_units,
      ccs_predictors = ccs_predictors,
      ccs_tolerance = ccs_tolerance
    )

    form_submit <- rvest::html_form_submit(form_set, submit = "commit")
    results_page <- httr::content(form_submit, as = "parsed")

    results_table <- tryCatch({
      results_page %>%
        rvest::html_element("table") %>%
        rvest::html_table(fill = TRUE) %>%
        as.data.frame() %>%
        `colnames<-`(., c("Compound", "Name", "Formula", "Monoisotopic Mass",
                          "Adduct", "Adduct M/Z", "Delta (ppm)", "CCS"))
    }, error = function(e) NULL)

    if (is.null(results_table)) return(NULL)

    # Filter excluded
    results_table_excluded <- if (!is.null(filter_adduct_type)) {
      dplyr::filter(results_table, !Adduct %in% filter_adduct_type)
    } else {
      results_table
    }

    # Scrape compound metadata
    compound_ids <- results_table$Compound
    compound_metadata <- purrr::map_df(compound_ids, get_compound_info)

    results <- list(compound_metadata = compound_metadata,
                    results_table_excluded = results_table_excluded)

    return(results)
  })

  # Combine results from multiple mz_rt values
  compound_metadata_all <- do.call(rbind, lapply(results_list, `[[`, "compound_metadata"))
  results_excluded_all <- do.call(rbind, lapply(results_list, `[[`, "results_table_excluded"))


  # # Number of rows of original results table
  nrow <- dim(compound_metadata_all)[1]

  message(sprintf("There are %d compound/s found before filtering of adducts was applied.", nrow))

  # # Number of rows of filtered results table
  nrow_filtered <- dim(results_excluded_all)[1]

  message(sprintf("There were %d compound/s excluded after filtering of adducts was applied. There are now only %d compounds remaining.", nrow_filtered, nrow - nrow_filtered))

  return(list(compound_metadata = compound_metadata_all,
              results_table_excluded = results_excluded_all))
}
