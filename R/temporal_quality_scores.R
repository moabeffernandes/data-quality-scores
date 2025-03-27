# Temporal Quality Scores


# Calculate Temporal Recency Quality Score (TQS) ----
#
# This function calculates a normalised "recency index" (TRQS) for each taxon,
# based on how recent the records are relative to a minimum year baseline.
#
# Arguments:
# - data: A data frame with observation years and taxon names
# - taxon_column: Name of the column with taxon names
# - year_column: Name of the column with year of observation
# - min_year_override: Optional fixed minimum year for normalisation
# - by_species: If TRUE, computes normalisation per species
#
# Returns:
# - A data frame with taxon and its TQS score (0–1)

score_temporal_recency <- function(data,
                                   taxon_column = "taxon_name",
                                   year_column = "year",
                                   min_year_override = NULL,
                                   by_species = FALSE) {
  
  # Validate required columns
  required_cols <- c(taxon_column, year_column)
  if (!all(required_cols %in% names(data))) {
    stop("Missing required columns: ", paste(required_cols, collapse = ", "))
  }
  
  current_year <- as.numeric(format(Sys.Date(), "%Y"))
  
  # Clean and prepare data
  cleaned_data <- data %>%
    dplyr::mutate(observed_year = as.numeric(.data[[year_column]])) %>%
    dplyr::filter(!is.na(observed_year) & observed_year <= current_year) %>%
    dplyr::mutate(observed_date = as.Date(paste0(observed_year, "-01-01")))
  
  if (nrow(cleaned_data) == 0) {
    stop("No valid records after filtering future or missing years.")
  }
  
  # Compute record age in years
  cleaned_data <- cleaned_data %>%
    dplyr::mutate(record_age_years = as.numeric(difftime(Sys.Date(), observed_date, units = "days")) / 365.25)
  
  ### User-specified minimum year ###
  if (!is.null(min_year_override)) {
    custom_min_year <- as.numeric(min_year_override)
    earliest_year <- min(cleaned_data$observed_year, na.rm = TRUE)
    
    if (custom_min_year > earliest_year) {
      stop("min_year_override (", custom_min_year, ") is later than the earliest year in the data (", earliest_year, ").")
    }
    
    max_possible_age <- as.numeric(difftime(Sys.Date(), as.Date(paste0(custom_min_year, "-01-01")), units = "days")) / 365.25
    
    result <- cleaned_data %>%
      dplyr::mutate(taxon = .data[[taxon_column]]) %>%
      dplyr::group_by(taxon) %>%
      dplyr::summarise(
        Recency_Index = mean(record_age_years, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(TQS = Recency_Index / max_possible_age)
    
    ### Normalise per species ###
  } else if (by_species) {
    min_years_by_taxon <- cleaned_data %>%
      dplyr::mutate(taxon = .data[[taxon_column]]) %>%
      dplyr::group_by(taxon) %>%
      dplyr::summarise(min_year = min(observed_year, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(max_possible_age = as.numeric(difftime(Sys.Date(), as.Date(paste0(min_year, "-01-01")), units = "days")) / 365.25)
    
    result <- cleaned_data %>%
      dplyr::mutate(taxon = .data[[taxon_column]]) %>%
      dplyr::group_by(taxon) %>%
      dplyr::summarise(
        Recency_Index = mean(record_age_years, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::left_join(min_years_by_taxon, by = "taxon") %>%
      dplyr::mutate(TQS = Recency_Index / max_possible_age)
    
    ### Global minimum year ###
  } else {
    global_min_year <- suppressWarnings(min(cleaned_data$observed_year, na.rm = TRUE))
    
    if (is.na(global_min_year) || is.infinite(global_min_year)) {
      stop("No valid minimum year found in the dataset.")
    }
    
    max_possible_age <- as.numeric(difftime(Sys.Date(), as.Date(paste0(global_min_year, "-01-01")), units = "days")) / 365.25
    
    result <- cleaned_data %>%
      dplyr::mutate(taxon = .data[[taxon_column]]) %>%
      dplyr::group_by(taxon) %>%
      dplyr::summarise(
        Recency_Index = mean(record_age_years, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(TRQS = Recency_Index / max_possible_age)
  }
  
  # Rename taxon column to match original input
  result <- result %>%
    dplyr::select(!!taxon_column := taxon, TRQS)
  
  # Return TRQS
  return(result)
  
}



# Calculate Logistic Temporal Quality Score (LTQS) ----
#
# Uses a logistic decay function to weight recent records more heavily than old ones.
#
# Arguments:
# - data: Data frame with observation year and taxon columns
# - taxon_column: Name of the taxon column (default: "taxon_name")
# - year_column: Name of the observation year column
# - max_year: Optional maximum year to use for age calculation (default = current year)
# - by_species: Normalise logistic curve per species (TRUE) or globally (FALSE)
# - by_record: Return score per record (TRUE) or species mean (FALSE)
#
# Returns:
# - A tibble with either per-record or per-species LTQS scores

score_logistic_temporal_quality <- function(data,
                                            taxon_column = "taxon_name",
                                            year_column = "year",
                                            max_year = NULL,
                                            by_species = FALSE,
                                            by_record = FALSE) {
  
  # Validate columns
  required_cols <- c(taxon_column, year_column)
  if (!all(required_cols %in% names(data))) {
    stop("Missing required columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Define year for recency calculation
  year_now <- if (!is.null(max_year)) {
    as.numeric(max_year)
  } else {
    as.numeric(format(Sys.Date(), "%Y"))
  }
  
  # Clean and prepare year data
  cleaned_data <- data %>%
    dplyr::mutate(
      observed_year = as.numeric(.data[[year_column]]),
      year_diff = year_now - observed_year
    ) %>%
    dplyr::filter(!is.na(year_diff) & year_diff >= 0)
  
  if (nrow(cleaned_data) == 0) {
    stop("No valid year data after filtering (check for future or missing years).")
  }
  
  if (by_species) {
    # Normalise per species (each species gets its own logistic curve)
    cleaned_data <- cleaned_data %>%
      dplyr::mutate(taxon = .data[[taxon_column]]) %>%
      dplyr::group_by(taxon) %>%
      dplyr::mutate(
        h = max(year_diff, na.rm = TRUE),  # species-specific max age
        LTQS = (1 - ((year_diff / h)^2))^2
      ) %>%
      dplyr::ungroup()
  } else {
    # Global normalisation
    h_global <- max(cleaned_data$year_diff, na.rm = TRUE)
    
    cleaned_data <- cleaned_data %>%
      dplyr::mutate(
        taxon = .data[[taxon_column]],
        LTQS = (1 - ((year_diff / h_global)^2))^2
      )
  }
  
  if (by_record) {
    # Return individual scores per record
    result <- cleaned_data %>%
      dplyr::select(!!taxon_column := taxon, LTQS)
  } else {
    # Return average score per species
    result <- cleaned_data %>%
      dplyr::group_by(taxon) %>%
      dplyr::summarise(
        LTQS = mean(LTQS, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::select(!!taxon_column := taxon, LTQS)
  }
  
  # Return LTQS
  return(result)
  
}



# Calculate Temporal Spread Score (TSS) ----
#
# This function evaluates how evenly distributed observations are across time.
# The TSS is based on normalised gaps between observation years (lower = more clumped).
#
# Arguments:
# - data: Data frame with observation years and taxon names
# - taxon_column: Name of the column with taxon names
# - year_column: Name of the column with observation years
# - override_max_gap: Optional override for normalization (max gap between years)
# - by_species: If TRUE, normalization is computed per species
#
# Returns:
# - A data frame with one row per taxon and its TSS score

score_temporal_spread <- function(data,
                                  taxon_column = "taxon_name",
                                  year_column = "year",
                                  override_max_gap = NULL,
                                  by_species = FALSE) {
  
  # Validate input columns
  required_cols <- c(taxon_column, year_column)
  if (!all(required_cols %in% names(data))) {
    stop("Missing required columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Clean and prepare year column
  cleaned_data <- data %>%
    dplyr::mutate(observed_year = as.numeric(.data[[year_column]])) %>%
    dplyr::filter(!is.na(observed_year))
  
  if (nrow(cleaned_data) == 0) {
    stop("No valid year values found in the dataset.")
  }
  
  ### Helper function to compute TSS from a vector of years and a max gap ###
  compute_tss <- function(years_vector, max_gap) {
    years <- sort(unique(years_vector))
    if (length(years) < 2) return(NA_real_)
    gaps <- diff(years)
    normalised_gaps <- gaps / max_gap
    1 - mean(normalised_gaps)
  }
  
  #### Global max gap across all taxa ###
  global_max_gap <- cleaned_data %>%
    dplyr::mutate(taxon = .data[[taxon_column]]) %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarise(
      year_list = list(sort(unique(observed_year))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      max_gap = purrr::map_dbl(year_list, function(yrs) {
        if (length(yrs) < 2) NA_real_ else max(diff(yrs))
      })
    ) %>%
    dplyr::pull(max_gap) %>%
    suppressWarnings() %>%
    max(na.rm = TRUE)
  
  ### Validate override ###
  if (!is.null(override_max_gap)) {
    override_max_gap <- as.numeric(override_max_gap)
    if (override_max_gap > global_max_gap) {
      stop("Provided override_max_gap (", override_max_gap, ") exceeds global max gap (", global_max_gap, ").")
    }
  }
  
  ### Use override max gap for all if provided ###
  if (!is.null(override_max_gap)) {
    result <- cleaned_data %>%
      dplyr::mutate(taxon = .data[[taxon_column]]) %>%
      dplyr::group_by(taxon) %>%
      dplyr::summarise(
        TSS = compute_tss(observed_year, override_max_gap),
        .groups = "drop"
      )
    
    ### Normalise max gap per species ###
  } else if (by_species) {
    result <- cleaned_data %>%
      dplyr::mutate(taxon = .data[[taxon_column]]) %>%
      dplyr::group_by(taxon) %>%
      dplyr::summarise(
        TSS = {
          years <- sort(unique(observed_year))
          if (length(years) < 2) NA_real_ else {
            max_gap <- max(diff(years))
            compute_tss(years, max_gap)
          }
        },
        .groups = "drop"
      )
    
    ### Use global max gap for all ###
  } else {
    result <- cleaned_data %>%
      dplyr::mutate(taxon = .data[[taxon_column]]) %>%
      dplyr::group_by(taxon) %>%
      dplyr::summarise(
        TSS = compute_tss(observed_year, global_max_gap),
        .groups = "drop"
      )
  }
  
  # Rename taxon back to original name
  result <- result %>%
    dplyr::select(!!taxon_column := taxon, TSS)
  
  # Return TSS
  return(result)
  
}
