# Taxonomic Quality Scores


# Calculate Basis of Record Quality Score (BRQS) ----
#
# This function assigns a score to each species based on the type of record
# (e.g., preserved specimen, observation). More reliable record types score higher.
#
# Arguments:
# - data: Data frame containing taxon and record type columns
# - taxon_column: Name of the column with taxon names
# - record_column: Name of the column with basis of record types
# - score_mapping: Optional named vector assigning scores to record types
#
# Returns:
# - A data frame with one row per taxon and its BRQS score

score_basis_of_record <- function(data,
                                  taxon_column = "taxon_name",
                                  record_column = "basisOfRecord",
                                  score_mapping = NULL) {
  
  # Default score mapping if none is provided
  default_scores <- c(
    "PRESERVED_SPECIMEN" = 1.0,
    "MATERIAL_SAMPLE" = 0.9,
    "MATERIAL_CITATION" = 0.8,
    "MACHINE_OBSERVATION" = 0.7,
    "LIVING_SPECIMEN" = 0.6,
    "HUMAN_OBSERVATION" = 0.5,
    "OBSERVATION" = 0.4,
    "OCCURRENCE" = 0.2
  )
  
  score_lookup <- if (is.null(score_mapping)) default_scores else score_mapping
  
  # Check that required columns exist
  required_cols <- c(taxon_column, record_column)
  if (!all(required_cols %in% names(data))) {
    stop("Missing required columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Select and rename relevant columns for internal use
  clean_data <- data %>%
    dplyr::mutate(
      taxon = .data[[taxon_column]],
      record_type = .data[[record_column]]
    ) %>%
    dplyr::select(taxon, record_type) %>%
    dplyr::mutate(score = score_lookup[record_type])
  
  # Compute BRQS per taxon
  result <- clean_data %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarise(
      BRQS = sum(score, na.rm = TRUE) / sum(!is.na(score)),
      .groups = "drop"
    ) %>%
    dplyr::rename(!!taxon_column := taxon)
  
  # Return BRQS
  return(result)
  
}




# Calculate Name Confidence Quality Score (NQS) ----
#
# This function scores taxonomic name confidence based on synonym types.
# Accepted names get the highest score; homotypic and heterotypic synonyms score lower.
#
# Arguments:
# - data: A data frame with taxon names and synonym types
# - taxon_column: Column name with taxon names
# - synonym_column: Column name with synonym type (e.g., "Accepted", "Synonym_homotypic")
# - score_mapping: Optional named vector of synonym type scores
#
# Returns:
# - A data frame with one row per taxon and its average NQS score

score_name_confidence <- function(data,
                                  taxon_column = "taxon_name",
                                  synonym_column = "synonym_type",
                                  score_mapping = NULL) {
  
  # Default score mapping for synonym types
  default_scores <- c(
    "Accepted" = 1.0,
    "Synonym_homotypic" = 0.7,
    "Synonym_heterotypic" = 0.4
  )
  
  scores <- if (is.null(score_mapping)) default_scores else score_mapping
  
  # Check required columns exist
  required_cols <- c(taxon_column, synonym_column)
  if (!all(required_cols %in% names(data))) {
    stop("Missing required columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Select and rename relevant columns for internal use
  clean_data <- data %>%
    dplyr::select(
      taxon = tidyselect::all_of(taxon_column),
      synonym_type = tidyselect::all_of(synonym_column)
    ) %>%
    dplyr::mutate(
      score = scores[synonym_type]  # Assign scores based on type; returns NA if not found
    )
  
  # Compute mean score per taxon
  result <- clean_data %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarise(
      NCQS = mean(score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(!!taxon_column := taxon)
  
  # Return NCQS
  return(result)
  
}