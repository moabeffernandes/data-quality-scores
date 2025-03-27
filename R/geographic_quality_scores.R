# Geographic Quality Scores

# Calculate Proportion of Records Within Known Native Range (WCVP) ----
#
# For each species, this function returns the proportion of occurrence records
# that fall within its known native range (based on WCVP level-3 regions).
#
# Arguments:
# - dataset: Data frame with coordinates and taxon info
# - longitude_col: Column name for longitude
# - latitude_col: Column name for latitude
# - taxon_column: Column name with taxon names (e.g., "taxon_name")
# - species_id_column: Column name for plant_name_id (matching WCVP)
# - unique_id: ID column for individual records (default: "gbifID")
# - return_by_record: If TRUE, returns match result for each record
#
# Returns:
# - A summary table with the proportion of points per taxon inside its native range,
#   or a record-level table if return_by_record = TRUE

score_range_match_proportion <- function(dataset,
                                         longitude_col = "decimalLongitude",
                                         latitude_col = "decimalLatitude",
                                         taxon_column = "taxon_name",
                                         species_id_column = "accepted_species_id",
                                         unique_id = "gbifID",
                                         return_by_record = FALSE) {
  
  # Filter native distribution records from WCVP
  native_distribution <- rWCVPdata::wcvp_distributions %>%
    dplyr::filter(plant_name_id %in% unique(dataset[[species_id_column]])) %>%
    dplyr::filter(introduced == 0, extinct == 0, location_doubtful == 0) %>%
    dplyr::select(plant_name_id, area_code_l3)
  
  # Prepare taxon ID and name lookup table
  taxon_lookup <- dataset %>%
    dplyr::select(
      species_id = tidyselect::all_of(species_id_column),
      taxon = tidyselect::all_of(taxon_column)
    ) %>%
    dplyr::distinct()
  
  # Join taxon names into native distribution data
  native_distribution <- native_distribution %>%
    dplyr::left_join(taxon_lookup, by = c("plant_name_id" = "species_id"))
  
  # Load Level-3 region map
  level3_map <- rWCVPdata::wgsrpd3 %>%
    dplyr::select(LEVEL3_COD)
  
  # Prepare occurrence data as spatial points
  spatial_input <- dataset %>%
    dplyr::filter(!is.na(.data[[longitude_col]]) & !is.na(.data[[latitude_col]]))
  
  points_sf <- sf::st_as_sf(spatial_input,
                            coords = c(longitude_col, latitude_col),
                            crs = 4326,
                            remove = FALSE)
  
  # Spatial join with Level-3 region polygons
  joined_points <- sf::st_join(points_sf, level3_map, left = TRUE)
  
  # Standardise column names for internal processing
  points_clean <- joined_points %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(taxon = .data[[taxon_column]]) %>%
    dplyr::select(
      !!unique_id,
      species_id = tidyselect::all_of(species_id_column),
      taxon,
      LEVEL3_COD
    )
  
  # Match occurrence points with native regions
  match_result <- points_clean %>%
    dplyr::left_join(
      native_distribution,
      by = c("taxon", "LEVEL3_COD" = "area_code_l3")
    ) %>%
    dplyr::mutate(in_native_range = !is.na(plant_name_id))
  
  # Return per-record if requested
  if (return_by_record) {
    return(
      match_result %>%
        dplyr::rename(
          !!taxon_column := taxon,
          !!species_id_column := species_id
        )
    )
  }
  
  # Otherwise, summarise per species
  summary_result <- match_result %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarise(
      proportion_points = mean(in_native_range, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(!!taxon_column := taxon)
  
  # Return summary
  return(summary_result)
  
}



# Calculate Distribution Gap Score ----
#
# This function measures what proportion of a species' expected native regions
# (from WCVP) are actually represented in your dataset by at least one record.
#
# Arguments:
# - dataset: Data frame with species, coordinates, and plant_name_id
# - longitude_col: Column name for longitude
# - latitude_col: Column name for latitude
# - taxon_column: Column name with taxon names
# - species_id_column: Column with WCVP plant_name_id
#
# Returns:
# - A tibble with one row per taxon, showing the proportion of native regions with records

score_range_gap <- function(dataset,
                            longitude_col = "decimalLongitude",
                            latitude_col = "decimalLatitude",
                            taxon_column = "taxon_name",
                            species_id_column = "accepted_species_id") {
  
  # Filter native distribution regions from WCVP
  native_distribution <- rWCVPdata::wcvp_distributions %>%
    dplyr::filter(plant_name_id %in% unique(dataset[[species_id_column]])) %>%
    dplyr::filter(introduced == 0, extinct == 0, location_doubtful == 0) %>%
    dplyr::select(plant_name_id, area_code_l3)
  
  # Create lookup of taxon names and join to distribution
  taxon_lookup <- dataset %>%
    dplyr::select(
      species_id = tidyselect::all_of(species_id_column),
      taxon = tidyselect::all_of(taxon_column)
    ) %>%
    dplyr::distinct()
  
  native_distribution <- native_distribution %>%
    dplyr::left_join(taxon_lookup, by = c("plant_name_id" = "species_id"))
  
  # Load WCVP Level 3 region map
  level3_map <- rWCVPdata::wgsrpd3 %>%
    dplyr::select(LEVEL3_COD)
  
  # Convert dataset to spatial points
  spatial_input <- dataset %>%
    dplyr::filter(!is.na(.data[[longitude_col]]) & !is.na(.data[[latitude_col]])) %>%
    dplyr::mutate(
      longitude = .data[[longitude_col]],
      latitude = .data[[latitude_col]]
    )
  
  points_sf <- sf::st_as_sf(spatial_input,
                            coords = c("longitude", "latitude"),
                            crs = 4326,
                            remove = FALSE)
  
  # Join points with Level-3 regions
  points_with_region <- sf::st_join(points_sf, level3_map, left = TRUE) %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(taxon = .data[[taxon_column]]) %>%
    dplyr::select(
      species_id = tidyselect::all_of(species_id_column),
      taxon,
      LEVEL3_COD
    )
  
  # Check which native regions have actual records
  region_matches <- native_distribution %>%
    dplyr::left_join(
      points_with_region,
      by = c("taxon", "area_code_l3" = "LEVEL3_COD")
    ) %>%
    dplyr::mutate(record_found = !is.na(species_id)) %>%
    dplyr::distinct(taxon, area_code_l3, record_found)
  
  # Summarise: proportion of native regions with data per species
  summary_result <- region_matches %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarise(
      proportion_gap = mean(record_found),
      .groups = "drop"
    ) %>%
    dplyr::rename(!!taxon_column := taxon)
  
  # Return proportion of native regions with data per species
  return(summary_result)
  
}



# Score Out-of-Range Records ----
#
# This function assigns lower scores to records falling outside the native range
# of each species, based on their distance from the native centroid using a logistic
# decay function
#
# Arguments:
# - dataset: Data frame with taxon, species_id, and coordinates
# - longitude_col: Name of longitude column
# - latitude_col: Name of latitude column
# - taxon_column: Column name for taxon names
# - species_id_column: Column name for WCVP plant_name_id
# - return_by_record: If TRUE, returns per-record scores; otherwise per-taxon averages
# - include_in_range: If FALSE, exclude in-range records from the output
# - k: Logistic function steepness
# - inflection_prop: Distance (as proportion of max) where logistic curve drops
#
# Returns:
# - A tibble with either one row per taxon or one row per record

score_out_of_range_points <- function(dataset,
                                      longitude_col = "decimalLongitude",
                                      latitude_col = "decimalLatitude",
                                      taxon_column = "taxon_name",
                                      species_id_column = "accepted_species_id",
                                      return_by_record = FALSE,
                                      include_in_range = TRUE,
                                      k = 10,
                                      inflection_prop = 0.3) {
  
  # Filter valid coordinates
  dataset_clean <- dataset %>%
    dplyr::filter(!is.na(.data[[longitude_col]]) & !is.na(.data[[latitude_col]]))
  
  # Convert to sf object
  points_sf <- sf::st_as_sf(dataset_clean,
                            coords = c(longitude_col, latitude_col),
                            crs = 4326,
                            remove = FALSE)
  
  # Load native distribution
  native_distribution <- rWCVPdata::wcvp_distributions %>%
    dplyr::filter(plant_name_id %in% unique(dataset[[species_id_column]])) %>%
    dplyr::filter(introduced == 0, extinct == 0, location_doubtful == 0) %>%
    dplyr::select(plant_name_id, area_code_l3) %>%
    dplyr::left_join(
      dataset %>%
        dplyr::select(
          species_id = tidyselect::all_of(species_id_column),
          taxon = tidyselect::all_of(taxon_column)
        ) %>%
        dplyr::distinct(),
      by = c("plant_name_id" = "species_id")
    )
  
  # Spatial join with Level 3 regions
  level3_map <- rWCVPdata::wgsrpd3 %>%
    dplyr::select(LEVEL3_COD)
  
  points_with_region <- sf::st_join(points_sf, level3_map, left = TRUE) %>%
    dplyr::mutate(taxon = .data[[taxon_column]]) %>%
    sf::st_transform(3857)  # Project to meters for distance calculations
  
  # Score per species
  species_list <- unique(points_with_region$taxon)
  record_results <- list()
  summary_results <- list()
  
  for (sp in species_list) {
    sp_points <- points_with_region %>% dplyr::filter(taxon == sp)
    
    native_codes <- native_distribution %>%
      dplyr::filter(taxon == sp) %>%
      dplyr::pull(area_code_l3) %>%
      unique()
    
    if (length(native_codes) == 0) next
    
    sp_points <- sp_points %>%
      dplyr::mutate(in_range = LEVEL3_COD %in% native_codes)
    
    in_range_pts <- sp_points %>% dplyr::filter(in_range)
    out_range_pts <- sp_points %>% dplyr::filter(!in_range)
    
    if (nrow(out_range_pts) == 0 || nrow(in_range_pts) < 2) {
      scored_points <- sp_points %>%
        dplyr::mutate(geo_score = ifelse(in_range, 1, NA)) %>%
        sf::st_drop_geometry()
    } else {
      in_coords <- sf::st_coordinates(in_range_pts)
      centroid <- colMeans(in_coords)
      out_coords <- sf::st_coordinates(out_range_pts)
      
      distances <- sqrt((out_coords[, 1] - centroid[1])^2 + (out_coords[, 2] - centroid[2])^2)
      max_dist <- max(distances)
      inflection_point <- inflection_prop * max_dist
      
      scores <- 1 / (1 + exp(k * (distances - inflection_point) / inflection_point))
      
      out_range_pts$geo_score <- scores
      in_range_pts$geo_score <- 1
      
      scored_points <- dplyr::bind_rows(in_range_pts, out_range_pts) %>%
        sf::st_drop_geometry()
    }
    
    # Filter if needed
    filtered_points <- if (include_in_range) scored_points else scored_points %>% dplyr::filter(!in_range)
    
    if (return_by_record) {
      record_results[[sp]] <- filtered_points
    } else {
      summary <- filtered_points %>%
        dplyr::group_by(taxon) %>%
        dplyr::summarise(
          mean_geo_score = mean(geo_score, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::rename(!!taxon_column := taxon)
      
      summary_results[[sp]] <- summary
    }
  }
  
  # Combine and return results
  if (return_by_record) {
    return(dplyr::bind_rows(record_results))
  } else {
    return(dplyr::bind_rows(summary_results))
  }
  
}
