# Visualising Quality Scores


# Libraries----
library(tidyverse)
library(sf)
library(rWCVP)
library(rWCVPdata)



# Load custom functions----
source("R/taxonomic_quality_scores.R")
source("R/temporal_quality_scores.R")
source("R/geographic_quality_scores.R")



# Data loading and Processing----
# Load data
dataset <- readRDS("data/test_dataset.rds") 

# Process data for calculating name quality score
dataset <- dataset %>%
  mutate(synonym_type = case_when(
    wcvp_status == "Synonym" & wcvp_homotypic == TRUE ~ "Synonym_homotypic",
    wcvp_status == "Synonym" & is.na(wcvp_homotypic) ~ "Synonym_heterotypic",
    TRUE ~ "Accepted"
  ))



# Define score functions----
# Taxonomic Quality Functions
taxonomic_funs <- list(
  brqs = function(df) score_basis_of_record(df, taxon_column = "taxon_name", record_column = "basisOfRecord"),
  ncqs  = function(df) score_name_confidence(df, taxon_column = "taxon_name", synonym_column = "synonym_type")
)

# Temporal Quality Functions
temporal_funs <- list(
  trqs  = function(df) score_temporal_recency(df, taxon_column = "taxon_name", year_column = "year"),
  ltqs = function(df) score_logistic_temporal_quality(df, taxon_column = "taxon_name", year_column = "year"),
  tss  = function(df) score_temporal_spread(df, taxon_column = "taxon_name", year_column = "year")
)

# Geographic Quality Functions
geographic_funs <- list(
  rmp = function(df) score_range_match_proportion(df, longitude_col = "decimalLongitude", latitude_col = "decimalLatitude", taxon_column = "taxon_name", species_id_column = "accepted_species_id", unique_id = "gbifID", return_by_record = FALSE),
  rgp = function(df) score_range_gap(df, longitude_col = "decimalLongitude", latitude_col = "decimalLatitude", taxon_column = "taxon_name", species_id_column = "accepted_species_id"),
  por = function(df) score_out_of_range_points(df, longitude_col = "decimalLongitude", latitude_col = "decimalLatitude", taxon_column = "taxon_name", species_id_column = "accepted_species_id", return_by_record = FALSE, include_in_range = TRUE, k = 10, inflection_prop = 0.3)
)



# Applying functions and summarising results----
# Apply functions to the dataset
taxonomic_results <- map(taxonomic_funs, ~ .x(dataset))
temporal_results <- map(temporal_funs, ~ .x(dataset))
geographic_results <- map(geographic_funs, ~ .x(dataset))

# List of data frames to join
score_tables <- list(
  reduce(taxonomic_results, full_join, by = "taxon_name"),
  reduce(temporal_results, full_join, by = "taxon_name"),
  reduce(geographic_results, full_join, by = "taxon_name"),
  distinct(dataset[, c("taxon_name", "area_category")])
)

# Join them all and pivot longer
results_complete <- score_tables %>%
  reduce(left_join, by = "taxon_name") %>%
  pivot_longer(cols = -c(taxon_name, area_category),
               names_to = "metric",
               values_to = "value")

# # Write to the disk
# write_csv(results_complete, "results/results_test.csv")

# Visualisation----
# Density plot with distribution of all metrics
ggplot(results_complete, aes(x = value, fill = area_category)) +
  geom_density(alpha = 0.6, color = NA) +
  facet_wrap(~ metric, scales = "free") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Score",
    y = "Density",
    fill = "Area Category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Save plot
ggsave(filename = "figures/density_quality_scores.jpg", plot = last_plot(),
       width = 8, height = 6, dpi = 300)
