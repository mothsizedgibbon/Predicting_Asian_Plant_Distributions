## MINI README ##
# This code models an SDM for a single species, as part of TEAMAP
# It was written by Matthew Shaun Grainger
# SETUP:
  # Make a single working directory for everything (which I have called "XTBGProjectCode)
  # Install java
  # Download maxent.jar and place it in the working directory
  # Install all of the packages
  # Place the LGM raster (cclgmbi3.tif) in the working directory
  # Place the variables data (TEA_climsoil_all_predictors_30s_20102023.tif) in the working directory
  # Place the occurrence data (TEA_30s_woody_successful_spp_occurrences_clean_noIAS_allsources_14062023.csv) in the working directory
  # Set the species_name to the name of the species that you want to model from the occurrence data file
# OUTPUTS:
  # _restricted_ENMeval_object.rds: Full ENMevaluate object for potential reanalyses 
  # _restricted_habitat_suitability_map.tif: Habitat suitability SpatRaster within accessible area as a .tif file for potential reanalyses
  # _restricted_habitat_suitability_map_vis.png: Image of habitat suitability map within accessible area overlayed upon a map of the land
  # _restricted_performance_metrics.csv: Performance metrics for the best model
  # _restricted_rmm_metadata.csv: Comprehensive range model metadata
  # _restricted_variable_importance_best_model.csv: All of the variables that were used in the model and their importance
# NOTES:
  # All outputs can be found within ParallelisedSDMs_Output which is a subdirectory that the code makes within the working directory
  # All output file names are prefixed by the species name with an underscore replacing the space
  # This code generates a habitat suitability map that is restricted to only the accessible area

## Preparing workspace
# Setting working directory
#setwd("E:/Matthew/XTBGProjectCode") # This directory contains this R script, the occurrence data file and the environmental data file
setwd("C:/Users/User/Documents/XTBGProjectCode") # This directory contains this R script, the occurrence data file and the environmental data file
# Loading packages
library(rJava) # Needed for Maxent
library(dismo) # This package is used for the MESS prediction
library(terra) # Used for handling raster files
library(sf) # This package is used to make a spatial object
library(ENMeval) # This package is used for partitioning data and making SDMs
library(corrplot) # For visualising correlations between explanatory variables
library(caret) # Used for the findCorrelation() function
library(dplyr) # Used for cleaning occurrence data
library(usdm) # Used for calculating VIF
library(CoordinateCleaner) # For cleaning coordinates (removing anomalous points)
library(rangeModelMetadata) # For getting comprehensive metadata for the SDMs
library(ggplot2) # For suitability map PNG
library(ggspatial)  # For adding scale bars and north arrows to suitability map PNG
library(viridis) # Colour-blind friendly colour pallete for suitability map PNG
library(grid) # For using ggplot for the suitability map PNG
library(ggnewscale) # For using ggplot for the suitability map PNG
library(raster) # Needed because dismo uses raster and not terra, and because need to use dismo for ENMevaluate() with maxent.jar
library(svglite) # For exporting visualisation as an SVG vector file
library(stringr) # For extracting list of variables from the best_model
library(tibble) # Involved in the process of extracting variable importance
# Setting a random seed to make it easy to reproduce this analysis
set.seed(88888888)
# Loading the Maxent algorithm (maxent.jar)
options(dismo.maxent = "E:/Matthew/XTBGProjectCode/maxent.jar") # Tells dismo where to find it
if (!file.exists("maxent.jar")) stop("maxent.jar not found in the working directory.") # Stops the script if cannot find maxent.jar
# Making a directory for the output files, if there is not one already
if (!dir.exists("ParallelisedSDMs_Output")) { # `If the directory does not already exist
  dir.create("ParallelisedSDMs_Output") # Makes it
}

## Make sure not parallel
terra::terraOptions(threads = 1)
if ("package:rJava" %in% search()){
  options(java.parameters = "-Xmx60G")
}
parallel <- FALSE
foreach::registerDoSEQ()


## Importing environmental data and preparing it
variables_data <- terra::rast("TEA_climsoil_all_predictors_30s_20102023.tif") # Explanatory variables data
names(variables_data) <- c("bio01", "bio02", "bio03", "bio04", "bio05","bio06","bio07","bio08","bio09","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19", "aridity", "bedrock_depth", "bulk_density", "cec_30cm", "clay_30cm", "orgC_30cm", "soilpH_30cm", "silt_30cm", "sand_30cm")
temp_vars <- c("bio01", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11") # All the temperature-related variables
precip_vars <- c("bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19") # All the precipitation-related variables
soil_vars <- c("aridity", "bedrock_depth", "bulk_density", "cec_30cm", "clay_30cm", "orgC_30cm", "soilpH_30cm", "silt_30cm", "sand_30cm") # All the soil-related variables
var_groups <- list(temperature = temp_vars, precipitation = precip_vars, soil = soil_vars) # List of the different categories of variables
# Importing map of LGM to identify which land masses were connected during the LGM
lgm_raster <- rast("cclgmbi3.tif") # Importing the map of the last glacial maximum, to find out which areas of land were connected
lgm_raster <- classify(lgm_raster, matrix(c(NA, NA, 0,
                                            0, Inf, 1), 
                                            ncol = 3, byrow = TRUE)) # Makes it so that water has a value of 0 and land has a value of one (binary)
# Making a map of modern land vs sea
land_sea_map <- variables_data[[1]] # Extracting first bioclim variable
land_sea_map <- classify(land_sea_map, matrix(c(NA, NA, 0,
                                              0, Inf, 1), 
                                              ncol = 3, byrow = TRUE)) # Makes it so that water has a value of 0 and land has a value of one (binary)


## Importing occurrence data
all_occurrence_data <- read.csv("TEA_30s_woody_successful_spp_occurrences_clean_noIAS_allsources_14062023.csv", stringsAsFactors = FALSE)  # Load all occurrence data once


## Choosing species for this single SDM
species_name <- "Pistacia weinmanniifolia" # Choosing the species
safe_name <- tolower(gsub("[^[:alnum:]]", "_", species_name))  # Replace non-alphanumeric characters with underscores to create safe filename (remove spaces and special characters)  


## Filter and clean occurrence data for this species
occurrence_data <- dplyr::filter(all_occurrence_data, Species == species_name) # Extracting the occurrence points for the chosen species
occurrence_data <- occurrence_data[!duplicated(occurrence_data[, c("longmid", "latmid")]), ] # Removing duplicated coordinates
cc_result  <- clean_coordinates(occurrence_data, lon = "longmid", lat = "latmid", # Flagging points that are spatial outliers
                                     species = "Species", tests = "outliers", 
                                     outliers_method = "quantile", outliers_mtp = 5, 
                                     outliers_size = 7,
                                     verbose = TRUE, report = FALSE)
occurrence_data <- as.data.frame(occurrence_data[cc_result$.summary, ]) # Removing spatial outliers and then converting back to a normal data frame
occurrence_data_coords <- occurrence_data[, c("longmid", "latmid")] # Keep only coordinate columns
occurrence_points_vect <- terra::vect(occurrence_data_coords, geom = c("longmid", "latmid"), crs = crs(variables_data)) # Spatial object of occurrence points
reference_variables <- terra::extract(variables_data, occurrence_points_vect, ID = FALSE)  # Extract environmental values at occurrence points
occs.na <- which(rowSums(is.na(reference_variables)) > 0)  # Find points with missing environmental data
if(length(occs.na) > 0) {  # If there are occurrence points with missing environmental data
  occurrence_data_coords <- occurrence_data_coords[-occs.na, ] # Remove from occurrence data coords, as this is later used
  occurrence_points_vect <- occurrence_points_vect[-occs.na, ] # Also remove from spatial points object of occurrence points, in case change code to use it later
  occurrence_data <- occurrence_data[-occs.na, ] # Also remove from the main occurrence_data df which has the species name, in case change code to use it later
}
if(nrow(occurrence_data_coords) < 5) { # Skip species that have less than 5 occurrence points as this is too few
  stop(paste0(species_name, " has less than 5 occurrence points, so skipping SDM.")) # Skip and tell the user that the species has less than 5 occurrence points
  }


## Getting the accessible area for the species by connecting 500km buffers around occurrence points
occs.sf <- sf::st_as_sf(occurrence_data_coords, coords = c("longmid","latmid"), crs = terra::crs(variables_data))  # Convert occurrence data to spatial object
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"  # Define equal-area projection CRS because WGS84 (lat/lon) isn't good for distance calculations in metres
occs.sf_proj <- sf::st_transform(occs.sf, crs = eckertIV)  # Transform the occurrence points to equal-area projection to enable accurate buffering in metres
buffer_radius <- 500000 # Radius of the buffer will be 500km (500,000m) 
accessible_area <- occs.sf_proj |> # Make species-specific accessible area
  sf::st_buffer(dist = buffer_radius) |>  # Buffer each point
  sf::st_union() |> # Merge all buffers into one polygon
  sf::st_sf() |> # Ensure that this is a valid sf object
  sf::st_transform(crs = terra::crs(variables_data))  # Transform back to match rasters
accessible_area_vect <- terra::vect(accessible_area)  # Convert to SpatVsector terra object
accessible_mask <- terra::rasterize(accessible_area_vect, variables_data, field = 1) # Create binary accessible area mask


## Restricting the accessible area to only include the areas that were part of the same landmass as areas with occurrence points, during the last glacial maximum
connected_lands <- terra::patches(lgm_raster == 1, directions = 8) # Labels clusters of land (landmasses)
occurrence_points <- terra::vect(occurrence_data_coords, geom = c("longmid", "latmid"), crs = terra::crs(variables_data)) # Makes occurrence points into spatial object
occurrence_raster <- terra::rasterize(occurrence_points, connected_lands, field = 1) # Turning occurrence points into a raster
connected_with_occurrences <- unique(terra::values(connected_lands)[!is.na(terra::values(occurrence_raster))]) # Landmasses that contain at least one occurrence point
landmass_mask <- connected_lands %in% connected_with_occurrences # Creates raster where cells belonging to the relevant land masses are marked as 1 and the others are NA
landmass_mask_aligned <- terra::project(landmass_mask, accessible_mask, method = "near") # 
final_accessible_mask <- terra::mask(accessible_mask, landmass_mask_aligned) # Mask combining 500KM buffer-based accessible area and LGM land connectivity restraints
if (all(is.na(terra::values(final_accessible_mask)))) { # If the mask has only NA values
  stop(paste0(species_name, " has an accessible area with only NA values.")) # Stop and tell the user that the accessible area only has NA values
}
final_accessible_vect <- terra::as.polygons(final_accessible_mask) # Converting the mask to a vector
accessible_area_km2 <- terra::expanse(final_accessible_vect, unit = "km") # Accessible area mask size
if (accessible_area_km2 < 10) { # If the size of the accessible area mask is too small
  stop(paste0(species_name, " has an accessible area too small (smaller than 10km2).")) # Stop and tell user that accessible area too small
}


## Removing variables that are multicollinear within the accessible area for this species (whilst ensuring that there is at least 1 soil variable, 1 temperature variable and 1 precipitation variable)
env_data_restricted <- terra::mask(variables_data, final_accessible_mask) # Restricted variable data for use in calculating VIF
sample_full <- terra::spatSample(variables_data, size = 5000, as.raster = TRUE) # Sampling from the full raster (to use for finding proportion of cells in accessible area)
sample_restricted <- terra::spatSample(env_data_restricted, size = 5000, as.raster = TRUE) # Sampling from the restricted raster (to use for finding proportion of cells in accessible area)
non_na_full <- terra::global(!is.na(sample_full), "sum", na.rm = TRUE)[1,1] # Count non-NA cells in whole geographic region (to use for finding proportion of cells in accessible area)
non_na_restricted <- terra::global(!is.na(sample_restricted), "sum", na.rm = TRUE)[1,1] # Count non-NA cells in restricted region (to use for finding proportion of cells in accessible area)
proportion <- non_na_restricted / non_na_full # Proportion of cells in restricted compared to full, based upon the samples (as the only difference between these two rasters is that there are more NAs in the restricted one and therefore less non-NA cells, so just need to look at the proportion of non-NA cells)
reference_sample_size <- 36200 # Reference sample size is 36,200 because this gives roughly 10,000 points if applied to the whole map
proportional_sample_size <- round(proportion * reference_sample_size) # This is the proportional sample size to use when calculating VIF
restricted_environmental_data_sample <- as.data.frame(terra::spatSample(env_data_restricted, size = proportional_sample_size, method = "regular", na.rm = TRUE)) # Sampling the correct amount of points to calculate VIF
if (nrow(restricted_environmental_data_sample) < 10) { # If the sample size is less than 10
  stop(paste0(species_name, " has an insufficient sample size for VIF calculation.")) # Stop and tell the user that the sample size is too small
}
final_vars <- list() # Initiating list of final variables
for (group in names(var_groups)) { # For every group of variables
  vars <- var_groups[[group]] # The variables in the current group
  group_data <- restricted_environmental_data_sample[, vars, drop = FALSE] # Keeping only the columns with the current variables
  group_data <- na.omit(group_data) # Get rid of rows with NAs
  if (nrow(group_data) < 2) { # If the group data has less than 2 rows after removing NAs
    print(paste0(species_name, " does not have enough data for VIF in group", group, ", so using first variable.")) # Tell the user that there isn't enough data calculate VIF, so will use the first variable
    final_vars[[group]] <- vars[1] # Use the first variable if cannot calculate VIF (COULD BE GOOD TO DECIDE WHICH VARIABLE TO DEFAULT TO BASED UPON ECOLOGICAL RELEVANCE)
    next
  }
  vif_result <- usdm::vifstep(group_data, th = 5) # Calculate VIF in a stepwise process
  kept_vars <- vif_result@variables # Only keeping variables below the VIF threshold
  if (length(kept_vars) == 0) { # If no variables are below the VIF threshold
    full_vif_result <- usdm::vif(group_data) # Get the full set of VIF results
    lowest_vif_var <- names(full_vif_result)[which.min(full_vif_result)] # Get the name of the variable in the group with the lowest VIF
    kept_vars <- lowest_vif_var # Keeping only this variable
  }
  final_vars[[group]] <- kept_vars # The final variables for this group
}
all_kept_vars <- unlist(final_vars) # Removing the categories so that the variables are in one vector
env_data_restricted <- env_data_restricted[[all_kept_vars]]  # For sampling background points
variables_data_for_restricted <- variables_data[[all_kept_vars]] # For full-extent prediction


## Creating the background points within only the accessible area, for the restricted suitability map
bg_density_per_km2 <- 1 / 2  # Defining the background point density to be 1 point per 2 km²
bg_points_n <- round(accessible_area_km2 * bg_density_per_km2) # Calculating number of bg points needed for this density
bg_points_n <- min(bg_points_n, 5000)  # Maximum number of background points is 5000 for performance (and there is no minimum number)
bg <- terra::spatSample(env_data_restricted, size = bg_points_n, values = FALSE, xy = TRUE, na.rm = TRUE) |> # Sample appropriate number of background points
  as.data.frame()
colnames(bg)[1:2] <- c("longmid", "latmid") # Changing to the correct column names
bg <- bg[!paste0(bg$longmid, "_", bg$latmid) %in% paste0(occurrence_data_coords$longmid, "_", occurrence_data_coords$latmid), ] # Get rid of background points with the same coordinates as occurrence points
bg_points_used <- nrow(bg) # Storing number of background points used so that can be added to performance summary
if (nrow(bg) < 10) stop("Too few background points for reliable modeling.") # Stopping and telling the user that there are not enough background points


## Run ENMeval to make the restricted suitability map
occurrence_data_coords_xy <- occurrence_data_coords # Formatting the occurrence data coordinates for ENMevaluate
colnames(occurrence_data_coords_xy) <- c("x", "y") # By changing the column names to x and y
colnames(bg)[1:2] <- c("x", "y") # Change bg columns to match occurrence data
tune.args <- list(fc = c("L", "LQ", "LQH", "LQHP", "LQHPT", "H"), rm = c(0.5, 1, 2))  # Define model settings to test
# Running ENMevaluate
e.mx <- tryCatch(ENMevaluate(occs = occurrence_data_coords_xy, envs = env_data_restricted, bg = bg,  # Run model evaluation upon the accessible area
                    algorithm = 'maxent.jar', partitions = 'block',  # Use MaxEnt with block partitioning
                    tune.args = tune.args),  # Test the specified settings
error = function(e) stop(paste0("ENMevaluate failed for ", species_name, ": ", e$message))) # Prints the error if there is one
# For Pistacia with these settings it took 213 mins


## Save full ENMeval object to a file
saveRDS(e.mx, file = file.path("ParallelisedSDMs_Output", paste0(safe_name, "_restricted_ENMeval_object.rds"))) # Saving the full ENMeval object in case want to re-analyse stuff later


## Model selection for the restricted suitability map
model_evaluation_results <- ENMeval::eval.results(e.mx)  # Extract model evaluation results
model_evaluation_results <- model_evaluation_results[order(model_evaluation_results$delta.AICc), ]  # Sort by AICc (best first)
selected_model <- NULL  # Initialize variable to store selected model
for(i in 1:nrow(model_evaluation_results)) { # For every model (each row is a model)
  AICc <- model_evaluation_results$delta.AICc[i]
  AUC <- model_evaluation_results$auc.val.avg[i]
  orMTP <- model_evaluation_results$or.mtp.avg[i]
  if (is.na(AICc) || is.na(AUC)){
    next
  }
  
  if(AUC >= 0.7 && # If average test AUC is greater than 0.7
     (is.na(orMTP) || orMTP <= 0.4)) { # If average test omission rate is less than 0.4
    selected_model <- model_evaluation_results[i, ] # Pick this model
    break # If find best model, stop the loop
  }
}
if(is.null(selected_model)) { # If none of the models were selected based upon criteria above
  for(i in 1:nrow(model_evaluation_results)) { # For every model
    AICc <- model_evaluation_results$delta.AICc[i]
    AUC <- model_evaluation_results$auc.val.avg[i]
    if (is.na(AICc) || is.na(AUC)){
      next
    }
    if (AUC >= 0.7){
      selected_model <- model_evaluation_results[i, ]
      break
    }
  }
}
if(is.null(selected_model)) { # If still no model selected
  selected_model <- model_evaluation_results[1, ] # Just pick model with best AICc
}
best_index <- which(model_evaluation_results$fc == selected_model$fc & # Extracting the index of the model using its fc and rm
                      model_evaluation_results$rm == selected_model$rm)[1]
best_model <- e.mx@models[[best_index]] # Extracting the model using this index
# Save best_model
saveRDS(best_model, file = file.path("ParallelisedSDMs_Output", paste0(safe_name, "_restricted_best_model.rds"))) # Saving the best model in case want to re-analyse stuff later


## Generate the restricted habitat suitability map and save it on its own to a .tif file
habitat_suitability_full <- ENMeval::maxnet.predictRaster( # Predict on FULL geographic extent
  best_model, 
  variables_data_for_restricted,  # Full-extent variables (already VIF-filtered)
  type = "cloglog", # Continuous probability output
  progress = "text" # Shows a progress bar
)
terra::crs(habitat_suitability_full) <- terra::crs(variables_data_for_restricted) # Assigning correct crs
if (all(is.na(values(habitat_suitability_full))) || terra::global(habitat_suitability_full, sd, na.rm = TRUE)[[1]] == 0) { # Check whether the suitability map worked or whether it did not work (all NAs)
  print("Predicted map is empty or constant.")
}
final_accessible_mask_aligned <- terra::resample(final_accessible_mask, habitat_suitability_full, method = "near")  # Resampling the fina_acceessible_mask to match the habitatu_suitability_full
habitat_suitability_final <- terra::mask( # Apply final accessible area mask
  habitat_suitability_full, 
  final_accessible_mask_aligned,  # Combined buffer + LGM mask
  inverse = FALSE  # Keep areas inside mask, set outside to NA
)
terra::writeRaster(habitat_suitability_final, filename = file.path("ParallelisedSDMs_Output", paste0(safe_name, "_restricted_habitat_suitability_map.tif")),  # Save habitat suitability map as .tif
                    overwrite = TRUE)


## Saving restricted habitat suitability map overlayed upon land-sea map as a PNG and as a SVG
hab_suit_df <- as.data.frame(habitat_suitability_final, xy = TRUE) # Converts habitat suitability map SpatRaster into data frame for ggplot
names(hab_suit_df)[3] <- "Suitability" # The third column contains raster values; renaming it to "Suitability" for clarity in the plot
land_poly_df <- as.data.frame(land_sea_map, xy = TRUE) # Converts land-sea map to a data frame for plotting in ggplot
names(land_poly_df)[3] <- "LandSea"
land_poly_df$color <- ifelse(land_poly_df$LandSea == 1, "#0D0887FF", "white") # Land is same colour as a habitat suitability of 0, sea is white
# Save visualisation as a SVG
svg( # Open PNG device for output file
  file.path("ParallelisedSDMs_Output", paste0(safe_name, "_restricted_habitat_suitability_map_vis.svg")),
  width = 8, height = 6
)
ggplot() + # New ggplot object
  geom_raster(data = land_poly_df, aes(x = x, y = y, fill = color)) + # Add land and sea mask
  scale_fill_identity() +  # Ensures grey/white are used directly, no legend for land/sea
  ggnewscale::new_scale_fill() +  # Allow a second fill scale
  geom_raster(data = hab_suit_df, aes(x = x, y = y, fill = Suitability), alpha = 0.9) + # Overlay habitat suitability
  scale_fill_viridis( # Uses colour-blind colour palette
    name = "Habitat suitability", # Legend title
    option = "C",
    limits = c(0, 1), # Colour scale consistent across species (for when this is automated)
    na.value = "transparent" # Hides cells with NA values
  ) +
  coord_equal() + # 1:1 aspect ratio to preserve spatial correctness
  theme_minimal(base_size = 14) + # Clean background
  theme( # fine-tunes fonts, removes gridlines, places legend, and styles the title
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
  ) +
  labs(title = bquote(bold(atop("Restricted habitat suitability map for ", italic(.(species_name)))))) # Adds title
dev.off() # Close PNG device to output file
# Save alternative visualisation as a SVG
land_poly_df$color <- ifelse(land_poly_df$LandSea == 1, "grey80", "white") # Land is grey, sea is white
svg( # Open SVG device for output file
  file.path("ParallelisedSDMs_Output", paste0(safe_name, "_restricted_habitat_suitability_map_altvis.svg")),
  width = 8, height = 6
)
ggplot() + # New ggplot object
  geom_raster(data = land_poly_df, aes(x = x, y = y, fill = color)) + # Add land and sea mask
  scale_fill_identity() +  # Ensures grey/white are used directly, no legend for land/sea
  ggnewscale::new_scale_fill() +  # Allow a second fill scale
  geom_raster(data = hab_suit_df, aes(x = x, y = y, fill = Suitability), alpha = 0.9) + # Overlay habitat suitability
  scale_fill_viridis( # Uses colour-blind colour palette
    name = "Habitat suitability", # Legend title
    option = "C",
    limits = c(0, 1), # Colour scale consistent across species (for when this is automated)
    na.value = "transparent" # Hides cells with NA values
  ) +
  coord_equal() + # 1:1 aspect ratio to preserve spatial correctness
  theme_minimal(base_size = 14) + # Clean background
  theme( # fine-tunes fonts, removes gridlines, places legend, and styles the title
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
  ) +
  labs(title = bquote(bold(atop("Restricted habitat suitability map for ", italic(.(species_name)))))) # Adds title
dev.off() # Close PNG device to output file

  
## Calculate TSS and Boyce Index for the best model
colnames(occurrence_data_coords)[colnames(occurrence_data_coords) == "longmid"] <- "x" # Making sure correct column names
colnames(occurrence_data_coords)[colnames(occurrence_data_coords) == "latmid"] <- "y"
occurrence_points <- terra::vect(occurrence_data_coords_xy, geom = c("x", "y"), crs = terra::crs(variables_data_for_restricted)) # aligning occurrence points to match raster
occ_data_for_pred <- terra::extract(variables_data_for_restricted, occurrence_points, ID = FALSE) # Variables data at occurrence points
occ_data_for_pred <- as.data.frame(occ_data_for_pred) # Converting to data frame to use for predicting
occ_prediction_df <- cbind(occurrence_data_coords_xy, occ_data_for_pred) # Combine coordinates and environmental values into a df
occ_pred <- dismo::predict(
  object = best_model,
  x = occ_prediction_df,
  type = "cloglog"
)
occ_pred <- as.numeric(occ_pred) # Suitability predictions from occurrence points
colnames(bg)[1:2] <- c("x", "y") # Making sure correct column names 
bg_points <- terra::vect(bg, geom = c("x", "y"), crs = terra::crs(variables_data_for_restricted)) # aligning bg points to match raster
bg_data_for_pred <- terra::extract(variables_data_for_restricted, bg_points, ID = FALSE) # Variables data at bg points
bg_data_for_pred <- as.data.frame(bg_data_for_pred) # Converting to data frame to use for predicting
bg_prediction_df <- cbind(bg, bg_data_for_pred) # Combine coordinates and environmental values into a df
bg_pred <- dismo::predict(
  object = best_model,
  x = bg_data_for_pred,
  type = "cloglog"
)
bg_pred <- as.numeric(bg_pred) # Suitability predictions from background points
if (length(occ_pred) == 0 || length(bg_pred) == 0) {
  best_tss <- NA
  best_boyce <- NA
  print(paste0("For ", species_name, " there are no predictions available for TSS or Boyce Index metrics."))
} else {
  # Calculate TSS for the best model
  thresholds <- seq(0, 1, by = 0.01) # Creates a sequence of 100 possible threshold values from 0 to 1 to maximise TSS
  tss_values <- numeric(length(thresholds)) # Creates empty numeric vector that stores TSS values for each threshold
  for(j in 1:length(thresholds)) { # For each threshold value
    thresh <- thresholds[j]
    true_pos <- sum(occ_pred >= thresh, na.rm = TRUE) # Calculates how many points correctly predicted positive
    false_neg <- sum(occ_pred < thresh, na.rm = TRUE) # Calculates how many points incorrectly predicted negative
    true_neg <- sum(bg_pred < thresh, na.rm = TRUE) # Correct negative
    false_pos <- sum(bg_pred >= thresh, na.rm = TRUE) # Incorrect positive
    sensitivity <- true_pos / (true_pos + false_neg) # Calculates sensitivity (true positives rate)
    specificity <- true_neg / (true_neg + false_pos) # Calculates specificity (true negatives rate)
    tss_values[j] <- sensitivity + specificity - 1 # Calculates TSS for this threshold: TSS = sensitivity + specificity - 1
  }
  best_tss <- max(tss_values, na.rm = TRUE) # Get the maximum TSS value
  # Calculate Boyce Index for the best model
  occ_pred <- occ_pred[!is.na(occ_pred)] # Remove NAs
  bg_pred <- bg_pred[!is.na(bg_pred)] # Remove NAs
  bins <- seq(0, 1, length.out = 11) # Make 10 equal-interval bins
  occ_bin <- findInterval(occ_pred, bins) # Assign each occurrence prediction to a bin
  bg_bin <- findInterval(bg_pred, bins) # Assign each background prediction to a bin
  occ_freq <- table(factor(occ_bin, levels = 1:10)) # Create a frequency table to count how many occurrence points fall into each bin
  bg_freq <- table(factor(bg_bin, levels = 1:10)) # Count how many background points fall into each bin
  predicted <- occ_freq / sum(occ_freq) # Calculate the proportion of occurrence points in each bin
  expected <- bg_freq / sum(bg_freq) # Proportion of background points in each bin
  pe_ratio <- predicted / expected # How much more frequent occurrences are than expected by chance
  pe_ratio[is.infinite(pe_ratio) | is.nan(pe_ratio)] <- NA # Replaces infinite/undefined values (from when expected is 0) with NA
  best_boyce <- cor(1:10, pe_ratio, use = "complete.obs") # Calculates Boyce Index which is correlation between bin number and P/E ratio
}


## Making performance summary and saving to file
performance_summary <- data.frame(  # Create data frame with key results for best model
species = safe_name, # Species
test_AUC = selected_model$auc.val.avg, # AUC
test_TSS = best_tss, # TSS
boyce_index = best_boyce, # Boyce index
delta_AICc = selected_model$delta.AICc, # AICc
omission_rate_MTP  = selected_model$or.mtp.avg, # orMTP
omission_rate_10p  = selected_model$or.10p.avg,  # Omission rate at 10% threshold
n_parameters = selected_model$ncoef,  # Number of parameters
n_occurrence_points = nrow(occurrence_data_coords), # Number of occurrence points
n_background_points = bg_points_used, # Number of background points used
fc = selected_model$fc, # FC
rm = selected_model$rm, # RM
stringsAsFactors = FALSE) # Prevent factor conversion for strings
write.csv(performance_summary, file.path("ParallelisedSDMs_Output", paste0(safe_name, "_restricted_performance_metrics.csv")), row.names = FALSE) # Write metrics to a CSV


## Save full metadata to a file for restricted suitability map
rmm <- eval.rmm(e.mx) # Generate a rangeModelMetadata object based on the information stored in the ENMevaluate object
rmm$model$selectionRules <- "best delta AICc with test AUC ≥ 0.7 and average test omission rate ≤ 0.4; fallback to best delta AICc with test AUC ≥ 0.7; final fallback to best delta AICc regardless of performance metrics" # Selection rules
rmm$model$selectionCriteria <- list( # Model selection criteria
  primary = "delta AICc",
  performanceThresholds = c("test AUC ≥ 0.7", "average test omission rate ≤ 0.4"),
  validationMethod = "jackknife"
)
rmm$model$finalModelSettings <- paste0("Features: ", selected_model$features, # FC and RM of best model
                                       " | RM: ", selected_model$rm)
rmm$model$performanceMetrics <- list(
  species = safe_name, # Species
  test_AUC = selected_model$auc.val.avg, # AUC
  test_TSS = best_tss, # TSS
  boyce_index = best_boyce, # Boyce index
  delta_AICc = selected_model$delta.AICc, # AICc
  omission_rate_MTP  = selected_model$or.mtp.avg, # orMTP
  omission_rate_10p  = selected_model$or.10p.avg,  # Omission rate at 10% threshold
  n_parameters = selected_model$ncoef,  # Number of parameters
  n_occurrence_points = nrow(occurrence_data_coords), # Number of occurrence points
  n_background_points = bg_points_used, # Number of background points used
  fc = selected_model$fc, # FC
  rm = selected_model$rm # RM
)
rmm$prediction$continuous$minVal <- terra::global(habitat_suitability_final, "min", na.rm = TRUE)[[1]] # Lowest suitability value
rmm$prediction$continuous$maxVal <- terra::global(habitat_suitability_final, "max", na.rm = TRUE)[[1]] # Highest suitability value
rmm$prediction$continuous$meanVal <- terra::global(habitat_suitability_final, "mean", na.rm = TRUE)[[1]] # Average suitability value
rmm$prediction$continuous$units <- "suitability (cloglog transformation)" # Unit of suitability score
rmm$data$environmental$variableSelection <- "VIF threshold 5 with a minimum of 1 temperature variable, 1 precipitation variable and 1 soil variable retained" # Variable selection
rangeModelMetadata::rmmToCSV(rmm, file = file.path("ParallelisedSDMs_Output", paste0(safe_name, "_restricted_rmm_metadata.csv"))) # Saving metadata to file


## Saving the coefficients of each variable to a file
lambdas <- best_model@lambdas # Has variable names
lambdas_list <- strsplit(lambdas, ",") # Split each string in the vector by commas
lambda_df <- as.data.frame(do.call(rbind, lambdas_list), stringsAsFactors = FALSE) # Data frame where first column has names of variables
colnames(lambda_df) <- c("Variable", "LambdaValue", "Min", "Max")
lambda_df <- lambda_df[order(as.numeric(lambda_df[, 2]), decreasing = TRUE), ]
lambda_df <- lambda_df %>%
  filter(!(Variable %in% c("numBackgroundPoints", "densityNormalizer", 
                           "entropy", "linearPredictorNormalizer")))
write.csv(lambda_df, file.path("ParallelisedSDMs_Output", paste0(safe_name, "_restricted_variable_coefficients.csv")), row.names = FALSE) # Write metrics to a CSV


## Saving importance of each variable to a file
best_model_results <- best_model@results %>%
  as.data.frame() %>%
  rownames_to_column(var = "model_information") # Creates a new column named 'model_setting'
var_importance <- best_model_results %>%
  filter(
    grepl("permutation.importance|contribution", model_information, ignore.case = TRUE)
  )
var_importance <- var_importance[order(as.numeric(var_importance[, 2]), decreasing = TRUE), ]
write.csv(var_importance, file.path("ParallelisedSDMs_Output", paste0(safe_name, "_restricted_variable_importance.csv")), row.names = FALSE) # Write metrics to a CSV

