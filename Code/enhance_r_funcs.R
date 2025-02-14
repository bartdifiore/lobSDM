####  Common Resources  ####
# Calculating average distance among raster layer centroids, used to interpolate NA values for points just outside covariate layer grid

# Function to calculate nearest distance, in space and time, between a point that has NA covariate value and one that doesn't
get_raster_centroid_dist<- function(r) {
  
  distance <- sqrt(res(terra::rast(r))[1]^2 + res(terra::rast(r))[2]^2)
  
  return(distance)
}

# Helper function to find and fill NA values
fill_na_spatiotemporal <- function(na_pts, non_na_pts, value_col, date_col) {
  if(FALSE){
    na_pts = na_points
    non_na_pts = non_na_points
    dist_mat = distances
    value_col = names(summ_df_list)[j]
    date_col = "Date"
  }
  # Ensure both datasets have the same CRS
  if (st_crs(na_pts) != st_crs(non_na_pts)) {
    na_pts <- st_transform(non_na_pts, st_crs(na_pts))
  }

  if(is.null(date_col)) {
    for (i in 1:nrow(na_pts)) {
      # Get the point and date for current NA value
      current_point <- na_pts[i,]
    
      # Calculate distance to all points in source dataset
      distances <- dist_mat[i,]
    
      # Find nearest non-NA neighbor
      valid_indices <- which(!is.na(non_na_pts[[value_col]]))
      nearest_idx <- valid_indices[which.min(distances[valid_indices])]
    
      # Fill NA with nearest neighbor value
      na_pts[[value_col]][i] <- non_na_pts[[value_col]][nearest_idx]
    }
  }
  
  if(!is.null(date_col)){
     for (i in 1:nrow(na_pts)) {
      # Get the point and date for current NA value
      current_point <- na_pts[i,]
      current_date <- na_pts[[date_col]][i]

      # Get just those non_na points with same date
      non_na_pts_sub <- non_na_pts |>
        dplyr::filter(Date == current_date)
      
      # Distances
      distances<- st_distance(current_point, non_na_pts_sub)

      # Calculate time differences (in days)
      time_diffs <- abs(as.numeric(difftime(current_date, 
                                        non_na_pts_sub[[date_col]], 
                                        units = "days")))
    
      # Combine spatial and temporal distances (normalized)
      space_weight <- 0.5
      time_weight <- 0.5
    
      norm_distances <- as.numeric(distances) / max(as.numeric(distances))
      norm_time_diffs <- time_diffs / max(time_diffs)

      if(all(is.na(norm_time_diffs))){
        combined_distance <- space_weight * norm_distances
      } else {
        combined_distance <- space_weight * norm_distances + 
                        time_weight * norm_time_diffs
      }
      
      # Find nearest non-NA neighbor
      valid_indices <- which(!is.na(non_na_pts_sub[[value_col]]))
      nearest_idx <- valid_indices[which.min(combined_distance[valid_indices])]
    
      # Fill NA with nearest neighbor value
      na_pts[[value_col]][i] <- non_na_pts_sub[[value_col]][nearest_idx]
    } 
  }
  return(na_pts)
}

####  Functions  ####
####
## Dynamic loading function -- each takes in a directory path and reads in files from it. 
static_covariates_read<- function(static_covariates_dir){
  
  # Get the files...
  rast_read_layers<- list.files(static_covariates_dir, pattern = ".grd")
  
  # Output storage
  static_covariate_rasters_out<- vector(mode = "list", length = length(rast_read_layers))
  
  # Loop through each, save em in the list
  for(i in seq_along(rast_read_layers)){
    static_covariate_rasters_out[[i]]<- raster::stack(paste(static_covariates_dir, rast_read_layers[i], sep = "/"))
    names(static_covariate_rasters_out)[i]<- file_path_sans_ext(rast_read_layers[i])
  }
  
  # Return it
  return(static_covariate_rasters_out)
}

dynamic_covariates_read<- function(dynamic_covariates_dir){
  
  # Get the files...
  rast_read_stacks<- list.files(dynamic_covariates_dir, pattern = ".grd")
  
  # Output storage
  dynamic_covariate_rasters_out<- vector(mode = "list", length = length(rast_read_stacks))
  
  # Loop through each, save em in the list
  for(i in seq_along(rast_read_stacks)){
    dynamic_covariate_rasters_out[[i]]<- raster::stack(paste(dynamic_covariates_dir, rast_read_stacks[i], sep = "/"))
    names(dynamic_covariate_rasters_out)[i]<- file_path_sans_ext(rast_read_stacks[i])
  }
  
  # Return it
  return(dynamic_covariate_rasters_out)
}

# Rescaling covariates
covariate_rescale_func<- function(x, type, center = NULL, scale = NULL){
  if(type == "JT"){
    x.out<- x/100
    return(x.out)
  } else {
    x.out<- as.numeric(scale(abs(x), center = center, scale = scale))
    return(x.out)
  }
}

## Processing functions
#' @title Extract static covariates at point locations
#' 
#' @description This function reads a raster layer and a sf spatial points object. The function then does a simple extraction by overlaying the spatial points onto the raster layer and gathering the value for each raster layer in the stack at all of the point locations. The function then returns the original sf_points data set as either an sf object or dataframe with columns for the raster layer as a tibble and also saves this file as "model_covs.rds".
#'
#' @param rast = A raster layer of one static covariate variable 
#' @param cov_name = The column name to use for the extracted value for the static covariate 
#' @param sf_points = SF spatial points object specifying the locations where we want to extract raster stack values
#' @param date_col_name = The column name of sf_points with the date information, formatted as YYYY-MM-DD
#' @param df_sf = Character string of "df" or "sf" signaling whether the returned object should be a dataframe or sf object 
#' 
#' @return Either an SF object or dataframe, with the covariate value appended as a new column to the point locations file. This file is also saved in out_dir. 
#' 
#' @export
static_extract<- function(rast, cov_name, sf_points, interp_nas, date_col_name, df_sf){
  # For debugging
  if(FALSE){
    rast = rast_use
    cov_name = names(static_covariates_list)[[i]]
    sf_points = sf_points_run
    date_col_name = date_col_name
    df_sf = "sf"
  }
  
  # A few checks...
  # Do the projections match?
  if(!st_crs(rast) == st_crs(sf_points)){
    print("Check that projection in raster stack and spatial points match")
    stop()
  }
  
  # # Are there duplicate records?
  # if(any(duplicated(sf_points, by = c(EST_DATE, geometry)))){
  #   print("Check `sf_points` and remove duplicated observations to reduce extraction time")
  #   stop()
  # }
  
  # Extraction -----------------------------------------------------------
  sf_extract<- data.frame(raster::extract(rast, sf_points))
  names(sf_extract)<- cov_name
  
  # Bind to sf_unique
  sf_points <- cbind(sf_points, sf_extract)
  
   # Deal with NAs
  if(interp_nas){
    print("WARNING! You have elected to fill in NA covariate values based on their nearest, non-NA match in space and time. We recommend first running this with `interp_na = FALSE` to make sure you know what you are doing.")
    na_row_ids <- which(is.na(sf_points[cov_name]))
    na_points <- st_transform(sf_points[na_row_ids, ], 32619)
    non_na_points <- st_transform(sf_points[-na_row_ids, ], 32619)
    
    # Get distance matrix, just do this once
    # Calculate distances to all points in source dataset
    # distances <- st_distance(na_points, non_na_points) # This takes a longg time
    distances <- as.matrix(dist(st_coordinates(na_points), st_coordinates(non_na_points), method = "euclidean"))
    
    filled_points <- fill_na_spatiotemporal(na_pts = na_points, non_na_pts = non_na_points, dist_mat = distances, value_col = cov_name, date_col = NULL)

    # Put back in
    sf_points[na_row_ids, cov_name] <- st_drop_geometry(filled_points[, cov_name])
  }
  
  # Return processed file -----------------------------------------------------------
  if(df_sf == "sf"){
    # Keep it as sf object
    out<- sf_points
    return(out)
  } else {
    # Drop the geometry and save the data frame
    out<- st_drop_geometry(sf_points)
    return(out)
  }
  
}

#' @title Static extract at point locations wrapper
#' 
#' @description This short wrapper function is used to process all potential covariates in the static_covariates_dir.
#'
#' @param static_covariates_list = List with static raster layers
#' @param all_tows = All tows dataframe, created by `bind_nmfs_dfo_tows()`
#' @param out_dir = Directory to save the combined dataframe as an .rds file
#' 
#' @return A datafame with information of all unique tows and appended columns for each of the static habitat covariates. This file is also saved in out_dir. 
#' 
#' @export
static_extract_wrapper<- function(static_covariates_list, sf_points, interp_nas, date_col_name, df_sf, out_dir){
  # For debugging
  if(FALSE){
    tar_load(static_covariates_stack)
    tar_load("all_tows_sf")
    static_covariates_list = static_covariates_stack
    sf_points = all_tows_sf
    date_col_name = "EST_DATE"
    df_sf = "sf"
    out_dir
  }
  
  # For each of the covariate layers in `static_covariates_list` we are going to run our `static_extract` function. I think this might require some creative manipulating of the `sf_points` so we don't lose covariates as we move through them...
  sf_points_run<- sf_points
  for(i in seq_along(static_covariates_list)){
    rast_use<- static_covariates_list[[i]]
    temp_tows_with_covs<- static_extract(rast = rast_use, cov_name = names(static_covariates_list)[[i]], sf_points = sf_points_run, interp_nas = interp_nas, date_col_name = date_col_name, df_sf = "sf")
    
    # If there are more files to go, update sf_points_run
    if(i < length(static_covariates_list)){
      sf_points_run<- temp_tows_with_covs
    }
    
    # If "last" one, good to return it
    if(i == length(static_covariates_list)){
      tows_with_covs_out<- temp_tows_with_covs
    }
  }
  
  if(df_sf == "sf"){
    # Keep it as sf object
    out<- tows_with_covs_out
    if(!is.null(out_dir)){
      saveRDS(out, file = paste(out_dir, "all_tows_with_static_covs_sf.rds", sep = "/"))
    }
    return(out)
  } else {
    # Drop the geometry and save the data frame
    out<- st_drop_geometry(tows_with_covs_out)
    if(!is.null(out_dir)){
      saveRDS(out, file = paste(out_dir, "all_tows_with_static_covs.rds", sep = "/"))
    }
    return(out)
  }
}

#' @title Extract dynamic covariates at point locations
#' 
#' @description This function summarizes spatio-temporal 2d gridded data at sample point locations. The function takes as the primary inputs a stack of one dynamic variable and a sf spatial points object. The function then does the extraction by overlaying the spatial points onto the raster stack and gathering the value for each raster layer in the stack at all of the point locations. Then, depending on the t_summ and t_position arguments, the function summarizes the fully extracted data (all points all locations) at the appropriate temporal scale and position (leading up to, saddling, or in the future) relative to the observation time. The function returns either the original sf_points data set with a new column for the dynamic raster layer values or a dataframe of the same structure and also saves this file either as "model_covs.rds" if new_file_name is NULL or as the name provided in new_file_name.
#'
#' @param rast_ts_stack = A raster stack of one covariate variable, where each layer represents a different time step
#' @param stack_name = The column name to use for the extracted value for the dynamic covariate 
#' @param t_summ = A numeric value or character string that indicates what temporal resolution should be used in summarizing the covariate values. If numeric, the function will simply summarize the values in the raster stack from the matching period back `t_summ` numeric time steps. For example, if the `rast_ts_stack` provided daily values and t_summ = 90, the function would calculate a 90 day average, where the 90-day window would either be leading up to and including the day of observation, saddled around the observation, or include the day of the observation and 89 days into the future. If a character string, should be one of "daily", monthly", "seasonal", "annual", or "match". These options are built into the function to provide a bit easier specification to quickly calculate the monthly/seasonal/annual summaries of the raster stack values. When used, this automatically defines t_position = saddle.
#' @param t_position = A character vector of either NULL, "past", "saddle", or "future". If NULL, then values are extracted based on matching up the observation point with the dynamic raster stack at the level specified in the t_summ character vector (e.g., "daily", "monthly", "seasonal", "annual", "match". If not, then summaries are calculated leading up to and including the observation time ("past"), saddled around around the observation time ("saddle"), or including and in the future of the observation time ("future").
#' @param interp_nas = TRUE/FALSE logical value specifying whether observations with "NA" covariate values should take on the value of the nearest observation (in space and time) that has a non-NA value. This is implemented to account, in particular, for observations along coast and island shorelines where we might not have environmental data from regular gridded products.
#' @param sf_points = SF spatial points object specifying the locations where we want to extract raster stack values
#' @param date_col_name = The column name of sf_points with the date information, formatted as YYYY-MM-DD
#' @param df_sf = Character string of "df" or "sf" signaling whether the returned object should be a dataframe or sf object 
#' 
#' @return Either an SF object or dataframe, with the covariate values appended as new columns to the point locations file. This file is also saved in out_dir. 
#' 
#' @export

dynamic_2d_extract<- function(rast_ts_stack, stack_name, t_summ, t_position, interp_nas, sf_points, date_col_name, df_sf){
  
  # Custom rowMeans like function
  rowMeans_cust<- function(row_id, start_col, end_col, full_data){
    
    # If daily, then start_col and end_col are equal, so just use one of them
    if(start_col == end_col){
      out<- full_data[row_id, start_col]
      return(out)
    } else {
      out<- rowMeans(full_data[row_id,start_col:end_col])
      return(out)
    }
  }
  
  # For debugging
  if(FALSE){
    rast_ts_stack = dynamic_stacks[[3]]
    stack_name = unlist(lapply(dynamic_files, names))[3]
    t_summ = "match"
    t_position = NULL
    interp_nas = TRUE
    sf_points = sf_points_run
    date_col_name = "Date"
    df_sf = "df"
  }
  
  # A few checks...
  # If t_summ is a character, does it match one of "daily", monthly", "seasonal", or "annual" AND is t_position set to NULL?
  if(is.character(t_summ)){
    t_summ_check<- t_summ %in% c("daily", "monthly", "seasonal", "annual", "match") & is.null(t_position)
    if(!t_summ_check){
      print("Check 't_summ' argument and 't_position'. 't_summ' must be one of 'daily', monthly', 'seasonal', 'annual', 'match' and 't_position' = NULL")
      stop()
    }
  }
  
  # Check t_position argument
  t_position_check<- is.null(t_position) || t_position %in% c("past", "saddle", "future") 
  if(!t_position_check){
    print("Check 't_position' argument. Must be one of 'past', 'saddle' or 'future'")
    stop()
  }
  
  # # Are there duplicate records?
  # if(any(duplicated(sf_points, by = c(EST_DATE, geometry)))){
  #   print("Check `sf_points` and remove duplicated observations to reduce extraction time")
  #   stop()
  # }
  
  # Finally, do the projections match?
  if(!st_crs(rast_ts_stack) == st_crs(sf_points)){
    print("Check that projection in raster stack and spatial points match")
    stop()
  }
  
  # To do the matching for seasons, need some type of look up table.
  month_season_table<- data.frame("Month" = str_pad(seq(from = 1, to = 12, by = 1), 2, "left", 0), "Season" = c("Winter", "Winter", "Spring", "Spring", "Spring", "Summer", "Summer", "Summer", "Fall", "Fall", "Fall", "Winter"))
  
  # Rename the date column to be "EST_DATE" to match the function set up, This gets a bit tricky if we are running things in sequence. In other words, maybe we ran things through for SST and then want to do Chl. When we pass in the dataframe, "EST_DATE" may already exist. Going to try something else for now and pass in {{date_col_name}} below whenever EST_DATE is used
  # sf_points<- sf_points %>%
  #   rename(., "EST_DATE" := {{date_col_name}})

  # Full extraction, all points and layers -----------------------------------------------------------
  sf_extract <- data.frame(raster::extract(rast_ts_stack, sf_points))
  
  # Getting the values we want from full extraction -------------------------
  # This leaves us with a data frame where the rows are the sf_points and then the columns are the covariate value for a given location for each of the layers (i.e., time steps) in the stack. To start with calculating any summary, we first need to know the which column provides an exact match between the time of the observation and the time in the raster stack. Then, depending on if it is a past, saddle, or future, and how many columns are specified by t_pos, we can average the necessary columns.
  colnames(sf_extract)<- gsub("[.]", "-", gsub("X", "", colnames(sf_extract)))
  
  # Now, accounting for the different windows and positions. As mentioned earlier, this is going to be completely different depending on if t_summ is a character sting ("daily", "monthly", "seasonal", "annual") or if it is a numeric integer. So, here comes that split...
  summ_df_list<- vector("list", length(t_summ))
  
  #####
  ## t_summm as a character string
  #####
  if(is.character(t_summ)){
    
    for(i in seq_along(t_summ)){
      # Get t_summ for this iteration
      t_summ_use<- t_summ[i]
      
      # Create a data frame to store the start and end...
      summ_df<- data.frame("Point" = 1:nrow(sf_extract), "Date_Match" = rep(NA, nrow(sf_extract)), "Start_Summ" = rep(NA, nrow(sf_extract)), "End_Summ" = rep(NA, nrow(sf_extract)))
      
      # Get the exact date match
      sf_points_date_col<- as.character(st_drop_geometry(sf_points[,{{date_col_name}}])[,1])
      sf_points_date_col<- switch(t_summ_use,
                                  "match" = ymd(sf_points_date_col),
                                  "daily" = ymd(sf_points_date_col),
                                  "monthly" = format(as.Date(sf_points_date_col), "%Y-%m"),
                                  "seasonal" = format(as.Date(sf_points_date_col), "%Y-%m"),
                                  "annual" = format(as.Date(sf_points_date_col), "%Y"))
      
      # May need to revisit this at some point, but for now, shouldn't influence the results..
      if(t_summ_use == "daily" | t_summ_use == "match"){
        summ_df$Date_Match<- match(sf_points_date_col, as.Date(colnames(sf_extract)))
      } 
      if(t_summ_use == "monthly" | t_summ_use == "seasonal"){
        summ_df$Date_Match<- match(sf_points_date_col, format(as.Date(colnames(sf_extract)), "%Y-%m"))
      }
      if(t_summ_use == "annual"){
        summ_df$Date_Match<- match(sf_points_date_col, format(as.Date(colnames(sf_extract)), "%Y"))
      }
      
      # Column index start based on t_summ_use. Seasonal needs more finagling.
      if(t_summ_use != "seasonal"){
        # Might need to come back and check this!
        rast_ts_dates<- as.Date(gsub("X", "", gsub("[.]", "-", names(rast_ts_stack))))
        if(t_summ_use == "match"){
          summ_df$Start_Summ<- summ_df$Date_Match
          summ_df$End_Summ<- summ_df$Date_Match
        } 
        if(t_summ_use == "daily"){
          summ_df$Start_Summ<- match(sf_points_date_col, as.Date(colnames(sf_extract)))
          summ_df$End_Summ<- match(sf_points_date_col, as.Date(colnames(sf_extract)))
        } 
        if(t_summ_use == "monthly" | t_summ_use == "seasonal"){
          summ_df$Start_Summ<- sapply(sf_points_date_col, FUN = function(x) min(which(x == format(rast_ts_dates, "%Y-%m"))))
          summ_df$End_Summ<- sapply(sf_points_date_col, FUN = function(x) max(which(x == format(rast_ts_dates, "%Y-%m"))))
        }
        if(t_summ_use == "annual"){
          summ_df$Start_Summ<- sapply(sf_points_date_col, FUN = function(x) min(which(x == format(rast_ts_dates, "%Y"))))
          summ_df$End_Summ<- sapply(sf_points_date_col, FUN = function(x) max(which(x == format(rast_ts_dates, "%Y"))))
        }
      } else {
        # Season bits, need a year - season combo for both the points AND columns of sf_extract.
        colnames_orig<- as.Date(gsub("X", "", gsub("[.]", "-", colnames(sf_extract))))
        colnames_season<- month_season_table$Season[match(format(colnames_orig, "%m"), month_season_table$Month)]
        colnames_season_match<- paste(format(colnames_orig, "%Y"), colnames_season, sep = "-")
        
        sf_points$Season_Match<- paste(format(as.Date(paste(sf_points_date_col, "-01", sep = "")), "%Y"), month_season_table$Season[match(format(as.Date(paste(sf_points_date_col, "-01", sep = "")), "%m"), month_season_table$Month)], sep = "-")
        
        # Start index
        summ_df$Start_Summ<- match(sf_points$Season_Match, colnames_season_match)
        
        # Now index end
        summ_df$End_Summ<- sapply(sf_points$Season_Match, FUN = function(x) max(which(x == colnames_season_match)))
      }
      
      # Deal with the "-Inf" indices, arising when there are no matches for "max"
      summ_df[summ_df == -Inf] <- NA
      
      # Store it
      summ_df_list[[i]]<- summ_df
      names(summ_df_list)[i]<- paste(stack_name, t_summ_use, sep = "_")
    }
  } 
  
  #####
  ## t_summ as a numeric vector
  #####
  if(is.numeric(t_summ)){
    for(i in seq_along(t_summ)){
      # Create a data frame to store the start and end...
      summ_df<- data.frame("Point" = 1:nrow(sf_extract), "Date_Match" = rep(NA, nrow(sf_extract)), "Start_Summ" = rep(NA, nrow(sf_extract)), "End_Summ" = rep(NA, nrow(sf_extract)))
      
      # Get the exact date match
      summ_df$Date_Match<- match(sf_points_date_col, as.Date(colnames(sf_extract)))
      
      # Get summary window time range
      t_summ_use<- t_summ[i]
      
      # Column index start based on the t_position and t_summ_use
      summ_df$Start_Summ<- switch(t_position,
                                  "past" = summ_df$Date_Match - (t_summ_use),
                                  "saddle" = summ_df$Date_Match - (t_summ_use/2),
                                  "future" = summ_df$Date_Match)
      # Now index end
      summ_df$End_Summ<- switch(t_position,
                                "past" = summ_df$Date_Match,
                                "saddle" = summ_df$Date_Match + (t_summ_use/2),
                                "future" = summ_df$Date_Match + (t_summ_use))
      
      # Store it
      summ_df_list[[i]]<- summ_df
      names(summ_df_list)[i]<- paste(stack_name, t_summ_use, sep = "_")
    }
  }
  
  # Subsetting full extraction, taking row means across the columns signaled by Start_Summ and End_Summ. First, need to drop any NAs as those are going to through an error in the indexing AND any negative index numbers (this would indicate a situation where the date of observation falls during the beginning of the dynamic variable time series, such that there isn't data for the prior or saddle option t_summ_use time steps prior to the observation)
  # Has to be a better option for how to do this, but a loop for now..
  # Adding Point ID column, used for merging
  sf_points$Point<- seq(from = 1, to = nrow(sf_points))
  
  for(j in seq_along(summ_df_list)){
    
    summ_df_comp<- summ_df_list[[j]] %>%
      drop_na(., Start_Summ, End_Summ) %>%
      dplyr::filter(., Start_Summ > 0)
    
    # Calculate mean given start and end column index of extraction file and row based on observation point ID
    sf_extract_mean <- summ_df_comp %>%
      mutate(., "Summ_Val" = pmap_dbl(list(row_id = Point, start_col = Start_Summ, end_col = End_Summ, full_data = list(sf_extract)), rowMeans_cust))
    
    # Don't need to keep all the columns, just point and summ_val. Rename to match stack and t_summ
    point_mean<- sf_extract_mean %>%
        dplyr::select(., Point, Summ_Val)
    names(point_mean)[2]<- names(summ_df_list)[j]

    # Join back to original trawl sf_points
    sf_points<- sf_points %>%
      left_join(., point_mean, by = c("Point" = "Point")) %>%
      dplyr::select(., -Point)
  }

  # Deal with NAs
  if(interp_nas){
    print("WARNING! You have elected to fill in NA covariate values based on their nearest, non-NA match in space and time. We recommend first running this with `interp_na = FALSE` to make sure you know what you are doing.")
    na_row_ids<- which(is.na(sf_points[names(summ_df_list)[j]]))
    na_points <- st_transform(sf_points[na_row_ids, ], 32619)
    non_na_points <- st_transform(sf_points[-na_row_ids, ], 32619)
    
    # Get distance matrix, just do this once
    # Calculate distances to all points in source dataset
    # distances <- st_distance(na_points, non_na_points) # This takes a longg time
    # distances <- as.matrix(dist(st_coordinates(na_points), st_coordinates(non_na_points), method = "euclidean"))
    # distances <- st_distance(na_points, non_na_points, tolerance = 200)

    filled_points <- fill_na_spatiotemporal(na_pts = na_points, non_na_pts = non_na_points, value_col = names(summ_df_list)[j], date_col = "Date")

    # Put back in
    sf_points[na_row_ids, names(summ_df_list)[j]] <- st_drop_geometry(filled_points[, names(summ_df_list)[j]])
  }
  
  # Save and return processed file -----------------------------------------------------------
  if(df_sf == "sf"){
    # Keep it as sf object
    out<- sf_points
    return(out)
  } else {
    # Drop the geometry and save the data frame
    out<- st_drop_geometry(sf_points)
    return(out)
  }
  #########
  ## End
  #########
}

#' @title Extract dynamic covariates at point locations wrapper
#' 
#' @description This short wrapper function is used to process all potential covariates in the dynamic_covariates_dir.
#'
#' @param dynamic_covariates_list = 
#' @param t_summ
#' @param t_position 
#' @param sf_points
#' @param date_col_name
#' @param df_sf
#' @param out_dir = Directory to save the combined dataframe as an .rds file
#' 
#' @return A datafame with information of all unique tows with dynamic covariate values as new columns. This file is also saved in out_dir. 
#' 
#' @export
dynamic_2d_extract_wrapper<- function(dynamic_covariates_list, t_summ, t_position, interp_nas, sf_points, date_col_name, df_sf, out_dir){
  # For debugging
  if(FALSE){
    tar_load(dynamic_covariates_stack)
    dynamic_covariates_list = dynamic_covariates_stack
    t_summ = "seasonal"
    t_position = NULL
    tar_load(all_tows_with_static_covs)
    sf_points = all_tows_with_static_covs
    date_col_name = "DATE"
    df_sf ="df"
    out_dir = here::here("data/combined")
    
    dynamic_covariates_list = dynamic_covariates_list_use
    t_summ = "monthly"
    t_position = NULL
    sf_points = tows_sf
    date_col_name = "Date"
    df_sf = "sf"

    dynamic_covariates_list = dynamic_stacks[c(2, 5, 8, 11)]
    t_summ = "match"
    t_position = NULL
    interp_nas = TRUE
    sf_points = all_tows_with_static_covs
    date_col_name = "Date"
    df_sf = "df"
    out_dir = NULL
  }
  
  # For each of the covariate layers in `dynamic_covariates_list` we are going to run our `dynamic_2d_extract` function. I think this might require some creative manipulating of the `sf_points` so we don't lose covariates as we move through them...
  sf_points_run<- sf_points
  
  for(i in seq_along(dynamic_covariates_list)){
    stack_use<- dynamic_covariates_list[[i]]
    temp_tows_with_covs<- dynamic_2d_extract(rast_ts_stack = stack_use, stack_name = names(dynamic_covariates_list)[[i]], t_summ = t_summ, t_position = t_position, interp_nas = interp_nas, sf_points = sf_points_run, date_col_name = date_col_name, df_sf = "sf")
    
    # If there are more files to go, update sf_points_run
    if(i < length(dynamic_covariates_list)){
      sf_points_run<- temp_tows_with_covs
    }
    
    # If "last" one, good to return it
    if (i == length(dynamic_covariates_list)) {
      tows_with_covs_out <- temp_tows_with_covs
    }
    print(i)
  }
  
  if(df_sf == "sf"){
    # Keep it as sf object
    out<- tows_with_covs_out
    if(!is.null(out_dir)){
      saveRDS(out, file = paste(out_dir, "all_tows_all_covs_sf.rds", sep = "/"))
    }
    return(out)
  } else {
    # Drop the geometry and save the data frame
    out<- st_drop_geometry(tows_with_covs_out)
    if(!is.null(out_dir)){
      saveRDS(out, file = paste(out_dir, "all_tows_all_covs.rds", sep = "/"))
    }
    return(out)
  }
}

# New rescale function
rescale_all_covs<- function(all_tows_with_all_covs, depth_cut = depth_cut, covs_rescale = c("Depth", "BS_seasonal", "BT_seasonal", "SS_seasonal", "SST_seasonal"), type = "AJA", scale = TRUE, center = TRUE, out_dir){
  
  if(FALSE){
    tar_load(all_tows_with_all_covs)
    type = "AJA"
    scale = TRUE
    center = TRUE
    out_dir
  }
  
  # Cut depth
  all_tows_with_all_covs_temp<- all_tows_with_all_covs %>%
    dplyr::filter(., Depth  <= depth_cut)
  
  all_tows_with_all_covs_rescale<- all_tows_with_all_covs_temp %>%
    mutate_at(., .vars = {{covs_rescale}}, .funs = covariate_rescale_func, type = type, scale = scale, center = center)
  
  saveRDS(all_tows_with_all_covs_rescale, file = paste(out_dir, "all_tows_all_covs_rescale.rds", sep = "/"))
  return(all_tows_with_all_covs_rescale)
}

get_rescale_params<- function(all_tows_with_all_covs, depth_cut, out_dir){
  
  if(FALSE){
    tar_load(all_tows_with_all_covs)
    type = "AJA"
    scale = TRUE
    center = TRUE
  }
  
  # Cut depth
  all_tows_with_all_covs_temp<- all_tows_with_all_covs %>%
    dplyr::filter(., Depth  <= depth_cut)
  
  cov_names<- c("Depth", "BS_seasonal", "BT_seasonal", "SS_seasonal", "SST_seasonal")
  
  rescale_params<- all_tows_with_all_covs %>%
    summarize_at(., .vars = {{cov_names}}, .funs = c("Mean" = mean, "SD" = sd), na.rm = TRUE)
  
  saveRDS(rescale_params, file = paste(out_dir, "rescale_cov_params.rds", sep = "/"))
  return(rescale_params)
}

#' @title Aggregate covariate raster stack list
#' 
#' @description This function aggregates a covariate raster stack to a new time scale. For example, with year-month raster stack of SST, we may instead need the covariates summarized at year-season time scale.
#'
#' @param predict_covariates_dir = predict_covariates_dir
#' @param ensemble_stat = Either the climate model ensemble statistic to use when working with climate model projections, or NULL. This is only used in naming the output file
#' @param resample_template = Raster layer template with target resolution
#' @param summarize = Currently, either "annual" or "seasonal" to indicate whether the each dynamic raster stack should be summarized to an annual or seasonal time scale
#' @param out_dir = Directory to save the aggregated raster stack
#' 
#' @return  The output directory path. This allows the function to integrate easily within the Targets workflow. Along with returning the output directory, the function also saves the new aggregated raster stack for each covariate.
#' 
#' @export
predict_covariates_stack_agg<- function(predict_covariates_dir, ensemble_stat, resample_template, summarize, out_dir){
  if(FALSE){
    predict_covariates_dir = "~/Box/RES_Data/CMIP6/BiasCorrected"
    ensemble_stat = "mean"
    resample_template = here::here("scratch/aja/TargetsSDM/data/supporting/", "Rast0.25grid.grd")
    summarize = "seasonal"
    out_dir = here::here("scratch/aja/TargetsSDM/data/predict")
  }
  
  # Quick check, ensemble stat needs to be one of mean, 5th or 95th
  if(!ensemble_stat %in% c("mean", "5th", "95th")){
    stop("Check `ensemble stat`. Must be one of `mean`, `5thpercentile` or `95thpercentile`.")
  }
  
  # To do the matching for seasons, need some type of look up table.
  month_season_table<- data.frame("Month" = str_pad(seq(from = 1, to = 12, by = 1), 2, "left", 0), "Season" = c("Winter", "Winter", "Spring", "Spring", "Spring", "Summer", "Summer", "Summer", "Fall", "Fall", "Fall", "Winter"))
  
  # First, get the files that match our ensemble statistic
  rast_files_load<- list.files(predict_covariates_dir, pattern = paste0(ensemble_stat, ".grd"), full.names = TRUE)
  
  for(i in seq_along(rast_files_load)){
    
    # Bring in the stack
    rast_stack_temp<- raster::stack(rast_files_load[i])
    
    covariate_name<- case_when(
      grepl("bot_sal", rast_files_load[i]) ~ "BS",
      grepl("bot_temp", rast_files_load[i]) ~ "BT",
      grepl("surf_sal", rast_files_load[i]) ~ "SS",
      grepl("surf_temp", rast_files_load[i]) ~ "SST"
    )
    
    # Check layers
    layers<- nlayers(rast_stack_temp)
    
    if(layers > 1){
      # Check resolution, resample if necessary
      if(any(res(rast_stack_temp) != c(0.25, 0.25))){
        rast_template<- raster::stack(resample_template)
        rast_stack_temp<- resample(rast_stack_temp, rast_template, method = "bilinear")
      }
      # Need to summarize
      # Get indices
      if(summarize == "annual"){
        stack_dates<- as.Date(gsub("X", "", gsub("[.]", "-", paste(names(rast_stack_temp), ".16", sep = ""))))
        stack_years<- format(stack_dates, format = "%Y")
        indices<- as.numeric(factor(stack_years, levels = unique(stack_years)))
        stack_names<- unique(stack_years)
      } else if(summarize == "seasonal"){
        stack_dates<- as.Date(gsub("X", "", gsub("[.]", "-", paste(names(rast_stack_temp), ".16", sep = ""))))
        stack_season<- month_season_table$Season[match(format(stack_dates, "%m"), month_season_table$Month)]
        stack_season_years<- paste(format(stack_dates, format = "%Y"), stack_season, sep = "-")
        indices<- as.numeric(factor(stack_season_years, levels = unique(stack_season_years)))
        stack_names<- unique(stack_season_years)
      }
      # Calculate summary and store it
      rast_agg_temp<- stackApply(rast_stack_temp, indices = indices, fun = mean, na.rm = TRUE)
      names(rast_agg_temp)<- stack_names
      writeRaster(rast_agg_temp, filename = paste(out_dir, "/predict_stack_", covariate_name, "_", summarize, "_", ensemble_stat, ".grd", sep = ""), format = "raster", overwrite = TRUE)
    } else {
      # No summary needed, just check resolution and then save it
      if(any(res(rast_stack_temp) != c(0.25, 0.25))){
        rast_template<- raster::stack(resample_template)
        rast_stack_temp<- resample(rast_stack_temp, rast_template, method = "bilinear")
      }
      rast_agg_temp<- rast_stack_temp
      writeRaster(rast_agg_temp, filename = paste(out_dir, "/predict_stack_", covariate_name, "_", summarize, "_", ensemble_stat, ".grd", sep = ""), format = "raster", overwrite = TRUE)
    }
  }
  # Return the file path
  return(out_dir)
}