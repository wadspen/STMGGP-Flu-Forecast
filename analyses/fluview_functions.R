#all functions here besides fetch_flu() were taken from the old cdcflu

.region_map <- c(
  national = 3,
  hhs = 1,
  census = 2,
  state = 5
)

.cdcfluview_ua <- "Mozilla/5.0 (compatible; R-cdcvluview Bot/2.0; https://github.com/hrbrmstr/cdcfluview)"


.mcga <- function (tbl) 
{
  x <- colnames(tbl)
  x <- tolower(x)
  x <- gsub("[[:punct:][:space:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("(^_|_$)", "", x)
  x <- gsub("^x_", "", x)
  x <- make.unique(x, sep = "_")
  colnames(tbl) <- x
  tbl
}


to_num <- function (x) 
{
  x <- gsub("%", "", x, fixed = TRUE)
  x <- gsub(">", "", x, fixed = TRUE)
  x <- gsub("<", "", x, fixed = TRUE)
  x <- gsub(",", "", x, fixed = TRUE)
  x <- gsub(" ", "", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}


ilinet <- function (region = c("national", "hhs", "census", "state"), years = NULL) 
{
  region <- match.arg(tolower(region), c("national", "hhs", 
                                         "census", "state"))
  meta <- jsonlite::fromJSON("https://gis.cdc.gov/grasp/flu2/GetPhase02InitApp?appVersion=Public")
  params <- list(AppVersion = "Public", DatasourceDT = list(list(ID = 1, 
                                                                 Name = "ILINet")), RegionTypeId = .region_map[region])
  params$SubRegionsDT <- switch(region, national = {
    list(list(ID = 0, Name = ""))
  }, hhs = {
    lapply(1:10, function(i) list(ID = i, Name = as.character(i)))
  }, census = {
    lapply(1:9, function(i) list(ID = i, Name = as.character(i)))
  }, state = {
    lapply(1:59, function(i) list(ID = i, Name = as.character(i)))
  })
  available_seasons <- sort(meta$seasons$seasonid)
  if (is.null(years)) {
    years <- available_seasons
  }
  else {
    years <- as.numeric(years)
    years <- ifelse(years > 1996, years - 1960, years)
    years <- sort(unique(years))
    years <- years[years %in% available_seasons]
    if (length(years) == 0) {
      years <- rev(sort(meta$seasons$seasonid))[1]
      curr_season_descr <- meta$seasons[meta$seasons$seasonid == 
                                          years, "description"]
      message(sprintf("No valid years specified, defaulting to this flu season => ID: %s [%s]", 
                      years, curr_season_descr))
    }
  }
  params$SeasonsDT <- lapply(years, function(i) list(ID = i, 
                                                     Name = as.character(i)))
  tf <- tempfile(fileext = ".zip")
  td <- tempdir()
  on.exit(unlink(tf), TRUE)
  res <- httr::POST(url = "https://gis.cdc.gov/grasp/flu2/PostPhase02DataDownload", 
                    httr::user_agent(.cdcfluview_ua), httr::add_headers(Origin = "https://gis.cdc.gov", 
                                                                        Accept = "application/json, text/plain, */*", Referer = "https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html"), 
                    encode = "json", body = params, httr::write_disk(tf))
  httr::stop_for_status(res)
  nm <- unzip(tf, overwrite = TRUE, exdir = td)
  xdf <- read.csv(nm, skip = 1, stringsAsFactors = FALSE)
  xdf <- .mcga(xdf)
  xdf$weighted_ili <- to_num(xdf$weighted_ili)
  xdf$unweighted_ili <- to_num(xdf$unweighted_ili)
  xdf$age_0_4 <- to_num(xdf$age_0_4)
  xdf$age_25_49 <- to_num(xdf$age_25_49)
  xdf$age_25_64 <- to_num(xdf$age_25_64)
  xdf$age_5_24 <- to_num(xdf$age_5_24)
  xdf$age_50_64 <- to_num(xdf$age_50_64)
  xdf$age_65 <- to_num(xdf$age_65)
  xdf$ilitotal <- to_num(xdf$ilitotal)
  xdf$num_of_providers <- to_num(xdf$num_of_providers)
  xdf$total_patients <- to_num(xdf$total_patients)
  xdf$week_start <- MMWRweek::MMWRweek2Date(xdf$year, xdf$week)
  if (region == "national") 
    xdf$region <- "National"
  if (region == "hhs") 
    xdf$region <- factor(xdf$region, levels = sprintf("Region %s", 
                                                      1:10))
  class(xdf) <- c("tbl_df", "tbl", "data.frame")
  arrange(suppressMessages(readr::type_convert(xdf)), week_start)
}







#cdcepi function used to collect NHSN weekly data
fetch_flu <- function(){
  require(dplyr)
  require(lubridate)
  require(jsonlite)
  require(purrr)
  require(stringr)

  locations <- read.csv(file = "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/locations.csv") %>% dplyr::select(1:4)
  #require(RSocrata)
  
  #read data from data.cdc.gov, filtering for when flu reporting became mandatory
  # health_data = RSocrata::read.socrata(url = "https://data.cdc.gov/resource/mpgq-jmmr.json") %>% 
  #   dplyr::filter(weekendingdate >= as.Date("2022-02-01"))
  
  # Function to fetch paginated data
  # fetch_socrata_data <- function(base_url, limit = 1000, max_rows = 30000) {
  #   offsets <- seq(0, max_rows, by = limit)
  #   
  #   all_data <- map_dfr(offsets, function(offset) {
  #     paged_url <- paste0(base_url, "?$limit=", limit, "&$offset=", offset)
  #     message("Fetching offset: ", offset)
  #     tryCatch(
  #       fromJSON(paged_url),
  #       error = function(e) {
  #         message("Error at offset ", offset, ": ", e$message)
  #         NULL
  #       }
  #     )
  #   })
  #   
  #   return(all_data)
  # }
  fetch_socrata_data <- function(base_url, limit = 1000, max_rows = 30000, filter_query = NULL) {
    offsets <- seq(0, max_rows, by = limit)
    
    all_data <- map_dfr(offsets, function(offset) {
      # Construct the paged URL with optional filtering
      paged_url <- paste0(base_url, "?$limit=", limit, "&$offset=", offset)
      
      # Add filter query if provided
      if (!is.null(filter_query)) {
        paged_url <- paste0(paged_url, "&", filter_query)
      }
      
      message("Fetching offset: ", offset)
      tryCatch(
        fromJSON(paged_url),
        error = function(e) {
          message("Error at offset ", offset, ": ", e$message)
          NULL
        }
      )
    })
    
    return(all_data)
  }
  
  # CDC Socrata endpoint
  base_url <- "https://data.cdc.gov/resource/mpgq-jmmr.json"
  #base_url <- "https://data.cdc.gov/resource/ua7e-t2fy.json" #Friday data
  ed_url <- "https://data.cdc.gov/resource/rdmq-nq56.json"
  
  # Define the filter query for the ED data
  filter_query <- "$where=county='All'"
  
  # Download all data (adjust `max_rows` if needed, as of 7/15 about 9,500 rows used)
  health_data_raw <- fetch_socrata_data(base_url, limit = 1000, max_rows = 30000)
  ed_data_raw <- fetch_socrata_data(ed_url, limit = 1000, max_rows = 80000, filter_query = filter_query)
  
  # Now filter it
  health_data <- health_data_raw %>%
    mutate(weekendingdate = as.Date(weekendingdate)) %>%
    filter(weekendingdate >= as.Date("2022-02-01"))
  
  
  #remove  VI, AS, and regions as they are not included for FluSight, keep only necessary vars and add epiweek and epiyear 
  recent_data = health_data %>% 
    dplyr::filter(!jurisdiction %in% c("VI", "AS", "GU", "MP")) %>%
    filter(!str_detect(jurisdiction, "Region")) %>%
    dplyr::select(jurisdiction, weekendingdate, totalconfflunewadm) %>% 
    dplyr::rename("value" = "totalconfflunewadm", "date"="weekendingdate", "state"="jurisdiction") %>% 
    dplyr::mutate(date = as.Date(date), 
                  value = as.numeric(value),
                  state = gsub("USA", "US", state))
                  #target="wk inc flu hosp")
  
  #filter
  ed_data <- ed_data_raw %>% 
    dplyr::select(week_end, geography, percent_visits_influenza) %>% 
    dplyr::rename("date" = "week_end", "state" = "geography", "value" = "percent_visits_influenza") %>% 
    dplyr::mutate(date = as.Date(date), 
                  value = as.numeric(value)/100,
                  state = gsub("United States", "US", state))
                  #target="wk inc flu prop ed visits")
    
  
  ed_full_data = dplyr::left_join(ed_data, locations, by = join_by("state" == "location_name")) %>% 
    rename(location_name = state, state = abbreviation) %>% select(-state, -population)
  
  #bind state population data
  full_data = dplyr::left_join(recent_data, locations, by = join_by("state" == "abbreviation"))
  
  
  #calculates weekly rate 
  final_dat = full_data %>% 
    dplyr::mutate(weekly_rate = (value*100000)/population) %>% 
    select(date, location, location_name, value, weekly_rate)
  
  #return(final_dat, ed_full_data)
    ed_full_data = dplyr::left_join(ed_data, locations, by = join_by("state" == "location_name")) %>% 
    rename(location_name = state, state = abbreviation) %>% select(-state, -population)
  
  #bind state population data
  full_data = dplyr::left_join(recent_data, locations, by = join_by("state" == "abbreviation"))
  
  
  #calculates weekly rate 
  final_dat = full_data %>% 
    dplyr::mutate(weekly_rate = (value*100000)/population) %>% 
    select(date, location, location_name, value, weekly_rate)
  
  #return(final_dat, ed_full_data)
  return(list(final_dat = final_dat, ed_full_data = ed_full_data))
  
}




score_func_param_checker <- function (quantiles, values, actual_value, id = "") 
{
    id_str = paste0(id, ": ")
    if (length(actual_value) > 1) {
        if (length(actual_value) != length(values)) {
            stop(paste0(id_str, "actual_value must be a scalar or the same length", 
                " as values"))
        }
        actual_value = unique(actual_value)
    }
    if (length(actual_value) != 1) {
        stop(paste0(id_str, "actual_value must have exactly 1 unique value"))
    }
    if (length(quantiles) != length(values)) {
        stop(paste0(id_str, "quantiles and values must be of the same length"))
    }
    if (anyDuplicated(quantiles)) {
        stop(paste0(id_str, "quantiles must be unique."))
    }
}



find_quantile_match <- function (quantiles, val_to_match, tol = 1e-08) 
{
    return(abs(quantiles - val_to_match) < tol & !is.na(quantiles))
}



weighted_interval_score <- function (quantile, value, actual_value) 
{
    score_func_param_checker(quantile, value, actual_value, "weighted_interval_score")
    if (all(is.na(actual_value))) 
        return(NA)
    actual_value <- actual_value[[1L]]
    value <- value[!is.na(quantile)]
    quantile <- quantile[!is.na(quantile)]
    med <- value[find_quantile_match(quantile, 0.5)]
    if (length(med) > 1L) 
        return(NA)
    wis <- 2 * mean(pmax(quantile * (actual_value - value), (1 - 
        quantile) * (value - actual_value), na.rm = TRUE))
    return(wis)
}
