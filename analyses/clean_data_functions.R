library(USpopcenters)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(usmap)
library(maps)




clean_data <- function(dat, seasons = 2022:2026, keep_week, first_week = 9, final_week = 44, k_states = NULL,
                       include_us = FALSE) {
  
                  max_seas <- max(seasons)
                  keep_states <- if (is.null(k_states) || length(k_states) < 2) {
                                          unique(dat$location_name)
                                        } else {
                                          k_states
                                        }
                  if (include_us == FALSE) {keep_states <- keep_states[keep_states != "US"]}
                  dat <- dat |> 
                    mutate(year = year(date),
                    week = week(date))

  
  
                  cl_dat <- dat %>% 
                  mutate(season = ifelse(month(date) <= 7 | 
                                          month(date) == 8 & mday(date) < 7, 
                                        year(date) - 1, year(date)),
                        place = 1) %>% 
                  arrange(desc(season), date) %>% 
                  group_by(season, location) %>% 
                  mutate(season_week = cumsum(place)) %>%
                  select(-place) %>% 
                  dplyr::select(date, location, location_name, season, season_week, 
                                value, weekly_rate) 

  

                states_dat <- cl_dat %>% 
                  arrange(location, season, season_week) %>% 
                  filter(season >= min(seasons)) %>% 
                  mutate(pop = round(value/weekly_rate, 1)) %>% 
                  mutate(seas_week_dum = ifelse(location_name %in% c("Minnesota", 
                                                                    "Massachusetts",
                                                                    "West Virginia") &
                                                  season == 2024 &
                                                  season_week <= 9, NA,
                                                ifelse(location_name %in% c("Minnesota",
                                                                            "Massachusetts") &
                                                        season == 2023 &
                                                        season_week >= 44, NA, 
                                                      season_week))) %>% 
                  mutate(rate = ifelse(is.na(value), 0,
                                      value)) %>% 
                  filter(location_name %in% keep_states) %>% 
                  arrange(location_name, season, season_week) %>% 
                  filter(season_week <= final_week & season_week >= first_week) %>% 
                  mutate(season_week = season_week - (first_week - 1)) %>% 
                  filter(location_name %in% keep_states) |> 
                  filter(season < max_seas | (season == max_seas & season_week <= keep_week))
  
  
                forc_grid <- expand.grid(season_week = unique(states_dat$season_week), season = seasons, 
                                        location_name = unique(states_dat$location_name)) %>% 
                  arrange(location_name, season, season_week)
  
                forc_grid <- rbind(forc_grid |> 
                                    filter(location_name != "US"),
                                   forc_grid |> 
                                      filter(location_name == "US"))


                fin_dat <- forc_grid %>% 
                  left_join(states_dat %>% 
                              ungroup() %>% 
                              dplyr::select(date, location_name, season, season_week, rate,
                                            seas_week_dum),
                            by = c("location_name", "season", "season_week")) %>%
                  mutate(rate = ifelse(is.na(rate), 0, rate))


                
                pop_data <- dat %>% 
                    mutate(pop = round(value/weekly_rate, 1)) %>% 
                    filter(!is.na(pop)) %>% 
                    dplyr::select(location_name, pop) %>% 
                    unique() 
                  rownames(pop_data) <- 1:nrow(pop_data)
  

                fin_dat <- fin_dat %>% 
                  left_join(pop_data, by = "location_name") |> 
                  left_join(dat |> 
                              dplyr::select(location_name, location) |> 
                              unique(),
                          by = "location_name")
  
                last_date <- fin_dat |> 
                  mutate(date = as.Date(date)) |> 
                  slice_max(date, with_ties = TRUE) |>
                  dplyr::select(date, season_week) |> 
                  unique()
  
                return(list(fin_dat, forc_grid, last_date$season_week))
}



get_stan_data <- function(dat, k_states = NULL,
                          include_us = FALSE) {
                        
                      keep_states <- if (is.null(k_states) || length(k_states) < 2) {
                                                              unique(dat$location_name)
                                                            } else {
                                                              k_states
                                                            }
                                      if (include_us == FALSE) {keep_states <- keep_states[keep_states != "US"]}
  
  
                        data("state2020")
                        pop_cent <- state2020
                        colnames(pop_cent) <- c("location", "location_name", "pop", "lat", "long")


                        inds <- dat %>%
                          ungroup() %>% 
                          mutate(row_index = row_number()) %>%
                          group_by(season, location_name) %>%
                          summarise(
                            first_index = row_index[which.min(seas_week_dum)],
                            last_index  = row_index[which.max(seas_week_dum)],
                            .groups = "drop"
                          ) %>% 
                          arrange(location_name, season) %>% 
                          dplyr::select(first_index, last_index)


                        time_mat <- dist(unique(dat$season_week)) %>% 
                          as.matrix()

                        seas_mat <- dist(unique(dat$season)) %>% 
                          as.matrix()

                        keep_pops <- pop_cent %>% 
                          filter(!(location_name %in% c("US"))) %>%
                          filter(location_name %in% keep_states) %>% 
                          arrange(location_name)

                        locs_mat <- dist(keep_pops[,c("lat", "long")]) %>% 
                          as.matrix()

                        stan_list <- list(
                          N = nrow(dat),
                          nweeks = length(unique(dat$season_week)),
                          nseas = length(unique(dat$season)),
                          nlocs = length(unique(dat$location_name)),
                          weeks = dat$season_week, #unique(dat$season_week),
                          seas = as.numeric(unique(dat$season)),
                          locs = as.numeric(factor(unique(dat$location_name))),
                          tdist = time_mat,
                          sdist = seas_mat,
                          ldist = locs_mat,
                          rate = dat$rate,
                          pop = dat$pop,
                          first_ind = inds$first_index,
                          last_ind = inds$last_index,
                          sig_alpha = 1
                        )
                      
                      return(stan_list)
}








get_forecasts <- function(fit, dat, forc_grid, ref_seas = 2026, ref_date, keep_week = 52,
                          probs = c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99),
                          pop_dat) {
                    
                    locations <- dat |> 
                      dplyr::select(location_name, location) |> 
                      unique()
                    base_date <- ref_date - 7
                    week_dates <- data.frame(season_week = 1:52, 
                      prev_seas_last_date = base_date - 7*keep_week) |> 
                      mutate(date = prev_seas_last_date + 7*season_week)

  
                    draws <- fit$draws(format = "df")


                    fs <- draws %>%
                      dplyr::select(contains("pred_f["))

                    
                    
                    # pred_quants <- apply(fs, MARGIN = 2, FUN = quantile, probs = probs) %>% t()
  
                    tfs <- fs |> 
                      t() |> 
                      as.data.frame()
  
                    rownames(tfs) <- 1:nrow(tfs)
                    colnames(tfs) <- paste0("draw", 1:ncol(tfs))
                    
                    us_fs <- cbind(forc_grid |> 
                                      filter(location_name != "US"), tfs) |> 
                      pivot_longer(contains("draw"), names_to = "rep", values_to = "draw") |> 
                      group_by(season_week, season, rep) |> 
                      summarise(draw = sum(draw)) |> 
                      mutate(rep = as.numeric(str_replace(rep, "draw", ""))) |> 
                      arrange(season, season_week, rep) |> 
                      mutate(rep = paste0("draw", rep)) |> 
                      pivot_wider(names_from = "rep", values_from = "draw") |> 
                      ungroup() |> 
                      dplyr::select(contains("draw")) |> 
                      t() 
  
                    # us_pred_quants <- apply(us_fs, MARGIN = 2, FUN = quantile, probs = probs) %>% t()
                    afs <- cbind(fs, us_fs)
                    pred_quants <- apply(afs, MARGIN = 2, FUN = quantile, probs = probs) %>% t()
                    forc_grid_w <- cbind(forc_grid, pred_quants)
                    rownames(forc_grid_w) <- 1:nrow(forc_grid_w)
  
                    
                    forecasts <- forc_grid_w %>% 
                      filter(season == ref_seas,
                            season_week >= keep_week &
                              season_week <= keep_week + 4) %>% 
                      pivot_longer(contains("%"), names_to = "output_type_id", 
                                  values_to = "quantile") %>% 
                      mutate(output_type_id = 
                              as.character(as.numeric(str_replace(output_type_id, "%", ""))/100)) %>% 
                      left_join(pop_dat, by = "location_name") %>% 
                      left_join(dat %>% 
                                  dplyr::select(location_name, location) %>% 
                                  unique(), by = "location_name") %>% 
                      mutate(value = quantile) %>% 
                      mutate(horizon = season_week - keep_week - 1,
                            target = "wk inc flu hosp",
                            output_type = "quantile") %>% 
                      mutate(reference_date = ref_date,
                            target_end_date = as.Date(ref_date) + horizon * 7) %>% 
                      dplyr::select(reference_date, target, horizon, target_end_date,
                                    location, output_type, output_type_id, value)
  
  
  
  
                    date_hoz <- forecasts |> 
                      dplyr::select(target_end_date, horizon) |> 
                      unique()
  
                    seas_draws <- cbind(forc_grid, afs |> t()) |> 
                      filter(season == 2026)|> 
                      pivot_longer(contains("draw"), names_to = "rep", values_to = "draw")

                    rate_trend_dist <- seas_draws |> 
                      filter(between(season_week, keep_week + 1, keep_week + 4)) |> 
                      left_join(dat |> 
                                  filter(date == base_date) |> 
                                  dplyr::select(location_name, value) |> 
                                  unique(),
                              by = "location_name") |> 
                      left_join(pop_dat, by = "location_name") |> 
                      mutate(horizon = season_week - (keep_week + 1)) |> #filter(location_name == "US")
                      rowwise() |> 
                      mutate(count_change = (draw - value),
                             rate_change = (draw - value)/pop) |> 
                      dplyr::select(location_name, horizon, count_change, rate_change) |> 
                      mutate(
                            output_type_id = case_when(
                              # Stable: horizon-specific rate threshold OR fewer than 10 admissions
                              horizon == 0 & (abs(rate_change) < 0.3 | abs(count_change) < 10) ~ "stable",
                              horizon == 1 & (abs(rate_change) < 0.5 | abs(count_change) < 10) ~ "stable",
                              horizon == 2 & (abs(rate_change) < 0.7 | abs(count_change) < 10) ~ "stable",
                              horizon == 3 & (abs(rate_change) < 1.0 | abs(count_change) < 10) ~ "stable",

                              # Large increase
                              horizon == 0 & rate_change >= 1.7 ~ "large_increase",
                              horizon == 1 & rate_change >= 3.0 ~ "large_increase",
                              horizon == 2 & rate_change >= 4.0 ~ "large_increase",
                              horizon == 3 & rate_change >= 5.0 ~ "large_increase",

                              # Increase
                              horizon == 0 & rate_change > 0 ~ "increase",
                              horizon == 1 & rate_change > 0 ~ "increase",
                              horizon == 2 & rate_change > 0 ~ "increase",
                              horizon == 3 & rate_change > 0 ~ "increase",

                              # Large decrease
                              horizon == 0 & rate_change <= -1.7 ~ "large_decrease",
                              horizon == 1 & rate_change <= -3.0 ~ "large_decrease",
                              horizon == 2 & rate_change <= -4.0 ~ "large_decrease",
                              horizon == 3 & rate_change <= -5.0 ~ "large_decrease",

                              # Decrease
                              horizon == 0 & rate_change < 0 ~ "decrease",
                              horizon == 1 & rate_change < 0 ~ "decrease",
                              horizon == 2 & rate_change < 0 ~ "decrease",
                              horizon == 3 & rate_change < 0 ~ "decrease",

                              TRUE ~ NA_character_
                            )
                          )
                    
                    rate_trend_pmf <- rate_trend_dist %>%
                      count(location_name, horizon, output_type_id) %>%
                      group_by(location_name, horizon) %>%
                      tidyr::complete(
                        output_type_id = c(
                          "large_increase",
                          "increase",
                          "stable",
                          "decrease",
                          "large_decrease"
                        ),
                        fill = list(n = 0)
                      ) %>%
                      mutate(value = n / sum(n)) %>%
                      ungroup() |> 
                      mutate(output_type = "pmf",
                             reference_date = ref_date,
                             target = "wk flu hosp rate change") |> 
                      left_join(date_hoz, by = "horizon") |> 
                      left_join(locations, by = "location_name") |> 
                      dplyr::select(reference_date, target, horizon, target_end_date,
                                    location, output_type, output_type_id, value)
  
                    max_hosp_week <- seas_draws |> 
                      group_by(location_name, rep) |>
                      slice_max(draw, n = 1, with_ties = FALSE) |> 
                      ungroup() |> 
                      arrange(location_name, season_week) |> 
                      group_by(location_name, season_week) |> 
                      summarise(n = n(), .groups = "drop_last") |> 
                      mutate(prop = n / sum(n)) |> 
                      ungroup()
  
                    max_hosp_week <- seas_draws |>
                      group_by(location_name, rep) |>
                      slice_max(draw, n = 1, with_ties = FALSE) |>
                      ungroup() |>
                      filter(season_week %in% 10:43) |> 
                      count(location_name, season_week, name = "n") |>
                      group_by(location_name) |>
                      complete(
                        season_week = 1:43,
                        fill = list(n = 0)
                      ) |>
                      mutate(value = n / sum(n)) |>
                      ungroup() |>
                      arrange(location_name, season_week) |> 
                      filter(season_week %in% 10:43) |> 
                      left_join(week_dates, by = "season_week") |> 
                      left_join(locations, by = "location_name") |> 
                      mutate(output_type_id = as.character(date)) |> 
                      mutate(target = "peak week inc flu hosp",
                             output_type = "pmf", horizon = NA,
                             reference_date = ref_date, target_end_date = NA) |> 
                      dplyr::select(reference_date, target, horizon, target_end_date,
                                    location, output_type, output_type_id, value)
  
                    replace <- ifelse(nrow(draws) < 100, TRUE, FALSE)
                    samp_draws <- sample(1:nrow(draws), 100, replace = replace)
                    samp_forecasts <- seas_draws |> 
                      mutate(rep = as.numeric(str_replace(rep, "draw", ""))) |> 
                      filter(rep %in% samp_draws) |> 
                      arrange(location_name, season_week, rep) |> 
                      mutate(ind = 1) |> 
                      group_by(season_week, season, location_name) |> 
                      mutate(label = cumsum(ind)) |> 
                      mutate(output_type_id = paste0("SA_", label)) |> 
                      left_join(locations, by = "location_name") |> 
                      left_join(week_dates, by = "season_week") |> 
                      mutate(target_end_date = date, value = draw) |> 
                      left_join(date_hoz, by = "target_end_date") |> 
                      filter(!is.na(horizon)) |> 
                      mutate(reference_date = ref_date, target = "wk inc flu hosp", 
                             output_type = "sample") |> 
                      ungroup() |> 
                      dplyr::select(reference_date, target, horizon, target_end_date,
                                    location, output_type, output_type_id, value)
  
                  peak_quantiles <- seas_draws |> 
                    filter(season_week %in% 10:43) |> 
                    group_by(location_name, rep) |> 
                    summarise(value = max(draw)) |> 
                    group_by(location_name) |> 
                    reframe(value = quantile(value, probs = probs))
  
                  peak_quantiles <- cbind(peak_quantiles, 
                    output_type_id = rep(probs, 
                                         length(unique(peak_quantiles$location_name)))) |> 
                    left_join(locations, by = "location_name") |> 
                    mutate(output_type = "quantile", horizon = NA, target_end_date = NA,
                           reference_date = ref_date, target = "peak inc flu hosp") |> 
                    dplyr::select(reference_date, target, horizon, target_end_date, location,
                                  output_type, output_type_id, value)
  
                  all_targets <- rbind(forecasts, samp_forecasts, peak_quantiles, max_hosp_week, rate_trend_pmf)
  
  return(list(all_targets, forc_grid_w))
}
