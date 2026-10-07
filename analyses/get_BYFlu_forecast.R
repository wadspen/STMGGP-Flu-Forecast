# setwd(here::here())
source("./fluview_functions.R")
source("./clean_data_functions.R")
library(USpopcenters)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(usmap)
library(cmdstanr)
library(maps)
library(patchwork)

team_name <- "BYFlu-STMGGP"

warm <- 1000
draw <- 1000
mod <- cmdstan_model("../stan_models/big_flu_pois_nug.stan")

#data comes in on Wednesday but represents the confirmed hospitalizations from the 
#last week ending on the previous Saturday
#this line reads updated data straight from cdcepi
dat <- read.csv("https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/refs/heads/main/target-data/target-hospital-admissions.csv")

#in case the data have not been updated on the cdcepi website, the following lines of code will get the same data
# results <- fetch_flu()
# dat <- results$final_dat  %>% arrange(desc(date)) 


#clean the data to be formatted properly
#the data are formatted to be by season and season week rather than year. For example the 2024-25 flu season
#begins in August of 2024 and ends in July 2025. It is the 2024 season. Week 1 is the first week of August, and
#week 52 is the last week of the following July. However, if for some reason you don't want to use all data,
#let first_week = w_f for w the week you want as week 1, and final_week = w_f as the final week of data 
data <- clean_data(dat, keep_week = 52, first_week = 1, final_week = 52, seasons = 2022:2026,
                   include_us = TRUE)
cl_dat <- data[[1]]

#leave the US out in the model. To get US forecasts, the predictions for all other 52 regions are summed
#fin_dat is the cleaned data from which a model may be fit. Date is missing where there isn't data observed
#yet and where we want to forecast. These times/locations are still included because the GP that is fit is
#a latent variable with balanced data
fin_dat <- cl_dat |> 
  filter(location_name != "US")
forc_grid <- data[[2]]
last_week_obs <- data[[3]]

pop_dat <- cl_dat %>% 
  dplyr::select(location_name, pop) |> 
  unique()

#everything below here is strictly for fitting in stan and prepping the forecast submission
ref_date <- as.Date(max(fin_dat$date, na.rm = TRUE)) + 7

#prepare data for modeling in stan
stan_data <- get_stan_data(fin_dat)

fit <- mod$sample(data =  stan_data,
                  chains = 1,
                  parallel_chains = 1,
                  iter_warmup = warm,
                  iter_sampling = draw,
                  refresh = 10,
                  show_messages = TRUE)

forcs <- get_forecasts(fit = fit, dat = dat, forc_grid = forc_grid, ref_seas = max(fin_dat$season), 
                       ref_date = ref_date, keep_week = last_week_obs, pop_dat = pop_dat)
forecasts <- forcs[[1]]
forc_grid_preds <- forcs[[2]]

locations <- dat |> 
  dplyr::select(location, location_name) |> 
  unique() |> 
  arrange(location)



rec_file <- paste0("./", team_name, "/", ref_date, "-", team_name, ".csv")
write.csv(forecasts, rec_file, row.names = FALSE)

draws_file <- paste0("./forecast_draws/", ref_date, "-", team_name, ".csv")
write.csv(forc_grid_preds, draws_file, row.names = FALSE)


#this will save the forecast submission to the forked repo. Be careful if you save directly to avoid
#saving over a previously pushed submission. Doing so will lead to problems submitting forecasts to github

submission_file <- paste0("../../FluSight-forecast-hub/model-output/", 
                          team_name, "/", ref_date, "-", team_name, ".csv")
write.csv(forecasts, submission_file, row.names = FALSE)



image_file <- paste0("./forecast_plots/", ref_date, "-", team_name, ".pdf")
pdf(
  file = image_file,
  width = 16,
  height = 8,
  onefile = TRUE
)

for (l in 1:nrow(locations)) {
  loc <- locations$location[l]
  loc_name <- locations$location_name[l]
  q_plot <- forecasts |>
    filter(target == "wk inc flu hosp", output_type == "quantile") |> 
    mutate(prop = paste0("p", output_type_id)) |> 
    dplyr::select(-output_type_id) |> 
    pivot_wider(names_from = "prop") |> 
    filter(location == loc) |> 
    ggplot() +
    geom_ribbon(aes(x = target_end_date, ymin = p0.01, ymax = p0.99), fill = "darkred") +
    geom_ribbon(aes(x = target_end_date, ymin = p0.025, ymax = p0.975), fill = "red") +
    geom_ribbon(aes(x = target_end_date, ymin = p0.25, ymax = p0.75), fill = "pink") +
    geom_line(aes(x = target_end_date, y = p0.5), colour = "purple", size = 2) +
    geom_point(data = cl_dat |> 
              filter(!is.na(date), 
                location == loc, season == 2026) |> 
              mutate(date = as.Date(date)),
              aes(x = date, y = rate), size = 3) +
    ggtitle(loc_name) +
    xlab("Date") +
    ylab("Hospitalizations") +
    theme_bw() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 22),
          title = element_text(size = 32))
  
  s_plot <- forecasts |>
    filter(target == "wk inc flu hosp", output_type == "sample", location == loc) |> 
    group_by(reference_date, target, horizon, target_end_date, location) |> 
    summarise(med = median(value),
              upp = quantile(value, probs = 0.975),
              low = quantile(value, probs = 0.025)) |> 
    mutate(target_end_date = as.Date(target_end_date)) |> 
    ggplot() +
    geom_ribbon(aes(x = target_end_date, ymin = low, ymax = upp), fill = "pink") +
    geom_line(aes(x = target_end_date, y = med), colour = "purple", size = 2) +
    geom_point(data = cl_dat |> 
                filter(!is.na(date), 
                  location == loc, season == 2026) |> 
                mutate(date = as.Date(date)),
                aes(x = date, y = rate), size = 3) +
    theme_bw()

  t_plot <- forecasts |> 
  filter(target == "wk flu hosp rate change", location == loc) |> 
  mutate(Trend = output_type_id) |> 
  ggplot() +
  geom_col(aes(x = target_end_date, y = value, fill = Trend)) +
  theme_bw() +
  xlab("Date") +
    ylab("Hospitalizations") +
    theme_bw() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 22),
          title = element_text(size = 22),
          legend.text = element_text(size = 15))
  
  p_plot <- forecasts |> 
  filter(target == "peak week inc flu hosp", location == loc) |> 
  # group_by(output_type_id) |> 
  # summarise(sum_probs = mean(value)) |> 
  mutate(date = as.Date(output_type_id)) |> 
  ggplot() +
  geom_col(aes(x = date, y = value)) +
  ggtitle("Peak Week Probabilities") +
  xlab("Date") +
  ylab("") +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        title = element_text(size = 27))


  print((q_plot + s_plot) / (p_plot + t_plot))
}

dev.off()


library(hubValidations)
hubValidations::validate_model_metadata(
    hub_path="../../FluSight-forecast-hub/",
    file_path="BYFlu-STMGGP.yml")

Sys.setenv(TZDIR = path.expand("~/zoneinfo"))
validations <- validate_submission(
    hub_path = "../../FluSight-forecast-hub/",
    file_path = "BYFlu-STMGGP/2026-10-10-BYFlu-STMGGP.csv")

check_for_errors(validations)
















  