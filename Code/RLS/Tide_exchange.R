# Script for tide exhange calculations
# Em Lim
# Sept 4, 2023

# Load packages
library(tidyverse)
library(lubridate)

# Load data
# Data specifications:
# 5 min intervals
# 2 weeks from first day of blitz
# meters
# 24 hour time
# 2021 April 26 start
# added May 20, 2021 for the two Dixon sites
# 2022 April 25 start
# 2023 May 8 start

tide <- read_csv("Data/RLS/tides_1_min.csv") %>%
  mutate(date_time = ymd_hms(paste(date, time)))


# build an empty dataframe
# Do this each time you run the for loop!!!
tide_new <- data.frame()

# then write the loop
for (x in 1:nrow(rls_survey_info)) {
  
  survey_start <- ymd_hms(rls_survey_info$date_time[x:x])
  survey_end <- survey_start + hours(1)
  
  output = tide %>%
    filter(between(date_time, survey_start, survey_end)) %>%
    mutate(rate = 100 * (tide_m - lag(tide_m))/lag(tide_m)) %>%
    slice(-1) %>%
    summarise(avg_rate = mean(rate)) %>%
    mutate(survey_id = rls_survey_info$survey_id[x:x])
  
  tide_new = rbind(tide_new, output)
  
}

# Now I have the average rate of change of tide height for each survey!!!