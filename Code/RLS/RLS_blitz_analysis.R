# Code to look at all the spatial data together
# June 2, 2023
# Em Lim

# Load packages and functions ----
library(tidyverse)
library(visreg)
library(ggplot2)
library(PNWColors)
library(ggeffects)
library(lubridate)
library(sf)
library(lmerTest)

# Source pretty functions
source("Code/theme_black.R")
source("Code/Functions.R") # Length to weight function here!

# Set default plotting
theme_set(theme_bw())

# Load data + manipulate RLS data ------
# RLS data from the website!

# Just the pelagic and cryptic fish
fish <- read_csv("Data/RLS/RLS_data/reef_fish_abundance_and_biomass.csv",
                 show_col_types = FALSE) %>%
  rbind(read_csv("Data/RLS/RLS_data/cryptobenthic_fish_abundance.csv",
                 show_col_types = FALSE)) %>%
  as.data.frame() %>%
  mutate(survey_date = ymd(survey_date),
         year = as.factor(year(survey_date)),
         site_code = ifelse(site_name == "Swiss Boy", "BMSC24", site_code)) %>%
  rename(site_ID = site_code) %>%
  filter(month(survey_date) == 4 | month(survey_date) == 5) # Just the RLS blitz data for now

# Just the mobile inverts
invert <- read_csv("Data/RLS/RLS_data/mobile_macroinvertebrate_abundance.csv",
                   show_col_types = FALSE) %>%
  as.data.frame() %>%
  mutate(survey_date = ymd(survey_date),
         year = as.factor(year(survey_date)),
         site_code = ifelse(site_name == "Swiss Boy", "BMSC24", site_code),
         species_name = case_when(
           species_name == "Montereina nobilis" ~ "Peltodoris nobilis",
           species_name == "Parastichopus californicus" ~ "Apostichopus californicus",
           species_name == "Berthella californica" ~ "Berthella chacei",
           species_name == "Henricia leviuscula" ~ "Henricia spp.",
                        TRUE ~ as.character(species_name))) %>%
  rename(site_ID = site_code) %>%
  filter(month(survey_date) == 4 | month(survey_date) == 5) # Just the RLS blitz data for now

# Join all rls data together
rls <- rbind(fish, invert)

# extract just one row per survey to join with the pee data and tide data
rls_survey_info <- rls %>%
  select(site_ID, year, depth, survey_id, survey_date, hour) %>%
  rename(survey_depth = depth) %>%
  unique() %>%
  mutate(date_time = ymd_hms(paste(survey_date, hour)))

# RLS Biomass calculations -------

# A note on length to weight relationships
# The overall relationship is W = a*L^b
# log form of the above formula is log(W) = log(a) + b*log(L)
# The Siegel paper presents the log transformed a value
# Eg if the Siegel paper says a = -4.60517 that corresponds to log(0.01)
# So some papers will present a already log transformed (a negative value) and some will present a as an untransformed value (positive decimal)
# log(W) = log(a) + b*log(L) is the same as log10(W) = log10(a) + b*log10(L) so don't worry about the base as long as it's the same across the formula

# Use the a and b parameters from FishBase
# W = a*L^b
# log form of the above formula is log(W) = log(a) + b*log(L)
# W = exp(log(a) + b*log(L))
# Cite Froese R, Thorson JT, Reyes Jr RB. A Bayesian approach for estimating length‐weight relationships in fishes. Journal of Applied Ichthyology. 2014 Feb;30(1):78-85.

# Use WL relationships from fishbase to estimate weight of fish in RLS surveys
fishes <- fish %>%
  mutate(size_class = if_else(size_class == 187.5, 87.5, size_class)) %>% # shrink that one wolf eel
  length_to_weight() %>% # Use nice length to weight function!
  filter(species_name != "Bolinopsis infundibulum") %>%
  filter(species_name != "Pleuronichthys coenosus") %>%
  filter(species_name != "Pleurobrachia bachei") %>%
  filter(species_name != "Polyorchis penicillatus") %>% # Filter inverts
  filter(species_name != "Phoca vitulina") %>% # Remove seal
  filter(species_name != "Actinopterygii spp.") %>% # Remove unidentif fish
  filter(species_name != "Myoxocephalus aenaeus") # Remove east coast fish

# That one huge wolf eel can't be right
# Fishbase: max size = 240 cm, max weight = 18.4 kg
# Ours is apparently 187 cm and 63 kg
# That formula must be off
# I also have size data for inverts "Haliotis kamtschatkana", "Crassadoma gigantea", "Pycnopodia helianthoides", "Polyorchis penicillatus", "Bolinopsis infundibulum", Pleurobrachia bachei, Pleuronichthys coenosus

# We also measured sea cucumbers a few times???? Where did that data go?
# We also have urchin size data from this year

# Some people counted M2 fishes on M1, but not always so I don't actually want goby and sculpin counts from M1
# Figure out what to do here

# Load + manipulate nh4+ data ----

# RLS nh4 from 2021
# This is the spring blitz, the june samples, and the July samples
# The may samples have estimated matrix effects, not great
# The June and July samples had a own matrix spike from each site but were still compared to DI
bottles2021 <- read_csv("Output/Output_data/RLS_nh4_2021.csv") 

# RLS from 2022
# These had a own matrix spike from each site but were still compared to DI
# They came out quite negative maybe because of SFU DI, and maybe because of temperature difference between the standards and samples
# I took the lowest nh4 reading and added it to everything else to "set" that sample to 0 and bump everything up
# Maybe an underestimation
bottles2022 <- read_csv("Output/Output_data/RLS_nh4_2022.csv")

# RLS from 2023
# Did the full "proper" Taylor protocol with standard bottles + BF from each site
bottles2023 <- read_csv("Output/Output_data/RLS_nh4_2023.csv")

# combine these three years into one!
rls_nh4_3years <- rbind(bottles2021, bottles2022, bottles2023) %>%
  mutate(year = as.factor(year)) %>%
  filter(month == "May") 

# just keep the surveys where I have an nh4 AND rls survey on the same transect
rls_nh4 <- rls_nh4_3years %>% 
  left_join(rls_survey_info, by = c("site_ID", "year")) %>%
  depth_function() # only keep the RLS survey from the transect where the pee is from


# Load + manipulate tide exchange data ----

# Downloaded tide height data (m) from http://tbone.biol.sc.edu/tide/tideshow.cgi in 1 min intervals over the two week spanning each RLS spring blitz, in 24 hour time
# 2021 April 26 start
# added May 20, 2021 for the two Dixon sites
# 2022 April 25 start
# 2023 May 8 start

tide <- read_csv("Data/RLS/tides_1_min.csv") %>%
  mutate(date_time = ymd_hms(paste(date, time)))

# build an empty dataframe
# Do this each time you run the for loop!!!
tide_exchange <- data.frame()

# then write the loop
for (x in 1:nrow(rls_survey_info)) {
  
  survey_start <- ymd_hms(rls_survey_info$date_time[x:x])
  survey_end <- survey_start + hours(1)
  
  output = tide %>%
    filter(between(date_time, survey_start, survey_end)) %>%
    mutate(rate = 100 * (tide_m - lag(tide_m))/lag(tide_m)) %>%
    slice(-1) %>%
    summarise(avg_exchange_rate = mean(rate)) %>%
    mutate(survey_id = rls_survey_info$survey_id[x:x])
  
  tide_exchange = rbind(tide_exchange, output)
}
# Now I have the average rate of change of tide height for each survey!!!


# Site level averaging -----

# Average nh4+ for each survey
nh4_avg <- rls_nh4 %>%
group_by(survey_id) %>%
  mutate(nh4_avg = mean(nh4_conc),
         year = as.factor(year)) %>%
  ungroup() %>%
  select(c(site, site_ID, survey_id, date_time, nh4_avg)) %>%
  unique()

# Calculate total fish biomass per survey
fish_biomass <- fishes %>%
  group_by(survey_id) %>%
  summarize(weight_sum = sum(weight_size_class_sum)) %>%
  left_join(rls %>% # also add overall site species richness
      group_by(survey_id) %>%
      summarize(species_richness = n_distinct(species_name))%>%
      arrange(desc(species_richness))) 
# One row per transect
  # Per site, per depth, per year
  # Same as grouping by survey_ID

# Join nh4 + fish biomass + tide exchange data
rls_final <- nh4_avg %>%
  left_join(fish_biomass, by = c("survey_id")) %>%
  left_join(tide_exchange, by = "survey_id") # add in the tide exchange data

# One row per survey
# Each survey has an average nh4 concentration, total fish biomass, and tide exchange


# Data exploration ------

# Data checks
goby <- fish %>%
  filter(species_name == "Rhinogobiops nicholsii") %>%
  filter(site_name == "Goby Town") %>%
  filter(depth == "5.7") %>%
  select(site_name, depth, method, block, species_name, total, size_class)

# New nudis?
nudi <- invert %>%
  filter(class == "Gastropoda") %>%
  filter(family != "Acmaeidae") %>%
  filter(family != "Calliostomatidae") %>%
  filter(family != "Epitoniidae" ) %>%
  filter(family != "Fissurellidae" ) %>%
  filter(family != "Haliotidae" ) %>%
  filter(family != "Lottiidae" ) %>%
  filter(family != "Muricidae" ) %>%
  filter(family != "Naticidae" ) %>%
  filter(family != "Tegulidae" ) %>%
  filter(family != "Turbinidae") %>%
  count(species_name, year)

# New RLS species?
rls_new <- rls %>%
  count(species_name, year) %>%
  select(-n)

prev <- rls_new %>%
  filter(year != "2023") %>%
  mutate(seen_before = "yes") %>%
  select(-year) %>%
  unique()

new <- rls_new %>%
  filter(year == "2023") %>%
  mutate(seen_this_year = "yes") %>%
  select(-year) %>%
  left_join(prev, by = "species_name") %>%
  replace_na(list(seen_before = "no")) %>%
  filter(seen_before == "no") 

#write_csv(new, "Output/Output_data/new_species.csv")

# missing size data
no_sizes <- fish %>%
  filter(size_class == "0") %>%
  filter(species_name != "Bolinopsis infundibulum") %>%
  filter(species_name != "Pleuronichthys coenosus") %>%
  filter(species_name != "Pleurobrachia bachei") %>%
  filter(species_name != "Polyorchis penicillatus") %>%
  select(-year)

#write_csv(no_sizes, "Output/Output_data/missing_fish_sizes.csv")

# What's the most abundant species?
rls_abundant_species <- rls %>% 
  count(species_name)

# where are the most urchins???
# Can swap for any species
rls_urchins <- rls %>%
  filter(species_name == "Mesocentrotus franciscanus") %>%
  group_by(site_name) %>%
  summarize(urchins = sum(total)) %>%
  arrange(desc(urchins))

# So basically you can manipulate the data to summarize whatever thing you're interested in, and then you can join that summarized data with another df of interest by the common site code/site name :)

# rank each site by the year
rank_2021 <- rls_final %>%
  filter(year(date_time) == "2021") %>%
  group_by(site) %>%
  summarise(nh4_avg2021 = mean(nh4_avg)) %>%
  arrange(desc(nh4_avg2021)) %>% 
  mutate(rank2021 = 1:22,
         grade2021 = 100 - (100*rank2021/22))

rank_2022 <- rls_final %>%
  filter(year(date_time) == "2022") %>%
  group_by(site) %>%
  summarise(nh4_avg2022 = mean(nh4_avg)) %>%
  arrange(desc(nh4_avg2022)) %>% 
  mutate(rank2022 = 1:19,
         grade2022 = 100 - (100*rank2022/19))

rank_2023 <- rls_final %>%
  filter(year(date_time) == "2023") %>%
  group_by(site) %>%
  summarise(nh4_avg2023 = mean(nh4_avg)) %>%
  arrange(desc(nh4_avg2023)) %>% 
  mutate(rank2023 = 1:20,
         grade2023 = 100 - (100*rank2023/20))

rank <- rank_2021 %>%
  left_join(rank_2022, by= "site") %>%
  left_join(rank_2023, by= "site") %>%
  select(grade2021, grade2022, grade2023, site) %>%
  rowwise() %>%
  mutate(avg_grade = mean(c(grade2021, grade2022, grade2023), na.rm = TRUE))


# What were the matrix effects like in 2022?
rls_nh4 %>%
  filter(year == "2022") %>%
  summarise(matrix = mean(matrix))
# 14.6

# how about the end of 2021?
rbind(bottles2021, bottles2022, bottles2023) %>%
  mutate(year = as.factor(year)) %>%
  filter(year == "2021") %>%
  filter(month != "May") %>%
  summarise(matrix = mean(matrix))
# 7.82
# Go back to the 2021 data and set the ME to 7.82


#### code checking only to here -----
# Plot these data??? ----
pal <- pnw_palette("Sailboat", 3)

ggplot(rls_nh4) +
  geom_point(aes(x = reorder(site, -nh4_avg), 
                 nh4_conc, colour = year, fill = year, pch = month),
             size = 3, alpha = 0.75) +
  geom_point(aes(x = reorder(site, -nh4_avg), 
                 nh4_avg),
                 size = 3, colour = "black", pch = 20) +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5)) +
  labs(x = "Site", y = "NH4+ Concentration (umol/L)") +
  scale_colour_manual(values = pal) 

#black
ggplot(rls_nh4) +
  geom_point(aes(y = reorder(site, -nh4_avg), 
                 x = nh4_conc, colour = year, fill = year),
             size = 3, alpha = 0.9) +
  geom_point(aes(y = reorder(site, -nh4_avg), 
                 x = nh4_avg),
             size = 4, colour = "white", pch = 20) +
  labs(y = "Site", x = "NH4+ Concentration (umol/L)") +
  scale_colour_manual(values = pal) + 
  theme_black()

 ggsave("Output/Figures/RLS_nh4_all_years_black.png", device = "png",
        height = 9, width = 16, dpi = 400)


# Plot the ranking of each site year to year
# Plot rank
ggplot(rank) +
  geom_point(aes(reorder(site, -avg_grade), grade2021, colour = "2021"), size = 3) +
  geom_point(aes(reorder(site, -avg_grade), grade2022, colour = "2022"), size = 3) +
  geom_point(aes(reorder(site, -avg_grade), grade2023, colour = "2023"), size = 3) +
  labs(x= "Site", y = "Rank", colour = "Year") +
  scale_colour_manual(values = pal) + 
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5)) 

# ggsave("Output/Figures/rank_all_years.png", device = "png",
#        height = 5, width = 8, dpi = 400)


# Plot biomass vs nh4
ggplot(rls_final, aes(weight_sum, nh4_avg, label = site)) +
  geom_point(aes(colour = site_ID)) +
  geom_smooth(method = lm) +
  geom_text(check_overlap = TRUE, hjust = 1,  size = 3, ) +
  labs(x = "Total fish biomass (g)", y = "Ammonium concentration (umol)") +
  theme_classic()



# Correlation analysis -----
# Make a DF where we only have overlap between sites
b_2021 <- rls_nh4 %>%
  filter(year == "2021") %>%
  transmute(site = site,
            nh4_2021 = nh4_avg) %>%
  unique()

b_2022 <- rls_nh4 %>%
  filter(year == "2022") %>%
  transmute(site = site,
            nh4_2022 = nh4_avg) %>%
  unique()

b_2023 <- rls_nh4 %>%
  filter(year == "2023") %>%
  transmute(site = site,
            nh4_2023 = nh4_avg) %>%
  unique()

# throw them together and drop the NAs
cor_data <- b_2021 %>%
  left_join(b_2022, by = "site") %>%
  left_join(b_2023, by = "site") %>%
  drop_na()

# OK so now do the correlations!

# 2021 vs 2022
cor.test(cor_data$nh4_2021, cor_data$nh4_2022, 
                method = "spearman")
# Yes correlated?

# 2022 vs 2023
cor.test(cor_data$nh4_2022, cor_data$nh4_2023, 
         method = "spearman")
# Yes correlated?

# 2023 vs 2021
cor.test(cor_data$nh4_2023, cor_data$nh4_2021, 
         method = "spearman")
# Yes correlated?

# Stats -------
# Does pee vary by site?
simple_model <- lm(nh4_conc ~ site_ID + year, data = rls_nh4)
summary(simple_model)
visreg(simple_model)

# does pee vary with biomass
model <- lmer(nh4_avg ~ weight_sum + (1|site_ID), rls_final)
summary(model)
visreg(model)

pee <- lm(nh4_conc ~ site * depth, rls_nh4)
anova(pee)
summary(pee)

visreg(pee, "depth", by = "site")



sum_stats_pee <- ggpredict(simple_model, terms = c("site_ID", "period")) %>% 
  #and then we'll just rename one of the columns so it's easier to plot
  rename(site_ID = x,
         nh4_conc = predicted,
         period = group)
#View(sum_stats_crabs)

# Graphing ----

# Dot and whisker?
ggplot() +
  geom_point(data = sum_stats_pee, 
             aes(y = site_ID, x = nh4_conc, colour = period),
             size = 4) +
  geom_errorbar(data = sum_stats_pee, 
                aes(y = site_ID,
                    x = nh4_conc,
                    colour = period,
                    # and you can decide which type of error to show here
                    # we're using 95% CI
                    xmin = conf.low,
                    xmax = conf.high),
                width = 0.2,
                size = 1.2)  +
  geom_point(data = rls_data, aes (y = site_ID, x = nh4_conc, colour = period), alpha = 0.5, height = 0, size = 2) +
  labs(x= "Ammonium concentration (umol/L)", y = "Site") +
  theme_black() +
  theme(legend.position="none") 

#ggsave("Output/Figures/RLS_pee_black.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# Site vs pee
ggplot(rls_data, aes(adj_nh4_conc, site, colour = period)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  labs(x= "Ammonium concentration (umol)", y = "Site") +
  theme_classic()

# black background
ggplot(rls_data, aes(nh4_conc, site, fill = site)) +
  geom_boxplot(colour = "white") +
  theme_black() +
  labs(x= "Ammonium concentration (umol/L)", y = "Site") +
  theme(legend.position = "none") +
  xlim(c(0, 5))

#ggsave("Output/Figures/RLS_sites_pee.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# Visualize salinity
ggplot(bottles_f_may2021, aes(x = sal_est, y = site)) +
  geom_point() +
  labs(x = "Salinity Estimate", y = "Site")

# histogram
ggplot(bottles, aes(x = sal_est)) +
  geom_histogram(colour="black", fill="white") +
  geom_vline(aes(xintercept = 31.8),
             color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept = 27.7),
             color="red", linetype="dashed", size=1)

# Sal vs pee
ggplot(bottles_f_may2021, aes(sal_est, nh4_conc)) + geom_point() +
  geom_smooth(method = lm)

# Visualize temperature
ggplot(bottles, aes(x = temp_est, y = site)) +
  geom_point() +
  labs(x = "Salinity Estimate", y = "Site")

ggplot(bottles, aes(x = temp_est)) +
  geom_histogram(colour="black", fill="white")


# Make a map -----
# Load not great shapefile
potato_map <- sf::st_read("Data/Shapefiles/eez.shp") %>%
  st_sf() %>%
  st_set_crs(4326)

# Define colours
blue <- paste("#b9d1df", sep="")

# Make map without pies, just scaling size of point to %
sf_use_s2(FALSE)

# coords
rls_coords <- rls_nh4 %>%
  group_by(site_ID) %>%
  summarise(nh4_mean = mean(nh4_conc)) %>%
  left_join(
    rls %>% select(site_ID, site_name, latitude, longitude, survey_latitude, survey_longitude,) %>% 
      unique()
    ) %>%
  left_join(fish_biomass %>%
              group_by(site_ID) %>%
              summarize(mean_weight = mean(weight_sum/1000)) ) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  mutate(xjit = ifelse(site_ID == "BMSC6" | 
                         site_ID == "BMSC3" |
                         site_ID == "BMSC21" |
                         site_ID == "BMSC15" |
                         site_ID == "BMSC11" |
                         site_ID == "BMSC19", -0.015, 0.015))

## Create legend
#site_names <- as.list(paste(coords$site_num, coords$`Site name`, sep = ": "))

# Make map
ggplot() +
  geom_sf(data = potato_map, fill = blue, colour = "white") +
  geom_sf(data = rls_coords, 
          colour = "black",
          pch = 21,
          alpha = 0.9,
          aes(size = mean_weight,
              fill = nh4_mean)) +
  coord_sf(xlim = c(-125.4, -125.0), ylim = c(48.80, 49), expand = FALSE)  +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 13, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        panel.grid.major = element_line(color = "white")) +
  viridis::scale_fill_viridis(option="magma", direction = -1) +
  labs(x = "Longitude", y = "Latitude",
       fill = expression(paste("NH"[4]^" +","(umol)")),
       size = "Fish biomass (kg)")  +
  geom_text(data = rls_coords,
            aes(x = survey_longitude, y = survey_latitude, 
                label = site_ID),
            size = 3,
            nudge_x = rls_coords$xjit)


            
#ggsave("Pub_figs/Fig.1.png", device = "png",
#       height = 150, width = 250, units = c("mm"), dpi = 600)