# New R file for making maps of the RLS and Kelp pee data
# Em Lim
# Oct 17, 2023

# Load packages -----
library(tidyverse)
library(ggplot2)
library(ggspatial)
library(sf)

# Load functions ----
source("Code/Functions.R") # Length to weight function here!

# Load data ----
# these are csv files created by the rls and kelp analysis files

data_map <- read_csv("Output/Output_data/kelp_final.csv")

rls_final <- read_csv("Output/Output_data/rls_final.csv")

kelp_rls <- read_csv("Output/Output_data/kelp_rls.csv")

# load all RLS data
rls_nh4_all <- rbind(read_csv("Output/Output_data/RLS_nh4_2021.csv"),
                 read_csv("Output/Output_data/RLS_nh4_2022.csv"),
                 read_csv("Output/Output_data/RLS_nh4_2023.csv")) %>%
  rename(site_code = site_ID) %>%
  group_by(site_code) %>%
  mutate(nh4_overall_avg_all = mean(nh4_conc)) %>%
  ungroup() %>%
  select(site_code, nh4_overall_avg_all) %>%
  unique()

# Load GREAT shapefile ----
hakai_map <- sf::st_read("Data/Hakaii_coast/COAST_TEST2.shp") %>%
  st_sf() %>%
  st_set_crs(4326)

# Load not great shapefile that renders quicker
potato_map <- sf::st_read("Data/Potato_shapefiles/eez.shp") %>%
  st_sf() %>%
  st_set_crs(4326)

# Make map without pies, just scaling size of point to %
sf_use_s2(FALSE)

# Define colours
blue <- paste("#b9d1df", sep="")


# RLS coords ------
rls_coords <- rls_final %>%
  select(site_code, survey_id, year, nh4_avg, abundance, tide_cat) %>%
  left_join(
    # use real coordinates not the RLS rounded ones
    read_csv("Data/RLS/RLS_data/true_coords.csv") 
  ) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  group_by(site_code) %>%
  mutate(nh4_overall_avg = mean(nh4_avg),
         nh4_min = min(nh4_avg),
         nh4_max = max(nh4_avg)) %>%
  ungroup() %>%
  mutate(Habitat = "Reef") %>%
  left_join(rls_nh4_all)


# Kelp data coords -----
kelp_coords <- data_map %>%
  select(site_code, nh4_in_avg, nh4_out_avg, nh4_avg, kelp_sp, tide_cat) %>%
  left_join(
    kelp_rls %>%
      mutate(Habitat = "Kelp")
  ) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  group_by(site_code)


# Both sets of coords ----
all_coords <- rls_coords %>%
  transmute(site_code = site_code,
            nh4_avg = nh4_overall_avg_all,
            geometry = geometry,
            Habitat = Habitat) %>%
  unique() %>%
  rbind(
    (kelp_coords %>% transmute(nh4_avg = nh4_out_avg,
                               geometry = geometry,
                               Habitat = Habitat)) ) %>%
  # shrink the two sites over 1 uM to 2 uM so the scale is nicer
  mutate(Habitat = factor(as.factor(Habitat), levels = c("Reef", "Kelp")),
         nh4_avg = ifelse(nh4_avg > 2, 2, nh4_avg),
         dummy = 2) 



# Make maps! ------

# Barkley Sound map
map_daddy_np(lat_min = -127,
          lat_max = -123, 
          long_min = 48.5, 
          long_max = 51, 
          map_file = potato_map,
          invert = FALSE) +
  # add rectangle for zoomed in part
  geom_rect(aes(xmin = -125.4, xmax = -125.0, ymin = 48.80, ymax = 49),
            color = "red", fill = NA, inherit.aes = FALSE) 

# ggsave("Output/Figures/barkley_sound_map.png", device = "png", height = 9, width = 16, units = "cm", dpi = 400)


# Making RLS maps -----

# RLS site map
map_daddy(lat_min = -125.4,
          lat_max = -125.0, 
          long_min = 48.75, 
          long_max = 49, 
          coord_data = all_coords %>% filter(Habitat == "Reef"), 
          nh4_var = dummy, 
          point_size = 4, 
          kelp_var = Habitat,
          map_file = hakai_map) +
  guides(pch = "none",
         fill = "none")

# ggsave("Output/Figures/rls_site_map.png", device = "png", height = 9, width = 16, dpi = 400)


# RLS nh4 data!
map_daddy(lat_min = -125.4,
          lat_max = -125.0, 
          long_min = 48.80, 
          long_max = 49, 
          coord_data = all_coords %>% filter(Habitat == "Reef"), 
          nh4_var = nh4_avg, 
          point_size = 9, 
          kelp_var = Habitat,
          map_file = hakai_map,
          invert = FALSE,
          white_background = FALSE) +
  guides(pch = "none")

# ggsave("Output/Figures/rls_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)


# Kelp maps! ----
# Kelp site map
map_daddy(lat_min = -125.4,
          lat_max = -125.0, 
          long_min = 48.80, 
          long_max = 49, 
          coord_data = all_coords %>% filter(Habitat == "Kelp"), 
          nh4_var = dummy, 
          point_size = 4, 
          kelp_var = Habitat,
          map_file = hakai_map,
          invert = FALSE,
          white_background = FALSE) +
  guides(pch = "none",
         fill = "none")

# ggsave("Output/Figures/kelp_site_map.png", device = "png", height = 9, width = 16, dpi = 400)


# Kelp nh4 data!
map_daddy(lat_min = -125.3,
          lat_max = -125.0, 
          long_min = 48.78, 
          long_max = 48.93, 
          coord_data = all_coords %>% filter(Habitat == "Kelp"), 
          nh4_var = nh4_avg, 
          point_size = 9, 
          kelp_var = Habitat,
          map_file = hakai_map,
          invert = FALSE,
          white_background = FALSE) +
  guides(pch = "none")

#ggsave("Output/Figures/kelp_nh4_map2.png", device = "png", height = 9, width = 16, dpi = 400)



# All sites
map_daddy(lat_min = -125.4,
          lat_max = -125.0, 
          long_min = 48.80, 
          long_max = 49, 
          coord_data = all_coords, 
          nh4_var = nh4_avg, 
          kelp_var = Habitat,
          point_size = 7, 
          map_file = hakai_map,
          invert = FALSE,
          white_background = TRUE)

# ggsave("Output/Pub_figs/Fig.all_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)


# calcs -----
# calculate distance between each RLS point
dist <- read_csv("Data/RLS/RLS_data/true_coords.csv") %>%
  mutate(k = 1) %>%
  rename(x = longitude,
         y = latitude)

dist2 <- dist %>%
  full_join(dist, by = "k") %>%
  mutate(distance = sqrt((x.x - x.y)^2 + (y.x - y.y)^2)) 

# greatest distance is between Hosie and Wouwer = 24.18 km
# smallest distance is Baeria Rocks North Island Southside to Baeria Rocks North Island Northside = 65.65 m
# next smallest is Dodger to Taylor = 194.05 m 



# play with making a kml -----
ggplot() +
  geom_sf(data = hakai_map, fill = "white", colour = blue)
