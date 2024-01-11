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

data_s <- read_csv("Output/Output_data/kelp_final.csv")

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
load("~/Documents/PhD/Ch2_Spatial_pee/Data/bc_map.Rdata")
bc_map <- slice # rename

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

# just slack and ebb
coords_slack <- rls_coords %>%
  filter(tide_cat != "Flood") %>%
  group_by(site_code) %>%
  mutate(nh4_avg = mean(nh4_avg)) %>%
  ungroup()


# Kelp data coords -----
kelp_coords <- data_s %>%
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
           map_file = bc_map,
           invert = FALSE) +
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
           map_file = bc_map,
           invert = FALSE) +
  guides(pch = "none")

 ggsave("Output/Figures/rls_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)


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
           map_file = bc_map,
           invert = FALSE) +
  guides(pch = "none",
         fill = "none")

# ggsave("Output/Figures/kelp_site_map.png", device = "png", height = 9, width = 16, dpi = 400)


# Kelp nh4 data!
map_daddy(lat_min = -125.4,
           lat_max = -125.0, 
           long_min = 48.80, 
           long_max = 49, 
           coord_data = all_coords %>% filter(Habitat == "Kelp"), 
           nh4_var = nh4_avg, 
           point_size = 9, 
           kelp_var = Habitat,
           map_file = potato_map,
           invert = TRUE) +
  guides(pch = "none")

#ggsave("Output/Figures/kelp_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)



# All sites
map_daddy(lat_min = -125.4,
           lat_max = -125.0, 
           long_min = 48.80, 
           long_max = 49, 
           coord_data = all_coords, 
           nh4_var = nh4_avg, 
           kelp_var = Habitat,
           point_size = 9, 
           map_file = bc_map,
           invert = FALSE)

ggsave("Output/Figures/all_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)





# Graveyard ------

# use site map! -----
site_map(
  lat_min = -125.3,
  lat_max = -125,
  long_min = 48.8,
  long_max = 48.9,
  coord_data = all_coords,
  map_data = potato_map,
  add_points = TRUE,
  add_annotate = TRUE
)


# just rls sites
barkley_rls <- site_map(lat_min = -125.6, lat_max = -124.95, long_min = 48.75, long_max = 49.1,
                        coord_data = all_coords %>%
                          filter(Habitat == "Reef"), map_data = potato_map,
                        add_points = TRUE, add_annotate = TRUE) +
  guides(pch = "none")

van_isle_rls <- inset_map(rect_xmin = -125.6, rect_xmax = -124.95, rect_ymin = 48.75, rect_ymax = 49.1,
                          map_data = potato_map)


barkley_rls + 
  inset_element(
    van_isle_rls, 
    left = 0, 
    bottom = 0.5, 
    right = 0.5, 
    top = 1.03,
    align_to = 'panel'
  )



# Cursed pimple map -----

# need to use full nh4 data, go run the rls analysis script to generate "rls_nh4"
coords_nh4 <- rls_nh4 %>%
  select(site_code, year, nh4_conc) %>%
  left_join(
    # use real coordinates not the RLS rounded ones
    read_csv("Data/RLS/RLS_data/true_coords.csv") 
  ) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  group_by(site_code) %>%
  mutate(nh4_avg = mean(nh4_conc),
         nh4_min = min(nh4_conc),
         nh4_max = max(nh4_conc)) %>%
  ungroup()


# Can I do something cursed? 
ggplot() +
  geom_sf(data = potato_map, fill = blue, colour = "white") +
  geom_sf(data = coords_nh4, 
          colour = "black",
          pch = 21,
          alpha = 0.9,
          size = 13.5,
          aes(fill = nh4_min)) +
  geom_sf(data = coords_nh4, 
          colour = "black",
          pch = 21,
          alpha = 0.9,
          size = 11,
          aes(fill = nh4_avg)) +
  geom_sf(data = coords_nh4, 
          colour = "black",
          pch = 21,
          alpha = 0.9,
          size = 4.5,
          aes(fill = nh4_max)) +
  coord_sf(xlim = c(-125.4, -125.0), ylim = c(48.80, 49), expand = FALSE)  +
  theme_black() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "white")) +
  viridis::scale_fill_viridis(option="magma", direction = -1,
                              guide = guide_colorbar(frame.colour = "white", ticks.colour = "white")) +
  labs(x = "Longitude", y = "Latitude",
       fill = expression(paste("NH"[4]^" +",(mu*M)))) +
  scale_x_continuous(breaks = seq(-125.4, -125.0, by = 0.1))

#ggsave("Output/Figures/cursed_nh4_map.png", device = "png", height = 9, width = 16, dpi = 400)
