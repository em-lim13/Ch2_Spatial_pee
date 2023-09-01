# Code to look at all the spatial data together
# June 2, 2023
# Em Lim

# Load packages and data ----
library(tidyverse)
library(visreg)
library(ggplot2)
theme_set(theme_bw())
library(PNWColors)
library(ggeffects)
library(lubridate)
library(sf)
library(lmerTest)
source("Code/theme_black.R")

# Load data

# RLS nh4 from 2021
# This is the spring blitz, the june samples, and the July samples
# The may samples have estimated matrix effects, not great
# The June and July samples had a own matrix spike from each site but were still compared to DI
bottles2021 <- read_csv("Output/Output_data/RLS_nh4_2021") 

# RLS from 2022
# These had a own matrix spike from each site but were still compared to DI
# They came out quite negative maybe because of SFU DI, and maybe because of temperature difference between the standards and samples
# I took the lowest nh4 reading and added it to everything else to "set" that sample to 0 and bump everything up
# Maybe an underestimation
bottles2022 <- read_csv("Output/Output_data/RLS_nh4_2022")

# RLS from 2023
# Did the full "proper" Taylor protocol with standard bottles + BF from each site
bottles2023 <- read_csv("Output/Output_data/RLS_nh4_2023")

# combine these three years into one!
rls_nh4_3years <- rbind(bottles2021, bottles2022, bottles2023) %>%
  mutate(year = as.factor(year)) %>%
  filter(month == "May") %>%
  select(-matrix) 


# Add the RLS data ----
# Load the RLS data from the website!

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

# Join all data together
rls <- rbind(fish, invert)

# extract just one row per survey to join with the pee data
rls_survey_depths <- rls %>%
  select(site_ID, year, depth, survey_id, hour) %>%
  rename(survey_depth = depth) %>%
  unique()

# Biomass calculations -------

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
  mutate(size_class = if_else(size_class == 187.5, 87.5, size_class), # shrink that one wolf eel
    weight_per_fish_g = case_when(
    # Gobies
    species_name == "Rhinogobiops nicholsii" ~ exp(log(0.01047) + 3.03*log(size_class)),
    
    # Greenlings
    species_name == "Hexagrammos decagrammus" ~ exp(log(0.00813) + 3.13*log(size_class)),
    species_name == "Hexagrammos stelleri" ~ exp(log(0.00692) + 3.16*log(size_class)),
    species_name == "Oxylebius pictus" ~ exp(log(0.01122) + 3.04*log(size_class)),
    species_name == "Ophiodon elongatus" ~ exp(log(0.00389) + 3.12*log(size_class)),
    species_name == "Hexagrammos spp." ~ exp(log(0.00813) + 3.13*log(size_class)), #
    
    # Rockfish
    species_name == "Sebastes melanops" ~ exp(log(0.01000) + 3.09*log(size_class)),
    species_name == "Sebastes caurinus" ~ exp(log(0.01000) + 3.09*log(size_class)),
    species_name == "Sebastes flavidus" ~ exp(log(0.01000) + 3.09*log(size_class)),
    species_name == "Sebastes maliger" ~ exp(log(0.01000) + 3.09*log(size_class)),
    species_name == "Sebastes nebulosus" ~ exp(log(0.01000) + 3.09*log(size_class)),
    species_name == "Sebastes spp." ~ exp(log(0.01000) + 3.09*log(size_class)),
    species_name == "Sebastes spp. juv" ~ exp(log(0.01000) + 3.09*log(size_class)),
    species_name == "Sebastes pinniger" ~ exp(log(0.01000) + 3.09*log(size_class)),
    
    # Sculpins
    species_name == "Jordania zonope" ~ exp(log(0.00389) + 3.12*log(size_class)),
    species_name == "Artedius harringtoni" ~ exp(log(0.00631) + 3.15*log(size_class)),
    species_name == "Artedius lateralis" ~ exp(log(0.00631) + 3.15*log(size_class)),
    species_name == "Artedius fenestralis" ~ exp(log(0.00631) + 3.15*log(size_class)),
    species_name == "Hemilepidotus hemilepidotus" ~ exp(log(0.00631) + 3.15*log(size_class)),
    species_name == "Cottidae spp." ~ exp(log(0.00631) + 3.15*log(size_class)),
    species_name == "Enophrys bison" ~ exp(log(0.00794) + 3.13*log(size_class)),
    species_name == "Rhamphocottus richardsonii" ~ exp(log(0.01995) + 3.01*log(size_class)),
    species_name == "Scorpaenichthys marmoratus" ~ exp(log(0.00389) + 3.12*log(size_class)),
    species_name == "Oligocottus maculosus" ~ exp(log(0.00631) + 3.15*log(size_class)),
    species_name == "Leptocottus armatus" ~ exp(log(0.01096) + 3.19*log(size_class)),
    species_name == "Blepsias cirrhosus" ~ exp(log(0.00631) + 3.14*log(size_class)),
    species_name == "Myoxocephalus polyacanthocephalus" ~ exp(log(0.00832) + 3.14*log(size_class)),
    species_name == "Myoxocephalus aenaeus" ~ exp(log(0.00832) + 3.14*log(size_class)),
    species_name == "Asemichthys taylori" ~ exp(log(0.00631) + 3.15*log(size_class)),
    
    #Perch
    species_name == "Embiotoca lateralis" ~ exp(log(0.01950) + 2.97*log(size_class)),
    species_name == "Rhacochilus vacca" ~ exp(log(0.01950) + 2.97*log(size_class)),
    species_name == "Brachyistius frenatus" ~ exp(log(0.01318) + 3.05*log(size_class)),
    species_name == "Cymatogaster aggregata" ~ exp(log(0.01950) + 2.97*log(size_class)),
    species_name == "Embiotocidae spp." ~ exp(log(0.01950) + 2.97*log(size_class)),
    species_name == "Percidae spp." ~ exp(log(0.01950) + 2.97*log(size_class)),
    
    # Gunnels + gunnel-like fish
    species_name == "Anarrhichthys ocellatus" ~ exp(log(0.00398) + 3.17*log(size_class)),
    species_name == "Apodichthys flavidus" ~ exp(log(0.00102) + 3.06*log(size_class)),
    species_name == "Pholis ornata" ~ exp(log(0.00162) + 3.19*log(size_class)),
    species_name == "Pholis laeta" ~ exp(log(0.00162) + 3.19*log(size_class)),
    species_name == "Pholis clemensi" ~ exp(log(0.00162) + 3.19*log(size_class)),
    species_name == "Pholis spp." ~ exp(log(0.00162) + 3.19*log(size_class)),
    species_name == "Lumpenus sagitta" ~ exp(log(0.00129) + 2.99*log(size_class)),
    species_name == "Chirolophis nugator" ~ exp(log(0.00372) + 3.16*log(size_class)),
    
    #Misc
    species_name == "Liparis florae" ~ exp(log(0.00525) + 3.15*log(size_class)), #snailfish
    species_name == "Aulorhynchus flavidus" ~ exp(log(0.00263) + 3.14*log(size_class)), #tubesnout
    species_name == "Syngnathus leptorhynchus" ~ exp(log(0.00028) + 3.18*log(size_class)), #pipefish
    species_name == "Clupea pallasii" ~ exp(log(0.00603) + 3.13*log(size_class)), #herring
    species_name == "Gasterosteus aculeatus" ~ exp(log(0.00977) + 3.09*log(size_class)), #stickleback
    species_name == "Porichthys notatus" ~ exp(log(0.00562) + 3.16*log(size_class)), #plainfin
    species_name == "Gibbonsia metzi" ~ exp(log(0.00513) + 3.06*log(size_class)),
    species_name == "Citharichthys stigmaeus" ~ exp(log(0.00759) + 3.15*log(size_class)),
    TRUE ~ as.numeric(NA)),
    biomass_per_indiv = biomass/total,
    weight_size_class_sum = weight_per_fish_g*total) %>%
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


# Join pee + survey data ------

# just keep the surveys where I have a nh4 AND rls survey on the same transect
rls_nh4 <- rls_nh4_3years %>% # remove BMSC6 for now bc it's hard
  left_join(rls_survey_depths, by = c("site_ID", "year")) %>%
  mutate(correct = case_when(
    site_ID== "BMSC6" &year =="2022" &depth == "8.5" &survey_depth== "6.5" ~ "no",
    site_ID== "BMSC6" &year =="2022" &depth == "5.5" &survey_depth== "9" ~ "no",
    site_ID== "BMSC6" &year =="2022" &depth == "6" &survey_depth== "9" ~ "no",
    site_ID== "BMSC1" &year == "2021" & survey_depth == "6.5" ~ "no",
    site_ID== "BMSC1" &year == "2022" & survey_depth == "4.7" ~ "no",
    site_ID== "BMSC5" &year == "2022" & survey_depth == "6.2" ~ "no",
    site_ID== "BMSC11" &year == "2022" & survey_depth == "8.5" ~ "no",
    site_ID== "BMSC12" &year == "2022" & survey_depth == "9" ~ "no",
    site_ID== "BMSC11" &year == "2023" & survey_depth == "5.5" ~ "no",
    site_ID== "BMSC12" &year == "2023" & survey_depth == "6.5" ~ "no",
    site_ID== "BMSC24" &year == "2023" & survey_depth == "7.5" ~ "no",
    site_ID== "BMSC25" &year == "2023" & survey_depth == "5.5" ~ "no",
    site_ID== "BMSC26" &year == "2023" & survey_depth == "9.5" ~ "no",
    site_ID== "BMSC27" &year == "2023" & survey_depth == "7" ~ "no",
    site_ID== "BMSC1" &year == "2023" & survey_depth == "8" ~ "no",
    site_ID== "BMSC5" &year == "2023" & survey_depth == "8" ~ "no",
    site_ID== "BMSC6" &year == "2023" & survey_depth == "5.5" ~ "no",
    site_ID== "BMSC8" &year == "2023" & survey_depth == "7.5" ~ "no",
    TRUE ~ as.character("yes")
    # cut out the rls transects I didn't directly measure nh4 on 
  )) %>%
  filter(correct == "yes") %>%
  select(-c(correct, hour))

# tricky ones:
### BMSC1 2021: two survey depths, one nh4 sample 
### BMSC6 2022: two nh4 samples at 2 depths
### BMSC1 2022: two survey depths, one nh4 sample 
### BMSC5 2022: two surveys at same depth, one nh4 sample
# I was on the later shallower survey
### BMSC6 2022: two surveys at same depth, two nh4 samples!
### BMSC11 2022: two survey depths, one nh4 sample 
### BMSC12 2022: two survey depths, one nh4 sample 
### BMSC11 2023: two survey depths, one nh4 sample 
# I tried to take nh4 samples between the shallower and deeper surveys
### BMSC12 2023: two survey depths, one nh4 sample 
# I tried to take nh4 samples between the shallower and deeper surveys
# BMSC24 2023: two survey depths, one nh4 sample
# I tried to take nh4 samples between the shallower and deeper surveys
# BMSC25 2023: two survey depths, one nh4 sample
# these two were back to back, don't average
# BMSC26 2023: two survey depths, one nh4  (only 2 pee reps tho)
# these two were back to back, don't average
# I was with the 10 m team that got in first
# BMSC27 2023: two survey depths, one nh4 sample
# also back to back, don't average
# BMSC1 2023: two survey depths, one nh4 sample
# BMSC5 2023: two survey depths, one nh4 sample
# BMSC6 2023: two survey depths, one nh4 sample
# the shallower team was way shallower
# BMSC8 2023: two survey depths, one nh4 sample

# So I need to decide what to do about repeated surveys
# Is the nh4 sample I took specific to the transect I took it on?
# Which means I should cut the RLS surveys that aren't the transects where the pee samples were taken
# Or would I want to average the two transects and take the mean biomass from the two to relate to the overall pee sample

# For simplicity I think I just want to keep the RLS survey from the transect where the pee is from

# so basically I need a better way to join up the pee samples with the transect they were taken on

# Site level averaging -----
nh4_avg <- rls_nh4 %>%
group_by(survey_id) %>%
  mutate(nh4_avg = mean(nh4_conc),
         year = as.factor(year)) %>%
  ungroup()  # just keep the samples from the annual RLS spring surveys


# How much fish biomass does each site have?
fish_biomass <- fishes %>%
  group_by(survey_id) %>%
  summarize(weight_sum = sum(weight_size_class_sum))
# One row per transect
  # Per site, per depth, per year
  # Same as grouping by survey_ID

# Join pee + fish biomass data
rls_final <- nh4_avg %>%
  left_join(fish_biomass, by = c("survey_id"))

### EVERYTHING RUNS TO HERE. I CHANGED RLS_FINAL SO MIGHT NEED TO CHANGE FILE NAMES FROM HERE DOWN


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

no_sizes <- fish %>%
  filter(size_class == "0") %>%
  filter(species_name != "Bolinopsis infundibulum") %>%
  filter(species_name != "Pleuronichthys coenosus") %>%
  filter(species_name != "Pleurobrachia bachei") %>%
  filter(species_name != "Polyorchis penicillatus") %>%
  select(-year)

#write_csv(no_sizes, "Output/Output_data/missing_fish_sizes.csv")

# New species
# Armina californica
# Berthella californica and Berthella chacei
# Double check these, I think only one is local 
# Berthella chacei is the local one
# Cadlina sylviaearleae
# Triopha catalinae (we've seen modesta before)


# check out goby town aka why it's sooooo fishy
goby <- fishes %>%
  filter(site_name == "Goby Town") %>%
  filter(depth == "5.7")


# Which species of fishes did we see
fishes_sp <- fish %>% 
  count(species_name, family)

# What's the most abundant species?
rls_abundant_species <- rls %>% 
  count(species_name)

# Which site has the most species?
rls_richness <- rls %>%
  group_by(site_name) %>%
  summarize(species_richness = n_distinct(species_name))%>%
  arrange(desc(species_richness))

# where are the most urchins???
# Can swap for any species
rls_urchins <- rls %>%
  filter(species_name == "Mesocentrotus franciscanus") %>%
  group_by(site_name) %>%
  summarize(urchins = sum(total)) %>%
  arrange(desc(urchins))

# So basically you can manipulate the data to summarize whatever thing you're interested in, and then you can join that summarized data with another df of interest by the common site code/site name :)

# Which sites have the lowest and highest pee
rls_pee_summary <- rls_nh4 %>%
  group_by(site) %>%
  summarize(mean_nh4 = mean(nh4_conc)) %>%
  arrange(desc(mean_nh4)) 

# rank each site by the year
rank_2021 <- bottles2021 %>%
  group_by(site) %>%
  summarise(nh4_avg2021 = mean(nh4_conc)) %>%
  arrange(desc(nh4_avg2021)) %>% 
  mutate(rank2021 = 1:22,
         grade2021 = 100 - (100*rank2021/22))

rank_2022 <- bottles2022 %>%
  group_by(site) %>%
  summarise(nh4_avg2022 = mean(nh4_conc)) %>%
  arrange(desc(nh4_avg2022)) %>% 
  mutate(rank2022 = 1:19,
         grade2022 = 100 - (100*rank2022/19))

rank_2023 <- bottles2023 %>%
  group_by(site) %>%
  summarise(nh4_avg2023 = mean(nh4_conc)) %>%
  arrange(desc(nh4_avg2023)) %>% 
  mutate(rank2023 = 1:20,
         grade2023 = 100 - (100*rank2023/20))

rank <- rank_2021 %>%
  left_join(rank_2022, by= "site") %>%
  left_join(rank_2023, by= "site") %>%
  select(grade2021, grade2022, grade2023, site) %>%
  rowwise() %>%
  mutate(avg_grade = mean(c(grade2021, grade2022, grade2023), na.rm = TRUE))


# What were the matrix effects like at the end of 2021 and in 2022?
rls_nh4 %>%
  filter(year == "2022") %>%
  summarise(matrix = mean(matrix))
# 14.6

rls_nh4 %>%
  filter(year == "2021") %>%
  filter(month != "May") %>%
  summarise(matrix = mean(matrix))
# 7.82
# Go back to the 2021 data and set the ME to 7.82

#What is the tide height at those sites
rls_tide_summary <- rls_nh4 %>%
  group_by(site) %>%
  summarize(depth = mean(depth)) %>%
  arrange(desc(depth))



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