# Code to play with RLS data
# Em Lim
# updated Oct 11, 2022

#### I'm not using this anymore, migrated following code to RLS_nh4_analysis.R

# Load packages
library(tidyverse)
library(visreg)
library(ggplot2)
library(RColorBrewer)
library(ggeffects)
library(lubridate)
theme_set(theme_bw())

# Load data -------------
# Load the RLS data from the website!

# Just the pelagic and cryptic fish
fish <- read_csv("Data/RLS/RLS_data/reef_fish_abundance_and_biomass.csv",
                   show_col_types = FALSE) %>%
  rbind(read_csv("Data/RLS/RLS_data/cryptobenthic_fish_abundance.csv",
                 show_col_types = FALSE)) %>%
  as.data.frame() %>%
  mutate(survey_date = ymd(survey_date),
         year = year(survey_date)) %>%
  filter(month(survey_date) == 4 | month(survey_date) == 5) # Just the RLS blitz data for now

# Just the mobile inverts
invert <- read_csv("Data/RLS/RLS_data/mobile_macroinvertebrate_abundance.csv",
                   show_col_types = FALSE) %>%
  as.data.frame() %>%
  mutate(survey_date = ymd(survey_date),
         year = year(survey_date)) %>%
  filter(month(survey_date) == 4 | month(survey_date) == 5) # Just the RLS blitz data for now

# Join all data together
rls <- rbind(fish, invert)

# Data checks

goby <- fish %>%
  filter(species_name == "Rhinogobiops nicholsii") %>%
  filter(site_name == "Goby Town") %>%
  filter(depth == "5.7") %>%
  select(site_name, depth, method, block, species_name, total, size_class)

# Neat data manipulations ------

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

# Mass estimates -------

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
  mutate(weight_g = case_when(
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
    biomass_per_indiv = biomass/total) %>%
  filter(species_name != "Bolinopsis infundibulum") %>%
  filter(species_name != "Pleuronichthys coenosus") %>%
  filter(species_name != "Pleurobrachia bachei") %>%
  filter(species_name != "Polyorchis penicillatus") %>% # Filter inverts
  filter(species_name != "Phoca vitulina") %>% # Remove seal
  filter(species_name != "Actinopterygii spp.") %>% # Remove unidentif fish
  filter(species_name != "3)	Myoxocephalus aenaeus") # Remove east coast fish


# I also have size data for inverts "Haliotis kamtschatkana", "Crassadoma gigantea", "Pycnopodia helianthoides", "Polyorchis penicillatus", "Bolinopsis infundibulum", Pleurobrachia bachei, Pleuronichthys coenosus

# We also measured sea cucumbers a few times???? Where did that data go?
# We also have urchin size data from this year

# Some people counted M2 fishes on M1, but not always so I don't actually want goby and sculpin counts from M1
# Figure out what to do here




# Old way of loading data -----

# updated_RLS_2021.csv is the most updated 2021 RLS file
# updated_RLS_KCCA.csv is the updated 2021 file + an incomplete record from Kieran and Claire's data

rls <- read_csv("Data/RLS/RLS_data/updated_RLS_2021.csv") %>%
  as.data.frame() %>%
  select(-1) %>% # cuts the first column which is blank
  select(-Inverts) %>% # cuts the inverts column which is just NAs
  filter(Method != 0) %>% # get rid of all method 0's
  slice(2:n()) %>% # cuts the first blank row
  rename(
    site_ID = `Site No.`,
    site_name = `Site Name`, 
    common_name = `Common name`
  )  %>% # Rename columns with spaces
  mutate(Species = str_to_sentence(Species),
         common_name = str_to_sentence(common_name),
         Date = dmy(Date)) %>% 
  filter(Species != "Debris - Metal") %>%
  filter(Species != "Debris - Other") %>%
  filter(Species != "Debris - Wood") %>%
  filter(Species != "Debris - Glass") %>%
  filter(Species != "Debris - Fishing Gear") %>%
  filter(Species != "Debris - Zero") %>%
  filter(Species != "Survey not completed") %>%
  filter(Species != "No species found") %>%
  filter(Species != "NOT PRESENT") %>% # Cut the non-animal species
  filter(Date < "2021-05-21") # just the May data

