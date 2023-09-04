# Code to quality control RLS data
# Em Lim
# updated Feb 6, 2023

# Load tidyverse
library(tidyverse)

# Load data -------------
rls <- read_csv("Data/RLS/RLS_2022_KDC_CMA.csv") %>%
  select(-1) %>% # cuts the first column which is blank
  select(-Inverts) %>% # cuts the inverts column which is just NAs
  slice(2:n()) %>% # cuts the first blank row
  rename(
    site_ID = `Site No.`,
    site_name = `Site Name`, 
    common_name = `Common name`
  )  # Rename columns with spaces

# Normally would filter non animal and M0 data but I want to make sure all the codes are consistent so let's leave these for now 

# Summarize all the combos of latin and common names and codes
# Each common name should only have ONE latin name and ONE code
rls_common_name_check <- rls  %>% 
  count(common_name, Species, Code)

# Comments:
# Cali cukes named inconsistently
# Blood stars can't be identified to species
    # Pacific blood star should be Henricia spp.
# No more separating juv fish from unidentified fish. 

# Black rockfish with latin name Sebastes flavidus at KCCA1
  # Looked at Kieran's sheet, it's supposed to be a black rockfish

# Add common names for the other unidentified species

#Questions:
# Not sure what the two NOT PRESENT species are at KCCA13, one with code "mfra" and one code "cry"
  # Claire entered them
  # We'll probably filtering out any Species == "NOT PRESENT" from the data but this should be followed up on
# All clown nudibranchs found were reported as Triopha catalinae, but Triopha modesta is supposed to be the more commonly Northern species. Pls confirm all the clown nudis you saw had round tubercules, no branchy ones
  # Paper here: https://link.springer.com/article/10.1007/s12526-020-01107-2

# Methods things
# Chitons aren't counted under the RLS method, were you both counting chitons all the time? The species name is Tonicella, were you both only counting lined chitons?
# Brittle stars also aren't counted under the RLS method, were you both counting them all the time? Or just at KCCA21?
# No debris data??? If no debris found, should enter "Debris - Zero" for each block under M2

# Tidy up the data!!!
rls_tidy <- rls %>%
  mutate(Code = str_to_lower(Code)) %>% # capitalize first letter of each code
  mutate(Species = case_when(
    common_name == "California sea cucumber" | common_name == "Californian sea cucumber" ~ "Apostichopus californicus",
    Species == "Sebastes flavidus" & common_name == "Black rockfish" ~ "Sebastes melanops",
    Species == "Henricia leviuscula" ~ "Henricia spp.",
    Species == "Pagurus spp." ~ "Paguroidea spp.",
    Species == "Oncorhynchus" | Species == "oncorhynchus" ~ "Oncorhynchus spp.",
    Species == "Sebastes spp. juv" ~ "Sebastes spp.",
    TRUE ~ as.character(Species))) %>% 
  mutate(common_name = case_when(
    Species == "Apostichopus californicus" ~ "California sea cucumber",
    Species == "Hexagrammos spp." ~ "Unidentified greenling",
    Species == "Sebastes spp." ~ "Unidentified rockfish",
    Species == "Henricia spp." ~ "Unidentified blood star",
    Species == "Tonicella spp." ~ "Unidentified lined chiton",
    Species == "Pandulus spp." ~ "Unidentified shrimp",
    Species == "Paguroidea spp." ~ "Unidentified hermit crab",
    Species == "Clupea pallasii" ~ "Pacific herring",
    Species == "Oncorhynchus spp." ~ "Unidentified pacific salmon",
    Species == "Stylasterias forreri" ~ "Velcro sea star",
    TRUE ~ as.character(common_name))) %>%
  mutate(Code = case_when(
    Species == "Apostichopus californicus" ~ "aca",
    Species == "Sebastes melanops" ~ "sme",
    Species == "Henricia spp." ~ "henricia",
    Species == "Hexagrammos spp." ~ "hexagrammos",
    Species == "Tonicella spp." ~ "tonicella",
    Species == "Pandulus spp." ~ "pandulus",
    Species == "Paguroidea spp." ~ "hermit",
    Species == "Sebastes spp." ~ "sebastes",
    Species == "Oncorhynchus spp." ~ "oncorhynchus",
    TRUE ~ as.character(Code))) %>%
  filter(Species != "Ophiopholis aculeata") # Filter brittle star

# check to make sure each species only has one name and one code! 
rls_last_check <- rls_tidy %>%
  count(common_name, Species, Code) # Looks good
