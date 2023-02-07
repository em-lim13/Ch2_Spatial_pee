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
# Capitalize first letter of each code
# Cali cukes named inconsistently
# Blood stars can't be identified to species
    # Pacific blood star should be Henricia spp.
# Hermit crabs. Ugh. 
  # Jasmin's notation: Paguroidea spp.	Unidentified Hermit crab	Hermit
# Black rockfish with latin name Sebastes flavidus at KCCA1
  # Looked at Kieran's sheet, it's supposed to be a black rockfish
# 2 unidentified greenlings on KCCA1
  # datasheet called them "J Green" and they're 5 cm
  # Juvenile greenlings
  # Jasmin has used: Hsp	Hexagrammos sp.	Greenling juv. for these
# 5 unidentified rockfish on KCCA1
  # datasheet calls them "J RF" 5 cm and 7.5 cm
  # Juv rockfish
  # Sebastes spp. juv	Rockfish YOY	Sjuv
# Include juv in species names for those
# Add common names for the other unidentified species
# All common names capital first letter, lower case for rest
# Purple sea star vs Ochre star? Leave for now

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
  mutate(Code = str_to_sentence(Code)) %>% # capitalize first letter of each code
  mutate(Species = case_when(
    common_name == "California sea cucumber" | common_name == "Californian sea cucumber" ~ "Apostichopus californicus",
    Species == "Sebastes flavidus" & common_name == "Black rockfish" ~ "Sebastes melanops",
    Species == "Henricia leviuscula" ~ "Henricia spp.",
    Species == "Hexagrammos spp." ~ "Hexagrammos spp. juv",
    Species == "Pagurus spp." ~ "Paguroidea spp.",
    Species == "Sebastes spp." ~ "Sebastes spp. juv",
    Species == "Oncorhynchus" | Species == "oncorhynchus" ~ "Oncorhynchus spp.",
    TRUE ~ as.character(Species))) %>% 
  mutate(common_name = case_when(
    Species == "Apostichopus californicus" ~ "California sea cucumber",
    Species == "Crassadoma gigantea" ~ "Giant rock scallop",
    Species == "Hexagrammos spp. juv" ~ "Greenling YOY",
    Species == "Sebastes spp. juv" ~ "Rockfish YOY",
    Species == "Henricia spp." ~ "Unidentified blood star",
    Species == "Tonicella spp." ~ "Unidentified lined chiton",
    Species == "Pandulus spp." ~ "Unidentified shrimp",
    Species == "Paguroidea spp." ~ "Unidentified hermit crab",
    Species == "Clupea pallasii" ~ "Pacific herring",
    Species == "Oncorhynchus spp." ~ "Pacific salmon",
    Species == "Stylasterias forreri" ~ "Velcro sea star",
    TRUE ~ as.character(common_name))) %>%
  mutate(Code = case_when(
    Species == "Apostichopus californicus" ~ "Aca",
    Species == "Sebastes melanops" ~ "Sme",
    Species == "Henricia spp." ~ "Hsp",
    Species == "Hexagrammos spp. juv" ~ "Hjuv",
    Species == "Tonicella spp." ~ "Tsp",
    Species == "Pandulus spp." ~ "Psp",
    Species == "Paguroidea spp." ~ "Hermit",
    Species == "Sebastes spp. juv" ~ "Sjuv",
    Species == "Mesocentrotus franciscanus" ~ "Mfr",
    Species == "Oncorhynchus spp." ~ "Osp",
    TRUE ~ as.character(Code)))

# check to make sure each species only has one name and one code! 
rls_last_check <- rls_tidy %>%
  count(common_name, Species, Code) # Looks good
