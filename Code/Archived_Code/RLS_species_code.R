# Code to play with RLS data
# Em Lim
# updated Oct 11, 2022

#### I'm not using this anymore

# Load packages
library(tidyverse)
library(visreg)
library(ggplot2)
library(RColorBrewer)
library(ggeffects)
library(lubridate)
theme_set(theme_bw())


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

