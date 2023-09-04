# Code to create functions to calculate concentration of ammonium via fluorometric methods
# Using the Taylor protocol 1 + 2 methods

# First load packages
library(tidyverse)
library(patchwork)

# function Nikola built for me!!!
pee_calc <- function(data_path) {
  data_bf <- read_csv(data_path, show_col_types = FALSE) %>%
    filter(sample == "bf") %>% # just the BF sample
    transmute(bf = (FLU1 + FLU2 + FLU3)/3,
              site_code = site_code) 
  
  # Load the standards and the sample data
  data <- read_csv(data_path, show_col_types = FALSE) %>%
    as.data.frame() %>%
    filter(sample != "bf") %>% # cut the bf sample
    left_join(data_bf, by = "site_code") %>% # Add the background FLU as a column
    mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
           mean_FLU = mean_FLU_unadj - bf) # correct the samples with this bf value
  
  # Build the standard curve for protocol 1!
  # We need to estimate how much ammonium is in the standard water
  standard <- data %>%
    filter(sample == "standard") %>% 
    mutate(nh4_added_umol = nh4_vol_uL/1e6 * 200, #amount of NH4
           total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
           nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
  #of NH4 in seawater sample
  
  #model for protocol 1
  mod <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard)
  summary(mod)
  print(summary(mod)$adj.r.squared) # Check Adjusted R-squared value, we're looking for at least 0.98
  
  #graph protocol 2 standard curve
  proto1plot<- (ggplot(standard, aes(nh4_conc_final_umol_L, mean_FLU)) +
                  geom_point() +
                  geom_smooth(method = lm, se = FALSE))
  # Inspect this curve to make sure nothing is wonky!
  
  #save coefficients
  int <- coef(mod)[1]
  slope <- coef(mod)[2]
  
  # calculate ammonium in standards using these coefficients
  standard_b <- standard %>%
    mutate(int = int, 
           slope = slope, 
           nh4_conc = abs(int/slope),
           true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)
  
  # Now move into protocol 2! Use the standards to calculate ammonium in samples
  # BF-corrected FLU of standard curve against calculated NH4 concentration
  mod2 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b)
  summary(mod2)
  
  #graph protocol 2 standard curve
  proto2plot<- (ggplot(standard_b, aes(true_nh4_conc, mean_FLU)) +
          geom_point() +
          geom_smooth(method = lm, se = FALSE))
  # Inspect this curve to make sure nothing is wonky!
  
  # inspect standard curves
  print(proto1plot + proto2plot)

  # pull coefficients from model
  int2 <- coef(mod2)[1]
  slope2 <- coef(mod2)[2]
  
  # calculate how much ammonium is in the samples
  samples <- data %>%
    filter(nh4_vol_uL == "0") %>%
    mutate(nh4_conc = (mean_FLU - int2)/slope2)
}








pee_calc2 <- function(data_path) {
  data_bf <- read_csv(data_path, show_col_types = FALSE) %>%
    filter(sample == "bf") %>% # just the BF sample
    transmute(bf = (FLU1 + FLU2 + FLU3)/3,
              site_code = site_code) 
  
  # Load the standards and the sample data
  data <- read_csv(data_path, show_col_types = FALSE) %>%
    as.data.frame() %>%
    filter(sample != "bf") %>% # cut the bf sample
    left_join(data_bf, by = "site_code") %>% # Add the background FLU as a column
    mutate(mean_FLU_unadj = (FLU1 + FLU2 + FLU3)/3,
           mean_FLU = mean_FLU_unadj - bf) # correct the samples with this bf value
  
  # Build the standard curve for protocol 1!
  # We need to estimate how much ammonium is in the standard water
  standard <- data %>%
    filter(sample == "standard") %>% 
    mutate(nh4_added_umol = nh4_vol_uL/1e6 * 2014, #amount of NH4
           total_vol_L = nh4_vol_uL/1e6 + 0.04, #new volume of sample + NH4
           nh4_conc_final_umol_L = nh4_added_umol / total_vol_L) #concentration
  #of NH4 in seawater sample
  
  #model for protocol 1
  mod <- lm(mean_FLU ~ nh4_conc_final_umol_L, data = standard)
  summary(mod)
  print(summary(mod)$adj.r.squared) # Check Adjusted R-squared value, we're looking for at least 0.98
  
  #graph protocol 2 standard curve
  proto1plot<- (ggplot(standard, aes(nh4_conc_final_umol_L, mean_FLU)) +
                  geom_point() +
                  geom_smooth(method = lm, se = FALSE))
  # Inspect this curve to make sure nothing is wonky!
  
  #save coefficients
  int <- coef(mod)[1]
  slope <- coef(mod)[2]
  
  # calculate ammonium in standards using these coefficients
  standard_b <- standard %>%
    mutate(int = int, 
           slope = slope, 
           nh4_conc = abs(int/slope),
           true_nh4_conc = nh4_conc + nh4_conc_final_umol_L)
  
  # Now move into protocol 2! Use the standards to calculate ammonium in samples
  # BF-corrected FLU of standard curve against calculated NH4 concentration
  mod2 <- lm(mean_FLU ~ true_nh4_conc, data = standard_b)
  summary(mod2)
  
  #graph protocol 2 standard curve
  proto2plot<- (ggplot(standard_b, aes(true_nh4_conc, mean_FLU)) +
                  geom_point() +
                  geom_smooth(method = lm, se = FALSE))
  # Inspect this curve to make sure nothing is wonky!
  
  # inspect standard curves
  print(proto1plot + proto2plot)
  
  # pull coefficients from model
  int2 <- coef(mod2)[1]
  slope2 <- coef(mod2)[2]
  
  # calculate how much ammonium is in the samples
  samples <- data %>%
    filter(sample != "standard") %>%
    mutate(nh4_conc = (mean_FLU - int2)/slope2)
}


# Biomass function ----

length_to_weight <- function(datafile){
  new_data <- datafile %>%
    mutate(weight_per_fish_g = case_when(
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
           weight_size_class_sum = weight_per_fish_g*total)
  
}


# Joining surveys up by depth -----

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

depth_function <- function(datafile){
  new_data <- datafile %>%
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
}
