# Code to create functions to calculate concentration of ammonium via fluorometric methods
# Using the Taylor protocol 1 + 2 methods

# First load packages
library(tidyverse)
library(patchwork)

# Theme black ----
theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "white"),  
      axis.text.x = element_text(size = base_size*2, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*2, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", linewidth  =  0.2),  
      axis.title.x = element_text(size = base_size*2.5, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size*2.5, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "black",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*2, color = "white"),  
      legend.title = element_text(size = base_size*2.5, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "black"),  
      panel.grid.minor = element_line(color = "black"),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*2.5, color = "white"),  
      strip.text.y = element_text(size = base_size*2.5, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}


theme_white = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "black"),  
      axis.text.x = element_text(size = base_size*2, color = "black", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*2, color = "black", lineheight = 0.9),  
      axis.ticks = element_line(color = "black", linewidth  =  0.2),  
      axis.title.x = element_text(size = base_size*2.5, color = "black", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size*2.5, color = "black", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "white"),  
      legend.key = element_rect(color = "white",  fill = "white"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*1.5, color = "black"),  
      legend.title = element_text(size = base_size*2.5, face = "bold", hjust = 0, color = "black"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "white", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "black"),  
      panel.grid.major = element_line(color = "white"),  
      panel.grid.minor = element_line(color = "white"),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*2.5, color = "black"),  
      strip.text.y = element_text(size = base_size*2.5, color = "black",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "white", fill = "white"),  
      plot.title = element_text(size = base_size*1.2, color = "black"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}


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


# Fish length to weight ----

length_to_weight <- function(datafile){
  new_data <- datafile %>%
    mutate(weight_per_indiv_g = case_when(
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
             TRUE ~ as.numeric(0.5)),
           # set the really big wolf eel weight manually to largest record weight
           # otherwise the calc thinks it's MASSSSSIVE
           weight_per_indiv_g = if_else(size_class == 187.5, 18400, weight_per_indiv_g),
           weight_per_indiv_kg = weight_per_indiv_g/1000,
           weight_size_class_sum = weight_per_indiv_kg*total)
  
}


# Invert length to weight ------

# Mean cuke wet weight = 829 g
# Mean cuke dry weight = 39

# scallops length to weight
# scallop <- read_csv("Data/RLS/scallop_l_w.csv")

# mod <- lm(log(body_weight) ~ log(shell_height), data = scallop)
# summary(mod)
# intercept = log(a) = -8.77259
# slope = b = 2.39442


invert_length_to_weight <- function(datafile){
  new_data <- datafile %>%
    mutate(weight_per_indiv_g = case_when(
      # Echinoderms
      species_name == "Apostichopus californicus" ~ 829, # Cuke avg from my measurements
      # Urchins
      species_name == "Mesocentrotus franciscanus" ~ 209, # 19.38 dry Peters
      species_name == "Strongylocentrotus purpuratus" ~ 104, # 10.71 dry Peters 
      species_name == "Strongylocentrotus droebachiensis" ~ 51, # 2.2 gonad weight Siikavuopio
      # stars
      species_name == "Patiria miniata" ~ 63, # 6.83 g dry weight Peters
      species_name == "Dermasterias imbricata" ~ 63, # 6.83 g dry weight Peters (bat star)
      species_name == "Orthasterias koehleri" ~ 5.5, # Menge 1975 (Leptastarias)
      species_name == "Pisaster ochraceus" ~ 73, # Sanford + p_ochraceus 2019
      species_name == "Evasterias troschelii" ~ 63, # guess Sanford + p_ochraceus 2019
      species_name == "Pycnopodia helianthoides" ~ size_class*0.5, # size data, guess for now!!
      # pycno size is a guess, careful
      species_name == "Stylasterias forreri" ~ 63, # Sanford + p_ochraceus 2019
      species_name == "Mediaster aequalis" ~ 5.5, # Menge 1975 (Leptastarias)
      species_name == "Leptasterias hexactis" ~ 5.5, # Menge 1975
      species_name == "Henricia pumila" ~ 5.5, # Menge 1975 (Leptastarias)
      species_name == "Henricia spp." ~ 5.5, # Menge 1975 (Leptastarias)
      
      # Snails
      species_name == "Pomaulax gibberosus" ~ 255, # 6.85 dry Peters (Wavy snail)
      species_name == "Haliotis kamtschatkana" ~ 0.0000578*(size_class*10)^3.2, # Zhang 2007
      species_name == "Ceratostoma foliatum" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Calliostoma ligatum" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Nucella lamellosa" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Tegula funebralis" ~ 1, # ROUGH from Palmer 1982 + 1988
      # Limpets + scallop
      species_name == "Lottia scutum" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Acmaea mitra" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Megathura crenulata" ~ 60, # 16.7 dry Peters
      species_name == "Crassadoma gigantea" ~ exp(-8.77259 + 2.39442*log(size_class)), #McDonald
      # Nudis
      species_name == "Hermissenda crassicornis" ~ 1, # Megina 2007, also AVILA 1997
      species_name == "Peltodoris nobilis" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Dirona albolineata" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Doris montereyensis" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Dendronotus iris" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Diaulula odonoghuei" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Polycera tricolor" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Triopha modesta" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Doris odhneri" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Acanthodoris hudsoni" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Diodora aspera" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Limacia cockerelli" ~ 1, # ROUGH from Palmer 1982 + 1988
      species_name == "Cadlina luteomarginata" ~ 1, # ROUGH from Palmer 1982 + 1988
      # Crabs
      species_name == "Paguroidea spp." ~ 5, # Rough from Griggiths and Gosselin
      species_name == "Pagurus hemphilli" ~ 5, # Rough from Griggiths and Gosselin
      species_name == "Scyra acutifrons" ~ 1, # Hines
      species_name == "Cancer productus" ~ 200, # rough from my data
      species_name == "Glebocarcinus oregonensis" ~ 5, # Hines H. nudus
      species_name == "Pugettia producta" ~ 46, # Hines
      species_name == "Pagurus beringanus" ~ 5, # Rough from Griggiths and Gosselin
      TRUE ~ as.numeric(0.5)),
      weight_per_indiv_kg = weight_per_indiv_g/1000,
      weight_size_class_sum = weight_per_indiv_kg*total)
}

# just the top 89-45 inverts.... anything that was seen more than 9 times



# Fish home ranges ----

home_range <- function(datafile){
  new_data <- datafile %>%
    mutate(range = case_when(
      # Gobies
      species_name == "Rhinogobiops nicholsii" ~ 0.05, # guessed from painted greenling
      # Greenlings
      species_name == "Hexagrammos decagrammus" ~ 0.3, # 0.1 - 0.5
      species_name == "Hexagrammos stelleri" ~ 0.3, # from kelp greenling 0.1 - 0.5
      species_name == "Oxylebius pictus" ~ 0.02,
      species_name == "Ophiodon elongatus" ~ 28.3, # 3.3 - 498
      species_name == "Hexagrammos spp." ~ 0.3, # from kelp greenling 0.1 - 0.5
      
      # Rockfish
      species_name == "Sebastes melanops" ~ 0.5,
      species_name == "Sebastes caurinus" ~ 0.15,
      species_name == "Sebastes flavidus" ~ 500, # 140-1400 range
      species_name == "Sebastes maliger" ~ 0.05,
      species_name == "Sebastes nebulosus" ~ 0.01,
      species_name == "Sebastes spp." ~ 0.18, # avg of small range rockfish
      species_name == "Sebastes spp. juv" ~ 0.18, # avg of small range rockfish
      species_name == "Sebastes pinniger" ~ 700, # up to 700
      
      # Sculpins, all guesses based on cabezon + red irish
      species_name == "Jordania zonope" ~ 0.05,
      species_name == "Artedius harringtoni" ~ 0.05,
      species_name == "Artedius lateralis" ~ 0.05,
      species_name == "Artedius fenestralis" ~ 0.05,
      species_name == "Hemilepidotus hemilepidotus" ~ 0.05,# actual estimate
      species_name == "Cottidae spp." ~ 0.05,
      species_name == "Enophrys bison" ~ 0.05,
      species_name == "Rhamphocottus richardsonii" ~ 0.05,
      species_name == "Scorpaenichthys marmoratus" ~ 0.05,# actual estimate
      species_name == "Oligocottus maculosus" ~ 0.05,
      species_name == "Leptocottus armatus" ~ 0.05,
      species_name == "Blepsias cirrhosus" ~0.05,
      species_name == "Myoxocephalus polyacanthocephalus" ~ 0.05,
      species_name == "Myoxocephalus aenaeus" ~ 0.05,
      species_name == "Asemichthys taylori" ~ 0.05,
      
      #Perch
      species_name == "Embiotoca lateralis" ~ 1, # inferred
      species_name == "Rhacochilus vacca" ~ 1, # unknown, inferred from other perch
      species_name == "Brachyistius frenatus" ~ 1, # inferred
      species_name == "Cymatogaster aggregata" ~ 1, # unknown, inferred from other perch
      species_name == "Embiotocidae spp." ~ 1, # unknown, inferred from other perch
      species_name == "Percidae spp." ~ 1, # unknown, inferred from other perch
      
      # Gunnels + gunnel-like fish
      species_name == "Anarrhichthys ocellatus" ~ 0.05,
      species_name == "Apodichthys flavidus" ~ 0.05, # inferred from wolf eel
      species_name == "Pholis ornata" ~ 0.05, # inferred from wolf eel
      species_name == "Pholis laeta" ~ 0.05, # inferred from wolf eel
      species_name == "Pholis clemensi" ~ 0.05, # inferred from wolf eel
      species_name == "Pholis spp." ~ 0.05, # inferred from wolf eel
      species_name == "Lumpenus sagitta" ~ 0.05, # inferred from wolf eel
      species_name == "Chirolophis nugator" ~ 0.05, # inferred from wolf eel
      
      #Misc
      species_name == "Liparis florae" ~ 0.05, # inferred from wolf eel
      species_name == "Aulorhynchus flavidus" ~ 1, # inferred from perch
      species_name == "Syngnathus leptorhynchus" ~ 1, # inferred from perch
      species_name == "Clupea pallasii" ~ 525, # 50 - 1000
      species_name == "Gasterosteus aculeatus" ~ 1, # inferred from perch
      species_name == "Porichthys notatus" ~ 300, # midshipman migrate at night 0 - 366 m deep
      species_name == "Gibbonsia metzi" ~ 1, # inferred from perch
      species_name == "Citharichthys stigmaeus" ~ 100, # flatfish seem to move a lot, avg from them
      TRUE ~ as.numeric(0.01)), # to give all the inverts a range weight of 1
      range_weight = 1/(range*100), # scale smallest range to a weight of 1, everything else is fraction
      weight_weighted = weight_size_class_sum*range_weight # weight weights by range
      )
  
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
    site_code== "BMSC6" &year =="2022" &depth == "8.5" &survey_depth== "6.5" ~ "no",
    site_code== "BMSC6" &year =="2022" &depth == "5.5" &survey_depth== "9" ~ "no",
    site_code== "BMSC6" &year =="2022" &depth == "6" &survey_depth== "9" ~ "no",
    site_code== "BMSC1" &year == "2021" & survey_depth == "6.5" ~ "no",
    site_code== "BMSC1" &year == "2022" & survey_depth == "4.7" ~ "no",
    site_code== "BMSC5" &year == "2022" & survey_depth == "6.2" ~ "no",
    site_code== "BMSC11" &year == "2022" & survey_depth == "8.5" ~ "no",
    site_code== "BMSC12" &year == "2022" & survey_depth == "9" ~ "no",
    site_code== "BMSC11" &year == "2023" & survey_depth == "5.5" ~ "no",
    site_code== "BMSC12" &year == "2023" & survey_depth == "6.5" ~ "no",
    site_code== "BMSC24" &year == "2023" & survey_depth == "7.5" ~ "no",
    site_code== "BMSC25" &year == "2023" & survey_depth == "5.5" ~ "no",
    site_code== "BMSC26" &year == "2023" & survey_depth == "9.5" ~ "no",
    site_code== "BMSC27" &year == "2023" & survey_depth == "7" ~ "no",
    site_code== "BMSC1" &year == "2023" & survey_depth == "8" ~ "no",
    site_code== "BMSC5" &year == "2023" & survey_depth == "8" ~ "no",
    site_code== "BMSC6" &year == "2023" & survey_depth == "5.5" ~ "no",
    site_code== "BMSC8" &year == "2023" & survey_depth == "7.5" ~ "no",
    TRUE ~ as.character("yes")
    # cut out the rls transects I didn't directly measure nh4 on 
  )) %>%
    filter(correct == "yes") %>%
    select(-c(correct, hour))
}


# Mapping -----

# Site map
old_site_map <- function(lat_min, lat_max, long_min, long_max, 
                     coord_data, map_data, 
                     add_points, add_annotate){

# make basic map plot
  gg <- ggplot() +
    geom_sf(data = {{map_data}}, fill = "white", colour = blue) +
    coord_sf(xlim = c({{lat_min}}, {{lat_max}}), ylim = c({{long_min}}, {{long_max}}), expand = FALSE) +
    # add aesthetic elements
    theme_black() +
    theme(panel.background = element_rect(fill = blue),
          panel.grid.major = element_line(color = blue)) +
    # axis things
    labs(x = "Longitude", y = "Latitude",
         fill = "Habitat") +
    scale_x_continuous(breaks = seq({{lat_min}}, {{lat_max}}, by = 0.1))
  
  # add bells and whistles
  if(add_points == TRUE){
    gg + geom_sf(data = {{coord_data}}, 
                 colour = "black",
                 fill = "black",
                 alpha = 0.5,
                 size = 8,
                 aes(pch = Habitat)) +
      scale_shape_manual(values = c(21, 25), drop = F) +
      guides(pch = guide_legend(override.aes = 
                                  list(colour = "black")))} 

  # add arrow + scale bar
  if(add_annotate == TRUE){
    gg +
    annotation_scale(location = "br", width_hint = 0.4) +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
}
}

# Site map from Nikola
site_map <- function(lat_min,
                     lat_max,
                     long_min,
                     long_max,
                     coord_data,
                     map_data,
                     add_points,
                     add_annotate) {
  # make basic map plot
  plot <- ggplot() +
    geom_sf(data = map_data,
            fill = "white",
            colour = blue) +
    coord_sf(
      xlim = c(lat_min, lat_max),
      ylim = c(long_min, long_max),
      expand = FALSE
    ) +
    # add aesthetic elements
    theme_black() +
    theme(
      panel.background = element_rect(fill = blue),
      panel.grid.major = element_line(color = blue)
    ) +
    # axis things
    labs(x = "Longitude", y = "Latitude",
         fill = "Habitat") +
    scale_x_continuous(breaks = seq(lat_min, lat_max, by = 0.1))
  
  if (add_points) {
    plot <- plot + geom_sf(
      data = coord_data,
      colour = "black",
      fill = "black",
      alpha = 0.5,
      size = 8,
      aes(pch = Habitat)
    ) +
      scale_shape_manual(values = c(21, 25), drop = F) +
      guides(pch = guide_legend(override.aes =
                                  list(colour = "black")))
  }
  if (add_annotate) {
    plot <- plot + annotation_scale(location = "br", width_hint = 0.4) +
      annotation_north_arrow(
        location = "br",
        which_north = "true",
        pad_x = unit(0.0, "in"),
        pad_y = unit(0.2, "in"),
        style = north_arrow_fancy_orienteering
      )
  }
  print(plot)
}


# big inset map
inset_map <- function(rect_xmin, rect_xmax, rect_ymin, rect_ymax, 
                      map_data){
  
  ggplot() +
    geom_sf(data = {{map_data}}, fill = "white", colour = blue) +
    coord_sf(xlim = c(-128.5, -123), ylim = c(48.25, 51), expand = FALSE) +
    # add rectangle for zoomed in part
    geom_rect(aes(xmin = {{rect_xmin}}, xmax = {{rect_xmax}}, ymin = {{rect_ymin}}, ymax = {{rect_ymax}}), color = "red", fill = NA, inherit.aes = FALSE) +
    # add aesthetic elements
    theme_bw() +
    theme(panel.background = element_rect(fill = blue),
          panel.grid.major = element_line(color = blue),
          panel.border = element_rect(fill = NA, colour = "black"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(0, "pt"),
          plot.title = NULL,
          plot.margin=grid::unit(c(0,0,0,0), "mm"))
}

# add label for Vancouver
# geom_text(aes(x = -123.120694, y = 49.282694, 
#               label = "Vancouver",
#               fontface = "bold")) 



map_daddy <- function(lat_min, lat_max, long_min, long_max, 
                      coord_data, nh4_var, kelp_var, point_size, map_file, invert) {
  
  if(invert == FALSE){
    sea <- blue
    land <- "white"
  }
  if(invert == TRUE){
    sea <- "white"
    land <- blue
  }
    
  ggplot() +
    geom_sf(data = map_file, fill = land, colour = sea) +
    # add points
    geom_sf(data = coord_data, 
            colour = "black",
            alpha = 0.9,
            size = point_size,
            aes(fill = {{nh4_var}},
                pch = {{kelp_var}})) +
    viridis::scale_fill_viridis(option="magma", direction = -1,
                                limits = c(0, 2),
                                guide = guide_colorbar(frame.colour = "white", ticks.colour = "white")) +
    coord_sf(xlim = c(lat_min, lat_max), ylim = c(long_min, long_max), expand = FALSE)  +
    labs(fill = expression(paste("NH"[4]^" +",(mu*M)))) +
    scale_shape_manual(values = c(21, 25), drop = F) +
    guides(pch = guide_legend(override.aes = 
                                list(colour = "white"))) +
    # Themes
    theme_bw() +
    theme(
          # panel stuff
          panel.background = element_rect(fill = sea),
          panel.grid.major = element_line(color = sea),
          panel.border = element_rect(fill = NA, colour = "black"),
          # remove axis
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(0, "pt"),
          plot.title = NULL,
          plot.margin=grid::unit(c(0,0,0,0), "mm"),
          # Specify legend options
          legend.background = element_rect(color = NA, fill = "black"),  
          legend.key = element_rect(color = "black",  fill = "black"),  
          legend.key.size = unit(1.2, "lines"),  
          legend.key.height = NULL,  
          legend.key.width = NULL,      
          legend.text = element_text(size = 24, color = "white"),  
          legend.title = element_text(size = 24, face = "bold", hjust = 0, color = "white"),  
          legend.position = "right",  
          legend.text.align = NULL,  
          legend.title.align = NULL,  
          legend.direction = "vertical",  
          legend.box = NULL,
          # Specify facetting options
          strip.background = element_rect(fill = "grey30", color = "grey10"),  
          strip.text.x = element_text(size = 30, color = "white"),  
          strip.text.y = element_text(size = 30, color = "white",angle = -90),  
          # Specify plot options
          plot.background = element_rect(color = "black", fill = "black"),  
          ) +
    annotation_scale(location = "br", width_hint = 0.4) +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
}


map_daddy_np <- function(lat_min, lat_max, long_min, long_max, 
                      map_file, invert) {
  
  if(invert == FALSE){
    sea <- blue
    land <- "white"
  }
  if(invert == TRUE){
    sea <- "white"
    land <- blue
  }
  
  ggplot() +
    geom_sf(data = map_file, fill = land, colour = sea) +
    coord_sf(xlim = c(lat_min, lat_max), ylim = c(long_min, long_max), expand = FALSE)  +
    # Themes
    theme_bw() +
    theme(
      # panel stuff
      panel.background = element_rect(fill = sea),
      panel.grid.major = element_line(color = sea),
      panel.border = element_rect(fill = NA, colour = "black"),
      # remove axis
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      plot.title = NULL,
      plot.margin=grid::unit(c(0,0,0,0), "mm"),
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "black",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = 20, color = "white"),  
      legend.title = element_text(size = 20, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL,
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = 30, color = "white"),  
      strip.text.y = element_text(size = 30, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
    ) +
    annotation_scale(location = "br", width_hint = 0.4) +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
}



# Plotting -----

# Dot Whisker Plot -----

dot_whisker <- function(sum_data, all_data, x_var, y_var){
  ggplot() +
    geom_point(data = {{sum_data}},
               aes(x = {{x_var}}, y = {{y_var}}, colour = {{x_var}}),
               size = 6) +
    geom_errorbar(data = {{sum_data}},
                  aes(x = {{x_var}},
                      y = {{y_var}},
                      ymin = conf.low,
                      ymax = conf.high, 
                      colour = {{x_var}}),
                  width = 0.4,
                  linewidth = 1.5) +
    geom_jitter(data = {{all_data}}, 
                aes(x = {{x_var}}, y = {{y_var}}, colour = {{x_var}}), 
                size = 3, alpha = 0.5, height=0) +
    theme_black() + 
    theme(legend.position = "none",
          plot.title = element_text(size = 30)) +
    scale_colour_manual(values = rev(csee_pal)) +
    labs(y = expression(paste("Ammonium"~(mu*M)))) +
    ylim(c(0, 3.8))
}


# Standardize variables -----

scale_vars <- function(datafile){
  {{datafile}} %>%
    mutate(
      # kelp forest level variables
      forest_biomass = BiomassM*Area_m2,
      den_scale = c(scale(DensityM)), # make sure density is right
      kelp_bio_scale = c(scale(BiomassM)),
      log_kelp = log(BiomassM + 0.001),
      log_kelp_scale = c(scale(log_kelp)),
      forest_bio_scale = c(scale(forest_biomass)),
      area_scale = c(scale(Area_m2)),
      # transect level variables
      bio_tran_scale = c(scale(Biomassm2kg)),
      log_kelp_tran = log(Biomassm2kg + 0.001),
      den_tran_scale = c(scale(kelp_den)), # make sure density is right
      # log the pee diff?
      log_pee_diff = log(in_minus_out + 1),
      # the biomass + abundance variables
      weight_sum_scale = c(scale(weight_sum)),
      all_weighted_scale = c(scale(all_weight_weighted)),
      abundance_scale = c(scale(abundance)),
      # biodiversity variables
      rich_scale = c(scale(species_richness)),
      shannon_scale = c(scale(shannon)),
      simpson_scale = c(scale(simpson)),
      # abiotic variables I should control for
      depth_scale = c(scale(depth_avg)),
      tide_scale = c(scale(avg_exchange_rate)),
      tide_cat = factor(as.factor(ifelse(avg_exchange_rate < -0.1897325, "Ebb",
                                         ifelse(avg_exchange_rate < 0.1897325, "Slack", "Flood"))),
                        levels = c("Ebb", "Slack", "Flood")),
      # only slack and flood
    ) 
}
