# Try using random forest??????
# Em Lim
# Jan 2024


# Fuck me I guess

# Load packages
library(tidyverse)
library(randomForest)

# Load data from RLS Blitz analysis
# rls_final is the df I've been analysing, with a row for each survey
rls <- read_csv("Output/Output_data/rls_final.csv", show_col_types = FALSE) %>%
  # rls_wider is just the species data. one column per species, one row per survey
  left_join(read_csv("Output/Output_data/rls_wider.csv", show_col_types = FALSE), by = "survey_id") %>%
  as.data.frame() %>%
  mutate(site_code = as.factor(site_code))

# add underscores to the stupid latin names
colnames(rls) <- colnames(rls) %>% 
  str_replace_all(" ","_")


# Load just the rls data so I can play with taxonomic level
rls_sp <- read_csv("Output/Output_data/rls_species.csv", show_col_types = FALSE) %>%
  select(-c(FID, 3:18, 29:38)) %>% # cut extra columns
  clean_phylo_names() %>% # function to fix naming errors
  mutate(match = reporting_name == species_name) %>% 
  separate(species_name, c("genus", "species")) %>% # the additional pieced discarded just means it dropped the period for all the spp. species
  mutate(species = ifelse(species == "spp" | species == "shrimp", NA, species),
         genus = ifelse(genus == "Cottidae" |
                          genus == "Paguroidea" |
                          genus == "Unidentified" |
                          genus ==  "Brachyura" |
                          genus == "Actinopterygii" |
                          genus == "Embiotocidae" |
                          genus == "Percidae", NA, genus)
         )

# use this to check the species with mismatching reporting vs scientific names
# false <- rls_sp %>%
#   filter(match == FALSE) %>%
#   select(-c(1:3, 11:13)) %>%
#   unique()

# try grouping by species and see how many we get
# omit indivs not identified to species? na.omit(species) ?
species <- rls_sp %>%
  group_by(genus, species) %>%
  summarise()
# 110 groups

# try grouping by genus 
genus <- rls_sp %>%
  group_by(genus) %>%
  summarise(n_species = n_distinct(species))
# 95 groups

# try grouping by family 
family <- rls_sp %>%
  group_by(family) %>%
  summarise(n_genus = n_distinct(genus),
            n_species = n_distinct(species))
# 65 groups

# try grouping by order 
order <- rls_sp %>%
  group_by(order) %>%
  summarise(n_family = n_distinct(family),
            n_genus = n_distinct(genus),
            n_species = n_distinct(species))
# 31 groups

# try grouping by class 
class <- rls_sp %>%
  group_by(class) %>%
  summarise(n_order = n_distinct(order),
            n_family = n_distinct(family),
            n_genus = n_distinct(genus),
            n_species = n_distinct(species))
# 12 groups

# Try to make a random forest I guess????

# let's make a simple one first
(RF.model = randomForest(nh4_avg ~ abundance_stand + tide_stand + shannon_stand + depth_avg_stand, data = rls))

plot(RF.model)
randomForest::importance(RF.model)
varImpPlot(RF.model)
partialPlot(RF.model, rls, "abundance_stand")

# try to get all species names so i can include them without typing each out?
species_cols <- paste(colnames(rls[30:160]), collapse = " + ")

# try to make the formula
my_formula <- as.formula(paste("nh4_avg ~ abundance_stand + tide_stand + shannon_stand + depth_avg_stand +", paste(colnames(rls[30:160]), collapse = " + ")))

my_formula <- as.formula(paste("nh4_avg ~ ", paste(colnames(rls[30:160]), collapse = " + ")))

# try to use that formula
(RF.model = randomForest(my_formula, data = rls))

# evaluate that?
plot(RF.model)
randomForest::importance(RF.model)
varImpPlot(RF.model)
partialPlot(RF.model, rls, "abundance_stand")



# I can also subset df so i can just call all columns
rls_reduced <- rls %>%
  select(c(nh4_avg, abundance_stand, tide_stand, shannon_stand, depth_avg_stand, 30:160))

(RF.model = randomForest(nh4_avg ~ ., data = rls_reduced))
