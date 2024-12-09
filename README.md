# Ch3_Spatial_pee
Everything for my spatial ammonium variation chapter

## Code
Analyses for each of the three main sections of this paper:
1) RLS_blitz_analysis.R contains code for the among site, meso-scale experiment
2) Kelp_pee_analysis.R contains code for the within site, small-scale experiment
3) Cage_pee_analysis.R contains code for the within-meters, smallest-scale experiment


Preliminary calculations:

Prelim_calculations folder contains an .R script for each of the above sections where I calculate ammonium concentrations from the seawater samples which we analysed using fluorometry. The final dataframes from these scripts are saved to the Output/Output_data folder. I also have an .R script to calculate kelp forest metrics for use in the Kelp_pee_analysis.R script.

Functions.R stores all of the functions I've written, this file is called at the beginning of most other scripts

Maps.R contains code to generate the map in Figure 1 and additional maps for presentations.

## Data
Data is subsectioned into logical folders
1) Cage_experiment contains the raw cage data and fluorometry data from the cuke cage experiments, and the unpublished cage data from the crab caging experiment which was initially cleaned in the repo for my crab chapter. This folder also contains the two excretion data .csv's (all_cuke_excretion.csv and red_crab_excretion.csv) generated in other projects which I use here to calculate nh4+ supply rates from the cages
2) Emily_flow contains all the raw data from Emily's flow rate experiments, which didn't make it into the final manuscript but is shown in the Supplement
3) Hakai_coast contains the shapefiles needed to make the map in Fig. 1a
4) RLS contains three folders (2021, 2022, 2023) which each contain all of the raw nh4 data from each sampling year. RLS_data contains the data for the biological communities surveyed, the true site coordinates (true_coords.csv), and a folder of outdated RLS data which I don't use. Tides_1_min.csv is the downloaded tide data used to calculate tide height rate of change
5) Size_data is a folder of published size data for invertebrates that I used to calculate mean wet weights for certain taxa
6) Team_kelp contains all of the raw ammonium data (one .csv per site, labelled with the date and site code), downloaded tide height data, the RLS survey data, the kelp forest metric data, and a folder (Output_data) that contains slightly processed kelp forest metric data.


## Datasheets
Self explanatory! Datasheets!

## Images
Images of a red rock crab and Cali cuke used in Fig. 4

## Output

Output_data has all of the data that's created in scripts (eg the pee calculations) and stored here so I can load it into an analysis script without going through all of the cleaning steps each time. This is also where dataframes I want to export to .csv to include as tables live

Pub_figs contains all the figures in the final paper

Pres_figs contains figures for presentations

Figures are just the misc and intermediate figures I've made over the years

## Protocols
Self explanatory

## renv
Package manager library

## Writing
Drafts of this manuscript!
