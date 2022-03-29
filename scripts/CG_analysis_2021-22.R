### Title: Common garden analysis
### Author: Alyssa Phillips
### Date: 2/17/22

library(googlesheets4)
library(ggplot2)

###
### Load datasets ----
###

# > biomass ----
# biomass <- googlesheets4::read_sheet(
#   "https://docs.google.com/spreadsheets/d/1itGdDp74-ZyEBcsbkfP1GFyH81npUqBm-4a75GpCtkE/edit?usp=sharing", 
#   sheet = "biomass", 
#   skip = 2
# )
biomass <- read.csv("~/Andropogon/Common garden/raw data/2021/2021 phenotype data - biomass_d03102022.csv", skip = 2, header = F)
colnames(biomass) <- c("position", "block", "population", "genotype","biomass_g")

biomass$fposition <- as.factor(biomass$position)
biomass$fpopulation <- as.factor(biomass$population)
biomass$fgenotype <- as.factor(biomass$genotype)
biomass$fblock <- as.factor(biomass$block)

str(biomass)
head(biomass)

###
### Data exploration ----
###

hist(biomass$biomass_g, breaks = 14) # spread of the data
dotchart(biomass$biomass_g) # y axis is order in the dataset, another way to show spread
dotchart(biomass$biomass_g, color = biomass$fblock) # can color by category to look for clusters/patterns
dotchart(biomass$biomass_g, color = biomass$fpopulation)

#* TO DO ----
# test assumptions of env variables of interest & check what is confounded
# pairs(loyn[,c("L.AREA", "L.DIST", "ABUND", "GRAZE", "L.LDIST", "YR.ISOL", "ALT")])
# hist(loyn$AREA)
# dotchart(loyn$AREA)

# what the spread across populations? (variance homogeneity?)
boxplot(biomass$biomass_g ~ biomass$population) # not even 
boxplot(biomass$biomass_g ~ biomass$block) # pretty even variance

# is there equal spread of area across all levels of pop?
ggplot(biomass, aes(x = fpopulation, y = biomass_g)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  facet_grid(~fblock)

# ggplot(biomass, aes(x = fpopulation, y = biomass_g)) + 
#   geom_point() +
#   stat_smooth(method = "lm")
