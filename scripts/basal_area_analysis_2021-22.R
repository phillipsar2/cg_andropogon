### Title: Common garden analysis - general
### Author: Alyssa Phillips
### Date: 2/17/22

#library(googlesheets4)
library(ggplot2)

###
### Load datasets ----
###

# bA <- googlesheets4::read_sheet(
#   "https://docs.google.com/spreadsheets/d/1itGdDp74-ZyEBcsbkfP1GFyH81npUqBm-4a75GpCtkE/edit?usp=sharing", 
#   sheet = "bA", 
#   skip = 2
# )
bA <- read.csv("~/Andropogon/Common garden/raw data/2021/2021 phenotype data - basal aread03162022.csv", skip = 2, header = F)
colnames(bA) <- c("position", "block", "population", "genotype","diameter_vertical", "diameter_horizontal")

head(bA)
str(bA)

bA$fposition <- as.factor(bA$position)
bA$fpopulation <- as.factor(bA$population)
bA$fgenotype <- as.factor(bA$genotype)
bA$fblock <- as.factor(bA$block)

bA$ndiameter_horizontal <- as.numeric(bA$diameter_horizontal)
bA$ndiameter_vertical <- as.numeric(bA$diameter_vertical)

str(bA)
head(bA)

# > Calculate basal area ----
# Calculated assuming an ellipse A = pi*ab

bA$basal_area <- pi * bA$ndiameter_horizontal * bA$ndiameter_vertical
head(bA$basal_area)

###
### Data exploration ----
###

hist(bA$basal_area, breaks = 20) # spread of the data
dotchart(bA$basal_area) # y axis is order in the dataset, another way to show spread

'Not a great spread of data. Mostly clustered between 0 and 1000'

dotchart(bA$basal_area, color = bA$fblock) # can color by category to look for clusters/patterns
dotchart(bA$basal_area, color = bA$fpopulation)

#* TO DO ----
# test assumptions of env variables of interest & check what is confounded
# pairs(loyn[,c("L.AREA", "L.DIST", "ABUND", "GRAZE", "L.LDIST", "YR.ISOL", "ALT")])
# hist(loyn$AREA)
# dotchart(loyn$AREA)

# what the spread across populations? (variance homogeneity?)
boxplot(bA$basal_area ~ bA$fpopulation) # not even 
boxplot(bA$basal_area ~ bA$fblock) # pretty even variance

'Unequal variation across populations and some difference in variation across blocks but not very substantial'

# is there equal spread of basal area across all levels of pop?
ggplot(bA, aes(x = fpopulation, y = basal_area)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  facet_grid(~fblock)

ggplot(bA, aes(x = fpopulation, y = basal_area)) +
  geom_point() +
  stat_smooth(method = "lm")

###
### Model building ----
###


###
### Testing assumptions ----
###

### this actually comes after model building

# (1) Variance homogeneity  ----
# Want everything to look similar -- starry night
# what we see is actually a swoop (a trumpet plot) - very indicative of a log relationship
plot(mod.cat, 1) 


#2 normality of errors
# almost always see a little bit of tails so it's a judgment call (but really is least important)
plot(mod.cat, 2)
hist(resid(mod.cat))

#3 variance homogeneity
# look across all different predictors
# tells us what's causing the problem in graph #1
plot(resid(mod.cat) ~ clams$LENGTH) # want to see a starry night - clearly length is the problem
plot(resid(mod.cat) ~ clams$fMONTH) # want variance of boxplots to be the same
