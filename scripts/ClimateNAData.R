### Title: Data exploration of ClimateNA env variables
### Author: Alyssa Phillips
### Date: 3/24/22

# ClimateNA: https://climatena.ca/

# Input data
loc <- read.csv("~/Andropogon/ClimateNA data/CG_all_ClimateData_03242022.csv")
str(loc)

loc$period <- as.factor(loc$period)
loc$Pop <- as.factor(loc$Pop)
loc$Subpop <- as.factor(loc$Subpop)
