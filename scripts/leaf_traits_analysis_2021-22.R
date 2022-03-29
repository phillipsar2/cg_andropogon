### Title: Common garden analysis
### Author: Alyssa Phillips
### Date: 2/17/22

#library(googlesheets4)
library(ggplot2)
library(dplyr)

###
### Load datasets ----
###

# > data ----
# data <- googlesheets4::read_sheet(
#   "https://docs.google.com/spreadsheets/d/1itGdDp74-ZyEBcsbkfP1GFyH81npUqBm-4a75GpCtkE/edit?usp=sharing", 
#   sheet = "data", 
#   skip = 2
# )
data <- read.csv("~/Andropogon/Common garden/raw data/2021/2021 phenotype data - leaf traits_d03172022.csv", skip = 1, header = F)
colnames(data) <- c("population", "genotype","leaf", "leaf_press","dry_mass_g","area_cm2", "sd_cm2", "width_cm", "length_cm", "leaf_thickness_mm")

str(data)

data$fpopulation <- as.factor(data$population)
data$fgenotype <- as.factor(data$genotype)
data$fleaf <- as.factor(data$leaf)
data$fleaf_press <- as.factor(data$leaf_press)

str(data)
head(data)
dim(data)

# > metadata ----
meta <- read.csv("~/Andropogon/Common garden/raw data/2021/2021 phenotype data - genotype metadata_d03172022.csv", header = F, skip = 1)
colnames(meta) <- c("position", "block", "population", "genotype")

meta[colnames(meta)] <- lapply(meta[colnames(meta)], factor)

str(meta)

# > conditionally add metadata to data ----
temp <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = T)
colnames(temp) <- c("position", "block")

for (i in 1:dim(data)[1]){
  temp[i, ] <- unlist( meta[ match(data$genotype[i], meta$genotype), 1:2], use.names = F)
}

# bind data.frames together
df <- cbind(temp, data)

df$block <- as.factor(df$block)
df$position <- as.factor(df$position)

# > calculate SLA (cm/g) ----
df$sla <- df$area_cm2 / df$dry_mass_g

###
### Data exploration ----
###

# > leaf thickness ----
df$leaf_thickness_mm[493] <- 0.09

hist(df$leaf_thickness_mm, breaks = 14) # spread of the data
dotchart(df$leaf_thickness_mm) 
dotchart(df$leaf_thickness_mm, color = df$block) 
dotchart(df$leaf_thickness_mm, color = df$fpopulation)

boxplot(df$leaf_thickness_mm ~ df$population) # equal spread across levels?
boxplot(df$leaf_thickness_mm ~ df$block) 
ggplot(df, aes(x = fpopulation, y = leaf_thickness_mm)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  facet_grid(~block)

'Overall distribution of data points looks like a nice spread; variance homogenity is likely not met'

# > leaf width ----
df$width_cm[1212] <- NA

hist(df$width_cm, breaks = 14) # spread of the data
dotchart(df$width_cm) 
dotchart(df$width_cm, color = df$block) 
dotchart(df$width_cm, color = df$fpopulation)

boxplot(df$width_cm ~ df$population) # equal spread across levels?
boxplot(df$width_cm ~ df$block) 
ggplot(df, aes(x = fpopulation, y = width_cm)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  facet_grid(~block)

# > SLA ----
outliers <- df[which(df$sla > 500),]
outliers %>% select(sla, area_cm2, dry_mass_g)
plot(outliers$dry_mass_g, outliers$area_cm2)

hist(df$sla, breaks = 15) # spread of the data
dotchart(df$sla) 
dotchart(df$sla, color = df$block) 
dotchart(df$sla, color = df$fpopulation)
dotchart(df$sla, color = df$fleaf_press)

# exclude outliers
sla_no_outs <- df[which(df$sla < 600),]

hist(sla_no_outs$sla, breaks = 30) # spread of the data
dotchart(sla_no_outs$sla) 
dotchart(sla_no_outs$sla, color = sla_no_outs$block) 
dotchart(sla_no_outs$sla, color = sla_no_outs$fpopulation)
dotchart(sla_no_outs$sla, color = sla_no_outs$fleaf_press)

boxplot(sla_no_outs$sla ~ sla_no_outs$population) # spread across levels
boxplot(sla_no_outs$sla ~ sla_no_outs$block) 

#* TO DO ----
# test assumptions of env variables of interest & check what is confounded
# pairs(loyn[,c("L.AREA", "L.DIST", "ABUND", "GRAZE", "L.LDIST", "YR.ISOL", "ALT")])
# hist(loyn$AREA)
# dotchart(loyn$AREA)
