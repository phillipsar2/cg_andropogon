### Title: Big bluestem kinship matrix
### Author: Alyssa Phillips
### Date: 1/3/2023

# Code adapted from script provided by Julianna Porter and Jeffrey Ross-Ibarra
# Create a kinship matrix using eqn 4.9 from Caballero (pg 74) - Van Raden (2008)

library(ggplot2)
library("pheatmap")
library(RColorBrewer)

# Load data ----

## Common garden files
# geno_tab <- read.table(gzfile("~/Andropogon/kinship/cg.andro.lowcov.nomiss.ibs.gz"), header= T) # no missing data
# geno_tab <- read.table(gzfile("~/Andropogon/pca/lowcov/cg.lowcov.50k.ibs.gz"), header= T) # 50k SNPs

## All Andropogon files
geno_tab <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.ibs", header = T)

head(geno_tab)
dim(geno_tab)
### rows = loci
### col = ind

# Convert dataframe to matrix ----
og_mat <- as.matrix(geno_tab[,5:dim(geno_tab)[2]])

# > Treat missing data ----
paste0("sites with missing data: ", sum(og_mat == -1))

og_mat[og_mat == -1] <- NA

## Sites with no missing data
# og_mat_complete <- og_mat[complete.cases(og_mat),]
# og_mat_complete[1:5,1:5]
# dim(og_mat_complete)
# paste0("sites with complete data: ", dim(og_mat_complete)[1])

# Calculate qs ----
'qs = average frequency of marker s for all individiuals (N) across populations'
# Calculate snp means accounting for missing data 
snp_means <- rowSums(og_mat, na.rm = T) / ( dim(og_mat)[2] - rowSums(is.na(og_mat)))
hist(snp_means)

# Calculate Kinship Matrix ----
# > Calculate sum of squares (denominator of fMij) ----
ssq = sum(snp_means * (1 - snp_means)) 

# > Calculate off-diagonal elements ----
n = dim(og_mat)[2] # number of genotypes

# > replace NAs with the SNP means ----
# calculate amount of missing data at each site
miss_sites <- rowSums(is.na(og_mat)) / n
hist(miss_sites)

miss_geno <- colSums(is.na(og_mat)) / dim(og_mat)[1]
hist(miss_geno)

# replace
for (i in 1:dim(og_mat)[1]){
  og_mat[i, is.na(og_mat[i,])] <- snp_means[i]
}

sum(is.na(og_mat))
og_mat[1:10,1:5]

# calculate molecular coancestry (fMij) between individuals i and j
kin <- matrix(nrow =  n, ncol = n)
for (ind_i in 1:n){
  for (ind_j in 1:n){
    # f <- sum( (og_mat_complete[,ind_i] - snp_means) * (og_mat_complete[,ind_j] - snp_means) ) / ssq
    f <- sum( (og_mat[,ind_i] - snp_means) * (og_mat[,ind_j] - snp_means) ) / ssq
    kin[ind_i, ind_j] <- f
  }
}

kin[1:5,1:5]

# > Diagonal ----
# # (A) Center mean on 1
# diag_name = "centereddiag"
# fin_kin <-  kin / mean(diag(kin))
# dim(fin_kin)
# fin_kin[1:5,1:5]

# # (B) Set diagonal to NA
diag_name = "NAdiag"
fin_kin <- kin
diag(fin_kin) <- NA

# (C) Calculate diagonal from two read draws from 100k sites
# diag_name = "inbreddiag"
# run1 <- read.table("~/Andropogon/kinship/all.andro.lowcov.miss20.min2.100k.run1.ibs.txt")
# run2 <- read.table("~/Andropogon/kinship/all.andro.lowcov.miss20.min2.100k.run2.ibs.txt")

# run1_mat <- as.matrix(run1[,5:dim(run1)[2]])
# run2_mat <- as.matrix(run2[,5:dim(run2)[2]])

# run1_mat[run1_mat == -1] <- NA
# run2_mat[run2_mat == -1] <- NA
# 
# s = run1_mat + run2_mat
# in_af <- s / 2
# in_af[1:5,1:5]
# hist(in_af)

## % sites that are not identical by state
# nibs <- colSums(s == 1, na.rm = T) / (dim(s)[1] - colSums(is.na(s)))

# inbred <- colMeans(in_af, na.rm = T) # 
# hist(inbred)

# fin_kin <- kin
# diag(fin_kin) <- inbred

## missing data versus inbreeding coef
# plot(colSums(is.na(in_af)) / 100000, inbred,
     # xlab = "Missing data per genotype",
     # ylab = "inbreeding coef")
# cor(colSums(is.na(in_af)) / 100000, inbred)

# > Set values < 0 to 0 ----
fin_kin[fin_kin < 0] <- 0

# Is mean associated with missing data?
# plot(miss_geno, diag(fin_kin))
# cor(miss_geno, diag(fin_kin))

# Clean up kinship mat ----
# > Load metadata ----
# meta <- read.csv("~/Andropogon/pca/lowcov/cg.geno_meta_02152023.csv", header = T)
meta <- read.csv("~/Andropogon/pca/lowcov/geno_meta_02262023.csv", header = T)
head(meta)
dim(meta)
str(meta)

meta$ploidy <- as.factor(meta$ploidy)
meta$population <- as.factor(meta$population)
# meta$geno_short <- sapply(meta$genotype, strsplit(sapply(strsplit(as.character(meta$genotype), split = "_"), `[`, 1)
meta$pop_ploidy <- paste0(meta$population, "_", meta$ploidy)
# meta$sub_population <- as.factor(meta$sub_population)

# Add genotypes to columns and rows
colnames(fin_kin) = as.factor(meta$pop_ploidy)
rownames(fin_kin) = as.factor(meta$pop_ploidy)

# Plot kinship mat ----
# pdf(file = paste0("~/Andropogon/kinship/cg.50k.kinship_mat.",Sys.Date(),".pdf"),
#     width = 10,
#     height = 9)

# pos_df = data.frame("Pos", meta$population)
p <- pheatmap(fin_kin, 
         fontsize = 6
         # annotation_row = pos_df
         )
p
# dev.off()

# ggsave(paste0("~/Andropogon/kinship/cg.50k.kinship_mat.",Sys.Date(),".jpeg"), plot = p )
ggsave(paste0("~/Andropogon/kinship/all.andro.100k.kinship_mat.",Sys.Date(),".jpeg"),
       plot = p,
       height = 15,
       width = 15,
       unit = "in")
# add in ploidy info

## >  W Subset plot ----

west <- c("Boulder_CO_3","Boulder_CO_2","ESL_2","Konza_3","Konza_2","NM_2","Austin_TX_3", "SAL_2","CDB_2","WEB_2","REL_2","Austin_TX_2","BB_2","Afton_TX_2", "NE_2","SAL_3")

fin_kin_sub <- fin_kin[colnames(fin_kin) %in% west,rownames(fin_kin) %in% west]

p <- pheatmap(fin_kin_sub, 
              fontsize = 10)
p

ggsave(paste0("~/Andropogon/kinship/west.andro.100k.kinship_mat.",Sys.Date(),".jpeg"),
       plot = p,
       height = 15,
       width = 15,
       unit = "in")

# pl_color <- as.matrix(x = NA, nrow = length(meta$ploidy), ncol = 1)
# pop_color <- as.matrix(x = NA, nrow = length(meta$population), ncol = 1)
# 
# for (i in 1:length(meta$ploidy)){
#   if (meta$ploidy[i] == 2){
#     pl_color[i] <- "white"
#   }
#   else{
#     pl_color[i] <- "black"
#   }
# }
# 
# 
# heatmap(fin_kin,
#         RowSideColors = pl_color,
#         # ColSideColors = pop_color,
#         xlab = "population"
# )

# Export kinship matrix ----
write.table(fin_kin, paste0("~/Andropogon/kinship/cg.50k.kinship_mat.", diag_name,".",Sys.Date(),".txt"), row.names = T, col.names = T, quote = F)

