############################
###Distance decay - genus###
############################

setwd("/home/gilda/all")
load(file = "all_dataset_ps_srs_bacteria_hell.RData")

library(vegan)
library(geosphere)
library(phyloseq)

#https://jkzorz.github.io/2019/07/08/mantel-test.html

metadata <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)[samples,]
metadata$rasValue_swe[is.na(metadata$rasValue_swe)] <- 0

community <- t(genera_proportions_tab_hellinger[-1,])[samples,]
dim(community)
#[1]  988 1445
community_tran <- community[rowSums(community)>0,]
dim(community_tran)
#[1]  988 1412

#genus taxonomy
genera_counts_tab_trim <- data.frame(t(genera_proportions_tab_hellinger[-1,row.names(community_tran)]))
dim(genera_counts_tab_trim)
#[1]  988 1445

#geographic coordinates
lat <- as.numeric(metadata[,1])
lon <- as.numeric(metadata[,2])
geo_formatted <- cbind(lon,lat)
row.names(geo_formatted) <- row.names(metadata)
dim(geo_formatted)
#[1] 988   2

#bioclimatic data
biocl_formatted <- metadata[,c(3:12,14)]
dim(biocl_formatted)
#[1] 988   11

################
#Entire dataset#
################

#geographical distance
dist_taxa <- vegdist(genera_counts_tab_trim, method="bray") #remove unclassified
d_geo <- distm(geo_formatted, fun = distHaversine) #this function works with decimal degrees
dist_geo = as.dist(d_geo)
dc_geo_genera <- mantel(dist_taxa, dist_geo, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_geo_genera
# Mantel statistic based on Spearman's rank correlation rho
# 
# Call:
# mantel(xdis = dist_taxa, ydis = dist_geo, method = "spearman",      permutations = 1000, na.rm = TRUE)
# 
# Mantel statistic r: 0.1499
#       Significance: 0.000999
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99%
# 0.0180 0.0230 0.0275 0.0308
# Permutation: free
# Number of permutations: 1000

#bioclimatic variables
dist_biocl <- vegdist(biocl_formatted, method="eucl")
dc_biocl_genera <- mantel(dist_taxa, dist_biocl, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_biocl_genera
# Mantel statistic based on Spearman's rank correlation rho
# 
# Call:
# mantel(xdis = dist_taxa, ydis = dist_biocl, method = "spearman", permutations = 1000, na.rm = TRUE)
# 
# Mantel statistic r: 0.1025
#       Significance: 0.000999
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99%
# 0.0248 0.0321 0.0405 0.0511
# Permutation: free
# Number of permutations: 1000

################
#Island dataset#
################

#get only islands
metadata_isl <- metadata[metadata[,"island"] == "y",]
dim(metadata_isl)
#[1] 142  20

#get community information - only islands
community_isl0 <- genera_counts_tab_trim[row.names(metadata_isl),]
dim(community_isl0)
#[1]  142 1445
community_isl <- community_isl0[rowSums(community_isl0)>0,]
dim(community_isl)
#[1]  142 1445

#geo and biocl
geo_formatted_isl <- geo_formatted[row.names(metadata_isl),]
dim(geo_formatted_isl)
#[1] 142  2
biocl_formatted_isl <- biocl_formatted[row.names(metadata_isl),]
dim(biocl_formatted_isl)
#[1] 142  11

#geographical distance
dist_taxa_isl <- vegdist(community_isl, method="bray") #remove unclassified
d_geo <- distm(geo_formatted_isl, fun = distHaversine) #this function works with decimal degrees
dist_geo_isl = as.dist(d_geo)
dc_geo_genera_isl <- mantel(dist_taxa_isl, dist_geo_isl, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_geo_genera_isl
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = dist_taxa_isl, ydis = dist_geo_isl, method = "spearman",      permutations = 1000, na.rm = TRUE) 
# 
# Mantel statistic r: 0.3776 
#       Significance: 0.000999 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0644 0.0829 0.0956 0.1092 
# Permutation: free
# Number of permutations: 1000

#bioclimatic variables
dist_biocl_isl <- vegdist(biocl_formatted_isl, method="eucl")
dc_biocl_genera_isl <- mantel(dist_taxa_isl, dist_biocl_isl, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_biocl_genera_isl
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = dist_taxa_isl, ydis = dist_biocl_isl, method = "spearman",      permutations = 1000, na.rm = TRUE) 
# 
# Mantel statistic r: 0.3001 
#       Significance: 0.000999 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0597 0.0783 0.0913 0.1198 
# Permutation: free
# Number of permutations: 1000

##################
#Mainland dataset#
##################

#get only mainland
metadata_ml <- metadata[metadata[,"island"] == "n",]
dim(metadata_ml)
#[1] 846  20

#get community information - only mainland
community_ml0 <-genera_counts_tab_trim[row.names(metadata_ml),]
dim(community_ml0)
#[1]  846 1445
community_ml <- community_ml0[rowSums(community_ml0)>0,]
dim(community_ml)
#[1]  846 1445

#geo and biocl
geo_formatted_ml <- geo_formatted[row.names(metadata_ml),]
dim(geo_formatted_ml)
#[1] 846  2
biocl_formatted_ml <- biocl_formatted[row.names(metadata_ml),]
dim(biocl_formatted_ml)
#[1] 846  11

#geographical distance
dist_taxa_ml <- vegdist(community_ml, method="bray") #remove unclassified
d_geo <- distm(geo_formatted_ml, fun = distHaversine) #this function works with decimal degrees
dist_geo_ml = as.dist(d_geo)
dc_geo_genera_ml <- mantel(dist_taxa_ml, dist_geo_ml, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_geo_genera_ml
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = dist_taxa_ml, ydis = dist_geo_ml, method = "spearman",      permutations = 1000, na.rm = TRUE) 
# 
# Mantel statistic r: 0.1611 
#       Significance: 0.000999 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0133 0.0172 0.0209 0.0240 
# Permutation: free
# Number of permutations: 1000

#bioclimatic variables
dist_biocl_ml <- vegdist(biocl_formatted_ml, method="eucl")
dc_biocl_genera_ml <- mantel(dist_taxa_ml, dist_biocl_ml, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_biocl_genera_ml
hy6# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = dist_taxa_ml, ydis = dist_biocl_ml, method = "spearman",      permutations = 1000, na.rm = TRUE) 
# 
# Mantel statistic r: 0.1102 
#       Significance: 0.000999 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0239 0.0313 0.0374 0.0422 
# Permutation: free
# Number of permutations: 1000

###############################################
###Mainland dataset - elevation and distance###
###############################################

#elevation_distance
metadata_el_dist <- read.csv("distance_altitude.csv", row.names = 1)[samples,]
dim(metadata_el_dist)
#[1] 988   5

metadata_el_dist_ml <- metadata_el_dist[row.names(metadata_ml),]
dim(metadata_el_dist_ml)
#[1] 846   5

#correlations
dist_el <- vegdist(metadata_el_dist_ml[,5], method="eucl")
dc_el_genera <- mantel(dist_taxa_ml, dist_el, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_el_genera
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = dist_taxa_ml, ydis = dist_el, method = "spearman",      permutations = 1000, na.rm = TRUE) 
# 
# Mantel statistic r: 0.07805 
#       Significance: 0.000999 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0134 0.0177 0.0209 0.0249 
# Permutation: free
# Number of permutations: 1000

dist_coast <- vegdist(metadata_el_dist_ml[,3], method="eucl")
dc_coast_genera <- mantel(dist_taxa_ml, dist_coast, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_coast_genera
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = dist_taxa, ydis = dist_coast, method = "spearman",      permutations = 1000, na.rm = TRUE) 
# 
# Mantel statistic r: 0.01824 
#       Significance: 0.15684 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0240 0.0305 0.0357 0.0409 
# Permutation: free
# Number of permutations: 1000

dist_ocean <- vegdist(metadata_el_dist_ml[,4], method="eucl")
dc_ocean_genera <- mantel(dist_taxa_ml, dist_ocean, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_ocean_genera
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = dist_taxa, ydis = dist_ocean, method = "spearman",      permutations = 1000, na.rm = TRUE) 
# 
# Mantel statistic r: 0.01556 
#       Significance: 0.16783 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0219 0.0283 0.0334 0.0384 
# Permutation: free
# Number of permutations: 1000

save.image(file = "all_dataset_ps_srs_mantel_bacteria.RData")
#load(file = "all_dataset_ps_srs_mantel_bacteria.RData")

#plot
par(mfrow=c(5,2))

plot(dist_geo, dist_taxa)
abline(lm(dist_taxa ~ dist_geo), col = "red", lwd=3)

plot(dist_biocl, dist_taxa)
abline(lm(dist_taxa ~ dist_biocl), col = "red", lwd=3)

plot(dist_geo_isl, dist_taxa_isl)
abline(lm(dist_taxa_isl ~ dist_geo_isl), col = "red", lwd=3)

plot(dist_biocl_isl, dist_taxa_isl)
abline(lm(dist_taxa_isl ~ dist_biocl_isl), col = "red", lwd=3)

plot(dist_geo_ml, dist_taxa_ml)
abline(lm(dist_taxa_ml ~ dist_geo_ml), col = "red", lwd=3)

plot(dist_biocl_ml, dist_taxa_ml)
abline(lm(dist_taxa_ml ~ dist_biocl_ml), col = "red", lwd=3)

plot(dist_el, dist_taxa_ml)
abline(lm(dist_taxa_ml ~ dist_el), col = "red", lwd=3)

plot(dist_coast, dist_taxa_ml)
abline(lm(dist_taxa_ml ~ dist_coast), col = "red", lwd=3)

plot(dist_ocean, dist_taxa_ml)
abline(lm(dist_taxa_ml ~ dist_ocean), col = "red", lwd=3)

#dist_dec.png: 900 x 1400
dev.off()

