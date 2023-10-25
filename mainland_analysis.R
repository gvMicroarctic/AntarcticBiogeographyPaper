##################################
###analyses on mainland samples###
##################################

setwd("/home/gilda/all")
load(file = "all_dataset_ps_srs_bacteria_hell.RData")

################################
###add elevation and distance###
################################

#import community information
community <- t(genera_proportions_tab_hellinger[-1,samples])
dim(community)
#[1]  988 1445

#import metadata
metadata <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)[samples,]
dim(metadata)
#[1] 988  20

#elevation_distance
metadata_el_dist <- read.csv("distance_altitude.csv", row.names = 1)[samples,]
dim(metadata_el_dist)
#[1] 988   5

#get only mainland
metadata_ml <- metadata[metadata[,"island"] == "n",]
dim(metadata_ml)
#[1] 846  20

metadata_el_dist_ml <- metadata_el_dist[row.names(metadata_ml),]
dim(metadata_el_dist_ml)
#[1] 846   5

metadata_all <- cbind(metadata_ml[,c(3:12,14)], metadata_el_dist_ml[,c(3:5)])
metadata_all$rasValue_swe[is.na(metadata_all$rasValue_swe)] <- 0

#get community information - only islands
community_ml0 <- community[row.names(metadata_ml),]
dim(community_ml0)
#[1]  846 1445
community_ml <- community_ml0[rowSums(community_ml0)>0,]
dim(community_ml)
#[1]  846 1445

#correlations
library("vegan")
dist_taxa <- vegdist(community_ml, method="bray")
dist_coast <- vegdist(metadata_el_dist_ml[,3], method="eucl")
dc_coast_genera <- mantel(dist_taxa, dist_coast, method = "spearman", permutations = 1000, na.rm = TRUE)
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
dc_ocean_genera <- mantel(dist_taxa, dist_ocean, method = "spearman", permutations = 1000, na.rm = TRUE)
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

dist_el <- vegdist(metadata_el_dist_ml[,5], method="eucl")
dc_el_genera <- mantel(dist_taxa, dist_el, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_el_genera
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = dist_taxa, ydis = dist_el, method = "spearman",      permutations = 1000, na.rm = TRUE) 
# 
# Mantel statistic r: 0.07805 
#       Significance: 0.000999 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0133 0.0176 0.0216 0.0255 
# Permutation: free
# Number of permutations: 1000

###########
###dbRDA###
###########

#hellinger transformed genera dataset
community_ml

#standardize environmental variables
library(vegan)
env.z <- decostand(metadata_all, method = "standardize")
round(apply(env.z, 2, mean), 1)
apply(env.z, 2, sd)
env.z <- data.frame(env.z)

#run RDA
dbrda.upr <- capscale(community_ml ~ ., data = data.frame(env.z), dist="bray") # upper model limit (the "full" model; it takes in account all variables)
dbrda.lwr <- capscale(community_ml ~ 1, data = data.frame(env.z), dist="bray") # lower model limit

#run the ordistep model to identify variables with explanatory power
fwd.sel <- ordiR2step(dbrda.lwr, #lower model limit
                      scope = formula(dbrda.upr), #upper model limit
                      direction = "forward",
                      R2scope = TRUE, #can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) #change to TRUE to see the selection process!

#check the new model with forward-selected variables
fwd.sel$call
# capscale(formula = community_ml ~ rasValue_bio10 + rasValue_bio5 +
#            rasValue_bio4 + rasValue_bio14 + rasValue_bio18 + rasValue_bio1 +
#            rasValue_bio2 + rasValue_bio17 + dist_ocean + rasValue_swe +
#            rasValue_bio15 + dist_coast + rasValue_bio12, data = data.frame(env.z),
#          distance = "bray")

#get statistical significance of the model
anova.cca(fwd.sel, permutations = 1000) #same results as anova()
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 1000
# 
# Model: capscale(formula = community_ml ~ rasValue_bio10 + rasValue_bio5 + rasValue_bio4 + rasValue_bio14 + rasValue_bio18 + rasValue_bio1 + rasValue_bio2 + rasValue_bio17 + dist_ocean + rasValue_swe + rasValue_bio15 + dist_coast + rasValue_bio12, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs      F   Pr(>F)
# Model     13   55.877 15.286 0.000999 ***
#   Residual 832  233.950
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#get statistical significance of the model for each variable
anova.cca(fwd.sel, step = 1000, by = "term")
# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = community_ml ~ rasValue_bio10 + rasValue_bio5 + rasValue_bio4 + rasValue_bio14 + rasValue_bio18 + rasValue_bio1 + rasValue_bio2 + rasValue_bio17 + dist_ocean + rasValue_swe + rasValue_bio15 + dist_coast + rasValue_bio12, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs       F Pr(>F)
# rasValue_bio10   1    9.771 34.7476  0.001 ***
#   rasValue_bio5    1    6.801 24.1858  0.001 ***
#   rasValue_bio4    1    6.259 22.2583  0.001 ***
#   rasValue_bio14   1    4.751 16.8949  0.001 ***
#   rasValue_bio18   1    5.804 20.6423  0.001 ***
#   rasValue_bio1    1    4.566 16.2385  0.001 ***
#   rasValue_bio2    1    4.218 15.0016  0.001 ***
#   rasValue_bio17   1    3.906 13.8916  0.001 ***
#   dist_ocean       1    3.582 12.7370  0.001 ***
#   rasValue_swe     1    1.937  6.8889  0.001 ***
#   rasValue_bio15   1    1.642  5.8407  0.001 ***
#   dist_coast       1    1.424  5.0641  0.001 ***
#   rasValue_bio12   1    1.216  4.3251  0.001 ***
#   Residual       832  233.950
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#check for collinearity; threshold should be 20 (above 20 removed); but above 10 should also be inspected
#Removing one by one (or a few at the times if sure that they would result collinear anyway)
vif.cca(fwd.sel)
# rasValue_bio10  rasValue_bio5  rasValue_bio4 rasValue_bio14 rasValue_bio18
# 3900.392957     391.442647     342.728792     168.541877      25.979318
# rasValue_bio1  rasValue_bio2 rasValue_bio17     dist_ocean   rasValue_swe
# 4394.573196       7.882142     122.771696       3.285020      34.218882
# rasValue_bio15     dist_coast rasValue_bio12
# 4.923530       6.756692     372.765064

#remove first variable
#a. removing bio1
dbrda <- capscale(formula = community_ml ~ rasValue_bio10 + rasValue_bio5 +
                    rasValue_bio4 + rasValue_bio14 + rasValue_bio18 +
                    rasValue_bio2 + rasValue_bio17 + dist_ocean + rasValue_swe +
                    rasValue_bio15 + dist_coast + rasValue_bio12, data = data.frame(env.z),
                  distance = "bray")
vif.cca(dbrda)
# rasValue_bio10  rasValue_bio5  rasValue_bio4 rasValue_bio14 rasValue_bio18
# 444.474705     389.357435       5.861416     155.830597      24.886881
# rasValue_bio2 rasValue_bio17     dist_ocean   rasValue_swe rasValue_bio15
# 7.364611     122.591650       3.218868      34.199743       4.560657
# dist_coast rasValue_bio12
# 6.532985     358.532817

#b. removing bio5
dbrda <- capscale(formula = community_ml ~ rasValue_bio10 + 
                    rasValue_bio4 + rasValue_bio14 + rasValue_bio18 +
                    rasValue_bio2 + rasValue_bio17 + dist_ocean + rasValue_swe +
                    rasValue_bio15 + dist_coast + rasValue_bio12, data = data.frame(env.z),
                  distance = "bray")
vif.cca(dbrda)
# rasValue_bio10  rasValue_bio4 rasValue_bio14 rasValue_bio18  rasValue_bio2
# 8.306777       5.757955     155.807564      23.701924       2.313421
# rasValue_bio17     dist_ocean   rasValue_swe rasValue_bio15     dist_coast
# 122.561797       2.496510      33.759705       4.507966       6.499919
# rasValue_bio12
# 358.260144

#c. removing bio14
dbrda <- capscale(formula = community_ml ~ rasValue_bio10 + 
                    rasValue_bio4 + rasValue_bio18 +
                    rasValue_bio2 + rasValue_bio17 + dist_ocean + rasValue_swe +
                    rasValue_bio15 + dist_coast + rasValue_bio12, data = data.frame(env.z),
                  distance = "bray")
vif.cca(dbrda)
# rasValue_bio10  rasValue_bio4 rasValue_bio18  rasValue_bio2 rasValue_bio17
# 7.976151       5.533822      23.322590       2.310000     121.007997
# dist_ocean   rasValue_swe rasValue_bio15     dist_coast rasValue_bio12
# 2.341010      33.056026       4.084768       6.480150     192.605117

#d. removing bio17
dbrda <- capscale(formula = community_ml ~ rasValue_bio10 + 
                    rasValue_bio4 + rasValue_bio18 +
                    rasValue_bio2 + dist_ocean + rasValue_swe +
                    rasValue_bio15 + dist_coast + rasValue_bio12, data = data.frame(env.z),
                  distance = "bray")
vif.cca(dbrda)
# rasValue_bio10  rasValue_bio4 rasValue_bio18  rasValue_bio2     dist_ocean
# 7.976139       5.393123      21.750716       2.219588       2.282378
# rasValue_swe rasValue_bio15     dist_coast rasValue_bio12
# 31.879778       2.015071       6.468961      75.230168

#e. removing bio12
dbrda <- capscale(formula = community_ml ~ rasValue_bio10 + 
                    rasValue_bio4 + rasValue_bio18 +
                    rasValue_bio2 + dist_ocean + rasValue_swe +
                    rasValue_bio15 + dist_coast, data = data.frame(env.z),
                  distance = "bray")
vif.cca(dbrda)
# rasValue_bio10  rasValue_bio4 rasValue_bio18  rasValue_bio2     dist_ocean
# 7.476183       3.825525       8.926693       1.765363       2.247160
# rasValue_swe rasValue_bio15     dist_coast
# 10.968904       1.658991       5.900301

#f. re-add swe
dbrda <- capscale(formula = community_ml ~ rasValue_bio10 + 
                    rasValue_bio4 + rasValue_bio18 + elevation +
                    rasValue_bio2 + dist_ocean + rasValue_swe +
                    rasValue_bio15 + dist_coast, data = data.frame(env.z),
                  distance = "bray")
vif.cca(dbrda)
# rasValue_bio10  rasValue_bio4 rasValue_bio18      elevation  rasValue_bio2
# 13.301574       4.045374       9.220438      11.913786       1.790155
# dist_ocean   rasValue_swe rasValue_bio15     dist_coast
# 2.247573      11.095954       1.678698       7.211495

#now I have my final variables
#retest
fwd.sel <- ordiR2step(dbrda.lwr, #lower model limit
                      scope = formula(dbrda), #upper model limit
                      direction = "forward",
                      R2scope = TRUE, #can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) #change to TRUE to see the selection process!

#check the new model with forward-selected variables
fwd.sel$call
# capscale(formula = community_ml ~ rasValue_bio10 + rasValue_bio4 +
#            rasValue_bio2 + rasValue_bio18 + rasValue_bio15 + dist_ocean +
#            elevation + dist_coast, data = data.frame(env.z), distance = "bray")

#I am not including snow variables then as no relevant! I tried to plot dbrda with also swe but does not look good anyway

#run final dbrda
dbrda_final <- capscale(formula = community_ml ~ rasValue_bio10 + rasValue_bio4 +
                          rasValue_bio2 + rasValue_bio18 + rasValue_bio15 + dist_ocean +
                          elevation + dist_coast, data = data.frame(env.z), distance = "bray")

summary(dbrda_final)
# Partitioning of squared Bray distance:
#   Inertia Proportion
# Total          289.83     1.0000
# Constrained     37.05     0.1278
# Unconstrained  252.78     0.8722

#check R2
RsquareAdj(dbrda_final)
# $r.squared
# [1] 0.1482267
# 
# $adj.r.squared
# [1] 0.1400855

#get statistical significance of the model
anova.cca(fwd.sel, permutations = 1000) #same results as anova()
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 1000
# 
# Model: capscale(formula = community_ml ~ rasValue_bio10 + rasValue_bio4 + rasValue_bio2 + rasValue_bio18 + rasValue_bio15 + dist_ocean + elevation + dist_coast, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs      F   Pr(>F)
# Model      8   37.046 15.333 0.000999 ***
#   Residual 837  252.781
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#get statistical significance of the model for each variable
anova.cca(fwd.sel, step = 1000, by = "term")
# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = community_ml ~ rasValue_bio10 + rasValue_bio4 + rasValue_bio2 + rasValue_bio18 + rasValue_bio15 + dist_ocean + elevation + dist_coast, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs       F Pr(>F)
# rasValue_bio10   1    9.771 32.3522  0.001 ***
#   rasValue_bio4    1    6.404 21.2048  0.001 ***
#   rasValue_bio2    1    5.433 17.9897  0.001 ***
#   rasValue_bio18   1    5.179 17.1500  0.001 ***
#   rasValue_bio15   1    4.335 14.3547  0.001 ***
#   dist_ocean       1    2.887  9.5592  0.001 ***
#   elevation        1    1.794  5.9413  0.001 ***
#   dist_coast       1    1.242  4.1119  0.001 ***
#   Residual       837  252.781
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#plot rda
perc <- round(100*(summary(dbrda_final)$cont$importance[2, 1:2]), 2)

## extract scores - these are coordinates in the RDA space
sc_si <- scores(dbrda_final, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(dbrda_final, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(dbrda_final, display="bp", choices=c(1, 2), scaling=1)

#--> https://www.researchgate.net/post/Variation-partitioning-analysis-in-R
#--> https://rdrr.io/rforge/vegan/man/varpart.html
#--> very cool: https://r.qcbs.ca/workshop10/book-en/variation-partitioning.html

#create color vector
bioregions <- metadata_ml$bioregion
bioregions <- replace(bioregions, bioregions==0, "black")
bioregions <- replace(bioregions, bioregions==1, "#8b0001ff")
bioregions <- replace(bioregions, bioregions==3, "#fe0000ff")
bioregions <- replace(bioregions, bioregions==4, "#ffff00ff")
bioregions <- replace(bioregions, bioregions==6, "#8ebc8bff")
bioregions <- replace(bioregions, bioregions==7, "#32cd33ff")
bioregions <- replace(bioregions, bioregions==8, "#006400ff")
bioregions <- replace(bioregions, bioregions==9, "#28a8f7ff")
bioregions <- replace(bioregions, bioregions==10, "#4682b4ff")
bioregions <- replace(bioregions, bioregions==12, "#808080ff")
bioregions <- replace(bioregions, bioregions==16, "#0000ccff")

#plot rda
library("Cairo")
CairoSVG(file="dbrda_bioclim_ml.svg",
         width=15,
         height=15,
         pointsize=12)

length(rownames(sc_si))
#[1] 846

plot(dbrda_final,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-0.5,0.5), 
     ylim = c(-0.6,0.5),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add sample clusters
pl <- ordihull(dbrda_final, bioregions, draw = c("polygon"), scaling = 1, col = c("black", "darkgreen", "darkolivegreen4", "darkred", "deepskyblue", "deepskyblue4", "dimgray", "green2", "midnightblue", "red", "yellow"), alpha = 40, lty = 0)
# centres and areas of the hulls
summary(pl)
# add points for samples
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = bioregions, # fill colour
       cex = 1.2) # size
# add points for species scores
# add text labels for species abbreviations
#text(sc_si + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
#     labels = rownames(sc_si), 
#     col = "grey40", 
#     font = 2, # bold
#     cex = 0.6)
## add arrows for effects of the explanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2], # end them at the score value
       col = "red",
       lwd = 3)
# add text labels for arrows
text(x = sc_bp[,1] -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)
dev.off()


##########
###pcoa###
##########

library(ape)

#pcoa
varespec.bray <- vegdist(community_ml, method = "bray")
pcoaVS <- pcoa(varespec.bray)

#envfit
env_envfit <- cbind(metadata_ml[,c("rasValue_bio4", "rasValue_bio15", "rasValue_bio10", "rasValue_bio18", "rasValue_bio2", "bioregion")], metadata_el_dist_ml[,c(3:5)])
env_envfit$bioregion <- factor(env_envfit$bioregion, levels = as.character(unique(env_envfit$bioregion)))
en <- envfit(pcoaVS$vectors, env_envfit, permutations = 1000, strata = NULL, choices=c(1,2),  display = "species", na.rm = FALSE)
en
# ***VECTORS
# 
# Axis.1   Axis.2     r2   Pr(>r)
# rasValue_bio4  -0.99043 -0.13805 0.0055 0.077922 .
# rasValue_bio15  0.52123  0.85341 0.0619 0.000999 ***
#   rasValue_bio10  0.69210  0.72180 0.1376 0.000999 ***
#   rasValue_bio18 -0.32598  0.94538 0.0347 0.000999 ***
#   rasValue_bio2  -0.98785  0.15538 0.1125 0.000999 ***
#   dist_coast     -0.71419 -0.69995 0.0697 0.000999 ***
#   dist_ocean     -0.28141 -0.95959 0.1192 0.000999 ***
#   elevation      -0.70428 -0.70992 0.1125 0.000999 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Permutation: free
# Number of permutations: 1000
# 
# ***FACTORS:
#   
#   Centroids:
#   Axis.1  Axis.2
# bioregion6  -0.1127 -0.0326
# bioregion10  0.0221 -0.1514
# bioregion9  -0.0350 -0.0103
# bioregion8   0.2513  0.2525
# bioregion3   0.0037  0.3075
# bioregion16 -0.1071  0.1156
# bioregion7   0.1173  0.0304
# 
# Goodness of fit:
#   r2   Pr(>r)
# bioregion 0.2797 0.000999 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Permutation: free
# Number of permutations: 1000

#colors
col <- env_envfit
col[,1] <- as.character(col[,1])
col[,1][col[,"bioregion"] == "0"] <- "black"
col[,1][col[,"bioregion"] == "1"] <- "#8b0001ff"
col[,1][col[,"bioregion"] == "3"] <- "#fe0000ff"
col[,1][col[,"bioregion"] == "4"] <- "#ffff00ff"
col[,1][col[,"bioregion"] == "6"] <- "#8fbc8bff"
col[,1][col[,"bioregion"] == "7"] <- "#32cd33ff"
col[,1][col[,"bioregion"] == "8"] <- "#006400ff"
col[,1][col[,"bioregion"] == "9"] <- "#26a8f6ff"
col[,1][col[,"bioregion"] == "10"] <- "#4682b3ff"
col[,1][col[,"bioregion"] == "12"] <- "#808080ff"
col[,1][col[,"bioregion"] == "16"] <- "#0000ccff"

#explained variance
pc1 <- pcoaVS$values[1,3]*100
pc1
#[1] 4.272615
pc2 <- pcoaVS$values[2,3]*100
pc2
#[1] 4.225458

#plot - only bioregions 
library("Cairo")
CairoSVG(file="pcoa_bioregion_ml.svg",
    width=10,
    height=9,
    pointsize=12)
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], pch = 21, bg = col[,1], col = "black", cex = 1.5,
     xlab = paste0("PCoA1 (", pc1, "%)"), 
     ylab = paste0("PCoA2 (", pc2, "%)"))
df <- en$factors[1]$centroids
text(df[,1], df[,2], labels = rownames(df), col = "black", cex = 2)
dev.off()

#plot - bioclimatic variables
library("Cairo")
CairoSVG(file="pcoa_bioclimatic_ml.svg",
    width=10,
    height=9,
    pointsize=12)
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], pch = 21, bg = col[,1], col = "black", cex = 1.5,
     xlim = c(-1.2,1), 
     ylim = c(-1.2,1),
     xlab = paste0("PCoA1 (", pc1, "%)"), 
     ylab = paste0("PCoA2 (", pc2, "%)"))
df <- en$factors[1]$centroids
#text(df[,1], df[,2], labels = rownames(df), col = "red")
df <- en$vectors[1]$arrows
arrows(0,0, # start them from (0,0)
       df[,1], df[,2], # end them at the score value
       col = "red",
       lwd = 2)
text(df[,1], df[,2], labels = rownames(df), col = "red")
dev.off()

#############################
###permanova only mainland###
#############################

adonis2(community_ml ~ bioregion, data=metadata_ml, distance="bray", permutations=1000)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = community_ml ~ bioregion, data = metadata_ml, permutations = 1000, distance = "bray")
# Df SumOfSqs      R2      F   Pr(>F)    
# bioregion   1    6.319 0.02528 21.894 0.000999 ***
#   Residual  844  243.606 0.97472                    
# Total     845  249.925 1.00000                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#########################################
###variation partitioning for mainland###
#########################################

#geographical data
geography <- data.frame(metadata_ml[row.names(community_ml),2:1])
geography_pcnm <- as.data.frame(scores(pcnm(dist(geography))))
#environmental data
env.z <- decostand(metadata_ml[,c(3:12,14)], method = "standardize")
round(apply(env.z, 2, mean), 1)
apply(env.z, 2, sd)
env.z <- data.frame(env.z)
env <- env.z[row.names(community_ml),c("rasValue_bio10", "rasValue_bio18",
                                       "rasValue_bio4", "rasValue_bio15", "rasValue_bio2")]

par(mfrow=c(2,1))
#environmental VS geography
var.part <- varpart(community_ml, env, geography_pcnm)
#print(var.part$part)
#plot venn diagram
plot(var.part, 
     Xnames = c("Bioclimatic dataset", "Geography"), # name the partitions
     bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
#test for significance
anova <- anova.cca(rda(community_ml, geography_pcnm))
print(paste("islands geo no ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
#[1] "islands geo no  5.71053789474968   0.001"
anova <- anova.cca(rda(community_ml, env))
print(paste("islands env no ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
#[1] "islands env no  3.90164577515153   0.001"
anova <- anova.cca(rda(community_ml, geography_pcnm, env))
print(paste("islands geo env ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
#[1] "islands geo env  4.47945440543137   0.001"
anova <- anova.cca(rda(community_ml, env, geography_pcnm))
print(paste("islands env geo ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
#[1] "islands env geo  2.67056228583321   0.001"

#distance and elevation
env.z <- decostand(metadata_el_dist_ml[,c(3:5)], method = "standardize")
round(apply(env.z, 2, mean), 1)
apply(env.z, 2, sd)
env.z <- data.frame(env.z)
env_ed <- env.z[row.names(community_ml),]
#environmental VS geography vs el/dist
var.part <- varpart(community_ml, env, geography_pcnm, env_ed)
plot(var.part, 
     Xnames = c("Bioclimatic dataset", "Geography", "Elevation/distance"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "red"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

#saved: variation_partitioning_ml.svg
#1200x1200

##################################
###distance-decay only mainland###
##################################

library(geosphere)

#geographical distance
dist_taxa <- vegdist(community_ml, method="bray") #remove unclassified
d_geo <- distm(metadata_ml[,c(2:1)], fun = distHaversine) #this function works with decimal degrees
dist_geo = as.dist(d_geo)
dc_geo_genera <- mantel(dist_taxa, dist_geo, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_geo_genera
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = dist_taxa, ydis = dist_geo, method = "spearman",      permutations = 1000, na.rm = TRUE) 
# 
# Mantel statistic r: 0.1611 
#       Significance: 0.000999 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0142 0.0179 0.0216 0.0255 
# Permutation: free
# Number of permutations: 1000

#bioclimatic variables
dist_taxa <- vegdist(community_ml, method="bray")
dist_biocl <- vegdist(metadata_ml[,c(3:12,14)], method="eucl")
dc_biocl_genera <- mantel(dist_taxa, dist_biocl, method = "spearman", permutations = 1000, na.rm = TRUE)
dc_biocl_genera
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = dist_taxa, ydis = dist_biocl, method = "spearman",      permutations = 1000, na.rm = TRUE) 
# 
# Mantel statistic r: 0.1102 
#       Significance: 0.000999 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0244 0.0314 0.0388 0.0466 
# Permutation: free
# Number of permutations: 1000