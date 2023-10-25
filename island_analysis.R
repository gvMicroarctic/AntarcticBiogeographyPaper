##########################
###Analyses for islands###
##########################

setwd("/home/gilda/all")

#import community information
load(file = "all_dataset_ps_srs_bacteria_hell.RData")
community <- t(genera_proportions_tab_hellinger[-1,samples])
dim(community)
#[1]  988 1445

#import metadata
metadata <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)[samples,]
dim(metadata)
#[1] 988  20

#get only islands
metadata_isl <- metadata[metadata[,"island"] == "y",]
metadata_isl$rasValue_swe[is.na(metadata_isl$rasValue_swe)] <- 0
dim(metadata_isl)
#[1] 142  20

#get community information - only islands
community_isl0 <- community[row.names(metadata_isl),]
dim(community_isl0)
#[1]  142 1445
community_isl <- community_isl0[rowSums(community_isl0)>0,]
dim(community_isl)
#[1]  142 1445

###########################################
###how many island samples per bioregion###
###########################################

unique(metadata_isl$bioregion)
#[1]  6  1  0 12  3  9  7 16  4

nrow(metadata_isl[metadata_isl$bioregion == 0,])
#[1] 24
nrow(metadata_isl[metadata_isl$bioregion == 1,])
#[1] 6
nrow(metadata_isl[metadata_isl$bioregion == 3,])
#[1] 78
nrow(metadata_isl[metadata_isl$bioregion == 4,])
#[1] 8
nrow(metadata_isl[metadata_isl$bioregion == 6,])
#[1] 1
nrow(metadata_isl[metadata_isl$bioregion == 7,])
#[1] 11
nrow(metadata_isl[metadata_isl$bioregion == 9,])
#[1] 2
nrow(metadata_isl[metadata_isl$bioregion == 12,])
#[1] 11
nrow(metadata_isl[metadata_isl$bioregion == 16,])
#[1] 1

########################################
###permanova for islands vs continent###
########################################

library("vegan")

#run permanova
adonis2(community ~ island, data=metadata, distance="bray", permutations=1000)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = community ~ island, data = metadata, permutations = 1000, distance = "bray")
# Df SumOfSqs      R2      F   Pr(>F)    
# island     1   10.645 0.03615 36.985 0.000999 ***
#   Residual 986  283.776 0.96385                    
# Total    987  294.421 1.00000                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

################################
###permanova only for islands###
################################

library("vegan")

#run permanova
adonis2(community_isl ~ island_sp, data=metadata_isl, distance="bray", permutations=1000)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = community_isl ~ island_sp, data = metadata_isl, permutations = 1000, distance = "bray")
# Df SumOfSqs      R2      F   Pr(>F)    
# island_sp  29   20.141 0.59499 5.6737 0.000999 ***
#   Residual  112   13.710 0.40501                    
# Total     141   33.851 1.00000                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


for (n in 1:length(unique(metadata_isl$island_sp))) {
  num <- nrow(metadata_isl[metadata_isl$island_sp == unique(metadata_isl$island_sp)[n],])
  print(paste0(n, " ", unique(metadata_isl$island_sp)[n], " ", num))
}
# [1] "1 no_name 1"
# [1] "2 james_ross 6"
# [1] "3 bouvetoya 4"
# [1] "4 possession 3"
# [1] "5 bartolome 2"
# [1] "6 kerguelen 2"
# [1] "7 maher 2"
# [1] "8 lauff 2"
# [1] "9 siple 7"
# [1] "10 marion 1"
# [1] "11 peter 1"
# [1] "12 scott 6"
# [1] "13 south_georgia 5"
# [1] "14 king_george 26"
# [1] "15 ross 2"
# [1] "16 useful 1"
# [1] "17 wiencke 4"
# [1] "18 hop 5"
# [1] "19 samson 2"
# [1] "20 alexander 8"
# [1] "21 spert 2"
# [1] "22 alectoria 2"
# [1] "23 deception 2"
# [1] "24 blaiklock 1"
# [1] "25 adelaide 1"
# [1] "26 robert 1"
# [1] "27 livingston 27"
# [1] "28 west_lagoon 1"
# [1] "29 anvers 10"
# [1] "30 herring 5"

for (n in 1:length(unique(metadata_isl$island_sp))) {
  num <- nrow(metadata_isl[metadata_isl$island_sp == unique(metadata_isl$island_sp)[n],])
  if (num >= 4) {
    print(paste0(n, " ", unique(metadata_isl$island_sp)[n], " ", num)) 
  }
}
# [1] "2 james_ross 6"
# [1] "3 bouvetoya 4"
# [1] "9 siple 7"
# [1] "12 scott 6"
# [1] "13 south_georgia 5"
# [1] "14 king_george 26"
# [1] "17 wiencke 4"
# [1] "18 hop 5"
# [1] "20 alexander 8"
# [1] "27 livingston 27"
# [1] "29 anvers 10"
# [1] "30 herring 5"

metadata_isl_mj <- metadata_isl[metadata_isl$island_sp == "james_ross" | metadata_isl$island_sp == "bouvetoya" |
                                  metadata_isl$island_sp == "siple" | metadata_isl$island_sp == "scott" |
                                  metadata_isl$island_sp == "south_georgia" | metadata_isl$island_sp == "king_george" |
                                  metadata_isl$island_sp == "wiencke" | metadata_isl$island_sp == "hop" |
                                  metadata_isl$island_sp == "alexander" | metadata_isl$island_sp == "livingston" |
                                  metadata_isl$island_sp == "anvers" | metadata_isl$island_sp == "herring",]
dim(metadata_isl_mj)
#[1] 113  20

community_isl_mj <- community_isl[row.names(metadata_isl_mj),]
dim(community_isl_mj)
#[1]  113 1445

#run permanova
adonis2(community_isl_mj ~ island_sp, data=metadata_isl_mj, distance="bray", permutations=1000)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = community_isl_mj ~ island_sp, data = metadata_isl_mj, permutations = 1000, distance = "bray")
# Df SumOfSqs      R2      F   Pr(>F)    
# island_sp  11   13.318 0.53031 10.367 0.000999 ***
#   Residual  101   11.796 0.46969                    
# Total     112   25.114 1.00000                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

adonis2(community_isl_mj ~ bioregion, data=metadata_isl_mj, distance="bray", permutations=1000)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = community_isl_mj ~ bioregion, data = metadata_isl_mj, permutations = 1000, distance = "bray")
# Df SumOfSqs      R2      F   Pr(>F)    
# bioregion   1   1.0809 0.04304 4.9924 0.000999 ***
#   Residual  111  24.0328 0.95696                    
# Total     112  25.1138 1.00000                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

###########
###dbRDA###
###########

#hellinger transformed genera dataset
community_isl

#standardize environmental variables
library(vegan)
env.z <- decostand(metadata_isl[,c(3:12,14)], method = "standardize")
round(apply(env.z, 2, mean), 1)
apply(env.z, 2, sd)
env.z <- data.frame(env.z)

#run RDA
dbrda.upr <- capscale(community_isl ~ ., data = data.frame(env.z), dist="bray") # upper model limit (the "full" model; it takes in account all variables)
dbrda.lwr <- capscale(community_isl ~ 1, data = data.frame(env.z), dist="bray") # lower model limit

#run the ordistep model to identify variables with explanatory power
fwd.sel <- ordiR2step(dbrda.lwr, #lower model limit
                      scope = formula(dbrda.upr), #upper model limit
                      direction = "forward",
                      R2scope = TRUE, #can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) #change to TRUE to see the selection process!

#check the new model with forward-selected variables
fwd.sel$call
# capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 + 
#            rasValue_bio10 + rasValue_swe + rasValue_bio14 + rasValue_bio18 + 
#            rasValue_bio5 + rasValue_bio1 + rasValue_bio2 + rasValue_bio17, 
#          data = data.frame(env.z), distance = "bray")

#get statistical significance of the model
anova.cca(fwd.sel, permutations = 1000) #same results as anova()
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 1000
# 
# Model: capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 + rasValue_bio10 + rasValue_swe + rasValue_bio14 + rasValue_bio18 + rasValue_bio5 + rasValue_bio1 + rasValue_bio2 + rasValue_bio17, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs      F   Pr(>F)    
# Model     10   10.443 5.6817 0.000999 ***
#   Residual 131   24.078                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#get statistical significance of the model for each variable
anova.cca(fwd.sel, step = 1000, by = "term")
# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 + rasValue_bio10 + rasValue_swe + rasValue_bio14 + rasValue_bio18 + rasValue_bio5 + rasValue_bio1 + rasValue_bio2 + rasValue_bio17, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs       F Pr(>F)    
# rasValue_bio4    1   2.3947 13.0288  0.001 ***
#   rasValue_bio15   1   1.1905  6.4770  0.001 ***
#   rasValue_bio10   1   0.9920  5.3972  0.001 ***
#   rasValue_swe     1   0.8405  4.5730  0.001 ***
#   rasValue_bio14   1   0.9174  4.9909  0.001 ***
#   rasValue_bio18   1   1.3681  7.4431  0.001 ***
#   rasValue_bio5    1   0.8230  4.4775  0.001 ***
#   rasValue_bio1    1   0.7019  3.8189  0.001 ***
#   rasValue_bio2    1   0.6779  3.6881  0.001 ***
#   rasValue_bio17   1   0.5371  2.9221  0.001 ***
#   Residual       131  24.0783                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#check for collinearity; threshold should be 20 (above 20 removed); but above 10 should also be inspected
#Removing one by one (or a few at the times if sure that they would result collinear anyway)
vif.cca(fwd.sel)
# rasValue_bio4 rasValue_bio15 rasValue_bio10   rasValue_swe rasValue_bio14 rasValue_bio18  rasValue_bio5  rasValue_bio1 
# 469.365696       5.412373     925.650513       4.657005     175.708092      61.392232     140.944143    1491.287438 
# rasValue_bio2 rasValue_bio17 
# 19.904102     188.470348 

#remove first variable
#a. removing bio1
dbrda <- capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 +
                    rasValue_bio10 + rasValue_swe + rasValue_bio14 + rasValue_bio18 +
                    rasValue_bio5 + rasValue_bio2 + rasValue_bio17,
                  data = data.frame(env.z), distance = "bray")
vif.cca(dbrda)
# rasValue_bio4 rasValue_bio15 rasValue_bio10   rasValue_swe rasValue_bio14 rasValue_bio18  rasValue_bio5  rasValue_bio2 
# 13.609772       4.741590      47.686772       4.405763     175.689335      61.266141      42.379249      11.708347 
# rasValue_bio17 
# 188.462813

#b. removing bio17
dbrda <- capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 +
                    rasValue_bio10 + rasValue_swe + rasValue_bio14 + rasValue_bio18 +
                    rasValue_bio5 + rasValue_bio2,
                  data = data.frame(env.z), distance = "bray")
vif.cca(dbrda)
# rasValue_bio4 rasValue_bio15 rasValue_bio10   rasValue_swe rasValue_bio14 rasValue_bio18  rasValue_bio5  rasValue_bio2 
# 13.566443       4.737125      46.954543       3.019945      43.462726      59.740920      39.443158      11.589380

#c. removing bio14
dbrda <- capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 +
                    rasValue_bio10 + rasValue_swe + rasValue_bio18 +
                    rasValue_bio5 + rasValue_bio2,
                  data = data.frame(env.z), distance = "bray")
vif.cca(dbrda)
# rasValue_bio4 rasValue_bio15 rasValue_bio10   rasValue_swe rasValue_bio18  rasValue_bio5  rasValue_bio2 
# 12.910386       4.530233      46.129141       2.998435       5.347852      39.376636      10.777145 

#d. removing bio5
dbrda <- capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 +
                    rasValue_bio10 + rasValue_swe + rasValue_bio18 +
                    rasValue_bio2,
                  data = data.frame(env.z), distance = "bray")
vif.cca(dbrda)
# rasValue_bio4 rasValue_bio15 rasValue_bio10   rasValue_swe rasValue_bio18  rasValue_bio2 
# 12.335388       4.033083       4.331366       2.121170       4.560644       5.184200

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
# capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 + 
#            rasValue_bio10 + rasValue_swe + rasValue_bio18 + rasValue_bio2, 
#          data = data.frame(env.z), distance = "bray")

#I am not including snow variables then as no relevant! I tried to plot dbrda with also swe but does not look good anyway

#run final dbrda
dbrda_final <- capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 + 
                          rasValue_bio10 + rasValue_swe + rasValue_bio18 + rasValue_bio2, 
                        data = data.frame(env.z), distance = "bray")

summary(dbrda_final)
# Partitioning of squared Bray distance:
#   Inertia Proportion
# Total          34.521     1.0000
# Constrained     6.884     0.1994
# Unconstrained  27.638     0.8006

#check R2
RsquareAdj(dbrda_final)
# $r.squared
# [1] 0.203353
# 
# $adj.r.squared
# [1] 0.1679465

#get statistical significance of the model
anova.cca(fwd.sel, permutations = 1000) #same results as anova()
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 1000
# 
# Model: capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 + rasValue_bio10 + rasValue_swe + rasValue_bio18 + rasValue_bio2, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs      F   Pr(>F)    
# Model      6   6.8838 5.6041 0.000999 ***
#   Residual 135  27.6377                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#get statistical significance of the model for each variable
anova.cca(fwd.sel, step = 1000, by = "term")
# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = community_isl ~ rasValue_bio4 + rasValue_bio15 + rasValue_bio10 + rasValue_swe + rasValue_bio18 + rasValue_bio2, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs       F Pr(>F)    
# rasValue_bio4    1   2.3947 11.6975  0.001 ***
#   rasValue_bio15   1   1.1905  5.8152  0.001 ***
#   rasValue_bio10   1   0.9920  4.8457  0.001 ***
#   rasValue_swe     1   0.8405  4.1058  0.001 ***
#   rasValue_bio18   1   0.9162  4.4752  0.001 ***
#   rasValue_bio2    1   0.5498  2.6854  0.002 ** 
#   Residual       135  27.6377                   
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
bioregions <- metadata_isl$bioregion
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
CairoSVG(file="dbrda_bioclim_isl.svg",
         width=15,
         height=15,
         pointsize=12)

length(rownames(sc_si))
#[1] 142

plot(dbrda_final,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-0.5,0.5), 
     ylim = c(-0.5,0.5),
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
###PCOA###
##########

#######library(ecodist)
library(vegan)
library(ape)

#pcoa
varespec.bray <- vegdist(community_isl, method = "bray")
pcoaVS <- pcoa(varespec.bray)

#envfit
env_envfit <- metadata_isl[row.names(community_isl),c("rasValue_bio4", "rasValue_bio15", "rasValue_bio10", "rasValue_swe", "rasValue_bio18", "rasValue_bio2", "bioregion")]
env_envfit$bioregion <- factor(env_envfit$bioregion, levels = as.character(unique(env_envfit$bioregion)))
en <- envfit(pcoaVS$vectors, env_envfit, permutations = 1000, strata = NULL, choices=c(1,2),  display = "species", na.rm = FALSE)
en
# ***VECTORS
# 
# Axis.1   Axis.2     r2   Pr(>r)    
# rasValue_bio4   0.35005 -0.93673 0.3853 0.000999 ***
#   rasValue_bio15  0.42129 -0.90693 0.4065 0.000999 ***
#   rasValue_bio10 -0.45571  0.89013 0.0999 0.001998 ** 
#   rasValue_swe    0.17192  0.98511 0.1811 0.000999 ***
#   rasValue_bio18 -0.38782  0.92173 0.2281 0.000999 ***
#   rasValue_bio2   0.29684 -0.95493 0.2382 0.000999 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Permutation: free
# Number of permutations: 1000
# 
# ***FACTORS:
#   
#   Centroids:
#   Axis.1  Axis.2
# bioregion6  -0.2058 -0.0388
# bioregion1   0.0224  0.0281
# bioregion0  -0.0295 -0.1167
# bioregion12  0.1951 -0.1718
# bioregion3  -0.0255  0.0966
# bioregion9  -0.0635 -0.3318
# bioregion7   0.0225 -0.2137
# bioregion16  0.0563 -0.2282
# bioregion4   0.0560  0.0332
# 
# Goodness of fit:
#   r2   Pr(>r)    
# bioregion 0.3327 0.000999 ***
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
col[,1][col[,"bioregion"] == "9"] <- "#26a8f6ff"
col[,1][col[,"bioregion"] == "12"] <- "#808080ff"
col[,1][col[,"bioregion"] == "16"] <- "#0000ccff"

#explained variance
pc1 <- pcoaVS$values[1,3]*100
pc1
#[1] 11.82147
pc2 <- pcoaVS$values[2,3]*100
pc2
#[1] 8.257236

#plot - only bioregions 
svg(filename="pcoa_isl_bioregion.svg",
    width=10,
    height=9,
    pointsize=12)
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], pch = 16, col = col[,1], ps = 1.5,
     xlab = paste0("PCoA1 (", pc1, "%)"), 
     ylab = paste0("PCoA2 (", pc2, "%)"), cex =2)
df <- en$factors[1]$centroids
text(df[,1], df[,2], labels = rownames(df), col = "black", cex = 1)
dev.off()

svg(filename="pcoa_isl_samples.svg",
    width=10,
    height=9,
    pointsize=12)
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], pch = 16, col = col[,1], ps = 1.5,
     xlab = paste0("PCoA1 (", pc1, "%)"), 
     ylab = paste0("PCoA2 (", pc2, "%)"), cex =2)
df <- en$factors[1]$centroids
text(pcoaVS$vectors[,1], pcoaVS$vectors[,2], labels = rownames(community_isl), pos = 2, col = "black", cex = 1)
dev.off()

#plot - bioclimatic variables
svg(filename="pcoa_isl_bioclimatic.svg",
    width=10,
    height=9,
    pointsize=12)
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], pch = 16, col = col[,1], 
     xlim = c(-1.2,1), 
     ylim = c(-1.2,1),
     xlab = paste0("PCoA1 (", pc1, "%)"), 
     ylab = paste0("PCoA2 (", pc2, "%)"))
df <- en$factors[1]$centroids
#text(df[,1], df[,2], labels = rownames(df), col = "red")
df <- en$vectors[1]$arrows
arrows(0,0, # start them from (0,0)
       df[,1], df[,2], # end them at the score value
       col = "blue",
       lwd = 2)
text(df[,1], df[,2], labels = rownames(df), col = "blue")
dev.off()

unique(metadata_isl$island_sp)

#get only islands from bioregion 3
metadata_bio3 <- metadata[metadata[,"bioregion"] == "3" & metadata[,"island"] == "y",]
dim(metadata_bio3)
#[1] 78  20
unique(metadata_bio3$island_sp)

for (n in 1:length(unique(metadata_bio3$island_sp))) {
  num <- nrow(metadata_bio3[metadata_bio3$island_sp == unique(metadata_bio3$island_sp)[n],])
  print(paste0(n, " ", unique(metadata_bio3$island_sp)[n], " ", num))
}
# [1] "1 king_george 26"
# [1] "2 useful 1"
# [1] "3 wiencke 4"
# [1] "4 spert 2"
# [1] "5 alectoria 2"
# [1] "6 deception 2"
# [1] "7 blaiklock 1"
# [1] "8 adelaide 1"
# [1] "9 robert 1"
# [1] "10 livingston 27"
# [1] "11 west_lagoon 1"
# [1] "12 anvers 10"

#changing shape only when more than 2 sammples per island
#1, 3, 10, 12

sh <- metadata_isl
sh[,1] <- as.character(sh[,1])
sh[,1][sh[,"bioregion"] == "0"] <- 16
sh[,1][sh[,"bioregion"] == "1"] <- 16
sh[,1][sh[,"bioregion"] == "3"] <- 16
sh[,1][sh[,"bioregion"] == "4"] <- 16
sh[,1][sh[,"bioregion"] == "6"] <- 16
sh[,1][sh[,"bioregion"] == "7"] <- 16
sh[,1][sh[,"bioregion"] == "9"] <- 16
sh[,1][sh[,"bioregion"] == "12"] <- 16
sh[,1][sh[,"bioregion"] == "16"] <- 16

sh[,1][sh[,"island_sp"] == unique(metadata_bio3$island_sp)[1]] <- 0
sh[,1][sh[,"island_sp"] == unique(metadata_bio3$island_sp)[3]] <- 2
sh[,1][sh[,"island_sp"] == unique(metadata_bio3$island_sp)[10]] <- 9
sh[,1][sh[,"island_sp"] == unique(metadata_bio3$island_sp)[12]] <- 11

sh[,1] <- as.numeric(sh[,1])

#plot - bioclimatic variables
svg(filename="pcoa_isl_bio3.svg",
    width=10,
    height=9,
    pointsize=12)
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], pch = sh[,1], col = col[,1], ps = 1.5,
     xlab = paste0("PCoA1 (", pc1, "%)"), 
     ylab = paste0("PCoA2 (", pc2, "%)"), cex =2)
df <- en$factors[1]$centroids
dev.off()

############################
###Variation Partitioning###
############################

#geographical data
geography <- data.frame(metadata_isl[,2:1])
geography_pcnm <- as.data.frame(scores(pcnm(dist(geography))))

library(vegan)
env.z <- decostand(metadata_isl[,c(3:12,14)], method = "standardize")
round(apply(env.z, 2, mean), 1)
apply(env.z, 2, sd)
env.z <- data.frame(env.z)

#environmental data
env <- env.z[row.names(community_isl),c("rasValue_bio4", "rasValue_bio15", "rasValue_bio10", "rasValue_swe", "rasValue_bio18", "rasValue_bio2")]

#environmental VS geography
var.part <- varpart(community_isl, env, geography_pcnm)
print(var.part$part)
# No. of explanatory tables: 2 
# Total variation (SS): 4669.7 
# Variance: 33.119 
# No. of observations: 142 
# 
# Partition table:
#   Df R.squared Adj.R.squared Testable
# [a+c] = X1            6   0.17497       0.13830     TRUE
# [b+c] = X2           13   0.30764       0.23732     TRUE
# [a+b+c] = X1+X2      19   0.41960       0.32922     TRUE
# Individual fractions                                    
# [a] = X1|X2           6                 0.09190     TRUE
# [b] = X2|X1          13                 0.19091     TRUE
# [c]                   0                 0.04641    FALSE
# [d] = Residuals                         0.67078    FALSE
# ---
#   Use function 'rda' to test significance of fractions of interest

#plot venn diagram
plot(var.part, 
     Xnames = c("Bioclimatic dataset", "Geography"), # name the partitions
     bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)