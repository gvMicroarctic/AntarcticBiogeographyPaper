##########
###PCOA###
##########

library(vegan)
library(ape)

setwd("/home/gilda/all")

load(file = "all_dataset_ps_srs_bacteria.RData")

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

#pcoa
varespec.bray <- vegdist(genera_counts_tab_trim, method = "bray")
pcoaVS <- pcoa(varespec.bray)

#envfit
env_envfit <- metadata[row.names(community_tran),c(3:12,14:15)]
env_envfit$bioregion <- factor(env_envfit$bioregion, levels = as.character(unique(env_envfit$bioregion)))
en <- envfit(pcoaVS$vectors, env_envfit, permutations = 1000, strata = NULL, choices=c(1,2),  display = "species", na.rm = FALSE)
en
# ***VECTORS
# 
# Axis.1   Axis.2     r2   Pr(>r)    
# rasValue_bio1   0.95466 -0.29769 0.2125 0.000999 ***
#   rasValue_bio2  -0.74799  0.66372 0.0980 0.000999 ***
#   rasValue_bio4  -0.99943 -0.03368 0.1034 0.000999 ***
#   rasValue_bio5   0.92686 -0.37540 0.2207 0.000999 ***
#   rasValue_bio10  0.91781 -0.39702 0.2265 0.000999 ***
#   rasValue_bio12  0.95957  0.28146 0.1516 0.000999 ***
#   rasValue_bio14  0.94370  0.33081 0.1475 0.000999 ***
#   rasValue_bio15 -0.31627 -0.94867 0.0289 0.000999 ***
#   rasValue_bio17  0.94409  0.32970 0.1409 0.000999 ***
#   rasValue_bio18  0.97659  0.21509 0.1648 0.000999 ***
#   rasValue_swe    0.96326  0.26858 0.1323 0.000999 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Permutation: free
# Number of permutations: 1000
# 
# ***FACTORS:
#   
#   Centroids:
#   Axis.1  Axis.2
# bioregion6  -0.0761  0.1083
# bioregion10 -0.1760 -0.0223
# bioregion9  -0.0339  0.0277
# bioregion1   0.0054  0.0483
# bioregion0   0.2219 -0.0242
# bioregion12  0.2739 -0.0554
# bioregion3   0.2271  0.0537
# bioregion8   0.2006 -0.2618
# bioregion16  0.0773  0.0990
# bioregion7   0.0065 -0.1124
# bioregion4   0.0583  0.0272
# 
# Goodness of fit:
#   r2   Pr(>r)    
# bioregion 0.333 0.000999 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Permutation: free
# Number of permutations: 1000

#colors
col <- env_envfit
col[,11] <- as.character(col[,11])
col[,11][col[,"bioregion"] == "0"] <- "black"
col[,11][col[,"bioregion"] == "1"] <- "#8b0001ff"
col[,11][col[,"bioregion"] == "3"] <- "#fe0000ff"
col[,11][col[,"bioregion"] == "4"] <- "#ffff00ff"
col[,11][col[,"bioregion"] == "6"] <- "#8fbc8bff"
col[,11][col[,"bioregion"] == "7"] <- "#32cd33ff"
col[,11][col[,"bioregion"] == "8"] <- "#006400ff"
col[,11][col[,"bioregion"] == "9"] <- "#26a8f6ff"
col[,11][col[,"bioregion"] == "10"] <- "#4682b3ff"
col[,11][col[,"bioregion"] == "12"] <- "#808080ff"
col[,11][col[,"bioregion"] == "16"] <- "#0000ccff"

#explained variance
pc1 <- pcoaVS$values[1,3]*100
pc1
#[1] 4.22761
pc2 <- pcoaVS$values[2,3]*100
pc2
#[1] 3.630098

#plot - only bioregions 
svg(filename="pcoa_bioregion.svg",
    width=10,
    height=9,
    pointsize=12)
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], pch = 21, bg = col[,11], col = "black", cex = 1.5,
     xlab = paste0("PCoA1 (", pc1, "%)"), 
     ylab = paste0("PCoA2 (", pc2, "%)"))
df <- en$factors[1]$centroids
text(df[,1], df[,2], labels = rownames(df), col = "black", cex = 2)
dev.off()

#plot - bioclimatic variables
svg(filename="pcoa_bioclimatic.svg",
    width=10,
    height=9,
    pointsize=12)
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], pch = 21, bg = col[,11], col = "black", cex = 1.5,
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

#color by island
col <- metadata[row.names(community_tran),]
col[,11] <- as.character(col[,11])
col[,11][col[,"island"] == "n"] <- "white"
col[,11][col[,"bioregion"] == "0" & col[,"island"] == "y"] <- "black"
col[,11][col[,"bioregion"] == "1" & col[,"island"] == "y"] <- "#8b0001ff"
col[,11][col[,"bioregion"] == "3" & col[,"island"] == "y"] <- "#fe0000ff"
col[,11][col[,"bioregion"] == "4" & col[,"island"] == "y"] <- "#ffff00ff"
col[,11][col[,"bioregion"] == "6" & col[,"island"] == "y"] <- "#8fbc8bff"
col[,11][col[,"bioregion"] == "7" & col[,"island"] == "y"] <- "#32cd33ff"
col[,11][col[,"bioregion"] == "8" & col[,"island"] == "y"] <- "#006400ff"
col[,11][col[,"bioregion"] == "9" & col[,"island"] == "y"] <- "#26a8f6ff"
col[,11][col[,"bioregion"] == "10" & col[,"island"] == "y"] <- "#4682b3ff"
col[,11][col[,"bioregion"] == "12" & col[,"island"] == "y"] <- "#808080ff"
col[,11][col[,"bioregion"] == "16" & col[,"island"] == "y"] <- "#0000ccff"

#explained variance
pc1 <- pcoaVS$values[1,3]*100
pc1
#[1] 4.22761
pc2 <- pcoaVS$values[2,3]*100
pc2
#[1] 3.630098

#plot - only bioregions 
svg(filename="pcoa_bioregion_col_islands.svg",
    width=10,
    height=9,
    pointsize=12)
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], pch = 21, bg = col[,11], col = "black", cex = 1.5,
     xlab = paste0("PCoA1 (", pc1, "%)"), 
     ylab = paste0("PCoA2 (", pc2, "%)"))
df <- en$factors[1]$centroids
text(df[,1], df[,2], labels = rownames(df), col = "black", cex = 2)
dev.off()

#save.image(file = "all_dataset_ps_srs_pcoa_bacteria.RData")

################################################
###PCOA highlighting samples from bioregion 7###
################################################

library(vegan)
library(ape)

setwd("/home/gilda/all")

load(file = "all_dataset_ps_srs_bacteria.RData")

metadata <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)[samples,]
metadata$rasValue_swe[is.na(metadata$rasValue_swe)] <- 0

bio7 <- read.csv("bio7.csv", row.names = 1)[samples,]
metadata$dataset <- bio7$dataset

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

#pcoa
varespec.bray <- vegdist(genera_counts_tab_trim, method = "bray")
pcoaVS <- pcoa(varespec.bray)

#colors
col <- metadata
col[,11] <- as.character(col[,11])
col[,11][col[,"bioregion"] == "0"] <- "white"
col[,11][col[,"bioregion"] == "1"] <- "white"
col[,11][col[,"bioregion"] == "3"] <- "white"
col[,11][col[,"bioregion"] == "4"] <- "white"
col[,11][col[,"bioregion"] == "6"] <- "white"
col[,11][col[,"bioregion"] == "7"] <- "#32cd33ff"
col[,11][col[,"dataset"] == "dt16_windmill"] <- "red"
col[,11][col[,"dataset"] == "dt16_vestfold"] <- "blue"
col[,11][col[,"bioregion"] == "8"] <- "white"
col[,11][col[,"bioregion"] == "9"] <- "white"
col[,11][col[,"bioregion"] == "10"] <- "white"
col[,11][col[,"bioregion"] == "12"] <- "white"
col[,11][col[,"bioregion"] == "16"] <- "white"

#explained variance
pc1 <- pcoaVS$values[1,3]*100
pc1
#[1] 4.22761
pc2 <- pcoaVS$values[2,3]*100
pc2
#[1] 3.630098

#plot - bio7 highlighted
svglite("pcoa_bio7.svg", width=10, height=9)
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], pch = 21, bg = col[,11], col = "black", cex = 1.5,
     xlab = paste0("PCoA1 (", pc1, "%)"), 
     ylab = paste0("PCoA2 (", pc2, "%)"))
dev.off()

#PERMANOVA
library("vegan")

metadata_bio7 <- metadata[metadata[,"dataset"] == "vestfold" | metadata[,"dataset"] == "windmill",] 
dim(metadata_bio7)
#[1] 191  20

community_bio7 <- community_tran[row.names(metadata_bio7),]
dim(community_bio7)

#run permanova
adonis2(community_bio7 ~ dataset, data=metadata_bio7, distance="bray", permutations=1000)
#adonis2(formula = community_bio7 ~ dataset, data = metadata_bio7, permutations = 1000, distance = "bray")
#Df SumOfSqs      R2      F   Pr(>F)    
#dataset    1    7.694 0.14418 31.841 0.000999 ***
#  Residual 189   45.671 0.85582                    
#Total    190   53.365 1.00000                         