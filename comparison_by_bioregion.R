#######################################
###Test specifics for each bioregion###
#######################################

setwd("/home/gilda/all")

#import community information
load(file = "all_dataset_ps_srs_analysis_bacteria.RData")
community <- t(genera_proportions_tab_hellinger[-1,])
dim(community)
#[1]  987 1412

#import information about bioregions and bioclimatic variables
bioclim <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)
dim(bioclim)
#[1] 987  20

#distance-decay
library(vegan)
library(geosphere)
library(phyloseq)

for (n in sort(unique(bioclim[,"bioregion"]))) {
  #trim datasets
  bioclim_trim <- bioclim[bioclim[,"bioregion"] == n,]
  community_trim <- community[row.names(bioclim_trim),]
  community_trim <- community_trim[rowSums(community_trim) > 0,]
  
  #geography
  geo_formatted <- data.frame(bioclim_trim[row.names(community_trim),2:1])
  dist_taxa <- vegdist(community_trim, method="bray") #remove unclassified
  d_geo <- distm(geo_formatted, fun = distHaversine) #this function works with decimal degrees
  dist_geo = as.dist(d_geo)
  dc_genera <- mantel(dist_taxa, dist_geo, method = "spearman", permutations = 1000, na.rm = TRUE)
  stat <- dc_genera$statistic
  pvalue <- dc_genera$signif
  print(paste("geo ", n, " ", nrow(geo_formatted), " ", stat, " ", pvalue))
  
  #bioclimatic variables
  geo_formatted <- data.frame(bioclim_trim[row.names(community_trim),3:12])
  dist_taxa <- vegdist(community_trim, method="bray") #remove unclassified
  dist_geo <- vegdist(geo_formatted, method="eucl")
  dc_genera <- mantel(dist_taxa, dist_geo, method = "spearman", permutations = 1000, na.rm = TRUE)
  stat <- dc_genera$statistic
  pvalue <- dc_genera$signif
  print(paste("env ", n, " ", nrow(geo_formatted), " ", stat, " ", pvalue))
}

################################
###Get variation partitioning###
################################

#https://r.qcbs.ca/workshop10/book-en/variation-partitioning.html
#https://rdrr.io/cran/vegan/man/varpart.html

setwd("/home/gilda/all")

library("vegan")

#import community information
load(file = "all_dataset_ps_srs_bacteria.RData")
community <- t(genera_proportions_tab_hellinger[-1,samples])
dim(community)
#[1]  988 1445

#bioclimatic variables
bioclim <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names=1)[samples,]
dim(bioclim)
#[1] 988  20
bioclim_trim <- bioclim[,c(3:12,14)]
bioclim_trim$rasValue_swe[is.na(bioclim_trim$rasValue_swe)] <- 0
dim(bioclim_trim)
#[1] 988  11

#standardize environmental variables
library(vegan)
env.z <- decostand(bioclim_trim, method = "standardize")
round(apply(env.z, 2, mean), 1)
apply(env.z, 2, sd)
env.z <- data.frame(env.z)

#geographical data
geography <- bioclim[,c(2,1)]
geography_pcnm <- as.data.frame(scores(pcnm(dist(geography))))

par(mfrow=c(5,2))

for (n in c(0, 3, 6, 7, 8, 9, 10, 16)) { #not processing 1, 4 and 12
  
  #trim datasets
  bioclim_trim <- bioclim[bioclim[,"bioregion"] == n,]
  community_trim <- community[row.names(bioclim_trim),]
  community_trim <- data.frame(community_trim[rowSums(community_trim) > 0,])
  
  #geographical data
  geography <- data.frame(bioclim_trim[row.names(community_trim),2:1])
  geography_pcnm <- as.data.frame(scores(pcnm(dist(geography))))
  
  #environmental data
  env <- env.z[row.names(community_trim),c("rasValue_bio10", "rasValue_bio18",
                  "rasValue_bio4", "rasValue_bio15", "rasValue_bio2")]
  
  #environmental VS geography
  var.part <- varpart(community_trim, env, geography_pcnm)
  #print(var.part$part)
  
  #plot venn diagram
  plot(var.part, 
       Xnames = c("Bioclimatic dataset", "Geography"), # name the partitions
       bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
       digits = 2, # only show 2 digits
       cex = 0.5)
  title(main = n)
  
  #test for significance
  anova <- anova.cca(rda(community_trim, geography_pcnm))
  print(paste(n, " geo no ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
  
  anova <- anova.cca(rda(community_trim, env))
  print(paste(n, " env no ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
  
  anova <- anova.cca(rda(community_trim, geography_pcnm, env))
  print(paste(n, " geo env ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
  
  anova <- anova.cca(rda(community_trim, env, geography_pcnm))
  print(paste(n, " env geo ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
}
# [1] "0  geo no  15.3994113632886   0.001"
# [1] "0  env no  18.1723632634114   0.001"
# [1] "0  geo env  9.55879427229223   0.002"
# [1] "0  env geo  12.331746172415   0.001"
# [1] "3  geo no  7.86706694013739   0.001"
# [1] "3  env no  5.74223663805503   0.001"
# [1] "3  geo env  6.04686856730887   0.001"
# [1] "3  env geo  3.92203826522651   0.001"
# [1] "6  geo no  5.25655148305462   0.001"
# [1] "6  env no  4.17331282169735   0.001"
# [1] "6  geo env  2.1010382585831   0.001"
# [1] "6  env geo  1.01779959722583   0.099"
# [1] "7  geo no  11.7093165851801   0.001"
# [1] "7  env no  9.3334903884564   0.001"
# [1] "7  geo env  4.38819232123465   0.001"
# [1] "7  env geo  2.01236612451092   0.001"
# [1] "8  geo no  26.9260281222129   0.001"
# [1] "8  env no  29.9629671419374   0.001"
# [1] "8  geo env  9.66153919722618   0.001"
# [1] "8  env geo  12.6984782169507   0.001"
# [1] "9  geo no  5.90767910083057   0.001"
# [1] "9  env no  4.02456583897281   0.001"
# [1] "9  geo env  5.43507809417663   0.001"
# [1] "9  env geo  3.55196483231888   0.001"
# [1] "10  geo no  10.004094361412   0.001"
# [1] "10  env no  6.65705911029103   0.001"
# [1] "10  geo env  5.64795918197857   0.001"
# [1] "10  env geo  2.30092393085763   0.001"
# [1] "16  geo no  5.2573806067532   0.001"
# [1] "16  env no  2.8258694467161   0.001"
# [1] "16  geo env  3.89164946432898   0.041"
# [1] "16  env geo  1.46013830429189   0.016"

#plot for islands
#trim datasets
bioclim_trim <- bioclim[bioclim[,"island"] == "y",]
community_trim <- community[row.names(bioclim_trim),]
community_trim <- data.frame(community_trim[rowSums(community_trim) > 0,])
#geographical data
geography <- data.frame(bioclim_trim[row.names(community_trim),2:1])
geography_pcnm <- as.data.frame(scores(pcnm(dist(geography))))
#environmental data
env <- env.z[row.names(community_trim),c("rasValue_bio10", "rasValue_bio18",
                                         "rasValue_bio4", "rasValue_bio15", "rasValue_bio2")]
#environmental VS geography
var.part <- varpart(community_trim, env, geography_pcnm)
#print(var.part$part)
#plot venn diagram
plot(var.part, 
     Xnames = c("Bioclimatic dataset", "Geography"), # name the partitions
     bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
title(main = "islands")
#test for significance
anova <- anova.cca(rda(community_trim, geography_pcnm))
print(paste("islands geo no ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
#[1] "islands geo no  10.1884726463826   0.001"
anova <- anova.cca(rda(community_trim, env))
print(paste("islands env no ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
#[1] "islands env no  4.71642319848897   0.001"
anova <- anova.cca(rda(community_trim, geography_pcnm, env))
print(paste("islands geo env ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
#[1] "islands geo env  4.71642319848897   0.001"
anova <- anova.cca(rda(community_trim, env, geography_pcnm))
print(paste("islands env geo ", anova$Variance[1], " ", anova$"Pr(>F)"[1]))
#[1] "islands env geo  3.03894206278714   0.001"

#saved as variation_partitioning_by_bioregion.svg

#########################
###Distance from coast###
#########################

#import community information
load(file = "all_dataset_ps_srs_analysis_bacteria.RData")
community <- t(genera_and_unidentified_counts_tab[-1,])
dim(community)
#[1]  987 1412

#import information about bioregions and bioclimatic variables
bioclim <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)
dim(bioclim)
#[1] 987  20

for (n in sort(unique(bioclim[,"bioregion"]))) {
  
  #trim datasets
  bioclim_trim <- bioclim[bioclim[,"bioregion"] == n,]
  community_trim <- community[row.names(bioclim_trim),]
  community_trim <- community_trim[rowSums(community_trim) > 0,]
  dim(community_trim)
  #[1]   80 1412
  
  #distance from coast
  dist_coast <- bioclim_trim[row.names(community_trim),18]
  
  #get top 100 genera
  top100 <- names(sort(colSums(community_trim), decreasing=TRUE)[1:100])
  
  #get Spearman correlation between taxa and distance from coast
  count = 0
  for (i in c(1:100)) {
    sp <- cor.test(community_trim[,top100[i]], dist_coast, method = "spearman", exact=FALSE)
    if (!(is.na(as.numeric(sp$estimate)))) {
      if (sp$p.value < 0.05) {
        count = count + 1
      }
    }
  }
  print(paste(n, " ", count))
}
# [1] "0   36"
# [1] "1   69"
# [1] "3   41"
# [1] "4   26"
# [1] "6   11"
# [1] "7   44"
# [1] "8   58"
# [1] "9   40"
# [1] "10   67"
# [1] "12   0"
# [1] "16   6"

##############################
###permanova for bioregions###
##############################

setwd("/home/gilda/all")

library("vegan")

#import community information
load(file = "all_dataset_ps_srs_bacteria.RData")
community <- t(genera_proportions_tab_hellinger[-1,samples])
dim(community)
#[1]  988 1445

#import information about bioregions and bioclimatic variables
bioclim <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)[samples,]
dim(bioclim)
#[1] 988  20

#remove bioregion represented by only few samples: keep only 3, 6, 7, 8, 9, 10 and 16
bioclim_trim <- bioclim[bioclim[,"bioregion"] == 0 | bioclim[,"bioregion"] == 3 | bioclim[,"bioregion"] == 6
                        | bioclim[,"bioregion"] == 7 | bioclim[,"bioregion"] == 8
                        | bioclim[,"bioregion"] == 9 | bioclim[,"bioregion"] == 10
                        | bioclim[,"bioregion"] == 16,]
community_trim <- community[row.names(bioclim_trim),]
community_trim <- data.frame(community_trim[rowSums(community_trim) > 0,])
dim(community_trim)
#[1]   963 1445

#distance from coast
bioclim_trim <- bioclim_trim[row.names(community_trim),]
bioclim_trim$bioregion <- as.character(bioclim_trim$bioregion)

#run permanova
adonis2(community_trim ~ bioregion, data=bioclim_trim, distance="bray", permutations=1000)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = community_trim ~ bioregion, data = bioclim_trim, permutations = 1000, distance = "bray")
# Df SumOfSqs      R2      F   Pr(>F)    
# bioregion   6   52.832 0.18972 36.369 0.000999 ***
#   Residual  932  225.644 0.81028                    
# Total     938  278.475 1.00000                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#When 0 is included
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = community_trim ~ bioregion, data = bioclim_trim, permutations = 1000, distance = "bray")
# Df SumOfSqs      R2     F   Pr(>F)    
# bioregion   7   55.995 0.19501 33.05 0.000999 ***
#   Residual  955  231.140 0.80499                   
# Total     962  287.135 1.00000                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#pairwise comparisons
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#pairwise.adonis(community_trim, factors=bioclim_trim$bioregion, p.adjust.m ='bonferroni', sim.method = "bray")
#pairwise comparisons are significant for pvalue < 0.05
pair <- pairwise.adonis(community_trim, factors=bioclim_trim$bioregion, p.adjust.m ='fdr', sim.method = "bray")
#all comparisons are significant for pvalue < 0.01

#setup dataframe
bioregion_R2 <- data.frame(rep(1,7),rep(1,7),rep(1,7),rep(1,7),rep(1,7),rep(1,7),rep(1,7))
row.names(bioregion_R2) <- c("3", "6", "7", "8", "9", "10", "16")
colnames(bioregion_R2) <- c("3", "6", "7", "8", "9", "10", "16")

for (i in c(1:length(pair$pairs))) {
  
  #get pairwise comparison bioregions 
  b1 <- strsplit(pair$pairs[i], " ")[[1]][1]
  b2 <- strsplit(pair$pairs[i], " ")[[1]][3]
  
  #get R2
  bioregion_R2[b1,b2] <- pair$R2[i]
  bioregion_R2[b2,b1] <- pair$R2[i]
}

#plot heatmap
library(gplots)
heatmap.2(as.matrix(bioregion_R2), Rowv=FALSE, Colv=FALSE, density.info="none",
          cellnote=round(as.matrix(bioregion_R2), digits = 2), notecol="white",
          trace="none", breaks=c(0,0.05,0.1,0.15,0.2,0.3,1), 
          col = c("orange", "darkorange1", "darkorange2", "darkorange3", "red2", "white"))
#saved as figure_permanova_pairwise_bioregions.svg and figure_permanova_pairwise_bioregions.png
#dimensions 1000x1000