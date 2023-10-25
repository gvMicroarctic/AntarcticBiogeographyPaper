#############################################################
###Plot tanglegrams for bioclimatic and community datasets###
#############################################################

setwd("/home/gilda/all")

load(file = "all_dataset_ps_srs_bacteria_hell.RData")
#https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html#tanglegram

#import information about bioregions and bioclimatic variables
bioclim <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)[samples,]
dim(bioclim)
#[1] 988  20

#community dataset
library(vegan)
community <- t(genera_proportions_tab_hellinger[-1,samples])
dim(community)
#[1]  988 1445
community_tran <- community[rowSums(community)>0,]
dim(community_tran)
#[1]  988 1445
dd <- vegdist(community_tran, method = "bray")
hc <- hclust(dd, method = "ward.D2")
hcd2 <- as.dendrogram(hc,hang=0.1)

#bioclimatic dataset
bioclim_tran <- bioclim[samples,c(3:12,14)]
bioclim_tran$rasValue_swe[is.na(bioclim_tran$rasValue_swe)] <- 0
dim(bioclim_tran)
#[1] 988  11
dd <- dist(scale(bioclim_tran), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
hcd1 <- as.dendrogram(hc,hang=0.1)

#geography dataset
geography_tran <- bioclim[samples,c(2,1)]
dim(geography_tran)
#[1] 988  2
library(geosphere)
dd0 <- distm(geography_tran, fun = distHaversine) #this function works with decimal degrees
dd1 <- as.dist(dd0)
library("usedist")
dd <- dist_setNames(dd1, samples)
hc <- hclust(dd, method = "ward.D2")
hcd3 <- as.dendrogram(hc,hang=0.1)

library(dendextend)

#calculate correlation between trees
cors <- cor.dendlist(dendlist(d1 = hcd1, d2 = hcd2, d3 = hcd3))
cors
# d1        d2        d3
# d1 1.0000000 0.1974212 0.4292164
# d2 0.1974212 1.0000000 0.1687299
# d3 0.4292164 0.1687299 1.0000000

cor.dendlist(dendlist(d1 = hcd1, d2 = hcd2, d3 = hcd3), method = "baker")
# d1         d2         d3
# d1 1.0000000 0.19274123 0.43825432
# d2 0.1927412 1.00000000 0.08322676
# d3 0.4382543 0.08322676 1.00000000

cor.dendlist(dendlist(d1 = hcd1, d2 = hcd2, d3 = hcd3), method = "common_nodes")
# d1        d2        d3
# d1 1.0000000 0.5331646 0.6992405
# d2 0.5331646 1.0000000 0.5422785
# d3 0.6992405 0.5422785 1.0000000

####create tanglegram

###Entire - bioclimatic vs community

#combine dendrograms
merged <- tanglegram(hcd1, hcd2, sort = TRUE)

#get order labels in dendogram1
order <- merged$dend1 %>% labels
#order metadata rows as in dendogram1
col1 <- data.frame(bioclim[order,])
#color
col1[,1] <- as.character(col1[,1])
col1[,1][col1[,"bioregion"] == "0"] <- "black"
col1[,1][col1[,"bioregion"] == "1"] <- "#8b0001ff"
col1[,1][col1[,"bioregion"] == "3"] <- "#fe0000ff"
col1[,1][col1[,"bioregion"] == "4"] <- "#ffff00ff"
col1[,1][col1[,"bioregion"] == "6"] <- "#8fbc8bff"
col1[,1][col1[,"bioregion"] == "7"] <- "#32cd33ff"
col1[,1][col1[,"bioregion"] == "8"] <- "#006400ff"
col1[,1][col1[,"bioregion"] == "9"] <- "#26a8f6ff"
col1[,1][col1[,"bioregion"] == "10"] <- "#4682b3ff"
col1[,1][col1[,"bioregion"] == "12"] <- "#808080ff"
col1[,1][col1[,"bioregion"] == "16"] <- "#0000ccff"

#get order labels in dendogram2
order <- merged$dend2 %>% labels
#order metadata rows as in dendogram2
col2 <- data.frame(bioclim[order,])
#color
col2[,1] <- as.character(col2[,1])
col2[,1][col2[,"bioregion"] == "0"] <- "black"
col2[,1][col2[,"bioregion"] == "1"] <- "#8b0001ff"
col2[,1][col2[,"bioregion"] == "3"] <- "#fe0000ff"
col2[,1][col2[,"bioregion"] == "4"] <- "#ffff00ff"
col2[,1][col2[,"bioregion"] == "6"] <- "#8fbc8bff"
col2[,1][col2[,"bioregion"] == "7"] <- "#32cd33ff"
col2[,1][col2[,"bioregion"] == "8"] <- "#006400ff"
col2[,1][col2[,"bioregion"] == "9"] <- "#26a8f6ff"
col2[,1][col2[,"bioregion"] == "10"] <- "#4682b3ff"
col2[,1][col2[,"bioregion"] == "12"] <- "#808080ff"
col2[,1][col2[,"bioregion"] == "16"] <- "#0000ccff"

#add shape and color to dendogram objects
hcd1_pl <- merged$dend1 %>% set("labels_cex", 10) %>% set("leaves_pch", 15) %>% set("leaves_cex", 2) %>% set("leaves_col", col1[,1]) %>% hang.dendrogram %>% hang.dendrogram(hang = -1) %>% set("branches_lwd", 1) %>% set("branches_lty", 1)
hcd2_pl <- merged$dend2 %>% set("labels_cex", 10) %>% set("leaves_pch", 15) %>% set("leaves_cex", 2) %>% set("leaves_col", col2[,1]) %>% hang.dendrogram %>% hang.dendrogram(hang = -1) %>% set("branches_lwd", 1) %>% set("branches_lty", 1)

#plot
png("figure_tanglegram_bio.png",
    pointsize=12)
tanglegram(hcd1_pl, hcd2_pl, sort = TRUE, color_lines=col1[,1], lwd = 0.5, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

#get entanglement value
#https://search.r-project.org/CRAN/refmans/dendextend/html/entanglement.html

entanglement(hcd1_pl, hcd2_pl, L =0) #e.g. panality for all the non horizontal lines! So very low score!!!!
#[1] 0.9979757

entanglement(hcd1_pl, hcd2_pl, L =0.5)
#[1] 0.6903678

entanglement(hcd1_pl, hcd2_pl, L =1)
#[1] 0.5134161

entanglement(hcd1_pl, hcd2_pl, L =1.5) #1.5 is the default value
#[1] 0.4001409

entanglement(hcd1_pl, hcd2_pl, L =2)
#[1] 0.320544

#Entanglement number: it ranges between 0 and 1! The higher it is, the worst it is (the lines connecting the two trees are more entangled)

all.equal(hcd1_pl, hcd2_pl, use.edge.length = TRUE)
#[1] "Difference in branch heights -  Mean relative difference: 1.000898"

all.equal(hcd1_pl, hcd2_pl, use.edge.length = FALSE)
#very long list! I am not reporting it here!

all.equal(hcd1_pl, hcd2_pl, use.edge.length = FALSE, use.topology = FALSE)
#[1] TRUE

###Entire - geography vs community

#combine dendrograms
merged <- tanglegram(hcd3, hcd2, sort = TRUE)

#get order labels in dendogram1
order <- merged$dend1 %>% labels
#order metadata rows as in dendogram1
col1 <- data.frame(bioclim[order,])
#color
col1[,1] <- as.character(col1[,1])
col1[,1][col1[,"bioregion"] == "0"] <- "black"
col1[,1][col1[,"bioregion"] == "1"] <- "#8b0001ff"
col1[,1][col1[,"bioregion"] == "3"] <- "#fe0000ff"
col1[,1][col1[,"bioregion"] == "4"] <- "#ffff00ff"
col1[,1][col1[,"bioregion"] == "6"] <- "#8fbc8bff"
col1[,1][col1[,"bioregion"] == "7"] <- "#32cd33ff"
col1[,1][col1[,"bioregion"] == "8"] <- "#006400ff"
col1[,1][col1[,"bioregion"] == "9"] <- "#26a8f6ff"
col1[,1][col1[,"bioregion"] == "10"] <- "#4682b3ff"
col1[,1][col1[,"bioregion"] == "12"] <- "#808080ff"
col1[,1][col1[,"bioregion"] == "16"] <- "#0000ccff"

#get order labels in dendogram2
order <- merged$dend2 %>% labels
#order metadata rows as in dendogram2
col2 <- data.frame(bioclim[order,])
#color
col2[,1] <- as.character(col2[,1])
col2[,1][col2[,"bioregion"] == "0"] <- "black"
col2[,1][col2[,"bioregion"] == "1"] <- "#8b0001ff"
col2[,1][col2[,"bioregion"] == "3"] <- "#fe0000ff"
col2[,1][col2[,"bioregion"] == "4"] <- "#ffff00ff"
col2[,1][col2[,"bioregion"] == "6"] <- "#8fbc8bff"
col2[,1][col2[,"bioregion"] == "7"] <- "#32cd33ff"
col2[,1][col2[,"bioregion"] == "8"] <- "#006400ff"
col2[,1][col2[,"bioregion"] == "9"] <- "#26a8f6ff"
col2[,1][col2[,"bioregion"] == "10"] <- "#4682b3ff"
col2[,1][col2[,"bioregion"] == "12"] <- "#808080ff"
col2[,1][col2[,"bioregion"] == "16"] <- "#0000ccff"

#add shape and color to dendogram objects
hcd1_pl <- merged$dend1 %>% set("labels_cex", 10) %>% set("leaves_pch", 15) %>% set("leaves_cex", 2) %>% set("leaves_col", col1[,1]) %>% hang.dendrogram %>% hang.dendrogram(hang = -1) %>% set("branches_lwd", 1) %>% set("branches_lty", 1)
hcd2_pl <- merged$dend2 %>% set("labels_cex", 10) %>% set("leaves_pch", 15) %>% set("leaves_cex", 2) %>% set("leaves_col", col2[,1]) %>% hang.dendrogram %>% hang.dendrogram(hang = -1) %>% set("branches_lwd", 1) %>% set("branches_lty", 1)

#plot
png("figure_tanglegram_geo.png",
    pointsize=12)
tanglegram(hcd1_pl, hcd2_pl, sort = TRUE, color_lines=col1[,1], lwd = 0.5, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

#get entanglement value
#https://search.r-project.org/CRAN/refmans/dendextend/html/entanglement.html

entanglement(hcd1_pl, hcd2_pl, L =0) #e.g. penalty for all the non horizontal lines! So very low score!!!!
#[1] 0.9979757

entanglement(hcd1_pl, hcd2_pl, L =0.5)
#[1] 0.7301359

entanglement(hcd1_pl, hcd2_pl, L =1)
#[1] 0.5631546

entanglement(hcd1_pl, hcd2_pl, L =1.5) #1.5 is the default value
#[1] 0.4493249

entanglement(hcd1_pl, hcd2_pl, L =2)
#[1] 0.3655436

#Entanglement number: it ranges between 0 and 1! The higher it is, the worst it is (the lines connecting the two trees are more entangled)

all.equal(hcd1_pl, hcd2_pl, use.edge.length = TRUE)
#[1] "Difference in branch heights -  Mean relative difference: 0.9999982"

all.equal(hcd1_pl, hcd2_pl, use.edge.length = FALSE)
#very long list! I am not reporting it here!

all.equal(hcd1_pl, hcd2_pl, use.edge.length = FALSE, use.topology = FALSE)
#[1] TRUE

###Single bioregions

##combine dendrograms
#bioclimatic
merged_bio <- tanglegram(hcd1, hcd2, sort = TRUE)
#get order labels in dendogram1
order <- merged_bio$dend1 %>% labels
#order metadata rows as in dendogram1
col1 <- data.frame(bioclim[order,])
#color
col1[,1] <- as.character(col1[,1])
col1[,1][col1[,"bioregion"] == "0"] <- "black"
col1[,1][col1[,"bioregion"] == "1"] <- "#8b0001ff"
col1[,1][col1[,"bioregion"] == "3"] <- "#fe0000ff"
col1[,1][col1[,"bioregion"] == "4"] <- "#ffff00ff"
col1[,1][col1[,"bioregion"] == "6"] <- "#8fbc8bff"
col1[,1][col1[,"bioregion"] == "7"] <- "#32cd33ff"
col1[,1][col1[,"bioregion"] == "8"] <- "#006400ff"
col1[,1][col1[,"bioregion"] == "9"] <- "#26a8f6ff"
col1[,1][col1[,"bioregion"] == "10"] <- "#4682b3ff"
col1[,1][col1[,"bioregion"] == "12"] <- "#808080ff"
col1[,1][col1[,"bioregion"] == "16"] <- "#0000ccff"
#get order labels in dendogram2
order <- merged_bio$dend2 %>% labels
#order metadata rows as in dendogram2
col2 <- data.frame(bioclim[order,])
#color
col2[,1] <- as.character(col2[,1])
col2[,1][col2[,"bioregion"] == "0"] <- "black"
col2[,1][col2[,"bioregion"] == "1"] <- "#8b0001ff"
col2[,1][col2[,"bioregion"] == "3"] <- "#fe0000ff"
col2[,1][col2[,"bioregion"] == "4"] <- "#ffff00ff"
col2[,1][col2[,"bioregion"] == "6"] <- "#8fbc8bff"
col2[,1][col2[,"bioregion"] == "7"] <- "#32cd33ff"
col2[,1][col2[,"bioregion"] == "8"] <- "#006400ff"
col2[,1][col2[,"bioregion"] == "9"] <- "#26a8f6ff"
col2[,1][col2[,"bioregion"] == "10"] <- "#4682b3ff"
col2[,1][col2[,"bioregion"] == "12"] <- "#808080ff"
col2[,1][col2[,"bioregion"] == "16"] <- "#0000ccff"
hcd1_bio_pl <- merged_bio$dend1 %>% set("labels_cex", 10) %>% set("leaves_pch", 15) %>% set("leaves_cex", 2) %>% set("leaves_col", col1[,1]) %>% hang.dendrogram %>% hang.dendrogram(hang = -1) %>% set("branches_lwd", 1) %>% set("branches_lty", 1)
hcd2_bio_pl <- merged_bio$dend2 %>% set("labels_cex", 10) %>% set("leaves_pch", 15) %>% set("leaves_cex", 2) %>% set("leaves_col", col2[,1]) %>% hang.dendrogram %>% hang.dendrogram(hang = -1) %>% set("branches_lwd", 1) %>% set("branches_lty", 1)
#geography
merged_geo <- tanglegram(hcd2, hcd3, sort = TRUE)
#get order labels in dendogram1
order <- merged_geo$dend1 %>% labels
#order metadata rows as in dendogram1
col1 <- data.frame(bioclim[order,])
#color
col1[,1] <- as.character(col1[,1])
col1[,1][col1[,"bioregion"] == "0"] <- "black"
col1[,1][col1[,"bioregion"] == "1"] <- "#8b0001ff"
col1[,1][col1[,"bioregion"] == "3"] <- "#fe0000ff"
col1[,1][col1[,"bioregion"] == "4"] <- "#ffff00ff"
col1[,1][col1[,"bioregion"] == "6"] <- "#8fbc8bff"
col1[,1][col1[,"bioregion"] == "7"] <- "#32cd33ff"
col1[,1][col1[,"bioregion"] == "8"] <- "#006400ff"
col1[,1][col1[,"bioregion"] == "9"] <- "#26a8f6ff"
col1[,1][col1[,"bioregion"] == "10"] <- "#4682b3ff"
col1[,1][col1[,"bioregion"] == "12"] <- "#808080ff"
col1[,1][col1[,"bioregion"] == "16"] <- "#0000ccff"
#get order labels in dendogram2
order <- merged_geo$dend2 %>% labels
#order metadata rows as in dendogram2
col2 <- data.frame(bioclim[order,])
#color
col2[,1] <- as.character(col2[,1])
col2[,1][col2[,"bioregion"] == "0"] <- "black"
col2[,1][col2[,"bioregion"] == "1"] <- "#8b0001ff"
col2[,1][col2[,"bioregion"] == "3"] <- "#fe0000ff"
col2[,1][col2[,"bioregion"] == "4"] <- "#ffff00ff"
col2[,1][col2[,"bioregion"] == "6"] <- "#8fbc8bff"
col2[,1][col2[,"bioregion"] == "7"] <- "#32cd33ff"
col2[,1][col2[,"bioregion"] == "8"] <- "#006400ff"
col2[,1][col2[,"bioregion"] == "9"] <- "#26a8f6ff"
col2[,1][col2[,"bioregion"] == "10"] <- "#4682b3ff"
col2[,1][col2[,"bioregion"] == "12"] <- "#808080ff"
col2[,1][col2[,"bioregion"] == "16"] <- "#0000ccff"
hcd1_geo_pl <- merged_geo$dend1 %>% set("labels_cex", 10) %>% set("leaves_pch", 15) %>% set("leaves_cex", 2) %>% set("leaves_col", col1[,1]) %>% hang.dendrogram %>% hang.dendrogram(hang = -1) %>% set("branches_lwd", 1) %>% set("branches_lty", 1)
hcd2_geo_pl <- merged_geo$dend2 %>% set("labels_cex", 10) %>% set("leaves_pch", 15) %>% set("leaves_cex", 2) %>% set("leaves_col", col2[,1]) %>% hang.dendrogram %>% hang.dendrogram(hang = -1) %>% set("branches_lwd", 1) %>% set("branches_lty", 1)

##islands

#get order labels in dendrogram for islands vs no_islands
order <- merged_bio$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"island"] == "n"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "0" & line_col[,"island"] == "y"] <- "black"
line_col[,1][line_col[,"bioregion"] == "1" & line_col[,"island"] == "y"] <- "#8b0001ff"
line_col[,1][line_col[,"bioregion"] == "3" & line_col[,"island"] == "y"] <- "#fe0000ff"
line_col[,1][line_col[,"bioregion"] == "4" & line_col[,"island"] == "y"] <- "#ffff00ff"
line_col[,1][line_col[,"bioregion"] == "6" & line_col[,"island"] == "y"] <- "#8fbc8bff"
line_col[,1][line_col[,"bioregion"] == "7" & line_col[,"island"] == "y"] <- "#32cd33ff"
line_col[,1][line_col[,"bioregion"] == "8" & line_col[,"island"] == "y"] <- "#006400ff"
line_col[,1][line_col[,"bioregion"] == "9" & line_col[,"island"] == "y"] <- "#26a8f6ff"
line_col[,1][line_col[,"bioregion"] == "10" & line_col[,"island"] == "y"] <- "#4682b3ff"
line_col[,1][line_col[,"bioregion"] == "12" & line_col[,"island"] == "y"] <- "#808080ff"
line_col[,1][line_col[,"bioregion"] == "16" & line_col[,"island"] == "y"] <- "#0000ccff"
#plot
png("figure_tanglegram_bio_island.png",
    pointsize=12)
tanglegram(hcd1_bio_pl, hcd2_bio_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

#get order labels in dendrogram for islands vs no_islands
order <- merged_geo$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"island"] == "n"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "0" & line_col[,"island"] == "y"] <- "black"
line_col[,1][line_col[,"bioregion"] == "1" & line_col[,"island"] == "y"] <- "#8b0001ff"
line_col[,1][line_col[,"bioregion"] == "3" & line_col[,"island"] == "y"] <- "#fe0000ff"
line_col[,1][line_col[,"bioregion"] == "4" & line_col[,"island"] == "y"] <- "#ffff00ff"
line_col[,1][line_col[,"bioregion"] == "6" & line_col[,"island"] == "y"] <- "#8fbc8bff"
line_col[,1][line_col[,"bioregion"] == "7" & line_col[,"island"] == "y"] <- "#32cd33ff"
line_col[,1][line_col[,"bioregion"] == "8" & line_col[,"island"] == "y"] <- "#006400ff"
line_col[,1][line_col[,"bioregion"] == "9" & line_col[,"island"] == "y"] <- "#26a8f6ff"
line_col[,1][line_col[,"bioregion"] == "10" & line_col[,"island"] == "y"] <- "#4682b3ff"
line_col[,1][line_col[,"bioregion"] == "12" & line_col[,"island"] == "y"] <- "#808080ff"
line_col[,1][line_col[,"bioregion"] == "16" & line_col[,"island"] == "y"] <- "#0000ccff"
#plot
png("figure_tanglegram_geo_island.png",
    pointsize=12)
tanglegram(hcd1_geo_pl, hcd2_geo_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

##bioregion 6

#get order labels in dendrogram for bioregion 6
order <- merged_bio$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#8fbc8bff"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#00000000"
#plot
png("figure_tanglegram_bio_bioregion6.png",
    pointsize=12)
tanglegram(hcd1_bio_pl, hcd2_bio_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

#get order labels in dendrogram for bioregion 6
order <- merged_geo$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#8fbc8bff"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#00000000"
#plot
png("figure_tanglegram_geo_bioregion6.png",
    pointsize=12)
tanglegram(hcd1_geo_pl, hcd2_geo_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

##bioregion 7

#get order labels in dendrogram for bioregion 7
order <- merged_bio$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#32cd33ff"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#00000000"
#plot
png("figure_tanglegram_bio_bioregion7.png",
    pointsize=12)
tanglegram(hcd1_bio_pl, hcd2_bio_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

#get order labels in dendrogram for bioregion 7
order <- merged_geo$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#32cd33ff"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#00000000"
#plot
png("figure_tanglegram_geo_bioregion7.png",
    pointsize=12)
tanglegram(hcd1_geo_pl, hcd2_geo_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

##bioregion 8

#get order labels in dendrogram for bioregion 8
order <- merged_bio$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#006400ff"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#00000000"
#plot
png("figure_tanglegram_bio_bioregion8.png",
    pointsize=12)
tanglegram(hcd1_bio_pl, hcd2_bio_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

#get order labels in dendrogram for bioregion 8
order <- merged_geo$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#006400ff"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#00000000"
#plot
png("figure_tanglegram_geo_bioregion8.png",
    pointsize=12)
tanglegram(hcd1_geo_pl, hcd2_geo_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()


##bioregion 9

#get order labels in dendrogram for bioregion 9
order <- merged_bio$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#26a8f6ff"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#00000000"
#plot
png("figure_tanglegram_bio_bioregion9.png",
    pointsize=12)
tanglegram(hcd1_bio_pl, hcd2_bio_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

#get order labels in dendrogram for bioregion 9
order <- merged_geo$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#26a8f6ff"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#00000000"
#plot
png("figure_tanglegram_geo_bioregion9.png",
    pointsize=12)
tanglegram(hcd1_geo_pl, hcd2_geo_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

##bioregion 10

#get order labels in dendrogram for bioregion 10
order <- merged_bio$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#4682b3ff"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#00000000"
#plot
png("figure_tanglegram_bio_bioregion10.png",
    pointsize=12)
tanglegram(hcd1_bio_pl, hcd2_bio_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

#get order labels in dendrogram for bioregion 10
order <- merged_geo$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#4682b3ff"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#00000000"
#plot
png("figure_tanglegram_geo_bioregion10.png",
    pointsize=12)
tanglegram(hcd1_geo_pl, hcd2_geo_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

##bioregion 16

#get order labels in dendrogram for bioregion 16
order <- merged_bio$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#0000ccff"
#plot
png("figure_tanglegram_bio_bioregion16.png",
    pointsize=12)
tanglegram(hcd1_bio_pl, hcd2_bio_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()

#get order labels in dendrogram for bioregion 16
order <- merged_geo$dend1 %>% labels
#order metadata rows as in dendogram1
line_col <- data.frame(bioclim[order,])
#color
line_col[,1] <- as.character(line_col[,1])
line_col[,1][line_col[,"bioregion"] == "0"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "1"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "3"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "4"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "6"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "7"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "8"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "9"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "10"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "12"] <- "#00000000"
line_col[,1][line_col[,"bioregion"] == "16"] <- "#0000ccff"
#plot
png("figure_tanglegram_geo_bioregion16.png",
    pointsize=12)
tanglegram(hcd1_geo_pl, hcd2_geo_pl, sort = TRUE, color_lines=line_col[,1], lwd=1, common_subtrees_color_lines = TRUE, fast = TRUE)
dev.off()
