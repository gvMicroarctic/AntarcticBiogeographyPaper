#########################################
###Table with mean, min and max values###
#########################################

setwd("/home/gilda/all")

#import information
bioclim0 <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)
dim(bioclim0)
#[1] 990  20

#elevation_distance
metadata_el_dist <- read.csv("distance_altitude.csv", row.names = 1)
dim(metadata_el_dist)
#[1] 990   5
ƒ
#merge
row.names(bioclim0) == row.names(metadata_el_dist)
bioclim1 <- cbind(bioclim0[,c(3:12,14)], metadata_el_dist[,c(3:5)], bioregion=bioclim0[,15])
bioclim <- bioclim1[row.names(bioclim1) != "ERX1598809" & row.names(bioclim1) != "ERX1598810",]
bioclim$rasValue_swe[is.na(bioclim$rasValue_swe)] <- 0
dim(bioclim)
#[1] 988  15
head(bioclim)

#table
bioclimatic_table <- data.frame()
row = 0
for (i in c(0,1,3,4,6,7,8,9,10,12,16)) {
  print(i)
  for (j in c(1:14)) {
    print(j)
    row = row + 1;
    mean <- mean(bioclim[bioclim[,"bioregion"] == i,j])
    sd <- sd(bioclim[bioclim[,"bioregion"] == i,j])
    min <- min(bioclim[bioclim[,"bioregion"] == i,j])
    max <- max(bioclim[bioclim[,"bioregion"] == i,j])
    
    bioclimatic_table[row,1] <- i
    bioclimatic_table[row,2] <- "0"
    bioclimatic_table[row,3] <- mean
    bioclimatic_table[row,4] <- sd
    bioclimatic_table[row,5] <- min
    bioclimatic_table[row,6] <- max
    
  }
}

colnames(bioclimatic_table) <- c("Bioregion", "Bioclimatic variable", "Mean", "Standard deviation", "Minimum", "Maximum")
bioclimatic_table[,2] <- rep(colnames(bioclim[,c(1:14)]),11)
bioclimatic_table
write.csv(bioclimatic_table, "bioclimatic_table.csv")

#save all for SI
#merge
row.names(bioclim0) == row.names(metadata_el_dist)
bioclim1 <- cbind(bioclim0[,c(3:12,14)], metadata_el_dist[,c(3:5)], bioregion=bioclim0[,c(15:17)])
bioclim_all <- bioclim1[row.names(bioclim1) != "ERX1598809" & row.names(bioclim1) != "ERX1598810",]
bioclim_all$rasValue_swe[is.na(bioclim_all$rasValue_swe)] <- 0
dim(bioclim_all)
#[1] 988  17
write.csv(bioclim_all, "bioclimatic_table_entire.csv")

bioclim_no_0 <- bioclim[,]
library("dplyr")
bioclim_no_0 <- bioclim %>%  filter(bioregion!='0')
dim(bioclim_no_0)
#[1] 964  15

#create boxplot
library("ggplot2")
bio1 <- ggplot(bioclim, aes(x=rasValue_bio1, y=factor(bioclim$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("bio1")

bio2 <- ggplot(bioclim, aes(x=rasValue_bio2, y=factor(bioclim$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("bio2")

bio4 <- ggplot(bioclim, aes(x=rasValue_bio4, y=factor(bioclim$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("bio4")

bio10 <- ggplot(bioclim, aes(x=rasValue_bio10, y=factor(bioclim$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("bio10")

bio12 <- ggplot(bioclim, aes(x=rasValue_bio12, y=factor(bioclim$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("bio12")

bio15 <- ggplot(bioclim, aes(x=rasValue_bio15, y=factor(bioclim$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("bio15")

bio18 <- ggplot(bioclim, aes(x=rasValue_bio18, y=factor(bioclim$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("bio18")

bio_swe <- ggplot(bioclim, aes(x=rasValue_swe, y=factor(bioclim$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("swe")

bio_el <- ggplot(bioclim_no_0, aes(x=elevation, y=factor(bioclim_no_0$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim_no_0$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("elevation")

bio_dc <- ggplot(bioclim_no_0, aes(x=dist_coast, y=factor(bioclim_no_0$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim_no_0$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("dist_coast")

bio_do <- ggplot(bioclim_no_0, aes(x=dist_ocean, y=factor(bioclim_no_0$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")), fill=factor(bioclim_no_0$bioregion, levels = c("1","3","4","6","7","8","9","10","12","16")))) + 
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Bioregion") + xlab("dist_ocean")


library(gridExtra)
svg(filename="bioclimatic_boxplot.svg", 
    width=11, 
    height=12, 
    pointsize=12)
grid.arrange(bio1, bio2, bio4, bio10, bio12, bio15, bio18, bio_swe, bio_el, bio_dc, bio_do,  nrow = 4)
dev.off()

#look at correlations
aov <- aov(rasValue_bio1 ~ as.character(bioregion), data=bioclim)
summary(aov)
aov <- aov(rasValue_bio2 ~ as.character(bioregion), data=bioclim)
summary(aov)
aov <- aov(rasValue_bio4 ~ as.character(bioregion), data=bioclim)
summary(aov)
aov <- aov(rasValue_bio10 ~ as.character(bioregion), data=bioclim)
summary(aov)
aov <- aov(rasValue_bio12 ~ as.character(bioregion), data=bioclim)
summary(aov)
aov <- aov(rasValue_bio15 ~ as.character(bioregion), data=bioclim)
summary(aov)
aov <- aov(rasValue_bio18 ~ as.character(bioregion), data=bioclim)
summary(aov)
aov <- aov(elevation ~ as.character(bioregion), data=bioclim)
summary(aov)
aov <- aov(dist_coast ~ as.character(bioregion), data=bioclim)
summary(aov)
aov <- aov(dist_ocean ~ as.character(bioregion), data=bioclim)
summary(aov)
##all significant

bioclim_trim <- bioclim[bioclim[,"bioregion"] == 3 | bioclim[,"bioregion"] == 6
                        | bioclim[,"bioregion"] == 7 | bioclim[,"bioregion"] == 8
                        | bioclim[,"bioregion"] == 9 | bioclim[,"bioregion"] == 10
                        | bioclim[,"bioregion"] == 16,]

#look at correlations
aov <- aov(rasValue_bio1 ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(rasValue_bio2 ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(rasValue_bio4 ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(rasValue_bio5 ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(rasValue_bio10 ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(rasValue_bio12 ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(rasValue_bio14 ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(rasValue_bio15 ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(rasValue_bio17 ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(rasValue_bio18 ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(rasValue_swe ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(elevation ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(dist_coast ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
aov <- aov(dist_ocean ~ as.character(bioregion), data=bioclim_trim)
summary(aov)
##all significant except from dist_coast