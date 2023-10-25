##############################
###Alpha diversity - genera###
##############################

setwd("/home/gilda/all")
load(file = "all_dataset_ps_srs_bacteria.RData")

library(vegan)

df <- t(genera_and_unidentified_counts_tab[-1,samples])
dim(df)
#[1]    988 1445

sha <- diversity(df, index = "shannon", MARGIN = 1, base = exp(1))
spe <- specnumber(df)

#get range for richness
min(spe)
#[1] 3
max(spe)
#[1] 278

metadata <- bioclim
dim(metadata) 
#[1] 988   20

#anosim
anosim(sha, metadata$bioregion, distance = "bray", permutations = 1000)
anosim(sha, metadata$bioregion, distance = "eucl", permutations = 1000)
anosim(spe, metadata$bioregion, distance = "bray", permutations = 1000)
anosim(spe, metadata$bioregion, distance = "eucl", permutations = 1000)
#Euclidean distances are better for alpha diversity indices

#Get the same data for unknowns
df_un <- as.matrix(read.csv("antarctic_soil_genus_taxonomy_rel_abundance.csv",check.names=FALSE))
un <- df_un[,samples]

#merge information
metadata_diversity <- cbind(sim, sha, spe, as.numeric(un[1,]), metadata)
colnames(metadata_diversity)[4] <- "un"

min(as.numeric(un[1,]))
#[1] 0
max(as.numeric(un[1,]))
#[1] 84.29579

#format
row.names(df) == row.names(metadata)
metadata$bioregion <- factor(metadata$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16"))

library("ggplot2")
#Shannon index
p_sha <- ggplot(metadata_diversity, aes(x=sha, y=factor(metadata$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")), fill=factor(metadata$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) +  
  geom_violin() + geom_boxplot(width=0.1) + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Bioregion") + xlab("Shannon index")

#Richness
p_spe <- ggplot(metadata_diversity, aes(x=spe, y=factor(metadata$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")), fill=factor(metadata$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) + 
  geom_violin() + geom_boxplot(width=0.1) + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Bioregion") + xlab("Richness")

#Unknowns
p_un <- ggplot(metadata_diversity, aes(x=un, y=factor(metadata$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")), fill=factor(metadata$bioregion, levels = c("0","1","3","4","6","7","8","9","10","12","16")))) + 
  geom_violin() + geom_boxplot(width=0.1) + theme_bw() +
  scale_fill_manual(values=c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8fbc8bff", "#32cd33ff", "#006400ff", "#26a8f6ff", "#4682b3ff", "#808080ff", "#0000ccff")) +
  coord_flip() + 
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Bioregion") + xlab("Unknown")

library(gridExtra)
library("Cairo")
CairoSVG(file="violin_plot_alpha_genus.svg", 
         width=11, 
         height=14, 
         pointsize=12)
grid.arrange(p_spe, p_sha, p_un, nrow = 3)
dev.off()

#test for correlation number of genera vs unknowns
cor.test(spe, as.numeric(un[1,]), method="pearson")
# Pearson's product-moment correlation
# data:  spe and as.numeric(un[1, ])
# t = -5.1684, df = 986, p-value = 2.858e-07
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2225239 -0.1010644
# sample estimates:
#        cor
# -0.1624093

#anova for Shannon diversity
sha.aov <- aov(sha ~ bioregion, data=metadata)
summary(sha.aov)
# Df Sum Sq Mean Sq F value Pr(>F)
# bioregion    10  294.0  29.404   66.71 <2e-16 ***
#   Residuals   977  430.6   0.441
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
sha.aov.t <- TukeyHSD(sha.aov)

#color dataset
dataset_col <- rbind(c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"))
row.names(dataset_col) <- c("0", "1", "3", "4", "6", "7", "8", "9", "10", "12", "16")
colnames(dataset_col) <- c("0", "1", "3", "4", "6", "7", "8", "9", "10", "12", "16")

#change colors in relation to pvalue
dt.div <- sha.aov.t$bioregion
for(i in 1:length(row.names(dt.div))) {
  pair <- row.names(dt.div)[i]
  n1 <- strsplit(pair, split = "-")[[1]][1]
  n2 <- strsplit(pair, split = "-")[[1]][2]
  #id p-value > 0.05 change color to "white"
  pvalue <- dt.div[i,4]
  if (pvalue >= 0.01) {
    dataset_col[n1,n2] <- "white"
    dataset_col[n2,n1] <- "white"
  } else {
    dataset_col[n1,n2] <- "black"
    dataset_col[n2,n1] <- "black"
  }
}

#figure reporting pvalue
dataset <- data.frame(rbind(c(0, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(1, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(3, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(4, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(6, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(7, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(8, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(9, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(10, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(12, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(16, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16)))
row.names(dataset) <- c("y0", "y1", "y3", "y4", "y6", "y7", "y8", "y9", "y10", "y12", "y16")
colnames(dataset) <- c("bio", "x0", "x1", "x3", "x4", "x6", "x7", "x8", "x9", "x10", "x12", "x16")

library("Cairo")
CairoSVG(file="alpha_sha_pairwise_anova.svg",
         width=8,
         height=6,
         pointsize=12)
dotchart(dataset$x0, pch = 21, labels = dataset$bio, bg = as.vector(dataset_col[1,]), pt.cex = 4, xlim = c(0, 16))
points(dataset$x1, 1:nrow(dataset), col = as.vector(dataset_col[2,]), pch = 19, cex = 4)
points(dataset$x3, 1:nrow(dataset), col = as.vector(dataset_col[3,]), pch = 19, cex = 4)
points(dataset$x4, 1:nrow(dataset), col = as.vector(dataset_col[4,]), pch = 19, cex = 4)
points(dataset$x6, 1:nrow(dataset), col = as.vector(dataset_col[5,]), pch = 19, cex = 4)
points(dataset$x7, 1:nrow(dataset), col = as.vector(dataset_col[6,]), pch = 19, cex = 4)
points(dataset$x8, 1:nrow(dataset), col = as.vector(dataset_col[7,]), pch = 19, cex = 4)
points(dataset$x9, 1:nrow(dataset), col = as.vector(dataset_col[8,]), pch = 19, cex = 4)
points(dataset$x10, 1:nrow(dataset), col = as.vector(dataset_col[9,]), pch = 19, cex = 4)
points(dataset$x12, 1:nrow(dataset), col = as.vector(dataset_col[10,]), pch = 19, cex = 4)
points(dataset$x16, 1:nrow(dataset), col = as.vector(dataset_col[11,]), pch = 19, cex = 4)
dev.off()

#anova for Richness diversity
spe.aov <- aov(spe ~ bioregion, data=metadata)
summary(spe.aov)
# Df  Sum Sq Mean Sq F value Pr(>F)
# bioregion    10 1131083  113108   56.84 <2e-16 ***
#   Residuals   977 1944115    1990
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
spe.aov.t <- TukeyHSD(spe.aov)

#color dataset
dataset_col <- rbind(c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"))
row.names(dataset_col) <- c("0", "1", "3", "4", "6", "7", "8", "9", "10", "12", "16")
colnames(dataset_col) <- c("0", "1", "3", "4", "6", "7", "8", "9", "10", "12", "16")

#change colors in relation to pvalue
dt.div <- spe.aov.t$bioregion
for(i in 1:length(row.names(dt.div))) {
  pair <- row.names(dt.div)[i]
  n1 <- strsplit(pair, split = "-")[[1]][1]
  n2 <- strsplit(pair, split = "-")[[1]][2]
  #id p-value > 0.05 change color to "white"
  pvalue <- dt.div[i,4]
  if (pvalue >= 0.01) {
    dataset_col[n1,n2] <- "white"
    dataset_col[n2,n1] <- "white"
  } else {
    dataset_col[n1,n2] <- "black"
    dataset_col[n2,n1] <- "black"
  }
}

#figure reporting pvalue
dataset <- data.frame(rbind(c(0, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(1, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(3, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(4, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(6, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(7, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(8, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(9, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(10, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(12, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(16, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16)))
row.names(dataset) <- c("y0", "y1", "y3", "y4", "y6", "y7", "y8", "y9", "y10", "y12", "y16")
colnames(dataset) <- c("bio", "x0", "x1", "x3", "x4", "x6", "x7", "x8", "x9", "x10", "x12", "x16")

library("Cairo")
CairoSVG(file="alpha_spe_pairwise_anova.svg",
         width=8,
         height=6,
         pointsize=12)
dotchart(dataset$x0, pch = 21, labels = dataset$bio, bg = as.vector(dataset_col[1,]), pt.cex = 4, xlim = c(0, 16))
points(dataset$x1, 1:nrow(dataset), col = as.vector(dataset_col[2,]), pch = 19, cex = 4)
points(dataset$x3, 1:nrow(dataset), col = as.vector(dataset_col[3,]), pch = 19, cex = 4)
points(dataset$x4, 1:nrow(dataset), col = as.vector(dataset_col[4,]), pch = 19, cex = 4)
points(dataset$x6, 1:nrow(dataset), col = as.vector(dataset_col[5,]), pch = 19, cex = 4)
points(dataset$x7, 1:nrow(dataset), col = as.vector(dataset_col[6,]), pch = 19, cex = 4)
points(dataset$x8, 1:nrow(dataset), col = as.vector(dataset_col[7,]), pch = 19, cex = 4)
points(dataset$x9, 1:nrow(dataset), col = as.vector(dataset_col[8,]), pch = 19, cex = 4)
points(dataset$x10, 1:nrow(dataset), col = as.vector(dataset_col[9,]), pch = 19, cex = 4)
points(dataset$x12, 1:nrow(dataset), col = as.vector(dataset_col[10,]), pch = 19, cex = 4)
points(dataset$x16, 1:nrow(dataset), col = as.vector(dataset_col[11,]), pch = 19, cex = 4)
dev.off()

#anova for unknown genera
un.aov <- aov(un[1,] ~ bioregion, data=metadata)
summary(un.aov)
# Df Sum Sq Mean Sq F value Pr(>F)
# bioregion    10  30245  3024.5   12.67 <2e-16 ***
#   Residuals   977 233281   238.8
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
un.aov.t <- TukeyHSD(un.aov)

#color dataset
dataset_col <- rbind(c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"), 
                     c("black", "#8b0001ff", "#fe0000ff", "#ffff00ff", "#8ebc8bff", "#32cd33ff", "#006400ff", "#28a8f7ff", "#4682b4ff", "#808080ff", "#0000ccff"))
row.names(dataset_col) <- c("0", "1", "3", "4", "6", "7", "8", "9", "10", "12", "16")
colnames(dataset_col) <- c("0", "1", "3", "4", "6", "7", "8", "9", "10", "12", "16")

#change colors in relation to pvalue
dt.div <- un.aov.t$bioregion
for(i in 1:length(row.names(dt.div))) {
  pair <- row.names(dt.div)[i]
  n1 <- strsplit(pair, split = "-")[[1]][1]
  n2 <- strsplit(pair, split = "-")[[1]][2]
  #id p-value > 0.05 change color to "white"
  pvalue <- dt.div[i,4]
  if (pvalue >= 0.01) {
    dataset_col[n1,n2] <- "white"
    dataset_col[n2,n1] <- "white"
  } else {
    dataset_col[n1,n2] <- "black"
    dataset_col[n2,n1] <- "black"
  }
}

#figure reporting pvalue
dataset <- data.frame(rbind(c(0, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(1, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(3, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(4, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(6, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(7, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(8, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(9, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(10, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16), c(12, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16),
                            c(16, 0, 1, 3, 4, 6, 7, 8, 9, 10, 12, 16)))
row.names(dataset) <- c("y0", "y1", "y3", "y4", "y6", "y7", "y8", "y9", "y10", "y12", "y16")
colnames(dataset) <- c("bio", "x0", "x1", "x3", "x4", "x6", "x7", "x8", "x9", "x10", "x12", "x16")

library("Cairo")
CairoSVG(file="alpha_un_pairwise_anova.svg",
         width=8,
         height=6,
         pointsize=12)
dotchart(dataset$x0, pch = 21, labels = dataset$bio, bg = as.vector(dataset_col[1,]), pt.cex = 4, xlim = c(0, 16))
points(dataset$x1, 1:nrow(dataset), col = as.vector(dataset_col[2,]), pch = 19, cex = 4)
points(dataset$x3, 1:nrow(dataset), col = as.vector(dataset_col[3,]), pch = 19, cex = 4)
points(dataset$x4, 1:nrow(dataset), col = as.vector(dataset_col[4,]), pch = 19, cex = 4)
points(dataset$x6, 1:nrow(dataset), col = as.vector(dataset_col[5,]), pch = 19, cex = 4)
points(dataset$x7, 1:nrow(dataset), col = as.vector(dataset_col[6,]), pch = 19, cex = 4)
points(dataset$x8, 1:nrow(dataset), col = as.vector(dataset_col[7,]), pch = 19, cex = 4)
points(dataset$x9, 1:nrow(dataset), col = as.vector(dataset_col[8,]), pch = 19, cex = 4)
points(dataset$x10, 1:nrow(dataset), col = as.vector(dataset_col[9,]), pch = 19, cex = 4)
points(dataset$x12, 1:nrow(dataset), col = as.vector(dataset_col[10,]), pch = 19, cex = 4)
points(dataset$x16, 1:nrow(dataset), col = as.vector(dataset_col[11,]), pch = 19, cex = 4)
dev.off()

############################################################
###Biplots with alpha diversity and bioclimatic variables###
############################################################

setwd("/home/gilda/all")
load(file = "all_dataset_ps_srs_bacteria.RData")

#alpha diversity
library(vegan)

df <- t(genera_and_unidentified_counts_tab[-1,samples])
dim(df)
#[1]    988 1445

sha <- diversity(df, index = "shannon", MARGIN = 1, base = exp(1))
spe <- specnumber(df)

#bioclimatic variables
#import information
bioclim0 <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)
dim(bioclim0)
#[1] 990  20
#elevation_distance
metadata_el_dist <- read.csv("distance_altitude.csv", row.names = 1)
dim(metadata_el_dist)
#[1] 990   5
#merge
row.names(bioclim0) == row.names(metadata_el_dist)
bioclim1 <- cbind(bioclim0[,c(3:12,14)], metadata_el_dist[,c(3:5)], bioregion=bioclim0[,15])
bioclim <- bioclim1[row.names(bioclim1) != "ERX1598809" & row.names(bioclim1) != "ERX1598810",]
bioclim$rasValue_swe[is.na(bioclim$rasValue_swe)] <- 0
dim(bioclim)
#[1] 988  15
head(bioclim)

#check data in same order
row.names(df) == row.names(bioclim)
dt <- cbind(bioclim, spe, sha)

#install.packages("ggpubr")
library("ggpubr")

#richness
p1 <- ggscatter(dt, x = "rasValue_bio1", y = "spe", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = TRUE, cor.method = "pearson",
                ylab = "richness", xlab = "BIO1")

p2 <- ggscatter(dt, x = "rasValue_bio2", y = "spe", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = TRUE, cor.method = "pearson",
                ylab = "richness", xlab = "BIO2")

p3 <- ggscatter(dt, x = "rasValue_bio4", y = "spe", 
                add = "reg.line", conf.int = TRUE, size = 1, 
                cor.coef = TRUE, cor.method = "pearson",
                ylab = "richness", xlab = "BIO4")

p4 <- ggscatter(dt, x = "rasValue_bio5", y = "spe", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = TRUE, cor.method = "pearson",
                ylab = "richness", xlab = "BIO5")

p5 <- ggscatter(dt, x = "rasValue_bio10", y = "spe", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = TRUE, cor.method = "pearson",
                ylab = "richness", xlab = "BIO10")

p6 <- ggscatter(dt, x = "rasValue_bio12", y = "spe", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = TRUE, cor.method = "pearson",
                ylab = "richness", xlab = "BIO12")

p7 <- ggscatter(dt, x = "rasValue_bio14", y = "spe", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = TRUE, cor.method = "pearson",
                ylab = "richness", xlab = "BIO14")

p8 <- ggscatter(dt, x = "rasValue_bio15", y = "spe", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = TRUE, cor.method = "pearson",
                ylab = "richness", xlab = "BIO15")

p9 <- ggscatter(dt, x = "rasValue_bio17", y = "spe", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = TRUE, cor.method = "pearson",
                ylab = "richness", xlab = "BIO17")

p10 <- ggscatter(dt, x = "rasValue_bio18", y = "spe", 
                 add = "reg.line", conf.int = TRUE, size = 1,
                 cor.coef = TRUE, cor.method = "pearson",
                 ylab = "richness", xlab = "BIO18")

p11 <- ggscatter(dt, x = "rasValue_swe", y = "spe", 
                 add = "reg.line", conf.int = TRUE, size = 1,
                 cor.coef = TRUE, cor.method = "pearson",
                 ylab = "richness", xlab = "SWE")

p12 <- ggscatter(dt, x = "dist_ocean", y = "spe", 
                 add = "reg.line", conf.int = TRUE, size = 1,
                 cor.coef = TRUE, cor.method = "pearson",
                 ylab = "richness", xlab = "Distance to ocean")

p13 <- ggscatter(dt, x = "elevation", y = "spe", 
                 add = "reg.line", conf.int = TRUE, size = 1,
                 cor.coef = TRUE, cor.method = "pearson",
                 ylab = "richness", xlab = "Elevation")

library("gridExtra")
library(svglite)
svglite("biplots_richness_biocl.svg", width=11, height=14)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8,
             p9, p10, p11, p12, p13, nrow = 5)
dev.off()

#shannon
p1 <- ggscatter(dt, x = "rasValue_bio1", y = "sha", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = FALSE, cor.method = "pearson",
                ylab = "Shannon diversity", xlab = "BIO1") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p2 <- ggscatter(dt, x = "rasValue_bio2", y = "sha", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = FALSE, cor.method = "pearson",
                ylab = "Shannon diversity", xlab = "BIO2") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p3 <- ggscatter(dt, x = "rasValue_bio4", y = "sha", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = FALSE, cor.method = "pearson",
                ylab = "Shannon diversity", xlab = "BIO4") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p4 <- ggscatter(dt, x = "rasValue_bio5", y = "sha", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = FALSE, cor.method = "pearson",
                ylab = "Shannon diversity", xlab = "BIO5") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p5 <- ggscatter(dt, x = "rasValue_bio10", y = "sha", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = FALSE, cor.method = "pearson",
                ylab = "Shannon diversity", xlab = "BIO10") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p6 <- ggscatter(dt, x = "rasValue_bio12", y = "sha", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = FALSE, cor.method = "pearson",
                ylab = "Shannon diversity", xlab = "BIO12") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p7 <- ggscatter(dt, x = "rasValue_bio14", y = "sha", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = FALSE, cor.method = "pearson",
                ylab = "Shannon diversity", xlab = "BIO14") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p8 <- ggscatter(dt, x = "rasValue_bio15", y = "sha", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = FALSE, cor.method = "pearson",
                ylab = "Shannon diversity", xlab = "BIO15") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p9 <- ggscatter(dt, x = "rasValue_bio17", y = "sha", 
                add = "reg.line", conf.int = TRUE, size = 1,
                cor.coef = FALSE, cor.method = "pearson",
                ylab = "Shannon diversity", xlab = "BIO17") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p10 <- ggscatter(dt, x = "rasValue_bio18", y = "sha", 
                 add = "reg.line", conf.int = TRUE, size = 1,
                 cor.coef = FALSE, cor.method = "pearson",
                 ylab = "Shannon diversity", xlab = "BIO18") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p11 <- ggscatter(dt, x = "rasValue_swe", y = "sha", 
                 add = "reg.line", conf.int = TRUE, size = 1,
                 cor.coef = FALSE, cor.method = "pearson",
                 ylab = "Shannon diversity", xlab = "SWE") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p12 <- ggscatter(dt, x = "dist_ocean", y = "sha", 
                 add = "reg.line", conf.int = TRUE, size = 1,
                 cor.coef = FALSE, cor.method = "pearson",
                 ylab = "Shannon diversity", xlab = "Distance to ocean") + stat_cor(label.x.npc = "middle", label.y = 0.2)

p13 <- ggscatter(dt, x = "elevation", y = "sha", 
                 add = "reg.line", conf.int = TRUE, size = 1,
                 cor.coef = FALSE, cor.method = "pearson",
                 ylab = "Shannon diversity", xlab = "Elevation") + stat_cor(label.x.npc = "middle", label.y = 0.2)

svglite("biplots_shannon_biocl.svg", width=11, height=14)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8,
             p9, p10, p11, p12, p13, nrow = 5)
dev.off()