###########
###dbRDA###
###########

setwd("/home/gilda/all/")
load(file = "all_dataset_ps_srs_bacteria.RData")

#hellinger transformed genera dataset
genera_counts_tab_trim <- data.frame(t(genera_proportions_tab_hellinger))[samples,-1]
dim(genera_counts_tab_trim)
#[1]  988 1445

#bioclimatic variables
bioclim <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names=1)[samples,]
dim(bioclim)
#[1] 988  20
#scd: Snow cover days; I am removing it because only two values 365 and 0 and it looks like it is not accurate as it gives 365 even in areas were the snow coverage is not all year around
#swe: Snow water equivalent; it has NA values in samples with scd = 0 so I am changing NA to 0
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

#run RDA
dbrda.upr <- capscale(genera_counts_tab_trim ~ ., data = data.frame(env.z), dist="bray") # upper model limit (the "full" model; it takes in account all variables)
dbrda.lwr <- capscale(genera_counts_tab_trim ~ 1, data = data.frame(env.z), dist="bray") # lower model limit

#run the ordistep model to identify variables with explanatory power
fwd.sel <- ordiR2step(dbrda.lwr, #lower model limit
                      scope = formula(dbrda.upr), #upper model limit
                      direction = "forward",
                      R2scope = TRUE, #can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) #change to TRUE to see the selection process!

#check the new model with forward-selected variables
fwd.sel$call
# capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 +
#            rasValue_bio18 + rasValue_bio5 + rasValue_bio15 + rasValue_bio2 +
#            rasValue_bio4 + rasValue_bio1 + rasValue_bio14 + rasValue_bio17 +
#            rasValue_bio12, data = data.frame(env.z), distance = "bray")

#get statistical significance of the model
anova.cca(fwd.sel, permutations = 1000) #same results as anova()
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 1000
# 
# Model: capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 + rasValue_bio18 + rasValue_bio5 + rasValue_bio15 + rasValue_bio2 + rasValue_bio4 + rasValue_bio1 + rasValue_bio14 + rasValue_bio17 + rasValue_bio12, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs      F   Pr(>F)
# Model     10   50.337 16.774 0.000999 ***
#   Residual 977  293.185
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#get statistical significance of the model for each variable
anova.cca(fwd.sel, step = 1000, by = "term")
# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 + rasValue_bio18 + rasValue_bio5 + rasValue_bio15 + rasValue_bio2 + rasValue_bio4 + rasValue_bio1 + rasValue_bio14 + rasValue_bio17 + rasValue_bio12, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs       F Pr(>F)
# rasValue_bio10   1   13.904 46.3321  0.001 ***
#   rasValue_bio18   1    7.815 26.0415  0.001 ***
#   rasValue_bio5    1    6.508 21.6865  0.001 ***
#   rasValue_bio15   1    4.902 16.3357  0.001 ***
#   rasValue_bio2    1    4.516 15.0506  0.001 ***
#   rasValue_bio4    1    4.356 14.5171  0.001 ***
#   rasValue_bio1    1    2.732  9.1025  0.001 ***
#   rasValue_bio14   1    1.622  5.4057  0.001 ***
#   rasValue_bio17   1    1.895  6.3144  0.001 ***
#   rasValue_bio12   1    2.087  6.9550  0.001 ***
#   Residual       977  293.185
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#check for collinearity; threshold should be 20 (above 20 removed); but above 10 should also be inspected
#Removing one by one (or a few at the times if sure that they would result collinear anyway)
vif.cca(fwd.sel)
# rasValue_bio10 rasValue_bio18  rasValue_bio5 rasValue_bio15  rasValue_bio2
# 2888.207232      58.222430     333.952722       4.543789       7.172540
# rasValue_bio4  rasValue_bio1 rasValue_bio14 rasValue_bio17 rasValue_bio12
# 404.662022    3606.816697     146.899863     321.006445      96.910726

#remove first variable
#a. removing bio14 and bio17
dbrda <- capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 +
           rasValue_bio18 + rasValue_bio5 + rasValue_bio15 + rasValue_bio2 +
           rasValue_bio4 + rasValue_bio1 +
           rasValue_bio12, data = data.frame(env.z), distance = "bray")
vif.cca(dbrda)
# rasValue_bio10 rasValue_bio18  rasValue_bio5 rasValue_bio15  rasValue_bio2
# 2638.262463      21.159725     324.577317       1.661119       6.421342
# rasValue_bio4  rasValue_bio1 rasValue_bio12
# 358.005578    3187.609342      20.821558

#b. removing bio5
dbrda <- capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 +
                    rasValue_bio18 + rasValue_bio15 + rasValue_bio2 +
                    rasValue_bio4 + rasValue_bio1 +
                    rasValue_bio12, data = data.frame(env.z), distance = "bray")
vif.cca(dbrda)
# rasValue_bio10 rasValue_bio18 rasValue_bio15  rasValue_bio2  rasValue_bio4
# 1599.786989      20.004945       1.629620       1.887979     326.539120
# rasValue_bio1 rasValue_bio12
# 2874.384754      20.114233

#c. removing bio1
dbrda <- capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 +
                    rasValue_bio18 + rasValue_bio15 + rasValue_bio2 +
                    rasValue_bio4 +
                    rasValue_bio12, data = data.frame(env.z), distance = "bray")
vif.cca(dbrda)
# rasValue_bio10 rasValue_bio18 rasValue_bio15  rasValue_bio2  rasValue_bio4
# 1.811938      19.020321       1.572339       1.844141       3.495063
# rasValue_bio12
# 20.058045

#d. removing bio12
dbrda <- capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 +
                    rasValue_bio18 + rasValue_bio15 + rasValue_bio2 +
                    rasValue_bio4, data = data.frame(env.z), distance = "bray")
vif.cca(dbrda)
# rasValue_bio10 rasValue_bio18 rasValue_bio15  rasValue_bio2  rasValue_bio4
# 1.805544       2.918990       1.462253       1.661183       3.285190

#e. re-add swe
dbrda <- capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 +
                    rasValue_bio18 + rasValue_bio15 + rasValue_bio2 +
                    rasValue_bio4 + rasValue_swe, data = data.frame(env.z), distance = "bray")
vif.cca(dbrda)
# rasValue_bio10 rasValue_bio18 rasValue_bio15  rasValue_bio2  rasValue_bio4
# 1.808832       5.368751       1.472217       1.682156       3.335474
# rasValue_swe
# 4.070677

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
# capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 +
#            rasValue_bio18 + rasValue_bio2 + rasValue_bio15 + rasValue_bio4,
#          data = data.frame(env.z), distance = "bray")

#I am not including snow variables then as no relevant! I tried to plot dbrda with also swe but does not look good anyway

#run final dbrda
dbrda_final <- capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 +
           rasValue_bio18 + rasValue_bio2 + rasValue_bio15 + rasValue_bio4,
           data = data.frame(env.z), distance = "bray")

summary(dbrda_final)
# Partitioning of squared Bray distance:
# Inertia Proportion
# Total          343.52     1.0000
# Constrained     36.78     0.1071
# Unconstrained  306.74     0.8929

#check R2
RsquareAdj(dbrda_final)
# $r.squared
# [1] 0.1249245
# 
# $adj.r.squared
# [1] 0.1204689

#get statistical significance of the model
anova.cca(fwd.sel, permutations = 1000) #same results as anova()
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 1000
# 
# Model: capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 + rasValue_bio18 + rasValue_bio2 + rasValue_bio15 + rasValue_bio4, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs     F   Pr(>F)
# Model      5    36.78 23.55 0.000999 ***
#   Residual 982   306.74
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#get statistical significance of the model for each variable
anova.cca(fwd.sel, step = 1000, by = "term")
# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = genera_counts_tab_trim ~ rasValue_bio10 + rasValue_bio18 + rasValue_bio2 + rasValue_bio15 + rasValue_bio4, data = data.frame(env.z), distance = "bray")
# Df SumOfSqs      F Pr(>F)
# rasValue_bio10   1   13.904 44.511  0.001 ***
#   rasValue_bio18   1    7.815 25.018  0.001 ***
#   rasValue_bio2    1    5.789 18.533  0.001 ***
#   rasValue_bio15   1    5.230 16.742  0.001 ***
#   rasValue_bio4    1    4.043 12.944  0.001 ***
#   Residual       982  306.741
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
bioregions <- bioclim$bioregion
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
CairoSVG(file="dbrda_bioclim.svg",
         width=15,
         height=15,
         pointsize=12)

length(rownames(sc_si))
#[1] 988

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
## add arrows for effects of the expanatory variables
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

############################
###variation partitioning###
############################

#collinear variables in explanatory tables must not be removed prior to partitioning
#my understanding is that environmental variables MUST not be trasformed or normalized -- but Pedro says yes to be consistent with dbRDA analyses
#coordinates: pcnm: https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/pcnm
#https://r.qcbs.ca/workshop10/book-en/variation-partitioning.html

#metadata
metadata <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)[samples,]

#geographical data
geography <- metadata[samples,c(2,1)]
geography_pcnm <- as.data.frame(scores(pcnm(dist(geography))))

#bioregions
#I cannot use it because it is a quanlitative variable

library(vegan)
#environmental VS geography
vp <- varpart(genera_counts_tab_trim, ~ rasValue_bio10 + rasValue_bio18 +
                rasValue_bio4 + rasValue_bio15 + rasValue_bio2, 
              geography_pcnm, data = data.frame(env.z))
vp$part
# No. of explanatory tables: 2
# Total variation (SS): 38920
# Variance: 39.433
# No. of observations: 988
# 
# Partition table:
#   Df R.squared Adj.R.squared Testable
# [a+c] = X1            5   0.09767       0.09308     TRUE
# [b+c] = X2           26   0.21978       0.19867     TRUE
# [a+b+c] = X1+X2      31   0.24971       0.22539     TRUE
# Individual fractions
# [a] = X1|X2           5                 0.02671     TRUE
# [b] = X2|X1          26                 0.13231     TRUE
# [c]                   0                 0.06637    FALSE
# [d] = Residuals                         0.77461    FALSE
# ---
#   Use function 'rda' to test significance of fractions of interest

#plot variation partitioning
library("Cairo")
CairoSVG(file="variation_partitioning.svg",
         width=8,
         height=6,
         pointsize=12)
plot(vp,
     Xnames = c("Bioclimatic dataset", "Geography"), # name the partitions
     bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
dev.off()

#test for significance

#geography without controlling for environmental variables
anova.cca(rda(genera_counts_tab_trim, geography_pcnm))
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(X = genera_counts_tab_trim, Y = geography_pcnm)
# Df Variance      F Pr(>F)
# Model     26   8.6666 10.412  0.001 ***
#   Residual 961  30.7659
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#environmental variables without controlling geography
anova.cca(rda(genera_counts_tab_trim ~ rasValue_bio10 + rasValue_bio18 +
                rasValue_bio4 + rasValue_bio15 + rasValue_bio2, 
              data = data.frame(env.z)))
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = genera_counts_tab_trim ~ rasValue_bio10 + rasValue_bio18 + rasValue_bio4 + rasValue_bio15 + rasValue_bio2, data = data.frame(env.z))
# Df Variance      F Pr(>F)
# Model      5    3.851 21.259  0.001 ***
#   Residual 982   35.581
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#geography alone
anova.cca(rda(genera_counts_tab_trim, geography_pcnm, data.frame(env.z)[,c("rasValue_bio10", "rasValue_bio18", 
                                                                          "rasValue_bio4", "rasValue_bio15", 
                                                                          "rasValue_bio2")]))
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(X = genera_counts_tab_trim, Y = geography_pcnm, Z = data.frame(env.z)[, c("rasValue_bio10", "rasValue_bio18", "rasValue_bio4", "rasValue_bio15", "rasValue_bio2", "rasValue_scd")])
# Df Variance      F Pr(>F)
# Model     26    6.150 7.7027  0.001 ***
#   Residual 955   29.327
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#environmental variables alone
anova.cca(rda(genera_counts_tab_trim, data.frame(env.z)[,c("rasValue_bio10", "rasValue_bio18", 
                                                           "rasValue_bio4", "rasValue_bio15", 
                                                           "rasValue_bio2")] , geography_pcnm))
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(X = genera_counts_tab_trim, Y = geography_pcnm, Z = data.frame(env.z)[, c("rasValue_bio10", "rasValue_bio18", "rasValue_bio4", "rasValue_bio15", "rasValue_bio2")])
# Df Variance      F Pr(>F)
# Model     26   5.9955 7.4512  0.001 ***
#   Residual 956  29.5856
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

save.image(file = "all_dataset_ps_srs_rda_vp_bacteria.RData")
