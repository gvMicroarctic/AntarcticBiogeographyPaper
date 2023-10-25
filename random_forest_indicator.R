######################################
###Format dataset for random forest###
######################################

#get dataset
setwd("/home/gilda/all/")

#import community information
load(file = "all_dataset_ps_srs_bacteria.RData")
community <- t(genera_proportions_tab[-1,])[samples,]
dim(community)
#[1]  988 1445

#trim dataset
unique(metadata$dataset)
# [1] "dt1"       "dt2"       "dt3"   "dt4"        "dt5"     
# [6] "dt6"        "dt71"      "dt72"         "dt8"      "dt9"      
# [11] "dt10"       "dt11"      "dt12"   "dt13"     "dt14"
# [16] "dt15"      "dt16"     "dt17"
metadata_trim <- metadata[metadata$dataset == "dt1" |
                            metadata$dataset == "dt3" |
                            metadata$dataset == "dt4" |
                            metadata$dataset == "dt5" |
                            metadata$dataset == "dt6" |
                            metadata$dataset == "dt71" |
                            metadata$dataset == "dt72" |
                            metadata$dataset == "dt8" |
                            metadata$dataset == "dt11" |
                            metadata$dataset == "dt12" |
                            metadata$dataset == "dt13" |
                            metadata$dataset == "dt14" |
                            metadata$dataset == "dt15",]
unique(metadata_trim$bioregion)
nrow(metadata_trim[metadata_trim$bioregion == 0,])
nrow(metadata_trim[metadata_trim$bioregion == 1,])
nrow(metadata_trim[metadata_trim$bioregion == 3,])
nrow(metadata_trim[metadata_trim$bioregion == 4,])
nrow(metadata_trim[metadata_trim$bioregion == 8,])
nrow(metadata_trim[metadata_trim$bioregion == 9,])
nrow(metadata_trim[metadata_trim$bioregion == 10,])

#metadata_trim <- metadata_trim[,c(3:12,14)]
metadata_trim <- metadata_trim[,c("rasValue_bio10", "rasValue_bio18", "rasValue_bio4", "rasValue_bio15", "rasValue_bio2", "rasValue_swe")]
dim(metadata_trim)
#[1] 503  6
samples_trim <- row.names(metadata_trim)
metadata_trim$rasValue_swe[is.na(metadata_trim$rasValue_swe)] <- 0

community_trim0 <- community[samples_trim,]
dim(community_trim0)
#[1]  503 1445

community_trim <- community_trim0[,colSums(community_trim0) > 0]
dim(community_trim)
#[1]  503 1303

#trim the community: I am keeping only genera present in more than 1/10 of samples
#dominant genera: genera present with a relative abundance higher or equal to 1 in at least one sample
community_dominant0 <- community_trim
community_dominant0[community_dominant0 > 0] = 1
genera_dominant_1 <- colnames(community_dominant0[,colSums(community_dominant0) > nrow(community_trim)/10])

community_dominant0 <- community_trim
community_dominant0[community_dominant0 < 1] = 0
community_dominant0[community_dominant0 >= 1] = 1
genera_dominant_2 <- colnames(community_dominant0[,colSums(community_dominant0) > 1])

genera_dominant <- intersect(genera_dominant_1, genera_dominant_2)
community_dominant <- community_trim[,genera_dominant]
dim(community_dominant)
#[1] 503 149

###########################
###Heatmap rel abundance###
###########################

community_abs_uncl <- data.frame(t(genera_and_unidentified_counts_tab)[row.names(community_dominant),])
dim(community_abs_uncl)
#[1]  457 1446

#metadata with sorted bioregion
metadata_bio <- metadata[row.names(community_dominant),]
metadata_rel <- metadata_bio[order(metadata_bio$bioregion),]
dim(metadata_rel)
#[1] 457  20

#I need to use community with also unclassified
df_sum <- data.frame(rep(0, ncol(community_abs_uncl)))
for(i in unique(metadata_rel$bioregion)) {
  print(i)
  df_sub <- community_abs_uncl[row.names(metadata_rel[metadata_rel$bioregion == i,]),]
  dim(df_sub)
  sum <- colSums(df_sub)
  df_sum <- cbind(df_sum, sum)
}
df_sum <- df_sum[,-1]
colnames(df_sum) <- unique(metadata_rel$bioregion)
row.names(df_sum) <- colnames(community_abs_uncl)
dim(df_sum)
#[1] 1446    8

#generate relative abundance table
df_sum_rel <- apply(df_sum, 2, function(x) x/sum(x)*100)
dim(df_sum_rel)
row.names(df_sum_rel) <- row.names(df_sum)
#[1] 1446    8

colSums(df_sum_rel)

################################
###Run random forest analysis###
################################

library(randomForest)
library(ppcor)
library(rfPermute)

row.names(metadata_trim) == row.names(community_dominant)

df_r <- data.frame()
df_mse <- data.frame(rep(0,ncol(metadata_trim)))
for (i in 1:ncol(community_dominant)) {
  print(i)
  data = cbind(response=community_dominant[,i], metadata_trim)
  classify <- randomForest(response~., data = data, ntree = 1000, do.trace=100, importance=TRUE)
  r <- classify$rsq[1000]*100
  if (r > 30) {
    df_r <- rbind(df_r, c(colnames(community_dominant)[i],r))
    df_mse <- cbind(df_mse, as.numeric(importance(classify)[,1]))
  }
}

colnames(df_r) <- c("genus", "R")
df_r
dim(df_r)
#[1] 51  2

df_mse <- df_mse[,-1]
row.names(df_mse) <- colnames(metadata_trim)
colnames(df_mse) <- df_r$genus
df_mse
dim(df_mse)
#[1] 6 51

#get statistics
#get the variable with highest mse for each genus
df_count <- data.frame(1)
for(i in 1:ncol(df_mse)) {
  df_count <-rbind(df_count, row.names(df_mse[df_mse[,i] == max(df_mse[,i]),]))
}
dim(df_count)
#[1] 52  1

#how many first hits per bioclimatic variable
count_mse <- c(length(df_count[df_count[,1] == "rasValue_swe",]),length(df_count[df_count[,1] == "rasValue_bio2",]),
  length(df_count[df_count[,1] == "rasValue_bio4",]), length(df_count[df_count[,1] == "rasValue_bio10",]),
  length(df_count[df_count[,1] == "rasValue_bio15",]), length(df_count[df_count[,1] == "rasValue_bio18",]))
count_mse

#averaged mse for each bioclimatic variable
average_mse <- c(mean(as.numeric(df_mse[6,])), mean(as.numeric(df_mse[5,])), mean(as.numeric(df_mse[3,])), mean(as.numeric(df_mse[1,])),
mean(as.numeric(df_mse[4,])), mean(as.numeric(df_mse[2,])))
average_mse

mse <- cbind(count_mse, average_mse, var = c("rasValue_swe","rasValue_bio2","rasValue_bio4","rasValue_bio10","rasValue_bio15","rasValue_bio18"))
mse <- as.data.frame(mse)
mse$var <- factor(mse$var, c("rasValue_swe","rasValue_bio2","rasValue_bio4","rasValue_bio10","rasValue_bio15","rasValue_bio18"))
mse$count_mse <- as.numeric(mse$count_mse)
mse$average_mse <- as.numeric(mse$average_mse)

#plot
library(ggplot2)
p1 <- ggplot(mse, aes(y=count_mse, x=var)) + 
  geom_bar(position='stack', stat='identity', fill="orange") +
  theme_bw()

p2 <- ggplot(mse, aes(y=average_mse, x=var)) + 
  geom_bar(position='stack', stat='identity', fill="blue") +
  theme_bw()

library(gridExtra)
grid.arrange(p1, p2, nrow = 2)
#400x800: barplot_random_forest.svg

###########################################################################
###Semi-partial correlation analyses on taxa selected from random forest###
###########################################################################

#set-up tables for results
n <- nrow(df_mse) + 1

table_e <- as.data.frame(matrix(0,nrow(df_mse),ncol(df_mse)))
colnames(table_e) <- colnames(df_mse)

table_p <- as.data.frame(matrix(0,nrow(df_mse),ncol(df_mse)))
colnames(table_p) <- colnames(df_mse)

#perform best-subset correlation with Spearman correlations
for(i in 1:length(colnames(df_mse))) {
  print(i)
  community_sub <- community_trim[,colnames(df_mse)[i]]
  combined <- cbind(metadata_trim, taxon=community_sub)
  table = spcor(combined, method = c("spearman"))
  table_e[,i] <- table$estimate[c(1:nrow(df_mse)),n]
  table_p[,i] <- table$p.value[c(1:nrow(df_mse)),n]
}
row.names(table_e) <- row.names(df_mse)
row.names(table_p) <- row.names(df_mse)

table_e
table_p

##heatmap to report random forest results
#keep only significant pvalues
table_p0 <- table_p
table_e0 <- table_e
table_p0[table_p0 >= 0.01] = 100
table_e0[table_p0 >= 0.01] = 100
table_e0 <- table_e0[,colSums(table_p0) != 600]
#set all statistics = 1 if value > 0
table_e0[table_e0 > 0 & table_e0 != 100] = 1
#set all statistics = -1 if value < 0
table_e0[table_e0 < 0] = -1
#format characters
colnames(table_e0) <- gsub("-", ".", colnames(table_e0))
colnames(table_e0) <- gsub(" ", ".", colnames(table_e0))
geo <- c("rasValue_swe","rasValue_bio2","rasValue_bio4","rasValue_bio10","rasValue_bio15","rasValue_bio18")
table_heat_rf <- t(table_e0[geo,])
#select and format genera selcted by random forest analysis
taxa0 <- row.names(table_heat_rf)
taxa0 <- gsub("-", ".", taxa0)
taxa_rf <- gsub(" ", ".", taxa0)

################################
###Bacterial indicator genera###
################################

library("ggplot2")
library("microeco")
library("phyloseq")
library("file2meco")

ps_srs1 = subset_samples(ps_srs, sample_names(ps_srs) %in% samples)
ps_srs1
ps_srs2 <- subset_samples(ps_srs1, dataset %in% c("dt1", "dt3", "dt4", "dt5", "dt6", "dt71", "dt72", "dt8", "dt11", "dt12",  "dt13", "dt14", "dt15"))
ps_srs2

#add bioregion information
metadata_trim <- metadata[row.names(sample_data(ps_srs2)),]
sample_data(ps_srs2)$bioregion <- metadata_trim$bioregion

#convert phyloseq object
ps_srs_meco <- phyloseq2meco(ps_srs2)
head(ps_srs_meco$sample_table)
ps_srs_meco$cal_abund()
row.names(ps_srs_meco$sample_table) == colnames(ps_srs_meco$otu_table)

ps_srs_meco$sample_table$bioregion
unique(ps_srs_meco$sample_table$bioregion)

#lefse
lefse <- trans_diff$new(dataset = ps_srs_meco, taxa_level = "Genus", method = "lefse", group = "bioregion", alpha = 0.01, lefse_subgroup = NULL, p_adjust_method = "holm")

#get genera selected from taxon indicator
genus_bi0 <- row.names(lefse$res_diff)
genus_bi <- gsub(".*_","",genus_bi0)
phylum_bi <- gsub(".*p__(.+)\\|c__.*", "\\1", genus_bi0)
taxa0 <- genus_bi
taxa0 <- gsub("-", ".", taxa0)
taxa <- gsub(" ", ".", taxa0)
taxa_int <- intersect(taxa, row.names(df_sum_rel))

#format table with all information
lr <- cbind(lefse$res_diff, genus=taxa, phylum=phylum_bi)
row.names(lr) <- lr$genus
lr

##heatmap to report indicator taxa results
setdiff(taxa_rf, row.names(lr))
#character(0)
#only genera present in both
int <- intersect(taxa_rf, row.names(lr))
int
#get information about group and phylum
table_heat_it0 <- as.data.frame(lr[int,c(4,10)])
table_heat_it <- table_heat_it0[taxa_rf,]

write.csv(lefse$res_diff, "lefse_results.csv")

###################################
###Heatmap rel abundance, subset###
###################################

table_heat_ra <- df_sum_rel[row.names(table_heat_it),]
dim(table_heat_ra)
#[1] 42  8

#######################
###Plot all heatmaps###
#######################

library(ComplexHeatmap)
library("cluster")

ha1 = rowAnnotation(table_heat_rf, annotation_name_side = "right")

#random forest
ht1 = Heatmap(as.matrix(table_heat_rf[,c("rasValue_bio2", "rasValue_bio4", "rasValue_bio10", "rasValue_bio15")]), name = "ht1", 
              col = c("100" = "white", "-1" = "darkorchid4", "1" = "orange"),
              column_names_gp = gpar(fontsize = 9), column_order = c("rasValue_bio2", "rasValue_bio4", "rasValue_bio10", "rasValue_bio15"))
ht1 #removed rasValue_swe and rasValue_bio18 because no significant correlations
dim(table_heat_rf)

#relative abundance
ht2 = Heatmap(as.matrix(table_heat_ra), name = "ht2", 
              col = circlize::colorRamp2(c(0, 1, 3), c("ivory1", "gray73", "black")),
              column_names_gp = gpar(fontsize = 9),
              column_order = c("0", "1", "3", "4", "8", "9", "10", "12"))
ht2

#indicator taxa
unique(table_heat_it$Group)
#[1] "0"  "12" "4"  "3"  "8"  "9"  "1"  "10"
ht3 = Heatmap(as.matrix(table_heat_it)[,1], name = "ht3", 
              col = c("10" = "#4682b3ff", "9" = "#26a8f6ff", "1" = "#8b0001ff", "0" = "black", "12" = "#808080ff", "3" = "#fe0000ff", "8" = "#006400ff", "4" = "#ffff00ff"),
              column_names_gp = gpar(fontsize = 9))
ht3

#phyla
unique(table_heat_it$phylum)
#[1] "Myxococcota"       "Actinobacteriota"  "Proteobacteria"    "Verrucomicrobiota" "Bacteroidota"      "Chloroflexi"       "Cyanobacteria"     "Nitrospirota"  "Acidobacteriota"
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(8, "Paired"))(length(unique(table_heat_it$phylum)))
ht4 = Heatmap(as.matrix(table_heat_it)[,2], name = "ht4",
              col = colors,
              column_names_gp = gpar(fontsize = 9))
ht4

#plot all
ht1 + ht2 + ht4 + ht3
#combined_heatmap: 1200 x 900

#####################################################
###Plot heatmap with the remaining dominant genera###
#####################################################

#dominant taxa
dominant_taxa0 <- colnames(community_dominant)
dominant_taxa0 <- gsub("-", ".", dominant_taxa0)
dominant_taxa <- gsub(" ", ".", dominant_taxa0)
length(dominant_taxa)
#[1] 149

#get dominant genera not represented in above analyses
genera_dominant_rest <- setdiff(dominant_taxa, row.names(table_heat_it))
length(genera_dominant_rest)
#[1] 107

#subset table for right genera
table_heat_ra1 <- df_sum_rel[genera_dominant_rest,]
dim(table_heat_ra1)
#[1] 107   8

#relative abundance
ht = Heatmap(as.matrix(table_heat_ra1), name = "ht", 
             col = circlize::colorRamp2(c(0, 1, 3), c("ivory1", "gray73", "black")),
             column_names_gp = gpar(fontsize = 9),
             column_order = c("0", "1", "3", "4", "8", "9", "10", "12"))
ht
#600x2500: heatmap_genus_si.svg