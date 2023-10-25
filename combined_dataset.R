########################
###Merge all datasets###
########################

setwd("/home/gilda/all")

library(dada2)
library(phyloseq)
library(ggplot2)
library(Biostrings)
library(decontam)
library(dplyr)
library(tidyr)
library(vegan)
library(zetadiv)
library(geosphere)

###############
###Dataset 1###
###############

load("dataset1_16S.RData")

metadata <- read.table("metadata_dataset1.txt", head=T, row.names = 1) #new
rownames(metadata) <- metadata$sample
dim(seqtab.nochim)
#[1]   175 12147
dim(metadata)
#[1] 175   9
dim(taxa)
#[1] 12147     7

ps_dt1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt1))
names(dna) <- taxa_names(ps_dt1)
ps_dt1 <- merge_phyloseq(ps_dt1, dna)
taxa_names(ps_dt1) <- paste0("ASV_N", seq(ntaxa(ps_dt1)))
ps_dt1

#remove blank samples
#get mock samples
blanks <- row.names(metadata[metadata$location == "blank", ])
blanks
#[1] "PCR_B1" "PCR_B3" "PCR_B4" "SGB1"   "SGB2"   "SGB3"   "SGB4"   "SGB5"
sample_data(ps_dt1)$is.neg <- sample_data(ps_dt1)$location == "blank"
#use and try out different decontam methods. Decided for "combined" at the end
blank_ASVs <- isContaminant(ps_dt1, method="prevalence", neg="is.neg", threshold=0.5)
blank_ASVs_id <- row.names(blank_ASVs[blank_ASVs[,"contaminant"] == TRUE,])
table(blank_ASVs$contaminant) #identified 17 as contaminants (vs 12130 no contaminants)
#FALSE  TRUE 
#12130    17 
which(blank_ASVs$contaminant)
# [1]  144  219  232  275  371  625  711  767  833  897  946  979 1080 1220 1843 4613 8378

#format asv table
asv <- t(otu_table(ps_dt1))
dim(asv)
#[1] 12147   175
asv_filt <- asv[!rownames(asv) %in% blank_ASVs_id,] #filter contaminant ASVs
dim(asv_filt)
#[1] 12130   175
asv_filt <- asv_filt[rowSums(asv_filt)> 0,]
dim(asv_filt)
#[1] 12130   175
asv_filt <- asv_filt[,!colnames(asv_filt) %in% blanks] #filter blank samples
dim(asv_filt)
#[1] 12130   167
asv_filt <- asv_filt[,!colnames(asv_filt) %in% "MSP1_3"] #removing this sample because NO geographical coordinates
dim(asv_filt)
#[1] 12130   166
asv_filt <- asv_filt[,!colnames(asv_filt) %in% "NP2_1"] #removing this sample because wrong geographical coordinates
dim(asv_filt)
#[1] 12130   165
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 12130   165

#format metadata table
retrieve <- as.vector(colnames(asv_filt))
metadata <- sample_data(ps_dt1)[retrieve,c(2:9)]
dim(metadata)
#[1] 165   8

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #165   8
dim(metadata_trim) #154   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #12130   154

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa <- tax_table(ps_dt1)[retrieve,]
dim(taxa)
#[1] 12130     7

#create new phyloseq object
ps_dt1_formatted <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                    sample_data(metadata_trim), 
                    tax_table(taxa))

#save phyloseq object
save(ps_dt1_formatted, file = "all_dataset.RData")
rm(list = ls())

###############
###Dataset 2###
###############

load("dataset2_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   120 15968
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 15968   120
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 15968   119

#format metadata table
metadata <- read.table("metadata_dataset2.txt", head=T, row.names = 1)
dim(metadata)
#[1] 120   8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 119   8

#format taxa table
dim(taxa)
#[1] 15968     7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #119   8
dim(metadata_trim) #116   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #15968   116

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 15968     7

ps_dt2 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                    sample_data(metadata_trim), 
                    tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt2))
names(dna) <- taxa_names(ps_dt2)
ps_dt2 <- merge_phyloseq(ps_dt2, dna)
taxa_names(ps_dt2) <- paste0("ASV_G", seq(ntaxa(ps_dt2)))
ps_dt2

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, file = "all_dataset.RData")
rm(list = ls())

###############
###Dataset 3###
###############

load("dataset3_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   41 15342
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 15342    41
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 15342    41

#format metadata table
metadata <- read.table("metadata_dataset3.txt", head=T, row.names = 1)
rownames(metadata) <- metadata$sample
dim(metadata)
#[1] 41   9
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,c(2:9)]
dim(metadata)
#[1] 41   8

#format taxa table
dim(taxa)
#[1] 15342     7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #41   8
dim(metadata_trim) #40   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #15342    40

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 15342     7

ps_dt3 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                    sample_data(metadata_trim), 
                    tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt3))
names(dna) <- taxa_names(ps_dt3)
ps_dt3 <- merge_phyloseq(ps_dt3, dna)
taxa_names(ps_dt3) <- paste0("ASV_L", seq(ntaxa(ps_dt3)))
ps_dt3

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, file = "all_dataset.RData")
rm(list = ls())

###############
###Dataset 4###
###############

load("dataset4_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   16 9117
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 9117   16
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 9117   16

#format metadata table
metadata <- read.table("metadata_dataset4.txt", head=T, row.names = 1)
dim(metadata)
#[1] 16   8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 16   8

#format taxa table
dim(taxa)
#[1] 9117    7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #16   8
dim(metadata_trim) #16   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #2027   16

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 9117    7

ps_dt4 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                    sample_data(metadata_trim), 
                    tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt4))
names(dna) <- taxa_names(ps_dt4)
ps_dt4 <- merge_phyloseq(ps_dt4, dna)
taxa_names(ps_dt4) <- paste0("ASV_M", seq(ntaxa(ps_dt4)))
ps_dt4

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, file = "all_dataset.RData")
rm(list = ls())

###############
###Dataset 5###
###############

load("dataset5_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]     6 12158
#remove underscore from sample name
metadata <- read.table("metadata_dataset5.txt", head=T, row.names = 1)
row.names(seqtab.nochim) <- row.names(metadata)
#format asv table
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 12158     6
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 12158     6

#format metadata table
metadata <- read.table("metadata_dt5.txt", head=T, row.names = 1)
dim(metadata)
#[1] 6 8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 6 8

#format taxa table
dim(taxa)
#[1] 12158     7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #6   8
dim(metadata_trim) #6   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #12158     6

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 12158     7

ps_dt5 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                    sample_data(metadata_trim), 
                    tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt5))
names(dna) <- taxa_names(ps_dt5)
ps_dt5 <- merge_phyloseq(ps_dt5, dna)
taxa_names(ps_dt5) <- paste0("ASV_ME", seq(ntaxa(ps_dt5)))
ps_dt5

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, file = "all_dataset.RData")
rm(list = ls())

#################
###Dataset 7.2###
#################

load("dataset7_2_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   16 4809
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 4809   16
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 4809   16

#format metadata table
metadata <- read.table("metadata_dataset7_2.txt", head=T, row.names = 1)
dim(metadata)
#[1] 16 8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 16 8

#format taxa table
dim(taxa)
#[1] 4809     7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #16   8
dim(metadata_trim) #16   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #4809   16

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 4809    7

ps_dt72 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                     sample_data(metadata_trim), 
                     tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt72))
names(dna) <- taxa_names(ps_dt72)
ps_dt72 <- merge_phyloseq(ps_dt72, dna)
taxa_names(ps_dt72) <- paste0("ASV_C", seq(ntaxa(ps_dt72)))
ps_dt72

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, file = "all_dataset.RData")
rm(list = ls())

###############
###Dataset 6###
###############

load("dataset6_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   35 15557
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 15557    35
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 15557    35

#format metadata table
metadata <- read.table("metadata_dataset6.txt", head=T, row.names = 1)
dim(metadata)
#[1] 35 8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 35 8

#format taxa table
dim(taxa)
#[1] 15557     7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #35   8
dim(metadata_trim) #35   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #15557    35

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 15557     7

ps_dt6<- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                   sample_data(metadata_trim), 
                   tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt6))
names(dna) <- taxa_names(ps_dt6)
ps_dt6 <- merge_phyloseq(ps_dt6, dna)
taxa_names(ps_dt6) <- paste0("ASV_P", seq(ntaxa(ps_dt6)))
ps_dt6
#plot_richness(ps_dt6, color="location", measures=c("Shannon", "Simpson"))

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, file = "all_dataset.RData")
rm(list = ls())

#################
###Dataset 7.1###
#################

load("dataset7_1_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   53 30986
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 30986    53
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 30986    53

#format metadata table
metadata <- read.table("metadata_dataset7_1.txt", head=T, row.names = 1)
dim(metadata)
#[1] 53 8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 53 8

#format taxa table
dim(taxa)
#[1] 30986     7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #53   8
dim(metadata_trim) #53   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #30986    53

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 30986     7

ps_dt71 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                    sample_data(metadata_trim), 
                    tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt71))
names(dna) <- taxa_names(ps_dt71)
ps_dt71 <- merge_phyloseq(ps_dt71, dna)
taxa_names(ps_dt71) <- paste0("ASV_AP", seq(ntaxa(ps_dt71)))
ps_dt71

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, file = "all_dataset.RData")
rm(list = ls())

###############
###Dataset 8###
###############

load("dataset8_16S.RData")

metadata <- read.table("metadata_dataset8.txt", head=T, row.names = 1) #new
dim(seqtab.nochim)
#[1]  114 7700
dim(metadata)
#[1] 57   8
dim(taxa)
#[1] 7700     7
seqtab.nochim_trim <- seqtab.nochim[row.names(metadata),]
dim(seqtab.nochim_trim)
#57 7700

ps_dt8 <- phyloseq(otu_table(seqtab.nochim_trim, taxa_are_rows=FALSE), 
                     sample_data(metadata), 
                     tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt8))
names(dna) <- taxa_names(ps_dt8)
ps_dt8 <- merge_phyloseq(ps_dt8, dna)
taxa_names(ps_dt8) <- paste0("ASV_B", seq(ntaxa(ps_dt8)))
ps_dt8

#remove blank samples

#get mock samples
blanks <- row.names(metadata[metadata$location == "blank", ])
blanks
#[1] "N1"      "N2"      "N3"      "N4"      "N5"      "REAGENT"
sample_data(ps_dt8)$is.neg <- sample_data(ps_dt8)$location == "blank"
#use and try out different decontam methods. Decided for "combined" at the end
blank_ASVs <- isContaminant(ps_dt8, method="prevalence", neg="is.neg", threshold=0.5)
blank_ASVs_id <- row.names(blank_ASVs[blank_ASVs[,"contaminant"] == TRUE,])
table(blank_ASVs$contaminant)
#FALSE  TRUE 
#7697     3
which(blank_ASVs$contaminant)
#[1] 1522 2567 3084

#format asv table
asv <- t(otu_table(ps_dt8))
dim(asv)
#[1] 7700   57
asv_filt <- asv[!rownames(asv) %in% blank_ASVs_id,] #filter contaminant ASVs
dim(asv_filt)
#[1] 7697   57
asv_filt <- asv_filt[rowSums(asv_filt)> 0,]
dim(asv_filt)
#[1] 7125   57
asv_filt <- asv_filt[,!colnames(asv_filt) %in% blanks] #filter blank samples
dim(asv_filt)
#[1] 7125   51
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 7125   51

#format metadata table
retrieve <- as.vector(colnames(asv_filt))
metadata <- sample_data(ps_dt8)[retrieve,c(1:8)]
dim(metadata)
#[1] 51  8

#format taxa table
retrieve <- as.vector(row.names(asv_filt))
taxa <- tax_table(ps_dt8)[retrieve,]
dim(taxa)
#[1] 7125    7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #51   8
dim(metadata_trim) #45   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #7125    45

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 7125     7

#create new phyloseq object
ps_dt8_formatted <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                               sample_data(metadata_trim), 
                               tax_table(taxa_trim))

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, ps_dt8_formatted, file = "all_dataset.RData")
rm(list = ls())

################
###Dataset 10###
################

load("dataset10_16S.RData")

metadata <- read.table("metadata_dataset10.txt", head=T, row.names = 1) #new
dim(seqtab.nochim)
#[1]  153 55961
dim(metadata)
#[1] 153   8
dim(taxa)
#[1] 55961     7

ps_dt10 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                    sample_data(metadata), 
                    tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt10))
names(dna) <- taxa_names(ps_dt10)
ps_dt10 <- merge_phyloseq(ps_dt10, dna)
taxa_names(ps_dt10) <- paste0("ASV_MARK", seq(ntaxa(ps_dt10)))
ps_dt10

#remove blank samples

#get mock samples
blanks <- row.names(metadata[metadata$location == "blank", ])
blanks
#[1] "water"
sample_data(ps_dt10)$is.neg <- sample_data(ps_dt10)$location == "blank"
#use and try out different decontam methods. Decided for "combined" at the end
blank_ASVs <- isContaminant(ps_dt10, method="prevalence", neg="is.neg", threshold=0.5)
blank_ASVs_id <- row.names(blank_ASVs[blank_ASVs[,"contaminant"] == TRUE,])
table(blank_ASVs$contaminant)
#FALSE  TRUE 
#55960     1
which(blank_ASVs$contaminant)
#[1] 21108

#format asv table
asv <- t(otu_table(ps_dt10))
dim(asv)
#[1] 55961   153
asv_filt <- asv[!rownames(asv) %in% blank_ASVs_id,] #filter contaminant ASVs
dim(asv_filt)
#[1] 55960   153
asv_filt <- asv_filt[rowSums(asv_filt)> 0,]
dim(asv_filt)
#[1] 55960   153
asv_filt <- asv_filt[,!colnames(asv_filt) %in% blanks] #filter blank samples
dim(asv_filt)
#[1] 55960   152
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 55960   152

#format metadata table
retrieve <- as.vector(colnames(asv_filt))
metadata <- sample_data(ps_dt10)[retrieve,c(1:8)]
dim(metadata)
#[1] 152   8

#format taxa table
retrieve <- as.vector(row.names(asv_filt))
taxa <- tax_table(ps_dt10)[retrieve,]
dim(taxa)
#[1] 55960     7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #152   8
dim(metadata_trim) #152   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #55960   152

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 55960     7

#create new phyloseq object
ps_dt10_formatted <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                              sample_data(metadata_trim), 
                              tax_table(taxa_trim))

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, ps_dt8_formatted, ps_dt10_formatted, file = "all_dataset.RData")
rm(list = ls())

###############
###Dataset 9###
###############

load("dataset9_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   37 7355
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 7355   37
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 7355   34

#format metadata table
metadata <- read.table("metadata_dataset9.txt", head=T, row.names = 1)
dim(metadata)
#[1] 37  9
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 34 9

#rename asv and metadata table
row.names(metadata) == colnames(asv_filt)
colnames(asv_filt) <- metadata$sample
row.names(metadata) <- metadata$sample
metadata <- metadata[,2:9]
dim(metadata)
#[1] 34 8

#format taxa table
dim(taxa)
#[1] 7355    7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #34   8
dim(metadata_trim) #19   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #7355   19

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 7355    7

ps_dt9 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                    sample_data(metadata_trim), 
                    tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt9))
names(dna) <- taxa_names(ps_dt9)
ps_dt9 <- merge_phyloseq(ps_dt9, dna)
taxa_names(ps_dt9) <- paste0("ASV_PAUL", seq(ntaxa(ps_dt9)))
ps_dt9

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, ps_dt8_formatted, ps_dt10_formatted, ps_dt9, file = "all_dataset.RData")
rm(list = ls())

################
###Dataset 14###
################

load("dataset14_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   84 5828
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 5828   84
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 5828   80

#format metadata table
metadata <- read.table("metadata_dataset14.txt", head=T, row.names = 1)
dim(metadata)
#[1] 84  8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 80 8

#format taxa table
dim(taxa)
#[1] 5828    7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #80   8
dim(metadata_trim) #75   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #5828   75

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 5828    7

ps_dt14 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                     sample_data(metadata_trim), 
                     tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt14))
names(dna) <- taxa_names(ps_dt14)
ps_dt14 <- merge_phyloseq(ps_dt14, dna)
taxa_names(ps_dt14) <- paste0("ASV_SOL", seq(ntaxa(ps_dt14)))
ps_dt14

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, ps_dt8_formatted, ps_dt10_formatted, ps_dt9, ps_dt14, file = "all_dataset.RData")
rm(list = ls())

################
###Dataset 12###
################

load("dataset12_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   9 1001
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 1001    9
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 1001    9

#format metadata table
metadata <- read.table("metadata_dataset12.txt", head=T, row.names = 1)
dim(metadata)
#[1] 9  8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 9 8

#format taxa table
dim(taxa)
#[1] 1001    7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #9   8
dim(metadata_trim) #5   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #1001    5

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 1001    7

ps_dt12 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                        sample_data(metadata_trim), 
                        tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt12))
names(dna) <- taxa_names(ps_dt12)
ps_dt12 <- merge_phyloseq(ps_dt12, dna)
taxa_names(ps_dt12) <- paste0("ASV_BOR", seq(ntaxa(ps_dt12)))
ps_dt12

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, ps_dt8_formatted, ps_dt10_formatted, ps_dt9, ps_dt14, ps_dt12, file = "all_dataset.RData")
rm(list = ls())

################
###Dataset 11###
################

load("dataset11_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   36 79785
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 79785    36
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 79785    36

#format metadata table
metadata <- read.table("metadata_dataset11.txt", head=T, row.names = 1)
dim(metadata)
#[1] 36  9
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 36 9

#rename asv table
row.names(metadata) == colnames(asv_filt)
colnames(asv_filt) <- metadata$sample
row.names(metadata) <- metadata$sample
metadata <- metadata[,2:9]
dim(metadata)
#[1] 36 8

#format taxa table
dim(taxa)
#[1] 79785     7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #36   8
dim(metadata_trim) #36   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #79785    36

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 79785     7

ps_dt11 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                      sample_data(metadata_trim), 
                      tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt11))
names(dna) <- taxa_names(ps_dt11)
ps_dt11 <- merge_phyloseq(ps_dt11, dna)
taxa_names(ps_dt11) <- paste0("ASV_ALM", seq(ntaxa(ps_dt11)))
ps_dt11

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, ps_dt8_formatted, ps_dt10_formatted, ps_dt9, ps_dt14, ps_dt12, ps_dt11, file = "all_dataset.RData")
rm(list = ls())

################
###Dataset 13###
################

load("dataset13_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   60 5370
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 5370   60
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 5370   60

#format metadata table
metadata <- read.table("metadata_dataset13.txt", head=T, row.names = 1)
dim(metadata)
#[1] 60  8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 60 8

#format taxa table
dim(taxa)
#[1] 5370    7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #60   8
dim(metadata_trim) #52   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #5370    52

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 5370     7

ps_dt13 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                          sample_data(metadata_trim), 
                          tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt13))
names(dna) <- taxa_names(ps_dt13)
ps_dt13 <- merge_phyloseq(ps_dt13, dna)
taxa_names(ps_dt13) <- paste0("ASV_SEV", seq(ntaxa(ps_dt13)))
ps_dt13

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, ps_dt8_formatted, ps_dt10_formatted, ps_dt9, ps_dt14, ps_dt12, ps_dt11, ps_dt13, file = "all_dataset.RData")
rm(list = ls())

################
###Dataset 15###
################

load("dataset15_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   2 8027
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 8027    2
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 8027    2

#format metadata table
metadata <- read.table("metadata_dataset15.txt", head=T, row.names = 1)
dim(metadata)
#[1] 2  8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 2 8

#format taxa table
dim(taxa)
#[1] 8027    7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #2   8
dim(metadata_trim) #2   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #8027    2

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 8027     7

ps_dt15 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                          sample_data(metadata_trim), 
                          tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt15))
names(dna) <- taxa_names(ps_dt15)
ps_dt15 <- merge_phyloseq(ps_dt15, dna)
taxa_names(ps_dt15) <- paste0("ASV_ZH", seq(ntaxa(ps_dt15)))
ps_dt15

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, ps_dt8_formatted, ps_dt10_formatted, ps_dt9, ps_dt14, ps_dt12, ps_dt11, ps_dt13, ps_dt15, file = "all_dataset.RData")
rm(list = ls())

################
###Dataset 17###
################

load("dataset17_16S.RData")

#format asv table
dim(seqtab.nochim)
#[1]   53 11037
asv <- t(seqtab.nochim)
asv_filt <- asv[rowSums(asv)> 0,]
dim(asv_filt)
#[1] 11037    53
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 11037    53

#format metadata table
metadata <- read.table("metadata_dataset17.txt", head=T, row.names = 1)
dim(metadata)
#[1] 53  8
retrieve <- as.vector(colnames(asv_filt))
metadata <- metadata[retrieve,]
dim(metadata)
#[1] 53 8

#format taxa table
dim(taxa)
#[1] 11037     7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #53   8
dim(metadata_trim) #53   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #11037    53

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 11037     7

ps_dt17 <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                      sample_data(metadata_trim), 
                      tax_table(taxa_trim))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt17))
names(dna) <- taxa_names(ps_dt17)
ps_dt17 <- merge_phyloseq(ps_dt17, dna)
taxa_names(ps_dt17) <- paste0("ASV_TYT", seq(ntaxa(ps_dt17)))
ps_dt17

#save phyloseq object
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, ps_dt8_formatted, ps_dt10_formatted, ps_dt9, ps_dt14, ps_dt12, ps_dt11, ps_dt13, ps_dt15, ps_dt17, file = "all_dataset.RData")
rm(list = ls())

#################
###Dataset  16###
#################

load("dataset16_16S.RData")

metadata <- read.table("metadata_dataset16.txt", head=T, row.names = 1) #new
dim(seqtab.nochim)
#[1]  154 29392
dim(metadata)
#[1] 154   8
dim(taxa)
#[1] 29392     7

ps_dt16 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                  sample_data(metadata), 
                  tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps_dt16))
names(dna) <- taxa_names(ps_dt16)
ps_dt16 <- merge_phyloseq(ps_dt16, dna)
taxa_names(ps_dt16) <- paste0("ASV_VW", seq(ntaxa(ps_dt16)))
ps_dt16

#remove blank samples

#get mock samples
blanks <- row.names(metadata[metadata$location == "blank", ])
blanks
#[1] "Badt72ck_S101_L001" "SoilDNA_S102_L001" "NTC_S21_L001"
sample_data(ps_dt16)$is.neg <- sample_data(ps_dt16)$location == "blank"
#use and try out different decontam methods. Decided for "combined" at the end
blank_ASVs <- isContaminant(ps_dt16, method="prevalence", neg="is.neg", threshold=0.5)
blank_ASVs_id <- row.names(blank_ASVs[blank_ASVs[,"contaminant"] == TRUE,])
table(blank_ASVs$contaminant)
#FALSE  TRUE 
#29384     8
which(blank_ASVs$contaminant)
#[1] 1585  1787  2246  2655  3475  3880  7453 12350

#format asv table
asv <- t(otu_table(ps_dt16))
dim(asv)
#[1] 29392   154
asv_filt <- asv[!rownames(asv) %in% blank_ASVs_id,] #filter contaminant ASVs
dim(asv_filt)
#[1] 29384   154
asv_filt <- asv_filt[rowSums(asv_filt)> 0,]
dim(asv_filt)
#[1] 29384   154
asv_filt <- asv_filt[,!colnames(asv_filt) %in% blanks] #filter blank samples
dim(asv_filt)
#[1] 29384   151
asv_filt <- asv_filt[,!colnames(asv_filt) %in% "TRT32_S63_L001"] #removing this sample because NO geographical coordinates
dim(asv_filt)
#[1] 29384   150
asv_filt <- asv_filt[,colSums(asv_filt)> 0]
dim(asv_filt)
#[1] 29384   150

#format metadata table
retrieve <- as.vector(colnames(asv_filt))
metadata <- sample_data(ps_dt16)[retrieve,c(1:8)]
dim(metadata)
#[1] 150   8

#format taxa table
retrieve <- as.vector(row.names(asv_filt))
taxa <- tax_table(ps_dt16)[retrieve,]
dim(taxa)
#[1] 29384     7

#trim samples with less than 5000 reads
metadata_trim <- metadata[metadata$reads >= 5000,]
dim(metadata) #150   8
dim(metadata_trim) #146   8

#trim samples from asv table
asv_trim <- asv_filt[,row.names(metadata_trim)]
asv_trim <- asv_trim[,colSums(asv_trim)> 0]
dim(asv_trim) #29384   146

#trim taxa from taxa table
retrieve <- as.vector(row.names(asv_trim))
taxa_trim <- taxa[retrieve,]
dim(taxa_trim)
#[1] 29384     7

#create new phyloseq object
ps_dt16_formatted <- phyloseq(otu_table(asv_trim, taxa_are_rows=TRUE), 
                            sample_data(metadata_trim), 
                            tax_table(taxa_trim))

#save phyloseq object from dt10 dataset
load("all_dataset.RData")
save(ps_dt1_formatted, ps_dt2, ps_dt3, ps_dt4, ps_dt5, ps_dt72, ps_dt6, ps_dt71, ps_dt8_formatted, ps_dt10_formatted, ps_dt9, ps_dt14, ps_dt12, ps_dt11, ps_dt13, ps_dt15, ps_dt17, ps_dt16_formatted, file = "all_dataset.RData")
rm(list = ls())

####################
###merge datasets###
####################

library(phyloseq)

setwd("/home/gilda/all")

#load datasets
load("all_dataset.RData")

#merge datasets
ps_asv <- merge_phyloseq(otu_table(ps_dt2), otu_table(ps_dt1_formatted), otu_table(ps_dt3), otu_table(ps_dt4), otu_table(ps_dt5), 
                                                       otu_table(ps_dt72), otu_table(ps_dt6), otu_table(ps_dt71), otu_table(ps_dt8_formatted), 
                                                       otu_table(ps_dt10_formatted), otu_table(ps_dt9), otu_table(ps_dt14), otu_table(ps_dt12), 
                                                       otu_table(ps_dt11), otu_table(ps_dt13), otu_table(ps_dt15), otu_table(ps_dt17), otu_table(ps_dt16_formatted))
ps_metadata <- merge_phyloseq(sample_data(ps_dt2), sample_data(ps_dt1_formatted), sample_data(ps_dt3), sample_data(ps_dt4), sample_data(ps_dt5), 
                              sample_data(ps_dt72), sample_data(ps_dt6), sample_data(ps_dt71), sample_data(ps_dt8_formatted), 
                              sample_data(ps_dt10_formatted), sample_data(ps_dt9), sample_data(ps_dt14), sample_data(ps_dt12), 
                              sample_data(ps_dt11), sample_data(ps_dt13), sample_data(ps_dt15),sample_data(ps_dt17), sample_data(ps_dt16_formatted))
ps_taxa <- merge_phyloseq(tax_table(ps_dt2), tax_table(ps_dt1_formatted), tax_table(ps_dt3), tax_table(ps_dt4), tax_table(ps_dt5), 
                          tax_table(ps_dt72), tax_table(ps_dt6), tax_table(ps_dt71), tax_table(ps_dt8_formatted), 
                          tax_table(ps_dt10_formatted), tax_table(ps_dt9), tax_table(ps_dt14), tax_table(ps_dt12), 
                          tax_table(ps_dt11), tax_table(ps_dt13), tax_table(ps_dt15), tax_table(ps_dt17), tax_table(ps_dt16_formatted))
#export taxonomy file
dim(data.frame(ps_taxa))
#[1] 328385      7
write.csv(data.frame(ps_taxa), "taxonomy.csv")

ps <- phyloseq(otu_table(ps_asv, taxa_are_rows=TRUE), 
               sample_data(ps_metadata), 
               tax_table(ps_taxa))

#check number of samples
dim(ps_metadata[ps_metadata$dataset == 'dt1',]) #[1] 154   8
dim(ps_metadata[ps_metadata$dataset == 'dt2',]) #[1] 116   8
dim(ps_metadata[ps_metadata$dataset == 'dt3',]) #[1] 40  8
dim(ps_metadata[ps_metadata$dataset == 'dt4',]) #[1] 16  8
dim(ps_metadata[ps_metadata$dataset == 'dt11',]) #36  8
dim(ps_metadata[ps_metadata$dataset == 'dt12',]) #[1] 5 8
dim(ps_metadata[ps_metadata$dataset == 'dt8',]) #[1] 45  8
dim(ps_metadata[ps_metadata$dataset == 'dt71',]) #[1] 53  8
dim(ps_metadata[ps_metadata$dataset == 'dt72',]) #[1] 16  8
dim(ps_metadata[ps_metadata$dataset == 'dt10',]) #[1] 152   8
dim(ps_metadata[ps_metadata$dataset == 'dt5',]) #[1] 6   8
dim(ps_metadata[ps_metadata$dataset == 'dt9',]) #[1] 19   8
dim(ps_metadata[ps_metadata$dataset == 'dt6',]) #[1] 35   8
dim(ps_metadata[ps_metadata$dataset == 'dt13',]) #[1] 52   8
dim(ps_metadata[ps_metadata$dataset == 'dt15',]) #[1] 2   8
dim(ps_metadata[ps_metadata$dataset == 'dt14',]) #[1] 75   8
dim(ps_metadata[ps_metadata$dataset == 'dt17',]) #[1] 53   8
#sum
dim(ps_metadata[ps_metadata$dataset == 'vestfold',]) #[1] 86   8
dim(ps_metadata[ps_metadata$dataset == 'windmill',]) #[1] 60   8

#save(ps, file = "all_dataset_ps.RData")
#rm(list = ls())

#load("all_dataset_ps.RData")

#remove reads associated to unwanted taxonomy
ps_b0 <- subset_taxa(ps, Kingdom == "Bacteria")
ps_b1 <- subset_taxa(ps_b0, Family != "Mitochondria" | is.na(Family))
ps_b <- subset_taxa(ps_b1, Order != "Chloroplast" | is.na(Order))
#N.B. https://github.com/joey711/phyloseq/issues/999

#check how many reads before and after taxonomy trimming
read_num <- cbind(all=colSums(otu_table(ps)), excluded=colSums(otu_table(ps_b)))
dim(read_num[read_num[,2] >= 5000,])
#[1] 990   2

ps_bacteria <- prune_samples(sample_sums(ps_b)>=5000, ps_b)

#check number of samples
ps_metadata <- sample_data(ps_bacteria)
dim(ps_metadata[ps_metadata$dataset == 'dt1',]) #[1] 154   8
dim(ps_metadata[ps_metadata$dataset == 'dt2',]) #[1] 116   8
dim(ps_metadata[ps_metadata$dataset == 'dt3',]) #[1] 40  8
dim(ps_metadata[ps_metadata$dataset == 'dt4',]) #[1] 16  8
dim(ps_metadata[ps_metadata$dataset == 'dt11',]) #36  8
dim(ps_metadata[ps_metadata$dataset == 'dt12',]) #[1] 5 8
dim(ps_metadata[ps_metadata$dataset == 'dt8',]) #[1] 45  8
dim(ps_metadata[ps_metadata$dataset == 'dt71',]) #[1] 53  8
dim(ps_metadata[ps_metadata$dataset == 'dt72',]) #[1] 16  8
dim(ps_metadata[ps_metadata$dataset == 'dt10',]) #[1] 152   8
dim(ps_metadata[ps_metadata$dataset == 'dt5',]) #[1] 6   8
dim(ps_metadata[ps_metadata$dataset == 'dt9',]) #[1] 18   8
dim(ps_metadata[ps_metadata$dataset == 'dt6',]) #[1] 35   8
dim(ps_metadata[ps_metadata$dataset == 'dt13',]) #[1] 46   8
dim(ps_metadata[ps_metadata$dataset == 'dt15',]) #[1] 2   8
dim(ps_metadata[ps_metadata$dataset == 'dt14',]) #[1] 51   8
dim(ps_metadata[ps_metadata$dataset == 'dt17',]) #[1] 53   8
#sum
dim(ps_metadata[ps_metadata$dataset == 'vestfold',]) #[1] 86   8
dim(ps_metadata[ps_metadata$dataset == 'windmill',]) #[1] 60   8

#get data for map
coordinates <- data.frame(sample_data(ps_bacteria))[,c(3,4)]
write.table(coordinates, "all_coordinates.txt")

#asv table
asv <- otu_table(ps_bacteria)
dim(asv) 
#[1] 312839    990
min <- min(colSums(asv))
max <- max(colSums(asv))
min
#[1] 5247
max
#[1] 988658
Cmin <- min(colSums(asv))
Cmin
#[1] 5247

library(SRS)
asv_srs <- SRS(data.frame(asv, check.names = FALSE), Cmin)
row.names(asv_srs) <- row.names(asv)
dim(asv_srs)
#[1] 312839 990
asv_srs_trim <- asv_srs[rowSums(asv_srs)> 0,]
dim(asv_srs_trim)
#[1] 270985    990

#metadata
metadata_trim <- sample_data(ps_bacteria)
dim(metadata_trim)
#[1] 990   8

#taxonomy table
retrieve <- as.vector(row.names(asv_srs_trim))
taxa_trim <- tax_table(ps_bacteria)[retrieve,]
dim(taxa_trim)
#[1] 270985      7

#check
colnames(asv_srs_trim) == row.names(metadata_trim)
row.names(asv_srs_trim) == row.names(taxa_trim)

#create phyloseq object
ps_srs_0 <- phyloseq(otu_table(asv_srs_trim, taxa_are_rows = TRUE), 
                     sample_data(metadata_trim), 
                     taxa_trim)

library("ape")
random_tree = rtree(ntaxa(ps_srs_0), rooted=TRUE, tip.label=taxa_names(ps_srs_0))

ps_srs <- phyloseq(otu_table(asv_srs_trim, taxa_are_rows = TRUE), 
                   sample_data(metadata_trim), 
                   taxa_trim, random_tree)

#save(ps_srs, file = "all_dataset_ps_srs_bacteria.RData")
#rm(list = ls())

##############################
##Create genus-level dataset##
##############################

setwd("/home/gilda/all")

load("all_dataset_ps_srs_bacteria.RData")

library(phyloseq)

#making a vector of genus names to set as row names
genera_counts_tab <- otu_table(tax_glom(ps_srs, taxrank="Genus"))
dim(genera_counts_tab)
#[1] 1445  990
genera_tax_vec <- as.vector(tax_table(tax_glom(ps_srs, taxrank="Genus"))[,6])
rownames(genera_counts_tab) <- as.vector(genera_tax_vec)
unclassified_tax_counts <- colSums(otu_table(ps_srs)) - colSums(genera_counts_tab)
# and we'll add this row to our genus count table:
genera_and_unidentified_counts_tab <- rbind("Unclassified"=unclassified_tax_counts, genera_counts_tab)
genera_and_unidentified_counts_tab0 <- genera_and_unidentified_counts_tab[rowSums(genera_and_unidentified_counts_tab) > 0,]
genera_and_unidentified_counts_tab <- genera_and_unidentified_counts_tab0[,colSums(genera_and_unidentified_counts_tab0) > 0]

#format and save genus table
paths <- tax_table(tax_glom(ps_srs, taxrank="Genus"))
dim(paths)
#[1] 1446    5
#check same order of genera
paths[,6] == genera_tax_vec
paths[,6] == row.names(genera_and_unidentified_counts_tab[-1,])
paths_new <- rbind(c("Unclassified", "Unclassified", "Unclassified", "Unclassified", "Unclassified"), paths[,c(1:5)])
dim(paths_new)
#[1] 1446    5
#table to save
genera_ntw <- cbind(paths_new,genera_and_unidentified_counts_tab)
row.names(genera_ntw) <- row.names(genera_and_unidentified_counts_tab)
row.names(genera_ntw)
colnames(genera_ntw)
#write
write.csv(genera_ntw, "antarctic_soil_genus_taxonomy.csv")

#generate relative abundance table
genera_proportions_tab <- apply(genera_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
dim(genera_proportions_tab)
#[1] 1446  990
#get Hellinger transformed data
genera_proportions_tab_hellinger <- sqrt(genera_proportions_tab)

#nmds
library(vegan)
nmds <- metaMDS(t(genera_proportions_tab_hellinger), distance = "bray") #Best solution was not repeated -- monoMDS stopping criteria:
#plot(nmds)

#calculate number of dominant genera
#dominant genera: genera present with a relative abundance higher or equal to 1 in at least one sample
genera_proportions_tab0 <- genera_proportions_tab
genera_proportions_tab0[genera_proportions_tab0 < 1] = 0
dominant_genera_proportions_tab <- genera_proportions_tab0[rowSums(genera_proportions_tab0) > 0,]
dim(dominant_genera_proportions_tab)
#[1] 492 990

#number of genera
nrow(genera_proportions_tab) - 1
#[1] 1445

#number of dominant genera
nrow(dominant_genera_proportions_tab) - 1
#[1] 491

#get sample names to analyses: all except from samples with unknown component = 100%
#remove the two samples (ERX1598809, ERX1598810) that had unknown component = 100%
samples <- colnames(genera_proportions_tab[,genera_proportions_tab[1,]<100])
length(samples)

#write relative abundance table
genera_proportions_tab_ntw <- cbind(genera_ntw[, c(1:5)], genera_proportions_tab[,samples])
write.csv(genera_proportions_tab_ntw, "antarctic_soil_genus_taxonomy_rel_abundance.csv")

#save.image(file = "all_dataset_ps_srs_bacteria.RData")
#rm(list = ls())

###############################
##Create family-level dataset##
###############################

#making a vector of family names to set as row names
family_counts_tab <- otu_table(tax_glom(ps_srs, taxrank="Family"))
dim(family_counts_tab)
#[1] 496  990
family_tax_vec <- as.vector(tax_table(tax_glom(ps_srs, taxrank="Family"))[,5])
rownames(family_counts_tab) <- as.vector(family_tax_vec)
unclassified_tax_counts <- colSums(otu_table(ps_srs)) - colSums(family_counts_tab)
# and we'll add this row to our family count table:
family_and_unidentified_counts_tab <- rbind("Unclassified"=unclassified_tax_counts, family_counts_tab)
family_and_unidentified_counts_tab0 <- family_and_unidentified_counts_tab[rowSums(family_and_unidentified_counts_tab) > 0,]
family_and_unidentified_counts_tab <- family_and_unidentified_counts_tab0[,colSums(family_and_unidentified_counts_tab0) > 0]

#generate relative abundance table
family_proportions_tab <- apply(family_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
dim(family_proportions_tab)
#[1] 497  990

#write relative abundance table
family_proportions_tab_trim0 <- family_proportions_tab[,samples]
family_proportions_tab_trim <- family_proportions_tab_trim0[rowSums(family_proportions_tab_trim0) > 0,]
dim(family_proportions_tab_trim)
#[1] 497  988

#number of genera
nrow(family_proportions_tab_trim) - 5
#[1] 495

write.csv(family_proportions_tab_trim, "antarctic_soil_family_taxonomy_rel_abundance.csv")

##############################
##Create order-level dataset##
##############################

#making a vector of order names to set as row names
order_counts_tab <- otu_table(tax_glom(ps_srs, taxrank="Order"))
dim(order_counts_tab)
#[1] 338  990
order_tax_vec <- as.vector(tax_table(tax_glom(ps_srs, taxrank="Order"))[,4])
rownames(order_counts_tab) <- as.vector(order_tax_vec)
unclassified_tax_counts <- colSums(otu_table(ps_srs)) - colSums(order_counts_tab)
# and we'll add this row to our order count table:
order_and_unidentified_counts_tab <- rbind("Unclassified"=unclassified_tax_counts, order_counts_tab)
order_and_unidentified_counts_tab0 <- order_and_unidentified_counts_tab[rowSums(order_and_unidentified_counts_tab) > 0,]
order_and_unidentified_counts_tab <- order_and_unidentified_counts_tab0[,colSums(order_and_unidentified_counts_tab0) > 0]

#generate relative abundance table
order_proportions_tab <- apply(order_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
dim(order_proportions_tab)
#[1] 339  990

#write relative abundance table
order_proportions_tab_trim0 <- order_proportions_tab[,samples]
order_proportions_tab_trim <- order_proportions_tab_trim0[rowSums(order_proportions_tab_trim0) > 0,]
dim(order_proportions_tab_trim)
#[1] 339  988

#number of genera
nrow(order_proportions_tab_trim) - 1
#[1] 338

write.csv(order_proportions_tab_trim, "antarctic_soil_order_taxonomy_rel_abundance.csv")

##############################
##Create class-level dataset##
##############################

#making a vector of class names to set as row names
class_counts_tab <- otu_table(tax_glom(ps_srs, taxrank="Class"))
dim(class_counts_tab)
#[1] 132  990
class_tax_vec <- as.vector(tax_table(tax_glom(ps_srs, taxrank="Class"))[,3])
rownames(class_counts_tab) <- as.vector(class_tax_vec)
unclassified_tax_counts <- colSums(otu_table(ps_srs)) - colSums(class_counts_tab)
# and we'll add this row to our class count table:
class_and_unidentified_counts_tab <- rbind("Unclassified"=unclassified_tax_counts, class_counts_tab)
class_and_unidentified_counts_tab0 <- class_and_unidentified_counts_tab[rowSums(class_and_unidentified_counts_tab) > 0,]
class_and_unidentified_counts_tab <- class_and_unidentified_counts_tab0[,colSums(class_and_unidentified_counts_tab0) > 0]

#generate relative abundance table
class_proportions_tab <- apply(class_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
dim(class_proportions_tab)
#[1] 133  990

#write relative abundance table
class_proportions_tab_trim0 <- class_proportions_tab[,samples]
class_proportions_tab_trim <- class_proportions_tab_trim0[rowSums(class_proportions_tab_trim0) > 0,]
dim(class_proportions_tab_trim)
#[1] 133  988

#number of genera
nrow(class_proportions_tab_trim) - 1
#[1] 132

write.csv(class_proportions_tab_trim, "antarctic_soil_class_taxonomy_rel_abundance.csv")

###############################
##Create phylum-level dataset##
###############################

#making a vector of phylum names to set as row names
phylum_counts_tab <- otu_table(tax_glom(ps_srs, taxrank="Phylum"))
dim(phylum_counts_tab)
#[1] 51  990
phylum_tax_vec <- as.vector(tax_table(tax_glom(ps_srs, taxrank="Phylum"))[,2])
rownames(phylum_counts_tab) <- as.vector(phylum_tax_vec)
unclassified_tax_counts <- colSums(otu_table(ps_srs)) - colSums(phylum_counts_tab)
# and we'll add this row to our phylum count table:
phylum_and_unidentified_counts_tab <- rbind("Unclassified"=unclassified_tax_counts, phylum_counts_tab)
phylum_and_unidentified_counts_tab0 <- phylum_and_unidentified_counts_tab[rowSums(phylum_and_unidentified_counts_tab) > 0,]
phylum_and_unidentified_counts_tab <- phylum_and_unidentified_counts_tab0[,colSums(phylum_and_unidentified_counts_tab0) > 0]

#generate relative abundance table
phylum_proportions_tab <- apply(phylum_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
dim(phylum_proportions_tab)
#[1] 52  990

#write relative abundance table
phylum_proportions_tab_trim0 <- phylum_proportions_tab[,samples]
phylum_proportions_tab_trim <- phylum_proportions_tab_trim0[rowSums(phylum_proportions_tab_trim0) > 0,]
dim(phylum_proportions_tab_trim)
#[1] 52  988

#number of genera
nrow(phylum_proportions_tab_trim) - 1
#[1] 51

write.csv(phylum_proportions_tab_trim, "antarctic_soil_phylum_taxonomy_rel_abundance.csv")

###########################################
###Check number of samples per bioregion###
###########################################

#import information
bioclim <- read.csv("bioclim_bioregion_antarctic_study_curated_final.csv", row.names = 1)[samples,]
dim(bioclim)
#[1] 988  20

#check number of samples
dim(bioclim[bioclim$bioregion == '0',]) #[1] 24   20
dim(bioclim[bioclim$bioregion == '1',]) #[1] 6   20
dim(bioclim[bioclim$bioregion == '3',]) #[1] 80  20
dim(bioclim[bioclim$bioregion == '4',]) #[1] 8  20
dim(bioclim[bioclim$bioregion == '6',]) #[1] 169  20
dim(bioclim[bioclim$bioregion == '7',]) #[1] 191 20
dim(bioclim[bioclim$bioregion == '8',]) #[1] 49  20
dim(bioclim[bioclim$bioregion == '9',]) #[1] 189  20
dim(bioclim[bioclim$bioregion == '10',]) #[1] 154  20
dim(bioclim[bioclim$bioregion == '12',]) #[1] 11   20
dim(bioclim[bioclim$bioregion == '16',]) #[1] 107   20

#check that all empty
dim(bioclim[bioclim$bioregion == '2',])
dim(bioclim[bioclim$bioregion == '5',])
dim(bioclim[bioclim$bioregion == '11',])
dim(bioclim[bioclim$bioregion == '13',])
dim(bioclim[bioclim$bioregion == '14',])
dim(bioclim[bioclim$bioregion == '15',])
#all these bioregions have 0 samples

#check how many islands
dim(bioclim[bioclim$island == 'y',]) #[1] 142  20
dim(bioclim[bioclim$island == 'n',]) #[1] 846  20

#how many island soils are attributed to bioregions
dim(bioclim[bioclim$island == 'y' & bioclim$bioregion != '0',]) #[1] 118  20