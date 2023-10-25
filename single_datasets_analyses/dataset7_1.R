#################
###Dataset 7.1###
#################

#check, import in R and process sequences from dataset 7.1

#set working directory
setwd("/home/gilda/dataset7_1")

#import/save R environment
#save.image(file='dataset7_1_16S.RData')
#load('dataset7_1_16S.RData')

library("dada2")
#https://benjjneb.github.io/dada2/tutorial.html

#path to fastq files
path <- "/home/gilda/dataset7_1"
list.files(path)

#upload files
fn <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fn), ".fastq"), `[`, 1)
sample.names

#trim sequences
filt <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
names(filt) <- sample.names
out <- filterAndTrim(fn, filt, 
                     maxN=0, trimRight= 10, trimLeft=10, maxEE=5, rm.phix=TRUE, truncQ=0,
                     compress=TRUE, multithread=TRUE)
head(out)

#remove sample with no reads
filt <- filt[file.exists(filt)]

#calculate error rates
err <- learnErrors(filt, multithread=TRUE)
#103485428 total bases in 361985 reads from 9 samples will be used for learning the error rates.

#sample inference
inf <- dada(filt, err=err, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

#construct sequence table
seqtab <- makeSequenceTable(inf)
dim(seqtab)
#[1]    53 39388

#chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 7793 bimeras out of 40578 input sequences.
dim(seqtab.nochim)
#[1]    53 32785
sum(seqtab.nochim)/sum(seqtab)
#[1] 0.7858551

#number of reads after each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(inf, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, file = "dataset7_1_read_count.txt", sep = "\t")

#download samples to assign taxonomy
#wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
#wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
taxa <- assignTaxonomy(seqtab.nochim, "/home/gildav/tax_silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/home/gildav/tax_silva/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa #removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
