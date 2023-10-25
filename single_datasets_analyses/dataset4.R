###############
###Dataset 4###
###############

#check, import in R and process sequences from dataset 4

#check how many reads
#for file in *1.fastq.gz 
#do
#count0=$(zcat ${file}|wc -l)
#count=$(( count0 / 4))
#combine="$file $count"
#echo $combine
#done

#run trimmomatic
#for file in *1.fastq.gz
#do
#filename="${file%_1.fastq.gz}"
#echo $filename
#java -jar /home/gilda/Trimmomatic-0.39/trimmomatic-0.39.jar PE $file ${filename}_2.fastq.gz ${filename}_R1_paired.fastq.gz ${filename}_R1_unpaired.fastq.gz ${filename}_R2_paired.fastq.gz ${filename}_R2_unpaired.fastq.gz CROP:285 ILLUMINACLIP:../Trimmomatic-0.39/adapters/all.fa:2:30:10
#done

#set working directory
setwd("/home/gilda/dataset4")

#import/save R environment
#save.image(file='dataset4_16S.RData')
#load('dataset4_16S.RData')

library("dada2", lib.loc="/home/gilda/R_library")
#https://benjjneb.github.io/dada2/tutorial.html

#path to fastq files
path <- "/home/gilda/dataset4"
list.files(path)

#upload files
fnFs <- sort(list.files(path, pattern="R1_paired", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_paired", full.names = TRUE))
#extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_R1_paired.fastq.gz"), `[`, 1)
sample.names

#trim sequences even if the quality is nice already
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, trimRight=c(20,50),maxEE=c(5,5), truncQ=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

#remove sample with no reads
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

#create quality check plots
dir.create("/home/gilda/dataset4/quality_plots")
library("gridExtra")
#raw forward
for(file in 1:length(filtFs)) {
  sample <- as.character(attributes(filtFs[file]))
  name <- paste("./quality_plots/quality_", sample, ".jpeg", sep="")
  print(name)
  jpeg(file=name)
  plot_Fs <- plotQualityProfile(fnFs[file])
  plot_Rs <- plotQualityProfile(fnRs[file])
  plot_Fs_filt <- plotQualityProfile(filtFs[file])
  plot_Rs_filt <- plotQualityProfile(filtRs[file])
  grid.arrange(plot_Fs, plot_Rs, plot_Fs_filt, plot_Rs_filt, ncol = 2, nrow = 2)
  dev.off()
}

#calculate error rates
errF <- learnErrors(filtFs, multithread=TRUE)
#103833377 total bases in 392209 reads from 8 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#106464903 total bases in 453561 reads from 9 samples will be used for learning the error rates.
#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

#sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#merge reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#[1] 16 25476

#chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 16359 bimeras out of 25476 input sequences.
dim(seqtab.nochim)
#[1]   16 9117
sum(seqtab.nochim)/sum(seqtab)
#[1] 0.5438023

#number of reads after each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
head(track)
write.table(track, file = "dataset4_read_count.txt", sep = "\t")

#download samples to assign taxonomy
#wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
#wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
taxa <- assignTaxonomy(seqtab.nochim, "/home/gilda/tax_silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/home/gilda/tax_silva/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa #removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)