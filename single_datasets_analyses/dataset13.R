################
###Dataset 13###
################

#check, import in R and process sequences from dataset 13

#check how many reads
#for file in *_1.fastq.gz
#do
#count0=$(zcat ${file}|wc -l)
#count=$(( count0 / 4))
#combine="$file $count"
#echo $combine
#done

#trimmomatic
#for file in *_1.fastq.gz
#do
#filename="${file%_1*}"
#echo $filename
#java -jar /home/gilda/Trimmomatic-0.39/trimmomatic-0.39.jar PE $file ${filename}_2.fastq.gz ${filename}_R1_paired.fastq.gz ${filename}_R1_unpaired.fastq.gz ${filename}_R2_paired.fastq.gz ${filename}_R2_unpaired.fastq.gz ILLUMINACLIP:../Trimmomatic-0.39/adapters/all.fa:2:30:10
#done

#set working directory
setwd("/home/gilda/dataset13")

#import/save R environment
#save.image(file='dataset13_16S.RData')
#load('dataset13_16S.RData')

library("dada2")
#https://benjjneb.github.io/dada2/tutorial.html

#path to fastq files
path <- "/home/gilda/dataset13"
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
                     maxN=0, trimRight=c(5,20), maxEE=c(5,5), truncQ=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

#remove sample with no reads
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

#calculate error rates
errF <- learnErrors(filtFs, multithread=TRUE)
#106161292 total bases in 436152 reads from 9 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#111297140 total bases in 487194 reads from 10 samples will be used for learning the error rates.

#sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#plotErrors(errR, nominalQ=TRUE)

#merge reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#[1] 60 18813

#chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 13443 bimeras out of 18813 input sequences.
dim(seqtab.nochim)
#[1]   60 5370
sum(seqtab.nochim)/sum(seqtab)
#[1] 0.5502161

#number of reads after each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#rownames(track) <- sample.names
head(track)
write.table(track, file = "dataset13_read_count.txt", sep = "\t")

#download samples to assign taxonomy
#wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
#wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
taxa <- assignTaxonomy(seqtab.nochim, "/home/gildav/tax_silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/home/gildav/tax_silva/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa #removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)