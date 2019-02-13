# -----------------------------------------------------------------------------#
# Singapore and Malaysia Fungi Processing Raw Reads
# Author: Geoffrey Zahn
# Requirements: ITSx v 1.1b1
# -----------------------------------------------------------------------------#

# Prepare demultiplexed reads in Bash... ####
# Run ITSxpress on all fasta files 
# for i in ./Fwd/*.fastq.gz; do itsxpress --fastq $i  --outfile $i.FungalITS1.fastq.gz --region ITS1 --taxa Fungi -s --threads 10 --log $i.ITSx.log; done
# rename the files
# rename -v -e 's/^(.{7}).*/$1.fastq.gz/' *.FungalITS1.fastq.gz
# No ITS1 found in REV reads


# Load packages ####
# library("devtools")
# devtools::install_github("benjjneb/dada2")
# devtools::install_github("joey711/phyloseq")
# devtools::install_github("bryandmartin/corncob")
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)
library(DECIPHER)
library(decontam)
library(phangorn)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(gridExtra)
library(ggpubr)
library(dada2); packageVersion("dada2")
library(stringr)
library(purrr)
library(corncob)

# Load dada2 and prep ####


# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "./fastq/ITS1" # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = ".fastq.gz"))

sample.names <- sapply(strsplit(basename(fns), ".fastq"), `[`, 1)


# visualize a couple of fwd read quality profiles to help you decide reasonable filtration parameters
plotQualityProfile(fns[3:4])

source("~/Desktop/GIT_REPOSITORIES/My_Scripts/qual_maxEE_plot.R")
qualMaxEEplot("~/Desktop/GIT_REPOSITORIES/Data_Course/Data/Fastq_16S/F3D0_S188_L001_R1_001.fastq")

# Filter and trim ####
filts <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
# filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))



out <- filterAndTrim(fns, filts, # fnRs, filtRs,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=16, trimLeft = 15) # On Windows set multithread=FALSE


plotQualityProfile(c(fns[100],filts[100]))
?plotQualityProfile

out[order(out[,2]),]
no.pass = which(as.data.frame(out)$reads.out == 0)
out[no.pass,]





# learn error rates ####
# Since some samples may have had zero reads pass QC, reassign filtFs and filtRs
filts <- sort(list.files(filtpath, full.names = TRUE))

errF <- learnErrors(filts, multithread=TRUE, MAX_CONSIST = 20, nbases = 1e10)

# sanity check
plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)


# Dereplication ####
derep <- derepFastq(filts, verbose=TRUE)
# derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names

# Since some samples were removed (no reads passed QC), reassign sample.names
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1)
# sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derep) <- sample.names
# names(derepRs) <- sample.names


# Sample inference (special settings for IonTorrent) ####
dadaFs <- dada(derep, err=errF, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
# dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


# Make a sequence table ####
seqtab <- makeSequenceTable(dadaFs)


# Remove Chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# Track Reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[,1], sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$filter.loss = (track[,1]-track[,2])/track[,1]
write.csv(track, file = "./Output/read_counts_at_each_step.csv", row.names = TRUE)
write.csv(sample.names, "./Output/sample_names.csv")

# Save intermediate seqtable object
saveRDS(seqtab.nochim, "./seqtab.nochim.RDS")
seqtab.nochim <- readRDS("./seqtab.nochim.RDS")

write.csv(sample.names, "./sample_names.csv", row.names = FALSE)


# import metadata ####
meta = read.csv("./metadata.csv", header = TRUE)
row.names(meta) <- as.character(meta$IDENT)
row.names(seqtab.nochim)

# reorder metadata
meta = meta[order(row.names(meta)),]
meta = as.data.frame(meta)

# Find controlsamples (extraction negatives) and clean ####
# meta$controls <- meta$Island == "Blank"

# remove missing samples excluded due to poor QC 

# row.names(seqtab.nochim) <- map(strsplit(row.names(seqtab.nochim), split = "_"),1)
# good.samples = (row.names(meta) %in%  row.names(seqtab.nochim))
# meta <- (meta[good.samples,])
# rm(good.samples)


# find contaminants
# seqtab.nochim = readRDS(file = "./Output/clean_dada2_seqtable.RDS")
# contams = isContaminant(seqtab.nochim, neg = meta$controls, normalize = TRUE)
# table(contams$contaminant)  # No control samples passed QC...only 1 sequence made it through


# remove them
# seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
# seqtab.nochim = seqtab.nochim[meta$controls == FALSE,]
# meta = meta[meta$controls == FALSE,]

# Remove all seqs with fewer than 100 nucleotides ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# Assign Taxonomy ####

# Save intermediate seqtab and re-load before assigning taxonomy to reduce virtual memory usage

taxa <- assignTaxonomy(seqtab.nochim, "./Taxonomy/RDP_UNITE_with_outgroups.fasta", multithread=20)

# Save intermediate files
write.csv(as.data.frame(seqtab.nochim), file = "./Output/SeqTable_no-chimera.csv", row.names = TRUE, quote = FALSE)
saveRDS(seqtab.nochim, file = "./Output/clean_dada2_seqtable.RDS")
saveRDS(taxa, file = "./Output/RDP_Taxonomy_from_dada2.RDS")
seqs <- getSequences(seqtab.nochim)

# Hand off to Phyloseq ####

seqtab.nochim = readRDS(file = "./Output/clean_dada2_seqtable.RDS")

taxa = readRDS(file = "./Output/RDP_Taxonomy_from_dada2.RDS")

unique(taxa[,2])

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))


# Save RDS object for Phyloseq
saveRDS(ps, file = "./Output/clean_phyloseq_object.RDS")
ps = readRDS(file = "./Output/clean_phyloseq_object.RDS")
