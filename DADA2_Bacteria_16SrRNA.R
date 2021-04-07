################################################################
                # DADA2 pipeline tutorial #
################################################################
#Source: https://benjjneb.github.io/dada2/tutorial.html
#Instalation: https://benjjneb.github.io/dada2/dada-installation.html
#Author: Michaela de Melo
#Date: September 2020 
#Goal: Generate an amplicon sequence variant (ASV) table for 16S rRNA data (bacterial sequence)
################################################################

# Install binaries from Bioconductor
# Binary installation is available through the Bioconductor package repository. Binary installation of the current release version (1.8) requires R 3.5.0 and Bioconductor version 3.7:
# try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("dada2")

# Load dada2 package
library(dada2); packageVersion("dada2")
# documents explaining the pipeline
#vignette("dada2")

# Set working directory
setwd("C:/Users/...") #change it

# Define the following path variable so that it points to the extracted directory on your machine:
path <- "C:/Users/..." #change it
list.files(path)

# Read fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.
# Identify the pattern sample names
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), '[', 1)
sample.names

# We start by visualizing the quality profiles of the forward reads:
plotQualityProfile(fnFs[1:5])
# In gray-scale is a heat map of the frequency of each quality score at each base position. 
# The median quality score at each position is shown by the green line, and the quartiles of the quality score distribution by the orange lines. 
# The red line shows the scaled proportion of reads that extend to at least that position (this is more useful for other sequencing technologies, as Illumina reads are typically all the same lenghth, hence the flat red line).
# The forward reads are normally good quality. We generally advise trimming the last few nucleotides to avoid less well-controlled errors that can arise there. 

# Now we visualize the quality profile of the reverse reads:
plotQualityProfile(fnRs[6:11]) #change intervals to check all quality profiles

################################################################
# Filter and trim #

# Assign the filenames for the filtered fastq.gz files.
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# We'll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of "expected errors" allowed in a read, which is a better filter than simply averaging quality scores.
# truncQ = Truncate reads at the first instance of a quality score less than or equal to truncQ.
# truncLen = Truncate reads after truncLen bases. Reads shorter than this are discarded. Chech plotQuality Profiles to decide the limits
# trimLeft = The number of nucleotides to remove from the start of each read. 
# If both truncLen and trimLeft are provided, filtered reads will have length truncLen-trimLeft.
# Trimleft for primers removal:
# Primers= 515F (YRYRGTGCCAGCMGCCGCGGTAA) (19+ 4 added to increase complexity) and 806R (GGACTACHVGGGTWTCTAAT) (20 bases)
# There are often 3 bases before them so I will remove number of bases in primers +3 bases.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200), trimLeft=c(26,23),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out) #Chec number of reads in and out, the parameters should be edit (e.g. truncLen) in case the number of reads drop a lot reads in--> out
#Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later.
summary(out)
out

#Plot filtered
plotQualityProfile(filtFs, aggregate = TRUE)
plotQualityProfile(filtRs, aggregate = TRUE)


## Split forward and reverse sequences into two folders FWD and REV respectively
### Create new folders FWD and REV and sort files:
setwd("C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered")
dir.create("C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered/FWD")
list.files("C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered", pattern="_F_filt.fastq.gz")
FWDfiles <- sort(list.files("C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered", pattern="_F_filt.fastq.gz"))
FWDfiles
dir.create("C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered/REV")
list.files("C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered", pattern="_R_filt.fastq.gz")
REVfiles <- sort(list.files("C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered", pattern="_R_filt.fastq.gz"))
REVfiles
### Move files:
install.packages("filesstrings")
library(filesstrings)
file.move(FWDfiles, "C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered/FWD")
list.files("C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered/FWD")
file.move(REVfiles, "REV")
list.files("C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered/REV")

# File parsing
filtpathF <- "C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered/FWD"
filtpathR <- "C:/Users/micha/Desktop/DNA_JB/bac/Max/filtered/REV"
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.names
sample.namesR
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100) # Initialize random number generator for reproducibility


################################################################
# Learn the Error Rates #

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
plotErrors(errF, nominalQ = TRUE) #visualize it 

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
plotErrors(errR, nominalQ = TRUE) #visualize it 

################################################################
# Sample inference algorithm #
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

################################################################
# Merge paired ends #
#We now merge the forward and reverse reads together to obtain the full denoised sequences. 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

################################################################
# Construct an amplicon sequence variant table (ASV) table #
# This is a sample-by-sequence feature table valued by the number of times each sequence was observed in each sample.
seqtab <- makeSequenceTable(mergers)
dim(seqtab) #dimensions
table(nchar(getSequences(seqtab))) #Inspect distribution of sequence lengths

################################################################
# Remove chimeras #
# We now remove chimeric sequences by comparing each inferred sequence to the others in the table, and removing those that can be reproduced by stitching together two more abundant sequences.
seqtabnm.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtabnm.nochim)
sum(seqtabnm.nochim)/sum(seqtabnm) #considered good if higher than 0.7
saveRDS(seqtabnm.nochim, "seqtab_final.rds") # Save final sequence table (no mismatch and no chimera)

################################################################
# Track reads through the pipeline #
#As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(ddF, getN), sapply(ddR, getN), sapply(mergers, getN), rowSums(seqtabnm.nochim)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) by getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track # we should lose more reads at the filtering step and much less in the following steps.
write.csv(track, "Track_reads.csv") #save the output as a table

#Option to track reads through the pipeline (in case the one above doesn't work)
#track2= cbind(out, rowSums(seqtab),rowSums(seqtabnm.nochim))
#colnames(track2) <- c("input", "filtered", "merged", "chimera removal")
#rownames(track2) <- sample.names
#track2
#write.csv(track2, "track2.csv")

################################################################
# Assign taxonomy #
#Download database and into same folder as fasta files: https://benjjneb.github.io/dada2/training.html
tax <- assignTaxonomy(seqtabnm.nochim, "silva_nr_v138_train_set.fa.gz", multithread=TRUE)
dim(tax) # check dimension of data frame tax 
unname(head(tax, n=20)) # check head of data frame tax

# Save final taxonomy table (until genera)
saveRDS(tax, "C:/Users...") #Change path

#Add species:
#taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v132.fa.gz")

#Inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#END!
#Next steps are optional and doesn't use DADA2
################################################################
# Rarefaction and normalization using PHYLOSEQ #

library(phyloseq)
packageVersion("phyloseq") #check for updates

## Create the object phyloseq ##
setwd("/Users...") #change it

# Read in bacterial asv table
asv <- readRDS("seqtab_final.rds")
asv <- t(asv) # transpose (rows <-> columns)
dim(asv); nrow(asv); ncol(asv) 
# Read bacterial taxonomy table 
taxbact <- readRDS("tax_final.rds")
taxa_names(taxbact) <- paste0("Seq", seq(ntaxa(taxbact))) #rename sequences by simpler names

# Transform the bacterial taxa table to matrix
dim(taxbact); ncol(taxbact); nrow(taxbact) # sequences are rows, taxonomic ranks are columns (6)
taxbact_mat <- as.matrix(taxbact)
dim(taxbact_mat); ncol(taxbact_mat); nrow(taxbact_mat) 

# Check if asv_ba and taxbact_mat have the same number of taxa
nrow(asv); nrow(taxbact_mat) 


## Create SAMPLE DATAFRAME ##
#Preparing samdf = sample dataframe with features of samples, i.e. site, season, depth, etc
setwd("/Users...") #change path

asv <- readRDS("seqtab_final.rds")
samples.out <- rownames(asv)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)

max=read.csv("Sample_Data.csv", header = TRUE )
saveRDS(max, file = "samdf.rds")
samdf_phy <- readRDS("samdf.rds")
nrow(samdf_phy) #check number of rows

samdf_final=cbind(subject, samdf_phy, row.names=1)

## Create phyloseq class object ##
#with asv table, taxa table and sample metadata

physeq <- phyloseq(otu_table(asv, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxbact))


## Remove Archaea, Chloroplasts and mitochondria ##
library(dplyr)
clean <- physeq %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Order  != "Chloroplast" &
      Family   != "Mitochondria"
  )
dim(otu_table(clean)) #check how many reads were removed from ASV tables


## RAREFY phyloseq object  at minimum reads per sample ## 
sample_sums(clean) # number of reads in each sample

random.seed <- .Random.seed #To save the random number generation seed for reproducibility
ps_rare <- rarefy_even_depth(clean, sample.size = 12849,#change it
                             rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
#Check how many ASVS were removed during the rarefaction
attr(ps_rare,"seed") <- random.seed # to save the random number generation seed for reproducibility

# Check if saving the random number generation seed worked, ps_2 should give the same result as ps_rare2
.Random.seed <- attr(ps_rare,"seed") 
ps_rare_2 <- rarefy_even_depth(clean, sample.size = 12849,
                               rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

dim(otu_table(ps_rare)) 
# Save B-Ai rarefied ASV table for later manipulations
saveRDS(otu_table(ps_rare), "seqtab_bact_rarefied.rds")

## Plot rarefaction curve ##
library(vegan)
rarecurve(otu_table(clean), step=22, cex=0.5)
#OR
raremax <- 12849
col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)
out <- with(pars[1:26, ],
            rarecurve(otu_table(clean), sample= raremax, step = 22, col = col,
                      lty = lty, label = FALSE))

##########################################################
# Normalize: transform rarefied ASV table to proportions #

# Import rarefied ASV table if not already in the environment
seqtab_rare <- readRDS("seqtab_bact_rarefied.rds")
dim(seqtab_rare); nrow(seqtab_rare); ncol(seqtab_rare) # sequences are rows, and samples are columns

seqtab_norm <- transform_sample_counts(seqtab_rare, function(x) x / sum(x) ) # The function provided is applied to the vector of counts from each sample
dim(seqtab_norm) 
head(seqtab_norm) #asvs abundances are ratios from 0 to 1
sample_sums(seqtab_norm) #all sample sizes are equal to 1

saveRDS(seqtab_norm,'seqtab_bac_raref_normaliz.rds')
