###### PREPARE WORKSPACE #####
library(dada2)
library(ShortRead)

setwd("~/caribou")

# This function (from the dada2 tutorial) creates all possible orientations of an input DNA sequence.
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# This function counts the number of reads in which a primer is found, in any orientation.
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# This function calculates the number of unique occurrences (i.e., reads) in an object
getN <- function(x) sum(getUniques(x))

# This function replaces ASV names (in sequence letter) with temporary taxon names (temp1, temp2, temp3...)
# I use it just to make taxon names easier to work with during the filtering steps.
# It also returns the number of sequences that lack a kingdom-level classification.
rename_temp <- function(pseq_object){
  dna <- Biostrings::DNAStringSet(taxa_names(pseq_object))
  names(dna) <- taxa_names(pseq_object)
  pseq_object <- merge_phyloseq(pseq_object, dna)
  rm(dna)
  taxa_names(pseq_object) <- paste0("temp", seq(ntaxa(pseq_object)))
  print(ntaxa(pseq_object) - ntaxa(subset_taxa(pseq_object, Kingdom !="NA")))
  return(pseq_object)
}

# This function separates messy BLAST results into a cleaner table #
# There is probably a better way to do this but I don't know it! 
organize_blast <- function(input_table){
  new_table <- input_table %>%
    tidyr::separate(., "sscinames", into=c("name1", "name2"),
                    sep="(?= internal transcribed)", remove=FALSE)
  
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= 5.8S ribosomal)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= 18S r)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= small subunit)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= genes for)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= mRNA)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= external transcribed)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= genom)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= chromosome)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= small ribosomal)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= 16S r)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= gene for)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= 16S-like)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= partial)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table <- new_table %>%
    tidyr::separate(., "name1", into=c("name1", "name2"),
                    sep="(?= DNA for)", remove=FALSE)
  new_table$name2[is.na(new_table$name2)] <- ""
  new_table$seqtype <- paste0(new_table$name2, new_table$seqtype)
  new_table$name2 <- NULL
  
  new_table$sscinames <- new_table$name1
  new_table$name1 <- NULL
  return(new_table)
}

# Import initial sample data
sample_data <- read.csv("raw_data/sample_data.csv", row.names="SampleID")

###### INITIAL SEQUENCE IMPORT AND PROCESSING ########
## [1] Confirm that primers have been removed from all the reads. ####
set.seed(100)

input_paths <- list()
fnFs <- list()
fnRs <- list()
sampleNames <- list()

FWD.orients <- list()
REV.orients <- list()
primer_check <- list()

primers <- data.frame(
  X16S = c("GTGYCAGCMGCCGCGGTAA", "CCGYCAATTYMTTTRAGTTT"),
  X18S = c("CCAGCASCYGCGGTAATTCC", "ACTTTCGTTCTTGATYRA"),
  FITS = c("GAACGCAGCRAANNGYGA", "TCCTCCGCTTATTGATATGC"),
  PITS = c("ATGCGATACTTGGTGTGAAT", "GACGCTTCTCCAGACTACAAT")
)

# Import sequence files and check for any residual primers in the data.
for(i in c("X16S", "X18S", "FITS", "PITS")){
  # Name the folder where the sequence files are stored. 
  input_paths[[i]] <- file.path("sequences", i, "cutadapt")
  
  # Identify forward (Fs) and reverse (Rs) reads.
  fnFs[[i]] <- sort(list.files(input_path[[i]], pattern="_R1.fastq"))
  fnRs[[i]] <- sort(list.files(input_path[[i]], pattern="_R2.fastq"))

  # Identify sample names (should be the beginning of the filename).
  sampleNames[[i]] <- sapply(strsplit(fnFs[[i]], "_"), '[', 2)
  fnFs[[i]] <- file.path(input_path, fnFs[[i]])
  fnRs[[i]] <- file.path(input_path, fnRs[[i]])
  
  # Determine orientations of all primers
  FWD.orients[[i]] <- allOrients(primers[1,i])
  REV.orients[[i]] <- allOrients(primers[2,i])
  
  primer_check[[i]] <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[i]][[1]]), 
                             FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[i]][[1]]), 
                             REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[i]][[1]]), 
                             REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[i]][[1]]))
}

rm(FWD.orients, REV.orients)

## [2] Truncate and quality-filter reads. ####
filt_paths <- list()
filtFs <- list()
filtRs <- list()
out <- list()

for(i in c("X16S", "X18S", "FITS", "PITS")){
  # Create a new directory where filtered reads will be stored.
  filt_paths[[i]] <- file.path(input_paths[[i]], "filtered")
  if(!file_test("-d", filt_paths[[i]])) dir.create(filt_paths[[i]])
  
  # Name filtered sequences with the same names as input sequences with _filt added.
  filtFs[[i]] <- file.path(filt_paths[[i]], paste0(sampleNames[[i]], "_F_filt.fastq.gz"))
  filtRs[[i]] <- file.path(filt_paths[[i]], paste0(sampleNames[[i]], "_R_filt.fastq.gz"))
  
  # 16S, 18S, and fungal ITS2 reads are truncated at 220 bp (forward read) and 180 bp (reverse read)
  if(i !="PITS"){
    out[[i]] <- filterAndTrim(fnFs[[i]], filtFs[[i]], fnRs[[i]], filtRs[[i]], 
                              truncLen=c(220,180), maxN=0, maxEE=c(2,2), truncQ=2,
                              rm.phix=TRUE, compress=TRUE, multithread=TRUE)}
  
  # Plant ITS2 reads are truncated at 260 bp (forward read) and 220 bp (reverse read)
  if(i=="PITS"){
    out[[i]] <- filterAndTrim(fnFs[[i]], filtFs[[i]], fnRs[[i]], filtRs[[i]], 
                              truncLen=c(260,220), maxN=0, maxEE=c(2,2), truncQ=2, 
                              rm.phix=TRUE, compress=TRUE, multithread=TRUE)}
  
  out[[i]] <- as.data.frame(out[[i]])
  rownames(out[[i]]) <- sampleNames
  colnames(out[[i]]) <- c("input", "filtered")
  
  # Remove samples that no longer have any reads after the filtering steps.
  bad.samples <- as.data.frame(out)
  bad.samples <- subset(bad.samples, filtered==0)
  bad.samples <- rownames(bad.samples)
  
  bad.samplesF <- file.path(filt_path, paste0(bad.samples, "_F_filt.fastq.gz"))
  bad.samplesR <- file.path(filt_path, paste0(bad.samples, "_R_filt.fastq.gz"))
  
  sampleNames <- sampleNames[!(sampleNames %in% bad.samples)]
  filtFs <- filtFs[!(filtFs %in% bad.samplesF)]
  filtRs <- filtRs[!(filtRs %in% bad.samplesR)]
  
  rm(bad.samples, bad.samplesF, bad.samplesR)
  
}  

## [3] Calculate error rates, merge paired-end reads, and create sequence tables. ####
derepFs <- list()
derepRs <- list()

errF <- list()
errR <- list()

dadaFs <- list()
dadaRs <- list()

mergers <- list()

seqtabAll <- list()

for(i in c("X16S", "X18S", "FITS", "PITS")){
  # Dereplicate sequences.
  derepFs[[i]] <- derepFastq(filtFs[[i]], verbose=TRUE)
  derepRs[[i]] <- derepFastq(filtRs[[i]], verbose=TRUE)
  
  # Assign sample names.
  names(derepFs[[i]]) <- sampleNames[[i]]
  names(derepRs[[i]]) <- sampleNames[[i]]
  
  # Calculate error probabilities.
  errF[[i]] <- learnErrors(filtFs[[i]], multithread=TRUE)
  errR[[i]] <- learnErrors(filtRs[[i]], multithread=TRUE)
  
  # Use error probabilities to determine exact sequence variants.
  dadaFs[[i]] <- dada(derepFs[[i]], err=errF[[i]], multithread=TRUE, pool=TRUE)
  dadaRs[[i]] <- dada(derepRs[[i]], err=errR[[i]], multithread=TRUE, pool=TRUE)
  
  # Merged paired-end reads and construct an OTU table.
  mergers[[i]] <- mergePairs(dadaFs[[i]], derepFs[[i]], dadaRs[[i]], derepRs[[i]])
  seqtabAll[[i]] <- makeSequenceTable(mergers[[i]])
}

## [4] Remove chimeras ####
# Examine the distribution of sequence lengths. Most sequences should be ~372bp.
table(nchar(getSequences(seqtabAll[["X16S"]])))

table(nchar(getSequences(seqtabAll[["X18S"]])))

table(nchar(getSequences(seqtabAll[["FITS"]])))

table(nchar(getSequences(seqtabAll[["PITS"]])))

# Remove chimeras and summarize the amount of reads lost in each processing step.
seqtabNoC <- list()
for(i in c("X16S", "X18S", "FITS", "PITS")){
  seqtabNoC[[i]] <- removeBimeraDenovo(seqtabAll[[i]], multithread=TRUE)  
  
  if(i=="X16S"){
    read_counts[[i]] <- cbind(out[[i]], sapply(dadaFs[[i]], getN), sapply(dadaRs[[i]], getN), sapply(mergers[[i]], getN),
                              rowSums(seqtabAll[[i]]), rowSums(seqtabNoC[[i]]))
    colnames(read_counts[[i]]) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "seqtab", "nochimera")
    rownames(read_counts[[i]]) <- sampleNames[[i]]
  }
  
  if(i !="X16S"){
    read_counts[[i]] <- cbind(sapply(dadaFs[[i]], getN), sapply(dadaRs[[i]], getN), sapply(mergers[[i]], getN),
                              rowSums(seqtabAll[[i]]), rowSums(seqtabNoC[[i]]))
    colnames(read_counts[[i]]) <- c("denoisedF", "denoisedR", "merged", "seqtab", "nochimera")
    rownames(read_counts[[i]]) <- sampleNames
    read_counts[[i]] <- merge(out[[i]], read_counts[[i]], by=0, all=TRUE)
    rownames(read_counts[[i]]) <- read_counts[[i]]$Row.names
    read_counts[[i]]$Row.names <- NULL
  }
  
}

###### ASSIGN AND CLEAN TAXONOMY #####
## [1] Bacterial 16S ####
# This was downloaded from https://benjjneb.github.io/dada2/training.html.
fastaRef <- "~/Documents/z_taxonomy_databases/rdp_train_set_18.fa.gz"

# Assign taxonomy to the sequence table.
taxtabNoC <- assignTaxonomy(seqtabNoC[["16S"]], refFasta=fastaRef, multithread=TRUE)

# Check to see that it worked.
unname(head(taxtabNoC[["16S"]]))

# Import into phyloseq object and clean workspace
pseq.16S.raw <- phyloseq(otu_table(seqtabNoC[["X16S"]], taxa_are_rows=FALSE),
                         sample_data(sample_data),
                         tax_table(taxtabNoC))

rm(taxtabNoC, fastaRef)

# a - Remove chloroplasts. ####
# Identify chloroplasts.
temp <- as.data.frame(tax_table(pseq.16S.raw))
temp$chloroplast <- temp$Class == "Chloroplast"
temp <- subset(temp, chloroplast=="TRUE") # n = 14

# Check the abundance of these sequences identified as chloroplasts.
df <- subset(as.data.frame(taxa_sums(pseq.16S.raw)),
             rownames(as.data.frame(taxa_sums(pseq.16S.raw))) %in% rownames(temp))

# Remove chloroplasts from the data frame.
badTaxa <- rownames(temp)
goodTaxa <- setdiff(taxa_names(pseq.16S.raw), badTaxa)
pseq.16S.filter <- prune_taxa(goodTaxa, pseq.16S.raw)

rm(temp, goodTaxa, badTaxa, df)

# b - Remove archaea ####
# Convert to relative abundance and subset to archaea.
archaea <- transform_sample_counts(pseq.16S.filter, function(x) 100*x/sum(x))
archaea <- subset_taxa(archaea, Kingdom=="Archaea")
archaea <- prune_samples(sample_sums(archaea) > 0, archaea)

# Check the presence and abundance of archaeal ASVs.
tax_table(archaea) # Only four archaeal ASVs, present in only 17 samples
rowMeans(otu_table(archaea)) # All at extremely low abundances

# Remove archaea from analysis
pseq.16S.filter <- subset_taxa(pseq.16S.filter, !(Kingdom %in% c("Archaea")))

# Remove sequences that are not in the correct length distribution.
temp <- otu_table(pseq.16S.raw)
table(nchar(getSequences(temp)))

temp <- temp[,nchar(colnames(temp)) %in% seq(369,375)]
pseq.16S.filter <- subset_taxa(pseq.16S.filter, taxa_names(pseq.16S.filter) %in% colnames(temp))
rm(temp)

# Create a list to start storing information about things that were removed from the data.
ctrl_data <- list()
ctrl_data$archaea <- archaea
rm(archaea)

# c - Remove taxa unclassified at the phylum level ####
# Rename sequences to temp1, temp2, temp3... and count number of ASVs unclassified at kingdom level
pseq.16S.filter <- rename_temp(pseq.16S.filter) # 0

# Check for taxa not classified at the PHYLUM level.
nophylum.16S <- subset_taxa(pseq.16S.filter, Phylum !="NA")
nophylum.16S <- as.data.frame(cbind(tax_table(nophylum.16S)))
nophylum.16S <- subset_taxa(pseq.16S.filter,
                                 !(taxa_names(pseq.16S.filter) %in% rownames(nophylum.16S)))

# Export those taxa for a BLAST search.
seqtab.nochim <- as.matrix(otu_table(nophylum.16S))
seqs <- refseq(nophylum.16S)
ids <- taxa_names(nophylum.16S)
db_out <- data.frame(ids=ids, seqs=seqs, count=colSums(seqtab.nochim))
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))

nophylum.16S <- list(db_out, fasta)
rm(db_out, fasta, ids, seqs, seqtab.nochim)

write.csv(nophylum.16S[[1]], "unclassified_seqs/nophylum_16s.csv") # n = 46
writeFasta(nophylum.16S[[2]], file="unclassified_seqs/nophylum_16S.fna")

# (BLAST search of the FASTA file - done remotely) #

# Calculate the percentage of reads that lack a phylum-level classification.
100*sum(nophylum.16S[[1]]$count)/sum(sample_sums(pseq.16S.filter)) # 0.0525% of reads
nrow(nophylum.16S[[1]]) # representing 46 ASVs

# Remove bacterial taxa that lack a phylum-level classification.
pseq.16S.filter <- subset_taxa(pseq.16S.filter, is.na(Phylum)=="FALSE")

ctrl_data$nophylum <- list()
ctrl_data$nophylum$X16S <- nophylum.16S
rm(nophylum.16S)

## [2] Eukaryotic 18S ####
# a - Run taxonomic annotation separately in QIIME ####
# Export sequences to assign taxonomy using the more updated QIIME database.
db_out <- data.frame(
  ids = paste0("original", c(1:nrow(seqtabNoC))),
  seqs = colnames(seqtabNoC)
)

fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
writeFasta(fasta, file = "original.18S.sequences.fna")

# (run analysis in QIIME) #

# b - Import and clean new taxonomy table ####
tax.qiime.cut.corrected <- read.csv("intermediate_data/X18S_raw_taxonomy.tsv", sep="\t")
tax.qiime.cut.corrected <- tidyr::separate(tax.qiime.cut.corrected,
                                           "Taxon",
                                           into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                           sep="; ")

taxtabNoC <- tax.qiime.cut.corrected
rownames(taxtabNoC) <- colnames(seqtabNoC[["X18S"]])
taxtabNoC$Feature.ID <- NULL
taxtabNoC$Confidence <- NULL

# Remove identifier prefixes
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("d__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("k__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("p__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("c__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("o__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("f__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("g__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("s__", "", x)}))

# c - Clean some kingdom- and phylum-level identifications and import to phyloseq ####
taxtabNoC$Kingdom <- as.character(taxtabNoC$Kingdom)
taxtabNoC$Kingdom[taxtabNoC$Phylum %in% c("Ascomycota", "Basidiomycota", "Chytridiomycota",
                                          "Cryptomycota", "Mucoromycota", "Neocallimastigomycota",
                                          "Zoopagomycota", "LKM15", "Aphelidea")] <- "Fungi"
taxtabNoC$Kingdom[taxtabNoC$Phylum %in% c("Chlorophyta", "Dinoflagellata", "Klebsormidiophyceae",
                                          "Charophyta", "Ochrophyta")] <- "Algae"
taxtabNoC$Kingdom[taxtabNoC$Phylum %in% c("Vertebrata", "Arthropoda", "Nematozoa",
                                          "Tardigrada", "Rotifera", "Platyhelminthes",
                                          "Holozoa")] <- "Animalia"
taxtabNoC$Kingdom[taxtabNoC$Phylum %in% c("Archamoebae", "Apicomplexa", "Parabasalia",
                                          "Cercozoa", "Amoebozoa", "Ciliophora",
                                          "Gracilipodida", "Protalveolata", "Cavosteliida",
                                          "SAR", "MAST-12", "Peronosporomycetes",
                                          "Bicosoecida")] <- "Protista"
taxtabNoC$Kingdom[taxtabNoC$Phylum %in% c("Phragmoplastophyta")] <- "Plant"

taxtabNoC$Phylum[taxtabNoC$Kingdom=="Eukaryota" & taxtabNoC$Phylum==""]

df <- unique(taxtabNoC$X18S[,c("Kingdom","Phylum")])
df <- unique(taxtabNoC[,c("Kingdom","Phylum", "Class")])
df <- subset(unique(taxtabNoC[,c("Kingdom","Phylum", "Class","Order","Family","Genus","Species")]),
             Kingdom=="Eukaryota")
rm(df, tax.qiime.cut.corrected)

rownames(taxtabNoC) <- colnames(seqtabNoC[["X18S"]])
taxtabNoC <- as.matrix(taxtabNoC)

# Check to see that it worked.
unname(head(taxtabNoC))

# Subset sample data to the samples with 18S reads
sample_data_temp <- subset(sample_data, rownames(sample_data) %in% rownames(seqtabNoC[["X18S"]]))

# Import into phyloseq object and clean workspace
pseq.18S.raw <- phyloseq(otu_table(seqtabNoC[["X18S"]], taxa_are_rows=FALSE),
                         sample_data(sample_data_temp),
                         tax_table(taxtabNoC))
rm(taxtabNoC, sample_data_temp)

# d - Use BLAST to identify taxa not assigned at the phylum level ####
pseq.18S.filter <- pseq.18S.raw

# Rename sequences to temp1, temp2, temp3... and count number of ASVs unclassified at kingdom level
pseq.18S.filter <- rename_temp(pseq.18S.filter) # 0

# Check for taxa unclassified at the phylum level
nophylum.18S <- subset_taxa(pseq.18S.filter, Phylum !="NA")
nophylum.18S <- as.data.frame(cbind(tax_table(nophylum.18S)))
nophylum.18S <- subset_taxa(pseq.18S.filter,
                            !(taxa_names(pseq.18S.filter) %in% rownames(nophylum.18S)))

# Export those taxa for a BLAST search.
seqtab.nochim <- as.matrix(otu_table(nophylum.18S))
seqs <- refseq(nophylum.18S)
ids <- taxa_names(nophylum.18S)
db_out <- data.frame(ids=ids, seqs=seqs, count=colSums(seqtab.nochim))
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))

nophylum.18S <- list(db_out, fasta)
rm(db_out, fasta, ids, seqs, seqtab.nochim)

write.csv(nophylum.18S[[1]], "unclassified_seqs/nophylum_18S.csv") # n = 46
writeFasta(nophylum.18S[[2]], file="unclassified_seqs/nophylum_18S.fna")

# (BLAST search of the FASTA file - done remotely) #
# TERMINAL: ~/Downloads/ncbi-blast-2.10.0+/bin/blastn -query nophylum_18S.fna -db nt -out nophylum_18S_annotations.txt -outfmt "10 std sscinames stitle" -max_target_seqs 10 -remote

# Calculate the percentage of reads that lack a phylum-level classification.
sum(nophylum.18S[[1]]$count)/sum(sample_sums(pseq.18S.filter)) # 18S = 4.77% of reads lack phylum classification

# Load and organize BLAST results
nophylum.18S[[3]] <- read.csv("intermediate_data/X18S_nophylum_blast.txt", header=FALSE,
                         col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                                     "qstart", "qend", "sstart", "send",
                                     "evalue", "bitscore", "staxids", "sscinames",
                                     "a1", "a2", "a3", "a4", "a5", "a6", "a7"))
qseqids.18S <- forcats::fct_inorder(unique(nophylum.18S[[3]]$qseqid))

# Remove taxa that had no significant BLAST similarity
df <- data.frame(problem_taxa = nophylum.18S[[1]]$ids, 
                 BLAST = nophylum.18S[[1]]$ids %in% qseqids.18S)
df <- subset(df, BLAST=="FALSE") # 11 taxa were missing from the BLAST search b/c no significant matches were found

pseq.18S.filter <- subset_taxa(pseq.18S.filter,
                               !(taxa_names(pseq.18S.filter) %in% df$problem_taxa))
rm(df)

# Clean BLAST results.
nophylum.18S[[3]]$seqtype <- paste0(nophylum.18S[[3]]$a1, nophylum.18S[[3]]$a2, 
                                    nophylum.18S[[3]]$a3, nophylum.18S[[3]]$a4,
                                    nophylum.18S[[3]]$a5, nophylum.18S[[3]]$a6)
nophylum.18S[[3]][,c("a1","a2","a3","a4","a5","a6","a7")] <- NULL

nophylum.18S[[3]] <- organize_blast(nophylum.18S[[3]])

nophylum.18S[[3]] <- nophylum.18S[[3]] %>%
  tidyr::separate(., "sscinames", into=c("genus", "species"),
                  sep=" ", remove=FALSE)

df <- unique(nophylum.18S[[3]]$genus)

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="D.scoparium"] <- "scoparium"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="D.scoparium"] <- "Dicranum"

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="L.pyriforme"] <- "pyriforme"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="L.pyriforme"] <- "Lycoperdon"

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="S.racemosum"] <- "racemosum"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="S.racemosum"] <- "Syncephalastrum"

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="D.albidus"] <- "albidus"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="D.albidus"] <- "Dryocosmus"

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="N.frontalis"] <- "frontalis"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="N.frontalis"] <- "Neocallimastix"

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="R.squarrosus"] <- "squarrosus"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="R.squarrosus"] <- "Rhytidiadelphus"

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="H.splendens"] <- "splendens"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="H.splendens"] <- "Hylocomium"

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="P.adiantoides"] <- "adiantoides"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="P.adiantoides"] <- "Piptadenia"

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="J.leiantha"] <- "leiantha"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="J.leiantha"] <- "Jungermannia"

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="F.pusilla"] <- "pusilla"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="F.pusilla"] <- "Fabronia"

nophylum.18S[[3]]$species[nophylum.18S[[3]]$genus=="P.undulatum"] <- "undulatum"
nophylum.18S[[3]]$genus[nophylum.18S[[3]]$genus=="P.undulatum"] <- "Pittosporum"

# Separate BLAST hits that have a clean genus/species name from BLAST results that have messy hits
unclass.help <- subset(nophylum.18S[[3]], genus %in% c("Uncultured", "uncultured", "Fish", "Marine",
                                                       "PREDICTED:", "Fungal", "Unicellular", "[Candida]",
                                                       "Soil", "Moss") |
                         is.na(species)==TRUE)

unclass.good <- subset(nophylum.18S[[3]], !(rownames(nophylum.18S[[3]]) %in% rownames(unclass.help)))

# Get full list of available taxon information for the "good" genera
genera <- unique(unclass.good$genus)
families <- taxize::tax_name(genera, get=c("superkingdom", "kingdom", "phylum",
                                           "subphylum", "class", "subclass",
                                           "order", "family", "genus"),
                             db="ncbi", message=FALSE)

families$db <- NULL
families$genus_final <- families$genus
families$genus <- NULL

unclass.good$query <- unclass.good$genus
unclass.good <- merge(unclass.good,
                      families, by="query", all=TRUE)

unclass.good[,c("genus","species","query")] <- NULL

# Manually clean the names of the "bad" genera
unclass.help$genus[unclass.help$genus=="Uncultured"] <- "uncultured"

unclass.help.uncultured <- subset(unclass.help, genus=="uncultured")
unclass.help.other <- subset(unclass.help, !(rownames(unclass.help) %in% rownames(unclass.help.uncultured)))

lookup <- unique(unclass.help.uncultured$species)
result <- taxize::tax_name(lookup, get=c("superkingdom", "kingdom", "phylum",
                                           "subphylum", "class", "subclass",
                                           "order", "family", "genus"),
                             db="ncbi", message=FALSE)
result$db <- NULL
result$genus_final <- result$genus
result$genus <- NULL
result.manual <- subset(result, is.na(superkingdom)=="TRUE")$query

# Rename "uncultured" to lowest taxonomic level #
unclass.help.uncultured$query <- unclass.help.uncultured$species
unclass.help.uncultured <- merge(unclass.help.uncultured,
                                 result, by="query", all=TRUE)

unclass.help.uncultured$phylum[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="Amoeboza"] <- "Amoebozoa"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="Amoeboza"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="Amoeboza"] <- NA

unclass.help.uncultured$superkingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="eukaryote"] <- "Eukaryote"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="eukaryote"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="eukaryote"] <- NA

unclass.help.uncultured$kingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="alveolate"] <- "Protista"
unclass.help.uncultured$phylum[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="alveolate"] <- "Alveolata"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="alveolate"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="alveolate"] <- NA

unclass.help.uncultured$kingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="dinoflagellate"] <- "Protista"
unclass.help.uncultured$phylum[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="dinoflagellate"] <- "Dinoflagellata"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="dinoflagellate"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="dinoflagellate"] <- NA

unclass.help.uncultured$superkingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="marine"] <- "Eukaryota"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="marine"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="marine"] <- NA

unclass.help.uncultured$superkingdom[unclass.help.uncultured$genus=="uncultured" & 
                                 unclass.help.uncultured$species=="soil" &
                                 grepl("Uncultured soil eukaryote", unclass.help.uncultured$sscinames)] <- "Eukaryota"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & 
                                unclass.help.uncultured$species=="soil" &
                                grepl("Uncultured soil eukaryote", unclass.help.uncultured$sscinames)] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & 
                                  unclass.help.uncultured$species=="soil" &
                                  grepl("Uncultured soil eukaryote", unclass.help.uncultured$sscinames)] <- NA

unclass.help.uncultured$kingdom[unclass.help.uncultured$genus=="uncultured" & 
                                       unclass.help.uncultured$species=="soil" &
                                       grepl("Uncultured soil fungus", unclass.help.uncultured$sscinames)] <- "Fungi"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & 
                                unclass.help.uncultured$species=="soil" &
                                grepl("Uncultured soil fungus", unclass.help.uncultured$sscinames)] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & 
                                  unclass.help.uncultured$species=="soil" &
                                  grepl("Uncultured soil fungus", unclass.help.uncultured$sscinames)] <- NA

unclass.help.uncultured$phylum[unclass.help.uncultured$genus=="uncultured" & 
                      unclass.help.uncultured$species %in% c("cercozoa", "Cercozoa", "cercozoan", "Cercozoan")] <- "Cercozoa"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & 
                     unclass.help.uncultured$species %in% c("cercozoa", "Cercozoa", "cercozoan", "Cercozoan")] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & 
                      unclass.help.uncultured$species %in% c("cercozoa", "Cercozoa", "cercozoan", "Cercozoan")] <- NA

unclass.help.uncultured$kingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="flagellate"] <- "Protista"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="flagellate"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="flagellate"] <- NA

unclass.help.uncultured$kingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="fungus"] <- "Fungi"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="fungus"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="fungus"] <- NA

unclass.help.uncultured$superkingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="freshwater"] <- "Eukaryote"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="freshwater"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="freshwater"] <- NA

unclass.help.uncultured$superkingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="microeukaryote"] <- "Eukaryote"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="microeukaryote"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="microeukaryote"] <- NA

unclass.help.uncultured$kingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="haptophyte"] <- "Algae"
unclass.help.uncultured$phylum[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="haptophyte"] <- "Haptophyta"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="haptophyte"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="haptophyte"] <- NA

unclass.help.uncultured$kingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="stramenopile"] <- "Protista"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="stramenopile"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="stramenopile"] <- NA

unclass.help.uncultured$superkingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="Banisveld"] <- "Eukaryote"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="Banisveld"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="Banisveld"] <- NA

unclass.help.uncultured$superkingdom[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="eukaryotic"] <- "Eukaryote"
unclass.help.uncultured$genus[unclass.help.uncultured$genus=="uncultured" & unclass.help.uncultured$species=="eukaryotic"] <- NA
unclass.help.uncultured$species[is.na(unclass.help.uncultured$genus)=="TRUE" & unclass.help.uncultured$species=="eukaryotic"] <- NA

# Rename the other random ones to lowest taxonomic level #
unclass.help.other[,c("superkingdom", "kingdom", "phylum", "subphylum", "class", "subclass",
                      "order", "family", "genus_final")] <- rep(NA)
unclass.help.other <- tibble::add_column(unclass.help.other, .before=1,
                                         query=unclass.help.other$genus)
unique(unclass.help.other[,c("genus","species")])

unclass.help.other$kingdom[unclass.help.other$genus=="Fish" & unclass.help.other$species=="environmental"] <- "Animalia"
unclass.help.other$phylum[unclass.help.other$genus=="Fish" & unclass.help.other$species=="environmental"] <- "Vertebrata"
unclass.help.other$class[unclass.help.other$genus=="Fish" & unclass.help.other$species=="environmental"] <- "Fish"
unclass.help.other$genus[unclass.help.other$genus=="Fish" & unclass.help.other$species=="environmental"] <- NA
unclass.help.other$species[is.na(unclass.help.other$genus)=="TRUE" & unclass.help.other$species=="environmental"] <- NA

unclass.help.other$kingdom[unclass.help.other$genus=="Marine" & unclass.help.other$species=="Metazoa"] <- "Animalia"
unclass.help.other$genus[unclass.help.other$genus=="Marine" & unclass.help.other$species=="Metazoa"] <- NA
unclass.help.other$species[is.na(unclass.help.other$genus)=="TRUE" & unclass.help.other$species=="Metazoa"] <- NA

unclass.help.other$kingdom[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Salmo"] <- "Animalia"
unclass.help.other$phylum[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Salmo"] <- "Vertebrata"
unclass.help.other$class[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Salmo"] <- "Actinopterygii"
unclass.help.other$order[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Salmo"] <- "Salmoniformes"
unclass.help.other$family[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Salmo"] <- "Salmonidae"
unclass.help.other$genus_final[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Salmo"] <- "Salmo"
unclass.help.other$species[unclass.help.other$species=="Salmo"] <- "salar"

unclass.help.other$kingdom[unclass.help.other$genus=="Fungal" & unclass.help.other$species=="sp."] <- "Fungi"
unclass.help.other$genus[unclass.help.other$genus=="Fungal" & unclass.help.other$species=="sp."] <- NA
unclass.help.other$species[is.na(unclass.help.other$genus)=="TRUE" & unclass.help.other$species=="sp."] <- NA

unclass.help.other$superkingdom[unclass.help.other$genus=="Unicellular" & unclass.help.other$species=="eukaryote"] <- "Eukaryote"
unclass.help.other$genus[unclass.help.other$genus=="Unicellular" & unclass.help.other$species=="eukaryote"] <- NA
unclass.help.other$species[is.na(unclass.help.other$genus)=="TRUE" & unclass.help.other$species=="eukaryote"] <- NA

unclass.help.other$phylum[unclass.help.other$genus=="Soil" & unclass.help.other$species=="amoeba"] <- "Amoebozoa"
unclass.help.other$genus[unclass.help.other$genus=="Soil" & unclass.help.other$species=="amoeba"] <- NA
unclass.help.other$species[is.na(unclass.help.other$genus)=="TRUE" & unclass.help.other$species=="sp."] <- NA

unclass.help.other$kingdom[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Hylobates"] <- "Animalia"
unclass.help.other$phylum[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Hylobates"] <- "Vertebrata"
unclass.help.other$class[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Hylobates"] <- "Mammalia"

unclass.help.other$kingdom[unclass.help.other$genus=="[Candida]"] <- "Fungi"
unclass.help.other$genus_final[unclass.help.other$genus=="[Candida]"] <- "Candida"

unclass.help.other$kingdom[unclass.help.other$genus=="Moss" & unclass.help.other$species=="environmental"] <- "Plant"
unclass.help.other$phylum[unclass.help.other$genus=="Moss" & unclass.help.other$species=="environmental"] <- "Embryophyta"
unclass.help.other$order[unclass.help.other$genus=="Moss" & unclass.help.other$species=="environmental"] <- "Bryophyta"
unclass.help.other$genus[unclass.help.other$genus=="Moss" & unclass.help.other$species=="environmental"] <- NA
unclass.help.other$species[is.na(unclass.help.other$genus)=="TRUE" & unclass.help.other$species=="environmental"] <- NA

unclass.help.other$kingdom[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Mangifera"] <- "Plant"
unclass.help.other$genus_final[unclass.help.other$genus=="PREDICTED:" & unclass.help.other$species=="Mangifera"] <- "Mangifera"

# Put the "bad" BLAST results back together again and standardize names #
unclass.help <- rbind(unclass.help.uncultured,
                      unclass.help.other)

unclass.help[,c("genus","species","query")] <- NULL

rm(unclass.help.uncultured, unclass.help.other)

unclass.help$superkingdom <- rep("Eukaryota")
unclass.help$phylum[unclass.help$phylum == "Perkinsozoa"] <- "Apicomplexa"
unclass.help$kingdom[unclass.help$family=="Colpodellaceae"] <- "Apicomplexa"

unclass.help$subphylum[unclass.help$phylum == "Evosea"] <- "Evosea"
unclass.help$phylum[unclass.help$phylum == "Evosea"] <- "Amoebozoa"
unclass.help$kingdom[unclass.help$phylum %in% c("Cercozoa", "Amoebozoa", "Apicomplexa")] <- "Protista"
unclass.help$kingdom[unclass.help$phylum %in% c("Haptophyta", "Bacillariophyta")] <- "Algae"

unclass.help$phylum[unclass.help$genus=="Candida"] <- "Ascomycota"
unclass.help$class[unclass.help$genus=="Candida"] <- "Saccharomycetes"
unclass.help$order[unclass.help$genus=="Candida"] <- "Saccharomycetales"
unclass.help$family[unclass.help$genus=="Candida"] <- "Saccharomycetaceae"

unclass.help$order[unclass.help$genus=="Mangifera"] <- "Sapindales"
unclass.help$family[unclass.help$genus=="Mangifera"] <- "Anacardioceae"

unclass.help$class[unclass.help$order=="Bryophyta"] <- "Bryophyta"

# Put ALL the unclassified taxa / BLAST results back together again
unclass.18S <- rbind(unclass.good, unclass.help)

nophylum.18S_tbl <- unclass.18S %>%
  tibble::as_tibble(.) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::group_split(.keep=TRUE)

for(i in c(1:length(nophylum.18S_tbl))){
  nophylum.18S_tbl[[i]] <- nophylum.18S_tbl[[i]][order(-nophylum.18S_tbl[[i]]$bitscore), ]
  nophylum.18S_tbl[[i]] <- rbind(nophylum.18S_tbl[[i]], c(NA, NA, NA, NA, NA, NA, NA, NA,
                                                           NA, NA, NA, NA, NA, NA, NA, NA, NA))
}

# e - Continue cleaning BLAST results to ensure accurate taxonomic classifications ####
# Write the BLAST results to a csv file
nophylum.18S_tbl <- dplyr::bind_rows(nophylum.18S_tbl)
write.csv(nophylum.18S_tbl, "nophylum.18S.tbl.csv")

# Manually inspect the CSV file to choose the best blast hit for each ASV.
nophylum.18S_audit <- read.csv("intermediate_data/X18S_nophylum_audited.csv")

orig.tax.table <- as.data.frame(cbind(tax_table(pseq.18S.filter)))
orig.tax.table <- subset(orig.tax.table, !rownames(orig.tax.table) %in% nophylum.18S_audit$qseqid)

rownames(nophylum.18S_audit) <- nophylum.18S_audit$qseqid
nophylum.18S_audit$qseqid <- NULL

orig.tax.table <- tibble::add_column(orig.tax.table, .before=1,
                                     superkingdom=rep("Eukaryota"))

for(i in c(1:ncol(nophylum.18S_audit))){
  nophylum.18S_audit[,i] <- as.character(nophylum.18S_audit[,i])
}

nophylum.18S_audit$kingdom[nophylum.18S_audit$phylum=="Chlorophyta"] <- "Algae"

nophylum.18S_audit$subclass <- NULL
nophylum.18S_audit$subphylum <- NULL

nophylum.18S_audit$phylum[nophylum.18S_audit$genus=="Entamoeba"] <- "Archamoebae"
nophylum.18S_audit$class[nophylum.18S_audit$genus=="Entamoeba"] <- "Entamoebida"
nophylum.18S_audit$order[nophylum.18S_audit$genus=="Entamoeba"] <- "Entamoebida"
nophylum.18S_audit$family[nophylum.18S_audit$genus=="Entamoeba"] <- "Entamoebida"

nophylum.18S_audit$family[nophylum.18S_audit$order=="Eugregarinorida"] <- "Eugregarinorida"
nophylum.18S_audit$order[nophylum.18S_audit$order=="Eugregarinorida"] <- "Gregarinasina"

nophylum.18S_audit$family[nophylum.18S_audit$genus=="Chloromonas"] <- "Chlorophyceae"
nophylum.18S_audit$order[nophylum.18S_audit$genus=="Chloromonas"] <- "Chlorophyceae"

nophylum.18S_audit$family[nophylum.18S_audit$genus=="Neocystis"] <- "Chlorophyceae"
nophylum.18S_audit$order[nophylum.18S_audit$genus=="Neocystis"] <- "Chlorophyceae"

nophylum.18S_audit$family[nophylum.18S_audit$genus=="Chlorellales"] <- "Trebouxiophyceae"
nophylum.18S_audit$order[nophylum.18S_audit$genus=="Chlorellales"] <- "Trebouxiophyceae"

nophylum.18S_audit$family[nophylum.18S_audit$genus=="Chloroidium"] <- "Trebouxiophyceae"
nophylum.18S_audit$order[nophylum.18S_audit$genus=="Chloroidium"] <- "Trebouxiophyceae"

nophylum.18S_audit$family[nophylum.18S_audit$order=="Sarcoptiformes"] <- "Acari"

nophylum.18S_audit$family[nophylum.18S_audit$family=="Aphelenchoididae"] <- "Rhabditida"

nophylum.18S_audit$order[nophylum.18S_audit$order=="Microbotryales"] <- NA

nophylum.18S_audit$family[nophylum.18S_audit$family=="Sainouridae"] <- NA
nophylum.18S_audit$genus[nophylum.18S_audit$genus=="Acantholus"] <- NA

nophylum.18S_audit$kingdom[nophylum.18S_audit$kingdom=="Viridiplantae"] <- "Plant"

nophylum.18S_audit$kingdom[nophylum.18S_audit$phylum=="Streptophyta"] <- "Phragmoplastophyta"
nophylum.18S_audit$kingdom[nophylum.18S_audit$phylum=="Streptophyta"] <- "Phragmoplastophyta"

nophylum.18S_audit$order[nophylum.18S_audit$class=="Bryopsida"] <- "Bryophyta"
nophylum.18S_audit$class[nophylum.18S_audit$order=="Bryophyta"] <- "Embryophyta"
nophylum.18S_audit$class[nophylum.18S_audit$order=="Polytrichales"] <- "Embryophyta"

nophylum.18S_audit$species <- rep(NA)
colnames(nophylum.18S_audit) <- colnames(orig.tax.table)

new.tax.table <- rbind(nophylum.18S_audit, orig.tax.table)
colnames(new.tax.table)[1] <- "Domain"
new.tax.table <- as.matrix(new.tax.table)

tax_table(pseq.18S.filter) <- new.tax.table

temp <- as.data.frame(cbind(tax_table(pseq.18S.filter)))
for(i in c(1:ncol(temp))){temp[,i] <- as.character(temp[,i])}
temp$Kingdom[temp$Kingdom=="Phragmoplastophyta"] <- "Plant"
temp$Phylum[temp$Order=="Bryophyta"] <- "Phragmoplastophyta"
temp$Phylum[temp$Class=="Embryophyta"] <- "Phragmoplastophyta"
temp$Kingdom[temp$Genus=="Blastocystis"] <- "Protista"
temp$Phylum[temp$Genus=="Blastocystis"] <- "Bigyra"
temp$Class[temp$Genus=="Blastocystis"] <- "Blastocystea"
temp$Order[temp$Genus=="Blastocystis"] <- "Blastocystida"
temp$Family[temp$Genus=="Blastocystis"] <- "Blastocystidae"
temp$Class[temp$Phylum=="Mucoromycota" & temp$Class=="Incertae_Sedis"] <- "Mucoromycota"
temp$Class[temp$Phylum=="Cryptomycota" & temp$Class=="Incertae_Sedis"] <- NA
temp$Order[temp$Phylum=="Cryptomycota" & temp$Order=="Incertae_Sedis"] <- NA
temp$Family[temp$Phylum=="Cryptomycota" & temp$Family=="Incertae_Sedis"] <- NA
temp$Class[temp$Genus=="Kraken"] <- "Imbricatea"
temp$Order[temp$Genus=="Kraken"] <- "Krakenida"
temp$Family[temp$Genus=="Kraken"] <- "Krakenidae"
temp$Order[temp$Order=="Incertae_Sedis" & is.na(temp$Genus)=="TRUE"] <- NA
temp$Order[temp$Order=="Incertae_Sedis" & temp$Genus=="Incertae_Sedis"] <- NA
temp$Order[temp$Genus=="Microsporomyces"] <- NA
temp$Family[temp$Genus=="Microsporomyces"] <- NA
temp$Family[temp$Family=="Incertae_Sedis" & temp$Order=="Saccharomycetales"] <- "Saccharomycetales"
temp$Family[temp$Family=="Incertae_Sedis" & is.na(temp$Order)=="TRUE"] <- NA
temp$Genus[temp$Genus=="Incertae_Sedis"] <- NA
temp$Order[temp$Family=="Acari"] <- "Acari"
temp$Class[temp$Genus=="Bracteacoccus"] <- "Chlorophyceae"
temp$Order[temp$Genus=="Bracteacoccus"] <- "Chlorophyceae"
temp$Family[temp$Genus=="Bracteacoccus"] <- "Chlorophyceae"
temp$Family[temp$Genus=="Leucosporidium"] <- "Luecosporidiaceae"
temp$Order[temp$Class=="Embryophyta" & temp$Order=="Embryophyta"] <- NA
temp$Family[temp$Class=="Embryophyta" & temp$Family=="Embryophyta"] <- NA
temp$Family[temp$Order=="Agaricales" & temp$Family=="Agaricales"] <- NA
temp$Kingdom[temp$Genus=="Dicksoniaceae"] <- "Plant"
temp$Phylum[temp$Genus=="Dicksoniaceae"] <- "Phragmoplastophyta"
temp$Class[temp$Genus=="Dicksoniaceae"] <- "Embryophyta"
temp$Order[temp$Genus=="Dicksoniaceae"] <- "Bryophyta"
temp$Family[temp$Genus=="Dicksoniaceae"] <- "Bryophyta"
temp$Phylum[temp$Genus=="Phaseoleae" & temp$Family=="Amb-18S-504"] <- "Chlorophyta"
temp$Class[temp$Genus=="Phaseoleae" & temp$Family=="Amb-18S-504"] <- "Chlorophyceae"
temp$Order[temp$Genus=="Phaseoleae" & temp$Family=="Amb-18S-504"] <- "Chlorophyceae"
temp$Family[temp$Genus=="Phaseoleae" & temp$Family=="Amb-18S-504"] <- "Chlorophyceae"
temp$Genus[temp$Genus=="Phaseoleae" & temp$Family=="Chlorophyceae"] <- "Neocystis"
temp$Genus[temp$Species=="uncultured_Trebouxia"] <- "Trebouxia"
temp$Genus[temp$Species=="uncultured_Chloromonas"] <- "Chloromonas"
temp$Genus[temp$Species=="uncultured_Coccomyxa"] <- "Coccomyxa"
temp$Genus[temp$Species=="uncultured_Apatococcus"] <- "Apatococcus"
temp$Genus[temp$Species=="uncultured_Stichococcus"] <- "Stichococcus"
temp$Genus[temp$Species=="uncultured_Xylochloris"] <- "Xylochloris"

temp$Order[temp$Genus=="Gregarina"] <- NA
temp$Family[temp$Genus=="Gregarina"] <- NA
temp$Species[temp$Genus=="Gregarina"] <- NA
temp$Genus[temp$Genus=="Gregarina"] <- NA
tax_table(pseq.18S.filter) <- as.matrix(temp)

rm(df, families, new.tax.table, nophylum.18S_audit, nophylum.18S_tbl, orig.tax.table, result,
   temp, unclass.18S, unclass.good, unclass.help, lookup, qseqids.18S, result.manual)

# f - Filter taxa w/o kingdom classification and taxa called as Mammalia ####
# Remove taxa that lack a kingdom-level classification (beyond "eukaryote")
# Calculate the percentage of reads that lack a kingdom-level classification.
sum(sample_sums(subset_taxa(pseq.18S.filter, is.na(Kingdom)=="TRUE"))) / sum(sample_sums(pseq.18S.filter)) # Only 0.1108375% of reads lack a kingdom classification.
ntaxa(subset_taxa(pseq.18S.filter, is.na(Kingdom)=="TRUE")) # Representing 65 ASVs.

# Goes from 3640 ASVs to 3575 ASVs
pseq.18S.filter <- subset_taxa(pseq.18S.filter, is.na(Kingdom)=="FALSE")

# Check for and remove taxa that were called as mammals.
temp.mammals <- subset_taxa(pseq.18S.filter, Class=="Mammalia")
ntaxa(temp.mammals) # 38 ASVs are called as mammals
sum(sample_sums(temp.mammals)) # 160,082 reads are classified as mammals
sum(sample_sums(temp.mammals)) / sum(sample_sums(pseq.18S.filter)) # Mammals = 8.16% of reads are classified as mammals

db_out <- data.frame(ids=taxa_names(temp.mammals),
                     seqs=refseq(temp.mammals),
                     abund=taxa_sums(temp.mammals))

fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
writeFasta(fasta, file="mammals.fna")
rm(fasta)

# TERMINAL: ~/Downloads/ncbi-blast-2.10.0+/bin/blastn -query mammals.fna -db nt -out mammals_18S_annotations.txt -outfmt "10 std sscinames stitle" -max_target_seqs 10 -remote

pseq.18S.filter <- subset_taxa(pseq.18S.filter, !(taxa_names(pseq.18S.filter) %in% taxa_names(temp.mammals)))

rm(temp.mammals, db_out)

saveRDS(pseq.18S.raw, "pseq_objects/pseq.18S.raw.rds")
saveRDS(pseq.18S.filter, "pseq_objects/pseq.18S.filter.rds")
saveRDS(nophylum.18S, "pseq_objects/nophylum.18S.rds")

## [3] Fungal ITS ####
# a - Run taxonomic annotation using DADA2. ####
# Identify the reference database you will use for assigning taxonomy. This was downloaded from https://benjjneb.github.io/dada2/training.html.
fastaRef <- "~/Documents/z_taxonomy_databases/sh_general_release_2021.10.05.tgz"

# Assign taxonomy to your sequence table.
taxtabNoC <- assignTaxonomy(seqtabNoC[["FITS"]], refFasta=fastaRef, multithread=TRUE)

# Check to see that it worked.
unname(head(taxtabNoC))

# Remove underscore taxon level identifiers
taxtabNoC <- as.data.frame(taxtabNoC)
taxtabNoC$seq <- rownames(taxtabNoC)
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("k__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("p__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("c__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("o__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("f__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("g__", "", x)}))
taxtabNoC <- data.frame(lapply(taxtabNoC, function(x) {gsub("s__", "", x)}))
rownames(taxtabNoC) <- taxtabNoC$seq
taxtabNoC$seq <- NULL
taxtabNoC <- as.matrix(taxtabNoC)

# Import initial sample data
sample_data_temp <- subset(sample_data, rownames(sample_data) %in% rownames(seqtabNoC))

# Import into phyloseq object and clean workspace
pseq.FITS.raw <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE),
                          sample_data(sample_data_temp),
                          tax_table(taxtabNoC))

# b - Identify taxa not assigned at the phylum level ####
pseq.FITS.filter <- pseq.FITS.raw

# Rename sequences to temp1, temp2, temp3... and count number of ASVs unclassified at kingdom level
pseq.FITS.filter <- rename_temp(pseq.FITS.filter) # 0

# Check taxa unclassified at the phylum level  #
nophylum.FITS <- subset_taxa(pseq.FITS.filter, Phylum !="NA")
nophylum.FITS <- as.data.frame(cbind(tax_table(nophylum.FITS)))
nophylum.FITS <- subset_taxa(pseq.FITS.filter,
                             !(taxa_names(pseq.FITS.filter) %in% rownames(nophylum.FITS)))

# Export those taxa for a BLAST search.
seqtab.nochim <- as.matrix(otu_table(nophylum.FITS))
seqs <- refseq(nophylum.FITS)
ids <- taxa_names(nophylum.FITS)
db_out <- data.frame(ids=ids, seqs=seqs, count=colSums(seqtab.nochim))
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))

nophylum.FITS <- list(db_out, fasta)
rm(db_out, fasta, ids, seqs, seqtab.nochim)

write.csv(nophylum.FITS[[1]], "unclassified_seqs/nophylum_FITS.csv") # n = 46
writeFasta(nophylum.FITS[[2]], file="unclassified_seqs/nophylum_FITS.fna")

# (BLAST search of the FASTA file - done remotely) #

# Calculate the percentage of reads that lack a phylum-level classification.
sum(nophylum.FITS[[1]]$count)/sum(sample_sums(pseq.FITS.filter)) # 4.34% of reads
nrow(nophylum.FITS[[1]]) # representing 620 ASVs

# c - Load and clean BLAST results ####
nophylum.FITS[[3]] <- read.csv("unclassified_seqs/FITS_nophylum_blast.txt", header=FALSE,
                         col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                                     "qstart", "qend", "sstart", "send",
                                     "evalue", "bitscore", "staxids", "sscinames",
                                     "a1", "a2", "a3", "a4", "a5", "a6", "a7"))
qseqids.FITS <- forcats::fct_inorder(unique(nophylum.FITS[[3]]$qseqid))

# See if any taxa were missing from the BLAST search.
df <- data.frame(problem_taxa = nophylum.FITS[[1]]$ids, BLAST = nophylum.FITS[[1]]$ids %in% qseqids.FITS)
subset(df, BLAST=="FALSE") # temp2712
rm(df)

# "No significant similarity found" for this sequence, so let's remove it
pseq.FITS.filter <- subset_taxa(pseq.FITS.filter, taxa_names(pseq.FITS.filter) != "temp2712")

nophylum.FITS[[3]]$seqtype <- paste0(nophylum.FITS[[3]]$a1, nophylum.FITS[[3]]$a2, nophylum.FITS[[3]]$a3, nophylum.FITS[[3]]$a4,
                                     nophylum.FITS[[3]]$a5, nophylum.FITS[[3]]$a6, nophylum.FITS[[3]]$a7)
nophylum.FITS[[3]][,c("a1","a2","a3","a4","a5","a6","a7")] <- NULL

nophylum.FITS[[3]] <- organize_blast(nophylum.FITS[[3]])

nophylum.FITS[[3]] <- nophylum.FITS[[3]] %>%
  tidyr::separate(., "sscinames", into=c("genus", "species"),
                  sep=" ", remove=FALSE)

# Separate BLAST hits that have a clean genus/species name from BLAST results that have messy hits
check.unclass.FITS <- subset(nophylum.FITS[[3]], genus %in% c("Uncultured", "uncultured", "Unclassified", "unclassified",
                                                        "Unidentified", "Viridiplantae", "T.jamesii", "[Myrmecia]") |
                               is.na(species)==TRUE)

# Remove these results from the broader taxon lookup. 
nophylum.FITS[[3]] <- subset(nophylum.FITS[[3]], !(genus %in% check.unclass.FITS$genus) &
                         is.na(species)==FALSE)

# Get full list of available taxon information for the "good" genera
genera <- unique(nophylum.FITS[[3]]$genus)

genera1 <- genera[c(1:50)]
genera2 <- genera[c(51:100)]
genera3 <- genera[c(101:length(genera))]

families <- taxize::tax_name(genera1, get=c("superkingdom", "kingdom", "phylum",
                                           "subphylum", "class", "subclass",
                                           "order", "family", "genus"),
                             db="ncbi", message=FALSE)
familes2 <- taxize::tax_name(genera2, get=c("superkingdom", "kingdom", "phylum",
                                            "subphylum", "class", "subclass",
                                            "order", "family", "genus"),
                             db="ncbi", message=FALSE)
familes3 <- taxize::tax_name(genera3, get=c("superkingdom", "kingdom", "phylum",
                                            "subphylum", "class", "subclass",
                                            "order", "family", "genus"),
                             db="ncbi", message=FALSE)

families <- rbind(families, familes2, familes3)
rm(familes2, familes3, genera1, genera2, genera3)

nophylum.FITS[[3]]$query <- nophylum.FITS[[3]]$genus
nophylum.FITS[[3]] <- merge(nophylum.FITS[[3]], families, by="query", all=TRUE)
rm(genera, families)

# Manually clean the names of the "bad" genera
check.unclass.FITS$genus[check.unclass.FITS$genus=="Uncultured"] <- "uncultured"
check.unclass.FITS$genus[check.unclass.FITS$genus=="Unclassified"] <- "unclassified"

#
manual.fix <- subset(check.unclass.FITS, species %in% c("fungus", "alga", "eukaryote", "soil", "Fungus", "Dikarya",
                                                        "isolate", "environmental", "bisecta", "bacterium", "basal", 
                                                        "basidiomycete", "organism") | genus %in% c("T.jamesii"))

check.unclass.FITS <- subset(check.unclass.FITS, !(species %in% manual.fix$species))
check.unclass.FITS <- subset(check.unclass.FITS, !(genus %in% "T.jamesii"))

genera2 <- unique(check.unclass.FITS$species)

families2 <- taxize::tax_name(genera2, get=c("superkingdom", "kingdom", "phylum",
                                             "subphylum", "class", "subclass",
                                             "order", "family", "genus"),
                              db="ncbi", message=FALSE)

check.unclass.FITS$query <- check.unclass.FITS$species
check.unclass.FITS <- merge(check.unclass.FITS, families2, by="query", all=TRUE)
nophylum.FITS[[3]] <- rbind(nophylum.FITS[[3]], check.unclass.FITS)
rm(check.unclass.FITS, families2, genera2)

# Rename the few that needed manual help...
manual.fix <- manual.fix %>%
  tibble::add_column(query = manual.fix$genus,
                     .before=1)
manual.fix[,c("db", "superkingdom", "kingdom", "phylum", "subphylum", "class", "subclass", "order", "family", "genus.y")] <- rep(NA)
colnames(manual.fix)[16] <- "genus.x"

# 2 - uncultured fungus
new_df <- subset(manual.fix, query=="uncultured" & species=="fungus")
new_df$superkingdom <-  "Eukaryota"
new_df$kingdom <-  "Fungi"

# 3 - uncultured alga
a <- subset(manual.fix, query=="uncultured" & species=="alga")
a$superkingdom <- "Eukaryota"
a$kingdom <- "Alga"
new_df <- rbind(new_df, a)

# 6 - uncultured eukaryote
a <- subset(manual.fix, query=="uncultured" & species=="eukaryote")
a$superkingdom <- "Eukaryota"
new_df <- rbind(new_df, a)

# 10 - uncultured soil
a <- subset(manual.fix, query=="uncultured" & species=="soil")
a$superkingdom <- "Eukaryota"
a$kingdom <- "Fungi"
new_df <- rbind(new_df, a)

# 18 - T. jamesii
a <- subset(manual.fix, query=="T.jamesii")
a$superkingdom <- "Eukaryota"
a$kingdom <- "Alga"
a$phylum <- "Chlorophyta"
a$class <- "Trebouxiophyceae"
a$order <- "Trebouxiales"
a$family <- "Trebouxiaceae"
a$genus.y <- "Trebouxia"
new_df <- rbind(new_df, a)

# 24 - uncultured Fungus
a <- subset(manual.fix, query=="uncultured" & species=="Fungus")
a$superkingdom <- "Eukaryota"
a$kingdom <- "Fungi"
new_df <- rbind(new_df, a)

# 25 - uncultured Dikarya
a <- subset(manual.fix, query=="uncultured" & species=="Dikarya")
a$superkingdom <- "Eukaryota"
a$kingdom <- "Fungi"
new_df <- rbind(new_df, a)

# 30 - Unidentified isolate
a <- subset(manual.fix, query=="Unidentified" & species=="isolate")
a$superkingdom <- "Eukaryota"
new_df <- rbind(new_df, a)

# 33 - Viridiplantae environmental
a <- subset(manual.fix, query=="Viridiplantae" & species=="environmental")
a$superkingdom <- "Eukaryota"
a$kingdom <- "Viridiplantae"
new_df <- rbind(new_df, a)

# 34 - [Myrmecia] bisecta
a <- subset(manual.fix, query=="[Myrmecia]" & species=="bisecta")
a$superkingdom <- "Eukaryota"
a$kingdom <- "Alga"
a$phylum <- "Chlorophyta"
a$class <- "Trebouxiophyceae"
a$order <- "Trebouxiales"
a$family <- "Trebouxiaceae"
a$genus.y <- "Myrmecia"
new_df <- rbind(new_df, a)

# 35 - uncultured bacterium
a <- subset(manual.fix, query=="uncultured" & species=="bacterium")
a$superkingdom <- "Bacteria"
a$kingdom <- "Bacteria"
new_df <- rbind(new_df, a)

# uncultured basal
a <- subset(manual.fix, query=="uncultured" & species=="basal")
a$superkingdom <- "Eukaryota"
a$kingdom <- "Fungi"
new_df <- rbind(new_df, a)

# uncultured basidiomycete
a <- subset(manual.fix, query=="uncultured" & species=="basidiomycete")
a$superkingdom <- "Eukaryota"
a$kingdom <- "Fungi"
a$phylum <- "Basidiomycota"
new_df <- rbind(new_df, a)

# uncultured organism
a <- subset(manual.fix, query=="uncultured" & species=="organism")
new_df <- rbind(new_df, a)

rm(a)

nophylum.FITS[[3]] <- rbind(nophylum.FITS[[3]], new_df)

rm(manual.fix, new_df)

# Split into a separate table for every ASV.
FITS_tbl <- nophylum.FITS[[3]] %>%
  tibble::as_tibble(.) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::group_split(.keep=TRUE)

# Create a taxonomy table for the fungal ASVs based on their BLAST results.
FITS_list <- list()

unclass.FITS_tax_table <- data.frame(matrix(ncol=8))
colnames(unclass.FITS_tax_table) <- c("qseqid", "kingdom", "phylum", "class", "order", "family", "genus", "species")

for(i in 1:length(FITS_tbl)){
  FITS_list[[i]] <- as.data.frame(FITS_tbl[[i]])
  
  FITS_list[[i]] <- subset(FITS_list[[i]],
                           bitscore >= max(FITS_list[[i]]$bitscore) & pident > 70)
  
  unclass.FITS_tax_table[i,1] <- as.character(unique(FITS_list[[i]]$qseqid))
  
  if(length(unique(FITS_list[[i]]$kingdom))==1){
    unclass.FITS_tax_table[i,2] <- unique(FITS_list[[i]]$kingdom)}
  
  if(length(unique(FITS_list[[i]]$phylum))==1){
    unclass.FITS_tax_table[i,3] <- unique(FITS_list[[i]]$phylum)}
  
  if(length(unique(FITS_list[[i]]$class))==1){
    unclass.FITS_tax_table[i,4] <- unique(FITS_list[[i]]$class)}
  
  if(length(unique(FITS_list[[i]]$order))==1){
    unclass.FITS_tax_table[i,5] <- unique(FITS_list[[i]]$order)}
  
  if(length(unique(FITS_list[[i]]$family))==1){
    unclass.FITS_tax_table[i,6] <- unique(FITS_list[[i]]$family)}
  
  if(length(unique(FITS_list[[i]]$genus))==1){
    unclass.FITS_tax_table[i,7] <- unique(FITS_list[[i]]$genus)}
  
  if(length(unique(FITS_list[[i]]$species))==1){
    unclass.FITS_tax_table[i,8] <- unique(FITS_list[[i]]$species)}
}

# d - Remove ASVs assigned to algae, plants, or bacteria ####
# Make a list of "bad taxa" (algae, plants) and "good taxa" (fungi)
badTaxa <- subset(unclass.FITS_tax_table, kingdom %in% c("Viridiplantae", "Alga"))
goodTaxa <- subset(unclass.FITS_tax_table, !(qseqid %in% badTaxa$qseqid))

# Determine the total number and abundance of sequence being removed
FITS.otu <- as.data.frame(otu_table(pseq.FITS.filter))
FITS.otu_sums <- as.data.frame(colSums(FITS.otu))
FITS.otu_sums$badTaxa <- rownames(FITS.otu_sums) %in% badTaxa$qseqid

# Check how much of the data is being removed
nrow(badTaxa) # 407 taxa assigned as plants or algae
sum(subset(FITS.otu_sums, badTaxa=="TRUE")$`colSums(FITS.otu)`) # 84,963 reads worth of "bad taxa"
sum(subset(FITS.otu_sums, badTaxa=="TRUE")$`colSums(FITS.otu)`) / sum(sample_sums(pseq.FITS.filter)) # Total 3.4% of reads

# Remove them from sample
pseq.FITS.filter <- subset_taxa(pseq.FITS.filter, !(taxa_names(pseq.FITS.filter) %in% badTaxa$qseqid))

# The BLAST search gave additional annotations to about 6 fungal taxa that were unclassified in the UNITE database.
# All of these taxa were relatively low abundance, so leaving them as "uncl fungi" seems sufficient.
FITS.goodtax <- subset(unclass.FITS_tax_table, qseqid %in% taxa_names(pseq.FITS.filter))
rownames(FITS.goodtax) <- FITS.goodtax$qseqid
FITS.goodtax <- merge(FITS.goodtax, taxa_sums(pseq.FITS.filter), by=0, all.x=TRUE)
rownames(FITS.goodtax) <- FITS.goodtax$Row.names
FITS.goodtax$Row.names <- NULL

x <- as.data.frame(cbind(tax_table(pseq.FITS.filter)))
FITS.goodtax <- merge(FITS.goodtax, x, by=0, all.x=TRUE)

rm(FITS.otu_sums, FITS.otu, badTaxa, goodTaxa, FITS_list, FITS_tbl, x, FITS.goodtax)

nophylum.FITS$classified_taxa <- unclass.FITS_tax_table
rm(unclass.FITS_tax_table)

# Reclassify a few taxa that UNITE appeared to label incorrectly.
tax.tab.fungi <- as.data.frame(cbind(tax_table(pseq.FITS.filter)))

# Change Austrolecia and Rhizoplaca to Lecanora
tax.tab.fungi$Genus[tax.tab.fungi$Genus=="Austrolecia"] <- "Lecanora"
tax.tab.fungi$Genus[tax.tab.fungi$Genus=="Rhizoplaca"] <- "Lecanora"

tax.tab.fungi$Family[tax.tab.fungi$Genus=="Lecanora"] <- "Lecanoraceae"
tax.tab.fungi$Order[tax.tab.fungi$Genus=="Lecanora"] <- "Lecanorales"
tax.tab.fungi$Class[tax.tab.fungi$Genus=="Lecanora"] <- "Lecanoromycetes"
tax.tab.fungi$Phylum[tax.tab.fungi$Genus=="Lecanora"] <- "Ascomycota"
tax.tab.fungi$Species[tax.tab.fungi$Genus=="Lecanora"] <- NA

# Change Phylliscum to Rhinocladiella
tax.tab.fungi$Genus[tax.tab.fungi$Genus=="Phylliscum"] <- "Rhinocladiella"
tax.tab.fungi$Family[tax.tab.fungi$Genus=="Rhinocladiella"] <- "Herpotrichiellaceae"
tax.tab.fungi$Order[tax.tab.fungi$Genus=="Rhinocladiella"] <- "Chaetothyriales"
tax.tab.fungi$Class[tax.tab.fungi$Genus=="Rhinocladiella"] <- "Eurotiomycetes"
tax.tab.fungi$Phylum[tax.tab.fungi$Genus=="Rhinocladiella"] <- "Ascomycota"

tax.tab.fungi <- as.matrix(tax.tab.fungi)
tax_table(pseq.FITS.filter) <- tax.tab.fungi

rm(tax.tab.fungi)

## [4] Plant ITS ####
# a - Export all sequences for a BLAST search. ####
db_out <- data.frame(
  ids = paste0("seq", seq(ncol(seqtabNoC))),
  seqs = colnames(seqtabNoC),
  abund = colSums(seqtabNoC)
)

fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
writeFasta(fasta, file = "plant_seqs.fna")

# Do the taxon analysis via BLAST in the terminal.

# b - Import and clean BLAST results. ####
plant_annotations <- read.csv("intermediate_data/PITS_blast_annotations.csv", header=FALSE,
                              col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                                          "qstart", "qend", "sstart", "send",
                                          "evalue", "bitscore", "staxids", "sscinames",
                                          "a1", "a2", "a3", "a4", "a5", "a6"))

qseqids <- forcats::fct_inorder(unique(plant_annotations$qseqid))

plant_annotations$seqtype <- paste0(plant_annotations$a1, plant_annotations$a2, plant_annotations$a3,
                                    plant_annotations$a4, plant_annotations$a5, plant_annotations$a6)
plant_annotations[,c("a1","a2","a3","a4","a5","a6")] <- NULL

plant_annotations <- organize_blast(plant_annotations)

plant_annotations <- plant_annotations %>%
  tidyr::separate(., "sscinames", into=c("genus", "species"),
                  sep=" ", remove=FALSE)

# Clean species names that were not clearly at the beginning of their cell.
plant_annotations$genus[plant_annotations$sscinames=="PREDICTED: Populus trichocarpa"] <- "Populus"
plant_annotations$species[plant_annotations$sscinames=="PREDICTED: Populus trichocarpa"] <- "trichocarpa"

plant_annotations$genus[plant_annotations$sscinames=="PREDICTED: Setaria italica uncharacterized LOC105914727 (LOC105914727)"] <- "Setaria"
plant_annotations$species[plant_annotations$sscinames=="PREDICTED: Setaria italica uncharacterized LOC105914727 (LOC105914727)"] <- "italica"

plant_annotations <- subset(plant_annotations, genus !="Embryophyte")

# Use 'taxize' package to obtain the kingdom, phylum, class, order, and family names for all the taxa.
genera <- unique(plant_annotations$genus)
genera <- genera[!(genera %in% c("Uncultured", "uncultured", "Unclassified"))]

families <- taxize::tax_name(genera, get=c("kingdom", "phylum", "class", "order", "family"),
                             db="ncbi", message=FALSE)

# Proof these annotations to ensure their accuracy.
plant_annotations$genus[plant_annotations$genus=="X"] <- "Elyhordeum"
plant_annotations$species[plant_annotations$species=="Elyhordeum"] <- "langei"
families <- subset(families, query !="X")
families <- rbind(families, c("ncbi", "Elyhordeum", "Viridiplantae", "Streptophyta", "Magnoliopsida", "Poales", "Poaceae"))

families$kingdom[families$query=="Bacillus"] <- "Bacteria"

families$db <- NULL
colnames(families)[1] <- "genus"


# Add full taxonomic information to the annotations file.
plant_annotations <- merge(plant_annotations, families,
                           by="genus", all=TRUE)

# Remove fungi and bacteria.
plant_annotations <- subset(plant_annotations, !(genus %in% c("uncultured", "Uncultured")))

plant_annotations <- subset(plant_annotations, !genus %in% c("Fritillaria", "Aspergillus", "Conyza"))

# Split into a separate table for every ASV.
plant_tbl <- plant_annotations %>%
  tibble::as_tibble(.) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::group_split(.keep=TRUE)

# Create a taxonomy table for the plant ASVs based on their BLAST results.
plant_list <- list()

plant_tax_table <- data.frame(matrix(ncol=8))
colnames(plant_tax_table) <- c("qseqid", "kingdom", "phylum", "class", "order", "family", "genus", "species")

for(i in 1:length(plant_tbl)){
  plant_list[[i]] <- as.data.frame(plant_tbl[[i]])
  
  plant_list[[i]] <- subset(plant_list[[i]],
                            bitscore >= max(plant_list[[i]]$bitscore) & pident > 70)
  
  plant_tax_table[i,1] <- as.character(unique(plant_list[[i]]$qseqid))
  
  if(length(unique(plant_list[[i]]$kingdom))==1){
    plant_tax_table[i,2] <- unique(plant_list[[i]]$kingdom)}
  
  if(length(unique(plant_list[[i]]$phylum))==1){
    plant_tax_table[i,3] <- unique(plant_list[[i]]$phylum)}
  
  if(length(unique(plant_list[[i]]$class))==1){
    plant_tax_table[i,4] <- unique(plant_list[[i]]$class)}
  
  if(length(unique(plant_list[[i]]$order))==1){
    plant_tax_table[i,5] <- unique(plant_list[[i]]$order)}
  
  if(length(unique(plant_list[[i]]$family))==1){
    plant_tax_table[i,6] <- unique(plant_list[[i]]$family)}
  
  if(length(unique(plant_list[[i]]$genus))==1){
    plant_tax_table[i,7] <- unique(plant_list[[i]]$genus)}
  
  if(length(unique(plant_list[[i]]$species))==1){
    plant_tax_table[i,8] <- unique(plant_list[[i]]$species)}
}

# Proof these assignments to make sure they're accurate!
temp <- subset(plant_tax_table, is.na(family)=="TRUE")
temp2 <- subset(plant_annotations, qseqid %in% temp$qseqid)

plant_tax_table[,c(2:8)][plant_tax_table$qseqid=="seq549", ] <- c("Viridiplantae", "Streptophyta", "Magnoliopsida", 
                                                                  "Ericales", "Ericaceae", "Arctostaphylos", "uva-ursi")
plant_tax_table[,c(2:8)][plant_tax_table$qseqid=="seq828", ]
plant_tax_table[,c(2:8)][plant_tax_table$qseqid=="seq857", ]
plant_tax_table[,c(2:8)][plant_tax_table$qseqid=="seq892", ]

rm(temp, temp2)

plant.BLAST.results <- plant_annotations %>%
  dplyr::relocate(genus, .after="family")
plant.BLAST.results$qseqid <- factor(plant.BLAST.results$qseqid,
                                     levels=qseqids)
plant.BLAST.results <- plant.BLAST.results[order(plant.BLAST.results$qseqid), ]

plant_tax_table$qseqid <- factor(plant_tax_table$qseqid, levels=qseqids)
plant_tax_table <- plant_tax_table[order(plant_tax_table$qseqid), ]

rm(plant_annotations, plant_list, plant_tbl, qseqids)

rownames(plant_tax_table) <- plant_tax_table$qseqid
plant_tax_table$qseqid <- NULL

# c - Merge the BLAST results with the original sequence names and create a phyloseq object. ####
plant_tax_table$ids <- rownames(plant_tax_table)

plant_tax_table <- subset(plant_tax_table, ids !="seq435")
plant_tax_table <- rbind(plant_tax_table,
                         c("Viridiplantae", "Streptophyta", "Magnoliopsida", "Caryophyllales",
                           "Chenopodiaceae", "Chenopodium", NA, "seq435"))

plant_tax_table[,8] <- NULL
colnames(plant_tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
plant_tax_table <- as.matrix(plant_tax_table)


# Import initial sample data
sample_data_temp <- subset(sample_data, rownames(sample_data) %in% rownames(seqtabNoC))

# Import into phyloseq object and clean workspace
pseq.PITS.raw <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE),
                          sample_data(sample_data_temp),
                          tax_table(plant_tax_table))

pseq.PITS.raw <- subset_taxa(pseq.PITS.raw, Kingdom=="Viridiplantae")

# d - Remove algae ####
pseq.PITS.filter <- pseq.PITS.raw

# Rename sequences to temp1, temp2, temp3... and count number of ASVs unclassified at kingdom level
pseq.PITS.filter <- rename_temp(pseq.PITS.filter) # 0

# Remove chlorophyta from plants #
temp <- subset_taxa(pseq.PITS.filter, Phylum=="Chlorophyta")
sum(sample_sums(temp)) # 30 reads
sum(sample_sums(temp)) / sum(sample_sums(pseq.PITS.filter)) # with <0.001% total abundance

pseq.PITS.filter <- subset_taxa(pseq.PITS.filter, Phylum !="Chlorophyta")
rm(temp)


###### EXTRACTION BLANKS ####
## [1] Identify sequences that were present in the negative controls and their abundance in other samples ####
pseq.raw <- list(pseq.16S.raw, pseq.18S.raw,
                      pseq.FITS.raw, pseq.PITS.raw)

pseq.filter.list <- list(pseq.16S.filter, pseq.18S.filter,
                         pseq.FITS.filter, pseq.PITS.filter)

rm(fnFs, fnRs, filtFs, filtRs, fastaRef, derepFs, derepRs,
   dadaFs, dadaRs, errF, errR, seqtabAll, seqtabNoC,
   pseq.16S.raw, pseq.18S.raw, pseq.FITS.raw, pseq.PITS.rawr)

ctrl_data <- list()
ctrl_data$extraction.blanks <- list()

for (i in c(1:2)){ # The negative controls were removed from the FITS and PITS data because they had no reads.
  # Create a data frame showing the taxa present in the negative controls, as well as their abundances.
  temp <- transform_sample_counts(pseq.filter.list[[i]], function(x) 100*x/sum(x))
  
  pseq.neg.ctrl <- subset_samples(temp, sample_names(pseq.filter.list[[i]]) %in% c("K1-NC", "NC-K2"))
  pseq.neg.ctrl <- subset_taxa(pseq.neg.ctrl, Class !="Mammalia")
  pseq.neg.ctrl <- prune_taxa(taxa_sums(pseq.neg.ctrl) > 0, pseq.neg.ctrl)
  
  # Save a data frame with the taxonomy information for the negative controls
  ctrl_data$extraction.blanks[[i]] <- list()
  ctrl_data$extraction.blanks[[i]]$tax_table <- as.data.frame(cbind(tax_table(pseq.neg.ctrl)))
  
  bad_pseq <- subset_taxa(temp, taxa_names(temp) %in% rownames(ctrl_data$extraction.blanks[[i]][[1]]))
  
  otu <- as.data.frame(otu_table(bad_pseq))
  otu <- merge(sample_data, otu, by=0, all=TRUE)
  
  ctrl_data$extraction.blanks[[i]]$otu_table <- otu
  
  rm(temp, pseq.neg.ctrl, bad_pseq, otu)
}

rm(pseq.filter.list)

## [2] Remove these sequences from the data, and remove negative controls from the analysis objects ####
# Remove taxa that were present in negative controls from the 16S and 18S data sets.
pseq.16S.filter <- subset_taxa(pseq.16S.filter, !(taxa_names(pseq.16S.filter) %in% ctrls.neg[[1]][[2]]$Row.names))
pseq.18S.filter <- subset_taxa(pseq.18S.filter, !(taxa_names(pseq.18S.filter) %in% ctrls.neg[[2]][[2]]$Row.names))

# Remove negative controls from sample data and phyloseq objects
sample_data <- subset(sample_data, !(rownames(sample_data) %in% c("K1-NC", "NC-K2")))
pseq.16S.filter <- subset_samples(pseq.16S.filter, !(sample_names(pseq.16S.filter) %in% c("K1-NC", "NC-K2")))
pseq.18S.filter <- subset_samples(pseq.18S.filter, !(sample_names(pseq.18S.filter) %in% c("K1-NC", "NC-K2")))

###### MOCK COMMUNITY ANALYSIS ####
## [1] Separate the mock community from the analysis object ####
# Subset to mock community and only the taxa present in the mock community.
ctrl_data$mock.commmunity <- list()
ctrl_data$mock.community$pseq.objects <- list(
  X16S = subset_samples(pseq.16S.filter, sample_names(pseq.16S.filter)=="MOCK"),
  # 18S not included because the mock community is not designed for these primers.
  FITS = subset_samples(pseq.FITS.filter, sample_names(pseq.FITS.filter)=="MOCK")
  # PITS not included because its mock community has no reads.
)

for(i in c(1:3)){
  ctrl_data$mock.community$pseq.objects[[i]] <- prune_taxa(taxa_sums(ctrl_data$mock.community$pseq.objects[[i]]) > 0, 
                                                           ctrl_data$mock.community$pseq.objects[[i]])
}

taxa_names(ctrl_data$mock.community$pseq.objects[[1]]) <-
  paste0("MC-16S-ASV", seq(1:ntaxa(ctrl_data$mock.community$pseq.objects[[1]])))
taxa_names(ctrl_data$mock.community$pseq.objects[[2]]) <-
  paste0("MC-18S-ASV", seq(1:ntaxa(ctrl_data$mock.community$pseq.objects[[2]])))
taxa_names(ctrl_data$mock.community$pseq.objects[[3]]) <-
  paste0("MC-FITS-ASV", seq(1:ntaxa(ctrl_data$mock.community$pseq.objects[[3]])))

# Extract OTU table and sequences
ctrl_data$mock.community$summary.data <- list()

for(i in c(1:3)){
  temp.pseq <- ctrl_data$mock.community$pseq.objects[[i]]
  temp.df <- 
    as.data.frame(t(otu_table(transform_sample_counts(temp.pseq,
                                                      function(x) 100*x/sum(x)))))
  temp.df <- cbind(as.data.frame(cbind(tax_table(temp.pseq))),
                   temp.df)
  
  temp.df <- temp.df[order(-temp.df$MOCK), ]
  
  levels(temp.df$Genus)[levels(temp.df$Genus)=="Escherichia/Shigella"] <- "Escherichia"
  
  ctrl_data$mock.community$summary.data[[i]] <- temp.df
  rm(temp.pseq, temp.df)
}

names(ctrl_data$mock.community$summary.data) <- c("X16S", "X18S", "FITS")

# Remove mock community from the main analysis objects
pseq.16S.filter <- subset_samples(pseq.16S.filter, sample_names(pseq.16S.filter) !="MOCK")
pseq.16S.filter <- prune_taxa(taxa_sums(pseq.16S.filter) > 0, pseq.16S.filter)

pseq.18S.filter <- subset_samples(pseq.18S.filter, sample_names(pseq.18S.filter) !="MOCK")
pseq.18S.filter <- prune_taxa(taxa_sums(pseq.18S.filter) > 0, pseq.18S.filter)

pseq.FITS.filter <- subset_samples(pseq.FITS.filter, sample_names(pseq.FITS.filter) !="MOCK")
pseq.FITS.filter <- prune_taxa(taxa_sums(pseq.FITS.filter) > 0, pseq.FITS.filter)

pseq.PITS.filter <- subset_samples(pseq.PITS.filter, sample_names(pseq.PITS.filter) !="MOCK")
pseq.PITS.filter <- prune_taxa(taxa_sums(pseq.PITS.filter) > 0, pseq.PITS.filter)

sample_data <- subset(sample_data, rownames(sample_data) !="MOCK")

# Keep track of who was lost (i.e., what ASVs were *only* present in the mock community?)
# Based on whether the mock community ASV is still present in the final sample...
ctrl_data$mock.community$summary.data[[1]]$unique.to.mock <- !rownames(ctrls.mock.comm.summary[[1]]) %in% taxa_names(pseq.16S.filter)
ctrl_data$mock.community$summary.data[[2]]$unique.to.mock <- !rownames(ctrls.mock.comm.summary[[2]]) %in% taxa_names(pseq.18S.filter)
ctrl_data$mock.community$summary.data[[3]]$unique.to.mock <- !rownames(ctrls.mock.comm.summary[[3]]) %in% taxa_names(pseq.FITS.filter)

## [2] Calculate Zymo's proprietary `MIQ score` for the 16S mock community ####
# Create a data frame with actual (expected) values for mock community abundances.
ctrl_data$mock.community$actual.values <- data.frame(
  Genus = c("Pseudomonas", "Escherichia/Shigella", "Salmonella", "Lactobacillus", "Enterococcus", "Staphylococcus", "Listeria", "Bacillus"),
  X16S = c(4.2, 10.1, 10.4, 18.4, 9.9, 15.5, 14.1, 17.4)
)

# Create a data frame showing the taxonomy table and relative abundances of the mock community
temp <- tax_glom(ctrl_data$mock.community$pseq.objects$X16S, taxrank="Genus")
temp <- transform_sample_counts(temp, function(x) 100* x/sum(x))
temp.tax <- as.data.frame(cbind(tax_table(temp)))
temp.otu <- as.data.frame(otu_table(temp))
temp.tax$Abundance <- as.vector(t(temp.otu))

# Simple rename to fix taxon names.
temp.tax$Genus <- as.character(temp.tax$Genus)
temp.tax$Genus[temp.tax$Genus=="Limosilactobacillus"] <- "Lactobacillus"

# Merge actual and expected values
temp.tax <- merge(temp.tax, ctrl_data$mock.community$actual.values, by="Genus", all.x=TRUE)

# Calculate relative difference above/below expected values
# (100 = exact match, >100 = MC higher than expected, <100 = MC lower than expected)
temp.tax$scale <- 100 * (temp.tax$Abundance / temp.tax$X16S)
temp.tax <- subset(temp.tax, scale !="NA")

# Define acceptable limits (+/- 15% of the expected percentage)
temp.tax$low.limit <- 85
temp.tax$high.limit <- 115

# Calculate the difference between expected value ('scale') and closest acceptable limit
temp.tax$diff <- rep(0)
for(j in c(1:nrow(temp.tax))){
  if(temp.tax[j,"scale"] > temp.tax[j,"high.limit"]){
    temp.tax[j,"diff"] <- temp.tax[j,"scale"] - temp.tax[j,"high.limit"]
  }
  if(temp.tax[j,"scale"] < temp.tax[j,"low.limit"]){
    temp.tax[j,"diff"] <- temp.tax[j,"low.limit"] - temp.tax[j,"scale"]
  }
}

# Take the square of the difference
temp.tax$diffsq <- temp.tax$diff * temp.tax$diff

# Calculate RMSE
100 - sqrt(sum(temp.tax$diffsq)/nrow(temp.tax)) # 76.93958

rm(temp.tax, temp.otu, temp)

## [3] Generate a phylogenetic tree for the 16S mock community ####
temp.df <- data.frame(
  taxid = taxa_names(ctrl_data$mock.community$pseq.objects$X16S),
  seqs = refseq(ctrl_data$mock.community$pseq.objects$X16S)
)

temp.df <- rbind(temp.df,
                 data.frame(
                   taxid = taxa_names(ctrl_data$mock.community$pseq.objects$X18S),
                   seqs = refseq(ctrl_data$mock.community$pseq.objects$X18S)))

temp.df <- rbind(temp.df,
                 data.frame(
                   taxid = taxa_names(ctrl_data$mock.community$pseq.objects$FITS),
                   seqs = refseq(ctrl_data$mock.community$pseq.objects$FITS)))

fasta <- ShortRead(sread = DNAStringSet(temp.df$seqs), id = BStringSet(temp.df$taxid))
writeFasta(fasta, file = "mock_comm_16S_pool.fasta")

# cat *.fasta > for_alignment.fasta
# mafft --treeout for_alignment.fasta > aligned.fasta

tree <- ggtree::read.tree("for_alignment.fasta.tree")

mock.tree <- ggtree(tree) + 
  geom_tiplab(size=3) + 
  xlim(0,0.9) + 
  labs(tag="a") + 
  plot_theme + 
  theme(panel.border=element_blank()) #+ geom_text(aes(label=node))

rm(tree, fasta, temp.df)

###### CLEAN REMAINING ANALYSIS OBJECTS ####
## [1] Examine read count distribution and sample completeness ####
# Create a single data frame showing the number of reads for every sample/amplicon pair.
read_counts.final <- transform(merge(as.data.frame(sample_sums(pseq.16S.filter)), as.data.frame(sample_sums(pseq.18S.filter)),
                                     by=0,all=TRUE),
                               row.names=Row.names, Row.names=NULL)
read_counts.final <- transform(merge(read_counts.final, as.data.frame(sample_sums(pseq.FITS.filter)),
                                     by=0,all=TRUE),
                               row.names=Row.names, Row.names=NULL)
read_counts.final <- transform(merge(read_counts.final, as.data.frame(sample_sums(pseq.PITS.filter)),
                                     by=0,all=TRUE),
                               row.names=Row.names, Row.names=NULL)
colnames(read_counts.final) <- c("X16S", "X18S", "FITS", "PITS")

# Calculate sample completeness for every sample/amplicon pair
iNext.output <- list(
  X16S = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.16S.filter))), q=0, datatype="abundance"),
  X18S = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.18S.filter))), q=0, datatype="abundance"),
  FITS = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.FITS.filter))), q=0, datatype="abundance"),
  PITS = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.PITS.filter))), q=0, datatype="abundance")
)

# Rename to "observed" (richness) and "SC" (sample completeness)
for(i in c(1:4)){
  rownames(iNext.output[[i]]$DataInfo) <- iNext.output[[i]]$DataInfo$site
  if(i==1){colnames(iNext.output[[i]]$DataInfo)[c(3:4)] <- c("Obs.16S", "SC.16S")}
  if(i==2){colnames(iNext.output[[i]]$DataInfo)[c(3:4)] <- c("Obs.18S", "SC.18S")}
  if(i==3){colnames(iNext.output[[i]]$DataInfo)[c(3:4)] <- c("Obs.FITS", "SC.FITS")}
  if(i==4){colnames(iNext.output[[i]]$DataInfo)[c(3:4)] <- c("Obs.PITS", "SC.PITS")}
}

# Add species richness and sample completeness information to the read counts object
for(i in c(1:4)){
  read_counts.final <- transform(merge(read_counts.final, iNext.output[[i]]$DataInfo[,c(3:4)],
                                       by=0,all=TRUE),
                                 row.names=Row.names, Row.names=NULL)
}

# Visualize read count distributions
read_counts.melt <- reshape2::melt(read_counts.final)

ggplot(read_counts.melt, aes(x=value)) +
  geom_histogram() +
  facet_wrap(~variable, scales="free_x")

rm(read_counts.melt)

## [2] Remove samples that produced fewer than 2000 reads. ####
pseq.16S.filter <- subset_samples(pseq.16S.filter, sample_sums(pseq.16S.filter) > 2000) # Removes 1 low-read sample
pseq.16S.filter <- prune_taxa(taxa_sums(pseq.16S.filter) > 0, pseq.16S.filter)

# Save the PITS low-read samples just for interest.
PITS.low.read.samples <- subset_samples(pseq.PITS.filter, sample_sums(pseq.PITS.filter) < 200)
PITS.low.read.samples <- prune_taxa(taxa_sums(PITS.low.read.samples) > 0, PITS.low.read.samples)

PITS.low.read.samples <- merge(as.data.frame(t(otu_table(PITS.low.read.samples))),
                               as.data.frame(cbind(tax_table(PITS.low.read.samples))),
                               by=0, all=TRUE)

pseq.PITS.filter <- subset_samples(pseq.PITS.filter, sample_sums(pseq.PITS.filter) > 10000) # Removes 4 low-read samples
pseq.PITS.filter <- prune_taxa(taxa_sums(pseq.PITS.filter) > 0, pseq.PITS.filter)
PITS.low.read.samples$unique <- !PITS.low.read.samples$Row.names %in% taxa_names(pseq.PITS.filter)

## [3] Remove low-abundance taxa (including singletons) ####
# Calculate the number of low-abundance taxa
ntaxa(pseq.16S.filter) # 11679 total
ntaxa(prune_taxa(taxa_sums(pseq.16S.filter)==1, pseq.16S.filter)) # 2952 singletons
x1 <- transform_sample_counts(pseq.16S.filter, function(x) 100*x/sum(x))
ntaxa(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE)) # 6362 taxa with mean relative abundance < 0.001
mean(sample_sums(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE))) # averaging 1.515% of each sample
sd(sample_sums(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE))) # +/- 1.306% of each sample

ntaxa(pseq.18S.filter) # 3534 total
ntaxa(prune_taxa(taxa_sums(pseq.18S.filter)==1, pseq.18S.filter)) # 1208 singletons
x1 <- transform_sample_counts(pseq.18S.filter, function(x) 100*x/sum(x))
ntaxa(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE)) # 2159 taxa with mean relative abundance < 0.001
mean(sample_sums(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE))) # averaging 0.442% of each sample
sd(sample_sums(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE))) # +/- 0.380% of each sample

ntaxa(pseq.FITS.filter) # 4635 total
ntaxa(prune_taxa(taxa_sums(pseq.FITS.filter)==1, pseq.FITS.filter)) # 36 singletons
x1 <- transform_sample_counts(pseq.FITS.filter, function(x) 100*x/sum(x))
ntaxa(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE)) # 1573 taxa with mean relative abundance < 0.001
mean(sample_sums(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE))) # averaging 0.803% of each sample
sd(sample_sums(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE))) # +/- 0.704% of each sample

ntaxa(pseq.PITS.filter) # 1228 total
ntaxa(prune_taxa(taxa_sums(pseq.PITS.filter)==1, pseq.PITS.filter)) # 199 singletons
x1 <- transform_sample_counts(pseq.PITS.filter, function(x) 100*x/sum(x))
ntaxa(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE)) # 359 taxa with mean relative abundance < 0.001
mean(sample_sums(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE))) # averaging 0.405% of each sample
sd(sample_sums(filter_taxa(x1, function(x) mean(x) < 0.001, TRUE))) # +/- 0.0737% of each sample

rm(x1)

# Remove low-abundance taxa from analysis
pseq.temp <- list(pseq.16S.filter, pseq.18S.filter, pseq.FITS.filter, pseq.PITS.filter)

pseq.post.trim.list <- list()

for(i in c(1:4)){
  # Transform by relative abundance
  temp <- transform_sample_counts(pseq.temp[[i]], function(x) 100*x/sum(x))
  
  # Save taxa with >0.001% abundance
  temp <- filter_taxa(temp, function(x) mean(x) > 0.001, TRUE)
  
  # Subset these taxa with absolute read counts
  pseq.post.trim.list[[i]] <- subset_taxa(pseq.temp[[i]], taxa_names(pseq.temp[[i]]) %in% taxa_names(temp))
  
  rm(temp)
}

rm(pseq.temp)

## [4] Rarefaction ####
taxa_names(pseq.post.trim.list[[1]]) <- paste0("ASV-16S-", seq(1:ntaxa(pseq.post.trim.list[[1]])))
taxa_names(pseq.post.trim.list[[2]]) <- paste0("ASV-18S-", seq(1:ntaxa(pseq.post.trim.list[[2]])))
taxa_names(pseq.post.trim.list[[3]]) <- paste0("ASV-FITS-", seq(1:ntaxa(pseq.post.trim.list[[3]])))
taxa_names(pseq.post.trim.list[[4]]) <- paste0("ASV-PITS-", seq(1:ntaxa(pseq.post.trim.list[[4]])))

pseq.clr <- pseq.post.trim.list
pseq.rarefied <- pseq.post.trim.list
rm(pseq.post.trim.list)

for(p in c(1:4)){
  rarefaction.average <- list()
  
  for(i in 1:1000){
    temp.rarefy <- rarefy_even_depth(pseq.rarefied[[p]], 
                                     sample.size = min(sample_sums(pseq.rarefied[[p]])),
                                     replace = FALSE,
                                     trimOTUs = FALSE,
                                     rngseed=i)
    rarefaction.average[[i]] <- as.data.frame(otu_table(temp.rarefy))
  }
  
  dfAvg <- Reduce("+", rarefaction.average)/length(rarefaction.average)
  dfAvg <- round(dfAvg, 0)
  dfAvg <- otu_table(dfAvg, taxa_are_rows=FALSE)
  
  # Replace feature table in phyloseq object with new, rarefied feature table.
  otu_table(pseq.rarefied[[p]]) <- dfAvg
  pseq.rarefied[[p]] <- prune_taxa(taxa_sums(pseq.rarefied[[p]]) > 0, pseq.rarefied[[p]])
  rm(temp.rarefy, rarefaction.average, dfAvg)
}

## [5] Construct a phylogenetic tree for both filtered and rarefied data #####
library(DECIPHER)
library(phangorn)

pseq.temp.list <- pseq.clr # Repeat with pseq.rarefied... just takes too much memory if I do them all at once.

for(i in c(1:4)){
  vector <- taxa_names(pseq.temp.list[[i]])
  taxa_names(pseq.temp.list[[i]]) <- refseq(pseq.temp.list[[i]])
  seqtabNoC <- as.matrix(as.data.frame(otu_table(pseq.temp.list[[i]])))
  seqs <- getSequences(seqtabNoC)
  
  names(seqs) <- seqs
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, verbose=FALSE)
  phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phangAlign)
  treeNJ <- NJ(dm)
  fit = pml(treeNJ, data=phangAlign)
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement="stochastic", control=pml.control(trace=0))
  
  phy_tree(pseq.temp.list[[i]]) <- fitGTR$tree
  taxa_names(pseq.temp.list[[i]]) <- vector
  
  rm(fitGTR, fit, treeNJ, dm, phangAlign, alignment, seqs, seqtabNoC, vector)
}

pseq.clr <- pseq.temp.list # Replace with pseq.rarefied afterwards
rm(pseq.temp.list)

# Check to make sure the phylogenetic trees are reasonably accurate.
subset <- prune_taxa(taxa_names(pseq.clr[[2]])[1:50], pseq.clr[[2]])
head(tax_table(subset))

taxa_names(subset) <- paste0(tax_table(subset)[,6], seq(c(1:50)))

plot_tree(subset, nodelabf=nodeplotblank, "treeonly", label.tips="taxa_names", ladderize="left")

## [6] Add guild information for the fungal taxa ####
devtools::install_github("brendanf/FUNGuildR")

taxtabNoC <- as.data.frame(cbind(tax_table(pseq.clr[[3]])))

taxtabNoC$Taxonomy <- paste0(taxtabNoC$Kingdom, ";",
                             taxtabNoC$Phylum, ";",
                             taxtabNoC$Class, ";",
                             taxtabNoC$Order, ";",
                             taxtabNoC$Family, ";",
                             taxtabNoC$Genus, ";",
                             taxtabNoC$Species)

otu_table <- as.data.frame(cbind(otu_table(pseq.clr[[3]])))
otu_table <- as.data.frame(t(otu_table))
otu_table$Taxonomy <- taxtabNoC$Taxonomy

db <- FUNGuildR::get_funguild_db()

fungal.guilds <- FUNGuildR::funguild_assign(otu_table, db = db, tax_col = "Taxonomy")

rownames(fungal.guilds) <- taxa_names(pseq.clr$FITS)
fungal.guilds <- fungal.guilds[,c((ncol(fungal.guilds)-10):ncol(fungal.guilds))]
fungal.guilds <- cbind(taxtabNoC, fungal.guilds)
 
rm(db, taxtabNoC)

# Change Dothiorella to Lichenoconium, but only for ASVs that BLAST as Lichenoconium.
# I only did this after I had named ASVs, so to ensure the code runs correctly, this code appears here
# rather than earlier, where the rest of the taxonomic cleaning was done.
# Once for rarefied data
tax.tab.fungi <- as.data.frame(cbind(tax_table(pseq.rarefied$FITS)))

temp <- subset(tax.tab.fungi, Genus=="Dothiorella")
temp <- subset(temp, !(rownames(temp) %in% c("ASV-FITS-194", "ASV-FITS-740")))
tax.tab.fungi$Genus[rownames(tax.tab.fungi) %in% rownames(temp)] <- "Lichenoconium"
tax.tab.fungi$Family[tax.tab.fungi$Genus=="Lichenoconium"] <- "Lichenoconiaceae"
tax.tab.fungi$Order[tax.tab.fungi$Genus=="Lichenoconium"] <- "Lichenoconiales"
tax.tab.fungi$Class[tax.tab.fungi$Genus=="Lichenoconium"] <- "Dothideomycetes"
tax.tab.fungi$Phylum[tax.tab.fungi$Genus=="Lichenoconium"] <- "Ascomycota"
tax.tab.fungi$Species[tax.tab.fungi$Genus=="Lichenoconium"] <- NA

tax.tab.fungi <- as.matrix(tax.tab.fungi)
tax_table(pseq.rarefied$FITS) <- tax.tab.fungi

# Repeat for CLR data
tax.tab.fungi <- as.data.frame(cbind(tax_table(pseq.clr$FITS)))

temp <- subset(tax.tab.fungi, Genus=="Dothiorella")
temp <- subset(temp, !(rownames(temp) %in% c("ASV-FITS-194", "ASV-FITS-740")))
tax.tab.fungi$Genus[rownames(tax.tab.fungi) %in% rownames(temp)] <- "Lichenoconium"
tax.tab.fungi$Family[tax.tab.fungi$Genus=="Lichenoconium"] <- "Lichenoconiaceae"
tax.tab.fungi$Order[tax.tab.fungi$Genus=="Lichenoconium"] <- "Lichenoconiales"
tax.tab.fungi$Class[tax.tab.fungi$Genus=="Lichenoconium"] <- "Dothideomycetes"
tax.tab.fungi$Phylum[tax.tab.fungi$Genus=="Lichenoconium"] <- "Ascomycota"
tax.tab.fungi$Species[tax.tab.fungi$Genus=="Lichenoconium"] <- NA

tax.tab.fungi <- as.matrix(tax.tab.fungi)
tax_table(pseq.clr$FITS) <- tax.tab.fungi
