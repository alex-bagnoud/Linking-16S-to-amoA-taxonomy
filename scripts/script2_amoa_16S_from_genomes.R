#!/usr/bin/env Rscript

### 1) Concatenate barrnap outputs while adding the path to the corresponding genome file at the end of each entry

## 1.1) List all barrnap files
files <- list.files(path="4-16S_genes/", pattern = "1-", full.names=T, recursive=FALSE)

## 1.2) Import their content in a dataframe together with the path of the corresponding genomes
vec <- character()
df <- NULL
for (f in files) {
    print(f)
    x <- read.table(f, sep = "\t")
    df <- rbind(df, x)
    y <- rep(f , nrow(x))
    vec <- c(vec, y)
}

df
vec2 <- sub("_barrnap.txt", "", sub(pattern = "4-16S_genes//1-", "", vec))
vec3 <- paste0("1-raw_data/arch_genomes_refseq_unarchived/", vec2, ".fna")
df <- cbind(df, vec3)

## 1.3) Export the file
write.table(df, "4-16S_genes/4-barrnap_list.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


### 2) Merge all data into a single dataframe

## 2.1) Import amoA contigs, amoA genomes, blast score

df.blast <- read.table("2-amoa_blast/3-cat_blast.txt", sep = "\t")
colnames(df.blast) <- c("amoa.contig", "best-hit", "pident", "length", "mismatch", "gapopen",
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore", "genome.path")

# remove low-quality hits
hist(df.blast$bitscore, breaks = 20)

df.blast2 <- df.blast[df.blast$bitscore > 900,]

## 2.2) Import amoA sequences 

library("Biostrings")

fastaToDf <- function(fastaFile){
    dnaSeq <- readBStringSet(fastaFile)
    fasta_df <- data.frame(header = names(dnaSeq), sequence = paste(dnaSeq))
}

df.amoa <- fastaToDf("3-amoa_seqs/1-amoa.fasta")
colnames(df.amoa) <- c("amoa.header", "amoa.seq")

# Extract the contig name
library("stringr")
df.amoa$amoa.contig <- str_split_fixed(df.amoa$amoa.header, ":", 2)[,1]

## 2.3) Import amoA annotation

df.annot <- read.table("3-amoa_seqs/2-uclust_annotations/1-amoa_tax_assignments.txt", sep = "\t")
colnames(df.annot) <- c("amoa.header", "amoa.annot", "V3", "V4")
df.annot$V3 <- NULL
df.annot$V4 <- NULL

# Extract the contig name
df.annot$amoa.contig <- str_split_fixed(df.annot$amoa.header, ":", 2)[,1]

## 2.4) Import barrnap list

df.barrnap <- read.table("4-16S_genes/4-barrnap_list.txt", sep = "\t")
colnames(df.barrnap) <- c("I6S.contig", "V2", "V3", "I6S.start", "I6S.end", "I16S.pval", "I16.dir", "V8", "I6S.product", "genome.path")
df.barrnap$V2 <- NULL
df.barrnap$V3 <- NULL
df.barrnap$V8 <- NULL

# Create a 16S header column
df.barrnap$I6S.header <- paste0(df.barrnap$I6S.contig, ":", df.barrnap$I6S.start, "-", df.barrnap$I6S.end)

## 2.5) Import all 16S sequences

df.16S <- fastaToDf("4-16S_genes/3-all_16S_seq.fasta")
colnames(df.16S) <- c("I6S.header", "I6S.seq")

# Extract the contig name
df.16S$I6S.contig <- str_split_fixed(df.16S$I6S.header, ":", 2)[,1]

## 2.6) Combine all dataframes

m1 <- merge(df.blast2, df.amoa, by.x = "amoa.contig", by.y = "amoa.contig")
m2 <- merge(m1, df.annot, by.x = "amoa.header", by.y = "amoa.header")
m3 <- merge(df.barrnap, df.16S, by.x = "I6S.header", by.y = "I6S.header")

m2$genome.path <- as.character(m2$genome.path)
m3$genome.path <- as.character(m3$genome.path)


m4 <- merge(m2, m3, by.x = "genome.path", by.y = "genome.path")


### 3) Export annotation files

dir.create("5-annotation_files")

## 3.1) 16S db fasta file

dfToFasta <- function(header, seq, file){
    sequence = BStringSet(seq)
    names(sequence) <- header
    writeXStringSet(sequence, file)
}

dfToFasta(m4$I6S.header, m4$I6S.seq, "5-annotation_files/1-16S_db.fasta")

## 3.2) Taxonomy file for QIIME

df.16Stax <- cbind(m4$I6S.header, as.character(m4$amoa.annot))
write.table(df.16Stax, "5-annotation_files/2-16S_amoa_tax_qiime.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(m4, "5-annotation_files/3-full_pairs_table.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

