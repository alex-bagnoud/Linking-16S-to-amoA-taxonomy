#!/usr/bin/env Rscript


### 2) Merge all data into a single dataframe

## 2.1) Import amoA contigs, amoA genomes, blast score

df.blast <- read.table("2-amoa_blast/1-blast.txt", sep = "\t")
colnames(df.blast) <- c("amoa.contig", "best-hit", "pident", "length", "mismatch", "gapopen",
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore")

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

df.barrnap <- read.table("4-16S_genes/4-barrnap.txt", sep = "\t")
colnames(df.barrnap) <- c("I6S.contig", "V2", "V3", "I6S.start", "I6S.end", "I16S.pval", "I16.dir", "V8", "I6S.product")
df.barrnap$V2 <- NULL
df.barrnap$V3 <- NULL
df.barrnap$V8 <- NULL

# Create a 16S header column
df.barrnap$I6S.header <- paste0(df.barrnap$I6S.contig, ":", df.barrnap$I6S.start, "-", df.barrnap$I6S.end)

## 2.5) Import all 16S sequences

df.16S <- fastaToDf("4-16S_genes/5-16S.fasta")
colnames(df.16S) <- c("I6S.header", "I6S.seq")

# Extract the contig name
df.16S$I6S.contig <- str_split_fixed(df.16S$I6S.header, ":", 2)[,1]

## 2.6) Combine all dataframes

m1 <- merge(df.blast2, df.amoa, by.x = "amoa.contig", by.y = "amoa.contig")
m2 <- merge(m1, df.annot, by.x = "amoa.header", by.y = "amoa.header")
m3 <- merge(df.barrnap, df.16S, by.x = "I6S.header", by.y = "I6S.header")

m4 <- merge(m2, m3, by.x = "amoa.contig.y", by.y = "I6S.contig.y")


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

