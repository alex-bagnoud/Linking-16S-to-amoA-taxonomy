# Linking-16S-to-amoA-taxonomy

Script written by alexandre.bagnoud@gmail.com in 2018.

## Introduction

This bioinformatic pipeline aims at linking 16S rRNA sequences of [Thaumarchaeota](https://en.wikipedia.org/wiki/Thaumarchaeota) (i.e. ammonia-oxidizing archaea, or AOA) to *amoA* phylogeny. *amoA* is the A subunit of the [ammonia-monooxygenase](https://en.wikipedia.org/wiki/Ammonia_monooxygenase) and it is used as a genetic marker to detect AOA in environments. It's phylogeny, unlike the one of 16S rRNA of Thaumarchaeota, is finely characterised and contains useful information about their distribution in environments, as described in [Alves et *al*., 2018](https://www.nature.com/articles/s41467-018-03861-1). By merging both phylogeny, we hope to transfer this useful information from the *amoA* phylogeny to the 16S rRNA phylogeny.

This work will be soon published in an article entitled *Linking 16S rRNA gene to amoA gene taxonomy reveals environmental distribution of ammonia-oxidizing archaea clades in peatland soils* by Haitao Wang, Alexandre Bagnoud, Rafael I. Ponce-Toledo, Melina Kerou, Micha Weil, Christa Schleper and Tim Urich.

## Approach

The approach used here is to screen genomic databases in order to find genomes, contigs or nucleotides that contain both a 16S rRNA and a *amoA* gene. Here are the main steps of the pipeline:

### Datasets and Database

Three different datasets were used here:

1. All [NCBI](https://www.ncbi.nlm.nih.gov/) sequences annotated as Archaea and that are at least 2000 bp long.
2. All the [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) archaeal genomes.
3. All the [GeneBank](https://www.ncbi.nlm.nih.gov/genbank/) archaeal genomes.

Each dataset was analysed separately.

For this study, dataset 1 was downloaded the 31st Oct. 2020, and datasets 2 and 3 were downloaded on the 7th of Nov. 2020.

The *amoA* database from [Alves et *al*., 2018](https://www.nature.com/articles/s41467-018-03861-1) (accessible from the supplementary information) was used to detect and annotate *amoA* genes from the 3 datasets.

### Finding and isolating *amoA* genes

The 3 nucleotides and genomes datasets were BLAST against Alves et *al*. database. Then, the sequence of positive matches were extracted using Samtools.

### Annotation of *amoA* genes

The newly found *amoA* sequences were then annotated based on Alves et *al*. database, using QIIME1, as recommanded by the authors.

### Extraction of 16S rRNA genes

For each sequence that contains an *amoA* gene, 16S rRNA gene was extracted (if present). 16S rRNA genes were detected with Barrnap and then the sequences were extraction with Samtools.

### Merging 16S rRNA and *amoA* sequences

Then, all 16S rRNA and *amoA* pairs (i.e. coming from the same original sequence) were reunited using an in-house R script.

### Merging the outcome of the 3 analyses

Finally, the results from all 3 runs were merged together using an in-house R script.

## Dependecies

Here is the list of the softwares used for this pipeline. The versions that were used for this study are indicated in parenthesis:

* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (v.2.10.1)
* [Samtools](http://www.htslib.org/) (v.1.8)
* [QIIME](http://qiime.org/) (v.1.9.1)
* [Barrnap](https://github.com/tseemann/barrnap) (v.0.9)
* [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) (v.0.2.6)
* [vsearch](https://github.com/torognes/vsearch) (v2.7.1)

## Database / Datasets download

### Alves *et al*. database

The *amoA* database was downoladed from the supplementary informatio of Alves *et al*. (https://doi.org/10.1038/s41467-018-03861-1). Supplementary files 1 et 3 were downloaded and unzipped. The files `AamoA.db_an96.aln_tax.annotated.fasta` (supplementary data 1), `AamoA.db_nr.aln.fasta`, and `AamoA.db_nr.aln_taxonomy_qiime.txt` (in the `AamoA.db_nr_qiime.mothur/` folder of supplementary data 3) were moved to `0-databases/`.

A BLAST database was then made using BLAST+ using this command line:
```sh
makeblastdb -dbtype nucl -in 0-databases/AamoA.db_an96.aln_tax.annotated.fasta -out 0-databases/AamoA.db_an96.aln_tax.annotated
```
(Unless specified, the code presented here is written in bash.)
### NCBI nucleotides sequences

Sequences were directly downloaded from this website: https://www.ncbi.nlm.nih.gov/nuccore, using as search query `(Archaea[Organism]) AND 2000:99999999999999[Sequence Length]`. The ouput was then dowloaded using `Send to file > FASTA format`. The sequences were downloaded on the 31st Oct. 2020. 

The fasta file was save as `1-ncbi_arch_2000/1-raw_data/ncbi_nucl_arch_min2000.fasta`. This file containes 335,202 sequences.

```sh
grep "^>" 1-ncbi_arch_2000/1-raw_data/ncbi_nucl_arch_min2000.fasta | wc -l
> 335202
```

### RefSeq archaeal genomes

The archaeal RefSeq genomes were downloaded with ncbi-genome-download on the 7th of November 2020, using this commande line:

```sh
mkdir 2-refseq_arch_genomes/
cd 2-refseq_arch_genomes/
ncbi-genome-download archaea --format fasta --assembly-level all --section refseq --output-folder 1-raw_data/arch_genomes_refseq archaea
```

Then, genomes were unarchived using this script:
```sh
mkdir 1-raw_data/arch_genomes_refseq_unarchived

for file in 1-raw_data/arch_genomes_refseq/refseq/archaea/GCF_*/*.fna.gz; do
	echo $file
	id=$(echo $file | cut -d "/" -f5)
	echo $id
	cp $file 1-raw_data/arch_genomes_refseq_unarchived/${id}.fna.gz
	gzip -d 1-raw_data/arch_genomes_refseq_unarchived/${id}.fna.gz
done
```
How many genomes were dowloaded?
```sh
ls -l 1-raw_data/arch_genomes_refseq_unarchived/ | wc -l
> 1106
```

### GeneBank archaeal genomes

GeneBank genomes were downloaded follwing the instructions  of this github repository : https://github.com/rprops/MetaG_analysis_workflow/wiki/09.-Download-genomes-NCBI-EDI. The genomes were downloaded on the 7th of November 2020 using these command lines:

Create folders:
```sh
mkdir 3-genbank_arch_genomes/
cd 3-genbank_arch_genomes/
mkdir 1-raw_data/
```

Download the list of archeal genomes:
```sh
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt
mv assembly_summary.txt 1-raw_data/1-assembly_summary.txt
```
How many genomes are there?
```sh
wc -l 1-raw_data/1-assembly_summary.txt
> 5738
```
Get the ftp links
```sh
less 1-raw_data/1-assembly_summary.txt | cut -f20 > 1-raw_data/2-ftp_links.txt
```
Download them all!
```sh
mkdir 1-raw_data/3-archaeal_genomes

for next in $(cat 1-raw_data/2-ftp_links.txt); do wget -P 1-raw_data/3-archael_genomes "$next"/*genomic.fna.gz; done
```
How many genomes were downloaded?
```sh
ls -lh 1-raw_data/3-archaeal_genomes/ | grep -v cds | grep -v rna |wc -l
> 5737
```

## Detailed script (and files)

Then, for each dataset, the following script was applied. Each command must be executed from the dataset folder `1-ncbi_arch_2000/`, `2-refseq_arch_genomes/`, or `3-genbank_arch_genomes/`.
The script presented below is the one that was specifically used to analyse RefSeq genomes. The two other dataset were analysed with script slightly different. All 3 scripts can be found in this [folder](scripts/).

### Blast amoa genes
```sh
mkdir 2-amoa_blast
```
This script blast each genome in the folder to the archaeal amoA db and keeps the best hit. It adds to the last line of the blase output file the path to the genome.

```sh
for file in 1-raw_data/arch_genomes_refseq_unarchived/*.fna; do
	id=$(echo ${file##*\/} | cut -d "." -f1)
	echo $id
	blastn -query $file -db ../0-databases/AamoA.db_an96.aln_tax.annotated \
	-outfmt 6 -max_target_seqs 1 -out 2-amoa_blast/1-${id}_blast.txt
	awk -v a="$file" -F"," 'BEGIN { OFS = "\t" } {$2=a; print}' 2-amoa_blast/1-${id}_blast.txt > 2-amoa_blast/2-${id}_blast2.txt
done

cat 2-amoa_blast/2-* > 2-amoa_blast/3-cat_blast.txt
rm 2-amoa_blast/1-*
rm 2-amoa_blast/2-*
```
Extract list of genomes that harbor an amoA:
```sh
less 2-amoa_blast/3-cat_blast.txt | cut -f13 | sort -u > 2-amoa_blast/4-amoa_genome_list.txt
```
### Get amoA sequences
```sh
mkdir 3-amoa_seqs
```
Extract amoA sequences with Samtools:
```sh
while read p; do
	echo $p
	start=$(echo $p | cut -d " " -f7)
	end=$(echo $p | cut -d " " -f8)
	seq=$(echo $p | cut -d " " -f1)
	file=$(echo $p | cut -d " " -f13)
	id=$(echo ${file##*\/} | cut -d "." -f1)
	echo "samtools faidx $file ${seq}:${start}-${end} > 2-amoa_seqs/1-${id}_amoa.fasta"
	samtools faidx $file ${seq}:${start}-${end} >> 3-amoa_seqs/1-amoa.fasta
done < 2-amoa_blast/3-cat_blast.txt

rm 1-raw_data/arch_genomes_refseq_unarchived/*.fna.fai
```
Annotate them using Alves *et al*. database and QIIME1:
```sh
db_seq="../0-databases/AamoA.db_nr.aln.fasta"
qiime_tax="../0-databases/AamoA.db_nr.aln_taxonomy_qiime.txt"

source activate qiime1
assign_taxonomy.py -i 3-amoa_seqs/1-amoa.fasta -t $qiime_tax -r $db_seq --similarity 0.9 \
-o 3-amoa_seqs/2-uclust_annotations/
source deactivate qiime1
```

### Extract 16S from amoA genomes
```sh
mkdir 4-16S_genes
```
Find 16S rRNA genes with Barrnap:
```sh
while read p; do
	id=$(echo ${p##*\/} | sed 's/.fna//')
	barrnap $p --kingdom arc| grep 16S > 4-16S_genes/1-${id}_barrnap.txt
done < 2-amoa_blast/4-amoa_genome_list.txt

find 4-16S_genes/ -size 0 -delete
```
Extract 16S rRNA sequences with Samtools:
```sh
for file in 4-16S_genes/1-*; do
	echo $file
	id=$(echo $file | cut -d "-" -f3 | sed 's/_barrnap.txt//')
	echo $id
	while read p; do
		seq=$(echo $p | cut -d " " -f1)
		start=$(echo $p | cut -d " " -f4)
		end=$(echo $p | cut -d " " -f5)
		echo "samtools faidx 0-genomes/${id}.fasta ${seq}:${start}-${end} >> 4-16S_genes/2-${id}_16S.fasta"
		samtools faidx 1-raw_data/arch_genomes_refseq_unarchived/${id}.fna ${seq}:${start}-${end} >> 4-16S_genes/2-${id}_16S.fasta
	done < $file
done 

rm 1-raw_data/arch_genomes_refseq_unarchived/*.fna.fai

cat 4-16S_genes/2* > 4-16S_genes/3-all_16S_seq.fasta
```


### Merging all the data and output annotations files
This exectuable in-house R script merge the sequence and annotation files.
```sh
Rscript ../scripts/script2_amoa_16S_from_genomes.R
```
One hase to run `script2_amoa_16S_from_contigs.R` for merging data from the NCBI nucleotide dataset.

### Concatenating outputs from the 3 datasets

Putting together all 16S from the 3 datasets:
```sh
cat ../*/5-annotation_files/1-16S_db.fasta > 1-16S_3runs_3.fasta
```

Finding unique 16S rRNA genes:
```sh
vsearch --derep_fulllength 1-16S_3runs_3.fasta --output 2-unique_16S.fasta
```

List of unique 16S sequences headers:
```sh
grep "^>" 2-unique_16S.fasta | sed 's/>//' > 3-unique_16S_header_list.txt
```

Putting together all QIIME annotations from the 3 datasets:
```sh
cat ../*/5-annotation_files/2-16S_amoa_tax_qiime.txt > 4-16S_amoa_tax_qiime_3runs.txt
```

Subset the unique annotations, using R:
```R
# Import files
header <- read.table("3-unique_16S_header_list.txt", header = FALSE)
annot <- read.table("4-16S_amoa_tax_qiime_3runs.txt", header = FALSE, sep = "\t")

# Remove duplicated entries in the annotations
annot2 <- annot[!duplicated(annot$V1),]

# merge df
m <- merge(header, annot2)

# Export new file

write.table(m, "5-unique_16S_amoa_tax_qiime_3runs.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```
Finally, concatenate all *amoA* sequences from the 3 dataset:
```sh
cat ../*/3-amoa_seqs/1-amoa.fasta > 6_amoa_seqs_3runs.fasta
```

## Relevant output files

## How to annotate your Thaumarchaeota 16S rRNA sequences with this database?
