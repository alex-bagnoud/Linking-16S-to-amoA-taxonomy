# Linking-16S-to-amoA-taxonomy

Script written by alexandre.bagnoud@gmail.com in 2018.

### Introduction

This bioinformatic pipeline aims at linking 16S rRNA sequences of [Thaumarchaeota](https://en.wikipedia.org/wiki/Thaumarchaeota) (i.e. ammonia-oxidizing archaea, or AOA) to *amoA* phylogeny. *amoA* is the A subunit of the [ammonia-monooxygenase](https://en.wikipedia.org/wiki/Ammonia_monooxygenase) and it is used as a genetic marker to detect AOA in environments. It's phylogeny, unlike the one of 16S rRNA of Thaumarchaeota, is finely characterised and contains useful information about their distribution in environments, as described in [Alves et *al*., 2018](https://www.nature.com/articles/s41467-018-03861-1). By merging both phylogeny, we hope to transfer this useful information from the *amoA* phylogeny to the 16S rRNA phylogeny.

This work will be soon published in an article entitled *Linking 16S rRNA gene to amoA gene taxonomy reveals environmental distribution of ammonia-oxidizing archaea clades in peatland soils* by Haitao Wang, Alexandre Bagnoud, Rafael I. Ponce-Toledo, Melina Kerou, Micha Weil, Christa Schleper and Tim Urich.

### Approach

The approach used here is to screen genomic databases in order to find genomes, contigs or nucleotides that contain both a 16S rRNA and a *amoA* gene. Here are the main steps of the pipeline:

#### Datasets and Database

Three different datasets were used here:

1. All [NCBI](https://www.ncbi.nlm.nih.gov/) sequences annotated as Archaea and that are at least 2000 bp long.
2. All the [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) archaeal genomes.
3. All the [GeneBank](https://www.ncbi.nlm.nih.gov/genbank/) archaeal genomes.

Each dataset was analysed separately.

For this study, dataset 1 was downloaded the 31st Oct. 2020, and datasets 2 and 3 were downloaded on the 7th of Nov. 2020.

The *amoA* database from [Alves et *al*., 2018](https://www.nature.com/articles/s41467-018-03861-1) (accessible from the supplementary information) was used to detect and annotate *amoA* genes from the 3 datasets.

#### Finding and isolating *amoA* genes

The 3 nucleotides and genomes datasets were BLAST against Alves et *al*. database. Then, the sequence of positive matches were extracted using Samtools.

#### Annotation of *amoA* genes

The newly found *amoA* sequences were then annotated based on Alves et *al*. database, using QIIME1, as recommanded by the authors.

#### Extraction of 16S rRNA genes

For each sequence that contains an *amoA* gene, 16S rRNA gene was extracted (if present). 16S rRNA genes were detected with Barrnap and then the sequences were extraction with Samtools.

#### Merging 16S rRNA and *amoA* sequences

Then, all 16S rRNA and *amoA* pairs (i.e. coming from the same original sequence) were reunited using an in-house R script.

#### Merging the outcome of the 3 analyses

Finally, the results from all 3 runs were merged together using an in-house R script.

### Dependecies

Here is the list of the softwares used for this pipeline. The versions that were used for this study are indicated in parenthesis:

* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (v.2.10.1)
* [Samtools](http://www.htslib.org/) (v.1.8)
* [QIIME](http://qiime.org/) (v.1.9.1)
* [Barrnap](https://github.com/tseemann/barrnap) (v.0.9)
* [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) (v.0.2.6)
* [vsearch](https://github.com/torognes/vsearch) (v2.7.1)

### Database / Datasets download

#### Alves *et al*. database

The *amoA* database was downoladed from the supplementary informatio of Alves *et al*. (https://doi.org/10.1038/s41467-018-03861-1). Supplementary files 1 et 3 were downloaded and unzipped. The files `AamoA.db_an96.aln_tax.annotated.fasta` (supplementary data 1), `AamoA.db_nr.aln.fasta`, and `AamoA.db_nr.aln_taxonomy_qiime.txt` (in the `AamoA.db_nr_qiime.mothur/` folder of supplementary data 3) were moved to `0-databases/`.

A BLAST database was then made using BLAST+ using this command line:
```
makeblastdb -dbtype nucl -in 0-databases/AamoA.db_an96.aln_tax.annotated.fasta -out 0-databases/AamoA.db_an96.aln_tax.annotated
```

#### NCBI nucleotides sequences

Sequences were directly downloaded from this website: https://www.ncbi.nlm.nih.gov/nuccore, using as search query `(Archaea[Organism]) AND 2000:99999999999999[Sequence Length]`. The ouput was then dowloaded using `Send to file > FASTA format`. The sequences were downloaded on the 31st Oct. 2020. 

The fasta file was save as `1-ncbi_arch_2000/1-raw_data/ncbi_nucl_arch_min2000.fasta`. This file containes 335,202 sequences.

```
grep "^>" 1-ncbi_arch_2000/1-raw_data/ncbi_nucl_arch_min2000.fasta | wc -l
> 335202
```

#### RefSeq archaeal genomes

The archaeal RefSeq genomes were downloaded with ncbi-genome-download on the 7th of November 2020, using this commande line:

```
mkdir 2-refseq_arch_genomes/
cd 2-refseq_arch_genomes/
ncbi-genome-download archaea --format fasta --assembly-level all --section refseq --output-folder 1-raw_data/arch_genomes_refseq archaea
```

Then, genomes were unarchived using this script:
```
mkdir 1-raw_data/arch_genomes_refseq_unarchived

for file in 1-raw_data/arch_genomes_refseq/refseq/archaea/GCF_*/*.fna.gz; do
	echo $file
	id=$(echo $file | cut -d "/" -f5)
	echo $id
	cp $file 1-raw_data/arch_genomes_refseq_unarchived/${id}.fna.gz
	gzip -d 1-raw_data/arch_genomes_refseq_unarchived/${id}.fna.gz
done
```

#### GeneBank archaeal genomes

GeneBank genomes were downloaded follwing the instructions  of this github repository : https://github.com/rprops/MetaG_analysis_workflow/wiki/09.-Download-genomes-NCBI-EDI. The genomes were downloaded on the 7th of November 2020 using these command lines:

```
# Create folders
mkdir 3-genbank_arch_genomes/
cd 3-genbank_arch_genomes/
mkdir 1-raw_data/

# Download the list of archeal genomes
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt
mv assembly_summary.txt 1-raw_data/1-assembly_summary.txt

# How many genomes are there?
wc -l 1-raw_data/1-assembly_summary.txt
5738

# Get the ftp links
less 1-raw_data/1-assembly_summary.txt | cut -f20 > 1-raw_data/2-ftp_links.txt

# Download them all!
mkdir 1-raw_data/3-archaeal_genomes

for next in $(cat 1-raw_data/2-ftp_links.txt); do wget -P 1-raw_data/3-archael_genomes "$next"/*genomic.fna.gz; done

# How many genomes were downloaded?
ls -lh 1-raw_data/3-archaeal_genomes/ | grep -v cds | grep -v rna |wc -l
> 5737
```

### Detailed script (and files)

### Relevant output files

### How to annotate your Thaumarchaeota 16S rRNA sequences with this database?
