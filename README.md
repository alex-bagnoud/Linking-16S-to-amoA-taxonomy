# linking-16S-to-amoA-taxonomy

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

Finally, the results from all 3 runs were merged together.

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

A BLAST database was then made using BLAST+ using this command line
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

```
ncbi-genome-download --format fasta --assembly-level all --section refseq --output-folder 1-raw_data/arch_genomes_refseq archaea
```



#### GeneBank archaeal genomes

### Detailed script (and files)

### Relevant output files

### How to annotate your Thaumarchaeota 16S rRNA sequences with this database?
