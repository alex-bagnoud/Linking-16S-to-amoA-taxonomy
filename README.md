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

Finally, the results from all 3 runs were merged together 
