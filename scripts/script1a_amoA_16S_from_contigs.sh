# This 'sh' file is not executable
# Please copy and paste each of the command in a UNIX terminal.

## 1) Download NCBI nucleotides sequences

# https://www.ncbi.nlm.nih.gov/nuccore
# NCBI query "(Archaea[Organism]) AND 2000:99999999999999[Sequence Length] "
# Send to file, FASTA format
# Downloaded the 31st Oct. 2020
# Save it as "1-raw_data/ncbi_nucl_arch_min2000.fasta"

# How many sequences were downloaded?
grep "^>" 1-raw_data/ncbi_nucl_arch_min2000.fasta | wc -l


## 2) blast seqs

mkdir 2-amoa_blast

# 2.1)

blastn -query 1-raw_data/ncbi_nucl_arch_min2000.fasta -db ../0-databases/AamoA.db_an96.aln_tax.annotated -outfmt 6 -max_target_seqs 1 -out 2-amoa_blast/1-blast.txt

# 2.2) Extract list of genomes that harbor an amoA

less 2-amoa_blast/1-blast.txt | cut -f1 | sort -u > 2-amoa_blast/2-amoa_ncbi_seq_list.txt

## 3) Get the amoA sequences

mkdir 3-amoa_seqs

## 3.1) Extract sequences

while read p; do
	echo $p
	start=$(echo $p | cut -d " " -f7)
	end=$(echo $p | cut -d " " -f8)
	seq=$(echo $p | cut -d " " -f1)
	echo "samtools faidx 1-raw_data/ncbi_nucl_arch_min2000.fasta ${seq}:${start}-${end} > 3-amoa_seqs/1-${seq}_amoa.fasta"
	samtools faidx 1-raw_data/ncbi_nucl_arch_min2000.fasta ${seq}:${start}-${end} >> 3-amoa_seqs/1-amoa.fasta
done < 2-amoa_blast/1-blast.txt

rm 1-raw_data/ncbi_nucl_arch_min2000.fasta.fai

## 3.2) Annotate them

db_seq="../0-databases/AamoA.db_nr.aln.fasta"
qiime_tax="../0-databases/AamoA.db_nr.aln_taxonomy_qiime.txt"


source activate qiime1
assign_taxonomy.py -i 3-amoa_seqs/1-amoa.fasta -t $qiime_tax -r $db_seq --similarity 0.9 \
-o 3-amoa_seqs/2-uclust_annotations/
source deactivate qiime1

## 4) Extract 16S from amoA genomes

mkdir 4-16S_genes

# 4.1) Maka big fasta file flat

less 1-raw_data/ncbi_nucl_arch_min2000.fasta | awk -v RS='>'     -v FS="\n"     -v OFS=""     -v ORS="" '{ if (NR > 1) { printf ">%s\n",$1; $1=""; printf "%s\n",$0 } }' > 1-raw_data/ncbi_nucl_arch_min2000_flat.fasta

# 4.2) List of all NCBI seq containing an amoA gene

less 3-amoa_seqs/2-uclust_annotations/1-amoa_tax_assignments.txt | cut -d ":" -f1 | sed s'/^/^>/' >  4-16S_genes/1-ncbi_seq_amoa_grep_list.txt
# less 3-amoa_seqs/2-uclust_annotations/1-amoa_tax_assignments.txt | cut -d ":" -f1 | sed s'/^/"^>/' | sed s'/$/"/' > 4-16S_genes/1-ncbi_seq_amoa_grep_list.txt

# 4.3) Extract sequences from original file

grep -f 4-16S_genes/1-ncbi_seq_amoa_grep_list.txt -A1 1-raw_data/ncbi_nucl_arch_min2000_flat.fasta > 4-16S_genes/2a-amoa_contigs.fasta
grep -v "^--" 3-16S_genes/2a-amoa_contigs.fasta > 3-16S_genes/2b-amoa_contigs_parsed.fasta

# 4.4) Barrnap

barrnap 4-16S_genes/2b-amoa_contigs_parsed.fasta --kingdom arc | grep 16S > 4-16S_genes/4-barrnap.txt


# 4.2) Get the sequences

while read p; do
	seq=$(echo $p | cut -d " " -f1)
	start=$(echo $p | cut -d " " -f4)
	end=$(echo $p | cut -d " " -f5)
	echo "samtools faidx ~/databases/ncbi_nucl_arch_min2000.fasta ${seq}:${start}-${end} >> 4-16S_genes/5-16S.fasta"
	samtools faidx 1-raw_data/ncbi_nucl_arch_min2000.fasta ${seq}:${start}-${end} >> 4-16S_genes/5-16S.fasta
done < 4-16S_genes/4-barrnap.txt


rm 1-raw_data/ncbi_nucl_arch_min2000.fasta.fai

## 5) Run 'script2_amoa_16S_from_contigs.R' to merge all these data and output annotations files

Rscript ../scripts/script2_amoa_16S_from_contigs.R






