# This 'sh' file is not executable
# Please copy and paste each of the command in a UNIX terminal.

# 1) Download genomes

# https://github.com/kblin/ncbi-genome-download
ncbi-genome-download archaea --format fasta --assembly-level all --section refseq --output-folder 1-raw_data/arch_genomes_refseq archaea

mkdir 1-raw_data/arch_genomes_refseq_unarchived

for file in 1-raw_data/arch_genomes_refseq/refseq/archaea/GCF_*/*.fna.gz; do
	echo $file
	id=$(echo $file | cut -d "/" -f5)
	echo $id
	cp $file 1-raw_data/arch_genomes_refseq_unarchived/${id}.fna.gz
	gzip -d 1-raw_data/arch_genomes_refseq_unarchived/${id}.fna.gz
done

# How many genomes are there?
ls -l 1-raw_data/arch_genomes_refseq_unarchived/ | wc -l

# 2) blast amoa genes

mkdir 2-amoa_blast

# This script blast each genome in the folder to the archaeal amoA db and keeps the best hit
# It adds to the last line the path to the genome

#for file in ~/databases/archaeal_refseq_genomes/*.fna; do

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


## 2.2) Extract list of genomes that harbor an amoA

less 2-amoa_blast/3-cat_blast.txt | cut -f13 | sort -u > 2-amoa_blast/4-amoa_genome_list.txt

# 3) Get the amoA sequences

mkdir 3-amoa_seqs

## 3.1) Extract sequences

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


## 3.2) Annotate them

db_seq="../0-databases/AamoA.db_nr.aln.fasta"
qiime_tax="../0-databases/AamoA.db_nr.aln_taxonomy_qiime.txt"


source activate qiime1
assign_taxonomy.py -i 3-amoa_seqs/1-amoa.fasta -t $qiime_tax -r $db_seq --similarity 0.9 \
-o 3-amoa_seqs/2-uclust_annotations/
source deactivate qiime1

# 4) Extract 16S from amoA genomes

mkdir 4-16S_genes

## 4.1) Barrnap

while read p; do
	id=$(echo ${p##*\/} | sed 's/.fna//')
	barrnap $p --kingdom arc| grep 16S > 4-16S_genes/1-${id}_barrnap.txt
done < 2-amoa_blast/4-amoa_genome_list.txt

find 4-16S_genes/ -size 0 -delete

## 4.2) Get the sequences

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

# 5) Run 'scripts/script2b_amoa_16S_from_genomes.R' to merge all these data and output annotations files

Rscript ./scripts/script2b_amoa_16S_from_genomes.R
