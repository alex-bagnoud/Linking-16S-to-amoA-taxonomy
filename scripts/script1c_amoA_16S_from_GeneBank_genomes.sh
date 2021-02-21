# This 'sh' file is not executable
# Please copy and paste each of the command in a UNIX terminal.

## 1.1) Download NCBI archaeal genomes

## Create folders:

mkdir 3-genbank_arch_genomes/
cd 3-genbank_arch_genomes/
mkdir 1-raw_data/

## Download the list of archeal genomes:

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt
mv assembly_summary.txt 1-raw_data/1-assembly_summary.txt

## How many genomes are there?

wc -l 1-raw_data/1-assembly_summary.txt

## Get the ftp links

less 1-raw_data/1-assembly_summary.txt | cut -f20 > 1-raw_data/2-ftp_links.txt

## Download them all!

mkdir 1-raw_data/3-archaeal_genomes

for next in $(cat 1-raw_data/2-ftp_links.txt); do wget -P 1-raw_data/3-archael_genomes "$next"/*genomic.fna.gz; done
## How many genomes were downloaded?

ls -lh 1-raw_data/3-archaeal_genomes/ | grep -v cds | grep -v rna |wc -l

## 1.2) Prepare genomes

## genomes downloaded on 16th June 2018

## Extract genome list

for file in 1-raw_data/3-archaeal_genomes/*.fna.gz; do
	echo $file | grep -v cds | grep -v rna >> 1-raw_data/4-genome_list.txt
done

## Unzip all genomes sequences

while read p; do
	echo $p
	gzip -d $p
done < 1-raw_data/4-genome_list.txt

## Rename them

for file in 1-raw_data/3-archaeal_genomes/*.fna; do
	new_file=$(echo $file | sed 's/\.[0-9]_.*genomic//')
	mv $file $new_file
done

## Remove gz files

for f in 1-raw_data/3-archaeal_genomes/*.fna.gz; do echo rm "$f"; done
for f in 1-raw_data/3-archaeal_genomes/*.fna.gz; do rm "$f"; done


# 2) blast genomes

mkdir 2-amoa_blast

## 2.1)

# This script blast each genome in the folder to the archaeal amoA db and keeps the best hit
# It adds to the last line the path to the genome

for file in 1-raw_data/3-archaeal_genomes/*.fna; do
	id=$(echo ${file##*\/} | cut -d "." -f1)
	echo $id
	blastn -query $file -db ../0-databases/AamoA.db_an96.aln_tax.annotated \
	-outfmt 6 -max_target_seqs 1 -out 2-amoa_blast/1-${id}_blast.txt
	awk -v a="$file" -F"," 'BEGIN { OFS = "\t" } {$2=a; print}' 2-amoa_blast/1-${id}_blast.txt > 2-amoa_blast/2-${id}_blast2.txt
done

cat 2-amoa_blast/2-* > 2-amoa_blast/3-cat_blast.txt

cd 2-amoa_blast
ls | grep blast2.txt | xargs cat >> 3-cat_blast.txt
cd ..

for f in 2-amoa_blast/[1-2]*.txt; do echo rm "$f"; done
for f in 2-amoa_blast/[1-2]*.txt; do rm "$f"; done


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
	echo "samtools faidx $file ${seq}:${start}-${end} > 3-amoa_seqs/1-${id}_amoa.fasta"
	samtools faidx $file ${seq}:${start}-${end} >> 3-amoa_seqs/1-amoa.fasta
done < 2-amoa_blast/3-cat_blast.txt

rm 1-raw_data/3-archaeal_genomes/*.fna.fai

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
	id=$(echo ${p##*\/} | cut -d "." -f1)
	barrnap $p --kingdom arc | grep 16S > 4-16S_genes/1-${id}_barrnap.txt
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
		echo "samtools faidx 1-raw_data/3-archaeal_genomes/${id}.fna ${seq}:${start}-${end} >> 4-16S_genes/2-${id}_16S.fasta"
		samtools faidx 1-raw_data/3-archaeal_genomes/${id}.fna ${seq}:${start}-${end} >> 4-16S_genes/2-${id}_16S.fasta
	done < $file
done 

rm 1-raw_data/3-archaeal_genomes/*.fna.fai

cat 4-16S_genes/2* > 4-16S_genes/3-all_16S_seq.fasta


# 5) Run 'script2b_amoa_16S_from_genomes.R' to merge all these data and output annotations files

Rscript ./scripts/script2b_amoa_16S_from_genomes.R






