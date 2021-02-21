# Download the amoA database :
# Download and unzip Supplementary Data 1 and 3 of this article https://doi.org/10.1038/s41467-018-03861-1
# Move files AamoA.db_an96.aln_tax.annotated.fasta (supplementary data 1),
# AamoA.db_nr.aln.fasta, and AamoA.db_nr.aln_taxonomy_qiime.txt
# (in the AamoA.db_nr_qiime.mothur/ folder of supplementary data 3) 0-databases/

# Create a BLAST database

makeblastdb -dbtype nucl -in 0-databases/AamoA.db_an96.aln_tax.annotated.fasta -out 0-databases/AamoA.db_an96.aln_tax.annotated
