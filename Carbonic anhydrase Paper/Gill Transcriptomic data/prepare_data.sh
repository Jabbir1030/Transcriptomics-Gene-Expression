#!/bin/bash
#WARNING: This script is for education purposes. You modify this script at your own risk

#we are interested in the raw-data(samples) and reference genome(fasta and gff) and phenotype information. 
#So we will put them in a directory called resources.All other files will be removed

resourcedir=resources
mkdir $resourcedir

mv data/genes/data.gtf $resourcedir
mv data/genome/data.fa $resourcedir
mv data/geuvadis_phenodata.csv $resourcedir
mv -v data/samples $resourcedir


#remove unwanted files
rm -fr data

echo "data preparating is completed"
echo "all data can be found in the directory : \"resources\" " 
