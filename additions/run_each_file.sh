#!/bin/bash
#SCRIPT FOR EACH FILE. It will run R and pl script for each file and create vcf gz and tabix.
DIR_R=r
DIR_VCF=vcf
DIR_VAR=var
FILE=$1
#Extract filename
FILENAME="$(basename -- $FILE)"
FILENAME="${FILENAME%.*}"
#Run R
cat $FILE | testsomatic.R > ../$DIR_R/$FILENAME.r.var
#Run pl scripts
cat ../$DIR_R/$FILENAME.r.var | var2vcf_paired.pl > ../$DIR_VCF/$FILENAME.vcf
#ZIP
bgzip -f ../$DIR_VCF/$FILENAME.vcf
tabix -f -p vcf ../$DIR_VCF/$FILENAME.vcf.gz 
#Delete raw, r and not zipped VCF files
rm ../$DIR_VAR/$FILENAME.var
rm ../$DIR_R/$FILENAME.r.var
rm ../$DIR_VCF/$FILENAME.vcf