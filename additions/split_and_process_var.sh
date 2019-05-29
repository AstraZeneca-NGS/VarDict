#!/bin/bash
#Script for run .R and .pl scripts on big VAR files (for example for gVCF or WGS calling).
#It will split VAR file by chromosomes and run the second script for each file.

#Usage: ./split_and_process_var.sh filename
#Requirements:
# * GNU parallel (can be installed from packages)
# * vcf-tools

INITIAL_VAR=$1
COMBINED="output"
DIR_R=r
DIR_VCF=vcf
DIR_VAR=var
mkdir $DIR_R
mkdir $DIR_VCF
mkdir $DIR_VAR

#Number of  concurrent jobs to run in one time. Can be changed by number of threads.
JOBS=4

#split var by chromosomes to var files in subfolder
awk '{print > "var/"$3".var"}' $INITIAL_VAR

#Run in parallel task for each var file
cd $DIR_VAR
find . -name "*.var" | parallel --jobs $JOBS -I% --max-args 1 ../run_each_file.sh %

#Concatecation of VCFs
cd ../$DIR_VCF
vcfs=(*.vcf.gz)
vcf-concat ${vcfs[@]} > ../$COMBINED.vcf
cd ..
bgzip -f $COMBINED.vcf 
tabix -f -p vcf $COMBINED.vcf.gz 