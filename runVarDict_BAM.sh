#!/bin/bash

# Usage: runVarDict_BAM.sh 1.fastq 2.fastq sample target_bed freq genome_seq GRCh37
PATH=$PATH:/group/cancer_informatics/tools_resources/NGS/bin
JBIN=/opt/az/oracle/java/jdk1.7.0_11/bin
SBIN=/group/cancer_informatics/tools_resources/NGS/bin
SNPEFF=/group/cancer_informatics/tools_resources/NGS/snpEff
REF=/ngs/reference_data/genomes/Hsapiens/hg19
BAM=$1
SAMPLE=$2
BED=$3
FREQ=$4     # Optional.  Default to 0.01
GENOME=$5   # Optional.  Default to hg19
GRCh37=$6   # Optional.
BEDBASE=`basename $BED`

if [ ! $FREQ ]
    then
        FREQ=0.01
fi

if [ ! $GENOME ]
    then
        GENOME=/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
fi

if [ $BED ]
    then
	$SBIN/splitBed.pl $BED 8
	for n in {1..8}; do
	    minivardict.sh $BAM ${SAMPLE} $BEDBASE $n $FREQ $GENOME &
	done
	waitVardict.pl vardict 8
	cat ${SAMPLE}_vars.txt.[1-9] > ${SAMPLE}_vars.txt
	$SBIN/teststrandbias.R ${SAMPLE}_vars.txt > ${SAMPLE}_vars.txt.t
	mv ${SAMPLE}_vars.txt.t ${SAMPLE}_vars.txt
	$SBIN/var2vcf_valid.pl -S -f $FREQ $GRCh37 ${SAMPLE}_vars.txt > ${SAMPLE}_vars.vcf
	$JBIN/java -Xmx4g -jar $SNPEFF/snpEff.jar eff -c $SNPEFF/snpEff.config -d -v -canon hg19 ${SAMPLE}_vars.vcf > ${SAMPLE}_vars.eff.vcf
	$JBIN/java -Xmx4g -jar $SNPEFF/SnpSift.jar annotate -v $REF/variation/dbsnp_latest.vcf ${SAMPLE}_vars.eff.vcf > ${SAMPLE}_vars.eff.dbsnp.vcf
	$JBIN/java -Xmx4g -jar $SNPEFF/SnpSift.jar annotate -v $REF/variation/clinvar_latest.vcf ${SAMPLE}_vars.eff.dbsnp.vcf > ${SAMPLE}_vars.eff.dbsnp.clin.vcf
	$JBIN/java -Xmx4g -jar $SNPEFF/SnpSift.jar annotate -v $REF/variation/CosmicCodingMuts_latest.vcf ${SAMPLE}_vars.eff.dbsnp.clin.vcf > ${SAMPLE}_vars.eff.dbsnp.clin.cosmic.vcf
	#rm ${SAMPLE}_vars.txt.[1-9]
	rm vardict.done.*
	touch vardict.done

	rm ${BEDBASE}.[1-8]
	rm snpEff_summary.html
	rm snpEff_genes.txt

fi

touch runVarDict_BAM.done

#samtools view $3.bam | mapping_stat_BWA.pl ${3}_bam > stat_bam.txt
#samtools view ${3}_sorted.bam | mapping_stat_BWA.pl ${3}_sorted > stat_sorted_bam.txt

# Remove duplicates by Picard.  Time: 64.63min, Mem: 22G (21,250,637,824)
#module load java
#java -jar ~/local/picard-tools-1.70/MarkDuplicates.jar INPUT=${3}_sorted.bam OUTPUT=${3}_sorted_NoDup.bam COMMENT="Remove Duplicates Using Picard" REMOVE_DUPLICATES=true METRICS_FILE=duplicates.txt ASSUME_SORTED=true
#samtools view ${3}_sorted_NoDup.bam | mapping_stat_BWA.pl ${3}_NoDup > stat_NoDup.txt
