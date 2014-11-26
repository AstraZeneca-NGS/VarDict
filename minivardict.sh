#!/bin/bash

# Usage: minivardict.sh bam sample bed N freq genome other_opt

BAM=$1
SAMPLE=$2
BED=$3
N=$4
FREQ=$5
GENOME=$6
OTHEROPT=$7

if [ ! $FREQ ]
    then
        FREQ=0.01
fi

if [ ! $GENOME ]
    then
        GENOME=/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
fi

if [ ! -e ${SAMPLE}_vars.txt.$N ]
    then
	vardict.pl -c 1 -S 2 -E 3 -g 4 -x 0 -G $GENOME -f $FREQ -k 3 -X 3 -N $SAMPLE -b $BAM $OTHEROPT ${BED}.$N > ${SAMPLE}_vars.txt.$N
fi

touch vardict.done.$N

