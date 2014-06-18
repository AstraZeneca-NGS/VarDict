#!/bin/bash

# Usage: minivardict.sh bam sample bed N freq

BAM=$1
SAMPLE=$2
BED=$3
N=$4
FREQ=$5
GENOME=$6

if [ ! $FREQ ]
    then
        FREQ=0.01
fi

if [ ! $GENOME ]
    then
        GENOME=/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
fi

vardict.pl -c 1 -S 2 -E 3 -g 4 -x 0 -G $GENOME -f $FREQ -k 3 -X 3 -N $SAMPLE -b $BAM ${BED}.$N > ${SAMPLE}_vars.txt.$N
touch vardict.done.$N

