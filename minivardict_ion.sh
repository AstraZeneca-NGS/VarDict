#!/bin/bash

# Usage: minivardict.sh bam sample bed N freq

BAM=$1
SAMPLE=$2
BED=$3
N=$4
FREQ=$5

if [ ! $FREQ ]
    then
        FREQ=0.01
fi

vardict.pl -q 20 -c 1 -S 2 -E 3 -g 4 -x 0 -f $FREQ -k 3 -X 3 -N $SAMPLE -b $BAM ${BED}.$N > ${SAMPLE}_vars.txt.$N
touch vardict.done.$N

