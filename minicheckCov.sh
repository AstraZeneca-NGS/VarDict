#!/bin/bash

# Usage: minicheckCov.sh bam sample bed N
SBIN=/group/cancer_informatics/tools_resources/NGS/bin
BAM=$1
SAMPLE=$2
BED=$3
N=$4

$SBIN/checkCov.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $SAMPLE -b $BAM -d 1:5:10:25:50:100:500:1000:5000:10000:50000 ${BED}.$N > ${SAMPLE}_cov.txt.$N
touch checkCov.done.$N

