#!/bin/bash

# Usage: minicheckSNV.sh bam sample bed N
SBIN=/group/cancer_informatics/tools_resources/NGS/bin
BAM=$1
SAMPLE=$2
BED=$3
N=$4

$SBIN/checkSNV.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $SAMPLE -b $BAM ${BED}.$N > ${SAMPLE}_alignment_summary.txt.$N
touch alignment_sum.done.$N

