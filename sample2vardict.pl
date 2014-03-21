#!/usr/bin/perl -w

# From a sample list to shell script for running variant detection using checkVar.pl

use Getopt::Std;
use strict;

our ($opt_p, $opt_b, $opt_H, $opt_a, $opt_A, $opt_o, $opt_f);
getopts("Hp:b:a:A:o:f:") || USAGE();
$opt_H && USAGE();

my $option = $opt_o ? $opt_o : "-l huge_ram=1 -pe smp 8";

my $prog = $opt_p ? $opt_p : "/group/cancer_informatics/tools_resources/NGS/bin/runVarDict_BAM.sh";
my $bed = $opt_b ? $opt_b : "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/human_hg19_1M_segs.txt";
my $freq = $opt_f ? $opt_f : 0.01;
while( <> ) {
    chomp;
    my ($sample, $bam) = split(/\t/);
    $bam =~ s/\@$//;
    $sample .= $opt_a if ( $opt_a );
    $sample = $opt_A . $sample if ( $opt_A );
    print "if [ ! -d \"$sample\" ]; then\n";
    print "    mkdir $sample\n";
    print "fi\n";
    print "cd $sample\n";
    print "qsub $option -cwd -V -N ${sample}_vardict -S /bin/bash $prog $bam $sample $bed $freq\n";
    print "cd ..\n\n";
}

sub USAGE {
    print <<USAGE;
    Usage: $0 -H -p alignment_program -b bed_file sample2bam

    The program will generate shell scripts for submitting vartect jobs to clusters.  When submitting big jobs, pipe the output
    to addQueue.pl, which will submit jobs only to chara, rask, espo, and orr

    -H Print this help page.
    -p alignment_program
       The ABSOLUTE path to the alignment script.  Default: /group/cancer_informatics/tools_resources/NGS/bin/vartect_hg19.sh,
       which will perform alignment using BWA 0.7.4 and call variants using checkVar.pl
    -b bed_file
       The ABSOLUTE path to the bed file.  If not supplied, it'll make variant calls to whole gene using 1Mb segments.
       Default: /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/human_hg19_1M_segs.txt
    -o options
       Options to be appended to qsub, such as "-l huge_ram=1", "-l mem_free=4G", etc.  Default: "-l huge_ram=1";
    -a string
       A string to append to the sample name, e.g. -mm10 when running against mouse genome
    -A string
       A string to suffix to the sample name, e.g. dis- when running after mouse disambiguation
    -f allele frequency
       The minimum allele frequency to call variants.  Default: 0.01

USAGE
    exit(0);
}
