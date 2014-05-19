#!/usr/bin/env perl
# Check the SNV rate in a sample
use warnings;
use lib "/users/kdld047/lib/perl5";
use lib "/users/kdld047/aris/lib";
use Getopt::Std;
use Stat::Basic;
use Fasta;
use strict;

our ($opt_h, $opt_H, $opt_b, $opt_D, $opt_M, $opt_d, $opt_s, $opt_c, $opt_S, $opt_E, $opt_n, $opt_N, $opt_e, $opt_g, $opt_x, $opt_f, $opt_r, $opt_B, $opt_z, $opt_v, $opt_p, $opt_F, $opt_C, $opt_m, $opt_Q, $opt_T, $opt_L);
unless( getopts( 'hHvzpDCFLMd:b:s:e:S:E:n:c:g:x:f:r:B:N:Q:m:T:' )) {
    USAGE();
}
USAGE() if ( $opt_H );
USAGE() unless ( $opt_b );
my $BAM = $opt_b; # the bam file
my $sample = $1 if ( $BAM =~ /([^\/\._]+?)_[^\/]*.bam/ );
if ( $opt_n ) {
    $sample = $1 if ( $BAM =~ /$opt_n/ );
}
$sample = $opt_N if ( $opt_N );
$opt_d = "\t" unless( $opt_d );
my $c_col = $opt_c ? $opt_c - 1 : 2;
my $S_col = $opt_S ? $opt_S - 1 : 6;
my $E_col = $opt_E ? $opt_E - 1 : 7;
my $s_col = $opt_s ? $opt_s - 1 : 9;
my $e_col = $opt_e ? $opt_e - 1 : 10;
my $g_col = $opt_g ? $opt_g - 1 : 12;

$s_col = $S_col if ( $opt_S && (!$opt_s) );
$e_col = $E_col if ( $opt_E && (!$opt_e) );

if ( $opt_L ) {
    $opt_p = 1;
    $c_col = 1 - 1;
    $S_col = 2 - 1;
    $E_col = 2 - 1;
    $s_col = 2 - 1;
    $e_col = 2 - 1;
    $g_col = 3 - 1;
    $opt_d = ":";
}
my $fasta = new Fasta( -fasta => "/users/kdld047/work/NGS/NGS/genomes/hg19_complete/hg19_complete.fa");
#899	NM_007300	chr17	-	41196311	41277500	41197694	41276113	24	41196311,41199659,41201137,41203079,41209068,41215349,41215890,41219624,41222944,41226347,41228504,41231350,41234420,41242960,41243451,41247862,41249260,41251791,41256138,41256884,41258472,41267742,41276033,41277287,	41197819,41199720,41201211,41203134,41209152,41215390,41215968,41219712,41223255,41226538,41228628,41231416,41234592,41243049,41246877,41247939,41249306,41251897,41256278,41256973,41258550,41267796,41276132,41277500,	0	BRCA1	cmpl	cmpl	1,0,1,0,0,1,1,0,1,2,1,1,0,1,1,2,1,0,1,2,2,2,0,-1,
my $SPLICE = defined($opt_x) ? $opt_x : 2;
my $FREQ = $opt_f ? $opt_f : 0.05;
my $BIAS = 0.05; # The cutoff to decide whether a positin has read strand bias
my $MINB = $opt_B ? $opt_B : 2; # The minimum reads for bias calculation
my $MINR = $opt_r ? $opt_r : 2; # The minimum reads for variance allele
if ( $opt_p ) {
    $FREQ = -1;
    $MINR = 0;
}
my $stat = new Stat::Basic;
if ( $opt_h ) {
    print join("\t", qw(Sample Gene Chr Start End Ref Alt Depth AltDepth RefFwdReads RefRevReads AltFwdReads AltRevReads Genotype AF Bias PMean PStd QMean QStd 5pFlankSeq 3pFlankSeq)), "\n";
}
my %snv;
my %all;
while( <> ) {
    chomp;
    next if ( /^#/ );
    next if ( /^browser/i );
    next if ( /^track/i );
    my @A = split(/$opt_d/);
    my ($chr, $cdss, $cdse, $gene) = @A[$c_col, $S_col, $E_col, $g_col];
    my @starts = split(/,/, $A[$s_col]);
    my @ends = split(/,/, $A[$e_col]);
    my @CDS = ();
    my @EXONSEQ = ();
    my %REF;
    $chr = "chr$chr" unless ($chr =~ /^chr/ );
    for(my $i = 0; $i < @starts; $i++) {
        my ($s, $e) = ($starts[$i], $ends[$i]);
	next if ( $cdss > $e ); # not a coding exon
	last if ( $cdse < $s ); # No more coding exon
	$s = $cdss if ( $s < $cdss );
	$e = $cdse if ( $e > $cdse );
	$s -= $SPLICE unless ( $s == $cdss );
	$e += $SPLICE unless ( $e == $cdse );
	$s++ if ( $opt_z );
	push(@CDS, [$s, $e]);
	my $s_start = $s - $SPLICE - 100 < 1 ? 1 : $s - $SPLICE - 100;
	my $exon = $fasta->getSeq( -id => $chr, -ori => 1, -start => $s_start, -end => $e + $SPLICE + 100);
	for(my $i = $s_start; $i <= $s_start + length($exon); $i++) {
	    $REF{ $i } = uc(substr( $exon, $i - ($s_start), 1 ));
	}
	push(@EXONSEQ, $exon);
    }
    #my %cov;
    my %var;
    #my %hash;
    for(my $i = 0; $i < @CDS; $i++) {
	my ($START, $END) = @{$CDS[$i]};
	$chr =~ s/^chr// if ( $opt_C );
	open(SAM, "samtools view $BAM $chr:$START-$END |");
	while( <SAM> ) {
	    chomp;
	    my @a = split(/\t/);
	    next if ( defined($opt_Q) and $a[4] < $opt_Q );
	    if ( $opt_m ) {
	        if ( /XM:i:(\d+)/ ) {  # number of mismatches.  Don't use NM since it includes gaps, which can be from indels
		    next if ( $1 > $opt_m );
		} else {
		    print STDERR "No XM tag for mismatches. $_\n"; #if ( $opt_D );
		}
	    }
	    my $start = $a[3];
	    my $n = 0;
	    my $p = 0;
	    my $dir = $a[1] & 0x10 ? "-" : "+";
	    my @segs = $a[5] =~ /(\d+)[MI]/g; # Only match and insertion counts toward read length
	    my @segs2 = $a[5] =~ /(\d+)[MIS]/g; # For total length, including soft-clipped bases
	    my $rlen = $stat->sum(\@segs); # The read length for matched bases
	    my $rlen2= $stat->sum(\@segs2); # The total length, including soft-clipped bases
	    while( $a[5] =~ /(\d+)([A-Z])/g ) {
		my $m = $1;
		my $C = $2;
		if ( $C eq "N" ) {
		    $start += $m;
		    next;
		} elsif ( $C eq "S" ) {
		    $n += $m;
		    next;
		} elsif ( $C eq "I" ) {
		    my $s = substr($a[9], $n, $m);
		    my $qstr = substr($a[10], $n, $m);
		    if ( $start - 1 >= $START && $start -1 <= $END ) {
			#$hash{ $start - 1 }->{ I }->{ "+$s" }->{ $dir }++;
			#push(@{ $hash{ $start - 1 }->{ I }->{ "+$s" }->{ p } }, $p < $rlen-$p ? $p + 1: $rlen-$p); 
			my $q = 0;
			for(my $i = 0; $i < length($qstr); $i++) {
			    #push(@{ $hash{ $start - 1 }->{ I }->{ "+$s" }->{ q } }, ord(substr($qstr, $i, 1))-33); 
			    $q += ord(substr($qstr, $i, 1))-33; 
			}
			$q /= length($qstr);
			$snv{ "+$s" }->{ cnt }++;
			$all{ ins }->{ cnt }++;
			$q > 25 ? $snv{ "+$s" }->{ hi }++ : $snv{ "+$s" }->{ lo }++;
			$q > 25 ? $all{ ins }->{ hi }++ : $all{ ins }->{ lo }++;
		    }
		    $n += $m;
		    $p += $m;
		    next;
		} elsif ( $C eq "D" ) {
		    my $s = "-$m";
		    if ( $start >= $START && $start <= $END ) {
			#$hash{ $start }->{ $s }->{ $dir }++; 
			#push(@{ $hash{ $start }->{ $s }->{ p } }, $p < $rlen-$p ? $p + 1: $rlen-$p); 
			#push(@{ $hash{ $start }->{ $s }->{ q } }, ord(substr($a[10], $n, 1))-33);
			#push(@{ $hash{ $start }->{ $s }->{ q } }, ord(substr($a[10], $n+1, 1))-33);
			my $k = "-";
			for(my $i = 0; $i < $m; $i++) {
			    #$cov{ $start+$i }->{ $s }++;
			    $k .= $REF{ $start + $i };
			}
			my $q = (ord(substr($a[10], $n, 1))-33 + (ord(substr($a[10], $n+1, 1))-33))/2;
			$snv{ $k }->{ cnt }++;
			$all{ del }->{ cnt }++;
			$q > 25 ? $snv{ $k }->{ hi }++ : $snv{ $k }->{ lo }++;
			$q > 25 ? $all{ del }->{ hi }++ : $all{ del }->{ lo }++;
		    }
		    $start += $m;
		    next;
		}
		for(my $i = 0; $i < $m; $i++) {
		    my $trim = 0;
		    if ( $opt_T ) {
		        if ( $dir eq "+" ) {
			    $trim = 1 if ( $n > $opt_T );
			} else {
			    $trim = 1 if ( $rlen2 - $n > $opt_T );
			}
		    }
		    my $s = substr($a[9], $n, 1);
		    my $q = substr($a[10], $n, 1);
		    unless( $trim ) {
			if ( $start >= $START && $start <= $END ) {
			    #$hash{ $start }->{ $s }->{ $dir }++;
			    #push(@{ $hash{ $start }->{ $s }->{ p } }, $p < $rlen-$p ? $p + 1: $rlen-$p);
			    #push(@{ $hash{ $start }->{ $s }->{ q } }, ord($q)-33);
			    my $k = "$REF{ $start }\t$s";
			    next if ( $s eq "N" || $REF{ $start } eq "N" );
			    $snv{ $k }->{ cnt }++;
			    (ord( $q )-33) > 25 ? $snv{ $k }->{ hi }++ : $snv{ $k }->{ lo }++;
			    if ( $REF{ $start } eq $s ) {
			        $all{ ref }->{cnt}++;
				(ord( $q )-33) > 25 ? $all{ ref }->{ hi }++ : $all{ ref }->{ lo }++;
			    } else {
			        $all{ snv }->{cnt}++;
				(ord( $q )-33) > 25 ? $all{ snv }->{ hi }++ : $all{ snv }->{ lo }++;
			    }
			    #$cov{ $start }->{ $s }++;
			}
		    }
		    $start++ unless( $C eq "I" );
		    $n++ unless( $C eq "D" );
		    $p++ unless( $C eq "D" );
		}
	    }
	}
	close( SAM );
    }
}
print join("\t", qw(Sample Ref Alt Type Count HQCnt LQCnt Hi_Lo Total% HQ% LQ%)), "\n";
while( my ($k, $v) = each %snv ) {
    my ($k1, $k2) = split(/\t/, $k);
    my $type = "ref";
    my $tk = $k;
    if ( $k2 ) {
        $type = "snv" if ( $k1 ne $k2 );
    } else {
        $type = $k1 =~ /^-/ ? "del" : "ins";
        $tk = $k1 =~ /^-/ ? "-\t$k" : "+\t$k";
    }
    $v->{ hi } = 0 unless ( $v->{ hi } );
    $v->{ lo } = 0 unless ( $v->{ lo } );
    $all{ $type }->{ hi } = 1 unless ( $all{ $type }->{ hi } );
    $all{ $type }->{ lo } = 1 unless ( $all{ $type }->{ lo } );
    my @tmp = map { sprintf("%.6f", $_); } ($v->{ cnt }/$all{ $type }->{ cnt }, $v->{ hi }/$all{ $type }->{ hi }, $v->{ lo }/$all{ $type }->{ lo });
    #my @tmp = $k1 eq $k2 ? ($v->{ cnt }/$all{ ref }->{ cnt }, $v->{ hi }/$all{ ref }->{ hi }, $v->{ lo }/$all{ ref }->{ lo }) : ($v->{ cnt }/$all{ diff }->{ cnt }, $v->{ hi }/$all{ diff }->{ hi }, $v->{ lo }/$all{ diff }->{ lo });
    #@tmp = map {sprintf("%.6f", $_);} @tmp;
    print join("\t", $sample, $tk, $type, $v->{ cnt }, $v->{ hi }, $v->{ lo }, sprintf( "%.6f", ($v->{ hi }+1)/($v->{ lo }+1)), @tmp), "\n";
}

my @types = ("ref", "snv", "ins", "del");
my $alnbases = $all{ ref }->{ cnt } + $all{ snv }->{ cnt }; # all aligned bases, excluding indels
my $alnhi = $all{ ref }->{ hi } + $all{ snv }->{ hi }; # all aligned bases, excluding indels
my $alnlo = $all{ ref }->{ lo } + $all{ snv }->{ lo }; # all aligned bases, excluding indels
foreach my $type (@types) {
    my $v = $all{ $type };
    $v->{ cnt } = $v->{ cnt } ? $v->{ cnt } : 0;
    $v->{ hi } = $v->{ hi } ? $v->{ hi } : 0;
    $v->{ lo } = $v->{ lo } ? $v->{ lo } : 0;
    my @tmp = map { sprintf("%.6f", $_); } ($v->{ cnt }/$alnbases, $v->{ hi }/$alnbases, $v->{ lo }/$alnbases);
    print join("\t", $sample, "All", $type, $type, $v->{ cnt }, $v->{ hi }, $v->{ lo }, sprintf( "%.6f", ($v->{ hi }+1)/($v->{ lo }+1)), @tmp), "\n";
}
#AZ01    chr6    106536253       G       38535   2       G/A
#FCB02N4ACXX:3:1106:11473:57062#GATGGTTC 163     chr3    38181980        50      80M188N10M      =       38182274        667     CTGAAGTTGTGTGTGTCTGACCGCGATGTCCTGCCTGGCACCTGTGTCTGGTGTATTGCTAGTGAGCTCATCGAAAAGAGGTGCCGCCGG ___\ceacgbe^cghghhhhhhfhifdhfhhhhefheb_agfe^MWWaegfa_MW\_S\Z\c^ddgece]acbURZ^b_```[bc^ca[a AS:i:-4 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:52C37      YT:Z:UU XS:A:+  NH:i:1  RG:Z:15
#FCB02N4ACXX:3:2206:20108:2526#GATGGTTC  163     chr3    38181981        50      79M188N11M      =       38182275        667     TGAAGTTGTGTGTGTCTGACCGCGATGTCCTGCCTGGCACCTGTGTCTGGTCTATTGCTAGTGAGCTCATCGTAAAGAGGTGCCGCCGGG \YY`c`\ZQPJ`e`b]e_Sbabc[^Ybfaega_^cafhR[U^ee[ec][R\Z\__ZZbZ\_\`Z`d^`Zb]bBBBBBBBBBBBBBBBBBB AS:i:-8 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:72A16A0    YT:Z:UU XS:A:+  NH:i:1  RG:Z:15
sub USAGE {
    print STDERR <<USAGE;
    $0 [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] [-x #_nu] [-g gene] [-f freq] [-r #_reads] [-B #_reads] region_info

    The program will calculate candidate variance for a given region(s) in an indexed BAM file.  The default
    input is IGV's one or more entries in refGene.txt, but can be any regions

    -H Print this help page
    -h Print a header row decribing columns
    -z Indicate wehther is zero-based cooridates, as IGV does.
    -v VCF format output
    -p Do pileup regarless the frequency
    -C Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
    -D Debug mode.  Will print some error messages and append full genotype at the end.
    -M Similar to -D, but will append individual quality and position data instead of mean
    -d delimiter
       The delimiter for split region_info, default to tab "\t"
    -n regular_expression
       The regular expression to extract sample name from bam filenames.  Default to: /([^\/\._]+?)_[^\/]*.bam/
    -N string
       The sample name to be used directly.  Will overwrite -n
    -b string
       The indexed BAM file
    -c INT
       The column for chr
    -S INT
       The column for region start, e.g. gene start
    -E INT
       The column for region end, e.g. gene end
    -s INT
       The column for segment starts in the region, e.g. exon starts
    -e INT
       The column for segment ends in the region, e.g. exon ends
    -g INT
       The column for gene name
    -x INT
       The number of nucleotide to extend for each segment, default: 2
    -f double
       The threshold for allele frequency, default: 0.05 or 5%
    -F Indicate to calculate the frequency as the sum of all non-reference variants, 
       instead of just the most frequent allele, which is the default
    -r minimum reads
       The minimum # of variance reads, default 3
    -B INT
       The minimum # of reads to determine strand bias, default 3
    -Q INT
       If set, reads with mapping quality less than INT will be filtered and ignored
    -m INT
       If set, reads with mismatches more than INT will be filtered and ignored.  Gaps are not counted as mismatches.  
       Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem currently doesn't have such SAM tag
    -T INT
       Trim bases after [INT] bases in the reads
    -L Used for command line pipe, such as "echo chr:pos:gene | checkVar.pl -L".  Will automatically set "-d : -p -c 1 -S 2 -E 2 -g 3"
USAGE
   exit(0);
}
