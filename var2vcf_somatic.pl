#!/usr/bin/env perl

use warnings;
use Getopt::Std;
use strict;

our ($opt_d, $opt_v, $opt_f, $opt_h, $opt_H, $opt_p, $opt_q, $opt_F, $opt_S, $opt_Q, $opt_o, $opt_C, $opt_M, $opt_P, $opt_N, $opt_I, $opt_m, $opt_c, $opt_D);
getopts('hHSCMd:v:f:p:q:F:Q:o:P:N:m:c:I:D:') || Usage();
($opt_h || $opt_H) && Usage();

my $MinDepth = $opt_d ? $opt_d : 7;
my $VarDepth = $opt_v ? $opt_v : 4;
my $Freq = defined($opt_f) ? $opt_f : 0.02;
my $Pmean = defined($opt_p) ? $opt_p : 5;
my $qmean = defined($opt_q) ? $opt_q : 23; # base quality
my $Qmean = defined($opt_Q) ? $opt_Q : 0; # mapping quality
my $GTFreq = defined($opt_F) ? $opt_F : 0.2; # Genotype frequency
my $SN = defined($opt_o) ? $opt_o : 1.5; # Signal to Noise
my $PVAL = defined($opt_P) ? $opt_P : 0.05; # the p-value from fisher test
my $DIFF = defined($opt_D) ? $opt_D : 0.2;
$opt_I = $opt_I ? $opt_I : 6;
$opt_m = $opt_m ? $opt_m : 4;
$opt_c = $opt_c ? $opt_c : 75;

my %hash;
my $sample;
while(<>) {
    chomp;
    next if (/R_HOME/);
    my @a = split(/\t/);
    $sample = $a[0];
    my $chr = $a[2];
    push( @{ $hash{ $chr }->{ $a[3] } }, \@a );
}
my $samplem = "${sample}-match";

if ( $opt_N ) {
    ($sample, $samplem) = split(/\|/, $opt_N);
    $samplem = "${sample}-match" unless( $samplem );
}
print <<VCFHEADER;
##fileformat=VCFv4.1
##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=END,Number=1,Type=Integer,Description="Chr End Position">
##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases to be shifted to 3 prime for deletions due to alternative alignment">
##INFO=<ID=MSI,Number=1,Type=Float,Description="MicroSattelite. > 1 indicates MSI">
##INFO=<ID=MSILEN,Number=1,Type=Float,Description="MSI unit repeat length in bp">
##INFO=<ID=SSF,Number=1,Type=Float,Description="P-value">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Odds ratio">
##INFO=<ID=LSEQ,Number=G,Type=String,Description="5' flanking seq">
##INFO=<ID=RSEQ,Number=G,Type=String,Description="3' flanking seq">
##INFO=<ID=STATUS,Number=1,Type=String,Description="Somatic or germline status">
##FILTER=<ID=q$qmean,Description="Mean Base Quality Below $qmean">
##FILTER=<ID=Q$Qmean,Description="Mean Mapping Quality Below $Qmean">
##FILTER=<ID=p$Pmean,Description="Mean Position in Reads Less than $Pmean">
##FILTER=<ID=SN$SN,Description="Signal to Noise Less than $SN">
##FILTER=<ID=Bias,Description="Strand Bias">
##FILTER=<ID=pSTD,Description="Position in Reads has STD of 0">
##FILTER=<ID=MAF0.05,Description="Matched sample has AF > 0.05, thus not somatic">
##FILTER=<ID=d$MinDepth,Description="Total Depth < $MinDepth">
##FILTER=<ID=v$VarDepth,Description="Var Depth < $VarDepth">
##FILTER=<ID=f$Freq,Description="Allele frequency < $Freq">
##FILTER=<ID=P$PVAL,Description="Not significant with p-value > $PVAL">
##FILTER=<ID=DIFF$DIFF,Description="Non-somatic or LOH and allele frequency difference < $DIFF">
##FILTER=<ID=P0.01Likely,Description="Likely candidate but p-value > 0.01/5**vd2">
##FILTER=<ID=InDelLikely,Description="Likely Indels are not considered somatic">
##FILTER=<ID=MSI$opt_I,Description="Variant in MSI region with $opt_I non-monomer MSI or 10 monomer MSI">
##FILTER=<ID=NM$opt_m,Description="Mean mismatches in reads >= $opt_m, thus likely false positive">
##FILTER=<ID=InGap,Description="The somatic variant is in the deletion gap, thus likely false positive">
##FILTER=<ID=InIns,Description="The somatic variant is adjacent to an insertion variant">
##FILTER=<ID=Cluster${opt_c}bp,Description="Two somatic variants are within $opt_c bp">
##FILTER=<ID=LongAT,Description="The somatic variant is flanked by long A/T (>=14)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Variant forward, reverse reads">
##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted AF for indels due to local realignment">
##FORMAT=<ID=BIAS,Number=1,Type=String,Description="Strand Bias Info">
##FORMAT=<ID=PMEAN,Number=1,Type=Float,Description="Mean position in reads">
##FORMAT=<ID=PSTD,Number=1,Type=Float,Description="Position STD in reads">
##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">
##FORMAT=<ID=QSTD,Number=1,Type=Float,Description="Quality score STD in reads">
##FORMAT=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">
##FORMAT=<ID=ODDRATIO,Number=1,Type=Float,Description="Strand Bias Oddratio">
##FORMAT=<ID=MQ,Number=1,Type=Float,Description="Mean Mapping Quality">
##FORMAT=<ID=SN,Number=1,Type=Float,Description="Signal to noise">
##FORMAT=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency using only high quality bases">
##FORMAT=<ID=NM,Number=1,Type=Float,Description="Mean mismatches in reads">
VCFHEADER

print join("\t", "#CHROM", qw(POS ID REF ALT QUAL FILTER INFO FORMAT), $sample, $samplem), "\n";
my @chrs = sort (keys %hash);
foreach my $chr (@chrs) {
    my @pos = sort { $a <=> $b } (keys %{ $hash{ $chr } });
    my ($pds, $pde) = (0, 0); # previous Deletion variant's start and end
    my ($pis, $pie) = (0, 0); # previous Insertion variant's start and end
    my ($pvs, $pve) = (0, 0); # previous SNV variant's start and end
    my ($pinfo1, $pfilter, $pinfo2) = ("", "", "");
    foreach my $p (@pos) {
	my @tmp = sort { $b->[14] <=> $a->[14] } @{ $hash{ $chr }->{ $p } };
	my $d = $tmp[0]; # Only the highest AF get represented
	my ($sample, $gene, $chrt, $start, $end, $ref, $alt, $dp1, $vd1, $rfwd1, $rrev1, $vfwd1, $vrev1, $gt1, $af1, $bias1, $pmean1, $pstd1, $qual1, $qstd1, $mapq1, $sn1, $hiaf1, $adjaf1, $nm1, $sbf1, $oddratio1, $dp2, $vd2, $rfwd2, $rrev2, $vfwd2, $vrev2, $gt2, $af2, $bias2, $pmean2, $pstd2, $qual2, $qstd2, $mapq2, $sn2, $hiaf2, $adjaf2, $nm2, $sbf2, $oddratio2, $shift3, $msi, $msilen, $lseq, $rseq, $seg, $status, $type, $pvalue, $oddratio)  = @$d;
	#$pvalue *= sqrt(60/($mapq1+length($ref)+length($alt)-1))*$af1;
	my @filters = ();
	if ( $oddratio1 eq "Inf" ) {
	    $oddratio1 = 0;
	} elsif ( $oddratio1 < 1 && $oddratio1 > 0 ) {
	    $oddratio1 = sprintf("%.2f", 1/$oddratio1);
	}
	if ( $oddratio2 eq "Inf" ) {
	    $oddratio2 = 0;
	} elsif ( $oddratio2 < 1 && $oddratio2 > 0 ) {
	    $oddratio2 = sprintf("%.2f", 1/$oddratio2);
	}
	if ($dp1 < $MinDepth) {
	    push( @filters, "d$MinDepth") unless ( $status eq "StrongSomatic" && $pvalue < 0.15 && $af1*$vd1 > 0.5);
	}
	if ($vd1 < $VarDepth) {
	    push( @filters, "v$VarDepth") unless ( $status eq "StrongSomatic" && $pvalue < 0.15 && $af1*$vd1 > 0.5);
	}
	if ( $status =~ /Somatic/ || $status =~ /SampleSpecific/ ) {
	    push( @filters, "f$Freq") if ($af1 < $Freq);
	    push( @filters, "MAF0.05") if ($qual2 >= $qmean && $pmean2 >= $Pmean && $mapq2 >= $Qmean && $sn2 >= $SN && $nm2 < $opt_m && $af2 > 0.05);
	    push( @filters, "p$Pmean") if ($pmean1 < $Pmean);
	    push( @filters, "pSTD") if ($pstd1 == 0 && $vd1 < $MinDepth);
	    push( @filters, "q$qmean") if ($qual1 < $qmean);
	    push( @filters, "Q$Qmean") if ($mapq1 < $Qmean);
	    push( @filters, "SN$SN") if ($sn1 < $SN);
	    push( @filters, "NM$opt_m") if ($nm1 >= $opt_m);
	    push( @filters, "Bias") if (($bias1 eq "2;1" || $bias1 eq "2;0") && $sbf1 < 0.01 && ($oddratio1 > 5 || $oddratio1 == 0));
	} elsif ( $status =~ /LOH/ || $status =~ /Deletion/ ) {
	    push( @filters, "f$Freq") if ($af2 < $Freq);
	    push( @filters, "p$Pmean") if ($pmean2 < $Pmean);
	    push( @filters, "pSTD") if ($pstd2 == 0 && $vd2 < $MinDepth);
	    push( @filters, "q$qmean") if ($qual2 < $qmean);
	    push( @filters, "Q$Qmean") if ($mapq2 < $Qmean);
	    push( @filters, "SN$SN") if ($sn2 < $SN);
	    push( @filters, "NM$opt_m") if ($nm2 >= $opt_m);
	    push( @filters, "Bias") if (($bias2 eq "2;1" || $bias2 eq "2;0") && $sbf2 < 0.01 && ($oddratio2 > 5 || $oddratio2 == 0));
	}
	if ( ($msi > $opt_I && $msilen > 1) || ($msi > 10 && $msilen == 1)) {
	    push( @filters, "MSI$opt_I") unless( $status eq "StrongSomatic" && $pvalue < 0.0005 && $msi < 10);
	}
	#push( @filters, "Bias") if (($a[15] eq "2;1" && $a[24] < 0.01) || ($a[15] eq "2;0" && $a[24] < 0.01) ); #|| ($a[9]+$a[10] > 0 && abs($a[9]/($a[9]+$a[10])-$a[11]/($a[11]+$a[12])) > 0.5));
	if ( $PVAL ) {
	    if ( $pvalue > $PVAL ) {
	        push(@filters, "P$PVAL") unless ($status eq "StrongSomatic" && (($pvalue < 0.15 && $af1 > 0.20) || ($pvalue < 0.10 && $af1 > 0.15 && ($vd1 <= $VarDepth || $type ne "SNV"))));
	        #push(@filters, "P$PVAL");
	    } elsif ( $opt_M && $status =~ /LikelySomatic/ && $pvalue > 0.05/5**$vd2 ) {
	        push(@filters, "P0.01Likely");
	    } elsif ( $opt_M && $status =~ /Likely/ && $type ne "SNV" ) {
	        push(@filters, "InDelLikely") unless(length($ref) <= 2 && length($alt) <= 2);
	    } elsif ( (!$opt_M) && abs($af1 - $af2) < $DIFF && $status !~ /Strong/ && $status !~ /Likely/) {
	        push(@filters, "DIFF$DIFF");
	    }
	}
	if ( @filters == 0 ) {
	    push( @filters, "InGap" ) if ( $pds && $type eq "SNV" && $start <= $pde && $end >= $pds && $status =~ /Somatic/ );
	    push( @filters, "InIns" ) if ( $pis && $type eq "SNV" && $start <= $pie && $end >= $pis && $status =~ /Somatic/ );
	    push( @filters, "LongAT") if (isLongAT($lseq) || isLongAT($rseq));
	}
	my $filter = @filters > 0 ? join(";", @filters) : "PASS";
	my $gt = (1-$af1 < $GTFreq) ? "1/1" : ($af1 >= 0.5 ? "1/0" : ($af1 >= $Freq ? "0/1" : "0/0"));
	my $gtm = (1-$af2 < $GTFreq) ? "1/1" : ($af2 >= 0.5 ? "1/0" : ($af2 >= $Freq ? "0/1" : "0/0"));
	$bias1 =~ s/;/,/;
	$bias2 =~ s/;/,/;
	my $qual = $vd1 > $vd2 ? int(log($vd1)/log(2) * $qual1) : int(log($vd2)/log(2) * $qual2);
	if ( $pfilter eq "PASS" && $pinfo2 =~ /Somatic/ && $pinfo2 =~ /TYPE=SNV/ && $filter eq "PASS" && $status =~ /Somatic/ && $type eq "SNV" && $start - $pvs < $opt_c ) {
	    $pfilter = "Cluster${opt_c}bp";
	    $filter = "Cluster${opt_c}bp";
	}
	if ( $pinfo1 ) {
	    print "$pinfo1\t$pfilter\t$pinfo2\n" unless ( ($opt_M && $pinfo2 !~ /Somatic/) || $opt_S && $pfilter ne "PASS" );
	}
	($pinfo1, $pfilter, $pinfo2) = (join("\t", $chr, $start, ".", $ref, $alt, $qual), $filter, join("\t", "STATUS=$status;SAMPLE=$sample;TYPE=$type;SHIFT3=$shift3;MSI=$msi;MSILEN=$msilen;SSF=$pvalue;SOR=$oddratio;LSEQ=$lseq;RSEQ=$rseq", "GT:DP:VD:AD:RD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM", "$gt:$dp1:$vd1:$vfwd1,$vrev1:$rfwd1,$rrev1:$af1:$bias1:$pmean1:$pstd1:$qual1:$qstd1:$sbf1:$oddratio1:$mapq1:$sn1:$hiaf1:$adjaf1:$nm1", "$gtm:$dp2:$vd2:$vfwd2,$vrev2:$rfwd2,$rrev2:$af2:$bias2:$pmean2:$pstd2:$qual2:$qstd2:$sbf2:$oddratio2:$mapq2:$sn2:$hiaf2:$adjaf2:$nm2"));
	($pds, $pde) = ($start+1, $end) if ($type eq "Deletion");
	($pis, $pie) = ($start-1, $end+1) if ($type eq "Insertion");
	($pvs, $pve) = ($start, $end) if ( $type eq "SNV" && $filter eq "PASS");
    }
    if ( $pinfo1 ) {
	print "$pinfo1\t$pfilter\t$pinfo2\n" unless ( ($opt_M && $pinfo2 !~ /Somatic/) || $opt_S && $pfilter ne "PASS" );
    }
}

sub isLongAT {
    my $seq = shift;
    return 1 if ( $seq =~ /T{14,}/ );
    return 1 if ( $seq =~ /A{14,}/ );
    return 0;
}

sub Usage {
print <<USAGE;
$0 [-hHS] [-p pos] [-q qual] [-d depth] [-v depth] [-f frequency] [-F frequency] vars.txt

The program will convert the variant output from checkVar.pl script into validated VCF file.

Options are:

    -h	Print this usage.
    -H	Print this usage.
    -C  If set, chrosomes will have names of 1,2,3,...,X,Y, instead of chr1, chr2, ..., chrX, chrY
    -S	If set, variants that didn't pass filters will not be present in VCF file
    -M  If set, output only candidate somatic
    -D  float (0-1)
        The minimum allele frequency difference between two samples required in addition to p-value.  Not compitable
	with -M option.  It's for interest of identifying variants with different AF, not just somatic.
    -c  int
        If two somatic candidates are within {int} bp, they're both filtered.  Default: 75
    -I  int
        The maximum non-monomer MSI allowed for a HT variant with AF < 0.6.  By default, 6, or any variants with AF < 0.6 in a region
	with >6 non-monomer MSI will be considered false positive.  For monomers, that number is 10.
    -m  int
        The maximum mean mismatches allowed.  Default: 4, or if a variant is supported by reads with more than 4 mismathes, it'll be considered
	false positive.  Mixmatches don't includes indels in the alignment.
    -N  Name(s)
        The sample name(s).  If only one name is given, the matched will be simply names as "name-match".  Two names
	are given separated by "|", such as "tumor|blood".
    -P  float
        The maximum p-value.  Default to 0.05.  If you want to keep all variants, set it to 0.
    -p	float
    	The minimum mean position of variants in the read.  Default: 5.
    -q	float
    	The minimum mean base quality.  Default to 23.0 for Illumina sequencing
    -Q	float
    	The minimum mapping quality.  Default to 0 for Illumina sequencing
    -d	integer
    	The minimum total depth.  Default to 7
    -v	integer
    	The minimum variant depth.  Default to 3
    -f	float
    	The minimum allele frequency.  Default to 0.02
    -o	signal/noise
    	The minimum signal to noise, or the ratio of hi/(lo+0.5).  Default to 1.5.  Set it higher for deep sequencing.
    -F	float
    	The minimum allele frequency to consider to be homozygous.  Default to 0.2.  Thus frequency > 0.8 (1-0.2) will 
	be considered homozygous "1/1", between 0.5 - (1-0.2) will be "1/0", between (-f) - 0.5 will be "0/1",
	below (-f) will be "0/0".

AUTHOR
       Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
       Report bugs to zhongwu\@yahoo.com

COPYRIGHT
       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

USAGE
exit(0);
}
