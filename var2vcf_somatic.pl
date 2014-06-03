#!/bin/env perl

use warnings;
use Getopt::Std;
use strict;

our ($opt_d, $opt_v, $opt_f, $opt_h, $opt_H, $opt_p, $opt_q, $opt_F, $opt_S, $opt_Q, $opt_o, $opt_C, $opt_M, $opt_P, $opt_N, $opt_m);
getopts('hHSCMd:v:f:p:q:F:Q:o:P:N:m:') || Usage();
($opt_h || $opt_H) && Usage();

my $TotalDepth = $opt_d ? $opt_d : 8;
my $VarDepth = $opt_v ? $opt_v : 4;
my $Freq = $opt_f ? $opt_f : 0.02;
my $Pmean = $opt_p ? $opt_p : 5;
my $qmean = $opt_q ? $opt_q : 23; # base quality
my $Qmean = $opt_Q ? $opt_Q : 0; # mapping quality
my $GTFreq = $opt_F ? $opt_F : 0.2; # Genotype frequency
my $SN = $opt_o ? $opt_o : 1.5; # Signal to Noise
my $PVAL = defined($opt_P) ? $opt_P : 0.05; # the p-value from fisher test
$opt_m = $opt_m ? $opt_m : 6;

my %hash;
my $sample;
while(<>) {
    chomp;
    my @a = split(/\t/);
    $sample = $a[0];
    my $chr = $a[2];
    #$chr = "chrX" if ( $chr eq "23" );
    #$chr = "chrY" if ( $chr eq "24" );
    #$chr = "chr$chr" if ( $chr !~ /^chr/ );
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
##INFO=<ID=BIAS,Number=1,Type=String,Description="Strand Bias Info">
##INFO=<ID=PMEAN,Number=1,Type=Float,Description="Mean position in reads">
##INFO=<ID=PSTD,Number=1,Type=Float,Description="Position STD in reads">
##INFO=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">
##INFO=<ID=QSTD,Number=1,Type=Float,Description="Quality score STD in reads">
##INFO=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">
##INFO=<ID=ODDRATIO,Number=1,Type=Float,Description="Strand Bias Oddratio">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mean Mapping Quality">
##INFO=<ID=SN,Number=1,Type=Float,Description="Signal to noise">
##INFO=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency using only high quality bases">
##INFO=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted AF for indels due to local realignment">
##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases to be shifted to 3 prime for deletions due to alternative alignment">
##INFO=<ID=MSI,Number=1,Type=Float,Description="MicroSattelite. > 1 indicates MSI">
##INFO=<ID=MSILEN,Number=1,Type=Float,Description="MSI unit repeat length in bp">
##INFO=<ID=LSEQ,Number=1,Type=Float,Description="5' flanking seq">
##INFO=<ID=RSEQ,Number=1,Type=Float,Description="3' flanking seq">
##FILTER=<ID=q$qmean,Description="Mean Base Quality Below $qmean">
##FILTER=<ID=Q$Qmean,Description="Mean Mapping Quality Below $Qmean">
##FILTER=<ID=p$Pmean,Description="Mean Position in Reads Less than $Pmean">
##FILTER=<ID=SN$SN,Description="Signal to Noise Less than $SN">
##FILTER=<ID=Bias,Description="Strand Bias">
##FILTER=<ID=pSTD,Description="Position in Reads has STD of 0">
##FILTER=<ID=d$TotalDepth,Description="Total Depth < $TotalDepth">
##FILTER=<ID=v$VarDepth,Description="Var Depth < $VarDepth">
##FILTER=<ID=f$Freq,Description="Allele frequency < $Freq">
##FILTER=<ID=P$PVAL,Description="Not significant with p-value > $PVAL">
##FILTER=<ID=P0.01Likely,Description="Likely candidate but p-value > 0.01">
##FILTER=<ID=MSIREN,Description="Variant in MSI region with $opt_m non-monomer MSI or 10 monomer MSI">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Variant forward, reverse reads">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference forward, reverse reads">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
VCFHEADER

print join("\t", "#CHROM", qw(POS ID REF ALT QUAL FILTER INFO FORMAT), $sample, $samplem), "\n";
#my @chrs = map { "chr$_"; } (1..22);
#push(@chrs, "chrX", "chrY", "chrM");
#if ( $opt_C ) {
#    @chrs = (1..22, "X", "Y", "MT");
#}
my @chrs = sort (keys %hash);
foreach my $chr (@chrs) {
    my @pos = sort { $a <=> $b } (keys %{ $hash{ $chr } });
    foreach my $p (@pos) {
	my @tmp = sort { $b->[14] <=> $a->[14] } @{ $hash{ $chr }->{ $p } };
	my $d = $tmp[0]; # Only the highest AF get represented
	#my @a = split(/\t/, $d);
	my @a = @$d;
	next if ( $opt_M && $a[51] !~ /Somatic/ );
	my @filters = ();
	if ( $PVAL ) {
	    if ( $a[53] > $PVAL ) {
	        push(@filters, "P$PVAL");
	    } elsif ( $a[51] =~ /Likely/ && $a[53] > 0.01 ) {
	        push(@filters, "P0.01Likely");
	    }
	}
	my $oddratio = $a[25];
	if ( $oddratio eq "Inf" ) {
	    $oddratio = 0;
	} elsif ( $oddratio < 1 && $oddratio > 0 ) {
	    $oddratio = sprintf("%.2f", 1/$oddratio);
	}
	push( @filters, "d$TotalDepth") if ($a[7] < $TotalDepth);
	push( @filters, "v$VarDepth") if ($a[8] < $VarDepth);
	push( @filters, "f$Freq") if ($a[14] < $Freq);
	push( @filters, "p$Pmean") if ($a[16] < $Pmean);
	#push( @filters, "pSTD") if ($a[17] == 0);
	push( @filters, "q$qmean") if ($a[18] < $qmean);
	push( @filters, "Q$Qmean") if ($a[20] < $Qmean);
	push( @filters, "SN$SN") if ($a[21] < $SN);
	push( @filters, "MSIREN") if ( ($a[46] > $opt_m && $a[47] > 1) || ($a[46] > 10 && $a[47] == 1));
	#push( @filters, "Bias") if (($a[15] eq "2;1" && $a[24] < 0.01) || ($a[15] eq "2;0" && $a[24] < 0.01) ); #|| ($a[9]+$a[10] > 0 && abs($a[9]/($a[9]+$a[10])-$a[11]/($a[11]+$a[12])) > 0.5));
	my $filter = @filters > 0 ? join(";", @filters) : "PASS";
	next if ( $opt_S && $filter ne "PASS" );
	my $gt = (1-$a[14] < $GTFreq) ? "1/1" : ($a[14] >= 0.5 ? "1/0" : ($a[14] >= $Freq ? "0/1" : "0/0"));
	my $gtm = (1-$a[33] < $GTFreq) ? "1/1" : ($a[33] >= 0.5 ? "1/0" : ($a[33] >= $Freq ? "0/1" : "0/0"));
	$a[15] =~ s/;/:/;
	#print STDERR join("\t", @a);
	my $qual = $a[8] > $a[27] ? int(log($a[8])/log(2) * $a[18]) : int(log($a[27])/log(2) * $a[37]);
	print  join("\t", $a[2], $a[3], ".", @a[5,6], $qual, $filter, "$a[51];SAMPLE=$a[0];TYPE=$a[52];DP=$a[7];VD=$a[8];AF=$a[14];BIAS=$a[15];PMEAN=$a[16];PSTD=$a[17];QUAL=$a[18];QSTD=$a[19];SBF=$a[24];ODDRATIO=$oddratio;MQ=$a[20];SN=$a[21];HIAF=$a[22];ADJAF=$a[23];SHIFT3=$a[45];MSI=$a[46];MSILEN=$a[47];SSF=$a[53];SOR=$a[54];LSEQ=$a[47];RSEQ=$a[48]", "GT:DP:VD:AD:RD:AF", "$gt:$a[7]:$a[8]:$a[11],$a[12]:$a[9],$a[10]:$a[14]", "$gtm:$a[26]:$a[27]:$a[30],$a[31]:$a[28],$a[29]:$a[33]"), "\n";
	#print  join("\t", $a[2], $a[3], ".", @a[5,6], $qual, $filter, "SOMATIC;DP=$a[7];VD=$a[8];AF=$a[14];ADJAF=$a[22]", "GT:DP:VD:AF", "$gt:$a[7]:$a[8]:$a[14]"), "\n";
    }
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
    -m  int
        The maximum non-monomer MSI allowed for a HT variant with AF < 0.6.  By default, 6, or any variants with AF < 0.6 in a region
	with >6 non-monomer MSI will be considered false positive.  For monomers, that number is 10.
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
    	The minimum total depth.  Default to 8
    -v	integer
    	The minimum variant depth.  Default to 4
    -f	float
    	The minimum allele frequency.  Default to 0.02
    -o	signal/noise
    	The minimum signal to noise, or the ratio of hi/(lo+0.5).  Default to 1.5.  Set it higher for deep sequencing.
    -F	float
    	The minimum allele frequency to consider to be homozygous.  Default to 0.2.  Thus frequency > 0.8 (1-0.2) will 
	be considered homozygous "1/1", between 0.5 - (1-0.2) will be "1/0", between (-f) - 0.5 will be "0/1",
	below (-f) will be "0/0".
USAGE
exit(0);
}
