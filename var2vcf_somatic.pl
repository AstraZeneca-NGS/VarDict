#!/usr/bin/perl -w

use Getopt::Std;
use strict;

our ($opt_d, $opt_v, $opt_f, $opt_h, $opt_H, $opt_p, $opt_q, $opt_F, $opt_S, $opt_Q, $opt_s, $opt_C);
getopts('hHSCd:v:f:p:q:F:Q:s:') || Usage();
($opt_h || $opt_H) && Usage();

my $TotalDepth = $opt_d ? $opt_d : 4;
my $VarDepth = $opt_v ? $opt_v : 2;
my $Freq = $opt_f ? $opt_f : 0.02;
my $Pmean = $opt_p ? $opt_p : 5;
my $qmean = $opt_q ? $opt_q : 23; # base quality
my $Qmean = $opt_Q ? $opt_Q : 0; # mapping quality
my $GTFreq = $opt_F ? $opt_F : 0.2; # Genotype frequency
my $SN = $opt_s ? $opt_s : 1.5; # Signal to Noise

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
    push( @{ $hash{ $chr }->{ $a[3] } }, $_ );
}

print <<VCFHEADER;
##fileformat=VCFv4.1
##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=END,Number=1,Type=Integer,Description="Chr End Position">
##INFO=<ID=VP,Number=1,Type=Integer,Description="Variant Depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=BIAS,Number=1,Type=String,Description="Strand Bias Info">
##INFO=<ID=REFBIAS,Number=1,Type=String,Description="Referece depth by strand">
##INFO=<ID=VARBIAS,Number=1,Type=String,Description="Variant depth by strand">
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
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=VP,Number=1,Type=Integer,Description="Variant Depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
VCFHEADER

print join("\t", "#CHROM", qw(POS ID REF ALT QUAL FILTER INFO FORMAT), $sample), "\n";
my @chrs = map { "chr$_"; } (1..22);
push(@chrs, "chrX", "chrY", "chrM");
if ( $opt_C ) {
    @chrs = (1..22, "X", "Y");
    #push(@chrs, "chrX", "chrY", "chrM");
}
foreach my $chr (@chrs) {
    my @pos = sort { $a <=> $b } (keys %{ $hash{ $chr } });
    foreach my $p (@pos) {
	foreach my $d (@{ $hash{ $chr }->{ $p } }) {
	    my @a = split(/\t/, $d);
	    my $oddratio = $a[25];
	    if ( $oddratio eq "Inf" ) {
	        $oddratio = 0;
	    } elsif ( $oddratio < 1 && $oddratio > 0 ) {
	        $oddratio = sprintf("%.2f", 1/$oddratio);
	    }
	    my @filters = ();
	    push( @filters, "d$TotalDepth") if ($a[7] < $TotalDepth);
	    push( @filters, "v$VarDepth") if ($a[8] < $VarDepth);
	    push( @filters, "f$Freq") if ($a[14] < $Freq);
	    push( @filters, "p$Pmean") if ($a[16] < $Pmean);
	    #push( @filters, "pSTD") if ($a[17] == 0);
	    push( @filters, "q$qmean") if ($a[18] < $qmean);
	    push( @filters, "Q$Qmean") if ($a[20] < $Qmean);
	    push( @filters, "SN$SN") if ($a[21] < $SN);
	    #push( @filters, "Bias") if (($a[15] eq "2;1" && $a[24] < 0.01) || ($a[15] eq "2;0" && $a[24] < 0.01) ); #|| ($a[9]+$a[10] > 0 && abs($a[9]/($a[9]+$a[10])-$a[11]/($a[11]+$a[12])) > 0.5));
	    my $filter = @filters > 0 ? join(";", @filters) : "PASS";
	    next if ( $opt_S && $filter ne "PASS" );
	    my $gt;
	    if (1 - $a[14] < $GTFreq) {
	        $gt = "1/1";
	    } elsif ($a[14] >= 0.5) {
	        $gt = "1/0";
	    } else {
	        $gt = "0/1";
	    }
	    $a[15] =~ s/;/:/;
	    my $qual = $a[50] =~ /Somatic/ ? int(log($a[8])/log(2) * $a[18]) : int(log($a[27])/log(2) * $a[37]);
	    print  join("\t", $a[2], $a[3], ".", @a[5,6], $qual, $filter, "$a[50];SAMPLE=$a[0];TYPE=$a[51];DP=$a[7];VD=$a[8];AF=$a[14];BIAS=$a[15];REFBIAS=$a[9]:$a[10];VARBIAS=$a[11]:$a[12];PMEAN=$a[16];PSTD=$a[17];QUAL=$a[18];QSTD=$a[19];SBF=$a[24];ODDRATIO=$oddratio;MQ=$a[20];SN=$a[21];HIAF=$a[22];ADJAF=$a[23];SHIFT3=$a[45];MSI=$a[46];SSF=$a[52];SOR=$a[53];LSEQ=$a[47];RSEQ=$a[48]", "GT:DP:VP:AF", "$gt:$a[7]:$a[8]:$a[14]"), "\n";
	    #print  join("\t", $a[2], $a[3], ".", @a[5,6], $qual, $filter, "SOMATIC;DP=$a[7];VD=$a[8];AF=$a[14];ADJAF=$a[22]", "GT:DP:VP:AF", "$gt:$a[7]:$a[8]:$a[14]"), "\n";
	}
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
    -S	If set, variants didn't pass filters will not be present in VCF file
    -p	float
    	The minimum mean position of variants in the read.  Default: 5.
    -q	float
    	The minimum mean base quality.  Default to 23.0 for Illumina sequencing
    -Q	float
    	The minimum mapping quality.  Default to 0 for Illumina sequencing
    -d	integer
    	The minimum total depth.  Default to 4
    -v	integer
    	The minimum variant depth.  Default to 2
    -f	float
    	The minimum allele frequency.  Default to 0.02
    -s	signal/noise
    	The minimum signal to noise, or the ratio of hi/(lo+0.5).  Default to 1.5.  Set it higher for deep sequencing.
    -F	float
    	The minimum allele frequency to consider to be homozygous.  Default to 0.2.  Thus frequency < 0.2 will 
	be considered homozygous REF, while frequency > 0.8 will be considered homozygous ALT.
USAGE
exit(0);
}
