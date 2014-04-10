#!/usr/bin/env perl
use warnings;
use Getopt::Std;
use strict;

our ($opt_d, $opt_v, $opt_f, $opt_h, $opt_H, $opt_p, $opt_q, $opt_F, $opt_S, $opt_Q, $opt_s, $opt_N, $opt_E);
getopts('hHSd:v:f:p:q:F:Q:s:N:E') || Usage();
($opt_h || $opt_H) && Usage();

my $TotalDepth = $opt_d ? $opt_d : 5;
my $VarDepth = $opt_v ? $opt_v : 2;
my $Freq = $opt_f ? $opt_f : 0.02;
my $Pmean = $opt_p ? $opt_p : 5;
my $qmean = $opt_q ? $opt_q : 25; # base quality
my $Qmean = $opt_Q ? $opt_Q : 15; # mapping quality
my $GTFreq = $opt_F ? $opt_F : 0.02; # Genotype frequency
my $SN = $opt_s ? $opt_s : 4; # Signal to Noise

my %hash;
my $sample;
my @chrs;
$sample = $opt_N if ( $opt_N );
while(<>) {
    chomp;
    my @a = split(/\t/);
    $sample = $a[0];
    my $chr = $a[2];
    #$chr = "chrX" if ( $chr eq "23" );
    #$chr = "chrY" if ( $chr eq "24" );
    #$chr = "chr$chr" if ( $chr !~ /^chr/ );
    push( @{ $hash{ $chr }->{ $a[3] } }, $_ );
    if (not grep /$chr/, @chrs) {
        push (@chrs, $chr);
    }
}

print <<VCFHEADER;
##fileformat=VCFv4.1
##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name">
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
#my @chrs = map { "chr$_"; } (1..22);
#push(@chrs, "chrX", "chrY", "chrM");

foreach my $chr (@chrs) {
    my @pos = sort { $a <=> $b } (keys %{ $hash{ $chr } });
    foreach my $p (@pos) {
        foreach my $d (@{ $hash{ $chr }->{ $p } }) {
            my @a = split(/\t/, $d);
            my $oddratio = $a[21];
            if ( $oddratio eq "Inf" ) {
                $oddratio = 0;
            } elsif ( $oddratio < 1 && $oddratio > 0 ) {
                $oddratio = 1/$oddratio;
            }
            my @filters = ();
            push( @filters, "d$TotalDepth") if ($a[7] < $TotalDepth);
            push( @filters, "v$VarDepth") if ($a[8] < $VarDepth);
            push( @filters, "f$Freq") if ($a[14] < $Freq);
            push( @filters, "p$Pmean") if ($a[16] < $Pmean);
            push( @filters, "pSTD") if ($a[17] == 0);
            push( @filters, "q$qmean") if ($a[18] < $qmean);
            push( @filters, "Q$Qmean") if ($a[22] < $Qmean);
            push( @filters, "SN$SN") if ($a[23] < $SN);
            push( @filters, "Bias") if (($a[15] eq "2;1" && $a[20] < 0.01) || ($a[15] eq "2;0" && $a[20] < 0.01) || ($a[9]+$a[10] > 0 && abs($a[9]/($a[9]+$a[10])-$a[11]/($a[11]+$a[12])) > 0.5));
            my $filter = @filters > 0 ? join(";", @filters) : "PASS";
            next if ( $opt_S && $filter ne "PASS" );
            my $gt;
            if (1 - $a[14] < $GTFreq) {
                $gt = "1/1";
            } elsif ($a[14] >= 0.5) {
                $gt = "1/0";
            } elsif ($a[14] > $GTFreq) {
                $gt = "0/1";
            } else {
                $gt = "0/0";
            }
            $a[15] =~ s/;/:/;
            my $qual = int(log($a[8])/log(2) * $a[18]);
            my $END = $opt_E ? "" :  ";END=$a[4]";
            print  join("\t", $a[2], $a[3], ".", @a[5,6], $qual, $filter, "SAMPLE=$a[0];DP=$a[7]$END;VP=$a[8];AF=$a[14];BIAS=$a[15];REFBIAS=$a[9]:$a[10];VARBIAS=$a[11]:$a[12];PMEAN=$a[16];PSTD=$a[17];QUAL=$a[18];QSTD=$a[19];SBF=$a[20];ODDRATIO=$oddratio;MQ=$a[22];SN=$a[23];HIAF=$a[24];ADJAF=$a[25];SHIFT3=$a[26];MSI=$a[27];LSEQ=$a[28];RSEQ=$a[29]", "GT:DP:VP:AF", "$gt:$a[7]:$a[8]:$a[14]"), "\n";
        }
    }
}
#AZ01	EZH2	chr7	148504716	148504717	AG	A	20852	17250	3495	0	17249	1	-1/G	0.827	1;1	41.5	2.44	CAGAGG	GGGGGA	A:28:F-28:R-0	T:12:F-12:R-0	C:17:F-17:R-0	-2:50:F-50:R-0	G:3495:F-3495:R-0	-1:17250:F-17249:R-1
#AZ01	EZH2	chr7	148506396	148506396	A	C	17774	15940	1801	1	15940	0	C/A	0.897	1;1	34.1	2.31	AAAGGT	CCTACC	A:1802:F-1801:R-1	T:17:F-17:R-0	C:15940:F-15940:R-0	G:9:F-9:R-0	-1:6:F-6:R-0
##fileformat=VCFv4.1
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#chr9    93606272        .       AG      A       .       PASS    DP=47549;END=93606273;AF=0.014;BIAS=2;2;SM=AZ47;PMEAN=33.0;PSTD=9.55
#chr9    93606299        .       T       G       .       PASS    DP=51402;END=93606299;AF=0.013;BIAS=2;1;SM=AZ47;PMEAN=1.2;PSTD=1.58
#chr9    93606397        .       A       C       .       PASS    DP=55670;END=93606397;AF=0.011;BIAS=2;1;SM=AZ47;PMEAN=14.3;PSTD=3.52

sub Usage {
print <<USAGE;
$0 [-hHS] [-p pos] [-q qual] [-d depth] [-v depth] [-f frequency] [-F frequency] vars.txt

The program will convert the variant output from checkVar.pl script into validated VCF file.

Options are:
    -h Print this usage.
    -H Print this usage.
    -S If set, variants that didn't pass filters will not be present in VCF file
    -p float
    	The minimum mean position of variants in the read.  Default: 5.
    -q float
    	The minimum mean base quality.  Default to 25.0 for Illumina sequencing
    -Q float
    	The minimum mapping quality.  Default to 15.0 for Illumina sequencing
    -d integer
    	The minimum total depth.  Default to 5
    -v integer
    	The minimum variant depth.  Default to 2
    -f float
    	The minimum allele frequency.  Default to 0.02
    -s signal/noise
    	The minimum signal to noise, or the ratio of hi/(lo+0.5).  Default to 4.0, that is both 2 variant reads are high quality.
    -F float
    	The minimum allele frequency to consider to be homozygous.  Default to 0.02.  Thus frequency < 0.02 will 
	   be considered homozygous REF, whilt frequency > 0.98 will be considered homozygous ALT.
    -N string
       The sample name to be used directly.
    -E If set, do not print END tag
USAGE
exit(0);
}
