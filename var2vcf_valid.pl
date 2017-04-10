#!/usr/bin/env perl
use warnings;
use Getopt::Std;
use strict;

our ($opt_d, $opt_v, $opt_f, $opt_h, $opt_H, $opt_p, $opt_q, $opt_F, $opt_S, $opt_Q, $opt_o, $opt_N, $opt_E, $opt_C, $opt_m, $opt_I, $opt_c, $opt_P, $opt_a, $opt_t, $opt_r, $opt_O, $opt_X, $opt_k, $opt_V, $opt_M, $opt_x, $opt_A, $opt_T);
getopts('htaHSCEAP:d:v:f:p:q:F:Q:s:N:m:I:c:r:O:X:k:V:M:x:T:') || Usage();
($opt_h || $opt_H) && Usage();

my $TotalDepth = $opt_d ? $opt_d : 3;
my $VarDepth = $opt_v ? $opt_v : 2;
my $Freq = $opt_f ? $opt_f : 0.02;
my $Pmean = $opt_p ? $opt_p : 8;
my $qmean = $opt_q ? $opt_q : 22.5; # base quality
my $Qmean = $opt_Q ? $opt_Q : 10; # mapping quality
my $GTFreq = $opt_F ? $opt_F : 0.2; # Genotype frequency
my $SN = $opt_o ? $opt_o : 1.5; # Signal to Noise
$opt_I = $opt_I ? $opt_I : 12;
$opt_m = $opt_m ? $opt_m : 5.25;
$opt_c = $opt_c ? $opt_c : 0;
$opt_P = defined($opt_P) ? $opt_P : 1; # Whether to filter pstd = 0 variant.
$opt_T = defined($opt_T) ? $opt_T : 1; # The minimum split read to support SV

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
$sample = $opt_N if ( $opt_N );
(my $sample_nowhitespace = $sample) =~ s/\s/_/g;

print <<VCFHEADER;
##fileformat=VCFv4.1
##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name (with whitespace translated to underscores)">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=END,Number=1,Type=Integer,Description="Chr End Position">
##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=BIAS,Number=1,Type=String,Description="Strand Bias Info">
##INFO=<ID=REFBIAS,Number=1,Type=String,Description="Reference depth by strand">
##INFO=<ID=VARBIAS,Number=1,Type=String,Description="Variant depth by strand">
##INFO=<ID=PMEAN,Number=1,Type=Float,Description="Mean position in reads">
##INFO=<ID=PSTD,Number=1,Type=Float,Description="Position STD in reads">
##INFO=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">
##INFO=<ID=QSTD,Number=1,Type=Float,Description="Quality score STD in reads">
##INFO=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">
##INFO=<ID=ODDRATIO,Number=1,Type=Float,Description="Strand Bias Odds ratio">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mean Mapping Quality">
##INFO=<ID=SN,Number=1,Type=Float,Description="Signal to noise">
##INFO=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency using only high quality bases">
##INFO=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted AF for indels due to local realignment">
##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases to be shifted to 3 prime for deletions due to alternative alignment">
##INFO=<ID=MSI,Number=1,Type=Float,Description="MicroSatellite. > 1 indicates MSI">
##INFO=<ID=MSILEN,Number=1,Type=Float,Description="MicroSatellite unit length in bp">
##INFO=<ID=NM,Number=1,Type=Float,Description="Mean mismatches in reads">
##INFO=<ID=LSEQ,Number=1,Type=String,Description="5' flanking seq">
##INFO=<ID=RSEQ,Number=1,Type=String,Description="3' flanking seq">
##INFO=<ID=GDAMP,Number=1,Type=Integer,Description="No. of amplicons supporting variant">
##INFO=<ID=TLAMP,Number=1,Type=Integer,Description="Total of amplicons covering variant">
##INFO=<ID=NCAMP,Number=1,Type=Integer,Description="No. of amplicons don't work">
##INFO=<ID=AMPFLAG,Number=1,Type=Integer,Description="Top variant in amplicons don't match">
##INFO=<ID=HICNT,Number=1,Type=Integer,Description="High quality variant reads">
##INFO=<ID=HICOV,Number=1,Type=Integer,Description="High quality total reads">
##INFO=<ID=SPLITREAD,Number=1,Type=Integer,Description="No. of split reads supporting SV">
##INFO=<ID=SPANPAIR,Number=1,Type=Integer,Description="No. of pairs supporting SV">
##INFO=<ID=DUPRATE,Number=1,Type=Float,Description="Duplication rate in fraction">
##FILTER=<ID=q$qmean,Description="Mean Base Quality Below $qmean">
##FILTER=<ID=Q$Qmean,Description="Mean Mapping Quality Below $Qmean">
##FILTER=<ID=p$Pmean,Description="Mean Position in Reads Less than $Pmean">
##FILTER=<ID=SN$SN,Description="Signal to Noise Less than $SN">
##FILTER=<ID=Bias,Description="Strand Bias">
##FILTER=<ID=pSTD,Description="Position in Reads has STD of 0">
##FILTER=<ID=d$TotalDepth,Description="Total Depth < $TotalDepth">
##FILTER=<ID=v$VarDepth,Description="Var Depth < $VarDepth">
##FILTER=<ID=f$Freq,Description="Allele frequency < $Freq">
##FILTER=<ID=MSI$opt_I,Description="Variant in MSI region with $opt_I non-monomer MSI or 13 monomer MSI">
##FILTER=<ID=NM$opt_m,Description="Mean mismatches in reads >= $opt_m, thus likely false positive">
##FILTER=<ID=InGap,Description="The variant is in the deletion gap, thus likely false positive">
##FILTER=<ID=InIns,Description="The variant is adjacent to an insertion variant">
##FILTER=<ID=Cluster${opt_c}bp,Description="Two variants are within $opt_c bp">
##FILTER=<ID=LongMSI,Description="The somatic variant is flanked by long A/T (>=14)">
##FILTER=<ID=AMPBIAS,Description="Indicate the variant has amplicon bias.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">
##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">
VCFHEADER

print join("\t", "#CHROM", qw(POS ID REF ALT QUAL FILTER INFO FORMAT), $sample), "\n";

# Exit if we don't have any variants to write
exit(0) unless( %hash );

#my @chrs = map { "chr$_"; } (1..22);
#push(@chrs, "chrX", "chrY", "chrM");
#if ( $opt_C ) {
#    @chrs = (1..22, "X", "Y", "MT");
#}
my @chrs = reorder(keys %hash);

foreach my $chr (@chrs) {
    my @pos = sort { $a <=> $b } (keys %{ $hash{ $chr } });
    my ($pds, $pde) = (0, 0); # previous Deletion variant's start and end
    my ($pis, $pie) = (0, 0); # previous Insertion variant's start and end
    my ($pvs, $pve) = (0, 0); # previous SNV variant's start and end
    my ($pinfo1, $pfilter, $pinfo2) = ("", "", "");
    foreach my $p (@pos) {
	my @tmp = sort { $b->[14] <=> $a->[14] } @{ $hash{ $chr }->{ $p } };
	#my @hds = qw(sp ep refallele varallele tcov cov rfc rrc fwd rev genotype freq bias pmean pstd qual qstd mapq qratio hifreq extrafreq shift3 msi msint nm leftseq rightseq);
	my $ALL = $opt_A ? @tmp + 0 : 1;
	my %seen = ();
	for(my $i = 0; $i < $ALL; $i++) {
	    my ($sample, $gene, $chrt, $start, $end, $ref, $alt, $dp, $vd, $rfwd, $rrev, $vfwd, $vrev, $genotype, $af, $bias, $pmean, $pstd, $qual, $qstd, $sbf, $oddratio, $mapq, $sn, $hiaf, $adjaf, $shift3, $msi, $msilen, $nm, $hicnt, $hicov, $lseq, $rseq, $seg, $type, $gamp, $tamp, $ncamp, $ampflag) = @{ $tmp[$i] };
	    next if ( $seen{ "$chrt-$start-$end-$ref-$alt" } );
	    $seen{ "$chrt-$start-$end-$ref-$alt" } = 1;
	    if ( not defined $type ) { $type = "REF"; }
	    my $isamp = 1 if ( defined($ampflag) );
	    my $rd = $rfwd + $rrev;
	    if ( $oddratio eq "Inf" ) {
		$oddratio = 0;
	    } elsif ( $oddratio < 1 && $oddratio > 0 ) {
		$oddratio = 1/$oddratio;
	    }
	    my @filters = ();
	    if ($dp < $TotalDepth) {
		push(@filters, "d$TotalDepth") unless($hicnt * $hiaf >= 0.5);
	    }
	    if ( $hicnt < $VarDepth ) {
		push(@filters, "v$VarDepth") unless ($hicnt * $hiaf >= 0.5);
	    }
	    push( @filters, "f$Freq") if ($af < $Freq);
	    push( @filters, "p$Pmean") if ($pmean < $Pmean);
	    push( @filters, "pSTD") if ($opt_P && $pstd == 0 && (! $isamp) && $af < 0.35);
	    push( @filters, "q$qmean") if ($qual < $qmean);
	    push( @filters, "Q$Qmean") if ($mapq < $Qmean && $af < 0.8);
	    push( @filters, "SN$SN") if ($sn < $SN);
	    #my $nmadj = length($ref) == length($alt) && length($ref) < 4 ? length($ref) - 1 : 0; # allow more mismatches more MNV, deprecated as it's adjusted in VarDict now.
	    push( @filters, "NM$opt_m") if ($nm >  $opt_m);
	    push( @filters, "MSI$opt_I") if ( ($msi > $opt_I && $msilen > 1 && $af < 0.2 && abs(length($ref) - length($alt)) == $msilen) || ($msi >= 13 && $msilen == 1 && $af <= 0.275 && abs(length($ref) - length($alt)) == $msilen) );


	    push( @filters, "Bias") if ($hiaf < 0.25 && $bias eq "2;1" && $sbf < 0.01 && ($oddratio > 5 || $oddratio == 0) && $end - $start < 100); #|| ($a[9]+$a[10] > 0 && abs($a[9]/($a[9]+$a[10])-$a[11]/($a[11]+$a[12])) > 0.5));
	    #push( @filters, "InGap" ) if ( $type eq "SNV" && (abs($start-$pds) <= 2 || abs($start-$pde) <= 2));
	    #push( @filters, "InIns" ) if ( $type eq "SNV" && (abs($start-$pis) <= 2 || abs($start-$pie) <= 2));
	    #push( @filters, "LongAT") if ($hiaf < 0.5 && (isLongAT($lseq) || isLongAT($rseq)));
	    if ( abs(length($ref)-length($alt)) == $msilen ) {
		if ($hiaf <= 0.275 && $msi >= 13) {
		    push( @filters, "LongMSI" );
		} elsif ( $hiaf <= 0.2 && $msi >= 8 && $msilen > 1 ) {
		    push( @filters, "LongMSI" );
		}
	    }
	    if ( $type eq "SNV" && @filters == 0 && $start - $pvs < $opt_c ) {
		$pfilter = "Cluster${opt_c}bp";
		push( @filters, "Cluster${opt_c}bp");
	    }
	    push( @filters, "AMPBIAS" ) if ( $isamp && (($gamp < $tamp-$ncamp) || $ampflag) );
	    my $filter = @filters > 0 ? join(";", @filters) : "PASS";
	    next if ( $opt_S && $filter ne "PASS" );
	    my $gt = (1-$af < $GTFreq) ? "1/1" : ($af >= 0.5 ? "1/0" : ($af >= $Freq ? "0/1" : "0/0"));
	    $bias =~ s/;/:/;
	    my $QUAL = ($vd le 1) ? 0 : int(log($vd)/log(2) * $qual);
	    my $END = $opt_E ? "" :  ";END=$end";
	    if ( $pinfo1 ) {
		print "$pinfo1\t$pfilter\t$pinfo2\n" unless ($opt_S && $pfilter ne "PASS");
	    }
	    my ($splitreads, $spanpairs, $cluster) = (0, 0, 0);
	    my $SVINFO = "";
	    unless( $isamp ) {
	        ($splitreads, $spanpairs, $cluster) = (split(/-/, $tamp)) if ( defined($tamp) && $tamp =~ /-/ ); 
		if ( $alt =~ /</ ) {
		    next unless ( $splitreads >= $opt_T ); # Ignore SV's without split read support for now until a better criteria
		    my $svlen = $end - $start;
		    $svlen++ if ( $alt =~ /INV/ );
		    $SVINFO = ";SVTYPE=$alt;SVLEN=$svlen";
		    $SVINFO =~ s/<|>//g;
		}
		$SVINFO .= ";SPLITREAD=$splitreads;SPANPAIR=$spanpairs" if ( defined($tamp) );
	    }
	    my $ampinfo = $isamp ? ";GDAMP=$gamp;TLAMP=$tamp;NCAMP=$ncamp;AMPFLAG=$ampflag" : "";
	    my $dupinfo = $isamp ? "" : (defined($gamp) ? ";DUPRATE=$gamp" : "");
	    my $crispr = $isamp ? "" : (defined($ncamp) ? ";CRISPR=$ncamp" : "");
	    ($pinfo1, $pfilter, $pinfo2) = (join("\t", $chr, $start, ".", $ref, $alt, $QUAL), $filter, join("\t", "SAMPLE=$sample_nowhitespace;TYPE=$type;DP=$dp$END;VD=$vd;AF=$af;BIAS=$bias;REFBIAS=$rfwd:$rrev;VARBIAS=$vfwd:$vrev;PMEAN=$pmean;PSTD=$pstd;QUAL=$qual;QSTD=$qstd;SBF=$sbf;ODDRATIO=$oddratio;MQ=$mapq;SN=$sn;HIAF=$hiaf;ADJAF=$adjaf;SHIFT3=$shift3;MSI=$msi;MSILEN=$msilen;NM=$nm;HICNT=$hicnt;HICOV=$hicov;LSEQ=$lseq;RSEQ=$rseq$ampinfo$dupinfo$crispr$SVINFO", "GT:DP:VD:AD:AF:RD:ALD", "$gt:$dp:$vd:$rd,$vd:$af:$rfwd,$rrev:$vfwd,$vrev"));
	    ($pds, $pde) = ($start+1, $end) if ($type eq "Deletion" && $filter eq "PASS" );
	    ($pis, $pie) = ($start-1, $end+1) if ($type eq "Insertion" && $filter eq "PASS" );
	    ($pvs, $pve) = ($start, $end) if ( $type eq "SNV" && $filter eq "PASS");
	}
    }
    if ( $pinfo1 ) {
	print "$pinfo1\t$pfilter\t$pinfo2\n" unless ($opt_S && $pfilter ne "PASS");
    }
}

sub isLongAT {
    my $seq = shift;
    return 1 if ( $seq =~ /T{14,}/ );
    return 1 if ( $seq =~ /A{14,}/ );
    return 0;
}

sub reorder {
    my @chr = @_;
    my @chrn = (); # numeric chromosomes
    my @nonchrn = (); # non-numeric chrosomes
    foreach my $c (@chr) {
	if ( $c =~ /\d/ && $c !~ /_/) {
	    my $t = $c;
	    $t =~ s/\D//g;
	    push(@chrn, [$t, $c]);
	} else {
	    next if ( $c =~ /X/ || $c =~ /Y/ );
	    next if ( $c eq "MT" || $c eq "chrM" );
	    push(@nonchrn, $c);
	}
    }
    @chrn = sort { $a->[0] <=> $b->[0]; } @chrn;
    @chr = map { $_->[1]; } @chrn;
    if ( $hash{ X } ) {
	push(@chr, 'X' );
    } elsif ( $hash{ chrX } ) {
	push(@chr, 'chrX' );
    } 
    if ( $hash{ Y } ) {
	push(@chr, 'Y' );
    } elsif ( $hash{ chrY } ) {
	push(@chr, 'chrY' );
    } 
    if ( $hash{ MT } ) {
	push(@chr, 'MT' );
    } elsif ( $hash{ chrM } ) {
	push(@chr, 'chrM' );
    }
    push ( @chr, @nonchrn );
    return (@chr);
}

sub Usage {
print <<USAGE;
$0 [-hHS] [-p pos] [-q qual] [-d depth] [-v depth] [-f frequency] [-F frequency] vars.txt

The program will convert the variant output from checkVar.pl script into validated VCF file.

Options are:
    -h Print this usage.
    -H Print this usage.
    -S If set, variants that didn't pass filters will not be present in VCF file
    -a For amplicon based variant calling.  Variant not supported by all amplicons will be considered false positve, with filter set to "AMPBIAS".
    -A Indicate to output all variants at the same position.  By default, only the variant with the highest allele frequency is converted to VCF.
    -c  int
        If two seemingly high quality SNV variants are within {int} bp, they're both filtered.  Default: 0, or no filtering
    -I  int
        The maximum non-monomer MSI allowed for a HT variant with AF < 0.5.  By default, 12, or any variants with AF < 0.5 in a region
        with >6 non-monomer MSI will be considered false positive.  For monomers, that number is 13.
    -m  double
        The maximum mean mismatches allowed.  Default: 5.25, or if a variant is supported by reads with more than 5.25 mean mismathes, it'll be considered
        false positive.  Mismatches don't includes indels in the alignment.
    -p  float
    	The minimum mean position of variants in the read.  Default: 8.
    -P  0 or 1
        Whehter to filter variants with pstd = 0.  Default: 1 or yes.  Set it to 0 for targeted PCR based sequencing, where pstd is expected.
    -q  float
    	The minimum mean base quality.  Default to 22.5 for Illumina sequencing
    -Q  float
    	The minimum mapping quality.  Default to 10.0 for Illumina sequencing.  Variants with lower quality (<10) will be kept only if af >= 0.8
    -d  integer
    	The minimum total depth.  Default to 3, might consider lower depth if (high quality variant count)*af >= 0.5
    -v  integer
    	The minimum high quality variant depth.  Default to 2.  Would consider lower depth (high quality variant count)*af >= 0.5
    -f  float
    	The minimum allele frequency.  Default to 0.02
    -o  signal/noise
    	The minimum signal to noise, or the ratio of hi/(lo+0.5).  Default to 1.5, that is both 2 variant reads are high quality.
    -F  float
    	The minimum allele frequency to consider to be homozygous.  Default to 0.2.  Thus frequency < 0.2 will 
	be considered homozygous REF, while frequency > 0.2 will be considered homozygous ALT.
    -N  string
        The sample name to be used directly.
    -E  If set, do not print END tag
    -T  integer
        The minimum number of split reads for SV.  Default: 1.  Change to 0 if you want SV called from discordant pairs only.

AUTHOR
       Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
       Report bugs to zhongwu\@yahoo.com

COPYRIGHT
       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

USAGE
exit(0);
}
