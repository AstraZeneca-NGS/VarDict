#!/usr/bin/perl -w

use lib '/users/kdld047/lib/perl5';
use Getopt::Std;
use Stat::Basic;
use strict;

our ($opt_F, $opt_R, $opt_r, $opt_n, $opt_N, $opt_H, $opt_q, $opt_p, $opt_b, $opt_f, $opt_c, $opt_u, $opt_D, $opt_Q, $opt_P, $opt_M, $opt_s, $opt_V, $opt_a, $opt_C);
getopts('uHabR:F:f:n:r:p:q:c:D:P:Q:M:s:V:C:') || USAGE();
USAGE() if ( $opt_H );
my $FRACTION = $opt_r ? $opt_r : 0.4;
my $MAXRATIO = $opt_R ? $opt_R : 1;
my $CNT = $opt_n ? $opt_n : 10;
my $AVEFREQ = $opt_f ? $opt_f : 0.15;
my $MINPMEAN = $opt_p ? $opt_p : 5;
my $MINQMEAN = $opt_q ? $opt_q : 25;
my $FILPMEAN = $opt_P ? $opt_P : 0; # will be filtered on the first place
my $FILQMEAN = $opt_Q ? $opt_Q : 0; # will be filtered on the first place
my $FILDEPTH = $opt_D ? $opt_D : 0; # will be filtered on the first place
my $MINFREQ = $opt_F ? $opt_F : 0.05;
my $MINMQ = $opt_M ? $opt_M : 20;
my $MINVD = $opt_V ? $opt_V : 2; # minimum variant depth
my $MAF = 0.0025;
my $SN = $opt_s ? $opt_s : 4;
my $control = $opt_c;


# SBF: Strand Bias Fisher Exact test
my @columns = qw(CDS AA END DP AF BIAS PMEAN PSTD QUAL QSTD SBF GMAF VD CLNSIG ODDRATIO HIAF MQ SN ADJAF SHIFT3 MSI dbSNPBuildID);
push(@columns, qw(GDAMP TLAMP NCAMP AMPFLAG)) if ( $opt_a );
my @amphdrs = $opt_a ? qw(GAmplicons TAmplicons NCAmplicons Ampflag) : ();
my @appcols = $opt_C ? split(/:/, $opt_C) : ();
push(@columns, @appcols) if ($opt_C);
print join("\t", qw(Sample Chr Start ID Ref Alt Type Effect Functional_Class Codon_Change Amino_Acid_Change Amino_Acid_Length Gene Transcript_bioType Gene_Coding Transcript Exon COSMIC_CDS_Change COSMIC_AA_Change End Depth AlleleFreq Bias Pmean Pstd Qual Qstd SBF GMAF VD CLNSIG ODDRATIO HIAF MQ SN AdjAF Shift3 MSI dbSNPBuildID), @appcols, @amphdrs, qw(N_samples N_Var Pcnt_sample Ave_AF PASS Var_Class)), "\n";
my @data;
my %sample;
my %var;
my $stat = new Stat::Basic;
my %CONTROL;
while( <> ) {
    next if ( /^#/ );
    chomp;
    my @a = split(/\t/);
    $a[7] .= ";";
    my %d;
    while( $a[7] =~ /([^=;]+)=([^=]+?);/g ) {
	$d{ $1 } = $2;
    }

    $d{ SBF } = $d{ SBF } < 0.0001 ? sprintf("%.1e", $d{ SBF }) : sprintf("%.4f", $d{ SBF }) if ( $d{ SBF } );
    $d{ ODDRATIO } = sprintf("%.3f", $d{ ODDRATIO }) if ( $d{ ODDRATIO } );
    my @effs = split(/,/, $d{ EFF });
    my $vark = join(":", @a[0,1,3,4]); # Chr Pos Ref Alt
    next if ( $FILDEPTH && $d{ DP } < $FILDEPTH );
    next if ( $FILPMEAN && $d{ PMEAN } < $FILPMEAN );
    next if ( $FILQMEAN && $d{ QUAL } < $FILQMEAN );
    if ( $control && $control eq $d{ SAMPLE } ) {
	my ($pmean, $qmean) = ($d{ PMEAN }, $d{ QUAL });
	my $pass = "TRUE";
	#$pass = "FALSE" unless ( $d{PSTD} > 0 );
	$pass = "FALSE" if ($qmean < $MINQMEAN );
	$pass = "FALSE" if ($pmean < $MINPMEAN );
	$pass = "FALSE" if ( $d{AF} < $MINFREQ );
	$pass = "FALSE" if ( $d{MQ} < $MINMQ );
	$pass = "FALSE" if ( $d{SN} < $SN );
	$pass = "FALSE" if ( $d{VD} && $d{VD} < $MINVD );
	my $class = $a[2] =~ /COSM/ ? "COSMIC" : ($a[2] =~ /^rs/ ? (checkCLNSIG($d{CLNSIG}) ? "ClnSNP" : "dbSNP") : "Novel");
        #$CONTROL{ $vark } = 1 if ( $pass eq "TRUE" && ($class ne "dbSNP" && $class ne "ClnSNP"));  # so that any novel or COSMIC variants showed up in control won't be filtered
        $CONTROL{ $vark } = 1 if ( $pass eq "TRUE" && $class eq "Novel");  # so that any novel or COSMIC variants showed up in control won't be filtered
    }
    unless( $opt_u && $d{ SAMPLE } =~ /Undetermined/i ) { # Undetermined won't count toward samples
	$sample{ $d{ SAMPLE } } = 1;
	push( @{ $var{ $vark } }, $d{ AF } );
    }
    foreach my $eff (@effs) {
        $eff =~ s/\)$//;
	my @e = split(/\|/, $eff, -1);

	# Move the aa position in multiple aa changes if they're silent.  e.g. GC796GS will become C797S
	if ( $e[3] && $e[3] =~ /^([A-Z]+)(\d+)([A-Z]+)$/ ) {
	    my ($aa1, $aap, $aa2) = ($1, $2, $3);
	    my $an = 0;
	    $an++ while($an < length($aa1)-1 && $an < length($aa2)-1 && substr($aa1, $an, 1) eq substr($aa2, $an, 1));
	    if ( $an ) {
	        $aa1 = substr($aa1, $an);
	        $aa2 = substr($aa2, $an);
		$aap += $an;
		$e[3] = "$aa1$aap$aa2";
	    }
	}
	my ($type, $effect) = split(/\(/, $e[0]);
	my @tmp = map { defined($d{ $_ }) ? $d{ $_ } : ""; } @columns;
	my @tmp2= map { defined($_) ? $_ : ""; } @e[1..9];
	push(@data, [$d{ SAMPLE }, @a[0..4], $type, $effect, @tmp2, @tmp]);
    }
}

my @samples = keys %sample;
my $sam_n = @samples + 0;
foreach my $d (@data) {
    my $vark = join(":", @$d[1, 2, 4, 5]); # Chr Pos Ref Alt
    next unless( $var{ $vark } ); # Likely just in Undetermined.
    my ($pmean, $qmean) = @$d[23,25];
    my $varn = @{ $var{ $vark } } + 0;
    my $ave_af = $stat->mean( $var{ $vark } );
    my $pass = ($varn/$sam_n > $FRACTION && $varn >= $CNT && $ave_af < $AVEFREQ && $d->[3] eq ".") ? "MULTI" : "TRUE"; # novel and present in $MAXRATIO samples
    #$pass = "FALSE" unless ( $d->[24] > 0 ); # all variants from one position in reads
    $pass = "DUP" if ( $d->[24] ==  0 && $d->[22] !~ /1$/ && $d->[22] !~ /0$/ && (!$opt_a) && $d->[21] < 0.35 ); # all variants from one position in reads
    $pass = "MAXRATE" if ( $varn/$sam_n >= $MAXRATIO && $d->[21] < 0.35 ); # present in $MAXRATIO samples, regardless of frequency
    $pass = "QMEAN" if ($qmean < $MINQMEAN );
    $pass = "PMEAN" if ($pmean < $MINPMEAN );
    $pass = "MQ" if ( $d->[33] < $MINMQ );
    $pass = "SN" if ( $d->[34] < $SN );
    $pass = "MINFREQ" if ( $d->[21] < $MINFREQ );
    $pass = "MINVD" if ( $d->[29] && $d->[29] < $MINVD );
    my $class = $d->[3] =~ /COSM/ ? "COSMIC" : ($d->[3] =~ /^rs/ ? (checkCLNSIG($d->[30]) == 1 ? "ClnSNP" : "dbSNP") : "Novel");

    # Rescue deleterious dbSNP, such as rs80357372 (BRCA1 Q139* that is in dbSNP, but not in ClnSNP or COSMIC
    if ( ($d->[6] eq "STOP_GAINED" || $d->[6] eq "FRAME_SHIFT") && $class eq "dbSNP" ) {
	my $pos = $1 if ( $d->[10] =~ /(\d+)/ );
	$class = "dbSNP_del" if ( $pos/$d->[11] < 0.95 );
    }

    $class = "dbSNP" if ( $d->[28] && $d->[28] > $MAF ); # if there's MAF with frequency, it'll be considered dbSNP regardless of COSMIC
    $pass = "CNTL" if ( $control && $CONTROL{ $vark } );
    $pass = "BIAS" if ( $opt_b && ($class eq "Novel"||$class eq "dbSNP") && ($d->[22] eq "2;1" || $d->[22] eq "2;0") && $d->[21] < 0.3 ); # Filter novel variants with strand bias.
    $pass = "NonClnSNP" if ( checkCLNSIG($d->[30]) == -1 && $class ne "COSMIC" );
    $pass = "AMPBIAS" if ( $opt_a && $d->[39] && $d->[39] < $d->[40] );
    print join("\t", @$d, $sam_n, $varn, sprintf("%.3f", $varn/$sam_n), $ave_af, $pass, $class), "\n";
}

sub checkCLNSIG {
    my $clnsig = shift;
    return 0 unless( $clnsig );
    my @cs = split(/\||,/, $clnsig );
    my $flag255 = 0;
    my $flagno = 0;
    foreach my $cs (@cs) {
        return 1 if ( $cs > 3 && $cs < 7 );
	$flagno++ if ( $cs < 3 );
	$flag255++ if ( $cs == 255 );
    }
    return -1 if ( $flagno > $flag255 );
    return 1 if ( $flag255 );
    return -1;
}

sub USAGE {
print <<USAGE;
    The program will convert an annotated vcf files by snfEFF using dbSNP and COSMIC back to txt format.  It also checks for quality
    and add "PASS" column.  It will not perform any filtering.
    
    Usage: $0 [-H] [-F var_fraction] [-n sample_cnt] [-f freq] [-p pos] [-q quality] vcf_files
    
    The program accepts more than one vcf files.

    Options:
    -H Print this help page
    -b Novel or dbSNP variants with strand bias "2;1" or "2;0" and AF < 0.3 will be considered as false positive
    -u Undeteremined won't be counted for the sample count.
    -a Indicate it's amplicon based calling and variants not supported by all amplicons will be filtered
    -r DOUBLE
        When a novel variant is present in more than [fraction] of samples and mean allele frequency is less than -f, it's 
	considered as likely false positive. Default 0.4.
	Used with -f and -n
    
    -f DOUBLE
        When the ave allele frequency is also below the [freq], the variant is considered likely false positive.  Default 0.15.
	Used with -r and -n

    -n INT
        When the variant is detected in greater or equal [sample_cnt] samples, the variant is considered likely false positive.  Default 10.
	Used with -r and -f

    -R DOUBLE
        When a variant is present in more than [fraction] of samples, and AF < 0.3, it's considered as 
	likely false positive, even if it's in COSMIC. Default 1.

    -F DOUBLE
        When indivisual allele frequency < feq for variants, it was considered likely false poitives. Default: 0.05 or 5%

    -p INT
        The minimum mean position in reads for variants.  Default: 5bp

    -q DOUBLE
        The minimum mean base quality phred score for variants.  Default: 25

    -P INT
        The filtering mean position in reads for variants.  The raw variant will be filtered on first place if the mean 
	posisiton is less then INT.  Default: 0bp

    -Q DOUBLE
        The filtering mean base quality phred score for variants.  The raw variant will be filtered on first place 
	if the mean quality is less then DOUBLE.  Default: 0

    -M DOUBLE
        The filtering mean mapping quality score for variants.  The raw variant will be filtered if the mean mapping quality score is less then specified.
	Default: 20

    -D INT
        The filtering total depth.  The raw variant will be filtered on first place if the total depth is less then INT.  Default: 0

    -V INT
        The filtering variant depth.  Variants with depth < INT will be considered false positive.  Default: 2, or at least 2 reads are needed for a variant

    -s signal
        The signal/noise value.  Default: 4

    -c Control
        The control sample name.  Any novel or COSMIC variants passing all above filters but also detected in Control sample will be deemed considered
	false positive.  Use only when there's control sample.

    -C additional_columns
        Add additional columns in VCF to be appended to the output.  Use : to separate multiple columns.  Only those defined in VCF are allowed.
	
    A novel variant (non-dbSNP, non-COSMIC) is considered false positive if all three conditions (-r -f -n) are met. Any variant meeting the -p
    or -q conditions are also considered likely false positive.  False positive variants are annotated "FALSE" in column PASS, "TRUE" otherwise.
USAGE
    exit(0);
}
