#!/usr/bin/env perl

use Getopt::Std;
use warnings;
use strict;

our ($opt_F, $opt_R, $opt_s, $opt_r, $opt_n, $opt_N, $opt_H, $opt_q, $opt_p, $opt_b, $opt_f, $opt_c, $opt_u, $opt_D, $opt_Q, $opt_P, $opt_M, $opt_o, $opt_V, $opt_a, $opt_C, $opt_G, $opt_A);
getopts('suHabR:F:f:n:r:p:q:c:D:P:Q:M:o:V:C:G:A:') || USAGE();

my %AA_code = (
    "ALA" => "A",  "ILE" => "I",  "LEU" => "L",  "VAL" => "V",
    "PHE" => "F",  "TRP" => "W",  "TYR" => "Y",  "ASN" => "N",
    "CYS" => "C",  "GLN" => "Q",  "MET" => "M",  "SER" => "S",
    "THR" => "T",  "ASP" => "D",  "GLU" => "E",  "ARG" => "R",
    "HIS" => "H",  "LYS" => "K",  "GLY" => "G",  "PRO" => "P", "TER" => "*"
);

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
my $MINMQ = $opt_M ? $opt_M : 10;
my $MINVD = $opt_V ? $opt_V : 2; # minimum variant depth
my $MAF = $opt_G ? $opt_G : 0.0025;
my $SN = $opt_o ? $opt_o : 1.5;
my @controls = $opt_c ? split(/:/, $opt_c) : ();
my %controls = map { ($_, 1); } @controls;
my %MultiMaf;
$opt_A = defined($opt_A) ? $opt_A : "/ngs/reference_data/genomes/Hsapiens/hg19/variation/dbSNP_multi_mafs_latest.txt";
setupMultiMaf($opt_A) if ( -e $opt_A );


# SBF: Strand Bias Fisher Exact test
my @columns = qw(CDS AA END DP AF BIAS PMEAN PSTD QUAL QSTD SBF CAF VD CLNSIG ODDRATIO HIAF MQ SN ADJAF SHIFT3 MSI dbSNPBuildID);
push(@columns, qw(GDAMP TLAMP NCAMP AMPFLAG)) if ( $opt_a );
my @amphdrs = $opt_a ? qw(GAmplicons TAmplicons NCAmplicons Ampflag) : ();
my @appcols = $opt_C ? split(/:/, $opt_C) : ();
push(@columns, @appcols) if ($opt_C);
print join("\t", qw(Sample Chr Start ID Ref Alt Type Effect Functional_Class Codon_Change Amino_Acid_Change cDNA_Change Amino_Acid_Length Gene Transcript_bioType Gene_Coding Transcript Exon COSMIC_CDS_Change COSMIC_AA_Change End Depth AlleleFreq Bias Pmean Pstd Qual Qstd SBF GMAF VD CLNSIG ODDRATIO HIAF MQ SN AdjAF Shift3 MSI dbSNPBuildID), @appcols, @amphdrs, qw(N_samples N_Var Pcnt_sample Ave_AF PASS Var_Class)), "\n";
my @data;
my %sample;
my %var;
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
#/    my @formats = split(/:/, $a[8]);
#    my @fdata = split(/:/, $a[9]);
#    for(my $i = 0; $i < @formats; $i++) {
#        $d{ $formats[$i] } = $fdata[$i];
#    }

    $d{ SBF } = $d{ SBF } < 0.0001 ? sprintf("%.1e", $d{ SBF }) : sprintf("%.4f", $d{ SBF }) if ( $d{ SBF } );
    $d{ ODDRATIO } = sprintf("%.3f", $d{ ODDRATIO }) if ( $d{ ODDRATIO } );
    my @effs = split(/,/, $d{ EFF });
    my $vark = join(":", @a[0,1,3,4]); # Chr Pos Ref Alt
    next if ( $FILDEPTH && $d{ DP } < $FILDEPTH );
    next if ( $FILPMEAN && $d{ PMEAN } < $FILPMEAN );
    next if ( $FILQMEAN && $d{ QUAL } < $FILQMEAN );
    if ( $controls{ $d{ SAMPLE } } ) {
	my ($pmean, $qmean) = ($d{ PMEAN }, $d{ QUAL });
	my $pass = "TRUE";
	#$pass = "FALSE" unless ( $d{PSTD} > 0 );
	$pass = "FALSE" if ($qmean < $MINQMEAN );
	$pass = "FALSE" if ($pmean < $MINPMEAN );
	$pass = "FALSE" if ( $d{AF} < $MINFREQ );
	$pass = "FALSE" if ( $d{MQ} < $MINMQ && $d{AF} < 0.8 );  # Keep low mapping quality but high allele frequency variants
	$pass = "FALSE" if ( $d{SN} < $SN );
	$pass = "FALSE" if ( $d{VD} && $d{VD} < $MINVD );
	my $class = $a[2] =~ /COSM/ ? "COSMIC" : ($a[2] =~ /^rs/ ? (checkCLNSIG($d{CLNSIG}) ? "ClnSNP" : "dbSNP") : "Novel");
        $CONTROL{ $vark } = 1 if ( $pass eq "TRUE" && $class eq "Novel");  # so that any novel variants showed up in control won't be filtered
    }
    unless( $opt_u && $d{ SAMPLE } =~ /Undetermined/i ) { # Undetermined won't count toward samples
	$sample{ $d{ SAMPLE } } = 1;
	push( @{ $var{ $vark } }, $d{ AF } );
    }
    my @alts = split(/,/, $a[4]);
    foreach my $eff (@effs) {
        $eff =~ s/\)$//;
	my @e = split(/\|/, $eff, -1);

	my ($type, $effect) = split(/\(/, $e[0]);
	my @tmp = map { defined($d{ $_ }) ? $d{ $_ } : ""; } @columns;
	my ($aachg, $cdnachg) = $e[3] ? split("/", $e[3]) : ("", "");
	($aachg, $cdnachg) = ("", $e[3]) if ( $e[3] =~ /^[cn]/ );
	if ( $aachg && $aachg =~ /^p\./ && (! $opt_s )) {
	    $aachg =~ s/^p\.//;
	    if ( $aachg =~ /^([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])$/ ) {
	        $aachg = "$AA_code{ uc($1) }$2$AA_code{ uc($3) }";
		print STDERR "$1 $3\n" unless( $AA_code{ uc($1) } && $AA_code{ uc($3) });
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(\d+)_([A-Z][a-z][a-z])(\d+)del$/ ) {
	        #$aachg = (length($a[3])-length($a[4]))/3 < $4 - $3 + 1 ? "$AA_code{$1}${2}del" : "$AA_code{$1}${2}_$AA_code{$3}${4}del";
	        $aachg = (length($a[3])-length($a[4]))/3 < $4 - $2 + 1 && $4 - $2 == 1 ? "$AA_code{uc($1)}${2}del" : "$AA_code{uc($1)}${2}_$AA_code{uc($3)}${4}del";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(\d+)_([A-Z][a-z][a-z])(\d+)ins([A-Z].*)$/ ) {
	        my $ins = "";
		for(my $i = 0; $i < (length($a[4])-length($a[3]))/3; $i += 3) {
		    $ins .= $AA_code{ uc(substr($5, $i, 3)) };
		}
		$aachg = "$AA_code{uc($1)}${2}_$AA_code{uc($3)}${4}ins$ins";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(\d+)(_.*)?fs$/ ) {
	        $aachg = "$AA_code{uc($1)}${2}fs";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(\d+)del$/ ) {
	        $aachg = "$AA_code{uc($1)}${2}del";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(.*)?(\d+)([\*\?])$/ ) {
	        $aachg = "$AA_code{uc($1)}${3}$4";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z][A-Z]\D*)(\d+)([A-Z][a-z][a-z][A-Z]\D*)$/ ) {
		my ($aa1, $aa2) = ("", "");
		for(my $i = 0; $i < length($1); $i += 3) {
		    $aa1 .= $AA_code{ uc(substr($1, $i, 3)) };
		}
		for(my $i = 0; $i < length($3); $i += 3) {
		    if ( substr($3, $i, 3) eq "ext" ) {
		        $aa2 .= "ext*?";
			last;
		    }
		    $aa2 .= $AA_code{ uc(substr($3, $i, 3)) };
		}
		$aachg = "$aa1$2$aa2";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(\d+)(_([A-Z][a-z][a-z])(\d+))?delins([A-Z].*)?$/ ) {
		my $insaa = "";
		$aachg = $AA_code{ uc($1) } . $2;
		$aachg .= $AA_code{ uc($4) } . $5 if ( $4 );
		if ( $6 ) {
		    for(my $i = 0; $i < length($6); $i += 3) {
			$insaa .= $AA_code{ uc(substr($6, $i, 3)) };
		    }
		}
		$aachg .= $insaa ? "delins$insaa" : "del";
	    } elsif ( $aachg =~ /^Ter(\d+)([A-Z][a-z][a-z])ext\*\?$/ ) {
	        $aachg = "*$1$AA_code{uc($2)}ext*?";
	    } else {
	        print STDERR "New format: $aachg\n";
	    }
	}
	# Move the aa position in multiple aa changes if they're silent.  e.g. GC796GS will become C797S
	if ( $aachg && $aachg =~ /^([A-Z]+)(\d+)([A-Z]+)$/ ) {
	    my ($aa1, $aap, $aa2) = ($1, $2, $3);
	    my $an = 0;
	    $an++ while($an < length($aa1)-1 && $an < length($aa2)-1 && substr($aa1, $an, 1) eq substr($aa2, $an, 1));
	    if ( $an ) {
	        $aa1 = substr($aa1, $an);
	        $aa2 = substr($aa2, $an);
		$aap += $an;
		$aachg = "$aa1$aap$aa2";
	    }
	}
	my @tmp2= map { defined($_) ? $_ : ""; } (@e[1, 2], $aachg, $cdnachg, @e[4..9]);
	push(@data, [$d{ SAMPLE }, @a[0..3], $alts[$e[10]-1], $type, $effect, @tmp2, @tmp]);
    }
}

my @samples = keys %sample;
my $sam_n = @samples + 0;
foreach my $d (@data) {
    my $vark = join(":", @$d[1, 2, 4, 5]); # Chr Pos Ref Alt
    next unless( $var{ $vark } ); # Likely just in Undetermined.
    my ($pmean, $qmean) = @$d[24,26];
    my $varn = @{ $var{ $vark } } + 0;
    my $ave_af = mean( $var{ $vark } );
    my $pass = ($varn/$sam_n > $FRACTION && $varn >= $CNT && $ave_af < $AVEFREQ && $d->[3] eq ".") ? "MULTI" : "TRUE"; # novel and present in $MAXRATIO samples
    #$pass = "FALSE" unless ( $d->[24] > 0 ); # all variants from one position in reads
    $pass = "DUP" if ( $d->[25] ==  0 && $d->[23] !~ /1$/ && $d->[23] !~ /0$/ && (!$opt_a) && $d->[22] < 0.35 ); # all variants from one position in reads
    $pass = "MAXRATE" if ( $varn/$sam_n >= $MAXRATIO && $d->[22] < 0.35 ); # present in $MAXRATIO samples, regardless of frequency
    $pass = "QMEAN" if ($qmean < $MINQMEAN );
    $pass = "PMEAN" if ($pmean < $MINPMEAN );
    $pass = "MQ" if ( $d->[34] < $MINMQ && $d->[22] < 0.8 ); # Keep low mapping quality but high allele frequency variants
    $pass = "SN" if ( $d->[35] < $SN );
    $pass = "MINFREQ" if ( $d->[22] < $MINFREQ );
    $pass = "MINVD" if ( $d->[30] && $d->[30] < $MINVD );
    my $class = $d->[3] =~ /COSM/ ? "COSMIC" : ($d->[3] =~ /^rs/ ? (checkCLNSIG($d->[31]) == 1 ? "ClnSNP" : "dbSNP") : "Novel");

    # Rescue deleterious dbSNP, such as rs80357372 (BRCA1 Q139* that is in dbSNP, but not in ClnSNP or COSMIC
    if ( ($d->[6] =~ /STOP_GAINED/i || $d->[6] =~ /FRAME_?SHIFT/i) && $class eq "dbSNP" ) {
	my $pos = $1 if ( $d->[10] =~ /(\d+)/ );
	$class = "dbSNP_del" if ( $pos/$d->[12] < 0.95 );
    }

    if ( $d->[6] =~ /SPLICE/i && $class eq "dbSNP" ) {
        $class = "dbSNP_del";
    }

    #$class = "dbSNP" if ( $d->[28] && $d->[28] > $MAF ); # if there's MAF with frequency, it'll be considered dbSNP regardless of COSMIC
    if ( $d->[29] ) {
        $d->[29] =~ s/^\[//; $d->[29] =~ s/\]$//;
	my @mafs = split(/,/, $d->[29]);
	if ( @mafs == 2 && $mafs[1] ne "." && $mafs[1] > $MAF ) {
	    $class = "dbSNP";
	} elsif ( @mafs > 2 ) {  # For dbSNP with multiple alleles in one position
	    my $mk = $1 if ( $d->[3] =~ /(rs\d+)/ );
	    $mk .= "-$d->[4]-$d->[5]";
	    $class = "dbSNP" if ($MultiMaf{ $mk } && $MultiMaf{ $mk } > $MAF );
	}
    }
    $pass = "CNTL" if ( $CONTROL{ $vark } );
    $pass = "BIAS" if ( $opt_b && ($class eq "Novel"||$class eq "dbSNP") && ($d->[23] eq "2;1" || $d->[23] eq "2;0") && $d->[22] < 0.3 ); # Filter novel variants with strand bias.
    $pass = "NonClnSNP" if ( checkCLNSIG($d->[31]) == -1 && $class ne "COSMIC" );
    $pass = "AMPBIAS" if ( $opt_a && $d->[40] && $d->[40] < $d->[41] );
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

sub mean {
    my $ref = shift;
    my ($sum, $n) = (0, 0);
    foreach( @$ref ) {
        $sum += $_;
	$n++;
    }
    return sprintf("%.3f", $sum/$n);
}

sub setupMultiMaf {
    my $in = shift;
    open(MMAF, $in);
    while( <MMAF> ) {
	chomp;
        my @a = split;
	$MultiMaf{ "$a[2]-$a[3]-$a[4]" } = $a[5];
    }
    close( MMAF );
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
    -s If set, it'll keep SNPEff's amino acid change as is.  Default: it'll change three letter code to one
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
        The filtering mean mapping quality score for variants.  The raw variant will be filtered if the mean mapping quality score is less then 
	specified unless the allele frequency is greater than 0.8.
	Default: 10

    -D INT
        The filtering total depth.  The raw variant will be filtered on first place if the total depth is less then INT.  Default: 0

    -V INT
        The filtering variant depth.  Variants with depth < INT will be considered false positive.  Default: 2, or at least 2 reads are needed for a variant

    -o signal
        The signal/noise value.  Default: 1.5 

    -c Control(s)
        The control sample name(s).  Any novel variants passing all above filters but also detected in Control sample will be deemed considered
	false positive.  Use only when there's control sample.  Multiple controls samples are separated using ":", e.g. s1:s2:s3.

    -C additional_columns
        Add additional columns in VCF to be appended to the output.  Use : to separate multiple columns.  Only those defined in VCF are allowed.
	
    -G DOUBLE
        The mininum GMAF value.  Any variants with GMAF above this value is deemed dbSNP, regardless whether it's in COSMIC or not.  Default: 0.0025
	
    -A file
        A file that contain GMAF when there're multiple alternative alleles.  It's not easy to be parsed from CAF as the order is not clear.
	Thus this extra file.  Use only if you have it available.  It should contain 6 columns, such as "chr1    907920  rs28430926      C       G       0.1107",
	where the last column is GMAF.  Default to: /ngs/reference_data/genomes/Hsapiens/hg19/variation/dbSNP_multi_mafs_latest.txt.  Use "", empty
	string to disable it if you don't have one.  If the default file doesn't exist, it'll be disabled.

    A novel variant (non-dbSNP, non-COSMIC) is considered false positive if all three conditions (-r -f -n) are met. Any variant meeting the -p
    or -q conditions are also considered likely false positive.  False positive variants are annotated "FALSE" in column PASS, "TRUE" otherwise.

AUTHOR
     Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
     Report bugs to zhongwu\@yahoo.com

COPYRIGHT
     This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

USAGE
    exit(0);
}
