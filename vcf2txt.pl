#!/usr/bin/env perl

use Getopt::Std;
use warnings;
use strict;

our ($opt_F, $opt_R, $opt_s, $opt_r, $opt_n, $opt_N, $opt_H, $opt_q, $opt_p, $opt_b, $opt_f, $opt_c, $opt_u, $opt_D, $opt_Q, $opt_P, $opt_M, $opt_o, $opt_V, $opt_C, $opt_G, $opt_A, $opt_L);
getopts('suHabR:F:f:n:r:p:q:c:D:P:Q:M:o:V:C:G:A:L') || USAGE();

my %AA_code = (
    "ALA" => "A",  "ILE" => "I",  "LEU" => "L",  "VAL" => "V",
    "PHE" => "F",  "TRP" => "W",  "TYR" => "Y",  "ASN" => "N",
    "CYS" => "C",  "GLN" => "Q",  "MET" => "M",  "SER" => "S",
    "THR" => "T",  "ASP" => "D",  "GLU" => "E",  "ARG" => "R",
    "HIS" => "H",  "LYS" => "K",  "GLY" => "G",  "PRO" => "P", "TER" => "*"
);

USAGE() if ( $opt_H );
my $FRACTION = $opt_r ? $opt_r : 0.4;
my $MAXRATIO = $opt_R ? $opt_R : 1.0;
my $CNT = $opt_n ? $opt_n : 10;
my $AVEFREQ = $opt_F ? $opt_F : 0.15;
my $MINPMEAN = $opt_p ? $opt_p : 5;
my $MINQMEAN = $opt_q ? $opt_q : 25;
my $FILPMEAN = $opt_P ? $opt_P : 0; # will be filtered on the first place
my $FILQMEAN = $opt_Q ? $opt_Q : 0; # will be filtered on the first place
my $FILDEPTH = $opt_D ? $opt_D : 0; # will be filtered on the first place
my $MINFREQ = $opt_f ? $opt_f : 0.02;
my $MINMQ = $opt_M ? $opt_M : 10;
my $MINVD = $opt_V ? $opt_V : 2; # minimum variant depth
my $MAF = $opt_G ? $opt_G : 0.0025;
my $SN = $opt_o ? $opt_o : 1.5;
my $PRINTLOF = $opt_L;
my @controls = $opt_c ? split(/:/, $opt_c) : ();
my %controls = map { ($_, 1); } @controls;
my %MultiMaf;
$opt_A = defined($opt_A) ? $opt_A : "/ngs/reference_data/genomes/Hsapiens/hg19/variation/dbSNP_multi_mafs_latest.txt";
setupMultiMaf($opt_A) if ( -e $opt_A );


# SBF: Strand Bias Fisher Exact test
my @columns = qw(CDS AA CNT END DP AF BIAS PMEAN PSTD QUAL QSTD SBF CAF VD CLNSIG ODDRATIO HIAF MQ SN ADJAF NM SHIFT3 MSI dbSNPBuildID);
my @ampcols = ();
my @paircols = ();
my @appcols = $opt_C ? split(/:/, $opt_C) : ();
push(@columns, @appcols) if ($opt_C);
my @data;
my %sample;
my %var;
my %CONTROL;
while( <> ) {
    next if ( /^#/ );
    chomp;
    my @a = split(/\t/);
    next unless( $a[6] eq "PASS" );
    $a[7] .= ";";
    my %d;
    while( $a[7] =~ /([^=;]+)=([^=]+?);/g ) {
	$d{ $1 } = $2;
    }
    @ampcols = qw(GDAMP TLAMP NCAMP AMPFLAG) if ( $d{GDAMP} && $d{TLAMP} );
    my @formats = split(/:/, $a[8]);
    my @fdata = split(/:/, $a[9]);
    my @mfdata = $a[10] ? split(/:/, $a[10]) : ();
    for(my $i = 0; $i < @formats; $i++) {
	$d{ $formats[$i] } = $fdata[$i];
	$d{ "M_$formats[$i]" } = $mfdata[$i] if ( $a[10] );
    }

    # Adapt for Mutect or FreeBayes
    unless( $d{ PMEAN } ) {  # Meaning not VarDict
    	delete $d{ AF };
    	# Use AD in Mutect
		if ( $d{ AD } ) {
			my @ads = split(/,/, $d{ AD });
			my $ads_sum = 0;
			$ads_sum += $_ foreach( @ads );
			if ($ads_sum > 0) {
				$d{ AF } = sprintf("%.3f", $ads[1]/$ads_sum);
			} else {
				$d{ AF } = '0';
			}
			$d{ VD } = $ads[1];
		}
		# In FreeBayes, we don't trust AD, but rather use AO and RO for allele freq calculation
		if ( $d{ AO } ) {
			my $ao_sum = 0;
			my @aos = split(/,/, $d{ AO });
			@aos = sort { $b <=> $a } @aos; # Just make sure the first is the most frequency
			$ao_sum += $_ foreach( @aos );
			$d{ AF } = sprintf("%.3f", $aos[0]/($ao_sum+$d{RO}));
			$d{ VD } = $aos[0];
		}
		next unless( $d{ AF } ); # No AD or AO found in Mutect/Freebayes - skipping this variant

		if ( $a[10] ) { # for somatic paired analysis
			if ( $d{ M_AD } ) {
			my @m_ads = split(/,/, $d{ M_AD });
			my $m_ads_sum = 0;
			foreach( @m_ads ) {
				$m_ads_sum += $_ if ( /\d/ );
			}
			$d{ M_AF } = $m_ads_sum ? sprintf("%.3f", $m_ads[1]/$m_ads_sum) : 0;
			$d{ M_VD } = $m_ads[1] && $m_ads[1] =~ /\d/ ? $m_ads[1] : 0;
			}
			# Use AO and RO for allele freq calculation for FreeBayes and overwrite AD even if it exists
			if ( $d{ M_AO } ) {
			my $m_ao_sum = 0;
			my @m_aos = split(/,/, $d{ M_AO });
			@m_aos = sort { $b <=> $a } @m_aos if ( $m_aos[0] =~ /\d/ ); # Just make sure the first is the most frequency
			foreach( @m_aos ) {
				$m_ao_sum += $_ if ( /\d/ );
			}
			my $m_ro = $d{M_RO} && $d{M_RO} =~ /\d/ ? $d{M_RO} : 0;
			$d{ M_AF } = $m_ao_sum + $m_ro > 0 ? sprintf("%.3f", $m_aos[0]/($m_ao_sum+$m_ro)) : 0;
			$d{ M_VD } = $m_aos[0] && $m_aos[0] =~ /\d/ ? $m_aos[0] : 0;
			}
		}
    }

    @paircols = qw(TYPE STATUS SSF SOR M_DP M_AF M_VD M_BIAS M_PMEAN M_PSTD M_QUAL M_QSTD M_HIAF M_MQ M_SN M_ADJAF M_NM) if ( $a[10] );

    $d{ SBF } = $d{ SBF } < 0.0001 ? sprintf("%.1e", $d{ SBF }) : sprintf("%.4f", $d{ SBF }) if ( $d{ SBF } );
    $d{ ODDRATIO } = sprintf("%.3f", $d{ ODDRATIO }) if ( $d{ ODDRATIO } );
    my @effs = ();
    if ( $d{ EFF } ) {
        @effs = split(/,/, $d{ EFF });
    } elsif ( $d{ ANN } ) {
        my @anns = split(/,/, $d{ ANN });
	foreach my $ann (@anns) { 
	    my ($allele, $ann, $impact, $gname, $gid, $ftype, $fid, $biotype, $rank, $hgvsc, $hgvsp, $cdnap, $cdsp, $protp, $dist) = split(/\|/, $ann);
	    my @ta = split(/&/, $ann);
	    $ann = $ta[0];
	    $protp = $1 if ( $protp =~ /\d+\/(\d+)/ );
	    #my @alts = split(/,/, $a[4]);
	    push(@effs, "$ann($impact|$ann|$hgvsc|$hgvsp/$hgvsc|$protp|$gname|$biotype|$ftype|$fid|$rank|1)");
	}
    } else {
        @effs = (" (||||||||||1)");
    }
    my $vark = join(":", @a[0,1,3,4]); # Chr Pos Ref Alt
    next if ( $FILDEPTH && $d{ DP } < $FILDEPTH );
    next if ( $FILPMEAN && $d{ PMEAN } && $d{ PMEAN } < $FILPMEAN );
    next if ( $FILQMEAN && $d{ QUAL } && $d{ QUAL } < $FILQMEAN );
    my ($pmean, $qmean) = ($d{ PMEAN }, $d{ QUAL });
    my $pass = "TRUE";
    #$pass = "FALSE" unless ( $d{PSTD} > 0 );
    $pass = "FALSE" if ( $qmean && $qmean < $MINQMEAN );
    $pass = "FALSE" if ( $pmean && $pmean < $MINPMEAN );
    $pass = "FALSE" if ( !$d{AF} || $d{AF} < $MINFREQ );
    $pass = "FALSE" if ( $d{MQ} && $d{MQ} < $MINMQ && $d{AF} < 0.5 );  # Keep low mapping quality but high allele frequency variants
    $pass = "FALSE" if ( $d{SN} && $d{SN} < $SN );
    $pass = "FALSE" if ( !$d{VD} || $d{VD} < $MINVD );
    if ( $d{ SAMPLE } && $controls{ $d{ SAMPLE } } ) {
	my $clncheck = checkCLNSIG($d{CLNSIG});
	my $class = $a[2] =~ /COSM/ ? "COSMIC" : ($a[2] =~ /^rs/ ? ($clncheck ? $clncheck : "dbSNP") : "Novel");
	$CONTROL{ $vark } = 1 if ( $pass eq "TRUE" && $class eq "Novel");  # so that any novel variants showed up in control won't be filtered
    }
    unless( $opt_u && $d{ SAMPLE } =~ /Undetermined/i ) { # Undetermined won't count toward samples
	$sample{ $d{ SAMPLE } } = 1;
	push( @{ $var{ $vark } }, $d{ AF } ) if ( $pass eq "TRUE" );
    }
    my @alts = split(/,/, $a[4]);
    for(my $i = 0; $i < @effs; $i++) {
        my $eff = $effs[$i];
	$eff =~ s/\)$//;
	my @e = split(/\|/, $eff, -1);

	my ($type, $effect) = split(/\(/, $e[0]);
	my @tmp = map { defined($d{ $_ }) ? $d{ $_ } : ""; } (@columns, @ampcols, @paircols);
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
		for(my $i = 0; $i < length($5); $i += 3) {
		    $ins .= $AA_code{ uc(substr($5, $i, 3)) };
		}
		$aachg = "$AA_code{uc($1)}${2}_$AA_code{uc($3)}${4}ins$ins";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(\d+)(_.*)?fs$/ ) {
		$aachg = "$AA_code{uc($1)}${2}fs";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(\d+)del$/ ) {
		$aachg = "$AA_code{uc($1)}${2}del";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(\D*)?(\d+)([\*\?])$/ ) {
		my $aa = $AA_code{uc($1)};
		my $ppos = $3;
		if ( $2 ) {
		    for(my $i = 0; $i < length($2); $i += 3) {
			$aa = $AA_code{ uc(substr($2, $i, 3)) };
			$ppos++;
		    }
		}
		$aachg = "$aa$ppos$4";
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
		$aachg .= "_" . $AA_code{ uc($4) } . $5 if ( $4 );
		if ( $6 ) {
		    for(my $i = 0; $i < length($6); $i += 3) {
			$insaa .= $AA_code{ uc(substr($6, $i, 3)) } ? $AA_code{ uc(substr($6, $i, 3)) } : "?";
		    }
		}
		$aachg .= $insaa ? "delins$insaa" : "del";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(\d+)(_([A-Z][a-z][a-z])(\d+))?delins\?+$/ ) {
		$aachg = $AA_code{ uc($1) } . $2;
		$aachg .= "_" . $AA_code{ uc($4) } . $5 if ( $4 );
		$aachg .= "?";
	    } elsif ( $aachg =~ /^([A-Z][a-z][a-z])(\d+)(_([A-Z][a-z][a-z])(\d+))?dup$/ ) {
		$aachg = $AA_code{ uc($1) } . $2;
		$aachg .= "_" . $AA_code{ uc($4) } . $5 if ( $4 );
		$aachg .= "dup";
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
	push(@data, [$d{ SAMPLE }, @a[0..3], $alts[$e[10]-1], $type, $effect, @tmp2, @tmp, $d{ LOF }]);
    }
}

my @amphdrs = @ampcols > 0 ? qw(GAmplicons TAmplicons NCAmplicons Ampflag) : ();
my @pairhdrs = @paircols > 0 ? qw(VType Status Paired-p_value Paired-OddRatio Matched_Depth Matched_AlleleFreq Matched_VD Matched_Bias Matched_Pmean Matched_Pstd Matched_Qual Matched_Qstd Matched_HIAF Matched_MQ Matched_SN Matched_AdjAF Matched_NM)  : ();
print join("\t", qw(Sample Chr Start ID Ref Alt Type Effect Functional_Class Codon_Change Amino_Acid_Change cDNA_Change Amino_Acid_Length Gene Transcript_bioType Gene_Coding Transcript Exon COSMIC_CDS_Change COSMIC_AA_Change COSMIC_Cnt End Depth AlleleFreq Bias Pmean Pstd Qual Qstd SBF GMAF VD CLNSIG ODDRATIO HIAF MQ SN AdjAF NM Shift3 MSI dbSNPBuildID), @appcols, @amphdrs, @pairhdrs, qw(N_samples N_Var Pcnt_sample Ave_AF PASS Var_Type Var_Class));
if ($PRINTLOF) {
	print "\tLOF";
} 
print "\n";

my @samples = keys %sample;
my $sam_n = @samples + 0;
foreach my $d (@data) {
    my $vark = join(":", @$d[1, 2, 4, 5]); # Chr Pos Ref Alt
    next unless( $var{ $vark } ); # Likely just in Undetermined.
    my $type = length($d->[4]) == length($d->[5]) ? (length($d->[4]) == 1 ? "SNV" : (length($d->[4]) <= 3 ? "MNV" : "Complex" )) : (substr($d->[4], 0, 1) ne substr($d->[5], 0, 1) ? "Complex" : (length($d->[4]) > length($d->[5]) ? "Deletion" : "Insertion" ));
    my ($af, $pmean, $qmean, $mq, $sn) = @$d[23, 25, 27, 35, 36];
    my $varn = @{ $var{ $vark } } + 0;
    my $ave_af = mean( $var{ $vark } );
    my $pass = ($varn/$sam_n > $FRACTION && $varn >= $CNT && $ave_af < $AVEFREQ && $d->[3] eq ".") ? "MULTI" : "TRUE"; # novel and present in $FRACTION samples
    my $clncheck = checkCLNSIG($d->[32]);
    my $class = $d->[3] =~ /COSM/ ? "COSMIC" : ($d->[3] =~ /^rs/ ? ($clncheck ? $clncheck : "dbSNP") : "Novel");
    #$pass = "FALSE" unless ( $d->[24] > 0 ); # all variants from one position in reads
    $pass = "DUP" if ( $pmean && $d->[26] ==  0 && $d->[24] !~ /1$/ && $d->[24] !~ /0$/ && (@amphdrs == 0) && $af < 0.35 ); # all variants from one position in reads
    $pass = "QMEAN" if (length($qmean) > 0 && $qmean < $MINQMEAN );
    $pass = "PMEAN" if ($pmean && $pmean < $MINPMEAN );
    $pass = "MQ" if ( length($mq) > 0 && $mq < $MINMQ && $af < 0.8 ); # Keep low mapping quality but high allele frequency variants
    $pass = "SN" if ( length($sn) > 0 && $sn < $SN );
    $pass = "MINFREQ" if ( $af < $MINFREQ );
    $pass = "MINVD" if ( $d->[31] && $d->[31] < $MINVD );

    # Rescue deleterious dbSNP, such as rs80357372 (BRCA1 Q139* that is in dbSNP, but not in ClnSNP or COSMIC
    if ( ($d->[6] =~ /STOP_GAINED/i || $d->[6] =~ /FRAME_?SHIFT/i) && $class eq "dbSNP" ) {
	my $pos = $1 if ( $d->[10] =~ /(\d+)/ );
	$class = "dbSNP_del" if ( $pos/$d->[12] < 0.95 );
    }

    # Consider splice variants deleterious
    if ( $d->[6] =~ /SPLICE/i && $d->[6] !~ /region/i && $class eq "dbSNP" ) {
	$class = "dbSNP_del";
    }

    #$class = "dbSNP" if ( $d->[28] && $d->[28] > $MAF ); # if there's MAF with frequency, it'll be considered dbSNP regardless of COSMIC
    if ( $d->[30] ) {  # GMAF
	$d->[30] =~ s/^\[//; $d->[30] = (split /\]/, $d->[30])[0];
	my @mafs = split(/,/, $d->[30]);
	if ( @mafs == 2 ) {
	    $class = "dbSNP" if ( $mafs[1] ne "." && $mafs[1] > $MAF );
	    $d->[30] = $mafs[1] ne "." ? $mafs[1] : "";
	} elsif ( @mafs > 2 ) {  # For dbSNP with multiple alleles in one position
	    my $mk = $d->[3] =~ /(rs\d+)/ ? $1 : "";
	    $mk .= "-$d->[4]-$d->[5]";
	    if ($MultiMaf{ $mk }) {
		$class = "dbSNP" if( $MultiMaf{ $mk } > $MAF );
		$d->[30] = $MultiMaf{ $mk };
	    #} else {
	        #print STDERR "$d->[3] with CAF $d->[29] not recorded\n";
		#$d->[29] = "";
	    } else {
		my @tmafs = ();
		foreach( @mafs ) {
		    push(@tmafs, $_) unless( $_ eq "." );
		}
		$d->[30] = join(",", @tmafs[1..$#tmafs]);
	    }
	}
    }
    $pass = "CNTL" if ( $CONTROL{ $vark } );
    $pass = "BIAS" if ( $opt_b && ($class eq "Novel"||$class eq "dbSNP") && ($d->[24] eq "2;1" || $d->[24] eq "2;0") && $d->[23] < 0.3 ); # Filter novel variants with strand bias.
    if ( $clncheck eq "dbSNP" && $class ne "COSMIC" && $class ne "dbSNP_del" ) {
	$class = "dbSNP"; 
    }
    $pass = "AMPBIAS" if ( @amphdrs > 0 && $d->[42+@appcols] && $d->[42+@appcols] < $d->[43+@appcols] );
    if ( $opt_R && $pass eq "TRUE" && $varn/$sam_n > $MAXRATIO && $varn > $CNT ) { # present in $MAXRATIO samples, regardless of frequency
        if ( max( $var{ $vark } ) > 0.35 ) { 
	    $class = "dbSNP";
	} else {
	    $pass = "MAXRATE" if ( $af < 0.35 );
	}
    }
    my $lof = pop(@$d);
    print join("\t", @$d, $sam_n, $varn, sprintf("%.3f", $varn/$sam_n), $ave_af, $pass, $type, $class) if ( $pass eq "TRUE" );
	if ( $PRINTLOF ) {
    	my $is_lof = "";
    	my $effect = @$d[7];
    	$is_lof = "YES" if ( $lof && $effect eq "HIGH" && index($lof, @$d[13]) != -1 );
    	print "\t$is_lof";
	}
	print "\n";
}

sub checkCLNSIG {
    my $clnsig = shift;
    return 0 unless( $clnsig );
    my @cs = split(/\||,/, $clnsig );
    my $flag255 = 0;
    my $flagno = 0;
    foreach my $cs (@cs) {
	return "ClnSNP_known" if ( $cs > 3 && $cs < 7 );
	$flagno++ if ( $cs <= 3 && $cs >= 2 );
	$flag255++ if ( $cs == 255 || $cs < 1 );
    }
    return "dbSNP" if ( $flagno > 1 && $flagno >= $flag255 );
    return "ClnSNP_unknown" if ( $flag255 ); # Keep unknown significant variants
    return "dbSNP";
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

sub max {
    my $ref = shift;
    my $max = 0;
    foreach( @$ref ) {
	$max = $_ if ( $max < $_ );
    }
    return $max;
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
	When a novel variant is present in more than [fraction] of samples and mean allele frequency is less than -F, it's 
	considered as likely false positive. Default 0.4.
	Used with -F and -n
    
    -F DOUBLE
	When the ave allele frequency is also below -F, the variant is considered likely false positive.  Default 0.15.
	Used with -r and -n

    -n INT
	When the variant is detected in greater or equal [sample_cnt] samples, the variant is considered likely false positive.  Default 10.
	Used with -r and -F

    -R DOUBLE
	When a passing variant is present in more than [fraction] of samples and at least -n samples , it's considered as 
	dbSNP, even if it's in COSMIC or apparent deleterious. Default 1.0. or no filtering.  Use with caution.  Don't use it for homogeneous 
	samples.  Use only for hetereogeneous samples, such as 0.5, or any variants present in 50% of samples are considered
	as dbSNP.

    -f DOUBLE
	When individual allele frequency < -f for variants, it was considered likely false poitives. Default: 0.02 or 2%

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

    A novel variant (non-dbSNP, non-COSMIC) is considered false positive if all three conditions (-r -F -n) are met. Any variant meeting the -p
    or -q conditions are also considered likely false positive.  False positive variants are annotated "FALSE" in column PASS, "TRUE" otherwise.

	-L 
	Report SNPEff's LOF annotation (will add additional column named LOF).

AUTHOR
     Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
     Report bugs to zhongwu\@yahoo.com

COPYRIGHT
     This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

USAGE
    exit(0);
}
