#!/usr/bin/perl -w

#use Getopt::Std;
use Getopt::Long qw(:config no_ignore_case);
use strict;

our ($opt_n, $opt_f, $opt_F, $opt_H, $opt_D, $opt_V, $opt_M, $opt_R, $opt_p, $opt_N, $opt_r, $opt_a, $opt_s, $opt_S);

#my $ruledir = "/ngs/reference_data/genomes/Hsapiens/hg19/variation/rules";  # "/users/kdld047/work/NGS/Wee1";
my $ruledir = "/users/kdld047/work/NGS/Wee1";
my $annotation_dir = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/cancer_informatics';  # "/ngs/reference_data/genomes/Hsapiens/hg19/variation/cancer_informatics";
my $filter_common_snp = "$annotation_dir/filter_common_snp.txt";
my $snpeffect_export_polymorphic = "$annotation_dir/snpeffect_export_Polymorphic.txt";
my $filter_common_artifacts = "$annotation_dir/filter_common_artifacts.txt";
my $actionable_hotspot = "$annotation_dir/actionable_hotspot.txt";
my $actionable = "$annotation_dir/actionable.txt";
my $compendia_ms7_hotspot = "$annotation_dir/Compendia.MS7.Hotspot.txt";
my $SPLICE = "$annotation_dir/SPLICE.txt";
my $LASTAA = "$annotation_dir/last_critical_aa.txt";
my $GMAF = 0.0025;

GetOptions ("n=s" => \$opt_n,
            "f=f" => \$opt_f,
            "F=f" => \$opt_F,
            "H|?" => \$opt_H,
            "D=i" => \$opt_D,
            "V=i" => \$opt_V,
            "M"   => \$opt_M,
            "S"   => \$opt_S,
            "N"   => \$opt_N,
            "r"   => \$opt_r,
            "R=f" => \$opt_R,
            "p=s" => \$opt_p,
            "a=s" => \$opt_a,
            "s=i" => \$opt_s,
            "GMAF=f" => \$GMAF,
            "ruledir=s" => \$ruledir,
            "annotdir=s" => \$opt_a,
            "filter_common_snp=s" => \$filter_common_snp,
            "snpeffect_export_polymorphic=s" => \$snpeffect_export_polymorphic,
            "filter_common_artifacts=s" => \$filter_common_artifacts,
            "actionable_hotspot=s" => \$actionable_hotspot,
            "actionable=s" => \$actionable,
            "compendia_ms7_hotspot=s" => \$compendia_ms7_hotspot) or USAGE();
if ( $opt_a ) {
    $annotation_dir = $opt_a eq "hg38" ? "/ngs/reference_data/genomes/Hsapiens/hg38/variation/cancer_informatics" : $opt_a;
    $filter_common_snp = "$annotation_dir/filter_common_snp.txt";
    $snpeffect_export_polymorphic = "$annotation_dir/snpeffect_export_Polymorphic.txt";
    $filter_common_artifacts = "$annotation_dir/filter_common_artifacts.txt";
    $actionable_hotspot = "$annotation_dir/actionable_hotspot.txt";
    $actionable = "$annotation_dir/actionable.txt";
    $compendia_ms7_hotspot = "$annotation_dir/Compendia.MS7.Hotspot.txt";
    $SPLICE = "$annotation_dir/SPLICE.txt";
}

#getopts( 'HMn:f:F:D:V:R:' ) || USAGE();
$opt_H && USAGE();

my $tp53group1 = parseMut_tp53( "$ruledir/Rules/DNE.txt" );
my $tp53group2 = parseMut_tp53( "$ruledir/Rules/TA0-25.txt" );
my $tp53group3 = parseMut_tp53( "$ruledir/Rules/TA25-50_SOM_10x.txt" );
my %tp53critical_aa_pos;
my %splice = ();
if ( -e $SPLICE ) {
    open(SP, $SPLICE);
    while( <SP> ) {
        chomp;
        my ($p, $g) = split(/\t/);
        $splice{ $g }->{ $p } = 1;
    }
    close( SP );
}

my %lastaa = ();
if ( -e $LASTAA ) {
    open(AA, $LASTAA);
    while( <AA> ) {
        chomp;
        my ($g, $ap, $acc) = split(/\t/);
        $lastaa{ $g } = $ap;
    }
    close( AA );
}

# Set up common SNP filter
my $MINAF = $opt_f ? $opt_f : 0.05;
my $MAXRATE = $opt_R ? $opt_R : 1.00;
my $MINCOMSAMPLE = $opt_s ? $opt_s : 5; # the minimum number of samples sharing the same variant
my $ACTMINAF = $opt_F ? $opt_F : ($MINAF/2 < 0.01 ? $MINAF/2 : 0.01);
my %filter_snp;
open( FSNP, $filter_common_snp );
while( <FSNP> ) {
    chomp;
    my @a = split(/\t/);
    my $key = join("-", @a[1..4]);
    $filter_snp{ $key } = 1;
}
close( FSNP );

my %snpeff_snp;
open( SESNP, $snpeffect_export_polymorphic );
while( <SESNP> ) {
    chomp;
    my @a = split(/\t/);
    my $key = $a[11] ? join("-", @a[11,2]) : $a[5];
    next unless($key);
    $snpeff_snp{ $key } = 1 unless( $key eq "-" );
}
close( SESNP );

# Set up artifact filter
my %filter_art;
open( ART, $filter_common_artifacts );
while( <ART> ) {
    chomp;
    my @a = split(/\t/);
    my $key = join("-", @a[1..4]);
    $filter_art{ $key } = 1;
}
close( ART );

# Set up actionable protein mutation
my %act_hot;
open( ACTHOT, $actionable_hotspot );
while( <ACTHOT> ) {
    chomp;
    my ($gene, $pchg, $sg) = split(/\t/);
    $act_hot{ "$gene-$pchg" } = 1;
}
close( ACTHOT );

# Set up actionable mutations
my %act_som;
my %act_germ;
my %rules;
open( ACT, $actionable );
while( <ACT> ) {
    chomp;
    my @a = split(/\t/);
    if ( $a[7] eq "germline" ) {
        my $key = join("-", @a[1..4]);
        $act_germ{ $key } = "germline";
    } elsif ( $a[7] eq "somatic" ) {
        if ( $a[6] eq "rule" ) {
            if ( $a[4] eq "*" && length($a[3]) == 1 ) {
                my $key = join("-", @a[1..3]);
                $act_som{ $key } = 1;
            } else {
                push(@{ $rules{ $a[5] }->{ $a[0] } }, [@a[1,2,3,4,8]]);
            }
        } else {
            my $key = join("-", @a[1..4]);
            $act_som{ $key } = $a[8] ? $a[8] : 1;
        }
    }
}
close( ACT );

# Set up Compendia Hotspot
my %hotspotnt;
my %hotspotprot;
open( HOT, $compendia_ms7_hotspot );
while( <HOT> ) {
    chomp;
    my @a = split(/\t/);
    next if ( $a[5] =~ /^g./ );
    my $key = join("-", @a[1..4]);
    $hotspotnt{ $key } = 1;
    next unless( $a[6] );
    $hotspotprot{ "$a[0]-$a[6]" } = 1;
}
close( HOT );

$opt_n = $opt_n ? qr/$opt_n/ : "";
my $hdr = <>; chomp $hdr;
my @hdrs = split(/\t/, $hdr);
my %hdrs;
for(my $i = 0; $i < @hdrs; $i++) {
    $hdrs{ $hdrs[$i] } = $i;
}
my $passcol = $hdrs{ PASS };
my $classcol = $hdrs{ Var_Class };
my $typecol = $hdrs{ Type };
my $funccol = $hdrs{ Functional_Class };
my $genecodecol = $hdrs{ Gene_Coding };
my $afcol = $hdrs{ AlleleFreq };
my $genecol = $hdrs{ Gene };
my $aachgcol = $hdrs{ Amino_Acid_Change };
my $cosmaachgcol = $hdrs{ COSMIC_AA_Change };
my $msicol = $hdrs{ MSI };
my $statuscol = $hdrs{ Status };
if ( $opt_M ) {
    print join("\t", "SAMPLE ID", "ANALYSIS FILE LOCATION", "VARIANT-TYPE", "GENE", "SOMATIC STATUS/FUNCTIONAL IMPACT", qw(SV-PROTEIN-CHANGE SV-CDS-CHANGE SV-GENOME-POSITION SV-COVERAGE SV-PERCENT-READS CNA-COPY-NUMBER CNA-EXONS CNA-RATIO CNA-TYPE REARR-GENE1 REARR-GENE2 REARR-DESCRIPTION REARR-IN-FRAME? REARR-POS1 REARR-POS2 REARR-NUMBER-OF-READS));
    print $statuscol ? "\tStatus\n" : "\n";
} else {
    print "$hdr\tStatus\n";
}
while( <> ) {
    chomp;
    my @a = split(/\t/);
    next unless( $a[$passcol] eq "TRUE" );
    my $sample = $a[0];
    my $chr = $a[1];
    $chr = "chr$chr" unless( $chr =~ /^chr/ );
    my $key = join("-", $chr, @a[2,4,5]);
    my $af = $a[$afcol];
    my $act = isActionable( $chr, @a[2,4,5], $a[$hdrs{Gene}], $a[$aachgcol], $a[$cosmaachgcol], \@a );
    unless( $act ) {
        next if ( $filter_snp{ $key } );
        next if ( $snpeff_snp{ "$a[$hdrs{Gene}]-$a[$hdrs{Amino_Acid_Change}]" } );
        next if ( $filter_art{ $key } && $af < 0.20 );
        my @gmaf = split(/,/, $a[$hdrs{GMAF}]);
        my $flag = @gmaf ? 0 : 1;
        foreach my $maf (@gmaf) {
            $flag = 1 if ( $maf && $maf <= $GMAF );
        }
        next unless( $flag );
    }
    next if ( $opt_D && $a[$hdrs{Depth}] < $opt_D );
    if ( $opt_V ) {
        if ( $a[$hdrs{VD}] < $opt_V ) {
            next unless( $af >= 0.5 );
        }
    }
    my @snps = $a[3] =~ /(rs\d+)/g;
    foreach my $rs (@snps) {
        next if ( $snpeff_snp{ $rs } );
    }
    my $aachg = $a[$aachgcol];
    my $gene = $a[$genecol];
    my $platform = $sample =~ /[_-]([^_\d]+?)$/ ? $1 : "";
    $platform = "" unless( $platform =~ /^WXS/i || $platform =~ /^RNA-Seq/i || $platform =~ /^VALIDATION/i || $platform =~ /^WGS/ );
    $platform = $opt_p if ( $opt_p );

    if ( $opt_n && $sample =~ /$opt_n/ ) {
        $sample = $1;
    } elsif ( $opt_n && $sample !~ /$opt_n/ ) {
        next;
    }
    my ($type, $fclass, $gene_coding) = @a[$typecol, $funccol, $genecodecol];

    # Filter low AF MSI
    if ( abs(length($a[4])-length($a[5])) == 1 ) {
        my $msi = $a[$msicol];
        if ( $msi <= 7 ) {
            next if ( $af < 0.05 );
        } elsif ( $msi == 8 ) {
            next if ( $af < 0.07 );
        } elsif ( $msi == 9 ) {
            next if ( $af < 0.125 );
        } elsif ( $msi == 10 ) {
            next if ( $af < 0.175 );
        } elsif ( $msi == 11 ) {
            next if ( $af < 0.25 );
        } elsif ( $msi == 12 ) {
            next if ( $af < 0.3 );
        } elsif ( $msi > 12 ) {
            next if ( $af < 0.35 );
        }
    }
    my $status = "unknown";
    if ( $a[$classcol] eq "ClnSNP" ) {
        $status = "likely";
    } elsif ( $a[$classcol] eq "dbSNP_del" ) {
        $status = "likely";
    } elsif ( $a[$classcol] eq "ClnSNP_known" ) {
        $status = "known";
    } elsif ( $a[$classcol] eq "ClnSNP_unknown" ) {
        $status = "unknown";
    } elsif ( $a[6] =~ /FRAME_?SHIFT/i ) {
        $status = "likely";
    } elsif ( $aachg =~ /^[A-Z]+\d+\*$/ ) {
        $status = "likely";
    }
    if ( $type =~ /splice/i && ($type =~ /acceptor/i || $type =~ /donor/i) ) {
        $status = "likely";
        $a[10] = "splice";
    } elsif ( length($a[10]) == 0 && ($type =~ /SPLICE/i && $type !~ /region_variant/) ) {
        if ( $a[$hdrs{cDNA_Change}] ) {
            if ($a[$hdrs{cDNA_Change}] =~ /\d+\+(\d+)/) {
                if ( $1 <= 2 ) {
                    $status = "likely";
                    $a[10] = "splice";
                }
            } elsif ( $a[$hdrs{cDNA_Change}] =~ /\d+-(\d+)[^_]\S+$/ ) {
                if ( $1 <= 2 ) {
                    $status = "likely";
                    $a[10] = "splice";
                }
            }
        } else {  # No cDNA_Change, for earlier version compatibility
            $status = "likely";
            $a[10] = "splice";
        }
    }
    if ( $a[$classcol] eq "COSMIC" ) {
        if ( $hdrs{ COSMIC_Cnt } ) {
            my @cnts = split(/,/, $a[$hdrs{ COSMIC_Cnt }]);
            foreach my $cnt (@cnts) {
                if ( $cnt >= 5 ) {
                    $status = "likely";
                }
            }
        }
        if ( $a[$cosmaachgcol] ) {
            $a[$cosmaachgcol] =~ s/^p\.//;
        }
    }
    if ( isHotspotNT($chr, @a[2,4,5]) ) {
        $status = "likely";
        #print STDERR "$a[6] $a[11] @a[2,4,5] $status\n";
    } elsif ( isHotspotProt( $a[$hdrs{Gene}], $a[$hdrs{Amino_Acid_Change}] ) ) {
        $status = "likely";
    }
    if ($act_som{ $key } || $act_germ{ $key }) {
        $status = "known";
    }
    if ( $act ) {
        $status = "known";
        next if ( $af < $ACTMINAF );
        #next if ( $af < 0.2 && $act eq "germline" );
    } else {
        #next if ( $type =~ /^INTRON/i || $type =~ /^SYNONYMOUS_/i || $fclass eq "SILENT" || ($type =~ /splice_region_variant/ && $a[10] eq "") );
        if ( $type =~ /^SYNONYMOUS_/i || $fclass eq "SILENT" ) {
            next if ( $a[$classcol] eq "dbSNP" || $a[$hdrs{ID}] =~ /rs/ );
            next if( $opt_S );
        }
        if ( ($type =~ /^INTRON/i || ($type =~ /splice/i && $a[10] eq "")) && $status eq "unknown" ) {
            if ( $opt_N ) {
                $status = "unknown";
            } else {
                next;
            }
        }
        next if ( $af < $MINAF );
    }
    #next if ( $status ne "known" && ($type =~ /^UPSTREAM/i || $type =~ /^DOWNSTREAM/i || $type =~ /^INTERGENIC/i || $type =~ /^INTRAGENIC/i || ($type =~ /UTR_/ && $type !~ /codon/i ) || $gene_coding =~ /NON_CODING/i || $fclass =~ /^NON_CODING/i) );
    next if ( $status ne "known" && ($type =~ /^UPSTREAM/i || $type =~ /^DOWNSTREAM/i || $type =~ /^INTERGENIC/i || $type =~ /^INTRAGENIC/i || $gene_coding =~ /NON_CODING/i || $fclass =~ /^NON_CODING/i) );
    if ( $status ne "known" && ($type =~ /UTR_/ && $type !~ /codon/i ) ) {
        if ( $opt_N ) {
            $status = "unknown";
        } else {
            next;
        }
    }
    next if ( $a[$classcol] eq "dbSNP" && $status ne "known" );

    # Ignore any variants that occur after last known critical amino acid
    if ( $aachg =~ /^[A-Z](\d+)/ ) {
        next if ( $lastaa{ $gene } && $1 >= $lastaa{ $gene } );
    }

    if ( $opt_R && $status ne "known" && $a[$hdrs{ Pcnt_sample }] > $MAXRATE && $a[$hdrs{ N_Var }] >= $MINCOMSAMPLE ) {
        if ( $opt_r ) {
            #print "$_\t$status\n";
            print join("\t", $sample, $platform, "short-variant", $gene, $status, $a[10], $a[$hdrs{cDNA_Change}], "$chr:$a[2]", $a[$hdrs{Depth}], $af*100, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-");
            print $statuscol ? "\t$a[$statuscol]\n" : "\n";
        } else {
            next;
        }
    }
    next if ( $opt_r );
    if ( $opt_M ) {
        print join("\t", $sample, $platform, "short-variant", $gene, $status, $a[10], $a[$hdrs{cDNA_Change}], "$chr:$a[2]", $a[$hdrs{Depth}], $af*100, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-");
        print $statuscol ? "\t$a[$statuscol]\n" : "\n";
    } else {
        print "$_\t$status\n";
    }
}

sub isActionable {
    my ($chr, $pos, $ref, $alt, $gene, $aachg, $cosmaachg, $ra) = @_;
    my $key = join("-", $chr, $pos, $ref, $alt);
    return $act_som{ $key } if ($act_som{ $key });
    return $act_germ{ $key } if ($act_germ{ $key });
    if ( length($ref) == 1 && length($ref) == length($alt) ) {
        return $act_som{ "$chr-$pos-$ref" } if ( $act_som{ "$chr-$pos-$ref" } );
    }
    return $act_hot{ "$gene-$aachg" } if ( $act_hot{ "$gene-$aachg" } && $aachg =~ /^([A-Z]\d+)[A-Z]$/ );
    if ( $aachg =~ /^([A-Z]\d+)[A-Z]$/ ) {
        return $act_hot{ "$gene-$1" } if ( $act_hot{ "$gene-$1" } );
    }
    if ( $gene eq "TP53" ) {
        my $tp53_group = classify_tp53($aachg, $pos, $ref, $alt);
        return "somatic" unless( $tp53_group eq "NA" );
    }
    if ( $rules{ "inframe-del" }->{ $gene } && length($ref) > length($alt) && (length($ref)-length($alt))%3 == 0 ) {
        foreach my $r ( @{ $rules{ "inframe-del" }->{ $gene } } ) {
            if ( $r->[0] eq $chr && $r->[1] <= $pos && $r->[2] >= $pos && (length($ref)-length($alt)) >= $r->[3] ) {
                $ra->[$hdrs{Amino_Acid_Change}] = $r->[4];
                return "somatic";
            }
        }
    } elsif ( $rules{ "inframe-ins" }->{ $gene } && length($ref) < length($alt) && (length($alt)-length($ref))%3 == 0 ) {
        foreach my $r ( @{ $rules{ "inframe-ins" }->{ $gene } } ) {
            if ( $r->[0] eq $chr && $r->[1] <= $pos && $r->[2] >= $pos && (length($alt)-length($ref)) >= $r->[3] ) {
                $ra->[$hdrs{Amino_Acid_Change}] = $r->[4];
                return "somatic";
            }
        }
    } elsif ( $rules{ "indel" }->{ $gene } && length($ref) != length($alt) ) {
        foreach my $r ( @{ $rules{ "indel" }->{ $gene } } ) {
            if ( $r->[0] eq $chr && $r->[1] <= $pos && $r->[2] >= $pos && abs(length($alt)-length($ref)) >= $r->[3]) {
                $ra->[$hdrs{Amino_Acid_Change}] = $r->[4];
                return "somatic";
            }
        }
    } elsif ( $rules{ "del" }->{ $gene } && length($ref) > length($alt) ) {
        foreach my $r ( @{ $rules{ "del" }->{ $gene } } ) {
            if ( $r->[0] eq $chr && $r->[1] <= $pos && $r->[2] >= $pos && (length($ref)-length($alt)) >= $r->[3]) {
                $ra->[$hdrs{Amino_Acid_Change}] = $r->[4];
                return "somatic";
            }
        }
    } elsif ( $rules{ "ins" }->{ $gene } && length($ref) < length($alt) ) {
        foreach my $r ( @{ $rules{ "ins" }->{ $gene } } ) {
            if ( $r->[0] eq $chr && $r->[1] <= $pos && $r->[2] >= $pos && (length($alt)-length($ref)) >= $r->[3]) {
                $ra->[$hdrs{Amino_Acid_Change}] = $r->[4];
                return "somatic";
            }
        }
    }
    return 0;
}

sub isHotspotNT {
    my ($chr, $pos, $ref, $alt) = @_;
    if ( length($ref) > length($alt) && $alt ne "-" ) {
        $ref = substr($ref, 1);
        $alt = length($alt) > 1 ? substr($alt, 1) : "-";
    } elsif ( length($alt) > length($ref) && $ref ne "-" ) {
        $alt = substr($alt, 1);
        $ref = length($ref) > 1 ? substr($ref, 1) : "-";
    }
    return $hotspotnt{ "$chr-$pos-$ref-$alt" } ? 1 : 0;
}

sub isHotspotProt {
    my ($gene, $pchg) = @_;
    $pchg =~ s/^p.//;
    return 0 unless( $pchg );
    return $hotspotprot{ "$gene-$pchg" } ? 1 : 0;
}

sub classify_tp53 {
    my ($aa, $pos, $ref, $alt) = @_;
    $ref =~ s/\s+//g;
    $alt =~ s/\s+//g;
    $pos =~ s/\s+//g;
    $aa =~ s/\s+//g;
    if ( $splice{ TP53 }->{ $pos } && length($ref) == 1 && length($alt) == 1 ) {
        return "Group 6";
    }
    $aa =~ s/^p.//;
    if ($aa =~ /^[A-Z]\d+[A-Z]$/ ) {
        return "Group 1" if ( $tp53group1->{ $aa } );
        return "Group 2" if ( $tp53group2->{ $aa } );
        return "Group 3" if ( $tp53group3->{ $aa } );
        return "NA";
    } elsif ( $aa =~ /^[A-Z](\d+)\*$/ ) {
        return $1 < 359 ? "Group 4" : "NA";
    } elsif ( $aa =~ /^[A-Z](\d+)fs/ ) {
        return $1 < 359 ? "Group 5" : "NA";
    }
    return "NA";
}

sub parseMut_tp53 {
    my $list = shift;
    open(MUT, $list);
    my %hash;
    while(<MUT>) {
        chomp;
        my @a = split(/\t/);
        next unless( $a[19] =~ /^p./ );
        $a[19] =~ s/^p.//;
        $hash{ $a[19] } = 1;
        $tp53critical_aa_pos{ $1 } = 1 if ( $a[19] =~ /^[A-Z](\d+)[A-Z]$/ );
    }
    close(MUT);
    return \%hash;
}
sub USAGE {
    print STDERR <<USAGE;
    Usage: $0 [-n reg_name] input
    The program will filter the VarDict output after vcf2txt.pl to candidate interpretable mutations, somatic or germline.

    Options:
    -H  Print this usage
    -M  Output in FM's format
    -S        Skip the silent mutation that are "unknown" in Status.  By default: all silent "Novel" and "COSMIC" only mutations will be kept.
        Silent mutation with entries in both COSMIC and dbSNP will be filtered.

    -N  If set, keep all intronic and UTR in the output, but will be set as "unknown".
    -D  int
        The minimum total depth
    -V  int
        The minimum reads supporting variant
    -n  reg
        The regular expression to extract sample name.  Default: Use as it is.  For TCGA, "(TCGA-..-....)-01" is preferred
        for tumor sample.

    -f  double
        The minimum allele frequency for regular variants. Default: 0.05

    -R  double
        If a variant is present in > [double] fraction of samples, it's deemed not a mutation.  Default: 1.0, or no filtering.
        Use with caution.  It'll filter even if it's in COSMIC, unless it's on actionable list. Don't use it if the cohort is homogeneous.
        Use only for heterogeneous cohorts.  Used together with -s option.

    -s  integer
        The minimum number of samples that sharing the same variant that are not on the known list to be filtered out.  Used together
        with -R option.  A variant has to satisfied both conditions, AND not on the knon list, to be considered too common to
        be functional and filtered out.  Use only for heterogeneous cohorts. Default: 5

    -r  If set, keep only those variants satisfying -R option.  The option is meant to find what are re-occuring variants or
        artifacts.

    -F  double
        The minimum allele frequency hotspot somatic mutations, typically lower then -f.  Default: 0.01 or half -f, whichever is less

    -p  string
        The platform, such as WXS, WGS, RNA-Seq, VALIDATION, etc.  No Default.  Used when output is in FM's format (-M option)

    -a  dirpath
        Directory path for containing annotation files.
        Default is /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation
        For hg38, use /group/cancer_informatics/tools_resources/NGS/genomes/hg38/Annotation
        If dirpath is "hg38", automatically set it to /group/cancer_informatics/tools_resources/NGS/genomes/hg38/Annotation

    --GMAF double
        When the GMAF is greater than specified, it's considered common SNP and filtered out.  Default: 0.0025 (0.25%)

    --ruledir  dirpath
        default is /users/kdld047/work/NGS/Wee1

    --annotdir  dirpath
        default is /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation

    --filter_common_snp  filepath
        default is /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/filter_common_snp.txt

    --snpeffect_export_polymorphic  filepath
        default is /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/snpeffect_export_Polymorphic.txt

    --filter_common_artifacts  filepath
        default is /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/filter_common_artifacts.txt

    --actionable_hotspot  filepath
        default is /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/actionable_hotspot.txt

    --actionable  filepath
        default is /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/actionable.txt

    --compendia_ms7_hotspot  filepath
        default is /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/Compendia.MS7.Hotspot.txt

USAGE
    exit(0);
}
