#!/usr/bin/perl -w

#use Getopt::Std;
use lib '/users/kdld047/lib/perl5';
use Getopt::Long qw(:config no_ignore_case);
use Stat::Basic;
use strict;

our ($opt_n, $opt_f, $opt_F, $opt_H, $opt_D, $opt_V, $opt_M, $opt_R, $opt_p, $opt_N, $opt_r, $opt_a, $opt_s, $opt_S, $opt_O, $opt_y, $opt_m, $opt_B);

#my $ruledir = "/ngs/reference_data/genomes/Hsapiens/hg19/variation/rules";  # "/users/kdld047/work/NGS/Wee1";
my $ruledir = "/users/kdld047/work/NGS/Wee1/Rules";
my $annotation_dir = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/cancer_informatics';  # "/ngs/reference_data/genomes/Hsapiens/hg19/variation/cancer_informatics";
my $filter_common_snp = "$annotation_dir/filter_common_snp.txt";
my $snpeffect_export_polymorphic = "$annotation_dir/snpeffect_export_Polymorphic.txt";
my $filter_common_artifacts = "$annotation_dir/filter_common_artifacts.txt";
my $actionable_hotspot = "$annotation_dir/actionable_hotspot.txt";
my $actionable = "$annotation_dir/actionable.txt";
my $compendia_ms7_hotspot = "$annotation_dir/Compendia.MS7.Hotspot.txt";
my $SPLICE = "$annotation_dir/SPLICE.txt";
my $LASTAA = "$annotation_dir/last_critical_aa.txt";
my $blackgenes = "$annotation_dir/blackgenes.txt";
my $GMAF = 0.0025;
my $recalfreq = '';
my @data = ();
my %samples;
my %mutcnt;
my %mut2sam;
my $stat = new Stat::Basic;

GetOptions ("n=s" => \$opt_n,
            "f=f" => \$opt_f,
            "F=f" => \$opt_F,
            "H|?" => \$opt_H,
            "D=i" => \$opt_D,
            "V=i" => \$opt_V,
            "M"   => \$opt_M,
            "O"   => \$opt_O,
            "S"   => \$opt_S,
            "N"   => \$opt_N,
            "r"   => \$opt_r,
            "y"   => \$opt_y,
            "R=f" => \$opt_R,
            "B=f" => \$opt_B,
            "p=s" => \$opt_p,
            "a=s" => \$opt_a,
            "s=i" => \$opt_s,
            "recalfreq" => \$recalfreq,
            "GMAF=f" => \$GMAF,
            "ruledir=s" => \$ruledir,
            "annotdir=s" => \$opt_a,
            "filter_common_snp=s" => \$filter_common_snp,
            "snpeffect_export_polymorphic=s" => \$snpeffect_export_polymorphic,
            "filter_common_artifacts=s" => \$filter_common_artifacts,
            "actionable_hotspot=s" => \$actionable_hotspot,
            "actionable=s" => \$actionable,
            "lastaa=s" => \$LASTAA,
            "blackgenes=s" => \$blackgenes,
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
    $LASTAA = "$annotation_dir/last_critical_aa.txt";
    $blackgenes = "$annotation_dir/blackgenes.txt";
}

#getopts( 'HMn:f:F:D:V:R:' ) || USAGE();
$opt_H && USAGE();

my $platform = $opt_p if ( $opt_p );
my $tp53group1 = parseMut_tp53( "$ruledir/DNE.txt" );
my $tp53group2 = parseMut_tp53( "$ruledir/TA0-25.txt" );
my $tp53group3 = parseMut_tp53( "$ruledir/TA25-50_SOM_10x.txt" );
my %tp53critical_aa_pos;
my %splice = ();
if ( -e $SPLICE ) {
    open(SP, $SPLICE);
    while( <SP> ) {
        chomp;
        next if ( /^##/ );
        my ($p, $g) = split(/\t/);
        $splice{ $g }->{ $p } = 1;
    }
    close( SP );
} else {
    print STDERR "Warning: No splice junction file '$SPLICE' supplied\nProceed...\n";
}

my %lastaa = ();
if ( -e $LASTAA ) {
    open(AA, $LASTAA);
    while( <AA> ) {
        chomp;
        next if ( /^##/ );
        my ($g, $ap, $acc) = split(/\t/);
        $lastaa{ $g } = $ap;
    }
    close( AA );
} else {
    print STDERR "Warning: No last known significant amino acid file '$LASTAA' supplied.\nProceed...\n";
}

my %blackgenes = ();
if ( -e $blackgenes ) {
    open(BLACK, $blackgenes);
    while( <BLACK> ) {
        chomp;
        next if ( /^##/ );
        my ($g, $reason) = split(/\t/, $_, 2);
        $blackgenes{ $g } = $reason;
    }
    close( BLACK );
} else {
    print STDERR "Warning: No last known significant amino acid file '$blackgenes' supplied.\nProceed...\n";
}

# Set up common SNP filter
my $MINAF = $opt_f ? $opt_f : 0.05;
my $MAXRATE = $opt_R ? $opt_R : 1.00; 
my $MINCOMSAMPLE = $opt_s ? $opt_s : 5; # the minimum number of samples sharing the same variant
my $ACTMINAF = $opt_F ? $opt_F : ($MINAF/2 < 0.01 ? $MINAF/2 : 0.01);
my %filter_snp;
if ( -e $filter_common_snp ) {
    open( FSNP, $filter_common_snp );
    while( <FSNP> ) {
        chomp;
        next if ( /^##/ );
        my @a = split(/\t/);
        my $key = join("-", @a[1..4]);
        $filter_snp{ $key } = 1;
    }
    close( FSNP );
} else {
    print STDERR "Warning: No curated common SNP file '$filter_common_snp' supplied.\nProceed...\n";
}

my %snpeff_snp;
if ( -e $snpeffect_export_polymorphic ) {
    open( SESNP, $snpeffect_export_polymorphic );
    while( <SESNP> ) {
        chomp;
        next if ( /^##/ );
        my @a = split(/\t/);
        my $key = $a[11] ? join("-", @a[11,2]) : $a[5];
        next unless($key);
        $snpeff_snp{ $key } = 1 unless( $key eq "-" );
    }
    close( SESNP );
} else {
    print STDERR "Warning: No SNPEffect polymorphic database file '$snpeffect_export_polymorphic' supplied!\nProceed...\n";
}

# Set up artifact filter
my %filter_art;
my %filter_rules;
if ( -e $filter_common_artifacts ) {
    open( ART, $filter_common_artifacts );
    while( <ART> ) {
        chomp;
        next if ( /^##/ );
        my @a = split(/\t/);
        if ( $a[5] eq "rule" ) {
            push(@{ $filter_rules{ $a[0] }->{ $a[4] } }, [@a[1,2,3]]);
        } else {
            my $key = join("-", @a[1..4]);
            $filter_art{ $key } = 1;
        }
    }
    close( ART );
} else {
    print STDERR "Warning: No common artifacts file '$filter_common_artifacts' supplied!\nProceed...\n";
}

# Set up actionable protein mutation
my %act_hot;
my %comm_snp;
if ( -e $actionable_hotspot ) {
    open( ACTHOT, $actionable_hotspot );
    while( <ACTHOT> ) {
        chomp;
        next if ( /^##/ ); # comment line
        my ($gene, $pchg, $sg) = split(/\t/);
        if ( $gene =~ /^#/ ) { # VUS, No special treatment for now
            $gene =~ s/^#//;
        } elsif ( $gene =~ /^\^/ ) { # common variant
            $gene =~ s/^\^//;
            $comm_snp{ "$gene-$pchg" } = 1;
        } else { # actionable
            $act_hot{ "$gene-$pchg" } = $sg eq "somatic" ? "somatic" : "germline";
        }
    }
    close( ACTHOT );
} else {
    print STDERR "Warning: No actionable hotspot file '$actionable_hotspot' supplied.\nProceed...\n";
}

# Set up actionable mutations
my %act_som;
my %act_germ;
my %rules;
if ( -e $actionable ) {
    open( ACT, $actionable );
    while( <ACT> ) {
        chomp;
        next if ( /^##/ );
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
} else {
    print STDERR "Warning: No actionable file '$actionable' supplied.\nProceed...\n";
}

# Set up Compendia Hotspot
my %hotspotnt;
my %hotspotprot;
if ( -e $compendia_ms7_hotspot ) {
    open( HOT, $compendia_ms7_hotspot );
    while( <HOT> ) {
        chomp;
        next if ( /^##/ );
        my @a = split(/\t/);
        next if ( $a[5] =~ /^g./ );
        my $key = join("-", @a[1..4]);
        $hotspotnt{ $key } = 1;
        next unless( $a[6] );
        $hotspotprot{ "$a[0]-$a[6]" } = 1;
    }
    close( HOT );
} else {
    print STDERR "Warning: No Compendia hotspot file '$compendia_ms7_hotspot' supplied.\nProceed...\n";
}

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
my $endcol = $hdrs{ End };
my $statuscol = $hdrs{ Status };
if ( $opt_M ) {
    print join("\t", "SAMPLE ID", "ANALYSIS FILE LOCATION", "VARIANT-TYPE", "GENE", "SOMATIC STATUS/FUNCTIONAL IMPACT", qw(SV-PROTEIN-CHANGE SV-CDS-CHANGE SV-GENOME-POSITION SV-COVERAGE SV-PERCENT-READS CNA-COPY-NUMBER CNA-EXONS CNA-RATIO CNA-TYPE REARR-GENE1 REARR-GENE2 REARR-DESCRIPTION REARR-IN-FRAME? REARR-POS1 REARR-POS2 REARR-NUMBER-OF-READS));
    print $statuscol ? "\tStatus\tSGZ\n" : "\tSGZ\n";
} else {
    print "$hdr\tStatus\tSGZ\n";
}
while( <> ) {
    chomp;
    my @a = split(/\t/, $_, -1);
    next unless( $a[$passcol] eq "TRUE" );
    my $sample = $a[0];
    my $chr = $a[1];
    $chr = "chr$chr" unless( $chr =~ /^chr/ );
    my $key = join("-", $chr, @a[2,4,5]);
    my $af = $a[$afcol];
    my $act = isActionable( $chr, @a[2,4,5], $a[$hdrs{Gene}], $a[$aachgcol], $a[$cosmaachgcol], \@a, $af, $a[$hdrs{CLNSIG}] );
    my $sgz = $act;
    if ( filterRule( $chr, @a[2,4,5], $a[$hdrs{Gene}], $a[$aachgcol], $a[$cosmaachgcol], \@a, $a[$hdrs{CLNSIG}] ) ) {
        print STDERR "$chr @a[2,4,5] $a[$hdrs{Gene}] filtered by rules.\n" if ( $opt_y );
        next;
    }
    if ( $comm_snp{ "$a[$hdrs{Gene}]-$a[$hdrs{Amino_Acid_Change}]" } ) {
        print STDERR "Filtered as it's curated as common SNP: $a[$hdrs{Gene}]-$a[$hdrs{Amino_Acid_Change}].\n" if ( $opt_y );
        next;
    }

    unless( $act ) {
        if ( $filter_snp{ $key } ) {
            print STDERR "Filtered as common SNP $key\n" if ( $opt_y );
            next;
        }
        if ( $filter_art{ $key } && $af < 0.35 ) {
            print STDERR "Filtered as likely artifact $key AF: $af < 0.35\n" if ( $opt_y );
            next;
        }
        my @gmaf = split(/,/, $a[$hdrs{GMAF}]);
        my $flag = @gmaf ? 0 : 1;
        foreach my $maf (@gmaf) {
            $sgz = "germline" unless( $sgz );
            $flag = 1 if ( $maf && $maf =~ /\d/ && $maf <= $GMAF );
        }
        unless( $flag ) {
            print STDERR "Filtered as GMAF is higher than $GMAF $key @gmaf\n" if ( $opt_y );
            next;
        }
        next if ( $a[$classcol] eq "dbSNP" );
        next if ( $snpeff_snp{ "$a[$hdrs{Gene}]-$a[$hdrs{Amino_Acid_Change}]" } && $a[$classcol] ne "ClnSNP_known" );
    }
    next if ( $opt_D && $a[$hdrs{Depth}] < $opt_D );

    # Filter strand biased variants
    if ( $opt_B && $a[$hdrs{Bias}] ne "2:2" ) {
        if ( $a[$hdrs{Bias}] eq "2:1" && $a[$hdrs{SBF}] < $opt_B * 10) {
            print STDERR "Filtered as it has only single strand support '2:1' while reference has both $key\n" if ( $opt_y );
            next;
        }
        if ( $a[$hdrs{Bias}] eq "2:0" && $a[$hdrs{SBF}] < $opt_B ) {
            print STDERR "Filtered as it has only single strand support '2:0' and p-value < $opt_B while reference has both $key\n" if ( $opt_y );
            next;
        }
        # Filter variants with opposite supporting strands for REF and ALT
        if ( $a[$hdrs{Bias}] eq "1:1" && $a[$hdrs{SBF}] < $opt_B/10 ) {
            print STDERR "Filtered as it has only single strand support '1:1' and p-value < $opt_B/10 while reference has both $key\n" if ( $opt_y );
            next;
        }

        if ( $hdrs{ ALD } ) {
            my ($frd, $rrd) = split(/,/, $a[$hdrs{ ALD }]);
            if ( $a[$hdrs{Bias}] =~ /^2/ && $frd * $rrd == 0 && ($frd + $rrd >= 4)) {
                print STDERR "Filtered as it has only single strand support ($frd, $rrd) while reference has both $key\n" if ( $opt_y );
                next;
            }
        }
    }

    if ( $opt_V ) {
        if ( $a[$hdrs{VD}] < $opt_V ) {
            next unless( $af >= 0.5 );
        }
    }
    my @snps = $a[3] =~ /(rs\d+)/g;
    my $gene = $a[$genecol];
    foreach my $rs (@snps) {
        if ( $snpeff_snp{ $rs } ) {
            print STDERR "Filtered as in SnpEff database. $rs $gene $key\n" if ( $opt_y );
            next;
        }
    }
    my $aachg = $a[$aachgcol];

    # Filter genes in black list
    if ( $blackgenes{ $gene } ) {
        unless( $act ) {
            print STDERR "Filtered due to $gene is on black gene list as $blackgenes{ $gene }.\n" if ( $opt_y );
            next;
        }
    }

    # Ignore HLA genes unless -O is specified
    unless( $opt_O ) {
        next if ( $gene =~ /^HLA-/ );
    }

    next if ( $gene =~ /^OR\d+[A-Z]/ ); # Ignore olfactory genes

    unless( $opt_p ) {
        $platform = $sample =~ /[_-]([^_\d]+?)$/ ? $1 : "";
        $platform = "" unless( $platform =~ /^WXS/i || $platform =~ /^RNA-Seq/i || $platform =~ /^VALIDATION/i || $platform =~ /^WGS/ );
    }

    if ( $opt_n && $sample =~ /$opt_n/ ) {
        $sample = $1;
    } elsif ( $opt_n && $sample !~ /$opt_n/ ) {
        next;
    }
    $samples{ $sample } = 1;
    my ($type, $fclass, $gene_coding) = @a[$typecol, $funccol, $genecodecol];
    my $cdna = $a[$hdrs{ cDNA_Change }];
    if ( $type =~ /upstream/i || $type =~ /downstream/i ) {
        next unless( $act );
    }
    if ( $type =~ /feature_ablation/i || $type =~ /transcript_ablation/i ) {
        next unless( $act );
    }
    if ( $cdna =~ /^c\.-\d+_\*/ ) {
        next unless( $act );
    }
    if ( $cdna =~ /^n\./ ) {
        next unless( $act || $type =~ /fusion/i );
    }

    # Filter low AF MSI
    my $msi = $a[$msicol];
    if ( abs(length($a[4])-length($a[5])) == 1 && $msi > 1 ) {  
        if ( $msi <= 2 ) {
            next if ( $af < 0.005 );
        } if ( $msi <= 4 ) {
            next if ( $af < 0.01 );
        } if ( $msi <= 7 ) {
            next if ( $af < 0.03 );
        } elsif ( $msi == 8 ) {
            next if ( $af < 0.06 );
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
    } elsif ( abs(length($a[4])-length($a[5])) == 3 && $msi >= 5 ) {
        next if ( $af < 0.10 ); # ignore low AF in 3nt MSI region
    }
    my $status = "unknown";
    if ( $a[$classcol] eq "ClnSNP" ) {
        $status = "likely";
    } elsif ( $a[$classcol] eq "dbSNP_del" ) {
        $status = "likely";
    } elsif ( $a[$classcol] eq "ClnSNP_known" ) {
        $status = "known";
    } elsif ( $a[6] =~ /FRAME_?SHIFT/i ) {
        $status = "likely";
    } elsif ( $aachg =~ /^[A-Z]+\d+\*$/ ) {
        $status = "likely";
    } elsif ( $aachg =~ /\*$/ || $aachg =~ /ins.*\*[A-Z]+$/ ) {
        $status = "likely";
    } elsif ( $a[$classcol] eq "ClnSNP_unknown" ) {
        $status = "unknown";
    } 
    if ( $type =~ /splice/i && ($type =~ /acceptor/i || $type =~ /donor/i) ) {
        my $spflag = 1;
        if ( $cdna =~ /[\+-](\d+)_-?\d+[\+-](\d+)/ ) {
            $spflag = 0 if ( $1 > 2 && $2 > 2 );
        } elsif ( $cdna =~ /^.\.\d+[\+-](\d+)[a-z]/ ) {
            $spflag = 0;
        }
        if ( $spflag ) {
            $status = "likely";
            $a[10] = "splice";
        }
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
                    $status = "likely" unless( $status eq "known" );
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
        #next if ( $af < 0.025 && $gene eq "TP53" ); # As TP53 is early event, it should be a little more stringent
        next if ( $af < 0.15 && $act eq "germline" );
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
        if ( $type =~ /^intron_variant/i && $status ne "known" ) {
            if ( $opt_N ) {
                $status = "unknown";
            } else {
                next;
            }
        }
        if ( $type =~ /^sequence_feature/i && $status ne "known" ) {
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
            print join("\t", $sample, $platform, "short-variant", $gene, $status, $a[10], $a[$hdrs{cDNA_Change}], "$chr:$a[2]", $a[$hdrs{Depth}], sprintf("%.2f", $af*100), "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-");
            print $statuscol ? "\t$a[$statuscol]\n" : "\n";
        } else {
            next;
        }
    }
    next if ( $opt_r && (!$recalfreq) );
    $sgz = "NA" unless( $sgz );
    if ( $recalfreq ) {
        push(@data, [@a, $status, "$chr-$a[2]-$a[4]-$a[5]", $sgz]);
        $mut2sam{ "$chr-$a[2]-$a[4]-$a[5]" }->{ $sample } = $af;
    } else {
        if ( $opt_M ) {
            unless( $type =~ /fusion/i ) {
                print join("\t", $sample, $platform, "short-variant", $gene, $status, $a[10], $hdrs{cDNA_Change} && $a[$hdrs{cDNA_Change}] ? $a[$hdrs{cDNA_Change}] : "", "$chr:$a[2]", $a[$hdrs{Depth}], sprintf("%.2f", $af*100), $a[3], "-", "-", "-", "-", "-", "-", "-", "-", "-", "-");
            } else {
                my ($g1, $g2) = split(/\&/, $gene, 2);
                $status = "likelY" if ( $status eq "unknown" );
                if ( $a[$hdrs{cDNA_Change}] && $a[$hdrs{cDNA_Change}] =~ /(\d+)_(\d+)/ ) {
                    ($g1, $g2) = ($g2, $g1) if ( $1 > $2 );
                }
                print join("\t", $sample, $platform, "fusion", $g1, $status, "$g1-$g2", $hdrs{cDNA_Change} && $a[$hdrs{cDNA_Change}] ? $a[$hdrs{cDNA_Change}] : "", "$chr:$a[2]", $a[$hdrs{Depth}], sprintf("%.2f", $af*100), $a[3], "-", "-", "-", $g1, $g2, "fusion", "-", "$chr:$a[2]", "$chr:$a[$endcol]", "-");
            }
            print $statuscol ? "\t$a[$statuscol]\t$sgz\n" : "\t$sgz\n";
        } else {
            print "$_\t$status\t$sgz\n";
        }
    }
}

if ( $recalfreq ) {
    my $samcnt = (keys %samples) + 0;
    while( my ($k, $r) = each %mut2sam ) {
        $mutcnt{ $k }->{ cnt } = (keys %$r) + 0;
        my @tmp = values %$r;
        $mutcnt{ $k }->{ pt75 } = $stat->prctile(\@tmp, 75);
        $mutcnt{ $k }->{ median } = $stat->median(\@tmp);
    }
    foreach my $d (@data) {
        my $sgz = pop( @$d );
        my $k = pop( @$d );
        $d->[$hdrs{ N_Var }] = $mutcnt{$k}->{ cnt };
        $d->[$hdrs{ N_samples }] = $samcnt;
        $d->[$hdrs{ Pcnt_sample }] = sprintf("%.2f", $mutcnt{ $k }->{ cnt }/$samcnt);
        $d->[$hdrs{ Ave_AF }] = $mutcnt{$k}->{ median };
        my $type = $d->[$hdrs{ Type }];
        my $status = pop( @$d );
        my $sample = $d->[0];
        $sample = $1 if ( $opt_n && $sample =~ /$opt_n/ );
        my $gene = $d->[$genecol];
        if ( $opt_R && $status ne "known" && $d->[$hdrs{ Pcnt_sample }] > $MAXRATE && $d->[$hdrs{ N_Var }] >= $MINCOMSAMPLE ) {
            if ( $opt_r ) {
                #print "$_\t$status\n";
                print join("\t", $sample, $platform, "short-variant", $gene, $status, $d->[10], $d->[$hdrs{cDNA_Change}], "$d->[1]:$d->[2]", $d->[$hdrs{Depth}], sprintf("%.2f", $d->[$afcol]*100), "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"), "\n";
                #print $statuscol ? "\t$a[$statuscol]\n" : "\n";
            } else {
                next;
            }
        }
        next if ( $opt_r );
        if ( $opt_M ) {
            unless( $type =~ /fusion/i ) {
                print join("\t", $sample, $platform, "short-variant", $gene, $status, $d->[10], $d->[$hdrs{cDNA_Change}], "$d->[1]:$d->[2]", $d->[$hdrs{Depth}], sprintf("%.2f", $d->[$afcol]*100), "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", $sgz), "\n";
            } else {
                my ($g1, $g2) = split(/\&/, $gene, 2);
                $status = "likelY" if ( $status eq "unknown" );
                if ( $d->[$hdrs{cDNA_Change}] && $d->[$hdrs{cDNA_Change}] =~ /(\d+)_(\d+)/ ) {
                    ($g1, $g2) = ($g2, $g1) if ( $1 > $2 );
                }
                print join("\t", $sample, $platform, "rearrangement", $g1, $status, "$g1-$g2", $d->[$hdrs{cDNA_Change}], "$d->[1]:$d->[2]", $d->[$hdrs{Depth}], sprintf("%.2f", $d->[$afcol]*100), $d->[3], "-", "-", "-", $g1, $g2, "fusion", "-", "$d->[1]:$d->[1]", "$d->[1]:$d->[$endcol]", "-", $sgz);
            }
        } else {
            print join("\t", @$d, $status, $sgz), "\n";
        }
    }
}

sub filterRule {
    my ($chr, $pos, $ref, $alt, $gene, $aachg, $cosmaachg, $ra, $clnsig) = @_;
    return 0 unless( $filter_rules{ $gene } );
    if ( $filter_rules{ $gene }->{ ignore } ) {
        foreach my $r (@{ $filter_rules{ $gene }->{ ignore } } ) {
            return 1 if ( $chr eq $r->[0] && $pos >= $r->[1] && $pos <= $r->[2] );
        }
    }
    return 0;
}

sub isActionable {
    my ($chr, $pos, $ref, $alt, $gene, $aachg, $cosmaachg, $ra, $af, $clnsig) = @_;
    my $key = join("-", $chr, $pos, $ref, $alt);
    return $act_som{ $key } if ($act_som{ $key });
    return $act_germ{ $key } if ($act_germ{ $key } && $af >= 0.15); # actionable germline need higher af
    if ( length($ref) == 1 && length($ref) == length($alt) ) {
        return $act_som{ "$chr-$pos-$ref" } if ( $act_som{ "$chr-$pos-$ref" } );
    }
    if ( $act_hot{ "$gene-$aachg" } && $aachg =~ /^([A-Z]\d+)[A-Z?]$/ ) {
        if ( $act_hot{ "$gene-$aachg" } eq "somatic" ) {
            return "somatic" if ( $af >= $ACTMINAF );
        } elsif ( $act_hot{ "$gene-$aachg" } eq "germline" ) {
            return "germline" if ( $af >= 0.15 );
        }
    }
    if ( $aachg =~ /^([A-Z]\d+)[A-Z]$/ || $aachg =~ /^(M1)\?$/ ) {
        return $act_hot{ "$gene-$1" } if ( $act_hot{ "$gene-$1" } );
    }
    if ( $gene eq "TP53" ) {
        my $tp53_group = classify_tp53($aachg, $pos, $ref, $alt);
        return "somatic" unless( $tp53_group eq "NA" );
    }
    if ( $rules{ "inframe-del" }->{ $gene } && length($ref) > length($alt) && (length($ref)-length($alt))%3 == 0 ) {
        foreach my $r ( @{ $rules{ "inframe-del" }->{ $gene } } ) {
            if ( $r->[0] eq $chr && $r->[1] <= $pos+length($ref)-1 && $r->[2] >= $pos && (length($ref)-length($alt)) >= $r->[3] ) {
                $ra->[$hdrs{Amino_Acid_Change}] = $r->[4] if ( $opt_M );
                return "somatic"; 
            }
        }
    } elsif ( $rules{ "inframe-ins" }->{ $gene } && length($ref) < length($alt) && (length($alt)-length($ref))%3 == 0 ) {
        foreach my $r ( @{ $rules{ "inframe-ins" }->{ $gene } } ) {
            if ( $r->[0] eq $chr && $r->[1] <= $pos+length($ref)-1 && $r->[2] >= $pos && (length($alt)-length($ref)) >= $r->[3] ) {
                $ra->[$hdrs{Amino_Acid_Change}] = $r->[4] if ( $opt_M );
                return "somatic";
            }
        }
    } elsif ( $rules{ "indel" }->{ $gene } && length($ref) != length($alt) ) {
        foreach my $r ( @{ $rules{ "indel" }->{ $gene } } ) {
            if ( $r->[0] eq $chr && $r->[1] <= $pos+length($ref)-1 && $r->[2] >= $pos && abs(length($alt)-length($ref)) >= $r->[3]) {
                $ra->[$hdrs{Amino_Acid_Change}] = $r->[4] if ( $opt_M );
                return "somatic"; 
            }
        }
    } elsif ( $rules{ "del" }->{ $gene } && length($ref) > length($alt) ) {
        foreach my $r ( @{ $rules{ "del" }->{ $gene } } ) {
            if ( $r->[0] eq $chr && $r->[1] <= $pos+length($ref)-1 && $r->[2] >= $pos && (length($ref)-length($alt)) >= $r->[3]) {
                $ra->[$hdrs{Amino_Acid_Change}] = $r->[4] if ( $opt_M );
                return "somatic";
            }
        }
    } elsif ( $rules{ "ins" }->{ $gene } && length($ref) < length($alt) ) {
        foreach my $r ( @{ $rules{ "ins" }->{ $gene } } ) {
            if ( $r->[0] eq $chr && $r->[1] <= $pos+length($ref)-1 && $r->[2] >= $pos && (length($alt)-length($ref)) >= $r->[3]) {
                $ra->[$hdrs{Amino_Acid_Change}] = $r->[4] if ( $opt_M );
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
    -m        Output only mutations that can be counted toward mutation burden calculation
    -S        Skip the silent mutation that are "unknown" in Status.  By default: all silent "Novel" and "COSMIC" only mutations will be kept.
        Silent mutation with entries in both COSMIC and dbSNP will be filtered.

    -N  If set, keep all intronic and UTR in the output, but will be set as "unknown".
    -B        double
        If set, filter all variants with strand bias "2:0" and p-value < double, useful for plasma ctDNA sequencing.  Suggest to start with 0.01
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

    -O  Indicate to keep IO genes.  Currently defined as "HLA-*".  By default, all HLA-* genes will be filtered out.

    --recalfreq  Indicate to re-calculate the variant frequency in the cohort.  Used when VCF is processed one at a time or combining
        two cohorts.

    --GMAF double
        When the GMAF is greater than specified, it's considered common SNP and filtered out.  Default: 0.0025 (0.25%)

    --ruledir  dirpath
        default is /users/kdld047/work/NGS/Wee1/Rules

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

    --lastaa  filepath
        default is /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/last_critical_aa.txt

    --blackgenes  filepath
        default is /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/blackgenes.txt

USAGE
    exit(0);
}
