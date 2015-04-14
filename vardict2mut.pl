#!/usr/bin/perl -w

use Getopt::Std;
use strict;

our ($opt_n, $opt_f, $opt_F, $opt_H, $opt_D, $opt_V, $opt_M, $opt_R);

getopts( 'HMn:f:F:D:V:R:' ) || USAGE();
$opt_H && USAGE();
my $RULEDIR = "/users/kdld047/work/NGS/Wee1";
my $tp53group1 = parseMut_tp53( "$RULEDIR/Rules/DNE.txt" );
my $tp53group2 = parseMut_tp53( "$RULEDIR/Rules/TA0-25.txt" );
my $tp53group3 = parseMut_tp53( "$RULEDIR/Rules/TA25-50_SOM_10x.txt" );
my %tp53critical_aa_pos;
my @tp53hg19splice = qw(7574034 7574035 7576851 7576852 7576927 7576928 7577017 7577018 7577156 7577157 7577497 7577498 7577609 7577610 7578175 7578176 7578290 7578291 7578369 7578370 7578555 7578556 7579310 7579311 7579591 7579592 7579698 7579699 7579722 7579723 7579837 7579838 7574033 7576926 7577155 7577608 7578289 7578554 7579590 7579721 7576853 7577019 7577499 7578177 7578371 7579312 7579700 7579839);
my %tp53hg19splice = map { ($_, 1); } @tp53hg19splice;
# Set up common SNP filter
my $MINAF = $opt_f ? $opt_f : 0.05;
my $MAXRATE = $opt_R ? $opt_R : 1.00;
my $ACTMINAF = $opt_F ? $opt_F : ($MINAF/2 < 0.01 ? $MINAF/2 : 0.01);
my %filter_snp;
open( FSNP, "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/filter_common_snp.txt" );
while( <FSNP> ) {
    chomp;
    my @a = split(/\t/);
    my $key = join("-", @a[1..4]);
    $filter_snp{ $key } = 1;
}
close( FSNP );

my %snpeff_snp;
open( SESNP, "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/snpeffect_export_Polymorphic.txt" );
while( <SESNP> ) {
    chomp;
    my @a = split(/\t/);
    my $key = $a[11] ? join("-", @a[11,2]) : $a[5];
    $snpeff_snp{ $key } = 1;
}
close( SESNP );

# Set up artifact filter
my %filter_art;
open( ART, "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/filter_common_artifacts.txt" );
while( <ART> ) {
    chomp;
    my @a = split(/\t/);
    my $key = join("-", @a[1..4]);
    $filter_art{ $key } = 1;
}
close( ART );

# Set up actionable protein mutation
my %act_hot;
open( ACTHOT, "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/actionable_hotspot.txt" );
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
open( ACT, "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/actionable.txt" );
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
	    } elsif ( $a[5] eq "inframe-del" ) {
	        push(@{ $rules{ "inframe-del" }->{ $a[0] } }, [@a[1,2,3,4]]);
	    } elsif ( $a[5] eq "inframe-ins" ) {
	        push(@{ $rules{ "inframe-ins" }->{ $a[0] } }, [@a[1,2,3,4]]);
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
open( HOT, "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/Compendia.MS7.Hotspot.txt" );
while( <HOT> ) {
    chomp;
    my @a = split(/\t/);
    my $key = join("-", @a[1..4]);
    $hotspotnt{ $key } = 1;
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
if ( $opt_M ) {
    print "SAMPLE ID	ANALYSIS FILE LOCATION	VARIANT-TYPE	GENE	SOMATIC STATUS/FUNCTIONAL IMPACT	SV-PROTEIN-CHANGE	SV-CDS-CHANGE	SV-GENOME-POSITION	SV-COVERAGE	SV-PERCENT-READS	CNA-COPY-NUMBER	CNA-EXONS	CNA-RATIO	CNA-TYPE	REARR-GENE1	REARR-GENE2	REARR-DESCRIPTION	REARR-IN-FRAME?	REARR-POS1	REARR-POS2	REARR-NUMBER-OF-READS\n";
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
    my $act = isActionable( $chr, @a[2,4,5], $a[$hdrs{Gene}], $a[$hdrs{Amino_Acid_Change}], $a[$hdrs{COSMIC_AA_Change}] );
    unless( $act ) {
	next if ( $filter_snp{ $key } );
	next if ( $snpeff_snp{ "$a[$hdrs{Gene}]-$a[$hdrs{Amino_Acid_Change}]" } );
	next if ( $filter_art{ $key } && $af < 0.3 );
    }
    next if ( $opt_D && $a[$hdrs{Depth}] < $opt_D );
    next if ( $opt_V && $a[$hdrs{VD}] < $opt_V );
    my @snps = $a[3] =~ /(rs\d+)/g;
    foreach my $rs (@snps) {
        next if ( $snpeff_snp{ $rs } );
    }
    my $aachg = $a[$aachgcol];
    my $gene = $a[$genecol];
    my $platform = $sample =~ /_([^_]+)$/ ? $1 : "";
    if ( $opt_n && $sample =~ /$opt_n/ ) {
        $sample = $1;
    } elsif ( $opt_n && $sample !~ /$opt_n/ ) {
        next;
    }
    my ($type, $fclass, $gene_coding) = @a[$typecol, $funccol, $genecodecol];
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
    } elsif ( $a[10] =~ /^[A-Z]+\d+\*$/ ) {
        $status = "likely";
    } 
    if ( length($a[10]) == 0 && ($type =~ /SPLICE/i && $type !~ /region_variant/) ) {
        $status = "likely";
	$a[10] = "splice";
    } elsif ( $type =~ /splice_donor/i || $type =~ /splice_acceptor/i ) {
        $status = "likely";
	$a[10] = "splice";
    }
    if ( $a[$classcol] eq "COSMIC" ) {
	if ( $hdrs{ COSMIC_Cnt } ) {
	    if ( $a[$hdrs{ COSMIC_Cnt }] >= 5 ) {
		$status = "likely";
	    }
	} else {
	    $status = "likely";
	}
	if ( $a[$cosmaachgcol] ) {
	    $a[$cosmaachgcol] =~ s/^p\.//;
	}
    }
    if ( isHotspotNT($chr, @a[2,4,5]) ) {
        $status = "likely";
    } elsif ( isHotspotProt( $a[$hdrs{Gene}], $a[$hdrs{Amino_Acid_Change}] ) ) {
        $status = "likely";
    }
    if ($act_som{ $key } || $act_germ{ $key }) {
	$status = "known";
    }
    if ( $act ) {
	$status = "known"; 
	next if ( $af < $ACTMINAF );
	next if ( $af < 0.2 && $act eq "germline" );
    } else {
        next if ( $type =~ /^INTRON/i || $type =~ /^SYNONYMOUS_/i || $fclass eq "SILENT" || $type =~ /splice_region_variant/ );
	next if ( $af < $MINAF );
    }
    next if ( $status ne "known" && ($type =~ /^UPSTREAM/i || $type =~ /^DOWNSTREAM/i || $type =~ /^INTERGENIC/i || $type =~ /^INTRAGENIC/i || ($type =~ /UTR_/ && $type !~ /codon/i ) || $gene_coding =~ /NON_CODING/i || $fclass =~ /^NON_CODING/i) );
    next if ( $a[$classcol] eq "dbSNP" && $status ne "known" );
    next if ( $opt_R && $status ne "known" && $a[$hdrs{ Pcnt_sample }] > $MAXRATE );
    if ( $opt_M ) {
	print join("\t", $sample, $platform, "short-variant", $gene, $status, $a[10], $a[$hdrs{cDNA_Change}], "$chr:$a[2]", $a[$hdrs{Depth}], $af*100, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"), "\n";
    } else {
	print "$_\t$status\n";
    }
}

sub isActionable {
    my ($chr, $pos, $ref, $alt, $gene, $aachg, $cosmaachg) = @_;
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
	    return "somatic" if ( $r->[0] eq $chr && $r->[1] <= $pos && $r->[2] >= $pos && (length($ref)-length($alt)) >= $r->[3] );
	}
    } elsif ( $rules{ "inframe-ins" }->{ $gene } && length($ref) < length($alt) && (length($alt)-length($ref))%3 == 0 ) {
        foreach my $r ( @{ $rules{ "inframe-ins" }->{ $gene } } ) {
	    return "somatic" if ( $r->[0] eq $chr && $r->[1] <= $pos && $r->[2] >= $pos && (length($alt)-length($ref)) >= $r->[3] );
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
    return $hotspotprot{ "$gene-$pchg" } ? 1 : 0;
}

sub classify_tp53 {
    my ($aa, $pos, $ref, $alt) = @_;
    $ref =~ s/\s+//g;
    $alt =~ s/\s+//g;
    $pos =~ s/\s+//g;
    $aa =~ s/\s+//g;
    if ( $tp53hg19splice{ $pos } && length($ref) == 1 && length($alt) == 1 ) {
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
	Use with caution.  It'll filter even if it's in COSMIC, unless if actionable. Don't use it if the sample is homogeneous.
	Use only in heterogeneous samples.
	
    -F  double
        The minimum allele frequency hotspot somatic mutations, typically lower then -f.  Default: 0.01 or half -f, whichever is less

USAGE
    exit(0);
}
