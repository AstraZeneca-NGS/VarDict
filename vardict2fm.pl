#!/usr/bin/perl -w

use Getopt::Std;
use strict;

our ($opt_n);

getopts( 'n:' );
# Set up common SNP filter
my %filter_snp;
open( FSNP, "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/filter_common_snp.txt" );
while( <FSNP> ) {
   chomp;
   my @a = split(/\t/);
   my $key = join("-", @a[1..4]);
   $filter_snp{ $key } = 1;
}
close( FSNP );

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

# Set up actionable rescue
my %rescue_act;
open( ACT, "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/rescue_actionable.txt" );
while( <ACT> ) {
   chomp;
   my @a = split(/\t/);
   my $key = join("-", @a[1..4]);
   $rescue_act{ $key } = $a[5] ? $a[5] : 1;
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

$opt_n = $opt_n ? qr/$opt_n/ : qr/(TCGA-..-....)-01/;
while( <> ) {
    chomp;
    my @a = split(/\t/);
    my $chr = $a[1];
    $chr = "chr$chr" unless( $chr =~ /^chr/ );
    my $key = join("-", $chr, @a[2,4,5]);
    next if ( $filter_snp{ $key } );
    next if ( $filter_art{ $key } && $a[21] < 0.3 );
    my $sample = "";
    if ( $a[0] =~ /$opt_n/ ) {
        $sample = $1;
    } else {
        next;
    }
    my $status = "unknown";
    if ( $a[44] eq "ClnSNP" ) {
        $status = "known";
    } elsif ( $a[44] eq "dbSNP_del" ) {
        $status = "likely";
    } elsif ( $a[10] =~ /^-/ && $a[6] =~ /FRAME_SHIFT/i ) {
        $status = "likely";
    } elsif ( $a[10] =~ /^[A-Z]+\d+\*$/ ) {
        $status = "likely";
    } 
    if ( length($a[10]) == 0 && $a[6] =~ /SPLICE/ ) {
        $status = "likely";
	$a[10] = "splice";
    }
    if ( $a[44] eq "COSMIC" ) {
        $status = "likely";
	if ( $a[18] ) {
	    $a[18] =~ s/^p\.//;
	    $a[10] = $a[18];
	}
    }
    if ( isHotspotNT($chr, @a[2,4,5]) ) {
        $status = "known";
    } elsif ( isHotspotProt( $a[12], $a[10] ) ) {
        $status = "known";
    }
    if ($rescue_act{ $key }) {
	$status = "known";
	$a[10] = $rescue_act{ $key } if ( $rescue_act{ $key } ne "1" );
    }
    print join("\t", $sample, "", "short-variant", $a[12], $status, $a[10], "$a[4]>$a[5]", "$chr:$a[2]", $a[20], $a[21]*100, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"), "\n";
}

sub isActionable {
    my ($chr, $pos, $ref, $alt) = @_;
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

