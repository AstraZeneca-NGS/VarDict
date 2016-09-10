#!/usr/bin/perl -w

use Getopt::Std;
use strict;

our ($opt_n);

getopts( 'n:' );

$opt_n = $opt_n ? qr/$opt_n/ : qr/(TCGA-..-....)-(\d\d)/;
my %mut;
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

my %samples;
while( <> ) {
    chomp;
    my @a = split(/\t/);
    my $chr = $a[1];
    $chr = "chr$chr" unless( $chr =~ /^chr/ );
    my $key = join("-", $chr, @a[2,4,5]);
    my ($sample, $stype) = ("", "");
    if ( length($a[10]) == 0 && $a[6] =~ /SPLICE/i ) {
        $a[10] = "splice";
    } elsif ( $a[6] =~ /splice_acceptor/i || $a[6] =~ /splice_donor/i ) {
        $a[10] = "splice";
    }
    if ( $a[0] =~ /$opt_n/ ) {
        $sample = $1;
	$stype = $2;
	$samples{ $sample }->{ $stype } = $a[0];
    } else {
        next;
    }
    my ($gene, $pchg, $cchg, $pos, $depth, $freq) = ($a[$genecol], $a[$aachgcol], $a[11], "$chr:$a[2]", $a[$hdrs{ Depth }], int($a[$afcol]*100));
    my $mutk = join("\t", $gene, $pchg, $cchg, $pos);
    $mut{ $sample }->{ $mutk }->{ $stype } = "$depth\t$freq";
}

while( my ($sample, $sv) = each %mut ) {
    while( my ($mut, $mv) = each %$sv ) {
        next unless( $mv->{ "01" } );
	my $sg = ( $mv->{ "10" } || $mv->{ "11" } ) ? "germline" : (($samples{ $sample }->{ 10 } || $samples{ $sample }->{ 11 } ) ? "somatic" : "NA");
	print join("\t", $sample, $mut, $mv->{ "01" }, $sg, "NA"), "\n";
    }
}
