#!/usr/bin/perl -w

# Split bed file
use strict;

my $bed = shift;
my $size = shift;

my ($lines, $t) = split(/\s+/, `wc -l $bed`);
my $seg = int($lines/$size);
$seg++ if ( $lines % $size > 0 && $lines % $size <= $seg );

my $n = 0;
open(BED, $bed);
my $base = `basename $bed`; chomp $base;
my %genes;
while(<BED>) {
    chomp;
    next if ( /^track/i || /^browser/i );
    my @a = split(/\t/);
    push(@{ $genes{ $a[3] } }, $_);
}
close( BED );
my @beds = ();
my @cur = ();
for(my $i = 1; $i <= $size; $i++) {
    my $out;
    open($out, ">$base.$i" );
    push(@beds, $out);
    $cur[$i-1] = 0;
}
while( my ($g, $v) = each %genes ) {
    my $N = 0;
    my $min = $cur[$N];
    for(my $i = 1; $i < @cur; $i++) {
        ($min, $N) = ($cur[$i], $i) if ( $cur[$i] < $min );
    }
    my $out = $beds[$N];
    foreach (@$v) {
        print $out "$_\n";
    }
    $cur[$N++] += @$v;
    $N = 0 if ( $N >= $size );
}

for(my $i = 0; $i <  $size; $i++) {
    close( $beds[$i] );
}
