#!/usr/bin/perl -w

# Split bed file
use strict;

my $bed = shift;
my $sub = shift;  # the sub number
my $size = shift; # the desired splits

my ($lines, $t) = split(/\s+/, `wc -l $bed`);
my $seg = int($lines/$size);
$seg++ if ( $lines % $size >= $seg );

my $head = $seg*$sub;
print `head -$head $bed | tail -$seg`;
if ($sub == 1 ) {
    my $tail = $lines % $seg;
    print `tail -$tail $bed`;
}
