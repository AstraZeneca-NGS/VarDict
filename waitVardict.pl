#!/usr/bin/env perl
use warnings;
use strict;

my $job = shift;
my $n = shift;

my $done = 0;
while(1) {
    for(my $i = 1; $i <= $n; $i++) {
        $done++ if ( -e "$job.done.$i" );
    }
    last if ( $done == $n );
    $done = 0;
    sleep(180);
}
