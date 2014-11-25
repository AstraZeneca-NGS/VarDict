#!/usr/bin/env perl

# USAGE: pickSeq -v delimiter column seq_ids_file sequence_file_fasta > outputfile
use vars qw($opt_u $opt_v $opt_p $opt_i $opt_I $opt_c $opt_d $opt_D $opt_P $opt_l $opt_L $opt_k $opt_K $opt_C $opt_m $opt_M $opt_a $opt_s $opt_S $opt_X);
use Getopt::Std;
use warnings;
use strict;

if ( $#ARGV < 0 ) {
    print "USAGE: $0 [-vuIC] [-d delimiter] [-p delimiter] [-i columns] [-c columns] [-D regexp] [-P regexp] [-l list] [-L char] [-k strings] [-K strings] ids_file files > outputfile\n";
    print "-v   If option is set, pick lines except those that are in ids_file.\n";
    print "-C   Indicating for cell line name picking, If option is set, set -I to true, -k and -K to '-& &_&;&.&:&/&(&)&,'.\n";
    print "-u   If option is set, pick only one/first line for one id, even if there are multiple.\n";
    print "-I   If option is set, IDs are case insensitive.\n";
    print "-d   The delimiter to split ids_file to get the id.  Default is '\\t'\n";
    print "-p   The delimiter to split files to get the id.  Default is '\\t'\n";
    print "-i   Specify which column(s) (starting from 1) in ids_file contains the uniq id.  Multiple columns are separated using ':'.  Can use range or combination, such as: 1..3:7 will get (1, 2, 3, 7).  1..5s2:9 will get (1, 3, 5, 9).\n";
    print "-c   Specify which column(s) (starting from 1) in files contains the uniq id. Multiple columns are separated using ':'.  Can use range or combination, such as: 1..3:7 will get (1, 2, 3, 7).  1..5s2:9 will get (1, 3, 5, 9).\n";
    print "-D   Regular expression to capture the ID in ids_file as \$1.\n";
    print "-P   Regular expression to capture the ID in file as \$1.\n";
    print "-l   A list of : delimited IDs.  No ids_file is expected if this option is set.\n";
    print "-L   Use in combination with -l.  The character to separate IDs specified by -l, default to ':'.\n";
    print "-k   A string to be stripped in IDs in ids_file.  Use ':' to separate multiple strings.  Should be used rarely.\n";
    print "-K   A string to be stripped in IDs in file.  Use ':' to separate multiple strings.  Should be used rarely.\n";
    print "-X   A character to seperate IDs for -k and -K.  Default to '&'.  Should be used rarely.\n";
    print "-m   A string to seperate IDs in ids_file, usually ':' or ';', if the fields contains multiple IDs.  Should be used rarely.\n";
    print "-M   A string to seperate IDs in file, usually ':' or ';', if the fields contains multiple IDs.  Should be used rarely.\n";
    print "-a   A list of : delimited IDs.  These IDs will be appended to the ID list.  Usually used to get the header line.  Should be used rarely.\n";
    print "-s   Strip the leading spaces in both files.\n";
    print "-S   Strip the trailing spaces in both files.\n";

    exit(1);
}

getopts( 'IvusSCd:p:i:c:D:P:l:L:k:K:m:M:a:X:' );
my ($del1, $delimiter);
my @col = (0);
my @col1 = (0);

if ( $opt_d ) {
    $del1 = qr/$opt_d/;
} else {
    $del1 = qr/\t/;
}

if ( $opt_p ) {
    $delimiter = qr/$opt_p/;
} else {
    $delimiter = qr/\t/;
}
if ( $opt_i ) {
    @col1 = indexes($opt_i);
}
if ( $opt_c ) {
    @col = indexes($opt_c);
}

my %ids;

$opt_X = $opt_X ? $opt_X : '&';
if ( $opt_C ) {
    $opt_I = 1;
    $opt_k = "-& &_&;&.&:&/&(&)&,";
    $opt_K = "-& &_&;&.&:&/&(&)&,";
}

my @k_str = $opt_k ? split( /[$opt_X]/, $opt_k) : ();
my @K_str = $opt_K ? split( /[$opt_X]/, $opt_K) : ();

# Automatically added any IDs if -a option is set
if ( $opt_a ) {
    my @a = split(/:/, $opt_a);
    foreach my $id (@a) {
        $id = lc($id) if ( $opt_I );
	if ( $opt_k ) {
	    foreach my $t (@k_str) {
	        $id =~ s/[$t]//g;
	    }
	}
	next unless( defined($id) );
	if ( $opt_m ) {
	    my @tmp = split( /$opt_m/, $id );
	    map { $ids{ $_ } = 1; } @tmp;
	} else {
	    $ids{ $id } = 1;
	}
    }
}
##############

if ( defined($opt_l) ) {
    my $d = $opt_L ? $opt_L : ":";
    my @a = split(/$d/, $opt_l);
    foreach my $id (@a) {
        $id = lc($id) if ( $opt_I );
	if ( $opt_k ) {
	    foreach my $t (@k_str) {
	        $id =~ s/[$t]//g;
	    }
	}
	next unless( defined($id) );
	if ( $opt_m ) {
	    my @tmp = split( /$opt_m/, $id );
	    map { $ids{ $_ } = 1; } @tmp;
	} else {
	    $ids{ $id } = 1;
	}
    }
    #%ids = map { $_ = lc($_) if ( $opt_I ); ($_, 1); } @a;
} else {
    my $id_file = shift;
    open( IDS, "$id_file" );
    while( <IDS> ) {
	chomp;
	next if ( /^#/ );
	s/^\s+// if ( $opt_s );
	s/\s+$// if ( $opt_S );
	my @line = split( /$del1/, $_, -1 );
	#next unless( length( join("", @line[@col1])) > 0 );
	my $id = join("\t", @line[@col1]);
	if ( $opt_D ) {
	    $id = $1 if ( /$opt_D/ );
	}
	$id = lc($id) if ( $id && $opt_I );
	if ( $opt_k ) {
	    foreach my $t (@k_str) {
		$id =~ s/[$t]//g;
	    }
	}
	next unless( defined($id) );
	if ( $opt_m ) {
	    my @tmp = split( /$opt_m/, $id );
	    map { $ids{ $_ } = 1; } @tmp;
	} else {
	    $ids{ $id } = 1;
	}
	#$ids{ $id } = 1 if ( $id );
    }
    close( IDS );
}

my $flag;
my %picked;
while( <> ) {
    $flag = 0;
    next if ( /^#/ );
    s/^\s+// if ( $opt_s );
    s/\s+$// if ( $opt_S );
    my $id;

    if ( $opt_P ) {
        $id = $1 if ( /$opt_P/ );
	next unless( defined($id) );
    } else {
	my @a = split( /$delimiter/, $_, -1 );
	#print STDERR "ID: '", join("", @a[@col]), "'\n";
	$id = join("\t", @a[ @col ]);
	chomp($id);
	$id =~ s/\n\t/\t/g;
	next unless( length( $id ) > 0  ); # Ignore if there's no ID
    }

    $id = lc($id) if ( $id && $opt_I );
    if ( $opt_K ) {
	foreach my $t (@K_str) {
	    $id =~ s/[$t]//g;
	}
    }
    if ( $opt_M ) {
        my @tmp = split( /$opt_M/, $id );
	map { $flag = 1 if ( $ids{ $_ } ); } @tmp;
    } else {
	$flag = 1 if ( $ids{ $id } );
    }
    if ( $flag ) {
	if ( $opt_u ) {
	    next if ( $picked{ $id } );
	    print if ( ! $opt_v );
	    $picked{ $id } = 1;
	} else {
	    print if ( ! $opt_v );
	}
    } else {
	if ( $opt_u ) {
	    next if ( $picked{ $id } );
	    print if ( $opt_v );
	    $picked{ $id } = 1;
	} else {
	    print if ( $opt_v );
	}
    }
}

sub indexes {
    my $c = shift;
    my @c = ();
    my @a = split( /:/, $c );
    foreach(@a) {
        if ( /^\d+$/ ) {
            push( @c, $_ );
        } elsif ( /^(\d+)\.\.(\d+)s(\d+)$/ ) {
            for(my $i = $1; $i <= $2; $i += $3) {
                push( @c, $i );
            }
        } elsif ( /^(\d+)\.\.(\d+)$/ ) {
            push( @c, ($1 .. $2) );
        } else {
            die "Illegal index range '$_' specified.";
        }
    }
    @c = map { $_-1; } @c;
    return @c;
}
