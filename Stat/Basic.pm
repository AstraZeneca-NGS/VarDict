package Stat::Basic;
use strict;

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 EXAMPLES

=head1 AUTHOR

  Zhongwu Lai

=cut

=head1 METHODS

=cut

use strict;

=head2 new

  Usage   :
  Function:
  Args    :
  Returns : An object reference

=cut

sub new {
    my $proto = shift;
    my $class = ref( $proto ) || $proto;

    my $self = {};

    my %p = @_;
    while( my ($k, $v) = each %p ) {
        $k = lc( $k );
        $k =~ s/^-/_/;
        $self->{ $k } = $v;
    }
    bless $self, $class;
    return $self;
}

=head2 _rm_null

  Usage   : $newarray = $obj->_rm_null($a);
  Function: Remove null values in the array
  Args    : $a	An array reference
  Returns : Reference to new array

=cut

sub _rm_null {
    my $this = shift;
    my $a = shift;
    my @na = ();
    foreach(@$a) { 
        push( @na, $_ ) if ( defined($_) && /\d/ );
    }
    return \@na;
}

=head2 mean

  Usage   : $mean = $obj->mean($a);
  Function: Calculate the mean for a given array, empty values are ignored.
  Args    : $a	An array reference or an array of values.  Any undef or missing values are ignored.

=cut

sub mean {
    my $this = shift;
    my $a = $_[0];
    my $n = 0;
    return undef unless( defined($a) );
    if ( ref($a) ) {  # it's a reference to array
	return 0 if ( @$a == 0 );
	foreach(@$a) { $n++ unless ( length($_) == 0 ); }
    } else {
	foreach(@_) { $n++ unless ( length($_) == 0 ); }
    }
    return $this->sum(@_)/$n;
}

=head2 geomean

  Usage   : $gmean = $obj->geomean($a);
  Function: Calculate the geomean for a given array
  Args    : $a	An array reference

=cut

sub geomean {
    my $this = shift;
    my $a = shift;
    return 0 if ( @$a == 0 );
    my $p = 1;
    my $c = 0;
    foreach(@$a) {
        next unless( $_ );
	$p *= $_;
	$c++;
    }
    return $p**(1/$c);
}

=head2 sum

  Usage   : $sum = $obj->sum($a);
  Function: Calculate the sum for a given array
  Args    : $a	An array reference or an actual array of values

=cut

sub sum {
    my $this = shift;
    my $a = $_[0];
    my $sum = 0;
    if ( ref( $a ) ) {
	foreach( @$a ) {
	    next unless ($_);
	    next unless (/\d/);
	    $sum += $_;
	}
    } else {
	foreach( @_ ) {
	    next unless ($_);
	    next unless (/\d/);
	    $sum += $_;
        }
    }
    return $sum; 
}

=head2 min

  Usage   : $min = $obj->min($a);
  Function: Find the minimum value for a given array
  Args    : $a	An array reference

=cut

sub min {
    my $this = shift;
    my $a = shift;
    my $min = $a->[0];
    foreach ( @$a ) {
        if ( $_ < $min ) {
            $min = $_;
        }
    }
    return $min;
}   
    
=head2 max

  Usage   : $max = $obj->max($a);
  Function: Find the maximum value for a given array
  Args    : $a	An array reference

=cut
    
sub max {
    my $this = shift;
    my $a = shift;
    my $max = $a->[0];
    foreach( @$a ) {
        if ( $_ > $max ) {
            $max = $_;
        }
    }
    return $max;
}

=head2 abmax

  Usage   : $abmax = $obj->abmax($a);
  Function: Find the absolute maximum value for a given array, so for [-5, 3], it will return -5, not 3.
  Args    : $a	An array reference

=cut
    
sub abmax {
    my $this = shift;
    my $a = shift;
    my $max = abs($a->[0]);
    my $abmax = $a->[0];
    foreach( @$a ) {
        if ( abs($_) > $max ) {
            $max = abs($_);
	    $abmax = $_;
        }
    }
    return $abmax;
}

=head2 abmin

  Usage   : $abmin = $obj->abmin($a);
  Function: Find the absolute minimum value for a given array, so for [-5, -2, 0.1, 3], it will return 0.1, not -5.
  Args    : $a	An array reference

=cut
    
sub abmin {
    my $this = shift;
    my $a = shift;
    my $min = abs($a->[0]);
    my $abmin = $a->[0];
    foreach( @$a ) {
        if ( abs($_) < $min ) {
            $min = abs($_);
	    $abmin = $_;
        }
    }
    return $abmin;
}

=head2 median

  Usage   : $median = $obj->median($a);
  Function: Calculate the median for a given array
  Args    : $a	An array reference

=cut

sub median {
    my $this = shift;
    my $aa = shift;
    $aa = $this->_rm_null($aa);

    return $aa->[0] if (@$aa == 1);
    my @sizes = sort { $a <=> $b } @$aa;
    if ( $#sizes%2 == 0 ) {
	return $sizes[ $#sizes/2 ];
    }

    return ($sizes[@sizes/2 - 1] + $sizes[@sizes/2])/2;
}

=head2 prctile

  Usage   : $prctile = $obj->prctile($a, $prctile, $flag);
  Function: Calculate the percentile values.  The same algorithm as in Matlab.  For n-element array:
  	    1. The sorted values in array are taken to be the 100(0.5/n), 100(1.5/n), ..., 
	       100([n-0.5]/n) percentiles.
	    2. Linear interpolation is used to compute percentiles for percent values between 
	       100(0.5/n) and 100([n-0.5]/n).
	    3. The minimum or maximum values in X are assigned to percentiles for percent values 
	       outside that range.
  Args    : $a	An array reference
  	    $prctile  A single percent value [0-100].  Negative values will be treated as 0.  Values
	    	      greater than 100 will be treated as 100.
	    $flag     Indicate whether the input array is already sorted.  Default 0, or not sorted.

=cut

sub prctile {
    my $this = shift;
    my $aa = shift;
    my $p  = shift;  $p = 0 if ($p < 0); $p = 100 if ( $p > 100 );
    my $flag = shift;  $flag = 0 unless ($flag);

    unless( $flag ) {
	my @sizes = sort { $a <=> $b } @$aa;
	$aa = \@sizes;
    }
    my $n = @$aa + 0;
    my $min = $aa->[0];
    my $max = $aa->[$n - 1];
    return $min if ( $p <= 100 * 0.5 / $n );
    return $max if ( $p >= 100 * ($n - 0.5)/$n );
    my $prev = 0;  my $next = $n - 1;
    while( $next - $prev != 1 ) {
        my $tmp = int( ($prev+$next)/2 );
	if ( $p > 100 * ($tmp + 0.5)/$n ) {
	    $prev = $tmp;
	} elsif ( $p < 100 * ($tmp + 0.5)/$n ) {
	    $next = $tmp;
	} else {
	    return $aa->[$tmp];
	}
    }

    my $tp = 100 * ($prev + 0.5)/$n;
    my $tn = 100 * ($next + 0.5)/$n;
    return $aa->[$prev] + ($aa->[$next] - $aa->[$prev]) * ($p - $tp)/($tn - $tp);
}

=head2 iqr

  Usage   : $var = $obj->iqr($a);
  Function: Calculate the IQR for a given array
  Args    : $a	An array reference

=cut

sub iqr {
    my $this = shift;
    my $a = shift;

    return $this->prctile($a, 75) - $this->prctile($a, 25);
}

=head2 var

  Usage   : $var = $obj->var($a, $mode);
  Function: Calculate the variance for a given array
  Args    : $a	An array reference
  	    $mode  0 or 1.  0 is the default.  Normalize by n-1
	    	   1 normalize by n

=cut

sub var {
    my $this = shift;
    my $a = shift;
    my $mode = shift;
    $mode = $mode ? $mode : 0;

    return 0 if ( @$a <= 1 );
    my $mean = $this->mean($a);
    my $var = 0;
    my $n = 0;
    foreach( @$a ) {
	next unless( /\d/ );
        $n++;
        $var += ($_ - $mean) * ($_ - $mean);
    }
    $var = $mode ? $var/$n : $var / ($n - 1);
    return $var;
}

=head2 std

  Usage   : $std = $obj->std($a, $mode);
  Function: Calculate the standard deviation for a given array
  Args    : $a	An array reference
  	    $mode  0 or 1.  0 is the default.  Normalize by n-1
	    	   1 normalize by n

=cut

sub std {
    my $this = shift;
    my $a = shift;
    my $mode = shift;
    $mode = $mode ? $mode : 0;

    my $var = $this->var($a, $mode);
    return sqrt($var);
}

=head2 mad

  Usage   : $mad = $obj->mad($a, $mode);
  Function: Calculate the mean/median absolute deviation for a given array
  Args    : $a	An array reference
  	    $mode  0 returns the mean absolute deviation.  Calculation is
	    	     based on mean
	    	   1 returns the median absolute deviation.  Calculation is
		     based on median

=cut

sub mad {
    my ($this, $a, $mode) = @_;
    $mode = $mode ? $mode : 0;

    $a = $this->_rm_null($a);
    my $m = $mode ? $this->median($a) : $this->mean($a);
    my @d = map { abs($_ - $m); } @$a;
    return $mode ? $this->median(\@d) : $this->mean(\@d);
}

=head2 rstd

  Usage   : $rstd = $obj->rstd($a);
  Function: Calculate the robust standard deviation for a given array
            rstd = MedianAbsoluteDeviation * 1.4826
  Args    : $a	An array reference

=cut

sub rstd {
    my $this = shift;
    my $a = shift;

    return $this->mad($a, 1) * 1.4826;
}

=head2 filter

  Usage   : $aref = $obj->filter($a, $cmp, $val);
  Function: Filter an array by a given value and criteria
  Args    : $a	An Array reference
  	    $cmp  Filtering method, options are: eq gt ge lt le
	    $val  The criteria
  Returns : An array reference

=cut

sub filter {
    my ($this, $a, $cmp, $val) = @_;
    my @tmp = ();
    die "Stat::Basic::filter expects three arguments" unless (defined($a) || $cmp || defined($val));
    foreach( @$a ) {
        if ( $cmp eq "eq" ) {
	    push(@tmp, $_) if ( $_ == $val );
	} elsif ( $cmp eq "gt" ) {
	    push(@tmp, $_) if ( $_ > $val );
	} elsif ( $cmp eq "ge" ) {
	    push(@tmp, $_) if ( $_ >= $val );
	} elsif ( $cmp eq "lt" ) {
	    push(@tmp, $_) if ( $_ < $val );
	} elsif ( $cmp eq "le" ) {
	    push(@tmp, $_) if ( $_ <= $val );
	}
    }
    return \@tmp;
}

=head2 standardize

  Usage   : $aref = $obj->standardize( $a, $mode );
  Function: Standardize an array
  Args    : $a	An Array reference
  	    $mode  0 or 1.  If 0 (the default), use standard deviation.
	    	   If 1, use median based mad.  s = x/std for 0.
		   s = x/mad for 1.  For constant arrays, where std and
		   mad will be 0, return array of 0s.
  Returns : An array reference

=cut

sub standardize {
    my $this = shift;
    my $a = shift;
    my $mode = shift;

    my @tmp = ();
    my $m = $mode ? $this->mad($a, 1) : $this->std($a);
    foreach( @$a ) {
        push(@tmp, $m == 0 ? 0 : $_/$m);
    }
    return \@tmp;
}

=head2 zscore

  Usage   : $aref = $obj->zscore( $a, $mode );
  Function: Standardize an array to z score
  Args    : $a	An Array reference
  	    $mode  0 or 1.  If 0 (the default), use standard deviation.
	           z = (x - mean)/std.
	    	   If 1, use median based mad. z = (x-median)/mad.
  Returns : An array reference

=cut

sub zscore {
    my $this = shift;
    my $a = shift;
    my $mode = shift;

    my @tmp = ();
    my $m = $mode ? $this->mad($a, 1) : $this->std($a);
    my $mid = $mode ? $this->median($a) : $this->mean($a);
    foreach( @$a ) {
	unless( /\d/ ) {
	    push(@tmp, $_); next;
	}
        push(@tmp, $m == 0 ? 0 : ($_ - $mid)/$m);
    }
    return \@tmp;
}

=head1 COPYRIGHT

  Copyright (c) 2007, AstraZeneca.  All Rights Reserved.

=cut

1;
