#!/usr/bin/env perl
# Parse a list of refseq and check CDS coverage

use warnings;
use Getopt::Std;
use strict;

our ($opt_h, $opt_H, $opt_b, $opt_D, $opt_d, $opt_s, $opt_c, $opt_S, $opt_E, $opt_n, $opt_N, $opt_e, $opt_g, $opt_x, $opt_f, $opt_r, $opt_B, $opt_z, $opt_v, $opt_p, $opt_F, $opt_C, $opt_m, $opt_Q, $opt_T, $opt_q, $opt_Z, $opt_X, $opt_P, $opt_3, $opt_k, $opt_R, $opt_G, $opt_a, $opt_o, $opt_O, $opt_V, $opt_t, $opt_y, $opt_I, $opt_i, $opt_M, $opt_J);
unless( getopts( 'hHvtzypDC3iF:d:b:s:e:S:E:n:c:g:x:f:r:B:N:Q:m:T:q:Z:X:P:k:R:G:a:o:O:V:I:M:J:' )) {
    USAGE();
}
USAGE() if ( $opt_H );
USAGE() unless ( $opt_b );
my $BAM = $opt_b; # the bam file
my $sample = $1 if ( $BAM =~ /([^\/\._]+).sorted[^\/]*.bam/ );
my $samplem = "";
if ( $opt_n ) {
    $sample = $1 if ( $BAM =~ /$opt_n/ );
}
$sample = $opt_N if ( $opt_N );
unless($sample) { $sample = $1 if ( $BAM =~ /([^\/]+)[_\.][^\/]*bam/ ); }

$opt_d = "\t" unless( $opt_d );
my $c_col = $opt_c ? $opt_c - 1 : 2;
my $S_col = $opt_S ? $opt_S - 1 : 6;
my $E_col = $opt_E ? $opt_E - 1 : 7;
my $s_col = $opt_s ? $opt_s - 1 : 9;
my $e_col = $opt_e ? $opt_e - 1 : 10;
my $g_col = $opt_g ? $opt_g - 1 : 12;
$opt_m = $opt_m ? $opt_m : 8;

$s_col = $S_col if ( $opt_S && (!$opt_s) );
$e_col = $E_col if ( $opt_E && (!$opt_e) );

$opt_F = defined($opt_F) ? $opt_F : "0x500";

my $VEXT = defined($opt_X) ? $opt_X : 3; # the extension of deletion and insertion for complex variants
$opt_P = defined($opt_P) ? $opt_P : 5;
$opt_k = defined($opt_k) ? $opt_k : 1; # Whether to perform local realignment.  Set to 0 to disable.
my $fasta = $opt_G ? $opt_G : "/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa";
my %CHRS; # Key: chr Value: chr_len
my $EXT = defined($opt_x) ? $opt_x : 0;
my $FREQ = $opt_f ? $opt_f : 0.05;
my $QRATIO = $opt_o ? $opt_o : 1.5; # The Qratio
my $BIAS = 0.05; # The cutoff to decide whether a positin has read strand bias
my $MINB = $opt_B ? $opt_B : 2; # The minimum reads for bias calculation
my $MINR = $opt_r ? $opt_r : 2; # The minimum reads for variance allele
my $GOODQ = defined($opt_q) ? $opt_q : 25; # The phred score in fastq to be considered as good base call
$opt_O = defined($opt_O) ? $opt_O : 0; # The minimun mean mapping quality to be considered
$opt_V = defined($opt_V) ? $opt_V : 0.05; # The minimun alelle frequency allowed in normal for a somatic mutation
my $RLEN = 0; # The read length
my $BUFFER = 200;
my $LOWQUAL = 10;
my $INDELSIZE = $opt_I ? $opt_I : 120;
my $MINMATCH = defined($opt_M) ? $opt_M : 25;
my %SPLICE;
if ( $opt_p ) {
    $FREQ = -1;
    $MINR = 0;
}
if ( $opt_h ) {
    print join("\t", qw(Sample Gene Chr Start End Ref Alt Depth AltDepth RefFwdReads RefRevReads AltFwdReads AltRevReads Genotype AF Bias PMean PStd QMean QStd 5pFlankSeq 3pFlankSeq)), "\n";
}

my ($BAM1, $BAM2) = split(/\|/, $BAM);
my @BAMS = split(/:/, $BAM1);
open(BAMH, "samtools view -H $BAMS[0] |");
while(<BAMH>) {
    if ( /^\@SQ/ ) {
	my $chr = $1 if ( /\s+SN:(\S+)/ );
	$CHRS{ $chr } = $1 if( /\sLN:(\d+)/ );
    }
}
close(BAMH);
$fasta = $CHRS{ 1 } ? "/ngs/reference_data/genomes/Hsapiens/GRCh37/seq/GRCh37.fa" : "/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa" unless( $opt_G );

my @SEGS = ();
if ( $opt_R ) {
    my ($chr, $reg, $gene) = split(/:/, $opt_R);
    $chr =~ /^chr/ ? ( $chr =~ s/^chr// ) : ($chr = "chr$chr") unless( $CHRS{ $chr } );
    $gene = $chr unless( $gene );
    my ($start, $end) = split(/-/, $reg);
    $end = $start unless( $end );
    $start =~ s/,//g;
    $end =~ s/,//g;
    $start -= $EXT;
    $end += $EXT;
    $start++ if ( $opt_z && $start < $end );
    $start = $end if ( $start > $end );
    push(@SEGS, [[$chr, $start, $end, $gene]]);
} else {
    my @segraw = ();
    while( <> ) {
	chomp;
	next if ( /^#/ || /^browser/ || /^track/ );
	unless($opt_a) {
	    my @a = split(/$opt_d/);
	    if ((! $opt_a) && @a == 8 && $a[6] =~ /^\d+$/ && $a[7] =~ /^\d+$/ && $a[6] >= $a[1] && $a[7] <= $a[2]) {
		$opt_a = "10:0.95"; 
		$opt_z = 1 unless( defined($opt_z) );
	    }
	}
	push(@segraw, $_);
    }
    if ($opt_a) {
	my ($pchr, $pend) = (0, 0);
	my $SI = 0;
	my %tsegs = ();
	foreach(@segraw) {
	    my ($chr, $start, $end, $gene, $score, $strand, $istart, $iend) = split(/$opt_d/);
	    $chr =~ /^chr/ ? ( $chr =~ s/^chr// ) : ($chr = "chr$chr") unless( $CHRS{ $chr } );
	    $start++ && $istart++ if ( $opt_z && $start < $end );
	    push( @{ $tsegs{ $chr } }, [$chr, $start, $end, $gene, $istart, $iend] );
	}
	while(my ($chr, $sv) = each %tsegs ) {
	    my @tmp = sort { $a->[4] <=> $b->[4]; } @$sv;
	    foreach my $tv (@tmp) {
		my ($chr, $start, $end, $gene, $istart, $iend) = @$tv;
		$SI++ if ( $pend && ($chr ne $pchr || $istart > $pend) );
		push(@{$SEGS[$SI]}, [$chr, $start, $end, $gene, $istart, $iend]);
		$pchr = $chr;
		$pend = $iend;
	    }
	}
	ampVardict();
	exit(0);
    } else {
	foreach(@segraw) {
	    my @A = split(/$opt_d/);
	    if ( (! $opt_c) && @A == 4 && $A[1] =~ /^\d+$/ && $A[2] =~ /^\d+$/ && $A[1] <= $A[2] ) {
		($c_col, $S_col, $E_col, $g_col, $s_col, $e_col) = (0, 1, 2, 3, 1, 2); 
		$opt_z = 1 unless( defined($opt_z) );  # only set if it's a 4 column BED file
	    }
	    my ($chr, $cdss, $cdse, $gene) = @A[$c_col, $S_col, $E_col, $g_col];
	    $chr =~ /^chr/ ? ( $chr =~ s/^chr// ) : ($chr = "chr$chr") unless( $CHRS{ $chr } );
	    my @starts = split(/,/, $A[$s_col]);
	    my @ends = split(/,/, $A[$e_col]);
	    my @CDS = ();
	    $gene = $chr unless( $gene );
	    for(my $i = 0; $i < @starts; $i++) {
		my ($s, $e) = ($starts[$i], $ends[$i]);
		next if ( $cdss > $e ); # not a coding exon
		last if ( $cdse < $s ); # No more coding exon
		$s = $cdss if ( $s < $cdss );
		$e = $cdse if ( $e > $cdse );
		$s -= $EXT; # unless ( $s == $cdss );
		$e += $EXT; # unless ( $e == $cdse );
		$s++ if ( $opt_z && $s < $e );
		push(@CDS, [$chr, $s, $e, $gene]);
	    }
	    push(@SEGS, \@CDS);
	}
    }
}

for(my $i = 0; $i < @SEGS; $i++) {
    for(my $j = 0; $j < @{ $SEGS[$i] }; $j++) {
	my $REF = getREF(@{ $SEGS[$i]->[$j] }[0..2]);
	if ( $BAM2 ) { # pair mode
	    if ( $opt_n ) {
		$sample = $1 if ( $BAM1 =~ /$opt_n/ );
		$samplem = $1 if ( $BAM2 =~ /$opt_n/ );
	    }
	    ($sample, $samplem) = split(/\|/, $opt_N) if ( $opt_N );
	    $samplem = $sample . "_match" unless( $samplem );
	    somdict($SEGS[$i]->[$j], toVars(@{ $SEGS[$i]->[$j] }[0..2], $BAM1, $REF), toVars(@{ $SEGS[$i]->[$j] }[0..2], $BAM2, $REF));
	} else {
	    vardict($SEGS[$i]->[$j], toVars(@{ $SEGS[$i]->[$j] }[0..2], $BAM1, $REF));
	}
    }
}

sub ampVardict {
    for(my $i = 0; $i < @SEGS; $i++) {
	my @vars = ();
	my %pos = ();
	my ($chr, $start, $end, $gene, $istart, $iend);
	for(my $j = 0; $j < @{$SEGS[$i]}; $j++) {
	    ($chr, $start, $end, $gene, $istart, $iend) = @{$SEGS[$i]->[$j]};
	    my $REF = getREF($chr, $start, $end);
	    push(@vars, toVars($chr, $start, $end, $BAM1, $REF));
	    for(my $p = $istart; $p <= $iend; $p++) {
		push(@{ $pos{ $p } }, [$j, $chr, $start, $end, $istart, $iend]);
	    }
	}
	while(my ($p, $v) = each %pos) {
	    my @gvs = (); # Good variants
	    my @ref = (); # reference
	    my $nt;
	    my $maxaf = 0;
	    my $vartype = "SNV";
	    my $flag = 0;
	    my $vref;
	    my $nocov=0; #
	    my $maxcov=0;
	    my %goodamp;
	    my @vcovs = ();
	    foreach my $amps (@$v) {
		my ($amp, $chr, $S, $E) = @$amps;
		if ( $vars[$amp]->{ $p }->{ VAR }->[0] ) {
		    my $tv = $vars[$amp]->{ $p }->{ VAR }->[0];
		    push(@vcovs, $tv->{ tcov });
		    $maxcov = $tv->{ tcov } if ( $tv->{ tcov } > $maxcov );
		    $vartype = varType($tv->{ refallele }, $tv->{ varallele });
		    if ( isGoodVar($tv, $vars[$amp]->{ $p }->{ REF }, $vartype) ) {
			push(@gvs, [$tv, "$chr:$S-$E"]);
			if ( $nt && $tv->{ n } ne $nt ) {
			    $flag = 1;
			}
			if ($tv->{ freq } > $maxaf ) {
			    ($maxaf, $nt, $vref) = ($tv->{ freq }, $tv->{ n }, $tv);
			}
			$goodamp{ "$amp-$tv->{ refallele }-$tv->{ varallele }" } = 1;
		    }
		} elsif ( $vars[$amp]->{ $p }->{ REF } ) {
		    push(@vcovs, $vars[$amp]->{ $p }->{ REF }->{ tcov });
		} else {
		    push(@vcovs, 0);
		}
		push(@ref, $vars[$amp]->{ $p }->{ REF }) if ( $vars[$amp]->{ $p }->{ REF } );
	    }
	    foreach my $t (@vcovs) {
		$nocov++ if ( $t < $maxcov/50 );  # The amplicon that has depth less than 1/50 of the max depth will be considered not working and thus not used.
	    }
	    @gvs = sort { $b->[0]->{ freq } <=> $a->[0]->{ freq } } @gvs if ( @gvs > 1 );
	    @ref = sort { $b->{ tcov } <=> $a->{ tcov } } @ref if ( @ref > 1 );
	    if ( @gvs < 1 ) { # Only referenece
		if ( $opt_p ) {
		    if ( @ref ) {
			$vref = $ref[0];
		    } else {
			print join("\t", $sample, $gene, $chr, $p, $p, "", "", 0, 0, 0, 0, 0, 0, "", 0, "0;0", 0, 0, 0, 0, 0, "", 0, 0, 0, 0, 0, 0, "", "", 0, 0, "$chr:$p-$p", "", 0, 0, 0, 0), "\n";
			next;
		    }
		} else {
		    next;
		}
	    } else {
		$vref = $gvs[0]->[0];
	    }
	    if ( $flag ) { # different good variants detected in different amplicons
		#next if ($gvs[0]->[0]->{ freq }/$gvs[1]->[0]->{ freq } < 10);  # need 10 times difference to overwrite it.
		my $gdnt = $gvs[0]->[0]->{ n };
		my @gcnt = ();
		foreach my $amps (@$v) {
		    my ($amp, $chr, $S, $E) = @$amps;
		    push(@gcnt, [$vars[$amp]->{ $p }->{ VARN }->{ $gdnt }, "$chr:$S-$E"]) if ( $vars[$amp]->{ $p }->{ VARN }->{ $gdnt } && isGoodVar( $vars[$amp]->{ $p }->{ VARN }->{ $gdnt }, $vars[$amp]->{ $p }->{ REF }) );
		}
		$flag = 0 if ( @gcnt+0 == @gvs+0 );
		@gvs = sort { $b->[0]->{ freq } <=> $a->[0]->{ freq } } @gcnt;
	    }
	    my @badv = ();
	    my $gvscnt = @gvs+0;
	    unless ( $gvscnt == @$v + 0 && (! $flag) ) {
		foreach my $amps (@$v) {
		    my ($amp, $chr, $S, $E, $iS, $iE) = @$amps;
		    next if ( $goodamp{"$amp-$vref->{ refallele }-$vref->{ varallele }"} );
		    my $tref = $vars[$amp]->{ $p }->{ VAR }->[0];
		    if ( $vref->{ sp } >= $iS && $vref->{ ep } <= $iE ) {
			if ( $vars[$amp]->{ $p }->{ VAR } ) {
			    push( @badv, [$vars[$amp]->{ $p }->{ VAR }->[0], "$chr:$S-$E"] );
			} elsif ( $vars[$amp]->{ $p }->{ REF } ) {
			    push( @badv, [$vars[$amp]->{ $p }->{ REF }, "$chr:$S-$E"] );
			} else {
			    push( @badv, [undef, "$chr:$S-$E"] );
			}
		    } elsif ( ($vref->{ sp } < $iE && $iE < $vref->{ ep }) || ($vref->{ sp } < $iS && $iS < $vref->{ ep }) ) { # the variant overlap with amplicon's primer
			#print STDERR "$iS $iE $vref->{sp} $vref->{ep} $vref->{n} $gvscnt\n";
			$gvscnt-- if ( $gvscnt > 1 );
		    }
		}
	    }
	    $flag = 0 if ( $flag && $gvscnt < @gvs+0 );
	    my @hds = qw(sp ep refallele varallele tcov cov rfc rrc fwd rev genotype freq bias pmean pstd qual qstd mapq qratio hifreq extrafreq shift3 msi msint nm hicnt hicov leftseq rightseq);
	    my @hds2 = qw(tcov cov rfc rrc fwd rev genotype freq bias pmean pstd qual qstd mapq qratio hifreq extrafreq);
	    $vartype = varType($vref->{ refallele }, $vref->{ varallele });
	    adjComplex($vref) if ( $vartype eq "Complex" );
	    print join("\t", $sample, $gene, $chr, (map { $vref->{ $_ } ? $vref->{ $_ } : 0; } @hds), $gvs[0]->[1] ? $gvs[0]->[1] : "", $vartype, $gvscnt, $gvscnt+@badv+0, $nocov, $flag);
	    print "\t", $vref->{ DEBUG } if ( $opt_D );
	    for(my $gvi = 1; $gvi < @gvs; $gvi++) {
		print "\tGood$gvi ", join(" ", (map {$gvs[$gvi]->[0]->{ $_ }; } @hds2), $gvs[$gvi]->[1]) if ( $opt_D );
	    }
	    for(my $bvi = 0; $bvi < @badv; $bvi++) {
		print "\tBad$bvi ", join(" ", (map {$badv[$bvi]->[0]->{ $_ } ? $badv[$bvi]->[0]->{ $_ } : 0; } @hds2), $badv[$bvi]->[1]) if ( $opt_D );
	    }
	    print "\n";
	}
    }
}

sub varType {
    my ($ref, $var) = @_;
    if ( length($ref) == 1 && length($var) == 1 ) {
	return "SNV";
    } elsif ( substr($ref, 0, 1) ne substr($var, 0, 1) ) {
	return "Complex";
    } elsif ( length($ref) == 1 && length($var) > 1 && $var =~ /^$ref/ ) {
	return "Insertion";
    } elsif (length($ref) > 1 && length($var) == 1 && $ref =~ /^$var/ ) {
	return "Deletion";
    }
    return "Complex";
}

sub vardict {
    my ($seg, $vars) = @_;
    my ($chr, $S, $E, $G) = @{ $seg };
    my @hds = qw(sp ep refallele varallele tcov cov rfc rrc fwd rev genotype freq bias pmean pstd qual qstd mapq qratio hifreq extrafreq shift3 msi msint nm hicnt hicov leftseq rightseq);
    push(@hds, "crispr") if ( $opt_J );
    for(my $p = $S; $p <= $E; $p++) {
	my @vts = ();
	my @vrefs = ();
	unless( $vars->{ $p }->{ VAR } ) {
	    next unless( $opt_p );
	    my $vref = $vars->{ $p }->{ REF };
	    unless($vref) {
		print join("\t", $sample, $G, $chr, $p, $p, "", "", 0, 0, 0, 0, 0, 0, "", 0, "0;0", 0, 0, 0, 0, 0, "", 0, 0, 0, 0, 0, 0, "", "", 0, 0, "$chr:$S-$E", ""), "\n";
		next;
	    }
	    push(@vts, "");
	    push(@vrefs, $vref);
	} else {
	    for(my $i = 0; $i < @{ $vars->{$p}->{VAR} }; $i++) {
		my $vref = $vars->{ $p }->{ VAR }->[$i];
		last if ( $vref->{refallele} =~ /N/ );
		my $vartype = varType($vref->{refallele}, $vref->{varallele});
		unless( isGoodVar( $vars->{ $p }->{ VAR }->[$i], $vars->{ $p }->{ REF }, $vartype ) ) {
		    next unless( $opt_p );
		}
		push(@vts, $vartype);
		push(@vrefs, $vref);
	    }
	}
	for(my $vi = 0; $vi < @vts; $vi++) {
	    my ($vartype, $vref) = ($vts[$vi], $vrefs[$vi]);
	    adjComplex($vref) if ( $vartype eq "Complex" );
	    print join("\t", $sample, $G, $chr, (map { $vref->{ $_ } ? $vref->{ $_ } : 0; } @hds), "$chr:$S-$E", $vartype);
	    print "\t", $vref->{ DEBUG } if ( $opt_D );
	    print "\n";
	}
    }
}

# Adjust the complex var
sub adjComplex {
    my $vref = shift;
    my ($refnt, $varnt) = ($vref->{ refallele }, $vref->{ varallele });
    my $n = 0;
    $n++ while( length($refnt) - $n > 1 && length($varnt) - $n > 1 && substr($refnt, $n, 1) eq substr($varnt, $n, 1) );
    if ( $n ) {
	$vref->{ sp } += $n;
	$vref->{ refallele } = substr($refnt, $n);
	$vref->{ varallele } = substr($varnt, $n);
	$vref->{ leftseq } .= substr($refnt, 0, $n);
	$vref->{ leftseq } = substr($vref->{ leftseq }, $n);
    }
    ($refnt, $varnt) = ($vref->{ refallele }, $vref->{ varallele });
    $n = 1;
    $n++ while( length($refnt) - $n > 0 && length($varnt) - $n > 0 && substr($refnt, -$n, 1) eq substr($varnt, -$n, 1) );
    if ( $n > 1 ) {
	$vref->{ ep } -= $n - 1;
	$vref->{ refallele } = substr($refnt, 0, -$n+1);
	$vref->{ varallele } = substr($varnt, 0, -$n+1);
	$vref->{ rightseq } = substr($refnt, -($n-1), $n-1) . substr($vref->{ rightseq }, 0, -$n+1);
    }
}

sub somdict {
    my ($seg, $vars1, $vars2) = @_;
    my ($chr, $S, $E, $G) = @{ $seg };
    my @hdrs = qw(tcov cov rfc rrc fwd rev genotype freq bias pmean pstd qual qstd mapq qratio hifreq extrafreq nm);
    my @hd1 = ("sp", "ep", "refallele", "varallele");
    my @hd2 = ("shift3", "msi", "msint", "leftseq", "rightseq");
    my $FISHERP = 0.01;
    for(my $p = $S; $p <= $E; $p++) {
	 next unless($vars1->{$p} || $vars2->{$p}); # both samples have no coverage
	 my $vartype = "";
	 if( ! $vars1->{ $p } ) { # no coverage for sample 1
	     $vartype = varType($vars2->{$p}->{ VAR }->[0]->{ refallele }, $vars2->{$p}->{ VAR }->[0]->{ varallele }) if ($vars2->{$p}->{ VAR });
	     next unless( $vars2->{$p}->{ VAR } && isGoodVar($vars2->{$p}->{ VAR }->[0], $vars2->{$p}->{REF}, $vartype) );
	     adjComplex($vars2->{$p}->{ VAR }->[0]) if ( $vartype eq "Complex" );
	     print join("\t", $sample, $G, $chr, (map{ $vars2->{$p}->{ VAR }->[0]->{ $_ }; } @hd1), (map{ 0; } @hdrs), (map { $vars2->{$p}->{ VAR }->[0]->{ $_ }; } (@hdrs, @hd2)), "$chr:$S-$E", "Deletion", $vartype), "\n";
	 } elsif( ! $vars2->{ $p } ) { # no coverage for sample 2
	     $vartype = varType($vars1->{$p}->{ VAR }->[0]->{ refallele }, $vars1->{$p}->{ VAR }->[0]->{ varallele }) if ($vars1->{$p}->{ VAR });
	     next unless( $vars1->{$p}->{ VAR } && isGoodVar($vars1->{$p}->{ VAR }->[0], $vars1->{$p}->{REF}, $vartype) );
	     adjComplex($vars1->{$p}->{ VAR }->[0]) if ( $vartype eq "Complex" );
	     print join("\t", $sample, $G, $chr, (map { $vars1->{$p}->{ VAR }->[0]->{ $_ }; } (@hd1, @hdrs)), (map { 0; } @hdrs), (map { $vars1->{$p}->{ VAR }->[0]->{$_};} @hd2), "$chr:$S-$E", "SampleSpecific", $vartype), "\n";
	 } else { # both samples have coverage
	     my ($v1, $v2) = ($vars1->{$p}, $vars2->{$p});
	     next unless ( $v1->{ VAR } || $v2->{ VAR } );
	     if ( $v1->{ VAR } ) {
		 my $N = 0;
		 while( $v1->{ VAR }->[$N] && isGoodVar($v1->{ VAR }->[$N], $v1->{ REF }, varType($v1->{ VAR }->[$N]->{ refallele }, $v1->{ VAR }->[$N]->{ varallele }))) {
		     my $VREF = $v1->{ VAR }->[$N];
		     my $nt = $VREF->{ n };
		     #if ( length($nt) > 1 && length($VREF->{ refallele }) == length($VREF->{ varallele }) && (! isGoodVar($v2->{ VARN }->{ $nt })) && $VREF->{ genotype } !~ /-/ && $VREF->{ genotype } !~ /m/ && $VREF->{ genotype } !~ /i/ ) {
		#	 my $fnt = substr($nt, 0, -1); $fnt =~ s/&$//;
		#	 my $lnt = substr($nt, 1); $lnt =~ s/^&//; substr($lnt, 1, 0) = "&" if ( length($lnt) > 1 );
		#	 if ( $v2->{ VARN }->{ $fnt } && isGoodVar($v2->{ VARN }->{ $fnt }, $v2->{ REF })) {
		#	     $VREF->{ sp } += length($VREF->{ refallele }) - 1;
		#	     $VREF->{ refallele } = substr($VREF->{ refallele }, -1 );
		#	     $VREF->{ varallele } = substr($VREF->{ varallele }, -1 );
		#	 } elsif ( $vars2->{ $p + length($nt) -2 }->{ VARN }->{ $lnt } && isGoodVar($vars2->{ $p + length($nt) -2 }->{ VARN }->{ $lnt }, $vars2->{$p+length($nt)-2}->{REF})) {
		#	     $VREF->{ ep } += length($VREF->{ refallele }) - 1;
		#	     $VREF->{ refallele } = substr($VREF->{ refallele }, 0, -1 );
		#	     $VREF->{ varallele } = substr($VREF->{ varallele }, 0, -1 );
		#	 }
		#     }
		     $vartype = varType($VREF->{ refallele }, $VREF->{ varallele });
		     adjComplex($VREF) if ( $vartype eq "Complex" );
		     if ( $v2->{ VARN }->{ $nt } ) {
			 #my $type = isGoodVar( $v2->{ VARN }->{ $nt }, $v2->{ REF } ) ? "Germline" : ($v2->{ VARN }->{ $nt }->{ freq } < $opt_V || $v2->{ VARN }->{ $nt }->{ cov } <= 1 ? "LikelySomatic" : "AFDiff");
			 my $type = isGoodVar( $v2->{ VARN }->{ $nt }, $v2->{ REF }, $vartype ) ? ($VREF->{ freq } > (1-$opt_V) && $v2->{ VARN }->{ $nt }->{ freq } < 0.8 && $v2->{ VARN }->{ $nt }->{ freq } > 0.2 ? "LikelyLOH" : ($v2->{ VARN }->{ $nt }->{ freq } < $opt_V || $v2->{ VARN }->{ $nt }->{ cov } <= 1 ? "LikelySomatic" : "Germline")) : ($v2->{ VARN }->{ $nt }->{ freq } < $opt_V || $v2->{ VARN }->{ $nt }->{ cov } <= 1 ? "LikelySomatic" : "AFDiff");
			 $type = "StrongSomatic" if ( $vartype eq "SNV" && isNoise($v2->{ VARN }->{ $nt }) );
			 #if ($type =~ /Somatic/) {
			 #    my $newtype = combineAnalysis($VREF, $v2->{ VARN }->{ $nt }, $chr, $p, $nt);
			 #    if ( $newtype eq "FALSE" ) {$N++; next;}
			 #    $type = $newtype if ( $newtype );
			 #}
			 #$type = "StrongSomatic" if ( isNoise($v2->{ VARN }->{ $nt }) );
			 print join("\t", $sample, $G, $chr, (map { $VREF->{ $_ }; } (@hd1, @hdrs)), (map { $v2->{ VARN }->{ $nt }->{ $_ }; } (@hdrs, @hd2)), "$chr:$S-$E", $type, $vartype), "\n";
		     } else { # sample 1 only, should be strong somatic
			 my @tvf = $v2->{ REF } ? (map { $v2->{ REF }->{ $_ } ? $v2->{ REF }->{ $_ } : 0; } @hdrs) : ($v2->{ VAR }->[0]->{ tcov } ? $v2->{ VAR }->[0]->{ tcov } : 0, map { 0; } (1..17));
			 my $type = "StrongSomatic";
			 #if ($vartype ne "SNV") {
			     $v2->{ VARN }->{ $nt }->{ cov } = 0;  # Ensure it's initialized before passing to combineAnalysis
			     my $newtype = $VREF->{ cov } < $MINR + 3 ? combineAnalysis($VREF, $v2->{ VARN }->{ $nt }, $chr, $p, $nt) : "";
			     if ( $newtype eq "FALSE" ) {$N++; next;}
			     $type = $newtype if ( $newtype );
			 #}
			 if ( $type ne "StrongSomatic" ) {
			     print join("\t", $sample, $G, $chr, (map { $VREF->{ $_ }; } (@hd1, @hdrs)), (map { $v2->{ VARN }->{ $nt }->{ $_ }; } @hdrs),(map { $VREF->{$_}; } @hd2), "$chr:$S-$E", $type, $vartype), "\n";
			 } else {
			     print join("\t", $sample, $G, $chr, (map { $VREF->{ $_ }; } (@hd1, @hdrs)), @tvf, (map { $VREF->{$_};} @hd2), "$chr:$S-$E", "StrongSomatic", $vartype), "\n";
			 }
		     }
		     $N++;
		 }
		 unless($N) {
		     next unless( $v2->{ VAR } );
		     $vartype = varType($v2->{ VAR }->[0]->{ refallele }, $v2->{ VAR }->[0]->{ varallele });
		     next unless( isGoodVar( $v2->{ VAR }->[0], $v2->{ REF }, $vartype ) );
		     # potentail LOH
		     my $nt = $v2->{ VAR }->[0]->{ n };
		     if ( $v1->{ VARN }->{ $nt } ) {
			 my $type = $v1->{ VARN }->{ $nt }->{ freq } < $opt_V ? "LikelyLOH" : "Germline";
			 adjComplex($v1->{ VARN }->{ $nt }) if ( $vartype eq "Complex" );
			 print join("\t", $sample, $G, $chr, (map { $v1->{ VARN }->{ $nt }->{ $_ }; } (@hd1, @hdrs)), (map { $v2->{ VAR }->[0]->{ $_ }; } (@hdrs, @hd2)), "$chr:$S-$E", $type, varType($v1->{ VARN }->{ $nt }->{ refallele }, $v1->{ VARN }->{ $nt }->{ varallele })), "\n";
		     } else {
			 my @th1 = $v1->{ REF } ? (map { $v1->{ REF }->{ $_ } } @hdrs) : ($v1->{ VAR }->[0]->{ tcov }, (map { 0; } @hdrs[1..$#hdrs]));
			 adjComplex($v2->{ VAR }->[0]) if ( $vartype eq "Complex" );
			 print join("\t", $sample, $G, $chr, (map { $v2->{ VAR }->[0]->{ $_ }; } @hd1), @th1, (map { $v2->{ VAR }->[0]->{ $_ }; } (@hdrs, @hd2)), "$chr:$S-$E", "StrongLOH", $vartype), "\n";
		     }
		 }
	     } elsif ( $v2->{ VAR } ) { # sample 1 has only reference
		 $vartype = varType($v2->{ VAR }->[0]->{ refallele }, $v2->{ VAR }->[0]->{ varallele });
		 next unless( isGoodVar( $v2->{ VAR }->[0], $v2->{ REF }, $vartype ) );
		 # potential LOH
		 my $nt = $v2->{ VAR }->[0]->{ n };
		 my $type = "StrongLOH";
		 $v1->{ VARN }->{ $nt }->{ cov } = 0;
		 my $newtype = $v2->{ VARN }->{ $nt }->{ cov } < $MINR + 3 ? combineAnalysis($v2->{ VARN }->{ $nt }, $v1->{ VARN }->{ $nt }, $chr, $p, $nt) : "";
		 next if ( $newtype eq "FALSE" );
		 $type = $newtype if ( $newtype );
		 my @th1 = $newtype ? (map { $v1->{ VARN }->{ $nt }->{ $_ } ? $v1->{ VARN }->{ $nt }->{ $_ } : 0; } @hdrs) : (map { $v1->{ REF }->{ $_ } ? $v1->{ REF }->{ $_ } : 0 } @hdrs);
		 adjComplex($v2->{ VAR }->[0]) if ( $vartype eq "Complex" );
		 print join("\t", $sample, $G, $chr, (map { $v2->{ VAR }->[0]->{ $_ } ? $v2->{ VAR }->[0]->{ $_ } : 0; } @hd1), @th1, (map { $v2->{ VAR }->[0]->{ $_ } ? $v2->{ VAR }->[0]->{ $_ } : 0; } (@hdrs, @hd2)), "$chr:$S-$E", $type, $vartype), "\n";
	     }
	 }
    }
}

# Taken a likely somatic indels and see whether combine two bam files still support somatic status.  This is mainly
# for Indels that softclipping overhang is too short to positively being called in one bam file, but can be called
# in the other bam file, thus creating false positives
sub combineAnalysis {
    my ($var1, $var2, $chr, $p, $nt) = @_;
    my ($bam1, $bam2) = split(/\|/, $BAM);
    print STDERR "Start Combine $p $nt\n" if ( $opt_y );
    my $REF = getREF($chr, $var1->{ sp } - $RLEN, $var1->{ ep } + $RLEN);
    my $vars = toVars($chr, $var1->{ sp } - $RLEN, $var1->{ ep } + $RLEN, "$bam1:$bam2", $REF);
    if ( $vars->{ $p }->{ VARN }->{ $nt } ) {
	my $vref = $vars->{ $p }->{ VARN }->{ $nt };
	print STDERR "Combine: 1: $var1->{ cov } comb: $vref->{ cov }\n" if ( $opt_y );
	if ( $vref->{ cov } - $var1->{ cov } >= $MINR ) {
	    map { $var2->{ $_ } = $vref->{ $_ } - $var1->{ $_ } > 0 ? $vref->{ $_ } - $var1->{ $_ } : 0; } qw(tcov cov rfc rrc fwd rev);
	    map { $var2->{ $_ } = sprintf("%.3f", ($vref->{ $_ }*$vref->{ cov } - $var1->{ $_ }*$var1->{ cov })/$var2->{ cov }); } qw(pmean qual mapq hifreq extrafreq nm);
	    map { $var2->{ $_ } = 1; } qw(pstd qstd);
	    return "FALSE" if ( $var2->{ tcov } <= 0 );
	    my $freq2 = sprintf("%.3f", $var2->{ cov }/$var2->{ tcov });
	    $var2->{ freq } = $freq2;
	    $var2->{ qratio } = $var1->{ qratio }; # Can't back calculate and should be inaccurate
	    $var2->{ genotype } = $vref->{ genotype };
	    $var2->{ bias } = strandBias($var2->{rfc}, $var2->{rrc}) . ";" . strandBias($var2->{fwd}, $var2->{rev});
	    return "Germline";
	} elsif ($vref->{ cov } < $var1->{ cov } - 2) {
	    print STDERR "Combine produce less: $chr $p $nt $vref->{ cov } $var1->{ cov }\n" if ( $opt_y );
	    return "FALSE";
	} else {
	    return "";
	}
    } else {
	return "FALSE";
    }
    return "";
}

# A variance is considered noise if the quality is below $GOODQ and there're no more than 3 reads
sub isNoise {
    my $vref = shift;
    if ( (($vref->{ qual } < 4.5 || ($vref->{ qual } < 12 && $vref->{ qstd } == 0)) && $vref->{ cov } <= 3) || ($vref->{ qual } < $GOODQ && $vref->{ freq } < 2*$opt_V && $vref->{ cov } <= 1) ) {
	$vref->{ tcov } -= $vref->{ cov };
	$vref->{ cov } = 0;
	$vref->{ fwd } = 0;
	$vref->{ rev } = 0;
	$vref->{ freq } = 0;
	$vref->{ hifreq } = 0;
	return 1;
    }
    return 0;
}
# Determine whether a variant meet specified criteria
sub isGoodVar {
    my ($vref, $rref, $type) = @_;
    return 0 unless($vref && $vref->{ refallele });
    $type = varType($vref->{ refallele }, $vref->{ varallele }) unless( $type );
    return 0 if ( $vref->{ freq } < $FREQ );
    return 0 if ( $vref->{ hicnt } < $MINR );
    return 0 if ( $vref->{ pmean } < $opt_P );
    return 0 if ( $vref->{ qual } < $GOODQ );
    #print STDERR "$vref->{sp} $vref->{ cov } $vref->{varallele} $type\n";
    #return 0 if ( $vref->{ pstd } == 0 );  # Leave it to var2vcf, especially for targeted PCR, where pstd = 0 might be expected.
    if ( $rref && $rref->{ hicnt } > $MINR && $vref->{ freq } < 0.25 ) {
	#The reference allele has much better mean mapq than var allele, thus likely false positives
	return 0 if ( ($vref->{ mapq } + (length($vref->{ refallele }) + length($vref->{ varallele })-2) < 5 && $rref->{ mapq } > 20) || ( 1 + ($vref->{ mapq } + (length($vref->{ refallele }) + length($vref->{ varallele }))))/($rref->{ mapq } + 1) < 0.25 );
	#return 0 if ( $type && $type eq "SNV" && (($vref->{ mapq } < 5 && $rref->{ mapq } > 20) || $rref->{ mapq } - $vref->{ mapq } > 30 ));
    }
    return 0 if ( $type eq "Deletion" && $SPLICE{ $vref->{sp} . "-" . $vref->{ep} } );
    return 0 if ( $vref->{ qratio } < $QRATIO );
    return 1 if ( $vref->{ freq } > 0.35 );
    return 0 if ( $vref->{ mapq } < $opt_O );
    return 0 if ( $vref->{ msi } >= 13 && $vref->{ freq } <= 0.275 && $vref->{ msint } == 1);
    return 0 if ( $vref->{ msi } >= 8 && $vref->{ freq } <= 0.2 && $vref->{ msint } > 1 );
    if ( $vref->{ bias } eq "2;1" && $vref->{ freq } < 0.25 ) {
	return 0 unless($type && $type ne "SNV" && (length($vref->{refallele}) > 3 || length($vref->{varallele}) > 3));
    }
    return 1;
}

sub addCnt {
    my ($vref, $dir, $rp, $q, $Q, $nm) = @_;  # ref dir read_position quality
    $vref->{ cnt }++;
    $vref->{ $dir }++;
    $vref->{ pmean } += $rp;
    $vref->{ qmean } += $q;
    $vref->{ Qmean } += $Q;
    $vref->{ nm } += $nm;
    $q >= $GOODQ ? $vref->{ hicnt }++ : $vref->{ locnt }++;
}

# Given a region, return the REF
sub getREF {
    my ($chr, $START, $END) = @_;
    my %REF;
    my $s_start = $START - $EXT - 700 < 1 ? 1 : $START - $EXT - 700;
    my $s_end = $END + $EXT + 700 > $CHRS{ $chr } ? $CHRS{ $chr } : $END + $EXT + 700;
    my ($header, $exon) = split(/\n/, `samtools faidx $fasta $chr:$s_start-$s_end`, 2);
    $exon =~ s/\s+//g;
    for(my $i = $s_start; $i <= $s_start + length($exon); $i++) {
	$REF{ $i } = uc(substr( $exon, $i - ($s_start), 1 ));
    }
    return \%REF;
}

# Add deletion structural variant
sub addSVDel {
    my ($sdref, $s, $e, $ms, $me, $dir, $rlen, $mlen) = @_;
    $sdref->{ cnt }++;
    $sdref->{ dir } = $dir;
    $sdref->{ start } = $s unless( $sdref->{ start } && $sdref->{ start } < $s );
    $sdref->{ end } = $e unless( $sdref->{ end } && $sdref->{ end } > $e );
    push(@{ $sdref->{ mates } }, [$ms, $me, $mlen]);
    $sdref->{ mstart } = $ms unless( $sdref->{ mstart } && $sdref->{ mstart } < $ms );
    $sdref->{ mend } = $ms+$rlen unless( $sdref->{ mend } && $sdref->{ mend } > $me );
}

# Construct a variant structure given a segment and BAM files.
sub parseSAM {
    my ($chr, $START, $END, $bam, $REF) = @_;
    my %cov;
    my %var;
    my %hash;
    my %dels5;
    my %ins;
    my %mnp; # Keep track of MNPs
    my %sclip3; # soft clipped at 3'
    my %sclip5; # soft clipped at 5'
    my @bams = split(/:/, $bam);
    $chr =~ s/^chr// if ( $opt_C );
    $REF = getREF($chr, $START, $END);
    foreach my $bami (@bams) {
	#my $tsamcnt = `samtools view $bami $chr:$START-$END | head -1`;
	#next unless( $tsamcnt || $opt_p ); # to avoid too much IO for exome and targeted while the BED is whole genome
	# Get the reference sequence
	my $SAMFILTER = $opt_F ? "-F $opt_F" : "";
	open(SAM, "samtools view $SAMFILTER $bami $chr:$START-$END |");
	my %DUP = ();
	my $DUPP = 0;
	while( <SAM> ) {
	    if ( $opt_Z ) {
		next if ( rand() <= $opt_Z );
	    }
	    my @a = split(/\t/);
	    next if ( defined($opt_Q) && $a[4] < $opt_Q ); # ignore low mapping quality reads
	    #next if ( $opt_F && ($a[1] & 0x200 || $a[1] & 0x100) ); # ignore "non-primary alignment", or quality failed reads.
	    #next if ( $opt_F && ($a[1] & 0x100) ); # ignore "non-primary alignment"
	    if ( $a[1] & 0x100 ) {
		next if ($opt_F); # -F 0 to turn it off
	    }
	    next if ( $a[9] eq "*" );
	    # filter duplicated reads if option -t is set
	    if ( $opt_t ) {
		%DUP = () if ( $a[3] ne $DUPP );
		if ( $a[7] =~ /\d/ ) {
		    next if ( $DUP{ "$a[3]-$a[6]-$a[7]" } );
		    $DUP{ "$a[3]-$a[6]-$a[7]" } = 1;
		    $DUPP = $a[3];
		} elsif ( $a[1] & 0x8 ) {  # mate not mapped
		    next if ( $DUP{ "$a[3]-$a[5]" } );
		    $DUP{ "$a[3]-$a[5]" } = 1;
		    $DUPP = $a[3];
		}
	    }
	    my $nm = 0;
	    my @segid = $a[5] =~ /(\d+)[ID]/g; # For total indels
	    my $idlen = 0; $idlen += $_ foreach(@segid);
	    if ( /NM:i:(\d+)/i ) {  # number of mismatches.  Don't use NM since it includes gaps, which can be from indels
		$nm = $1 - $idlen;
		next if ( $opt_m && $nm > $opt_m ); # edit distance - indels is the # of mismatches
	    } else {
		print STDERR "No NM tag for mismatches. $_\n" if ( $opt_y && $a[5] ne "*" );
		#next;
	    }
	    my $n = 0; # keep track the read position, including softclipped
	    my $p = 0; # keep track the position in the alignment, excluding softclipped
	    my $dir = $a[1] & 0x10 ? "-" : "+";

	    # Amplicon based calling
	    if ( $opt_a ) {
		my ($dis, $ovlp) = split( /:/, $opt_a );
		($dis, $ovlp) = (10, 0.95) unless($dis && $ovlp);
		my $rlen3= 0; $rlen3 += $1 while( $a[5] =~ /(\d+)[MD]/g ); # The total aligned length, excluding soft-clipped bases and insertions
		my ($segstart, $segend) = ($a[3], $a[3]+$rlen3-1);
		if ( $a[5] =~ /^(\d+)S/ ) {
		    my $ts1 = $segstart > $START ? $segstart : $START;
		    my $te1 = $segend < $END ? $segend : $END;
		    next unless( abs($ts1-$te1)/($segend-$segstart) > $ovlp);
		} elsif ($a[5] =~ /(\d+)S$/ ) {
		    my $ts1 = $segstart > $START ? $segstart : $START;
		    my $te1 = $segend < $END ? $segend : $END;
		    next unless( abs($te1-$ts1)/($segend-$segstart) > $ovlp);
		} else {
		    if ($a[6] eq "=" && $a[8]) {
			($segstart, $segend) = $a[8] > 0 ? ($segstart, $segstart+$a[8]-1) : ($a[7], $a[7]-$a[8]-1);
		    }
		    # No segment overlapping test since samtools should take care of it
		    my $ts1 = $segstart > $START ? $segstart : $START;
		    my $te1 = $segend < $END ? $segend : $END;
		    next unless( (abs($segstart - $START) <= $dis && abs($segend - $END) <= $dis ) && abs(($ts1-$te1)/($segend-$segstart)) > $ovlp);
		}
	    }
	    if ( $a[1] & 0x8 ) { # Mate unmapped, potential insertion
		# to be implemented
	    } else {
		if ( $a[6] eq "=" ) {
		    if ( /\tSA:Z:(\S+)/ ) {
			if ( $a[1] & 0x800 ) { # the supplementary alignment
			    next; # Ignore the supplmentary for now so that it won't skew the coverage
			}
		    }
		}
	    }

	    # Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's algorithm to figure out indels
	    my $offset = 0;
	    if( $a[5] =~ s/^(\d+)D// ) { $a[3] += $1; }
	    $a[5] =~ s/\d+D$//;
	    $a[5] =~ s/^(\d+)I/$1S/;
	    $a[5] =~ s/(\d+)I$/$1S/;
	    while( $idlen > 0 && $opt_k ) {
		my $flag = 0;
		if ($a[5] =~ /^(\d+)S(\d+)([ID])/) {
		    my $tslen = $1 + ($3 eq "I" ? $2 : 0);
		    $tslen .= "S";
		    $a[3] += $3 eq "D" ? $2 : 0;
		    $a[5] =~ s/^\d+S\d+[ID]/$tslen/;
		    $flag = 1;
		}
		if ($a[5] =~ /(\d+)([ID])(\d+)S$/) {
		    my $tslen = $3 + ($2 eq "I" ? $1 : 0);
		    $tslen .= "S";
		    $a[5] =~ s/\d+[ID]\d+S$/$tslen/;
		    $flag = 1;
		}
		if ($a[5] =~ /^(\d+)S(\d+)M(\d+)([ID])/) {
		    my $tslen = $1 + $2 + ($4 eq "I" ? $3 : 0);
		    my $tmid = $2;
		    $tslen .= "S";
		    if ( $tmid <= 10 ) {
			$a[3] += $2 + ($4 eq "D" ? $3 : 0);
			$a[5] =~ s/^\d+S\d+M\d+[ID]/$tslen/;
			$flag = 1;
		    }
		}
		if ($a[5] =~ /(\d+)([ID])(\d+)M(\d+)S$/) {
		    my $tslen = $4 + $3 + ($2 eq "I" ? $1 : 0);
		    my $tmid = $3;
		    $tslen .= "S";
		    if ( $tmid <= 10 ) {
			$a[5] =~ s/\d+[ID]\d+M\d+S$/$tslen/;
			$flag = 1;
		    }
		}

		# The following two clauses to make indels at the end of reads as softly
		# clipped reads and let VarDict's algorithm identify indels
		if ($a[5] =~ /^(\d)M(\d+)([ID])(\d+)M/) {
		    my $tslen = $1 + ($3 eq "I" ? $2 : 0);
		    my $mlen = $4;
		    if ( $1 <= 8 ) {
			$a[3] += $1 + ($3 eq "D" ? $2 : 0);
			my $n = 0;
			$n++ while( $n < $mlen && $REF->{ $a[3]+$n } && $REF->{ $a[3]+$n } ne substr($a[9], $tslen+$n, 1));
			$tslen += $n;
			$mlen -= $n;
			$a[3] += $n;
			$a[5] =~ s/^\dM\d+[ID]\d+M/${tslen}S${mlen}M/;
			$flag = 1;
		    }
		}
		if ($a[5] =~ /(\d+)([ID])(\d)M$/) {
		    my $tslen = $3 + ($2 eq "I" ? $1 : 0);
		    $tslen .= "S";
		    if ( $3 <= 8 ) {
			$a[5] =~ s/\d+[ID]\d+M$/$tslen/;
			$flag = 1;
		    }
		}

		# Combine two deletions and insertion into one complex if they are close
		if ($a[5] =~ /^(.*?)(\d+)M(\d+)D(\d)M(\d+)I(\d)M(\d+)D(\d+)M/) {
		    my $mid = $4 + $6;
		    if ( $mid <= 15 ) {
			my $tslen = $4 + $5 + $6;
			my $dlen = $3 + $4 + $6 + $7;
			my $ov5 = $1;
			my $refoff = $a[3] + $2;
			my $rdoff = $2;
			my $RDOFF = $2;
			my $rm = $8;
			if ( $ov5 ) {
			    my @rdp = $ov5 =~ /(\d+)[MIS]/g;  # read position
			    my @rfp = $ov5 =~ /(\d+)[MD]/g;  # reference position
			    $refoff += $_ foreach(@rfp); 
			    $rdoff += $_ foreach(@rdp); 
			}
			my $rn = 0;
			$rn++ while( $REF->{ $refoff + $rn } && $REF->{ $refoff + $rn } eq substr($a[9], $rdoff + $rn, 1) );
			$RDOFF += $rn;
			$dlen -= $rn;
			$tslen -= $rn;
			if ( $tslen <= 0 ) {
			    $dlen -= $tslen;
			    $rm += $tslen;
			    $tslen = $dlen ."D" . $rm ."M";
			} else {
			    $tslen = "${dlen}D${tslen}I${rm}M";
			}
			#print STDERR "B: $rn $RDOFF $dlen $tslen $a[5]\n";
			$a[5] =~ s/\d+M\d+D\dM\d+I\dM\d+D\d+M/${RDOFF}M$tslen/;
			$flag = 1;
		    }
		}

		# Combine two close deletions (<10bp) into one
		if ($a[5] =~ /(\d+)D(\d)M(\d+)([DI])(\d+I)?/ ) {
		    my $dlen = $1 + $2 + ($4 eq "D" ? $3 : 0);
		    my $ilen = $2 + ($4 eq "I" ? $3 : 0);
		    my $istr = $5;
		    if ( $istr && $4 eq "D") {
			$ilen += $1 if ( $istr =~ /(\d+)I/ );
		    }
		    $a[5] =~ s/\d+D\dM\d+[DI](\d+I)?/${dlen}D${ilen}I/;
		    $flag = 1;
		}

		# Combine two close indels (<10bp) into one
		if ($a[5] =~ /(\d+)I(\d)M(\d+)([DI])(\d+I)?/ ) {
		    my $dlen = $2 + ($4 eq "D" ? $3 : 0);
		    my $ilen = $1 + $2 + ($4 eq "I" ? $3 : 0);
		    my $istr = $5;
		    if ( $istr && $4 eq "D" ) {
			$ilen += $1 if ( $istr =~ /(\d+)I/ );
		    }
		    $a[5] =~ s/\d+I\dM\d+[DI](\d+I)?/${dlen}D${ilen}I/;
		    $flag = 1;
		}

		last unless ( $flag );
	    }
	    # The following two clauses to capture sometimes mis-softly clipped reads by aligner
	    if ($a[5] =~ /^(.*?)(\d+)M(\d+)S$/) {
		my $ov5 = $1;
		my $mch = $2;
		my $soft = $3;
		my $refoff = $a[3] + $2;
		my $rdoff = $2;
		if ( $ov5 ) {
		    my @rdp = $ov5 =~ /(\d+)[MIS]/g;  # read position
		    my @rfp = $ov5 =~ /(\d+)[MD]/g;  # reference position
		    $refoff += $_ foreach(@rfp); 
		    $rdoff += $_ foreach(@rdp); 
		}
		my $rn = 0;
		my %RN = ();
		while( $rn + 1 < $soft && $REF->{ $refoff + $rn + 1} && $REF->{ $refoff + $rn + 1} eq substr($a[9], $rdoff + $rn + 1, 1) && ord(substr($a[10], $rdoff + $rn + 1, 1))-33 > $LOWQUAL ) {
		    $rn++; 
		    $RN{ $REF->{ $refoff + $rn + 1} }++;
		}
		my @rn_nt = keys %RN; # don't adjust if homopolymer
		if ( ($rn > 3 && @rn_nt > 1) || ($REF->{ $refoff } && $REF->{ $refoff } eq substr($a[9], $rdoff, 1) )) {
		    $mch += $rn + 1;
		    $soft -= $rn + 1;
		    if ( $soft > 0 ) {
			$a[5] =~ s/\d+M\d+S$/${mch}M${soft}S/;
		    } else {
			$a[5] =~ s/\d+M\d+S$/${mch}M/;
		    }
		}
		if ( $rn == 0 ) {
		    $rn++ while( $rn < $mch && $REF->{ $refoff - $rn - 1 } && $REF->{ $refoff - $rn - 1 } ne substr($a[9], $rdoff - $rn - 1, 1) );
		    if ( $rn > 0 && $rn < $mch ) {
			$soft += $rn;
			$mch -= $rn;
			$a[5] =~ s/\d+M\d+S$/${mch}M${soft}S/;
		    }
		}
	    }
	    if ($a[5] =~ /^(\d+)S(\d+)M/ ) {
		my $mch = $2;
		my $soft = $1;
		my $rn = 0;
		my %RN = ();
		while( $rn + 1 < $soft && $REF->{ $a[3] - $rn - 2} && $REF->{ $a[3] - $rn - 2} eq substr($a[9], $soft - $rn - 2, 1) && ord(substr($a[10], $soft - $rn - 2, 1))-33 > $LOWQUAL ) {
		    $rn++;
		    $RN{ $REF->{ $a[3] - $rn - 2} }++;
		}
		my @rn_nt = keys %RN; # don't adjust if homopolymer
		if ( ($rn > 3 && @rn_nt > 1) || ($REF->{ $a[3] - 1 } && $REF->{ $a[3] - 1 } eq substr($a[9], $soft - 1, 1))) {
		    $mch += $rn + 1;
		    $soft -= $rn + 1;
		    if ( $soft > 0 ) {
			$a[5] =~ s/^\d+S\d+M/${soft}S${mch}M/;
		    } else {
			$a[5] =~ s/^\d+S\d+M/${mch}M/;
		    }
		    $a[3] -= $rn + 1;
		    $rn++;
		}
		if ( $rn == 0 ) {
		    $rn++ while( $rn < $mch && $REF->{ $a[3] + $rn } && $REF->{ $a[3] + $rn } ne substr($a[9], $soft + $rn, 1));
		    if( $rn > 0 && $rn < $mch ) {
			$soft += $rn;
			$mch -= $rn;
			$a[5] =~ s/^\d+S\d+M/${soft}S${mch}M/;
			$a[3] += $rn;
		    }
		}
	    }

	    my $start = $a[3];
	    my @cigar = $a[5] =~ /(\d+)([A-Z])/g;
	    my @segs = $a[5] =~ /(\d+)[MI]/g; # Only match and insertion counts toward read length
	    my @segs2 = $a[5] =~ /(\d+)[MIS]/g; # For total length, including soft-clipped bases
	    my $rlen = 0; $rlen += $_ foreach(@segs); #$stat->sum(\@segs); # The read length for matched bases
	    next if ( $opt_M && $rlen < $MINMATCH );
	    my $rlen2= 0; $rlen2 += $_ foreach(@segs2); #$stat->sum(\@segs2); # The total length, including soft-clipped bases
	    $RLEN = $rlen2 if ($rlen2 > $RLEN); # Determine the read length

	    for(my $ci = 0; $ci < @cigar; $ci += 2) {
		my $m = $cigar[$ci];
		my $C = $cigar[$ci+1];
		$C = "S" if ( ($ci == 0 || $ci == @cigar - 2) && $C eq "I" ); # Treat insertions at the edge as softclipping
		if ( $C eq "N" ) {
		    $SPLICE{ ($start-1) . "-" . ($start+$m-1) }++;
		    $start += $m;
		    $offset = 0;
		    next;
		} elsif ( $C eq "S" ) {
		    if ( $ci == 0 ) { # 5' soft clipped
			# align softclipped but matched sequences due to mis-softclipping
			while( $m-1 >= 0 && $start - 1 > 0 && $start - 1 <= $CHRS{ $chr } && $REF->{ $start-1 } && $REF->{ $start-1 } eq substr($a[9], $m-1, 1) && ord(substr($a[10],$m-1, 1))-33 > 10) {
			    $hash{ $start - 1 }->{ $REF->{ $start - 1 } }->{ cnt } = 0 unless( $hash{ $start - 1 }->{ $REF->{ $start - 1 } }->{ cnt } );
			    addCnt($hash{ $start - 1 }->{ $REF->{ $start - 1 } }, $dir, $m, ord(substr($a[10],$m-1, 1))-33, $a[4], $nm);
			    $cov{ $start - 1 }++;
			    $start--; $m--;
			}
			if ( $m > 0 ) {
			    my $q = 0;
			    my $qn = 0;
			    my $lowqcnt = 0;
			    for(my $si = $m-1; $si >= 0; $si--) {
				last if ( substr($a[9], $si, 1) eq "N" );
				my $tq = ord(substr($a[10], $si, 1))-33;
				$lowqcnt++ if ($tq <= 12);
				last if ( $lowqcnt > 1 );
				$q += $tq;
				$qn++;
			    }
			    if ( $qn >= 1 && $qn > $lowqcnt && $start >= $START - $BUFFER && $start <= $END + $BUFFER ) {
				for(my $si = $m-1; $m - $si <= $qn; $si--) {
				    $sclip5{ $start }->{ nt }->[$m-1-$si]->{ substr($a[9], $si, 1) }++;
				    $sclip5{ $start }->{ seq }->[$m-1-$si]->{ substr($a[9], $si, 1) }->{ cnt } = 0 unless( $sclip5{ $start }->{ seq }->[$m-1-$si]->{ substr($a[9], $si, 1) }->{ cnt });
				    addCnt($sclip5{ $start }->{ seq }->[$m-1-$si]->{ substr($a[9], $si, 1) }, $dir, $si - ($m-$qn), ord(substr($a[10], $si, 1))-33, $a[4], $nm);
				}
				$sclip5{ $start }->{ cnt } = 0 unless( $sclip5{ $start }->{ cnt } );
				addCnt( $sclip5{ $start }, $dir, $m, $q/$qn, $a[4], $nm);
			    }
			}
			$m = $cigar[$ci];
		    } elsif ( $ci == @cigar - 2 ) { # 3' soft clipped
			while( $n < length($a[9]) && $REF->{ $start } && $REF->{ $start } eq substr($a[9], $n, 1) && ord(substr($a[10], $n, 1))-33 > 10) {
			    $hash{ $start }->{ $REF->{ $start } }->{ cnt } = 0 unless( $hash{ $start }->{ $REF->{ $start } }->{ cnt } );
			    addCnt($hash{ $start }->{ $REF->{ $start } }, $dir, $rlen2-$p, ord(substr($a[10], $n, 1))-33, $a[4], $nm);
			    $cov{$start}++;
			    $n++; $start++; $m--; $p++;
			}
			if ( length($a[9]) - $n > 0 ) {
			    my $q = 0;
			    my $qn = 0;
			    my $lowqcnt = 0;
			    for(my $si = 0; $si < $m; $si++) {
				last if ( substr($a[9], $n+$si, 1) eq "N" );
				my $tq = ord(substr($a[10], $n+$si, 1))-33;
				$lowqcnt++ if ($tq <= 12);
				last if ( $lowqcnt > 1 );
				$q += $tq;
				$qn++;
			    }
			    if ( $qn >= 1 && $qn > $lowqcnt && $start >= $START - $BUFFER && $start <= $END + $BUFFER ) {
				for(my $si = 0; $si < $qn; $si++) {
				    $sclip3{ $start }->{ nt }->[$si]->{ substr($a[9], $n+$si, 1) }++;
				    $sclip3{ $start }->{ seq }->[$si]->{ substr($a[9], $n+$si, 1) }->{ cnt } = 0 unless( $sclip3{ $start }->{ seq }->[$si]->{ substr($a[9], $n+$si, 1) }->{ cnt } );
				    addCnt($sclip3{ $start }->{ seq }->[$si]->{ substr($a[9], $n+$si, 1) }, $dir, $qn - $si, ord(substr($a[10], $n+$si, 1))-33, $a[4], $nm);
				}
				$sclip3{ $start }->{ cnt } = 0 unless( $sclip3{ $start }->{ cnt } );
				addCnt( $sclip3{ $start }, $dir, $m, $q/$qn, $a[4], $nm);
			    }
			}
		    }
		    $n += $m;
		    $offset = 0;
		    $start = $a[3];  # had to reset the start due to softclipping adjustment
		    next;
		} elsif ( $C eq "H" ) {
		    $offset = 0;
		    next;
		} elsif ( $C eq "I" ) {
		    $offset = 0;
		    my $s = substr($a[9], $n, $m);
		    my $q = substr($a[10], $n, $m);
		    my $ss = "";
		    my ($multoffs, $multoffp, $nmoff) = (0, 0, 0); # For multiple indels within 10bp
		    if ( $opt_k && $cigar[$ci+2] && $cigar[$ci+2] <= $VEXT && $cigar[$ci+3] =~ /M/ && $cigar[$ci+5] && $cigar[$ci+5] =~ /[ID]/ ) {
			$s .= "#" . substr($a[9], $n+$m, $cigar[$ci+2]);
			$q .= substr($a[10], $n+$m, $cigar[$ci+2]); 
			$s .= $cigar[$ci+5] eq "I" ? ("^" . substr($a[9], $n+$m+$cigar[$ci+2], $cigar[$ci+4])) : ("^" . $cigar[$ci+4]);
			$q .= $cigar[$ci+5] eq "I" ? substr($a[10], $n+$m+$cigar[$ci+2], $cigar[$ci+4]) : substr($a[10], $n+$m+$cigar[$ci+2], 1);
			$multoffs += $cigar[$ci+2] + ($cigar[$ci+5] eq "D" ? $cigar[$ci+4] : 0);
			$multoffp += $cigar[$ci+2] + ($cigar[$ci+5] eq "I" ? $cigar[$ci+4] : 0);
			if ( $cigar[$ci+6] && $cigar[$ci+7] eq "M" ) {
			    my ($toffset, $tss, $tq) = findOffset($start+$multoffs, $n+$m+$multoffp, $cigar[$ci+6], \$a[9], \$a[10], $REF, \%cov);
			    ($offset, $ss) = ($toffset, $tss);
			    $q .= $tq;
			}
			$ci += 4;
		    } else {
			if ( $opt_k && $cigar[$ci+3] && $cigar[$ci+3] =~ /M/ ) {
			    my $vsn = 0;
			    for(my $vi = 0; $vsn <= $VEXT && $vi < $cigar[$ci+2]; $vi++) {
				last if ( substr($a[9], $n+$m+$vi, 1) eq "N" );
				last if ( ord(substr($a[10], $n+$m+$vi, 1))-33 < $GOODQ );
				if ($REF->{ $start+$vi } && substr($a[9], $n+$m+$vi, 1) ne $REF->{ $start+$vi }) {
				    $offset = $vi+1;
				    $nmoff++;
				    $vsn = 0;
				} elsif ($REF->{ $start+$vi } && substr($a[9], $n+$m+$vi, 1) eq $REF->{ $start+$vi }) {
				    $vsn++;
				}
			    }
			    if ($offset) {
				$ss .= substr($a[9], $n+$m, $offset);
				$q .= substr($a[10], $n+$m, $offset);
				for( my $osi = 0; $osi < $offset; $osi++ ) {
				    $cov{ $start + $osi }++;
				}
			    }
			}
		    }
		    $s = "$s&$ss" if ( $offset > 0);
		    if ( $start - 1 >= $START && $start -1 <= $END && $s !~ /N/ ) {
			$ins{ $start - 1 }->{ "+$s" }++;
			$hash{ $start - 1 }->{ I }->{ "+$s" }->{ $dir }++;
			my $hv = $hash{ $start - 1 }->{ I }->{ "+$s" };
			$hv->{ cnt }++;
			my $tp = $p < $rlen-$p ? $p + 1: $rlen-$p;
			my $tmpq = 0;
			for(my $i = 0; $i < length($q); $i++) {
			    $tmpq += ord(substr($q, $i, 1))-33; 
			}
			$tmpq /= length($q);
			unless( $hv->{ pstd } ) {
			    $hv->{ pstd } = 0;
			    $hv->{ pstd } = 1 if ($hv->{ pp } && $tp != $hv->{ pp });
			}
			unless( $hv->{ qstd } ) {
			    $hv->{ qstd } = 0;
			    $hv->{ qstd } = 1 if ($hv->{ pq } && $tmpq != $hv->{ pq });
			}
			$hv->{ pmean } += $tp;
			$hv->{ qmean } += $tmpq;
			$hv->{ Qmean } += $a[4];
			$hv->{ pp } = $tp;
			$hv->{ pq } = $tmpq;
			$tmpq >= $GOODQ ? $hv->{ hicnt }++ : $hv->{ locnt }++;
			$hv->{ nm } += $nm - $nmoff;

			# Adjust the reference count for insertion reads
			if ( $REF->{ $start - 1 } && $hash{ $start - 1 }->{ $REF->{ $start - 1 } } && substr($a[9], $n-1, 1) eq $REF->{ $start - 1 } ) {
			    #subCnt($hash{ $start - 1 }->{ $REF->{ $start - 1 } }, $dir, $tp, $tmpq, $a[4], $nm);
			    subCnt($hash{ $start - 1 }->{ substr($a[9], $n-1, 1) }, $dir, $tp, ord(substr($a[10], $n-1, 1))-33, $a[4], $nm - $nmoff);
			}
			# Adjust count if the insertion is at the edge so that the AF won't > 1
			if ( $ci == 2 && $cigar[1] =~ /S|H/ ) {
			    my $ttref = $hash{ $start - 1 }->{ $REF->{ $start - 1 } };
			    $ttref->{ $dir }++;
			    $ttref->{ cnt }++;
			    $ttref->{ pstd } = $hv->{ pstd };
			    $ttref->{ qstd } = $hv->{ qstd };
			    $ttref->{ pmean } += $tp;
			    $ttref->{ qmean } += $tmpq;
			    $ttref->{ Qmean } += $a[4];
			    $ttref->{ pp } = $tp;
			    $ttref->{ pq } = $tmpq;
			    $ttref->{ nm } += $nm - $nmoff;
			    #$cov{ $start - 1 }->{ $REF->{ $start - 1 } }++;
			    $cov{ $start - 1 }++;
			}
		    }
		    $n += $m+$offset+$multoffp;
		    $p += $m+$offset+$multoffp;
		    $start += $offset+$multoffs;
		    next;
		} elsif ( $C eq "D" ) {
		    $offset = 0;
		    my $s = "-$m";
		    my $ss = "";
		    my $q1 = substr($a[10], $n-1, 1);
		    my $q = "";
		    my ($multoffs, $multoffp, $nmoff) = (0, 0, 0); # For multiple indels within $VEXT bp 
		    if ( $opt_k && $cigar[$ci+2] && $cigar[$ci+2] <= $VEXT && $cigar[$ci+3] =~ /M/ && $cigar[$ci+5] && $cigar[$ci+5] =~ /[ID]/ ) {
			$s .= "#" . substr($a[9], $n, $cigar[$ci+2]);
			$q .= substr($a[10], $n, $cigar[$ci+2]);
			$s .= $cigar[$ci+5] eq "I" ? ("^" . substr($a[9], $n+$cigar[$ci+2], $cigar[$ci+4])) : ("^" . $cigar[$ci+4]);
			$q .= $cigar[$ci+5] eq "I" ? substr($a[10], $n+$cigar[$ci+2], $cigar[$ci+4]) : "";
			$multoffs += $cigar[$ci+2] + ($cigar[$ci+5] eq "D" ? $cigar[$ci+4] : 0);
			$multoffp += $cigar[$ci+2] + ($cigar[$ci+5] eq "I" ? $cigar[$ci+4] : 0);
			if ( $cigar[$ci+6] && $cigar[$ci+7] eq "M" ) {
			    my $vsn = 0;
			    my $tn = $n + $multoffp;
			    my $ts = $start + $multoffs + $m;
			    for(my $vi = 0; $vsn <= $VEXT && $vi < $cigar[$ci+6]; $vi++) {
				last if ( substr($a[9], $tn+$vi, 1) eq "N" );
				last if ( ord(substr($a[10], $tn+$vi, 1))-33 < $GOODQ );
				last if ( $REF->{ $ts+$vi } && $REF->{ $ts+$vi } eq 'N' );
				if ($REF->{ $ts+$vi } && substr($a[9], $tn+$vi, 1) ne $REF->{ $ts+$vi }) {
				    $offset = $vi+1;
				    $nmoff++;
				    $vsn = 0;
				} elsif ( $REF->{ $ts+$vi } && substr($a[9], $tn+$vi, 1) eq $REF->{ $ts+$vi }) {
				    $vsn++;
				}
			    }
			    if ($offset) {
				$ss .= substr($a[9], $tn, $offset);
				$q .= substr($a[10], $tn, $offset);
			    }
			}
			$ci += 4;
		    } elsif ( $opt_k && $cigar[$ci+2] && $cigar[$ci+3] eq "I" ) {
			$s .= "^" . substr($a[9], $n, $cigar[$ci+2]);
			$q .= substr($a[10], $n, $cigar[$ci+2]);
			$multoffp += $cigar[$ci+2];
			if ( $cigar[$ci+4] && $cigar[$ci+5] eq "M" ) {
			    my $vsn = 0;
			    my $tn = $n + $multoffp;
			    my $ts = $start + $m;
			    for(my $vi = 0; $vsn <= $VEXT && $vi < $cigar[$ci+4]; $vi++) {
				last if ( substr($a[9], $tn+$vi, 1) eq "N" );
				last if ( ord(substr($a[10], $tn+$vi, 1))-33 < $GOODQ );
				last if ( $REF->{ $ts+$vi } && $REF->{ $ts+$vi } eq 'N' );
				if ($REF->{ $ts+$vi } && substr($a[9], $tn+$vi, 1) ne $REF->{ $ts+$vi }) {
				    $offset = $vi+1;
				    $nmoff++;
				    $vsn = 0;
				} elsif ( $REF->{ $ts+$vi } && substr($a[9], $tn+$vi, 1) eq $REF->{ $ts+$vi }) {
				    $vsn++;
				}
			    }
			    if ($offset) {
				$ss .= substr($a[9], $tn, $offset);
				$q .= substr($a[10], $tn, $offset);
			    }
			}
			$ci += 2;
		    } else {
			if ( $opt_k && $cigar[$ci+3] && $cigar[$ci+3] =~ /M/ ) {
			    my $vsn = 0;
			    for(my $vi = 0; $vsn <= $VEXT && $vi < $cigar[$ci+2]; $vi++) {
				last if ( substr($a[9], $n+$vi, 1) eq "N" );
				last if ( ord(substr($a[10], $n+$vi, 1))-33 < $GOODQ );
				last if ( $REF->{ $start+$m+$vi } && $REF->{ $start+$m+$vi } eq 'N' );
				if ($REF->{ $start+$m+$vi } && substr($a[9], $n+$vi, 1) ne $REF->{ $start+$m+$vi }) {
				    $offset = $vi+1;
				    $nmoff++;
				    $vsn = 0;
				} elsif ( $REF->{ $start+$m+$vi } && substr($a[9], $n+$vi, 1) eq $REF->{ $start+$m+$vi }) {
				    $vsn++;
				}
			    }
			    if ($offset) {
				$ss .= substr($a[9], $n, $offset);
				$q .= substr($a[10], $n, $offset);
			    }
			}
		    }
		    $s = "$s&$ss" if ( $offset > 0 );
		    my $q2 = substr($a[10], $n+$offset, 1);
		    $q .= ord($q1) > ord($q2) ? $q1 : $q2;
		    if ( $start >= $START && $start <= $END ) {
			$hash{ $start }->{ $s }->{ $dir }++;
			$dels5{ $start }->{ $s }++;
			my $hv = $hash{ $start }->{ $s };
			$hv->{ cnt }++;
			my $tp = $p < $rlen-$p ? $p + 1: $rlen-$p;
			my $tmpq = 0;
			for(my $i = 0; $i < length($q); $i++) {
			    $tmpq += ord(substr($q, $i, 1))-33; 
			}
			$tmpq /= length($q);
			unless( $hv->{ pstd } ) {
			    $hv->{ pstd } = 0;
			    $hv->{ pstd } = 1 if ($hv->{ pp } && $tp != $hv->{ pp });
			}
			unless( $hv->{ qstd } ) {
			    $hv->{ qstd } = 0;
			    $hv->{ qstd } = 1 if ($hv->{ pq } && $tmpq != $hv->{ pq });
			}
			$hv->{ pmean } += $tp;
			$hv->{ qmean } += $tmpq;
			$hv->{ Qmean } += $a[4];
			$hv->{ pp } = $tp;
			$hv->{ pq } = $tmpq;
			$hv->{ nm } += $nm - $nmoff;
			$tmpq >= $GOODQ ? $hv->{ hicnt }++ : $hv->{ locnt }++;
			for(my $i = 0; $i < $m; $i++) {
			    #$cov{ $start+$i }->{ $s }++;
			    $cov{ $start+$i }++;
			}
		    }
		    $start += $m+$offset+$multoffs;
		    $n += $offset+$multoffp;
		    $p += $offset+$multoffp;
		    next;
		}

		# Now dealing with matching part
		my $nmoff = 0;
		my $moffset = 0;
		for(my $i = $offset; $i < $m; $i++) {
		    my $trim = 0;
		    if ( $opt_T ) {
			if ( $dir eq "+" ) {
			    $trim = 1 if ( $n > $opt_T );
			} else {
			    $trim = 1 if ( $rlen2 - $n > $opt_T );
			}
		    }
		    my $s = substr($a[9], $n, 1);
		    if ( $s eq "N" ) {
			$start++; $n++; $p++;
			next;
		    }
		    my $q = ord(substr($a[10], $n, 1))-33;
		    my $qbases = 1;
		    my $qibases = 0;
		    # for more than one nucleotide mismatch
		    my $ss = "";
		    # More than one mismatches will only perform when all nucleotides have quality > $GOODQ
		    # Update: Forgo the quality check.  Will recover later
		    while(($start + 1) >= $START && ($start + 1) <= $END && ($i + 1) < $m && $REF->{$start} && $REF->{$start} ne 'N' && substr($a[9], $n, 1) ne $REF->{$start} && $q >= $GOODQ) {
			#last if (ord(substr($a[10], $n+1, 1))-33 < $GOODQ);
			last if (substr($a[9], $n+1, 1) eq "N" );
			last if ($REF->{$start+1} && $REF->{$start+1} eq 'N');
			if ( substr($a[9], $n+1, 1) ne $REF->{ $start + 1 } ) {
			    $ss .= substr($a[9], $n+1, 1);
			    $q += ord(substr($a[10], $n+1, 1))-33;
			    $qbases++;
			    $n++;
			    $p++;
			    $i++;
			    $start++;
			    $nmoff++;
			} else {
			    last;
			}
		    }
		    $s .= "&$ss" if ( $ss );
		    my $ddlen = 0;
		    if ( $opt_k && $m-$i <= $VEXT && $cigar[$ci+2] && $cigar[$ci+3] eq "D" && $REF->{$start} && ($ss || substr($a[9], $n, 1) ne $REF->{$start}) && ord(substr($a[10], $n, 1))-33 >= $GOODQ ) {
			while($i+1 < $m) {
			    $s .= substr($a[9], $n+1, 1);
			    $q += ord(substr($a[10], $n+1, 1))-33;
			    $qbases++;
			    $i++; $n++; $p++; $start++;
			}
			$s =~ s/&//;
			$s = "-$cigar[$ci+2]&$s";
			$ddlen = $cigar[$ci+2];
			$ci += 2;
			if ( $cigar[$ci+3] && $cigar[$ci+3] eq "I" ) {
			    $s .= "^" . substr($a[9], $n+1, $cigar[$ci+2]);
			    for(my $qi = 1; $qi <= $cigar[$ci+2]; $qi++) {
				$q += ord(substr($a[10], $n+1+$qi, 1))-33;
				$qibases++;
			    }
			    $n += $cigar[$ci+2];
			    $p += $cigar[$ci+2];
			    $ci += 2;
			}
			if ( $cigar[$ci+2] && $cigar[$ci+3] eq "M" ) {
			    my ($toffset, $tss, $tq, $tnmoff) = findOffset($start+$ddlen+1, $n+1, $cigar[$ci+2], \$a[9], \$a[10], $REF, \%cov);
			    if ( $toffset ) {
				$moffset = $toffset;
				$nmoff += $tnmoff;
				$s .= "&$tss";
				for(my $qi = 0; $qi < length($tq); $qi++) {
				    $q += ord(substr($tq, $qi, 1)) - 33;
				    $qibases++;
				}
			    }
			}
		    }
		    unless( $trim ) {
			if ( $start - $qbases + 1 >= $START && $start - $qbases + 1 <= $END ) {
			    $hash{ $start - $qbases + 1 }->{ $s }->{ $dir }++;
			    $mnp{ $start - $qbases + 1 }->{ $s }++ if ( $s =~ /^[ATGC]&[ATGC]+$/ );
			    my $hv = $hash{ $start - $qbases + 1 }->{ $s };
			    $hv->{ cnt }++;
			    my $tp = $p < $rlen-$p ? $p + 1: $rlen-$p;
			    $q = $q/($qbases+$qibases);
			    unless( $hv->{ pstd } ) {
				$hv->{ pstd } = 0;
				$hv->{ pstd } = 1 if ($hv->{ pp } && $tp != $hv->{ pp });
			    }
			    unless( $hv->{ qstd } ) {
				$hv->{ qstd } = 0;
				$hv->{ qstd } = 1 if ($hv->{ pq } && $q != $hv->{ pq });
			    }
			    $hv->{ pmean } += $tp;
			    $hv->{ qmean } += $q;
			    $hv->{ Qmean } += $a[4];
			    $hv->{ pp } = $tp;
			    $hv->{ pq } = $q;
			    $hv->{ nm } += $nm - $nmoff;
			    $q >= $GOODQ ? $hv->{ hicnt }++ : $hv->{ locnt }++;
			    for(my $qi = 1; $qi <= $qbases; $qi++) {
				$cov{ $start - $qi + 1 }++;
			    }
			    if ( $s =~ /^-/ ) {
				$dels5{ $start - $qbases + 1 }->{ $s }++;
				for(my $qi = 1; $qi < $ddlen; $qi++) {
				    $cov{ $start + $qi }++;
				}
			    }
			}
		    }
		    if ( $s =~ /^-/ ) {
			$start += $ddlen;
		    }
		    $start++ unless( $C eq "I" );
		    $n++ unless( $C eq "D" );
		    $p++ unless( $C eq "D" );
		}
		if ( $moffset ) {
		    $offset = $moffset;
		    $n += $moffset;
		    $start += $moffset;
		    $p += $moffset;
		}
		last if ( $start > $END );
	    }
	}
	close( SAM );
    }

    if ( $opt_i ) {
	while( my ($intron, $icnt) = each %SPLICE ) {
	    print join("\t", $sample, $chr, $intron, $icnt), "\n";
	}
	return;
    }

    if ( $opt_k ) {
	print STDERR "Start Realigndel\n" if ( $opt_y );
	realigndel(\%hash, \%dels5, \%cov, \%sclip5, \%sclip3, $REF, $chr, \@bams);
	print STDERR "Start Realignins\n" if ( $opt_y );
	realignins(\%hash, \%ins, \%cov, \%sclip5, \%sclip3, $REF, $chr);
	print STDERR "Start Realignlgdel\n" if ( $opt_y );
	realignlgdel(\%hash, \%cov, \%sclip5, \%sclip3, $REF, $chr);
	print STDERR "Start Realignlgins\n" if ( $opt_y );
	realignlgins(\%hash, \%cov, \%sclip5, \%sclip3, $REF, $chr, \@bams);
	print STDERR "Start Realignlgins30\n" if ( $opt_y );
	realignlgins30(\%hash, \%cov, \%sclip5, \%sclip3, $REF, $chr, \@bams);
    }
    adjMNP(\%hash, \%mnp, \%cov, $chr, $REF, \%sclip5, \%sclip3);
    return (\%hash, \%cov, $REF);
}

sub toVars {
    my ($chr, $START, $END, $bam, $REF) = @_;
    my ($hash, $cov) = parseSAM($chr, $START, $END, $bam, $REF);
    my %vars; # the variant structure
    while( my ($p, $v) = each %$hash ) {
	next unless( %$v );
	next unless( $p >= $START && $p <= $END );
	my @tmp = ();
	my $vcov = 0; #the variance coverage
	my @var = ();
	my @vn = keys %$v;
	if ( @vn == 1 && $vn[0] eq $REF->{ $p } ) {
	    next unless( $opt_p || $BAM2 || $opt_a ); # ignore if only reference were seen and no pileup to avoid computation
	}
	#my @v = values %{ $cov->{ $p } };
	#my $tcov = 0; $tcov += $_ foreach(@v); #$stat->sum(\@v);
	next unless ( $cov->{ $p } );
	my $tcov = $cov->{ $p };
	my $hicov = 0;
	#$hicov += $_->{ hicnt } ? $_->{ hicnt } : 0 foreach( values %$v );
	next if ( $tcov == 0 ); # ignore when there's no coverage
	while( my ($n, $cnt) = each %$v ) {
	    if ( $n eq "I" ) {
		while( my ($in, $icnt) = each %$cnt ) {
		    $hicov += $icnt->{ hicnt } ? $icnt->{ hicnt } : 0;
		}
	    } else {
		$hicov += $cnt->{ hicnt } ? $cnt->{ hicnt } : 0;
	    }
	}
	while( my ($n, $cnt) = each %$v ) {
	    unless( $n eq "I") {
		next unless( $cnt->{ cnt } );
		my $fwd = $cnt->{ "+" } ? $cnt->{ "+" } : 0;
		my $rev = $cnt->{ "-" } ? $cnt->{ "-" } : 0;
		my $bias = strandBias($fwd, $rev);
		my $vqual = sprintf("%.1f", $cnt->{ qmean }/$cnt->{ cnt }); # base quality
		my $MQ = sprintf("%.1f", $cnt->{ Qmean }/$cnt->{ cnt }); # mapping quality
		my ($hicnt, $locnt) = ($cnt->{ hicnt } ? $cnt->{ hicnt } : 0, $cnt->{ locnt } ? $cnt->{ locnt } : 0);
		my $ttcov = ( $cnt->{ cnt } > $tcov && $cnt->{ cnt } - $tcov < $cnt->{ extracnt } ) ? $cnt->{ cnt } : $tcov;
		my $tvref = {n => $n, cov => $cnt->{ cnt }, fwd => $fwd, rev => $rev, bias => $bias, freq => $cnt->{ cnt }/$ttcov, pmean => sprintf("%.1f", $cnt->{ pmean }/$cnt->{ cnt } ), pstd => $cnt->{ pstd }, qual => $vqual, qstd => $cnt->{ qstd }, mapq => $MQ, qratio => sprintf("%.3f", $hicnt/($locnt ? $locnt : $locnt+0.5)), hifreq => ($hicov > 0 ? $hicnt/$hicov : 0), extrafreq => $cnt->{ extracnt } ? $cnt->{ extracnt }/$ttcov : 0, shift3 => 0, msi => 0, nm => sprintf("%.1f", $cnt->{ nm }/$cnt->{ cnt } ), hicnt => $hicnt, hicov => $hicov };
		push(@var, $tvref);
		if ( $opt_D ) {
		    push( @tmp, "$n:" . ($fwd + $rev) . ":F-$fwd:R-$rev:" . sprintf("%.3f", $tvref->{freq}) . ":$tvref->{bias}:$tvref->{pmean}:$tvref->{pstd}:$vqual:$tvref->{qstd}:" . sprintf("%.3f", $tvref->{hifreq}) . ":$tvref->{mapq}:$tvref->{qratio}");
		}
	    }
	}
	if ( $v->{ I } ) {
	    while( my ($n, $cnt) = each %{ $v->{ I } } ) {
		my $fwd = $cnt->{ "+" } ? $cnt->{ "+" } : 0;
		my $rev = $cnt->{ "-" } ? $cnt->{ "-" } : 0;
		my $bias = strandBias($fwd, $rev);
		my $vqual = sprintf("%.1f", $cnt->{ qmean }/$cnt->{ cnt }); # base quality
		my $MQ = sprintf("%.1f", $cnt->{ Qmean }/$cnt->{ cnt }); # mapping quality
		my ($hicnt, $locnt) = ($cnt->{ hicnt } ? $cnt->{ hicnt } : 0, $cnt->{ locnt } ? $cnt->{ locnt } : 0);
		#$hicov += $hicnt ? $hicnt : 0;
		my $ttcov = ( $cnt->{ cnt } > $tcov && $cnt->{ cnt } - $tcov < $cnt->{ extracnt } ) ? $cnt->{ cnt } : $tcov;
		my $tvref = {n => $n, cov => $cnt->{ cnt }, fwd => $fwd, rev => $rev, bias => $bias, freq => $cnt->{ cnt }/$ttcov, pmean => sprintf("%.1f", $cnt->{ pmean }/$cnt->{ cnt } ), pstd => $cnt->{ pstd }, qual => $vqual, qstd => $cnt->{ qstd }, mapq => $MQ, qratio => sprintf("%.3f", $hicnt/($locnt ? $locnt : $locnt+0.5)), hifreq => ($hicov > 0 ? $hicnt/$hicov : 0), extrafreq => $cnt->{ extracnt } ? $cnt->{ extracnt }/$ttcov : 0, shift3 => 0, msi => 0, nm => sprintf("%.1f", $cnt->{ nm }/$cnt->{ cnt } ), hicnt => $hicnt, hicov => $hicov };
		push(@var, $tvref);
		if ( $opt_D ) {
		    push( @tmp, "I$n:" . ($fwd + $rev) . ":F-$fwd:R-$rev:" . sprintf("%.3f", $tvref->{freq}) . ":$tvref->{bias}:$tvref->{pmean}:$tvref->{pstd}:$vqual:$tvref->{qstd}:" . sprintf("%.3f", $tvref->{hifreq}) . ":$tvref->{mapq}:$tvref->{qratio}" );
		}
	    }
	}
	@var = sort { $b->{ qual } * $b->{cov} <=> $a->{ qual } * $a->{cov}; } @var;
	my $maxfreq = 0;
	foreach my $tvar (@var) {
	    if ( $tvar->{ n } eq $REF->{ $p } ) {
		$vars{ $p }->{ REF } = $tvar;
	    } else {
		push( @{ $vars{ $p }->{ VAR } }, $tvar );
		$vars{ $p }->{ VARN }->{ $tvar->{ n } } = $tvar;
		$maxfreq = $tvar->{freq} if ($tvar->{freq} > $maxfreq );
	    }
	}
	unless( $opt_p || $maxfreq > $FREQ || $opt_a ) {
	    unless( $BAM2 ) {
		delete $vars{ $p };
		next;
	    }
	}
	# Make sure the first bias is always for the reference nucleotide
	my ($pmean, $pstd, $qual, $qstd, $mapq, $qratio);
	my ($rfc, $rrc) = (0, 0); # coverage for referece forward and reverse strands
	my $genotype1 = $vars{ $p }->{ REF } && $vars{ $p }->{ REF }->{ freq } >= $FREQ ? $vars{ $p }->{ REF }->{ n } : ($vars{ $p }->{ VAR } ? $vars{ $p }->{ VAR }->[0]->{ n } : $vars{ $p }->{ REF }->{ n });
	$genotype1 = "+" . (length($genotype1) - 1 ) if ( $genotype1 =~ /^\+/ );
	my $genotype2 = "";
	my $vn;

	if ( $vars{ $p }->{ REF } ) {
	    ($rfc, $rrc) = ($vars{ $p }->{ REF }->{ fwd }, $vars{ $p }->{ REF }->{ rev });
	}
	# only reference reads are observed.
	if ( $vars{ $p }->{ VAR } ) {
	    for(my $vi = 0; $vi < @{ $vars{ $p }->{ VAR } }; $vi++) {
		my $vref = $vars{ $p }->{ VAR }->[$vi];
		$genotype2 = $vref->{ n };
		$genotype2 = "+" . (length($genotype2) - 1 ) if ( $genotype2 =~ /^\+/ );
		$vn = $vref->{ n };
		my $dellen = $vn =~ /^-(\d+)/ ? $1 : 0;
		my $ep = $vn =~ /^\+/ ?  $p : ($vn =~ /^-/ ? $p + $dellen - 1: $p);
		my ($refallele, $varallele) = ("", "");
		my ($shift3, $msi, $msint) = (0, 0, "");  # how many bp can a deletion be shifted to 3 prime
		my $sp = $p;
		if ( $vn =~ /^\+/ ) {
		    unless( $vn =~ /&/ || $vn =~ /#/ ) {
			my $tseq1 = $vn;
			$tseq1 =~ s/^\+//;
			my $leftseq = join( "", (map { $REF->{ $_ }; } (($p-50 > 1 ? $p-50 : 1) .. $p)) ); # left 10 nt
			my $tseq2 = join("", (map { $REF->{ $_ }; } (($p+1) .. ($p+70 > $CHRS{ $chr } ? $CHRS{ $chr } : ($p+70)))));
			($msi, $shift3, $msint) = findMSI($tseq1, $tseq2, $leftseq);
			my ($tmsi, $tshift3, $tmsint) = findMSI($leftseq, $tseq2);
			($msi, $msint) = ($tmsi, $tmsint) if ( $msi < $tmsi );
			$msi = $shift3/length($tseq1) unless( $msi > $shift3/length($tseq1) );
			#print STDERR "$maxmsi $tseq1 $tseq2\n";
			#my $tseq2 = join("", (map { $REF->{ $_ }; } (($p+1) .. ($p+length($tseq1)))));
			#while( $tseq1 eq $tseq2 ) {
			#    $shift3 += length($tseq1);
			#    $p += length($tseq1);
			#    $tseq2 = join("", (map { $REF->{ $_ }; } (($p+1) .. ($p+length($tseq1)))));
			#}
			#$msi = $shift3/length($tseq1);
		    }
		    if ( $opt_3 ) {
			$sp += $shift3;
			$ep += $shift3;
		    }
		    ($refallele, $varallele) = ($REF->{$p}, $vn);
		    $varallele =~ s/^\+//;
		    $varallele = $REF->{$p} . $varallele;
		} elsif ( $vn =~ /^-/ ) {
		    $varallele = $vn;
		    $varallele =~ s/^-\d+//;
		    my $leftseq = join( "", (map { $REF->{ $_ }; } (($p-70 > 1 ? $p - 70 : 1) .. ($p - 1))) ); # left 10 nt
		    my $tseq = join("", (map { $REF->{ $_ }; } (($p) .. ($p+$dellen+70 > $CHRS{ $chr } ? $CHRS{ $chr } : ($p+$dellen+70)) )));
		    ($msi, $shift3, $msint) = findMSI(substr($tseq, 0, $dellen), substr($tseq, $dellen), $leftseq);
		    my ($tmsi, $tshift3, $tmsint) = findMSI($leftseq, substr($tseq, $dellen));
		    ($msi, $msint) = ($tmsi, $tmsint) if ( $msi < $tmsi );
		    $msi = $shift3/$dellen unless( $msi > $shift3/$dellen ); # if ( $shift3%$dellen == 0 );
		    unless( $vn =~ /&/ || $vn =~ /#/ || $vn =~ /\^/ ) {
			#while($REF->{ $sp } eq $REF->{ $ep+1 } ) {
			#    $sp++;
			#    $ep++;
			#    $p++;
			#    $shift3++;
			#}
			if ( $opt_3 ) {
			    $sp += $shift3;
			}
			$varallele = $REF->{ $p - 1 };
			$refallele = $REF->{ $p - 1 };
			$sp--;
		    }
		    for(my $i = 0; $i < $dellen; $i++) {
			$refallele .= $REF->{ $p + $i };
		    }
		    #print STDERR "$msi, $shift3, $msint $varallele $refallele $vn $p\n";
		} else {
		    my $tseq1 = join("", (map { $REF->{ $_ }; } (($p-30 > 1 ? ($p-30) : 1) .. ($p+1))));
		    my $tseq2 = join("", (map { $REF->{ $_ }; } (($p+2) .. ($p+70 > $CHRS{ $chr } ? $CHRS{ $chr } : ($p+70)))));
		    ($msi, $shift3, $msint) = findMSI($tseq1, $tseq2);
		    ($refallele, $varallele) = ($REF->{ $p }, $vn);
		}
		if ( $vn =~ /&([ATGC]+)/ ) {
		    my $extra = $1;
		    $varallele =~ s/&//;
		    for(my $m = 0; $m < length($extra); $m++) {
			$refallele .= $REF->{ $ep + $m + 1 };
			$genotype1 .= $REF->{ $ep + $m + 1 }; # unless( $genotype1 =~ /^-/ || $genotype1 =~ /^\+/ || length($genotype1) > 1 );
		    }
		    $ep += length($extra);
		    if ( $varallele =~ /&([ATGC]+)/ ) {
			my $vextra = $1;
			$varallele =~ s/&//;
			for(my $m = 0; $m < length($vextra); $m++) {
			    $refallele .= $REF->{ $ep + $m + 1 };
			    $genotype1 .= $REF->{ $ep + $m + 1 }; # unless( $genotype1 =~ /^-/ || $genotype1 =~ /^\+/ || length($genotype1) > 1 );
			}
			$ep += length($vextra);
		    }
		    if ( $vn =~ /^\+/ ) {
			substr($refallele, 0, 1) = "";
			substr($varallele, 0, 1) = "";
			$sp++;
		    }
		}
		if ($vn =~ /#(.+)\^(.+)/) {
		    my $mseq = $1;
		    my $tail = $2;
		    $ep += length($mseq);
		    $refallele .= join("", (map { $REF->{ $_ }; } (($ep-length($mseq)+1) .. $ep)));
		    if ( $tail =~ /^(\d+)/ ) {
			for(my $ti = 0; $ti < $1; $ti++) {
			    $refallele .= $REF->{ $ep + $ti + 1 };
			}
			$ep += $1;
		    }
		    $varallele =~ s/#//;
		    $varallele =~ s/\^(\d+)?//;
		    $genotype1 =~ s/#/m/;
		    $genotype1 =~ s/\^/i/;
		    $genotype2 =~ s/#/m/;
		    $genotype2 =~ s/\^/i/;
		} 
		if ($vn =~ /\^([ATGNC]+)/) { # for deletion followed directly by insertion in novolign
		    $varallele =~ s/\^//;
		    $genotype1 =~ s/\^/i/;
		    $genotype2 =~ s/\^/i/;
		}
		my $genotype = "$genotype1/$genotype2";
		$genotype =~ s/&//g;
		$genotype =~ s/#//g;
		$genotype =~ s/\^/i/g;

		# Perform adjustment to get as close to CRISPR site as possible
                if ( $opt_J && $shift3 && (length($refallele) != length($varallele) && substr($refallele, 0, 1) eq substr($varallele, 0, 1) ) ) {
                    unless( $sp == $opt_J || $ep == $opt_J ) {
			my $n = 0;
			my $dis = abs($opt_J - $sp) < abs($opt_J - $ep) ? abs($opt_J - $sp) : abs($opt_J - $ep);
			if ( $sp < $opt_J ) {
			    $n++ while( $sp + $n < $opt_J && $n < $shift3 && $ep + $n != $opt_J);
			    $n = 0 if ( abs($sp + $n - $opt_J) > $dis && abs($ep + $n - $opt_J) > $dis ); # Don't move if it makes it worse
			}
			if ( $ep < $opt_J && $n == 0 ) {
			    if ( abs($ep - $opt_J) <= abs($sp - $opt_J) ) {
                                $n++ while( $ep + $n < $opt_J && $n < $shift3 );
                            }
                        }
			if ( $n > 0 ) {
			    $sp += $n;
			    $ep += $n;
			    $refallele = join("", (map { $REF->{ $_ } } ($sp .. $ep)));
			    my $tva = "";
			    if ( length($refallele) < length($varallele) ) { # Insertion
				$tva = substr($varallele, 1);
				if ( length( $tva ) > 1 ) {
				    my $ttn = $n % length($tva);
				    $tva = substr($tva, $ttn) . substr($tva, 0, $ttn) unless( $ttn == 0 );
				}
			    }
			    $varallele = $REF->{ $sp } . $tva;
			}
			$vref->{ crispr } = $n;
                    }
                }
		$vref->{ leftseq } = join( "", (map { $REF->{ $_ }; } (($sp-20 < 1 ? 1 : $sp-20) .. ($sp - 1))) ); # left 20 nt
		$vref->{ rightseq } = join( "", (map { $REF->{ $_ }; } (($ep+1) .. ($ep+20 > $CHRS{ $chr } ? $CHRS{ $chr } : $ep+20))) ); # right 20 nt
		$vref->{ extrafreq } = sprintf("%.4f", $vref->{ extrafreq } ) if ($vref->{ extrafreq });
		$vref->{ freq } = sprintf("%.4f", $vref->{ freq }) if ($vref->{ freq });
		$vref->{ hifreq } = sprintf("%.4f", $vref->{ hifreq }) if ($vref->{ hifreq });
		$vref->{ msi } = $msi ? sprintf("%.3f", $msi) : $msi;
		$vref->{ msint } = $msint ? length($msint) : 0;
		$vref->{ shift3 } = $shift3;
		$vref->{ sp } = $sp;
		$vref->{ ep } = $genotype =~ /\+(\d+)/ ? $sp + $1 : $ep;
		$vref->{ refallele } = $refallele;
		$vref->{ varallele } = $varallele;
		$vref->{ genotype } = $genotype;
		$vref->{ tcov } = $tcov;
		$vref->{ rfc } = $rfc;
		$vref->{ rrc } = $rrc;
		$vref->{ bias } = $vars{$p}->{ REF } ? $vars{$p}->{ REF }->{ bias } . ";" . $vref->{ bias } : "0;" . $vref->{ bias };
		$vref->{ DEBUG } = join(" & ", @tmp) if ( $opt_D );
	    }
	}
	if ( $vars{$p}->{ REF } ) {
	    my $vref = $vars{$p}->{ REF };  # no variant reads are detected.
	    $vref->{ tcov } = $tcov;
	    $vref->{ cov } = 0;
	    $vref->{ freq } = 0;
	    $vref->{ rfc } = $rfc;
	    $vref->{ rrc } = $rrc;
	    $vref->{ fwd } = 0;
	    $vref->{ rev } = 0;
	    $vref->{ msi } = 0;
	    $vref->{ msint } = 0;
	    $vref->{ bias } .= ";0";
	    $vref->{ shift3 } = 0;
	    $vref->{ sp } = $p;
	    $vref->{ ep } = $p;
	    $vref->{ hifreq } = sprintf("%.4f", $vref->{ hifreq }) if ($vref->{ hifreq });
	    $vref->{ refallele } = $REF->{$p};
	    $vref->{ varallele } = $REF->{$p};
	    $vref->{ genotype } = "$REF->{$p}/$REF->{$p}";
	    $vref->{ leftseq } = "";
	    $vref->{ rightseq } = "";
	    $vref->{ DEBUG } = join(" & ", @tmp) if ( $opt_D );
	}
    }
    return \%vars;
}

# Find closest mismatches to combine with indels
sub findOffset {
    my ($refp, $readp, $mlen, $rdseq, $qstr, $REF, $cov) = @_;
    my ($offset, $ss, $q, $tnm) = (0, "", "", 0);
    my $vsn = 0;
    for(my $vi = 0; $vsn <= $VEXT && $vi < $mlen; $vi++) {
	last if ( substr($$rdseq, $readp+$vi, 1) eq "N" );
	last if ( ord(substr($$qstr, $readp+$vi, 1))-33 < $GOODQ );
	if ($REF->{ $refp+$vi } && substr($$rdseq, $readp+$vi, 1) ne $REF->{ $refp+$vi }) {
	    $offset = $vi+1;
	    $tnm++;
	    $vsn = 0;
	} elsif ($REF->{ $refp+$vi } && substr($$rdseq, $readp+$vi, 1) eq $REF->{ $refp+$vi }) {
	    $vsn++;
	}
    }
    if ($offset) {
	$ss .= substr($$rdseq, $readp, $offset);
	$q .= substr($$qstr, $readp, $offset);
	for( my $osi = 0; $osi < $offset; $osi++ ) {
	    $cov->{ $refp + $osi }++;
	}
    }
    return ($offset, $ss, $q, $tnm);
}

# Adjust MNP when there're breakpoints within MNP
sub adjMNP {
    my ($hash, $mnp, $cov, $chr, $REF, $sclip5, $sclip3) = @_;
    while( my ($p, $v) = each %$mnp ) {
	while( my ($vn, $mv) = each %$v ) {
	    next unless( $hash->{ $p }->{ $vn } ); # The variant is likely already been used by indel realignment
	    my $mnt = $vn; $mnt =~ s/&//;
	    for(my $i = 0; $i < length($mnt)-1; $i++) {
		my $left = substr($mnt, 0, $i+1);
		my $right = substr($mnt, -(length($mnt)-$i-1)); 
		if ( $hash->{ $p }->{ $left } ) {
		    my $tref = $hash->{ $p }->{ $left };
		    if ($tref->{ cnt } < $hash->{ $p }->{ $vn }->{ cnt } && $tref->{ pmean }/$tref->{ cnt } <= $i + 1) {
			print STDERR "  AdjMnt Left: $p $vn $tref->{ cnt }\n" if ( $opt_y );
			adjCnt($hash->{ $p }->{ $vn }, $tref);
			delete $hash->{ $p }->{ $left };
		    }
		}
		if ( $hash->{ $p+$i+1 }->{ $right } ) {
		    my $tref = $hash->{ $p+$i+1 }->{ $right };
		    if ($tref->{ cnt } < $hash->{ $p }->{ $vn }->{ cnt } && $tref->{ pmean }/$tref->{ cnt } <= length($mnt)-$i-1) {
			print STDERR "  AdjMnt Right: $p $vn $tref->{ cnt }\n" if ( $opt_y );
			adjCnt($hash->{ $p }->{ $vn }, $tref);
			$cov->{ $p } += $tref->{ cnt };
			delete $hash->{ $p+$i+1 }->{ $right };
		    }
		}
	    }
            if ( $sclip3->{ $p } ) {
                my $sc3v = $sclip3->{ $p };
                unless ($sc3v->{ used }) {
                    my $seq = findconseq( $sc3v );
                    if( $seq && length($seq) >= length($mnt) ) {
                        if ( $seq =~ /^$mnt/ ) {
                            if ( length($seq) == length($mnt) || ismatchref(substr($seq, length($mnt)), $REF, $p + length($mnt), 1 ) ) {
                                adjCnt($hash->{ $p }->{ $vn }, $sc3v);
                                $cov->{ $p } += $sc3v->{ cnt };
                                $sc3v->{ used } = 1;
                            }
                        }
                    }
                }
            }
            if ( $sclip5->{ $p + length($mnt) } ) {
                my $sc5v = $sclip5->{ $p + length($mnt) };
                unless ($sc5v->{ used }) {
                    my $seq = findconseq( $sc5v );
                    if( $seq && length($seq) >= length($mnt) ) {
                        $seq = reverse($seq);
                        if ( $seq =~ /$mnt$/ ) {
                            if ( length($seq) == length($mnt) || ismatchref(substr($seq, 0, length($seq) - length($mnt)), $REF, $p - 1, -1 ) ) {
                                adjCnt($hash->{ $p }->{ $vn }, $sc5v);
                                $cov->{ $p } += $sc5v->{ cnt };
                                $sc5v->{ used } = 1;
                            }
                        }
                    }
                }
            }
	}
    }
}

sub findMSI {
    my ($tseq1, $tseq2, $left) = @_;
    my $nmsi = 1;
    my $shift3 = 0;
    my ($maxmsi, $msicnt) = ("", 0);
    #print STDERR "$tseq1\t$tseq2\n";
    while( $nmsi <= length($tseq1) && $nmsi <= 6) {
	my $msint = substr($tseq1, -$nmsi, $nmsi);
	my $msimatch = $1 if ( $tseq1 =~ /(($msint)+)$/ );
	$msimatch = $1 if ( $left && "$left$tseq1" =~ /(($msint)+)$/ );
	my $curmsi += length($msimatch)/$nmsi;
	my $curshift3 = 0;
	if ( $tseq2 =~ /^(($msint)+)/ ) {
	    $curmsi += length($1)/$nmsi;
	}
	($maxmsi, $msicnt) = ($msint, $curmsi) if ( $curmsi > $msicnt );
	$nmsi++;
    }
    my $tseq = "$tseq1$tseq2";
    $shift3++ while(substr($tseq, $shift3, 1) eq substr($tseq2, $shift3, 1) && $shift3 < length($tseq2) );
    #print STDERR "$tseq1 $tseq2 $maxmsi $msicnt $shift3\n";
    return ($msicnt, $shift3, $maxmsi);
}

sub islowqual {
    my $q = shift;
    for(my $i = 0; $i < length($q); $i++) {
       return 1 if(ord(substr($q,$i,1))-33 < 7 && length($q) - $i > 3);
    }
    return 0;
}

# if any of the base is represented by 75%, it's a low complexity seq
sub islowcomplexseq {
    my ($seq) = @_;
    my $len = length($seq);
    return 1 if ( $len == 0 );
    my $a = $seq =~ s/A/A/g;
    return 1 if ( $a/$len > 0.75 );
    my $t = $seq =~ s/T/T/g;
    return 1 if ( $t/$len > 0.75 );
    my $g = $seq =~ s/G/G/g;
    return 1 if ( $g/$len > 0.75 );
    my $c = $seq =~ s/C/C/g;
    return 1 if ( $c/$len > 0.75 );
    return 0;
}

# this will try to realign large insertions (typically larger than 30bp)
sub realignlgins30 {
    my ($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams) = @_;
    my @tmp5;
    while(my($p, $sc5v) = each %$sclip5) {
	push(@tmp5, [$p, $sc5v, $sc5v->{cnt}]);
    }
    @tmp5 = sort {$b->[2] <=> $a->[2];} @tmp5;
    my @tmp3;
    while(my($p, $sc3v) = each %$sclip3) {
	push(@tmp3, [$p, $sc3v, $sc3v->{cnt}]);
    }
    @tmp3 = sort {$b->[2] <=> $a->[2];} @tmp3;
    for(my $i = 0; $i < @tmp5; $i++) {
	my ($p5, $sc5v, $cnt5) = @{ $tmp5[$i] };
	#print STDERR "Ins30: p5 $p5\n" if ( $opt_y );
	next if ($sc5v->{ used });
	for( my $j = 0; $j < @tmp3; $j++ ) {
	    my ($p3, $sc3v, $cnt3) = @{ $tmp3[$j] };
	    #print STDERR "Ins30: p3 $p3\n" if ( $opt_y );
	    next if ($sc3v->{ used });
	    #next unless($cnt3/$cnt5 < 7 && $cnt3/$cnt5 > 0.15 ); # they have to be somewhat balanced to prevent exception
	    next if ( $p5 - $p3 > $RLEN/1.5 );
	    next if ( $p3 - $p5 > $RLEN-10 );  # if they're too far away, don't even try
	    my $seq5 = findconseq($sc5v);
	    my $seq3 = findconseq($sc3v);
	    next unless(length($seq5) > 10 && length($seq3) > 10);
	    print STDERR "  Working lgins30: $p3 $p5 3: $seq3 $cnt3 5: ", scalar reverse($seq5), " $cnt5\n" if ( $opt_y );
	    my ($bp5, $bp3, $score) = find35match($seq5, $seq3, $p5, $p3, $REF);
	    next unless($score);
	    my $ins = $bp3 > 1 ? substr($seq3, 0, -$bp3 + 1) : $seq3;
	    $ins .= reverse(substr($seq5, 0, $bp5)) if ( $bp5 > 0 );
	    if ( islowcomplexseq($ins) ) {
		print STDERR "  Discard low complex insertion found $ins.\n" if ( $opt_y );
		next;
	    }
	    my $bi;
	    my $vref;
	    print STDERR "  Found candidate lgins30: $p3 $p5 $ins\n" if ( $opt_y );
	    if ( $p5 > $p3 ) {
		next if ( length($seq3) > length($ins) && (! ismatch(substr($seq3, length($ins)), join("", (map { $REF->{ $_ }; } ($p5 .. ($p5 + length($seq3) - length($ins)+2))) ), 1)));
		next if ( length($seq5) > length($ins) && (! ismatch(substr($seq5, length($ins)), join("", (map { $REF->{ $_ }; } ( ($p3 - (length($seq5)-length($ins)) - 2) .. ($p3 - 1)))), -1)));
		print STDERR "  Found lgins30 complex: $p3 $p5 ", length($ins), " $ins\n" if ( $opt_y );
		my $tmp = join("", (map { $REF->{ $_ }; } ($p3 .. ($p5-1))));
		if ( length($tmp) > length($ins) ) { # deletion is longer
		    $ins = -($p5 - $p3) . "^$ins";
		    $bi = $p3;
		    $hash->{ $p3 }->{ $ins }->{ cnt } = 0 unless( $hash->{ $p3 }->{ $ins }->{ cnt } );
		    $vref = $hash->{ $p3 }->{ $ins };
		} elsif (length($tmp) < length($ins) ) {
		    my $p3s = $p3 + length($tmp);
		    my $p3e = $p3s + length($seq3) - length($ins) + 2;
		    $ins = substr($ins, 0, length($ins)-length($tmp)) . "&" . substr($ins, -($p5-$p3));
		    $ins = "+$ins";
		    $bi = $p3 - 1;
		    $hash->{ $bi }->{ I }->{ $ins }->{ cnt } = 0 unless( $hash->{ $bi }->{ I }->{ $ins }->{ cnt } );
		    $vref = $hash->{ $bi }->{ I }->{ $ins };
		} else { # long MNP
		    $ins = "-" . length($ins) . "^$ins"; #substr($ins, 0, 1) . "&" . substr($ins, -(length($ins)-1));
		    $bi = $p3;
		    $hash->{ $p3 }->{ $ins }->{ cnt } = 0 unless( $hash->{ $p3 }->{ $ins }->{ cnt } );
		    $vref = $hash->{ $p3 }->{ $ins };
		}
	    } else {
		next if ( length($seq3) > length($ins) && (! ismatch(substr($seq3, length($ins)), join("", (map { $REF->{ $_ }; } ($p5 .. ($p5 + length($seq3) - length($ins)+2))) ), 1)));
		next if ( length($seq5) > length($ins) && (! ismatch(substr($seq5, length($ins)), join("", (map { $REF->{ $_ }; } ( ($p3 - (length($seq5)-length($ins)) - 2) .. ($p3 - 1)))), -1)));
		my $tmp = "";
		if ( length($ins) <= $p3 - $p5 ) { # Tandem duplication
		    for(my $p = $p5; $p < $p5 + ($p3 - $p5 - length($ins))/2; $p++) {
			$tmp .= $REF->{ $p };
		    }
		} else {
		    for(my $p = $p5; $p < $p3; $p++) {
			$tmp .= $REF->{ $p };
		    }
		}
		print STDERR "  Found lgins30: $p3 $p5 ", length($ins), " $tmp + $ins\n" if ( $opt_y );
		$ins = "+$tmp$ins";
		$bi = $p5 - 1;
		$hash->{ $bi }->{ I }->{ $ins }->{ cnt } = 0 unless( $hash->{ $bi }->{ I }->{ $ins }->{ cnt } );
		$vref = $hash->{ $bi }->{ I }->{ $ins };
	    }
	    $sc3v->{ used } = 1;
	    $sc5v->{ used } = 1;
	    $vref->{ pstd } = 1;
	    $vref->{ qstd } = 1;
	    $cov->{ $bi } += $sc5v->{ cnt };
	    print STDERR "  lgins30 Found: '$ins' $bi $bp3 $bp5\n" if ( $opt_y );
	    #adjCnt($hash->{ $bi }->{ I }->{ "+$ins" }, $hash->{$bi}->{ $REF->{$bi} }, $hash->{$bi}->{ $REF->{$bi} }) if (length($ins) <= $p3 - $p5 && $bams && $p3 - $p5 >= 5 && $p3 - $p5 < 75 && $hash->{$bi}->{ $REF->{$bi} } && noPassingReads($chr, $p5, $p3, $bams) && $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } > 2 * $hash->{$bi}->{ $REF->{$bi} }->{ cnt });
	    if ( $ins =~ /^\+/ ) {
		adjCnt($vref, $sc3v, $hash->{ $bi }->{ $REF->{ $bi } });
		adjCnt($vref, $sc5v);
		adjCnt($vref, $hash->{$bi}->{ $REF->{$bi} }, $hash->{$bi}->{ $REF->{$bi} }) if ( $bams && $p3 - $p5 >= 5 && $p3 - $p5 < $RLEN - 10 && $hash->{$bi}->{ $REF->{$bi} } && $hash->{$bi}->{ $REF->{$bi} }->{ cnt } && noPassingReads($chr, $p5, $p3, $bams) && $vref->{ cnt } > 2 * $hash->{$bi}->{ $REF->{$bi} }->{ cnt });
		my %tins = ();
		$tins{ $bi }->{ $ins } = $vref->{ cnt };
		realignins($hash, \%tins, $cov, $sclip5, $sclip3, $REF, $chr);
	    } elsif ( $ins =~ /^-/ ) {
		adjCnt($vref, $sc3v, $hash->{ $bi }->{ $REF->{ $bi } });
		adjCnt($vref, $sc5v);
		my %tdel = ();
		$tdel{ $bi }->{ $ins } = $vref->{ cnt };
		realigndel($hash, \%tdel, $cov, $sclip5, $sclip3, $REF, $chr);
	    } else {
		adjCnt($vref, $sc3v);
		adjCnt($vref, $sc5v);
	    }
	    last;
	}
    }
    print STDERR "Done: lgins30\n\n" if ( $opt_y );
}

# test whether the two soft-clipped reads match
# Returns the breakpoints in 5 and 3 prime soft-clipped reads
sub find35match {
    my ($seq5, $seq3, $p5, $p3, $REF) = @_;
    my $LONGMM = 3;
    my $max = 0;
    my ($B3, $B5) = (0, 0);
    for(my $i = 0; $i < length($seq5) - 8; $i++) {
	for(my $j = 1; $j < length($seq3) - 8; $j++) {
	    my $nm = 0;
	    my $n = 0;
	    while($n+$j <= length($seq3) && $i+$n <= length($seq5)) {
		$nm++ if (substr($seq3, -($j+$n), 1) ne substr($seq5, $i+$n, 1));
		last if ( $nm > $LONGMM );
		$n++;
	    }

	    #print STDERR "find35match: $i $j $n $nm $seq5 $seq3\n" if ( $opt_y );
	    if ( ($n+$j >= length($seq3) || $i+$n >= length($seq5)) && $n-$nm > $max && $n-$nm > 8 && $nm/$n < 0.1) {
		$max = $n - $nm;
		$B3 = $j;
		$B5 = $i;
		#print STDERR join("\t", $n+$j, length($seq3), $i+$n, length($seq5), $n, $nm, $i, $j), "\n";
		return ($B5, $B3, $max);
	    }
	}
    }
    return ($B5, $B3, $max);
}

# Realign large insertions that are not present in alignment
sub realignlgins {
    my ($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams) = @_;
    my @tmp;
    while(my($p, $sc5v) = each %$sclip5) {
	push(@tmp, [$p, $sc5v, $sc5v->{cnt}]);
    }
    @tmp = sort {$b->[2] <=> $a->[2];} @tmp;
    foreach my $t (@tmp) {
	my ($p, $sc5v, $cnt) = @$t;
	next if ($sc5v->{ used }); # already been used in 
	my $seq = findconseq($sc5v);
	next unless( $seq );
	print STDERR "  Working lgins: 5: $p $seq\n" if ( $opt_y );
	next if ($seq =~ /^.AAAAAAAA/ || $seq =~ /^.TTTTTTTT/ );
	next if ( length($seq) < 12 );
	next if ( islowcomplexseq($seq) );
	my ($bi, $ins) = findbi($seq, $p, $REF, -1, $chr);
	next unless( $bi );
	print STDERR "  Found candidate lgins from 5: $bi +$ins $p $seq Cnt: $cnt\n" if ($opt_y);
	$hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } = 0 unless( $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } );
	$hash->{ $bi }->{ I }->{ "+$ins" }->{ pstd } = 1;
	$hash->{ $bi }->{ I }->{ "+$ins" }->{ qstd } = 1;
	adjCnt($hash->{ $bi }->{ I }->{ "+$ins" }, $sc5v);
	my $rpflag = 1;  # A flag to indicate whether an insertion is a repeat
	for(my $i = 0; $i < length($ins); $i++) {
	    if ( $REF->{ $bi + 1 + $i } ne substr($ins, $i, 1) ) {
		$rpflag = 0;
		last;
	    }
	}
	$cov->{ $bi } += $sc5v->{ cnt };
	my $len = $ins =~ /&/ ? length($ins) - 1 : length($ins);
	for(my $ii = $len+1; $ii < @{ $sc5v->{ seq } }; $ii++) {
	    my $pii = $bi - $ii + $len;
	    while(my( $tnt, $tv) = each %{ $sc5v->{ seq }->[$ii] }) {
		$hash->{ $pii }->{ $tnt }->{ cnt } = 0 unless( $hash->{ $pii }->{ $tnt }->{ cnt } );
		adjCnt( $hash->{ $pii }->{ $tnt }, $tv );
		$hash->{ $pii }->{ $tnt }->{ pstd } = 1;
		$hash->{ $pii }->{ $tnt }->{ qstd } = 1;
		$cov->{ $pii } += $tv->{ cnt };
	    }
	}
	$sc5v->{ used } = 1;
	my %tins = ();
	$tins{ $bi }->{ "+$ins" } = $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt };
	realignins($hash, \%tins, $cov, $sclip5, $sclip3, $REF, $chr);
	adjCnt($hash->{ $bi }->{ I }->{ "+$ins" }, $hash->{$bi}->{ $REF->{$bi} }, $hash->{$bi}->{ $REF->{$bi} }) if ($rpflag && $bams && length($ins) >= 5 && length($ins) < $RLEN-10 && $hash->{$bi}->{ $REF->{$bi} } && $hash->{$bi}->{ $REF->{$bi} }->{ cnt } && noPassingReads($chr, $bi, $bi+length($ins), $bams) && $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } > 2 * $hash->{$bi}->{ $REF->{$bi} }->{ cnt });
    }
    @tmp = ();
    while(my($p, $sc3v) = each %$sclip3) {
	push(@tmp, [$p, $sc3v, $sc3v->{cnt}]);
    }
    @tmp = sort {$b->[2] <=> $a->[2];} @tmp;
    foreach my $t (@tmp) {
	my ($p, $sc3v, $cnt) = @$t;
	#print STDERR "Do 3: $p $cnt $sc3v->{ used }\n";
	next if ($sc3v->{ used }); # already been used in 
	my $seq = findconseq($sc3v);
	next unless( $seq );
	print STDERR "  Working lgins 3: $seq\n" if ( $opt_y );
	next if ($seq =~ /^.AAAAAAA/ || $seq =~ /^.TTTTTTT/ );
	next if ( length($seq) < 12 );
	next if ( islowcomplexseq($seq) );
	my ($bi, $ins, $be) = findbi($seq, $p, $REF, 1, $chr);
	#print STDERR "Here $seq $p $sc3v '$bi' $cnt\n";
	next unless( $bi );
	print STDERR "  Found candidate lgins from 3: $bi $ins $p $seq Cnt: $cnt\n" if ( $opt_y );
	$hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } = 0 unless( $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } );
	$hash->{ $bi }->{ I }->{ "+$ins" }->{ pstd } = 1;
	$hash->{ $bi }->{ I }->{ "+$ins" }->{ qstd } = 1;
	    #print STDERR "3B: $bi $cov->{ $bi } $ins $sc3v->{cnt}\n";
	adjCnt($hash->{ $bi }->{ I }->{ "+$ins" }, $sc3v, $hash->{ $bi }->{ $REF->{ $bi } });
	my $rpflag = 1;
	for(my $i = 0; $i < length($ins); $i++) {
	    if ( $REF->{ $bi + 1 + $i } ne substr($ins, $i, 1) ) {
		$rpflag = 0;
		last;
	    }
	}
	my $offset = $bi == $be ? ($p - $bi - 1) : -($p + $be - $bi);
	my $len = $ins =~ /&/ ? length($ins) - 1 : length($ins);
	#for(my $ii = $len-$offset; $ii < @{ $sc3v->{ seq } }; $ii++)
	for(my $ii = $len; $ii < @{ $sc3v->{ seq } }; $ii++) {
	    my $pii = $p + $ii - $len;
	    while(my( $tnt, $tv) = each %{ $sc3v->{ seq }->[$ii] }) {
		$hash->{ $pii }->{ $tnt }->{ cnt } = 0 unless( $hash->{ $pii }->{ $tnt }->{ cnt } );
		adjCnt( $hash->{ $pii }->{ $tnt }, $tv );
		$hash->{ $pii }->{ $tnt }->{ pstd } = 1;
		$hash->{ $pii }->{ $tnt }->{ qstd } = 1;
		$cov->{ $pii } += $tv->{ cnt };
	    }
	}
	$sc3v->{ used } = 1;
	my %tins = ();
	$tins{ $bi }->{ "+$ins" } = $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt };
	realignins($hash, \%tins, $cov, $sclip5, $sclip3, $REF, $chr);
	adjCnt($hash->{ $bi }->{ I }->{ "+$ins" }, $hash->{$bi}->{ $REF->{$bi} }, $hash->{$bi}->{ $REF->{$bi} }) if ($rpflag && $bams && length($ins) >= 5 && length($ins) < $RLEN-10 && $hash->{$bi}->{ $REF->{$bi} } && $hash->{$bi}->{ $REF->{$bi} }->{ cnt } && noPassingReads($chr, $bi, $bi+length($ins), $bams) && $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } > 2 * $hash->{$bi}->{ $REF->{$bi} }->{ cnt });
    }
}

# Realign large deletions that are not present in alignment
sub realignlgdel {
    my ($hash, $cov, $sclip5, $sclip3, $REF, $chr) = @_;
    my $LONGMM = 3;
    my @tmp = ();
    while(my($p, $sc5v) = each %$sclip5) {
	push(@tmp, [$p, $sc5v, $sc5v->{cnt}]);
    }
    @tmp = sort {$b->[2] <=> $a->[2];} @tmp;
    foreach my $t (@tmp) {
	my ($p, $sc5v, $cnt) = @$t;
	next if ($sc5v->{ used }); # already been used in 
	my $seq = findconseq($sc5v);
	next unless( $seq );
	next if ($seq =~ /^.AAAAAAA/ || $seq =~ /^.TTTTTTT/ );
	next if ( length($seq) < 7 );
	next if ( islowcomplexseq($seq) );
	my $bp = findbp($seq, $p - 5, $REF, $INDELSIZE, -1, $chr);
	my $dellen = $p - $bp;
	next unless( $bp );
	my $extra = "";
	my $en = 0;
	while( substr($seq, $en, 1) ne $REF->{ $bp - $en - 1 } && $en < length($seq) ) {
	    $extra .= substr($seq, $en, 1);
	    $en++;
	}
	my $gt = -$dellen;
	if ( $extra ) {
	    $extra = reverse($extra);
	    $gt = "-$dellen&$extra";
	    $bp -= length($extra);
	}
	print STDERR "  Found Realignlgdel: $bp $gt 5' $p $seq $cnt\n" if ( $opt_y );
	$hash->{$bp}->{ $gt }->{ cnt } = 0 unless( $hash->{$bp}->{ $gt }->{ cnt } );
	$hash->{$bp}->{ $gt }->{ qstd } = 1; # more accurate implementation later
	$hash->{$bp}->{ $gt }->{ pstd } = 1; # more accurate implementation later
	for(my $tp = $bp; $tp < $bp + $dellen; $tp++) {
	    $cov->{$tp} += $sc5v->{ cnt };
	}
	adjCnt( $hash->{$bp}->{ $gt }, $sc5v );
	$sc5v->{ used } = 1;

	# Work on the softclipped read at 3'
	my $n = 0;
	$n++ while( $REF->{ $bp+$n } && $REF->{ $bp + $dellen + $n} && $REF->{ $bp+$n } eq $REF->{ $bp + $dellen + $n} );
	my $sc3p = $bp + $n;
	my $str = "";
	my $mcnt = 0;
	while( $REF->{ $bp+$n } && $REF->{ $bp+$dellen+$n } && $REF->{ $bp+$n } ne $REF->{ $bp + $dellen + $n} && $mcnt <= $LONGMM) {
	    $str .= $REF->{ $bp + $dellen + $n};
	    $n++; $mcnt++;
	}
	if ( length($str) == 1 ) {
	    $n++ while( $REF->{ $bp+$n } && $REF->{ $bp+$dellen+$n } && $REF->{ $bp+$n } eq $REF->{ $bp + $dellen + $n} );
	    $sc3p = $bp + $n;
	}
	if ( $sclip3->{ $sc3p } && (! $sclip3->{ $sc3p }->{ used }) ) {
	    $sc3p > $bp ? adjCnt($hash->{$bp}->{ $gt }, $sclip3->{ $sc3p }, $hash->{$bp}->{ $REF->{ $bp } }) : adjCnt($hash->{$bp}->{ $gt }, $sclip3->{ $sc3p });
	    if ( $sc3p == $bp ) {
		for(my $tp = $bp; $tp < $bp + $dellen; $tp++) {
		    $cov->{$tp} += $sclip3->{ $sc3p }->{ cnt };
		}
	    }
	    for(my $ip = $bp+1; $ip < $sc3p; $ip++) {
		rmCnt($hash->{$ip}->{ $REF->{$dellen + $ip} }, $sclip3->{$sc3p});
		delete $hash->{$ip}->{ $REF->{$dellen + $ip} } if ( $hash->{$ip}->{ $REF->{$dellen + $ip} } && $hash->{$ip}->{ $REF->{$dellen + $ip} }->{ cnt } == 0);
		delete $hash->{$ip} if ( (keys %{$hash->{$ip}} ) == 0);
	    }
	    $sclip3->{ $sc3p }->{ used } = 1;
	}
	my %dels5;
	$dels5{ $bp }->{$gt} = $hash->{$bp}->{ $gt }->{cnt};
	realigndel($hash, \%dels5, $cov, $sclip5, $sclip3, $REF, $chr);
	print STDERR "  Found lgdel done: $bp $gt $p 5' $seq $hash->{$bp}->{ $gt }->{cnt}\n\n" if ( $opt_y );
    }

    # Work on 3' clipped reads
    @tmp = ();
    while(my($p, $sc3v) = each %$sclip3) {
	push(@tmp, [$p, $sc3v, $sc3v->{cnt}]);
    }
    @tmp = sort {$b->[2] <=> $a->[2];} @tmp;
    foreach my $t (@tmp) {
	my ($p, $sc3v, $cnt) = @$t;
	#next if ($cnt < $MINR);
	#use Object; print STDERR Object::Perl($sc3v);
	next if ($sc3v->{ used }); # already been used in 
	my $seq = findconseq($sc3v);
	#print STDERR "LG3: $p $cnt $seq\t\n";
	next unless( $seq );
	next if ($seq =~ /^.AAAAAAA/ || $seq =~ /^.TTTTTTT/ );
	next if ( length($seq) < 7 );
	next if ( islowcomplexseq($seq) );
	my $bp = findbp($seq, $p + 5, $REF, $INDELSIZE, 1, $chr);
	my $dellen = $bp - $p;
	next unless( $bp );
	my $extra = "";
	my $en = 0;
	while( substr($seq, $en, 1) ne $REF->{ $bp + $en } && $en < length($seq) ) {
	    $extra .= substr($seq, $en, 1);
	    $en++;
	}
	my $gt = -$dellen;
	$bp = $p; # Set it to 5'
	if ( $extra ) {
	    #$gt = -($dellen - length($extra)+1) . "&$extra";
	    $gt = "-$dellen&$extra";
	} else { # 5' adjustment
	    $bp-- while($REF->{ $bp - 1 } eq $REF->{ $bp + $dellen - 1 });
	}
	print STDERR "  Found Realignlgdel: $bp $gt 3' $p $seq $cnt\n" if ( $opt_y );
	$hash->{$bp}->{ $gt }->{ cnt } = 0 unless( $hash->{$bp}->{ $gt }->{ cnt } );
	$hash->{$bp}->{ $gt }->{ qstd } = 1; # more accurate implementation later
	$hash->{$bp}->{ $gt }->{ pstd } = 1; # more accurate implementation later
	for(my $tp = $bp; $tp < $bp + $dellen; $tp++) {
	    $cov->{$tp} += $sc3v->{ cnt };
	}
	$sc3v->{ pmean } += $dellen*$sc3v->{ cnt };
	adjCnt( $hash->{$bp}->{ $gt }, $sc3v );
	$sc3v->{ used } = 1;

	my %dels5;
	$dels5{ $bp }->{$gt} = $hash->{$bp}->{ $gt }->{cnt};
	realigndel($hash, \%dels5, $cov, $sclip5, $sclip3, $REF, $chr);
	print STDERR "  Found lgdel: $bp $gt $p 3' $seq $hash->{$bp}->{ $gt }->{cnt}\n\n" if ( $opt_y );
    }
    print STDERR "  Done: Realignlgdel\n\n" if ( $opt_y );
}

sub subCnt {
    my ($vref, $dir, $rp, $q, $Q, $nm) = @_;  # ref dir read_position quality
    $vref->{ cnt }--;
    $vref->{ $dir }--;
    $vref->{ pmean } -= $rp;
    $vref->{ qmean } -= $q;
    $vref->{ Qmean } -= $Q;
    $vref->{ nm } -= $nm;
    $q >= $GOODQ ? $vref->{ hicnt }-- : $vref->{ locnt }--;
}

sub rmCnt {
    my ($vref, $tv) = @_;
    $vref->{ cnt } -= $tv->{ cnt };
    $vref->{ hicnt } -= $tv->{ hicnt } ? $tv->{ hicnt } : 0;
    $vref->{ locnt } -= $tv->{ locnt } ? $tv->{ locnt } : 0;
    $vref->{ pmean } -= $tv->{ pmean };
    $vref->{ qmean } -= $tv->{ qmean };
    $vref->{ Qmean } -= $tv->{ Qmean };
    $vref->{ "+" } -= $tv->{ "+" } ? $tv->{ "+" } : 0;
    $vref->{ "-" } -= $tv->{ "-" } ? $tv->{ "-" } : 0;
    foreach my $k (qw(cnt hicnt locnt pmean qmean Qmean + -)) {
	$vref->{ $k } = 0 if ( $vref->{ $k } && $vref->{ $k } < 0 );
    }
}

# Find the consensus sequence in soft-clipped reads.  Consensus is called if
# the matched nucleotides are >90% of all softly clipped nucleotides.
sub findconseq {
    my ($scv) = @_;
    return $scv->{ SEQ } if ( defined($scv->{ SEQ }) );
    my $total = 0;
    my $match = 0;
    my $seq = "";
    my $flag = 0;
    for(my $i = 0; $i < @{ $scv->{ nt } }; $i++ ) {
	my $nv = $scv->{ nt }->[$i];
	my $max = 0;
	my $maxq = 0;
	my $mnt = "";
	my $tt = 0;
	while(my ($nt, $ncnt) = each %$nv) {
	    $tt += $ncnt;
	    if ( $scv->{ seq }->[$i]->{ $nt }->{ qmean } > $maxq) {
		$max = $ncnt;
		$mnt = $nt;
		$maxq = $scv->{ seq }->[$i]->{ $nt }->{ qmean };
	    }
	    #print STDERR "Qmean: $maxq $nt $ncnt\n";
	}
	unless ( ($tt-$max <= 2 && $max > $tt - $max) || $max/$tt >= 0.8) {
	    last if ( $flag );
	    $flag = 1;
	}
	$total += $tt;
	$match += $max;
	$seq .= $mnt;
    }
    my $SEQ = ($total && $match/$total > 0.9 && length($seq)/1.5 > @{ $scv->{ nt } } - length($seq) && (length($seq)/(@{ $scv->{ nt } } +0) > 0.8 || @{ $scv->{ nt } } - length($seq) < 10 || length($seq) > 25)) ? $seq : "";
    print STDERR "  candidate consensus: $seq M: $match T: $total Final: $SEQ\n" if ( $opt_y );
    $scv->{ SEQ } = $SEQ;
    return $SEQ;
}

# Find the insertion
sub findbi {
    my ($seq, $p, $REF, $dir, $chr) = @_;
    my $MAXMM = 3; # maximum mismatches allowed
    my $score = 0;
    my ($BI, $INS, $BI2) = (0, "", 0);
    for(my $n = 6; $n < length($seq); $n++) {
	my $mm = 0;
	my $i = 0;
	my %m = ();
	last if ( $p + 6 >= $CHRS{ $chr } );
	for($i = 0; $i + $n < length($seq); $i++) {
	    last if ( $p + $dir*$i - ($dir == -1 ? 1 : 0) < 1 );
	    last if ( $p + $dir*$i - ($dir == -1 ? 1 : 0) > $CHRS{ $chr } );
	    if ( substr($seq, $i + $n, 1) ne $REF->{ $p + $dir*$i - ($dir == -1 ? 1 : 0) } ) {
		$mm++; 
	    } else {
		$m{ substr($seq, $i + $n, 1) }++;
	    }
	    last if ( $mm > $MAXMM );
	}
	my @mnt = keys %m;
	next unless( @mnt >= 2 ); # at least three different NT for overhang sequences, weeding out low complexity seq
	#print STDERR "bi: $n $i ", substr($seq, $n), " $p $seq $mm\n";
	if ( (@mnt >= 3 && $i + $n >= length($seq) - 1 && $i >= 8  && $mm/$i < 0.15) || (@mnt >= 2 && $mm == 0 && $i + $n == length($seq) && $n >= 20 && $i >= 8)) {
	    my $ins = substr($seq, 0, $n);
	    my $extra = "";
	    my $ept = 0;
	    while(substr($seq, $n + $ept, 1) ne $REF->{ $p + $ept*$dir - ($dir == -1 ? 1 : 0) } || substr($seq, $n + $ept + 1, 1) ne $REF->{ $p + ($ept+1)*$dir - ($dir == -1 ? 1 : 0) }) {
		$extra .= substr($seq, $n + $ept, 1);
		$ept++;
	    }
#	print STDERR "bi: $n $i ", substr($seq, $n), " $p $seq $extra $ept $mm\n";
	    if ($dir == -1) {
		$ins .= $extra;
		$ins = reverse($ins);
		if ( $extra ) {
		    substr($ins, -length($extra), 0) = "&";
		}
		if ( $mm == 0 && $i + $n == length($seq)) {
		    ($BI, $INS, $BI2) = ($p-1-length($extra), $ins, $p-1);
		    ($BI, $INS, $BI2) = adjInsPos($BI, $INS, $REF) unless( $extra );
		    return ($BI, $INS, $BI2);
		} elsif ($i-$mm > $score) {
		    ($BI, $INS, $BI2) = ($p-1-length($extra), $ins, $p-1);
		    $score = $i-$mm;
		}
	    } else {
		my $s = -1;
		if ( $extra ) {
		    $ins .= "&$extra";
		} else {
		    while( substr($ins, $s, 1) eq $REF->{ $p + $s } && $s >= -$n ) {
			$s--;
		    }
		    if ( $s < -1 ) {
			my $tins = substr($ins, $s + 1, 1 - $s);
			substr($ins, $s+1) = "";
			substr($ins, 0, 0) = $tins;
		    }
		}
		if( $mm == 0 && $i + $n == length($seq)) {
		    ($BI, $INS, $BI2) = ($p+$s, $ins, $p+$s+length($extra));
		    ($BI, $INS, $BI2) = adjInsPos($BI, $INS, $REF) unless( $extra );
		    return ($BI, $INS, $BI2);
		} elsif ($i-$mm > $score) {
		    ($BI, $INS, $BI2) = ($p+$s, $ins, $p+$s+length($extra));
		    $score = $i-$mm;
		}
	    }
	}
    }
    unless( $BI2 != $BI ) {
	($BI, $INS) = adjInsPos( $BI, $INS, $REF ) if ( $INS && $BI);
    }
    return ($BI, $INS, $BI2);
}

# Adjust the insertion position if necessary
sub adjInsPos {
    my ($bi, $ins, $REF) = @_;
    my $n = 1;
    my $len = length($ins);
    while( $REF->{ $bi } eq substr($ins, -$n, 1) ) {
	$n++;
	$n = 1 if ( $n > $len );
	$bi--;
    }
    $ins = substr($ins, -($n-1)) . substr($ins, 0, -($n-1)) if ( $n > 1 );
    return( $bi, $ins, $bi );
}

# Find breakpoint
sub findbp {
    my ($seq, $sp, $REF, $dis, $dir, $chr) = @_;
    my $MAXMM = 3; # maximum mismatches allowed
    my $BP = 0;
    my $score = 0;
    for(my $n = 0; $n < $dis; $n++) {
	my $mm = 0;
	my %m = ();
	my $i = 0;
	for($i = 0; $i < length($seq); $i++) {
	    last if ( $sp + $dir*$n + $dir*$i < 1 );
	    last if ( $sp + $dir*$n + $dir*$i > $CHRS{ $chr } );
	    if ( substr($seq, $i, 1) ne $REF->{ $sp + $dir*$n + $dir*$i } ) {
		$mm++; 
	    } else {
		$m{ substr($seq, $i, 1) }++;
	    }
	    last if ( $mm >  $MAXMM - int($n/100) );
	}
	my @mnt = keys %m;
	next unless( @mnt >= 3 );
	if ( $mm <= $MAXMM - int($n/100) && $i >= length($seq) - 2 && $i >= 8 + int($n/10) && $mm/$i < 0.12) {
	    my $bp = $sp + $dir*$n - ($dir < 0 ? $dir : 0);
	    if ( $mm == 0 && $i == length($seq) ) {
		print STDERR "  Findbp: $seq $sp $bp $mm $i\n" if ( $opt_y );
		return $bp;
	    } elsif ($i-$mm > $score) {
		$BP = $bp;
		$score = $i-$mm;
	    }
	}
    }
    print STDERR "  Findbp with mismatches: $seq $sp $BP $dir $score\n" if ( $opt_y && $BP );
    return $BP;
}

# Find breakpoint for complex variants
sub findbp2 {
    my ($seq, $sp, $REF, $dis, $dir, $chr) = @_;
    my $MAXMM = 3; # maximum mismatches allowed
    my $BP = 0;
    my $BI = 0;
    my $SI = 0;
    my $score = 0;
    for(my $n = 0; $n < $dis; $n++) {
	my $mm = 0;
	my %m = ();
	my $i = 0;
	my $pmm = 0;
	for($i = length($seq) - 1 ; $i >= 0; $i--) {
	    next if ( $sp + $dir*$n + $dir*$i < 1 );
	    next if ( $sp + $dir*$n + $dir*$i > $CHRS{ $chr } );
	    if ( substr($seq, $i, 1) ne $REF->{ $sp + $dir*$n + $dir*$i } ) {
		$mm++;
		$pmm++;
		if ( $pmm > 1 ) { # two consecutive mismatches
		    $mm -= 2;
		    $i += 2;
		    last;
		}
	    } else {
		$pmm = 0;
		$m{ substr($seq, $i, 1) }++;
	    }
	    last if ( $mm >  $MAXMM - int($n/100) );
	}
	my @mnt = keys %m;
	next unless( @mnt >= 3 );
	#print STDERR "$mm $dir $n $sp $i\n";
	if ( $mm <= $MAXMM - int($n/100) && $i <= 2 && length($seq) - $i >= 8 && $mm/(length($seq) - $i) < 0.12) {
	    #print STDERR "$mm $dir $n $sp $i $seq\n";
	    my $bp = $sp + $dir*$n + ($dir < 0 ? -$i*$dir : $i);
	    if ( $mm == 0 && $i == 0 ) {
		return ($bp, 0);
	    } elsif (length($seq) - $i - $mm > $score) {
		$BP = $bp;
		$BI = $i;
		$score = $i-$mm;
	    }
	}
    }
    return ($BP, $BI);
}

# Realign deletions if already present in alignment
sub realigndel {
    my ($hash, $dels5, $cov, $sclip5, $sclip3, $REF, $chr, $bams) = @_;
    my $LONGMM = 3; # Longest continued mismatches typical aligned at the end
    my @tmp = ();
    while(my($p, $dv) = each %$dels5) {
	while(my($vn, $dcnt) = each %$dv) {
	    my $ecnt = 0;
	    $ecnt = length($1) if ( $vn =~ /([ATGC&]+)$/ );
	    push(@tmp, [$p, $vn, $dcnt, $ecnt]);
	}
    }
    #@tmp = sort {$b->[2] - $b->[3] <=> $a->[2] - $a->[3]} @tmp;
    @tmp = sort {$b->[2] <=> $a->[2]} @tmp;
    foreach my $tmpv (@tmp) {
	my ($p, $vn, $dcnt) = @$tmpv;
	print STDERR "  Realigndel for: $p $vn $dcnt cov: $cov->{ $p }\n" if ( $opt_y );
	my $vref = $hash->{ $p }->{ $vn };
	my $dellen = $vn =~ /^-(\d+)/ ? $1 : 0;
	$dellen += $1 if ( $vn =~ /\^(\d+)$/ );
	my $extrains = ($vn =~ /\^([ATGNC]+)/) ? $1 : "";
	my $extra = ($vn =~ /^-\d+(.*)/) ? $1 : ""; $extra =~ s/\^|&|#//g;
	my $wustart = ($p - $dellen - 100) > 1 ? ($p - $dellen - 100) : 1;
	my $wupseq = join( "", (map { $REF->{ $_ }; } ($wustart .. ($p-1)))) . $extra; # . $extra; # 5' flanking seq
	my $sanend = ($p + 2*$dellen + 100) > $CHRS{$chr} ? $CHRS{ $chr } : ($p + 2*$dellen + 100);
	my $sanpseq = $extra . join( "", (map { $REF->{ $_ }; } (($p + $dellen + length($extra) - length($extrains)) .. $sanend))); # 3' flanking seq
	my ($mm3, $sc3p, $nm3, $misp3, $misnt3) = findMM3($REF, $p, $sanpseq, $dellen, $sclip3); # mismatches, mismatch positions, 5 or 3 ends
	my ($mm5, $sc5p, $nm5, $misp5, $misnt5) = findMM5($REF, $p+$dellen+length($extra)-length($extrains)-1, $wupseq, $dellen, $sclip5);
	my @mm = (@$mm3, @$mm5);
	print STDERR "  Mismatches: misp3: $misp3-$misnt3 misp5: $misp5-$misnt5 sclip3: @$sc3p sclip5: @$sc5p\n" if ( $opt_y );
	for(my $mi = 0; $mi < @mm; $mi++) {
	    my ($mm, $mp, $me) = @{$mm[$mi]};
	    substr($mm, 1, 0) = "&" if (length($mm) > 1);
	    next unless( $hash->{ $mp } );
	    next unless( $hash->{ $mp }->{ $mm } );
	    my $tv = $hash->{ $mp }->{ $mm };
	    next unless( $tv->{ cnt } );
	    next if( $tv->{ qmean }/$tv->{ cnt } < $GOODQ );
	    next unless ( $tv->{ pmean }/$tv->{ cnt } <= ($me == 3 ? $nm3 + 4 : $nm5 + 4)); # || ($tv->{ pmean }/$tv->{ cnt } < ($sc3p[1] ? $sc3p[1] : $sc3p[0]) - $mp && $me == 3); #$opt_k;
	    next unless( $tv->{ cnt } < $dcnt + $dellen && $tv->{ cnt } / $dcnt < 8);
	    print STDERR "  Realigndel Adj: $mm $mp $me $nm3 $nm5 $p $tv->{ cnt } $tv->{ qmean } cov: $cov->{ $p }\n" if ( $opt_y );
	    # Adjust ref cnt so that AF won't > 1
	    if ( $mp > $p && $me == 5 ) {
		my $f = $tv->{pmean} ? ($mp-$p)/($tv->{pmean}/$tv->{cnt}) : 1;
		$f = 1 if ( $f > 1 );
		$cov->{ $p } += int($tv->{ cnt } * $f);
		adjRefCnt($tv, $hash->{$p}->{ $REF->{$p} }, $dellen);
	    }

	    my $ref = ($mp > $p && $me == 3) ? ($hash->{$p}->{ $REF->{$p} } ? $hash->{$p}->{ $REF->{$p} } : "") : "";
	    adjCnt($vref, $tv, $ref);
	    delete $hash->{ $mp }->{ $mm };
	    delete $hash->{ $mp } if ( (keys %{$hash->{ $mp }} < 1) );
	    print STDERR "  Realigndel AdjA: $mm $mp $me $nm3 $nm5 $p $tv->{ cnt } $tv->{ qmean } cov: $cov->{ $p }\n" if ( $opt_y );
	}
	delete $hash->{ $misp3 }->{ $misnt3 } if ( $misp3 && $hash->{ $misp3 }->{ $misnt3 } && @$mm3 == 1);
	delete $hash->{ $misp5 }->{ $misnt5 } if ( $misp5 && $hash->{ $misp5 }->{ $misnt5 } && @$mm5 == 1);

	#next if ( $dellen < 2 ); # soft-clipping only happens when deletion is at least 2 bp
	foreach my $sc5pp (@$sc5p) {
	    if ( $sclip5->{ $sc5pp } && (! $sclip5->{ $sc5pp }->{ used }) ) {
		my $tv = $sclip5->{ $sc5pp };
		my $seq = findconseq( $tv );
		print STDERR "  Realigndel 5: $p $sc5pp Seq: '$seq' ", scalar reverse($wupseq), " $tv->{ cnt } $dcnt $vn $p cov: $cov->{ $p }\n" if ( $opt_y );
		if ( $seq && ismatch($seq, $wupseq, -1) ) {
		    $cov->{ $p } += $tv->{ cnt } if ( $sc5pp > $p );
		    adjCnt($vref, $tv);
		    $sclip5->{ $sc5pp }->{ used } = 1;
		    print STDERR "  Realigndel 5: $p $sc5pp $seq ", scalar reverse($wupseq), " $tv->{ cnt } $dcnt $vn $p used cov: $cov->{ $p }\n" if ( $opt_y && $seq );
		}
	    }
	}
	#use Object; print STDERR Object::Perl($sclip3);
	#print STDERR keys(%$sclip3), "\n";
	foreach my $sc3pp (@$sc3p) {
	    if ( $sclip3->{ $sc3pp } && (! $sclip3->{ $sc3pp }->{ used }) ) {
		my $tv = $sclip3->{ $sc3pp };
		my $seq = findconseq( $tv );
		#if ( $seq && findbp($seq, $sc3pp + $dellen, $REF, 1, 1, $chr) ) 
		#print STDERR "$seq $sanpseq $sc3pp $p\n";
		print STDERR "  Realigndel 3: $p $sc3pp seq '$seq' $sanpseq $tv->{ cnt } $dcnt $vn $p $dellen ", substr($sanpseq, $sc3pp-$p), "\n" if ( $opt_y );
		if ( $seq && ismatch($seq, substr($sanpseq, $sc3pp-$p), 1) ) {
		    print STDERR "  Realigndel 3: $p $sc3pp $seq $sanpseq $tv->{ cnt } $dcnt $vn $p used\n" if ( $opt_y && $seq );
		    $cov->{ $p } += $tv->{ cnt } if ( $sc3pp <= $p );
		    my $ref = $sc3pp <= $p ? "" : $hash->{$p}->{ $REF->{$p} };
		    adjCnt($vref, $tv, $ref);
		    $sclip3->{ $sc3pp }->{ used } = 1;
		}
	    }
	}
	#my $pe = $p + $dellen + length($extra) + length($compm);
	my $pe = $p + $dellen + length($extra) - length($extrains);
	adjCnt($vref, $hash->{$p}->{ $REF->{$p} }, $hash->{$p}->{ $REF->{$p} }) if ($bams && $pe - $p >= 5 && $pe - $p < $RLEN-10 && $hash->{$p}->{ $REF->{$p} } &&  $hash->{$p}->{ $REF->{$p} }->{ cnt } && noPassingReads($chr, $p, $pe, $bams) && $vref->{ cnt } > 2 * $hash->{$p}->{ $REF->{$p} }->{ cnt } * (1 - ($pe-$p)/$RLEN));  # taking the size of gap into account
    }
    for(my $i = $#tmp; $i > 0; $i--) {
	my ($p, $vn, $icnt) = @{$tmp[$i]};
	next unless ($hash->{ $p }->{ $vn });
	my $vref = $hash->{ $p }->{ $vn };
	if ($vn =~ /(-\d+)&[ATGC]+$/ ) {
	    my $tn = $1;
	    if ( $hash->{ $p }->{ $tn } ) {
		my $tref = $hash->{ $p }->{ $tn };
		if ( $vref->{ cnt } < $tref->{ cnt } ) {
		    adjCnt($tref, $vref);
		    delete $hash->{ $p }->{ $vn };
		}
	    }
	}
    }
}

# check whether there're reads supporting wild type in deletions
# Only for indels that have micro-homology
sub noPassingReads {
    my ($chr, $s, $e, $bams) = @_;
    my $cnt = 0;
    my $midcnt = 0; # Reads end in the middle
    my $dlen = $e - $s;
    my $dlenqr = qr/${dlen}D/;
    foreach my $bam (@$bams) {
	open(SMB, "samtools view $bam $chr:$s-$e |");
	while(<SMB>) {
	    my @a = split(/\t/);
	    my $rs = $a[3];
	    my $rlen = 0; $rlen += $1 while( $a[5] =~ /(\d+)[MD]/g ); # The total aligned length, excluding soft-clipped bases and insertions
	    next if ( $a[5] =~ /$dlenqr/ );
	    my $re = $rs + $rlen;
	    $cnt++ if ( $re > $e+2 && $rs < $s-2 );
	    $midcnt++ if ( $rs < $s-2 && $re > $s && $re < $e );
	}
	close(SMB);
    }
    print STDERR "    Passing Read CNT: $cnt $chr $s $e $midcnt\n" if ( $opt_y );
    return $cnt > 0 ? 0 : $midcnt + 1;
}

# Adjust the count,  If $ref is given, the count for reference is also adjusted.
sub adjCnt {
    my ($vref, $tv, $ref) = @_;
    $vref->{ cnt } += $tv->{ cnt };
    $vref->{ extracnt } += $tv->{ cnt };
    $vref->{ hicnt } += $tv->{ hicnt } ? $tv->{ hicnt } : 0;
    $vref->{ locnt } += $tv->{ locnt } ? $tv->{ locnt } : 0;
    $vref->{ pmean } += $tv->{ pmean };
    $vref->{ qmean } += $tv->{ qmean };
    $vref->{ Qmean } += $tv->{ Qmean };
    $vref->{ nm } += $tv->{ nm };
    $vref->{ pstd } = 1;
    $vref->{ qstd } = 1;
    $vref->{ "+" } += $tv->{ "+" } ? $tv->{ "+" } : 0;
    $vref->{ "-" } += $tv->{ "-" } ? $tv->{ "-" } : 0;
    print STDERR "    AdjCnt: '+' $vref->{ cnt } $tv->{ cnt } ", $vref->{ "+" }, " ", $tv->{ "+" } ? $tv->{ "+" } : 0, " Ref: ", $ref && $ref->{ cnt } ? $ref->{ cnt } : "NA", "\n" if ( $opt_y );
    print STDERR "    AdjCnt: '-' $vref->{ cnt } $tv->{ cnt } ", $vref->{ "-" }, " ", $tv->{ "-" } ? $tv->{ "-" } : 0, " Ref: ", $ref && $ref->{ cnt } ? $ref->{ cnt } : "NA", "\n" if ( $opt_y );
    return unless($ref);
    $ref->{ cnt } -= $tv->{ cnt };
    $ref->{ hicnt } -= $tv->{ hicnt } ? $tv->{ hicnt } : 0;
    $ref->{ locnt } -= $tv->{ locnt } ? $tv->{ locnt } : 0;
    $ref->{ pmean } -= $tv->{ pmean };
    $ref->{ qmean } -= $tv->{ qmean };
    $ref->{ Qmean } -= $tv->{ Qmean };
    $ref->{ nm } -= $tv->{ nm };
    $ref->{ "+" } -= $tv->{ "+" } ? $tv->{ "+" } : 0;
    $ref->{ "-" } -= $tv->{ "-" } ? $tv->{ "-" } : 0;
    foreach my $k (qw(cnt hicnt locnt pmean qmean Qmean + -)) {
	$ref->{ $k } = 0 if ( $ref->{ $k } && $ref->{ $k } < 0 );
    }
}

# Adjust the reference count.
sub adjRefCnt {
    my ($tv, $ref, $len) = @_;
    return unless($ref);
    my $f = $tv->{ pmean } ? ($tv->{ pmean }/$tv->{ cnt }-$len+1)/($tv->{ pmean }/$tv->{ cnt }) : 0; # the adjustment factor
    $f = 1 if ( $f > 1 );
    return if ( $f < 0 );
    $ref->{ cnt } -= int($f*$tv->{ cnt });
    $ref->{ hicnt } -= $tv->{ hicnt } ? int($f*$tv->{ hicnt }) : 0;
    $ref->{ locnt } -= $tv->{ locnt } ? int($f*$tv->{ locnt }) : 0;
    $ref->{ pmean } -= $f * $tv->{ pmean };
    $ref->{ qmean } -= $f * $tv->{ qmean };
    $ref->{ Qmean } -= $f * $tv->{ Qmean };
    $ref->{ nm } -= $f * $tv->{ nm };
    $ref->{ "+" } -= $tv->{ "+" } ? int($f * $tv->{ "+" }) : 0;
    $ref->{ "-" } -= $tv->{ "-" } ? int($f * $tv->{ "-" }) : 0;
    foreach my $k (qw(cnt hicnt locnt pmean qmean Qmean + -)) {
	$ref->{ $k } = 0 if ( $ref->{ $k } && $ref->{ $k } < 0 );
    }
}

# Given a variant sequence, find the mismatches and potential softclipping positions
sub findMM3 {
    my ($REF, $p, $seq, $len, $sclip3) = @_;
    $seq =~ s/#|\^//g;
    my $LONGMM = 3;
    my @mm = (); # mismatches, mismatch positions, 5 or 3 ends
    my $n = 0;
    my $mn = 0;
    my $mcnt = 0;
    my $str = "";
    my @sc3p = ();
    $n++ while( $REF->{ $p+$n } && $REF->{ $p+$n } eq substr($seq, $n, 1) );
    push(@sc3p, $p + $n);
    my $Tbp = $p + $n;
    while( $REF->{ $p+$n } ne substr($seq, $n, 1) && $mcnt <= $LONGMM && $n < length($seq)) {
	$str .= substr($seq, $n, 1);
	push(@mm, [$str, $Tbp, 3]);
	$n++; $mcnt++;
    }
    # Adjust clipping position if only one mismatch
    my ($misp, $misnt) = (0, "");
    if ( length($str) == 1 ) {
	$n++ && $mn++ while( $REF->{ $p+$n } && $REF->{ $p+$n } eq substr($seq, $n, 1) );
	if ( $mn > 1 ) {
	    my $n2 = 0;
	    $n2++ while( $REF->{ $p+$n+1+$n2 } && $n+$n2+1 < length($seq) && $REF->{ $p+$n+1+$n2 } eq substr($seq, $n+$n2+1, 1) );
	    if ( $n2 > 2 && $n+$n2+1 < length($seq) ) {
		push(@sc3p, $p + $n + $n2);
		($misp, $misnt) = ($p + $n, substr($seq, $n, 1));
		$sclip3->{ $p+$n+$n2 }->{ used } = 1 if ( $sclip3->{ $p+$n+$n2 } ); #( && $len >= 4);
		$mn += $n2;
	    } else {
		push(@sc3p, $p + $n);
		$sclip3->{ $p+$n }->{ used } = 1 if ( $sclip3->{ $p+$n } ); #( && $len >= 4);
	    }
	}
    }
    #print STDERR "MM3: $seq $len $p '@sc3p'\n";
    return (\@mm, \@sc3p, $mn, $misp, $misnt);
}

sub findMM5 {
    my ($REF, $p, $seq, $len, $sclip5) = @_;
    $seq =~ s/#|\^//g;
    my $LONGMM = 3;
    my @mm = (); # mismatches, mismatch positions, 5 or 3 ends
    my $n = 0;
    my $mn = 0;
    my $mcnt = 0;
    my $str = "";
    my @sc5p = ();
    while( $REF->{ $p-$n } ne substr($seq, -1-$n, 1) && $mcnt < $LONGMM ) {
	$str = substr($seq, -1-$n, 1) . $str;
	push(@mm, [$str, $p-$n, 5]);
	$n++; $mcnt++;
    }
    push(@sc5p, $p + 1);
    my $Tbp = $p+1;
    # Adjust clipping position if only one mismatch
    my ($misp, $misnt) = (0, "");
    if ( length($str) == 1 ) {
	$n++ && $mn++ while( $REF->{ $p-$n } && $REF->{ $p-$n } eq substr($seq, -1-$n, 1) );
	if ( $mn > 1 ) {
	    my $n2 = 0;
	    $n2++ while( $REF->{ $p-$n-1-$n2 } && -1-$n-1-$n2 >= 0 && $REF->{ $p-$n-1-$n2 } eq substr($seq, -1-$n-1-$n2, 1) );
	    if ( $n2 > 2 ) {
		push(@sc5p, $p - $n - $n2);
		($misp, $misnt) = ($p - $n, substr($seq, -1-$n, 1));
		$sclip5->{ $p-$n-$n2 }->{ used } = 1 if ( $sclip5->{ $p-$n-$n2 } ); #(&& $len >= 4);
		$mn += $n2;
	    } else {
		push(@sc5p, $p - $n);
		$sclip5->{ $p-$n }->{ used } = 1 if ( $sclip5->{ $p-$n } ); #(&& $len >= 4);
	    }
	}
    }
    return (\@mm, \@sc5p, $mn, $misp, $misnt);
}

sub realignins {
    my ($hash, $ins, $cov, $sclip5, $sclip3, $REF, $chr) = @_;
    my @tmp = ();
    while(my($p, $iv) = each %$ins) {
	while(my($vn, $icnt) = each %$iv) {
	    my $ecnt = 0;
	    $ecnt = length($1) if ( $vn =~ /([ATGC&]+)$/ );
	    push(@tmp, [$p, $vn, $icnt, $ecnt]);
	}
    }
    #@tmp = sort {$b->[2] - $b->[3] <=> $a->[2] - $a->[3]} @tmp;
    @tmp = sort {$b->[2] <=> $a->[2]} @tmp;
    foreach my $tmpv (@tmp) {
	my ($p, $vn, $icnt) = @$tmpv;
	print STDERR "  Realign Ins: $p $vn $icnt\n" if ( $opt_y );
	my $vref = $hash->{ $p }->{ I }->{ $vn };
	my $ins = $1 if ( $vn =~ /^\+([ATGC]+)/ );
	next unless($ins);
	my $extra = ($vn =~ /&([ATGC]+)/) ? $1 : "";
	my $compm = ($vn =~ /#([ATGC]+)/) ? $1 : "";  # the match part for a complex variant
	my $newins = ($vn =~ /\^([ATGC]+)$/) ? $1 : "";  # the adjacent insertion
	my $newdel = ($vn =~ /\^(\d+)$/) ? $1 : 0;  # the adjacent deletion
	my $tn = $vn;  $tn =~ s/^\+//; $tn =~ s/&//; $tn =~ s/#//; $tn =~ s/\^\d+$//; $tn =~ s/\^//;
	my $wustart = ($p - 100 - length($vn)+1) > 1 ? ($p - 100 - length($vn)+1) : 1;
	my $wupseq = join( "", (map { $REF->{ $_ }; } ($wustart .. $p))) . $tn; # 5prime flanking seq
	my $sanend = $CHRS{ $chr } < ($p + length($vn) + 100) ? $CHRS{ $chr } : ($p + length($vn) + 100);
	my $sanpseq = $tn . join( "", (map { $REF->{ $_ }; } (($p+length($extra)+1+length($compm)+$newdel) .. $sanend))); # 3prime flanking seq
	my ($mm3, $sc3p, $nm3, $misp3, $misnt3) = findMM3($REF, $p+1, $sanpseq, length($ins)+length($compm), $sclip3); # mismatches, mismatch positions, 5 or 3 ends
	my ($mm5, $sc5p, $nm5, $misp5, $misnt5) = findMM5($REF, $p+length($extra)+length($compm)+$newdel, $wupseq, length($ins)+length($compm), $sclip5);
	my @mm = (@$mm3, @$mm5);
	for(my $mi = 0; $mi < @mm; $mi++) {
	    my ($mm, $mp, $me) = @{$mm[$mi]};
	    substr($mm, 1, 0) = "&" if (length($mm) > 1);
	    next unless( $hash->{ $mp } );
	    next unless( $hash->{ $mp }->{ $mm } );
	    my $tv = $hash->{ $mp }->{ $mm };
	    next unless( $tv->{ cnt } );
	    next if( $tv->{ qmean }/$tv->{ cnt } < $GOODQ );
	    next if ( $tv->{ pmean }/$tv->{ cnt } > ($me == 3 ? $nm3 + 4 : $nm5 + 4) ); #$opt_k;
	    next unless( $tv->{ cnt } < $icnt+length($ins) && $tv->{ cnt }/$icnt < 8);
	    print STDERR "    insMM: $mm\t$mp\t$me\t$nm3\t$nm5\t$vn\t$icnt\t$tv->{cnt}\t$tv->{qmean}\t$tv->{pmean}\t$cov->{ $p }\n" if ( $opt_y );
	    # Adjust ref cnt so that AF won't > 1
	    if ( $mp > $p && $me == 5 ) {
		$cov->{ $p } += $tv->{ cnt };
	    }

	    my $ref = ($mp > $p && $me == 3) ? ($hash->{$p}->{ $REF->{$p} } ? $hash->{$p}->{ $REF->{$p} } : "") : "";
	    adjCnt($vref, $tv, $ref);
	    delete $hash->{ $mp }->{ $mm };
	    delete $hash->{ $mp } if ( (keys %{$hash->{ $mp }} < 1) );
	}
	delete $hash->{ $misp3 }->{ $misnt3 } if ( $misp3 && $hash->{ $misp3 }->{ $misnt3 } && @$mm3 == 1);
	delete $hash->{ $misp5 }->{ $misnt5 } if ( $misp5 && $hash->{ $misp5 }->{ $misnt5 } && @$mm5 == 1);

	foreach my $sc5pp (@{ $sc5p }) {
	    if ( $sclip5->{ $sc5pp } && (! $sclip5->{ $sc5pp }->{ used }) ) {
		my $tv = $sclip5->{ $sc5pp };
		my $seq = findconseq( $tv );
		print STDERR "    ins5: $p $sc5pp $seq $wupseq $icnt $tv->{ cnt }\n" if ( $opt_y );
		if ( $seq && ismatch($seq, $wupseq, -1) ) {
		    print STDERR "    ins5: $p $sc5pp $seq $wupseq $icnt $tv->{ cnt } used\n" if ( $opt_y );
		    $cov->{ $p } += $tv->{ cnt } if ( $sc5pp > $p );
		    adjCnt($vref, $tv);
		    $sclip5->{ $sc5pp }->{ used } = 1;
		}
	    }
	}

	foreach my $sc3pp (@$sc3p) {
	    print STDERR "    33: $p $sc3pp $vn $sanpseq\n" if ( $opt_y );
	    if ( $sclip3->{ $sc3pp } && (! $sclip3->{ $sc3pp }->{ used }) ) {
		my $tv = $sclip3->{ $sc3pp };
		my $seq = findconseq( $tv );
		print STDERR "    ins3: $p $sc3pp $seq $sanpseq $vn $icnt $tv->{cnt}\n" if ($opt_y);
		if ( $seq && ismatch($seq, substr($sanpseq, $sc3pp - $p - 1), 1) ) {
		    print STDERR "    ins3: $p $sc3pp $seq $vn $icnt $tv->{cnt} used\n" if ( $opt_y );
		    $cov->{ $p } += $tv->{ cnt } if ( $sc3pp <= $p );
		    my $ref = $sc3pp <= $p ? "" : $hash->{$p}->{ $REF->{$p} };
		    adjCnt($vref, $tv, $ref);
		    $sclip3->{ $sc3pp }->{ used } = 1;
		}
	    }
	}
    }
    for(my $i = $#tmp; $i > 0; $i--) {
	my ($p, $vn, $icnt) = @{$tmp[$i]};
	next unless ($hash->{ $p }->{ I }->{ $vn });
	my $vref = $hash->{ $p }->{ I }->{ $vn };
	if ($vn =~ /(\+[ATGC]+)&[ATGC]+$/ ) {
	    my $tn = $1;
	    if ( $hash->{ $p }->{ I }->{ $tn } ) {
		my $tref = $hash->{ $p }->{ I }->{ $tn };
		if ( $vref->{ cnt } < $tref->{ cnt } ) {
		    adjCnt($tref, $vref, $hash->{ $p }->{ $REF->{ $p } });
		    delete $hash->{ $p }->{ I }->{ $vn };
		}
	    }
	}
    }
}

sub ismatch {
    my ($seq1, $seq2, $dir) = @_;
    print STDERR "    Matching $seq1 $seq2 $dir\n" if ( $opt_y );
    $seq2 =~ s/#|\^//g;
    my ($mm, $n) = (0, 0);
    my %nts = ();
    for(my $n = 0; $n < length($seq1) && $n < length($seq2); $n++) {
	$nts{ substr($seq1, $n, 1) } = 1;
	$mm++ if ( substr($seq1, $n, 1) ne substr($seq2, $dir*$n - ($dir == -1 ? 1 : 0), 1) );
    }
    my @nts = keys %nts;
    return ($mm <= 2 && $mm/length($seq1) < 0.15) ? 1 : 0;
}

sub ismatchref {
    my ($seq, $REF, $p, $dir) = @_;
    print STDERR "      Matching REF $seq $p $dir\n" if ( $opt_y );
    my ($mm, $n) = (0, 0);
    for(my $n = 0; $n < length($seq); $n++) {
        $mm++ if ( substr($seq, $dir == 1 ? $n : $dir * $n - 1, 1) ne $REF->{ $p + $dir*$n } );
    }
    return ($mm <= 3 && $mm/length($seq) < 0.15) ? 1 : 0;
}

sub strandBias {
    my ($fwd, $rev) = @_;
    if ( $fwd + $rev <= 12 ) {  # using p=0.01, because prop.test(1,12) = 0.01
	return $fwd*$rev > 0 ? 2 : 0;
    }
    return ($fwd/($fwd+$rev) >= $BIAS && $rev/($fwd+$rev) >= $BIAS && $fwd >= $MINB && $rev >= $MINB) ? 2 : 1;
}

#FCB02N4ACXX:3:2206:20108:2526#GATGGTTC  163     chr3    38181981	50      79M188N11M      =       38182275	667     TGAAGTTGTGTGTGTCTGACCGCGATGTCCTGCCTGGCACCTGTGTCTGGTCTATTGCTAGTGAGCTCATCGTAAAGAGGTGCCGCCGGG \YY`c`\ZQPJ`e`b]e_Sbabc[^Ybfaega_^cafhR[U^ee[ec][R\Z\__ZZbZ\_\`Z`d^`Zb]bBBBBBBBBBBBBBBBBBB AS:i:-8 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:72A16A0    YT:Z:UU XS:A:+  NH:i:1  RG:Z:15
sub USAGE {
    print STDERR <<USAGE;
    $0 [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] [-x #_nu] [-g gene] [-f freq] [-r #_reads] [-B #_reads] region_info

    VarDict is a variant calling program for SNV, MNV, indels (<120 bp default, but can be set using -I option), and complex variants.  It accepts any BAM format, either
    from DNA-seq or RNA-seq.  There're several distinct features over other variant callers.  First, it can perform local
    realignment over indels on the fly for more accurate allele frequencies of indels.  Second, it rescues softly clipped reads
    to identify indels not present in the alignments or support existing indels.  Third, when given the PCR amplicon information,
    it'll perform amplicon-based variant calling and filter out variants that show amplicon bias, a common false positive in PCR
    based targeted deep sequencing.  Forth, it has very efficient memory management and memory usage is linear to the region of
    interest, not the depth.  Five, it can handle ultra-deep sequencing and the performance is only linear to the depth.  It has
    been tested on depth over 2M reads.  Finally, it has a build-in capability to perform paired sample analysis, intended for
    somatic mutation identification, comparing DNA-seq and RNA-seq, or resistant vs sensitive in cancer research.  By default,
    the region_info is an entry of refGene.txt from IGV, but can be any region or bed files.

    -H Print this help page
    -h Print a header row decribing columns
    -v VCF format output
    -i Output splicing read counts
    -p Do pileup regarless the frequency
    -C Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
    -D Debug mode.  Will print some error messages and append full genotype at the end.
    -M Similar to -D, but will append individual quality and position data instead of mean
    -t Indicate to remove duplicated reads.  Only one pair with same start positions will be kept
    -3 Indicate to move indels to 3-prime if alternative alignment can be achieved.
    -F bit
       The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and duplicates).  Use -F 0 to turn it off.
    -z 0/1
       Indicate wehther is zero-based cooridates, as IGV does.  Default: 1 for BED file or amplicon BED file.  Use 0 to turn it off.
       When use -R option, it's set to 0
    -a int:float
       Indicate it's amplicon based calling.  Reads don't map to the amplicon will be skipped.  A read pair is considered belonging
       the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: 10:0.95
    -k 0/1
       Indicate whether to perform local realignment.  Default: 1 or yes.  Set to 0 to disable it.  For Ion or PacBio, 0 is recommended.
    -G Genome fasta
       The the reference fasta.  Should be indexed (.fai).  Default to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
    -R Region
       The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED is needed.
    -d delimiter
       The delimiter for split region_info, default to tab "\t"
    -n regular_expression
       The regular expression to extract sample name from bam filenames.  Default to: /([^\/\._]+?)_[^\/]*.bam/
    -N string
       The sample name to be used directly.  Will overwrite -n option
    -b string
       The indexed BAM file
    -c INT
       The column for chromosome
    -S INT
       The column for region start, e.g. gene start
    -E INT
       The column for region end, e.g. gene end
    -s INT
       The column for segment starts in the region, e.g. exon starts
    -e INT
       The column for segment ends in the region, e.g. exon ends
    -g INT
       The column for gene name, or segment annotation
    -x INT
       The number of nucleotide to extend for each segment, default: 0
    -f double
       The threshold for allele frequency, default: 0.05 or 5%
    -r minimum reads
       The minimum # of variance reads, default 2
    -B INT
       The minimum # of reads to determine strand bias, default 2
    -Q INT
       If set, reads with mapping quality less than INT will be filtered and ignored
    -q INT
       The phred score for a base to be considered a good call.  Default: 25 (for Illumina)
       For PGM, set it to ~15, as PGM tends to under estimate base quality.
    -m INT
       If set, reads with mismatches more than INT will be filtered and ignored.  Gaps are not counted as mismatches.  
       Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem is calculated as NM - Indels.  Default: 8,
       or reads with more than 8 mismatches will not be used.
    -T INT
       Trim bases after [INT] bases in the reads
    -X INT
       Extension of bp to look for mismatches after insersion or deletion.  Default to 3 bp, or only calls when they're within 3 bp.
    -P number
       The read position filter.  If the mean variants position is less that specified, it's considered false positive.  Default: 5
    -Z double
       For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with caution.  The
       downsampling will be random and non-reproducible.
    -o Qratio
       The Qratio of (good_quality_reads)/(bad_quality_reads+0.5).  The quality is defined by -q option.  Default: 1.5
    -O MapQ
       The reads should have at least mean MapQ to be considered a valid variant.  Default: no filtering
    -V freq
       The lowest frequency in normal sample allowed for a putative somatic mutations.  Default to 0.05
    -I INT
       The indel size.  Default: 120bp
    -M INT
       The minimum matches for a read to be considered.  If, after soft-clipping, the matched bp is less than INT, then the 
       read is discarded.  It's meant for PCR based targeted sequencing where there's no insert and the matching is only the primers.
       Default: 0, or no filtering 
    -J CRISPR_cutting_site
       The genomic position that CRISPR/Cas9 suppose to cut, typically 3bp from the PAM NGG site and within the guide.  For
       CRISPR mode only.  It will adjust the variants (mostly In-Del) start and end sites to as close to this location as possible,
       if there are alternatives.  The option should only be used for CRISPR mode.

AUTHOR
       Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
       Report bugs to zhongwu\@yahoo.com

COPYRIGHT
       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

USAGE
   exit(0);
}
