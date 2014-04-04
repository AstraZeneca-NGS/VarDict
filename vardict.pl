#!/usr/bin/env perl
# Parse a list of refseq and check CDS coverage
use warnings;
use Getopt::Std;
use strict;

our ($opt_h, $opt_H, $opt_b, $opt_D, $opt_M, $opt_d, $opt_s, $opt_c, $opt_S, $opt_E, $opt_n, $opt_N, $opt_e, $opt_g, $opt_x, $opt_f, $opt_r, $opt_B, $opt_z, $opt_v, $opt_p, $opt_F, $opt_C, $opt_m, $opt_Q, $opt_T, $opt_L, $opt_q, $opt_Z, $opt_X, $opt_P, $opt_3, $opt_k, $opt_R, $opt_G);
unless( getopts( 'hHvzpDCFLM3d:b:s:e:S:E:n:c:g:x:f:r:B:N:Q:m:T:q:Z:X:P:k:R:G:' )) {
    USAGE();
}
USAGE() if ( $opt_H );
USAGE() unless ( $opt_b );
my $BAM = $opt_b; # the bam file
my $sample = $1 if ( $BAM =~ /([^\/\._]+?).sorted[^\/]*.bam/ );
if ( $opt_n ) {
    $sample = $1 if ( $BAM =~ /$opt_n/ );
}
$sample = $opt_N if ( $opt_N );
unless($sample) { $sample = $1 if ( $BAM =~ /([^\/]+?)[_\.][^\/]*bam/ ); }

$opt_d = "\t" unless( $opt_d );
my $c_col = $opt_c ? $opt_c - 1 : 2;
my $S_col = $opt_S ? $opt_S - 1 : 6;
my $E_col = $opt_E ? $opt_E - 1 : 7;
my $s_col = $opt_s ? $opt_s - 1 : 9;
my $e_col = $opt_e ? $opt_e - 1 : 10;
my $g_col = $opt_g ? $opt_g - 1 : 12;
my $VEXT = $opt_X ? $opt_X : 5; # the extension of deletion and insertion for mismatches
$opt_P = $opt_P ? $opt_P : 5;
$opt_k = $opt_k ? $opt_k : 0; # The extension for indels for realignments

$s_col = $S_col if ( $opt_S && (!$opt_s) );
$e_col = $E_col if ( $opt_E && (!$opt_e) );

if ( $opt_L ) {
    $opt_p = 1;
    $c_col = 1 - 1;
    $S_col = 2 - 1;
    $E_col = 2 - 1;
    $s_col = 2 - 1;
    $e_col = 2 - 1;
    $g_col = 3 - 1;
    $opt_d = ":";
}
my $fasta = $opt_G ? $opt_G : "/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa";
my $EXT = defined($opt_x) ? $opt_x : 0;
my $FREQ = $opt_f ? $opt_f : 0.05;
my $BIAS = 0.05; # The cutoff to decide whether a position has read strand bias
my $MINB = $opt_B ? $opt_B : 2; # The minimum reads for bias calculation
my $MINR = $opt_r ? $opt_r : 2; # The minimum reads for variance allele
my $GOODQ = $opt_q ? $opt_q : 25; # The phred score in fastq to be considered as good base call
if ( $opt_p ) {
    $FREQ = -1;
    $MINR = 0;
}
if ( $opt_h ) {
    print join("\t", qw(Sample Gene Chr Start End Ref Alt Depth AltDepth RefFwdReads RefRevReads AltFwdReads AltRevReads Genotype AF Bias PMean PStd QMean QStd 5pFlankSeq 3pFlankSeq)), "\n";
}
my @SEGS = ();
if ( $opt_R ) {
    my ($chr, $reg, $gene) = split(/:/, $opt_R);
    $chr = "chr$chr" unless( $chr =~ /^chr/ );
    $gene = $chr unless( $gene );
    my ($start, $end) = split(/-/, $reg);
    $end = $start unless( $end );
    $start =~ s/,//g;
    $end =~ s/,//g;
    $start -= $EXT;
    $end += $EXT;
    push(@SEGS, [[$chr, $start, $end, $gene]]);
} else {
    while( <> ) {
        chomp;
        next if ( /^#/ );
        next if ( /^browser/i );
        next if ( /^track/i );
        my @A = split(/$opt_d/);
        my ($chr, $cdss, $cdse, $gene) = @A[$c_col, $S_col, $E_col, $g_col];
        my @starts = split(/,/, $A[$s_col]);
        my @ends = split(/,/, $A[$e_col]);
        my @CDS = ();
        $chr = "chr$chr" unless ($chr =~ /^chr/ );
        $gene = $chr unless( $gene );
        for(my $i = 0; $i < @starts; $i++) {
            my ($s, $e) = ($starts[$i], $ends[$i]);
            next if ( $cdss > $e ); # not a coding exon
            last if ( $cdse < $s ); # No more coding exon
            $s = $cdss if ( $s < $cdss );
            $e = $cdse if ( $e > $cdse );
            $s -= $EXT; # unless ( $s == $cdss );
            $e += $EXT; # unless ( $e == $cdse );
            $s++ if ( $opt_z );
            push(@CDS, [$chr, $s, $e, $gene]);
        }
        push(@SEGS, \@CDS);
    }
}

foreach my $segs (@SEGS) {
    my %REF;
    my %cov;
    my %var;
    my %hash;
    my ($chr, $START, $END, $gene);
    for(my $i = 0; $i < @$segs; $i++) {
	($chr, $START, $END, $gene) = @{$segs->[$i]};
	my $tsamcnt = `samtools view $BAM $chr:$START-$END | head`;
	next unless( $tsamcnt || $opt_p ); # to avoid too much IO for exome and targeted while the BED is whole genome
	my $s_start = $START - $EXT - 100 < 1 ? 1 : $START - $EXT - 100;
	my $s_end = $END + $EXT + 100;
	my ($header, $exon) = split(/\n/, `samtools faidx $fasta $chr:$s_start-$s_end`, 2);
	$exon =~ s/\s+//g;
	for(my $i = $s_start; $i <= $s_start + length($exon); $i++) {
	    $REF{ $i } = uc(substr( $exon, $i - ($s_start), 1 ));
	}
	$chr =~ s/^chr// if ( $opt_C );
	open(SAM, "samtools view $BAM $chr:$START-$END |");
	while( <SAM> ) {
	    if ( $opt_Z ) {
	        next if ( rand() <= $opt_Z );
	    }
	    my @a = split(/\t/);
	    next if ( defined($opt_Q) && $a[4] < $opt_Q );
	    if ( $opt_m ) {
	        if ( /XM:i:(\d+)/ ) {  # number of mismatches.  Don't use NM since it includes gaps, which can be from indels
		    next if ( $1 > $opt_m );
		} else {
		    print STDERR "No XM tag for mismatches. $_\n"; #if ( $opt_D );
		}
	    }
	    my $start = $a[3];
	    my $n = 0;
	    my $p = 0;
	    my $dir = $a[1] & 0x10 ? "-" : "+";
	    my @segs = $a[5] =~ /(\d+)[MI]/g; # Only match and insertion counts toward read length
	    my @segs2 = $a[5] =~ /(\d+)[MIS]/g; # For total length, including soft-clipped bases
	    my $rlen = 0; $rlen += $_ foreach(@segs); #$stat->sum(\@segs); # The read length for matched bases
	    my $rlen2= 0; $rlen2 += $_ foreach(@segs2); #$stat->sum(\@segs2); # The total length, including soft-clipped bases
	    my $offset = 0;
	    #while( $a[5] =~ /(\d+)([A-Z])/g ) {
	    my @cigar = $a[5] =~ /(\d+)([A-Z])/g;
	    for(my $ci = 0; $ci < @cigar; $ci += 2) {
		my $m = $cigar[$ci];
		my $C = $cigar[$ci+1];
		if ( $C eq "N" ) {
		    $start += $m;
		    $offset = 0;
		    next;
		} elsif ( $C eq "S" ) {
		    $n += $m;
		    $offset = 0;
		    next;
		} elsif ( $C eq "H" ) {
		    next;
		} elsif ( $C eq "I" ) {
		    my $s = substr($a[9], $n, $m);
		    my $q = substr($a[10], $n, $m);
		    my $ss = "";
                    my ($multoffs, $multoffp) = (0, 0); # For multiple indels within 10bp
                    if ( $cigar[$ci+2] && $cigar[$ci+2] <= $VEXT && $cigar[$ci+3] =~ /M/ && $cigar[$ci+5] && $cigar[$ci+5] =~ /[ID]/ ) {
                        $s .= "#" . substr($a[9], $n+$m, $cigar[$ci+2]);
                        $q .= substr($a[10], $n+$m, $cigar[$ci+2]); 
                        $s .= $cigar[$ci+5] eq "I" ? ("^" . substr($a[9], $n+$m+$cigar[$ci+2], $cigar[$ci+4])) : ("^" . $cigar[$ci+4]);
                        $q .= $cigar[$ci+5] eq "I" ? substr($a[10], $n+$m+$cigar[$ci+2], $cigar[$ci+4]) : substr($a[10], $n+$m+$cigar[$ci+2], 1);
                        $multoffs += $cigar[$ci+2] + ($cigar[$ci+5] eq "D" ? $cigar[$ci+4] : 0);
                        $multoffp += $cigar[$ci+2] + ($cigar[$ci+5] eq "I" ? $cigar[$ci+4] : 0);
                        $ci += 4;
		    } else {
			if ( $cigar[$ci+3] && $cigar[$ci+3] =~ /M/ ) {
			    for(my $vi = 0; $vi < $VEXT && $vi < $cigar[$ci+2]; $vi++) {
			        $offset = $vi+1 if ($REF{ $start+$vi } && substr($a[9], $n+$m+$vi, 1) ne $REF{ $start+$vi });
			    }
			    if ($offset) {
			        $ss .= substr($a[9], $n, $offset);
				$q .= substr($a[10], $n, $offset);
			    }
			}
		    }
		    $s = "$s&$ss" if ( $offset > 0);
		    if ( $start - 1 >= $START && $start -1 <= $END ) {
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
		    }
		    $n += $m+$offset+$multoffp;
		    $p += $m+$offset+$multoffp;
		    $start += $offset+$multoffs;
		    next;
		} elsif ( $C eq "D" ) {
		    my $s = "-$m";
		    my $ss = "";
                    my $q = substr($a[10], $n, 1);
                    my ($multoffs, $multoffp) = (0, 0); # For multiple indels within 10bp
                    if ( $cigar[$ci+2] && $cigar[$ci+2] <= $VEXT && $cigar[$ci+3] =~ /M/ && $cigar[$ci+5] && $cigar[$ci+5] =~ /[ID]/ ) {
                        $s .= "#" . substr($a[9], $n, $cigar[$ci+2]);
                        $q .= substr($a[10], $n, $cigar[$ci+2]);
                        $s .= $cigar[$ci+5] eq "I" ? ("^" . substr($a[9], $n+$cigar[$ci+2], $cigar[$ci+4])) : ("^" . $cigar[$ci+4]);
                        $q .= $cigar[$ci+5] eq "I" ? substr($a[10], $n+$cigar[$ci+2], $cigar[$ci+4]) : "";
                        $multoffs += $cigar[$ci+2] + ($cigar[$ci+5] eq "D" ? $cigar[$ci+4] : 0);
                        $multoffp += $cigar[$ci+2] + ($cigar[$ci+5] eq "I" ? $cigar[$ci+4] : 0);
                        $ci += 4;
                    } else {
			if ( $cigar[$ci+3] && $cigar[$ci+3] =~ /M/ ) {
			    for(my $vi = 0; $vi < $VEXT && $vi < $cigar[$ci+2]; $vi++) {
			        $offset = $vi+1 if ($REF{ $start+$m+$vi } && substr($a[9], $n+$vi, 1) ne $REF{ $start+$m+$vi });
			    }
			    if ($offset) {
			        $ss .= substr($a[9], $n, $offset);
				$q .= substr($a[10], $n, $offset);
			    }
			}
		    }
		    $s = "$s&$ss" if ( $offset > 0 );
		    $q .= substr($a[10], $n+$offset, 1);
		    if ( $start >= $START && $start <= $END ) {
			$hash{ $start }->{ $s }->{ $dir }++;
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
			$tmpq >= $GOODQ ? $hv->{ hicnt }++ : $hv->{ locnt }++;
			for(my $i = 0; $i < $m; $i++) {
			    $cov{ $start+$i }->{ $s }++;
			}
		    }
		    $start += $m+$offset+$multoffs;
		    $n += $offset+$multoffp;
		    $p += $offset+$multoffp;
		    next;
		}
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
                    my $q = ord(substr($a[10], $n, 1))-33;
                    # for more than one nucleotide mismatch
		    my $ss = "";
                    while(($start + 1) >= $START && ($start + 1) <= $END && ($i + 1) < $m && substr($a[9], $n, 1) ne $REF{$start} ) {
                        if ( substr($a[9], $n+1, 1) ne $REF{ $start + 1 } ) {
                            $ss .= substr($a[9], $n+1, 1);
                            $q += ord(substr($a[10], $n+1, 1))-33;
                            $n++;
                            $p++;
                            $i++;
                            $start++;
                        } else {
                            last;
                        }
                    }
		    $s .= "&$ss" if ( $ss );
		    unless( $trim ) {
			if ( $start - length($ss) >= $START && $start - length($ss) <= $END ) {
			    $hash{ $start - length($ss) }->{ $s }->{ $dir }++;
			    my $hv = $hash{ $start - length($ss) }->{ $s };
			    $hv->{ cnt }++;
			    my $tp = $p < $rlen-$p ? $p + 1: $rlen-$p;
			    $q = $ss? $q/(length($s)-1) : $q;
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
			    $q >= $GOODQ ? $hv->{ hicnt }++ : $hv->{ locnt }++;
			    $cov{ $start - length($ss) }->{ $s }++;
			}
		    }
		    $start++ unless( $C eq "I" );
		    $n++ unless( $C eq "D" );
		    $p++ unless( $C eq "D" );
		}
		$offset = 0;
	    }
	}
	close( SAM );
    }
    while( my ($p, $v) = each %hash ) {
	my @tmp = ();
	my $vcov = 0; #the variance coverage
	my @var = ();
	my @v = values %{ $cov{ $p } };
	my @vn = keys %$v;
	if ( @vn == 1 && $vn[0] eq $REF{ $p } ) {
	    next unless( $opt_p ); # ignore if only reference were seen and no pileup to avoid computation
	}
	my $tcov = 0; $tcov += $_ foreach(@v); #$stat->sum(\@v);
	my $hicov = 0;
	$hicov += $_->{ hicnt } ? $_->{ hicnt } : 0 foreach( values %$v );
	next if ( $tcov == 0 ); # ignore when there's no coverage
	# ignore when the ref allele is dominant.
	if ( $v->{ $REF{ $p } } ) {
	    my $refcnt = $v->{ $REF{ $p } }->{ cnt } ? $v->{ $REF{ $p } }->{ cnt } : 0;
	    next if ( (! $v->{ I } ) && $refcnt/$tcov > 1 - $FREQ );
	}
	while( my ($n, $cnt) = each %$v ) {
	    unless( $n eq "I") {
		my $fwd = $cnt->{ "+" } ? $cnt->{ "+" } : 0;
		my $rev = $cnt->{ "-" } ? $cnt->{ "-" } : 0;
		##my $bias = ($fwd/($fwd+$rev) >= $BIAS && $rev/($fwd+$rev) >= $BIAS && $fwd >= $MINB && $rev >= $MINB) ? 2 : 1;
		my $bias = 0;
		if ( $fwd + $rev <= 20 ) {
		    $bias = $fwd*$rev > 0 ? 2 : 0;
		} else {
		    $bias = ($fwd/($fwd+$rev) >= $BIAS && $rev/($fwd+$rev) >= $BIAS && $fwd >= $MINB && $rev >= $MINB) ? 2 : 1;
		}
		my $vqual = sprintf("%.1f", $cnt->{ qmean }/$cnt->{ cnt }); # base quality
		my $MQ = sprintf("%.1f", $cnt->{ Qmean }/$cnt->{ cnt }); # mapping quality
		my ($hicnt, $locnt) = ($cnt->{ hicnt } ? $cnt->{ hicnt } : 0, $cnt->{ locnt } ? $cnt->{ locnt } : 0);
		my $tvref = {n => $n, cov => $cnt->{ cnt }, fwd => $fwd, rev => $rev, bias => $bias, freq => ($fwd+$rev)/$tcov, pmean => sprintf("%.1f", $cnt->{ pmean }/$cnt->{ cnt } ), pstd => $cnt->{ pstd }, qual => $vqual, qstd => $cnt->{ qstd }, mapq => $MQ, qratio => sprintf("%.3f", $hicnt/($locnt+0.5)), hifreq => ($hicov > 0 ? $hicnt/$hicov : 0) };
		push(@var, $tvref);
		if ( $opt_M ) {
		    push( @tmp, "$n:" . ($fwd + $rev) . ":F-$fwd:R-$rev:$tvref->{freq}:$tvref->{bias}:" . join(",", @{$cnt->{p}}) . ":$tvref->{pstd}:" . join(",", @{$cnt->{q}}) . ":$tvref->{qstd}");
		} elsif ( $opt_D ) {
		    push( @tmp, "$n:" . ($fwd + $rev) . ":F-$fwd:R-$rev:" . sprintf("%.5f", $tvref->{freq}) . ":$tvref->{bias}:$tvref->{pmean}:$tvref->{pstd}:$vqual:$tvref->{qstd}:" . sprintf("%.5f", $tvref->{hifreq}) . ":$tvref->{mapq}:$tvref->{qratio}");
		}
	    }
	}
	if ( $v->{ I } ) {
	    while( my ($n, $cnt) = each %{ $v->{ I } } ) {
		my $fwd = $cnt->{ "+" } ? $cnt->{ "+" } : 0;
		my $rev = $cnt->{ "-" } ? $cnt->{ "-" } : 0;
		#my $bias = ($fwd/($fwd+$rev) >= $BIAS && $rev/($fwd+$rev) >= $BIAS && $fwd >= $MINB && $rev >= $MINB) ? 2 : 1;
		my $bias = 0;
		if ( $fwd + $rev <= 20 ) {
		    $bias = $fwd*$rev > 0 ? 2 : 0;
		} else {
		    $bias = ($fwd/($fwd+$rev) >= $BIAS && $rev/($fwd+$rev) >= $BIAS && $fwd >= $MINB && $rev >= $MINB) ? 2 : 1;
		}
		my $vqual = sprintf("%.1f", $cnt->{ qmean }/$cnt->{ cnt }); # base quality
		my $MQ = sprintf("%.1f", $cnt->{ Qmean }/$cnt->{ cnt }); # mapping quality
		my ($hicnt, $locnt) = ($cnt->{ hicnt } ? $cnt->{ hicnt } : 0, $cnt->{ locnt } ? $cnt->{ locnt } : 0);
		my $tvref = {n => $n, cov => $cnt->{ cnt }, fwd => $fwd, rev => $rev, bias => $bias, freq => ($fwd+$rev)/$tcov, pmean => sprintf("%.1f", $cnt->{ pmean }/$cnt->{ cnt } ), pstd => $cnt->{ pstd }, qual => $vqual, qstd => $cnt->{ qstd }, mapq => $MQ, qratio => sprintf("%.3f", $hicnt/($locnt+0.5)), hifreq => ($hicov > 0 ? $hicnt/$hicov : 0) };
		push(@var, $tvref);
		if( $opt_M ) {
		    push( @tmp, "I$n:" . ($fwd + $rev) . ":F-$fwd:R-$rev:$tvref->{freq}:$tvref->{bias}:" . join(",", @{$cnt->{p}}) . ":$tvref->{pstd}:" . join(",", @{$cnt->{q}}) . ":$tvref->{qstd}");
		} elsif ( $opt_D ) {
		    push( @tmp, "I$n:" . ($fwd + $rev) . ":F-$fwd:R-$rev:" . sprintf("%.5f", $tvref->{freq}) . ":$tvref->{bias}:$tvref->{pmean}:$tvref->{pstd}:$vqual:$tvref->{qstd}:" . sprintf("%.5f", $tvref->{hifreq}) . ":$tvref->{mapq}:$tvref->{qratio}" );
		}
	    }
	}
	@var = sort { $b->{ qual } * $b->{cov} <=> $a->{ qual } * $a->{cov}; } @var;
	my ($freq, $bias, $hifreq);

	# Make sure the first bias is always for the reference nucleotide
	my ($pmean, $pstd, $qual, $qstd, $mapq, $qratio);
	my ($rfc, $rrc) = (0, 0); # coverage for referece forward and reverse strands
	my ($vfc, $vrc) = (0, 0); # coverage for variance forward and reverse strands
	my $genotype1 = "$var[0]->{n}";
	my $genotype2 = "";
	my $vn;
	if ( $var[0]->{ n } eq $REF{ $p } ) {
	    unless( $var[1] ) {
	        next unless ($opt_p ); # ignore no or low frequency variance unless pileup is needed
		# When pileup is needed 
		$freq = 0;
		$hifreq = 0;
		$vcov = 0;
		$vn = $REF{ $p }; # Same as reference
		($vfc, $vrc) = (0, 0);
		$genotype2 = $vn;
		$bias = $var[0]->{ bias } . ";0";
		($pmean, $pstd, $qual, $qstd, $mapq, $qratio) = ($var[0]->{ pmean }, $var[0]->{ pstd }, $var[0]->{ qual }, $var[0]->{ qstd }, $var[0]->{ mapq }, $var[0]->{ qratio });
		push(@tmp, join("|", "Ref", $pmean, $pstd, $qual, $qstd, sprintf("%.5f", $hifreq), $mapq, $qratio));
	    } else { # when there're non-reference reads
		$freq = $var[1]->{ freq };
		$hifreq = $var[1]->{ hifreq };
		$vcov = $var[1]->{ cov };
		$vn = $var[1]->{ n };
		$genotype2 = $vn;
		($vfc, $vrc) = ($var[1]->{ fwd }, $var[1]->{ rev });
		$bias = $var[0]->{ bias } . ";" . $var[1]->{ bias };
		($pmean, $pstd, $qual, $qstd, $mapq, $qratio) = ($var[1]->{ pmean }, $var[1]->{ pstd }, $var[1]->{ qual }, $var[1]->{ qstd }, $var[1]->{ mapq }, $var[1]->{ qratio });
		push(@tmp, join("|", "Ref", $var[0]->{ pmean }, $var[0]->{ pstd }, $var[0]->{ qual }, $var[0]->{ qstd }, sprintf("%.5f", $var[0]->{ hifreq }), $var[0]->{ mapq }, $var[0]->{ qratio }));
		push(@tmp, join("|", "Alt", $pmean, $pstd, $qual, $qstd, sprintf("%.5f", $hifreq), $mapq, $qratio));
	    }
	    ($rfc, $rrc) = ($var[0]->{ fwd }, $var[0]->{ rev });
	} else {
	    $freq = $var[0]->{ freq };
	    $hifreq = $var[0]->{ hifreq };
	    $vcov = $var[0]->{ cov };
	    $vn = $var[0]->{ n };
	    ($vfc, $vrc) = ($var[0]->{ fwd }, $var[0]->{ rev });
	    ($pmean, $pstd, $qual, $qstd, $mapq, $qratio) = ($var[0]->{ pmean }, $var[0]->{ pstd }, $var[0]->{ qual }, $var[0]->{ qstd }, $var[0]->{ mapq }, $var[0]->{ qratio });
	    push(@tmp, join("|", "Alt", $pmean, $pstd, $qual, $qstd, sprintf("%.5f", $hifreq), $mapq, $qratio));
	    if ( $var[1] && $var[1]->{ n } eq $REF{ $p } ) {
		$bias = $var[1]->{ bias } . ";" . $var[0]->{ bias };
		$genotype2 = $var[1]->{ n };
		($rfc, $rrc) = ($var[1]->{ fwd }, $var[1]->{ rev });
		push(@tmp, join("|", "Ref", $var[1]->{ pmean }, $var[1]->{ pstd }, $var[1]->{ qual }, $var[1]->{ qstd }, sprintf("%.5f", $var[1]->{ hifreq }), $var[1]->{ mapq }, $var[1]->{ qratio }));
	    } else {
	        $bias = "0;" . $var[0]->{ bias };
		$genotype2 = $var[1]->{ n } && $var[1]->{ freq } >= $FREQ ? $var[1]->{ n } : $var[0]->{ n };
	    }
	}
	if ( $opt_F ) {
	    $freq = 0;
	    $hifreq = 0;
	    foreach my $v (@var) {
	        $freq += $v->{ freq } unless( $v->{n} eq $REF{ $p } );
	        $hifreq += $v->{ hifreq } unless( $v->{n} eq $REF{ $p } );
	    }
	}
	my $dellen = $vn =~ /^-(\d+)/ ? $1 : 0;
	my $ep = $vn =~ /^\+/ ?  $p : ($vn =~ /^-/ ? $p + $dellen - 1: $p);
	next unless ( $vcov >= $MINR );
	#next if ($bias eq "2;1"); # Ignore strand biased variances
	#next if ($pstd == 0); # Ignore if variances are all in the same position of reads
	my ($refallele, $varallele) = ("", "");
	my ($shift3, $msi) = (0, 0);  # how many bp can a deletion be shifted to 3'
	my $sp = $p;
	my $extrafreq = 0;
	if ( $vn =~ /^\+/ ) {
	    if ( $opt_k ) {
		my $extra = ($vn =~ /&([ATGC]+)/) ? $1 : "";
		my $wupseq = join( "", (map { $REF{ $_ }; } (($p - $opt_k - length($vn)+1) .. $p))) . $vn; # 5' flanking seq
		$wupseq =~ s/\+//;
		my $sanpseq = $vn . $extra . join( "", (map { $REF{ $_ }; } (($p + length($extra)+1) .. ($p + length($vn) + $opt_k)))); # 3' flanking seq
		$sanpseq =~ s/^\+//;
		for( my $tp5 = $p - $opt_k - length($vn); $tp5 <= $p + length($vn)+$opt_k; $tp5++) {
		    next unless( $hash{ $tp5 } );
		    while( my ($tn, $tv) = each %{ $hash{ $tp5 } }) {
		        next if ($tn eq $REF{ $tp5 });
			next if ( $tn eq "I" );
			next if ( $tn =~ /^-/ );
			next if ( $tv->{ qmean }/$tv->{ cnt } < $GOODQ );
			next if ( $tv->{ pmean }/$tv->{ cnt } > $dellen/2 + $opt_k);
			my $head = $tn;
			$head =~ s/&//;
			$head =~ s/^\+//;
			my $flag = 0;
			$flag = 1 if ( $tp5 == $p + 1 && $sanpseq =~ /^$head/ );
			$flag = 1 if ( $tp5 + length($head) == $p + 1 + length($extra) && $wupseq =~ /$head$/);
			#print STDERR join("\t", $tn, $head, $p, $tp5, $extra, $flag, $wupseq, $sanpseq), "\n";
			if ( $flag ) {
			    $freq += $tv->{ cnt }/$tcov;
			    $extrafreq += $tv->{ cnt }/$tcov;
			    $hifreq += $hicov > 0 ? $tv->{ hicnt }/$hicov : 0;
			    $vfc += $tv->{ "+" } ? $tv->{ "+" } : 0;
			    $vrc += $tv->{ "-" } ? $tv->{ "-" } : 0;
			    $vcov += ($tv->{ "+" } ? $tv->{ "+" } : 0) + ($tv->{ "-" } ? $tv->{ "-" } : 0);
			}
		    }
		}
	    }
	    unless( $vn =~ /&/ || $vn =~ /#/ ) {
		my $tseq1 = $vn;
		$tseq1 =~ s/^\+//;
		my $tseq2 = join("", (map { $REF{ $_ }; } (($p+1) .. ($p+length($tseq1)))));
		while( $tseq1 eq $tseq2 ) {
		    $shift3 += length($tseq1);
		    $p += length($tseq1);
		    $tseq2 = join("", (map { $REF{ $_ }; } (($p+1) .. ($p+length($tseq1)))));
		}
		$msi = $shift3/length($tseq1);
	    }
	    if ( $opt_3 ) {
		$sp += $shift3;
		$ep += $shift3;
	    } else {
	        $p -= $shift3;
	    }
	    ($refallele, $varallele) = ($REF{$p}, $vn);
	    $varallele =~ s/^\+//;
	    $varallele = $REF{$p} . $varallele;
	} elsif ( $vn =~ /^-/ ) {
	    $varallele = $vn;
	    $varallele =~ s/^-\d+//;
	    if ( $opt_k ) {
		my $extra = ($vn =~ /&([ATGC]+)/) ? $1 : "";
		my $wupseq = join( "", (map { $REF{ $_ }; } (($p - $dellen - $opt_k) .. ($p-1)))) . $extra; # 5' flanking seq
		my $sanpseq = $extra . join( "", (map { $REF{ $_ }; } (($p + $dellen + length($extra)) .. ($p + 2*$dellen + $opt_k)))); # 3' flanking seq
		for( my $tp5 = $p - $dellen - $opt_k; $tp5 <= $p + $dellen; $tp5++) {
		    next unless( $hash{ $tp5 } );
		    while( my ($tn, $tv) = each %{ $hash{ $tp5 } }) {
		        next if ($tn eq $REF{ $tp5 });
			next if ( $tn eq "I" );
			next if ( $tn =~ /^-/ );
			next if ( $tv->{ qmean }/$tv->{ cnt } < $GOODQ );
			next if ( $tv->{ pmean }/$tv->{ cnt } > $dellen/2 + $opt_k);
			my $head = $tn;
			$head =~ s/&//;
			my $tseq = "";
			if ( $tp5 > $p ) {
			    my $inseq1 = join( "", (map { $REF{ $_ }; } ($p .. ($tp5-1))) );
			    my $inseq2 = join( "", (map { $REF{ $_ }; } (($p+$dellen) .. ($dellen+$tp5-1))) );
			    $tseq = $inseq1 if ( $inseq1 eq $inseq2 );
			}
			my $flag = 0;
			$flag = 1 if ( $tp5 - length($tseq) == $p && $sanpseq =~ /^$tseq$head/ );
			$flag = 1 if ( $tp5 + length($head) == $p + $dellen + length($extra) && $wupseq =~ /$head$/);
			#print STDERR join("\t", $tn, $head, $p, $tp5, $tseq, $extra, $flag, $wupseq, $sanpseq), "\n";
			if ( $flag ) {
			    $freq += $tv->{ cnt }/$tcov;
			    $extrafreq += $tv->{ cnt }/$tcov;
			    $hifreq += $hicov > 0 ? $tv->{ hicnt }/$hicov : 0;
			    $vfc += $tv->{ "+" } ? $tv->{ "+" } : 0;
			    $vrc += $tv->{ "-" } ? $tv->{ "-" } : 0;
			    $vcov += ($tv->{ "+" } ? $tv->{ "+" } : 0) + ($tv->{ "-" } ? $tv->{ "-" } : 0);
			}
		    }
		}
	    }
	    unless( $vn =~ /&/ || $vn =~ /#/ ) {
		while($REF{ $sp } eq $REF{ $ep+1 } ) {
		    $sp++;
		    $ep++;
		    $p++;
		    $shift3++;
		}
		$msi = $shift3/$dellen; # if ( $shift3%$dellen == 0 );
		unless ( $opt_3 ) {
		    $sp -= $shift3;
		    $ep -= $shift3;
		    $p -= $shift3;
		}
		$varallele = $REF{ $p - 1 };
		$refallele = $REF{ $p - 1 };
		$sp--;
	    }
	    for(my $i = 0; $i < $dellen; $i++) {
	        $refallele .= $REF{ $p + $i };
	    }
	} else {
	    ($refallele, $varallele) = ($REF{ $p }, $vn);
	}
	unless( $opt_p ) {
	    next if ( $varallele =~ /N/ );
	    next if ($refallele =~ /N/ );
	    next if ($pmean < $opt_P);
	    next if ($qual < $GOODQ );
	}
	if ( $vn =~ /&(.*)$/ ) {
	    my $extra = $1;
	    $varallele =~ s/&//;
	    for(my $m = 0; $m < length($extra); $m++) {
	        $refallele .= $REF{ $ep + $m + 1 };
		$genotype1 .= $REF{ $ep + $m + 1 } unless( $genotype1 =~ /^-/ || $genotype1 =~ /^\+/ || length($genotype1) > 1 );
	    }
	    $ep += length($extra);
	}
	if ($vn =~ /#(.+)\^(.+)/) {
	    my $mseq = $1;
	    my $tail = $2;
	    $ep += length($mseq);
	    $refallele .= $mseq;
	    if ( $tail =~ /^(\d+)/ ) {
		for(my $ti = 0; $ti < $1; $ti++) {
		    $refallele .= $REF{ $ep + $ti + 1 };
		}
	        $ep += $1;
	    }
	    $varallele =~ s/#//;
	    $varallele =~ s/\^(\d+)?//;
	    $genotype1 =~ s/#//;
	    $genotype1 =~ s/\^/-/;
	    $genotype2 =~ s/#//;
	    $genotype2 =~ s/\^/-/;
	}
	next if ( $freq < $FREQ );
	my $left5 = join( "", (map { $REF{ $_ }; } (($sp - 10) .. ($sp - 1))) ); # left 10 nt
	my $right5 = join( "", (map { $REF{ $_ }; } (($ep + 1) .. ($ep + 10))) ); # right 10 nt
	my $genotype = "$genotype1/$genotype2";
	$genotype =~ s/&//g;
	$genotype =~ s/#//g;
	$genotype =~ s/\^/-/g;
	$freq = sprintf("%.5f", $freq);
	$hifreq = sprintf("%.5f", $hifreq);
	$extrafreq = sprintf("%.5f", $extrafreq) if ($extrafreq);
	if ( $opt_v ) {
	    print join("\t", $chr, $sp, ".", $refallele, $varallele, ".", "PASS", "DP=$tcov;END=$ep;AF=$freq;BIAS=$bias;SM=$sample;PMEAN=$pmean;PSTD=$pstd;QUAL=$qual;QSTD=$qstd"), "\n";
	} else {
	    print join("\t", $sample, $gene, $chr, $sp, $ep, $refallele, $varallele, $tcov, $vcov, $rfc, $rrc, $vfc, $vrc, $genotype, $freq, $bias, $pmean, $pstd, $qual, $qstd, $mapq, $qratio, $hifreq, $extrafreq, $shift3, $msi, $left5, $right5);
	    print "\t" . join(" & ", @tmp) if ( $opt_D || $opt_M );
	    print "\n";
        }
    }
    if ( $opt_p ) {
	for(my $i = 0; $i < @$segs; $i++) {
	    my ($tchr, $START, $END, $tgene) = @{$segs->[$i]};
	    for(my $j = $START; $j <= $END; $j++) {
		#print STDERR "Position $chr:$j has no coverage!\n" unless( $hash{ $j } );
		print join("\t", $sample, $tgene, $tchr, $j, $j, $REF{ $j }, "", 0, 0, 0, 0, 0, 0, "", 0, "0;0", 0, 0, 0, 0, 0, "", 0, 0, 0, 0, "", ""), "\n" unless( $hash{ $j } );
	    }
	}
    }
}
#FCB02N4ACXX:3:2206:20108:2526#GATGGTTC  163     chr3    38181981        50      79M188N11M      =       38182275        667     TGAAGTTGTGTGTGTCTGACCGCGATGTCCTGCCTGGCACCTGTGTCTGGTCTATTGCTAGTGAGCTCATCGTAAAGAGGTGCCGCCGGG \YY`c`\ZQPJ`e`b]e_Sbabc[^Ybfaega_^cafhR[U^ee[ec][R\Z\__ZZbZ\_\`Z`d^`Zb]bBBBBBBBBBBBBBBBBBB AS:i:-8 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:72A16A0    YT:Z:UU XS:A:+  NH:i:1  RG:Z:15
sub USAGE {
    print STDERR <<USAGE;
    $0 [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] [-x #_nu] [-g gene] [-f freq] [-r #_reads] [-B #_reads] region_info

    The program will calculate candidate variance for a given region(s) in an indexed BAM file. The default
    input is IGV's one or more entries in refGene.txt, but can be any regions in 1-based end-inclusive coordinates.

    -H Print this help page
    -h Print a header row decribing columns
    -z Indicate whether zero-based coordinates, as IGV does (and BED). Affects a given BED file, not option R below.
    -v VCF format output
    -p Do pileup regardless of frequency
    -C Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
    -D Debug mode.  Will print some error messages and append full genotype at the end.
    -M Similar to -D, but will append individual quality and position data instead of mean
    -3 Indicate to move deletions to 3-prime if alternative alignment can be achieved.
    -k Indel extension
       Indicate the number of bp to rescue forcely aligned reads in deletions and insertions to better represent frequency.  Use with caution.
    -G Genome fasta
       The the reference fasta. Should be indexed (.fai). Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
    -R Region
       The region of interest. In the format of chr:start-end, in 1-based end-inclusive coordinates. 
       If end is omitted, then a single position.  No BED is needed.
    -d delimiter
       The delimiter for split region_info, default to tab "\t"
    -n regular_expression
       The regular expression to extract sample name from bam filenames.  Default to: /([^\/\._]+?)_[^\/]*.bam/
    -N string
       The sample name to be used directly.  Will overwrite -n
    -b string
       The indexed BAM file
    -c INT
       The column for chr
    -S INT
       The column for region start, e.g. gene start
    -E INT
       The column for region end, e.g. gene end
    -s INT
       The column for segment starts in the region, e.g. exon starts
    -e INT
       The column for segment ends in the region, e.g. exon ends
    -g INT
       The column for gene name
    -x INT
       The number of nucleotide to extend for each segment, default: 0
    -f double
       The threshold for allele frequency, default: 0.05 or 5%
    -F Indicate to calculate the frequency as the sum of all non-reference variants, 
       instead of just the most frequent allele, which is the default
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
       Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem currently doesn't have such SAM tag
    -T INT
       Trim bases after [INT] bases in the reads
    -X INT
       Extension of bp to look for mismatches after insertion or deletion.  Default to 5 bp.
    -P number
       The read position filter.  If the mean variants position is less that specified, it's considered false positive.  Default: 5
    -Z double
       For downsampling fraction. .g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with caution.  The
       downsampling will be random and non-reproducible.
    -L Used for command line pipe, such as "echo chr:pos:gene | checkVar.pl -L".  Will automatically set "-d : -p -c 1 -S 2 -E 2 -g 3"
USAGE
   exit(0);
}
