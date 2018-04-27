#!/usr/bin/env perl
# Parse a list of refseq and check CDS coverage

use warnings;
#use Getopt::Std;
use Getopt::Long qw{ :config bundling no_ignore_case no_auto_abbrev};
use strict;

our ($opt_h, $opt_H, $opt_b, $opt_D, $opt_d, $opt_s, $opt_c, $opt_S, $opt_E, $opt_n, $opt_N, $opt_e, $opt_g, $opt_x, $opt_f, $opt_r, $opt_B, $opt_z, $opt_p, $opt_F, $opt_m, $opt_Q, $opt_T, $opt_q, $opt_Z, $opt_X, $opt_P, $opt_3, $opt_k, $opt_R, $opt_G, $opt_a, $opt_o, $opt_O, $opt_V, $opt_t, $opt_y, $opt_I, $opt_i, $opt_M, $opt_L, $opt_U, $opt_w, $opt_W, $opt_A, $opt_J, $opt_j, $opt_u);
#unless( getopts( 'hHtzypD3iUF:d:b:s:e:S:E:n:c:g:x:f:r:B:N:Q:m:T:q:Z:X:P:k:R:G:a:o:O:V:I:M:L:w:W:A:J:j:' )) {
#    USAGE();
#}
my @adaptor;
my $CHIMERIC;
GetOptions( "h|header" => \$opt_h,
	    "H|?" => \$opt_H,
	    "dedup|t" => \$opt_t,
	    "z" => \$opt_z,
	    "verbose|y" => \$opt_y,
	    "chimeric" => \$CHIMERIC,
	    "p" => \$opt_p,
	    "u" => \$opt_u,
	    "debug|D" => \$opt_D,
	    "3" => \$opt_3,
	    "splice|i" => \$opt_i,
	    "nosv|U" => \$opt_U,
	    "F=s" => \$opt_F,
	    "d=s" => \$opt_d,
	    "b=s" => \$opt_b,
	    "s=i" => \$opt_s,
	    "e=i" => \$opt_e,
	    "S=i" => \$opt_S,
	    "E=i" => \$opt_E,
	    "c=i" => \$opt_c,
	    "g=i" => \$opt_g,
	    "n=s" => \$opt_n,
	    "N=s" => \$opt_N,
	    "x=i" => \$opt_x,
	    "X=i" => \$opt_X,
	    "f=f" => \$opt_f,
	    "r=i" => \$opt_r,
	    "B=i" => \$opt_B,
	    "Q=f" => \$opt_Q,
	    "m=i" => \$opt_m,
	    "T|trim=i" => \$opt_T,
	    "q=f" => \$opt_q,
	    "Z|downsample=f" => \$opt_Z,
	    "P=i" => \$opt_P,
	    "k=i" => \$opt_k,
	    "R=s" => \$opt_R,
	    "a|amplicon=s" => \$opt_a,
	    "o=f" => \$opt_o,
	    "O=f" => \$opt_O,
	    "G=s" => \$opt_G,
	    "V=f" => \$opt_V,
	    "I=i" => \$opt_I,
	    "M=i" => \$opt_M,
	    "L=i" => \$opt_L,
	    "w|insert-size=i" => \$opt_w,
	    "W|insert-std=i" => \$opt_W,
	    "A=f" => \$opt_A,
	    "adaptor=s" => \@adaptor,
	    "J|crispr=i" => \$opt_J,
	    "j=i" => \$opt_j ) || Usage();
USAGE() if ( $opt_H );
USAGE() unless ( $opt_b );
my %adaptor;
my %adaptor_rev;
my $ADSEED = 6; # adaptor seed
if ( @adaptor ) {
    foreach my $seq (@adaptor) {
	for(my $i = 0; $i <= 6 && $i + $ADSEED < length($seq); $i++) {
	    my $seed = substr($seq, $i, $ADSEED);
	    my $rseed = reverse($seed);
	    $rseed =~ y/ATGC/TACG/;
	    $adaptor{ $seed } = $i + 1;
	    $adaptor_rev{ $rseed } = $i + 1;
	}
    }
}
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
$opt_m = defined($opt_m) ? $opt_m : 8;

$s_col = $S_col if ( $opt_S && (!$opt_s) );
$e_col = $E_col if ( $opt_E && (!$opt_e) );

$opt_F = defined($opt_F) ? $opt_F : "0x500";

my $VEXT = defined($opt_X) ? $opt_X : 3; # the extension of deletion and insertion for complex variants
$opt_P = defined($opt_P) ? $opt_P : 5;
$opt_k = defined($opt_k) ? $opt_k : 1; # Whether to perform local realignment.  Set to 0 to disable.
my %GENOMES = ( hg19 => "/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa",
		hg38 => "/ngs/reference_data/genomes/Hsapiens/hg38/seq/hg38.fa",
		mm10 => "/ngs/reference_data/genomes/Mmusculus/mm10/seq/mm10.fa" );
my $fasta = $opt_G ? ($GENOMES{ $opt_G } ? $GENOMES{ $opt_G } : $opt_G) : "/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa";
my %CHRS; # Key: chr Value: chr_len
my $EXT = defined($opt_x) ? $opt_x : 0;
my $FREQ = $opt_f ? $opt_f : 0.01;
my $QRATIO = $opt_o ? $opt_o : 1.5; # The Qratio
my $BIAS = 0.05; # The cutoff to decide whether a positin has read strand bias
my $MINB = $opt_B ? $opt_B : 2; # The minimum reads for bias calculation
my $MINR = $opt_r ? $opt_r : 2; # The minimum reads for variance allele
my $GOODQ = defined($opt_q) ? $opt_q : 22.5; # The phred score in fastq to be considered as good base call
$opt_O = defined($opt_O) ? $opt_O : 0; # The minimun mean mapping quality to be considered
$opt_V = defined($opt_V) ? $opt_V : 0.05; # The minimun alelle frequency allowed in normal for a somatic mutation
my $SVMINLEN = $opt_L ? $opt_L : 1000; # The minimum structural variant length to be presented using <DEL> <DUP> <INV> <INS>, etc.
my $RLEN = 0; # The read length
my $LOWQUAL = 10; # any base with quality <= 10 will be consider low quality in soft-clipped seq and extention will stop
my $SEED1 = 17; # The larger seed size
my $SEED2 = 12; # The smaller seed size
my $DISCPAIRQUAL = 35; # The minimum mapping quality when structural variant is only supported by discordant pairs
my $INDELSIZE = $opt_I ? $opt_I : 50;
my $MINMATCH = defined($opt_M) ? $opt_M : 25;
my %SPLICE;
my $INSSIZE = $opt_w ? $opt_w : 300; # Mean Insert size
my $INSSTD = $opt_W ? $opt_W : 100; # Insert std
my $INSSTDAMT = $opt_A ? $opt_A : 4; # Insert std amount
my $SVMAXLEN = 150000; # Max Structure variant size to be called in realignment step
my %REVCOMP = ( A => "T", T => "A", C => "G", G => "C" );
my $SVFLANK = 50; # the flanking sequence length for SV
my $MINSVCDIST = 1.5; # the minimum distance between two SV clusters in term of read length
my %SOFTP2SV;
if ( $opt_p ) {
    $FREQ = -1;
    $MINR = 0;
}
if ( $opt_h ) {
    print join("\t", qw(Sample Gene Chr Start End Ref Alt Depth AltDepth RefFwdReads RefRevReads AltFwdReads AltRevReads Genotype AF Bias PMean PStd QMean QStd MQ Sig_Noise HiAF ExtraAF shift3 MSI MSI_NT NM HiCnt HiCov 5pFlankSeq 3pFlankSeq Seg VarType Duprate SV_info));
    print $opt_J ? "\tCRISPR\n" : "\n";
}

my ($BAM1, $BAM2) = split(/\|/, $BAM);
my @BAMS = $BAM1 =~ /^http/ ? ($BAM1) : split(/:/, $BAM1);
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
	    s/\r//g;
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
	    my @A = split(/$opt_d/, $_, -1);
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
    } elsif ( $var =~ /<(...)>/ ) {
	return $1;
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
    while( my ($p, $pv) = each %$vars ) {
	unless( $pv->{ SV } ) {
	    next unless( $p >= $S && $p <= $E );
	}
	my @vts = ();
	my @vrefs = ();
	unless( $pv->{ VAR } ) {
	    next unless( $opt_p );
	    my $vref = $pv->{ REF };
	    unless($vref) {
		print join("\t", $sample, $G, $chr, $p, $p, "", "", 0, 0, 0, 0, 0, 0, "", 0, "0;0", 0, 0, 0, 0, 0, "", 0, 0, 0, 0, 0, 0, "", "", 0, 0, "$chr:$S-$E", "", 0, 0), "\n";
		next;
	    }
	    push(@vts, "");
	    push(@vrefs, $vref);
	} else {
	    for(my $i = 0; $i < @{ $pv->{VAR} }; $i++) {
		my $vref = $pv->{ VAR }->[$i];
		last if ( $vref->{refallele} =~ /N/ );
		my $vartype = varType($vref->{refallele}, $vref->{varallele});
		unless( isGoodVar( $pv->{ VAR }->[$i], $pv->{ REF }, $vartype ) ) {
		    next unless( $opt_p );
		}
		push(@vts, $vartype);
		push(@vrefs, $vref);
	    }
	}
	for(my $vi = 0; $vi < @vts; $vi++) {
	    my ($vartype, $vref) = ($vts[$vi], $vrefs[$vi]);
	    adjComplex($vref) if ( $vartype eq "Complex" );
	    my $duprate = $vref->{ duprate } ? $vref->{ duprate } : 0;
	    my $crispr = $opt_J ? ($vref->{ crispr } ? $vref->{ crispr } : 0) : "";
	    print join("\t", $sample, $G, $chr, (map { $vref->{ $_ } ? $vref->{ $_ } : 0; } @hds), "$chr:$S-$E", $vartype, $duprate, $pv->{ SV } ? $pv->{ SV } : 0);
	    print "\t$crispr" if ( $opt_J );
	    print "\t", $vref->{ DEBUG } if ( $opt_D );
	    print "\n";
	}
    }
    
}

# Adjust the complex var
sub adjComplex {
    my $vref = shift;
    my ($refnt, $varnt) = ($vref->{ refallele }, $vref->{ varallele });
    next if ( $varnt =~ /^</ );
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
    #print STDERR "$chr:$S-$E\n";
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
	     my $vref2 = $vars2->{$p}->{ VAR }->[0];
	     adjComplex($vars2->{$p}->{ VAR }->[0]) if ( $vartype eq "Complex" );
	     print join("\t", $sample, $G, $chr, (map{ $vref2->{ $_ }; } @hd1), (map{ 0; } @hdrs), (map { $vref2->{ $_ }; } (@hdrs, @hd2)), "$chr:$S-$E", "Deletion", $vartype, "0", "0", $vref2->{ duprate } ? $vref2->{ duprate } : 0, $vars2->{$p}->{ SV } ? $vars2->{$p}->{ SV } : 0 ), "\n";
	 } elsif( ! $vars2->{ $p } ) { # no coverage for sample 2
	     $vartype = varType($vars1->{$p}->{ VAR }->[0]->{ refallele }, $vars1->{$p}->{ VAR }->[0]->{ varallele }) if ($vars1->{$p}->{ VAR });
	     next unless( $vars1->{$p}->{ VAR } && isGoodVar($vars1->{$p}->{ VAR }->[0], $vars1->{$p}->{REF}, $vartype) );
	     my $vref1 = $vars1->{$p}->{ VAR }->[0];
	     adjComplex($vars1->{$p}->{ VAR }->[0]) if ( $vartype eq "Complex" );
	     print join("\t", $sample, $G, $chr, (map { $vref1->{ $_ }; } (@hd1, @hdrs)), (map { 0; } @hdrs), (map { $vref1->{$_};} @hd2), "$chr:$S-$E", "SampleSpecific", $vartype, $vref1->{ duprate } ? $vref1->{ duprate } : 0, $vars1->{$p}->{ SV } ? $vars1->{$p}->{ SV } : 0, 0, 0 ), "\n";
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
			 $type = "StrongSomatic" if ( isNoise($v2->{ VARN }->{ $nt }) && $vartype eq "SNV" );
			 #if ($type =~ /Somatic/) {
			 #    my $newtype = combineAnalysis($VREF, $v2->{ VARN }->{ $nt }, $chr, $p, $nt);
			 #    if ( $newtype eq "FALSE" ) {$N++; next;}
			 #    $type = $newtype if ( $newtype );
			 #}
			 #$type = "StrongSomatic" if ( isNoise($v2->{ VARN }->{ $nt }) );
			 print join("\t", $sample, $G, $chr, (map { $VREF->{ $_ }; } (@hd1, @hdrs)), (map { $v2->{ VARN }->{ $nt }->{ $_ }; } (@hdrs, @hd2)), "$chr:$S-$E", $type, $vartype, $VREF->{ duprate } ? $VREF->{ duprate } : 0, $v1->{ SV } ? $v1->{ SV } : 0, $v2->{ VARN }->{ $nt }->{ duprate } ? $v2->{ VARN }->{ $nt }->{ duprate } : 0, $v2->{ SV } ? $v2->{ SV } : 0), "\n";
		     } else { # sample 1 only, should be strong somatic
			 my @tvf = $v2->{ REF } ? (map { $v2->{ REF }->{ $_ } ? $v2->{ REF }->{ $_ } : 0; } @hdrs) : ($v2->{ VAR }->[0]->{ tcov } ? $v2->{ VAR }->[0]->{ tcov } : 0, map { 0; } (1..17));
			 my $type = "StrongSomatic";
			 if ($vartype ne "SNV" && (length($nt) > 10 || $nt =~ /-\d\d/)) {
			     $v2->{ VARN }->{ $nt }->{ cov } = 0;  # Ensure it's initialized before passing to combineAnalysis
			     my $newtype = $VREF->{ cov } < $MINR + 3 && $nt !~ /</ ? combineAnalysis($VREF, $v2->{ VARN }->{ $nt }, $chr, $p, $nt) : "";
			     if ( $newtype eq "FALSE" ) {$N++; next;}
			     $type = $newtype if ( $newtype );
			 }
			 if ( $type ne "StrongSomatic" ) {
			     print join("\t", $sample, $G, $chr, (map { $VREF->{ $_ }; } (@hd1, @hdrs)), (map { $v2->{ VARN }->{ $nt }->{ $_ }; } @hdrs),(map { $VREF->{$_}; } @hd2), "$chr:$S-$E", $type, $vartype, $VREF->{ duprate } ? $VREF->{ duprate } : 0, $v1->{ SV } ? $v1->{ SV } : 0, $v2->{ VARN }->{ $nt }->{ duprate } ? $v2->{ VARN }->{ $nt }->{ duprate } : 0, $v2->{ SV } ? $v2->{ SV } : 0), "\n";
			 } else {
			     print join("\t", $sample, $G, $chr, (map { $VREF->{ $_ }; } (@hd1, @hdrs)), @tvf, (map { $VREF->{$_};} @hd2), "$chr:$S-$E", "StrongSomatic", $vartype, $VREF->{ duprate } ? $VREF->{ duprate } : 0, $v1->{ SV } ? $v1->{ SV } : 0, $v2->{ VARN }->{ $nt }->{ duprate } ? $v2->{ VARN }->{ $nt }->{ duprate } : 0, $v2->{ SV } ? $v2->{ SV } : 0), "\n";
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
			 print join("\t", $sample, $G, $chr, (map { $v1->{ VARN }->{ $nt }->{ $_ }; } (@hd1, @hdrs)), (map { $v2->{ VAR }->[0]->{ $_ }; } (@hdrs, @hd2)), "$chr:$S-$E", $type, varType($v1->{ VARN }->{ $nt }->{ refallele }, $v1->{ VARN }->{ $nt }->{ varallele }),  $v1->{ VARN }->{ $nt }->{ duprate } ? $v1->{ VARN }->{ $nt }->{ duprate } : 0, $v1->{ SV } ? $v1->{ SV } : 0, $v2->{ VAR }->[0]->{ duprate } ? $v2->{ VAR }->[0]->{ duprate } : 0, $v2->{ SV } ? $v2->{ SV } : 0), "\n";
		     } else {
			 my @th1 = $v1->{ REF } ? (map { $v1->{ REF }->{ $_ } } @hdrs) : ($v1->{ VAR }->[0]->{ tcov }, (map { 0; } @hdrs[1..$#hdrs]));
			 adjComplex($v2->{ VAR }->[0]) if ( $vartype eq "Complex" );
			 my $vref2 = $v2->{ VAR }->[0];
			 print join("\t", $sample, $G, $chr, (map { $v2->{ VAR }->[0]->{ $_ }; } @hd1), @th1, (map { $v2->{ VAR }->[0]->{ $_ }; } (@hdrs, @hd2)), "$chr:$S-$E", "StrongLOH", $vartype, 0, 0, $vref2->{ duprate } ? $vref2->{ duprate } : 0, $v2->{ SV } ? $v2->{ SV } : 0), "\n";
		     }
		 }
	     } elsif ( $v2->{ VAR } ) { # sample 1 has only reference
		 $vartype = varType($v2->{ VAR }->[0]->{ refallele }, $v2->{ VAR }->[0]->{ varallele });
		 next unless( isGoodVar( $v2->{ VAR }->[0], $v2->{ REF }, $vartype ) );
		 # potential LOH
		 my $nt = $v2->{ VAR }->[0]->{ n };
		 my $type = "StrongLOH";
		 $v1->{ VARN }->{ $nt }->{ cov } = 0;
		 my $newtype = $v2->{ VARN }->{ $nt }->{ cov } < $MINR + 3 && $nt !~ /</ && (length($nt) > 10 || $nt =~ /-\d\d/ ) ? combineAnalysis($v2->{ VARN }->{ $nt }, $v1->{ VARN }->{ $nt }, $chr, $p, $nt) : "";
		 next if ( $newtype eq "FALSE" );
		 $type = $newtype if ( $newtype );
		 my @th1 = $newtype ? (map { $v1->{ VARN }->{ $nt }->{ $_ } ? $v1->{ VARN }->{ $nt }->{ $_ } : 0; } @hdrs) : (map { $v1->{ REF }->{ $_ } ? $v1->{ REF }->{ $_ } : 0 } @hdrs);
		 adjComplex($v2->{ VAR }->[0]) if ( $vartype eq "Complex" );
		 my $vref2 = $v2->{ VAR }->[0];
		 print join("\t", $sample, $G, $chr, (map { $v2->{ VAR }->[0]->{ $_ } ? $v2->{ VAR }->[0]->{ $_ } : 0; } @hd1), @th1, (map { $v2->{ VAR }->[0]->{ $_ } ? $v2->{ VAR }->[0]->{ $_ } : 0; } (@hdrs, @hd2)), "$chr:$S-$E", $type, $vartype, 0, 0, $vref2->{ duprate } ? $vref2->{ duprate } : 0, $v2->{ SV } ? $v2->{ SV } : 0), "\n";
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
    return "" if ( $var1->{ ep } - $var1->{ sp } > $SVMINLEN ); # Don't do it for structural variants
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
    return 1 if ( $vref->{ freq } > 0.30 );
    return 0 if ( $vref->{ mapq } < $opt_O );
    return 0 if ( $vref->{ msi } >= 15 && $vref->{ freq } <= 0.25 && $vref->{ msint } == 1);
    return 0 if ( $vref->{ msi } >= 12 && $vref->{ freq } <= 0.1 && $vref->{ msint } > 1 );
    if ( $vref->{ bias } eq "2;1" && $vref->{ freq } < 0.20 ) {
	return 0 unless($type && $type ne "SNV" && (length($vref->{refallele}) >= 3 || length($vref->{varallele}) >= 3));
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
    my ($chr, $START, $END, $ref, $ext) = @_;
    my %REF = ();
    $ref = $ref ? $ref : \%REF;
    my $extension = $ext ? $ext : 15000;
    my $s_start = $START - $EXT - $extension < 1 ? 1 : $START - $EXT - $extension;
    my $s_end = $END + $EXT + $extension > $CHRS{ $chr } ? $CHRS{ $chr } : $END + $EXT + $extension;
    print STDERR "TIME: Getting REF: ", time(), "\n" if ( $opt_y );
    my ($header, $exon) = split(/\n/, `samtools faidx $fasta $chr:$s_start-$s_end`, 2);
    return $ref if ( isLoaded($chr, $s_start, $s_end, $ref) );
    push(@{ $ref->{ loaded } }, [$chr, $s_start, $s_end]);
    $exon =~ s/\s+//g;
    for(my $i = 0; $i <= length($exon); $i++) {
	next if ( $ref->{ $i + $s_start } ); # don't process it more than once
	$ref->{ $i + $s_start } = uc(substr( $exon, $i, 1 ));
	push( @{ $ref->{ uc(substr( $exon, $i, $SEED1)) } }, $i + $s_start ) if ( length($exon) - $i > $SEED1 );
	push( @{ $ref->{ uc(substr( $exon, $i, $SEED2)) } }, $i + $s_start ) if ( length($exon) - $i > $SEED2 );
    }
    print STDERR "TIME: Got REF.", time(), "\n" if ( $opt_y );
    return $ref;
}

# Check whether a region is alreday loaded in reference
sub isLoaded {
    my ($chr, $s, $e, $ref) = @_;
    return 0 unless( $ref->{ loaded } );
    foreach my $r (@{ $ref->{ loaded } }) {
	return 1 if ( $chr eq $r->[0] && $s >= $r->[1] && $e <= $r->[2] );
    }
    return 0;
}

# Add deletion structural variant
sub addSV {
    my ($sdref, $s, $e, $ms, $me, $dir, $rlen, $mlen, $softp, $rp, $q, $Q, $nm) = @_;
    $sdref->{ cnt }++;
    $sdref->{ dir } = $dir;
    $sdref->{ $dir }++;
    #$sdref->{ pmean } += $rp;
    #$sdref->{ qmean } += $q;
    #$sdref->{ Qmean } += $Q;
    #$sdref->{ nm } += $nm;
    $q >= $GOODQ ? $sdref->{ hicnt }++ : $sdref->{ locnt }++;
    $sdref->{ start } = $s unless( $sdref->{ start } && $sdref->{ start } < $s );
    $sdref->{ end } = $e unless( $sdref->{ end } && $sdref->{ end } > $e );
    push(@{ $sdref->{ mates } }, [$ms, $me, $mlen, $s, $e, $rp, $q, $Q, $nm]);
    $sdref->{ mstart } = $ms unless( $sdref->{ mstart } && $sdref->{ mstart } < $ms );
    $sdref->{ mend } = $ms+$rlen unless( $sdref->{ mend } && $sdref->{ mend } > $me );
    if ($softp) {
	if ( $dir == 1 ) {
	    $sdref->{ soft }->{ $softp }++ if ( abs($softp - $sdref->{end}) < 10 );
	} else {
	    $sdref->{ soft }->{ $softp }++ if ( abs($softp - $sdref->{start}) < 10 );
	}
    }
}

# Construct a variant structure given a segment and BAM files.
sub parseSAM {
    my ($chr, $START, $END, $bams, $REF, $hash, $cov, $sclip5, $sclip3, $svflag) = @_;
    my %COV = ();
    $cov = $cov ? $cov : \%COV;
    my %HASH = ();
    $hash = $hash ? $hash : \%HASH;
    my %dels5;
    my %ins;
    my %mnp; # Keep track of MNPs
    my %SCLIP3 = (); # soft clipped at 3'
    my %SCLIP5 = (); # soft clipped at 5'
    $sclip3 = $sclip3 ? $sclip3 : \%SCLIP3;
    $sclip5 = $sclip5 ? $sclip5 : \%SCLIP5;
    %SPLICE = ();
    my @svfdel;
    my @svrdel;
    my ($svdelfend, $svdelrend) = (0, 0); # for strutral variant: Deletion
    my @svfdup;
    my @svrdup;
    my ($svdupfend, $svduprend) = (0, 0); # for strutral variant: Dupplication
    my @svfinv3;
    my @svfinv5;
    my @svrinv3;
    my @svrinv5;
    my ($svinvfend3, $svinvrend3) = (0, 0); # for strutral variant: Inversion
    my ($svinvfend5, $svinvrend5) = (0, 0); # for strutral variant: Inversion
    my %svffus;
    my %svrfus;
    my %svfusfend;
    my %svfusrend;
    my @svfins;
    my @svrins;
    my ($svinsfend, $svinsrend) = (0, 0); # for strutral variant: Insertion
    #$REF = getREF($chr, $START, $END);
    my ($totalreads, $dupreads) = (0, 0);
    foreach my $bami (@$bams) {
	#my $tsamcnt = `samtools view $bami $chr:$START-$END | head -1`;
	#next unless( $tsamcnt || $opt_p ); # to avoid too much IO for exome and targeted while the BED is whole genome
	# Get the reference sequence
	my $SAMFILTER = $opt_F ? "-F $opt_F" : "";
	print STDERR "TIME: Start parsing SAM: ", time(), " $bami $chr:$START-$END\n" if ( $opt_y );
	open(SAM, "samtools view $SAMFILTER $bami $chr:$START-$END |");
	my %DUP = ();
	my $DUPP = 0;
	while( <SAM> ) {
	    if ( $opt_Z ) {
		next if ( rand() <= $opt_Z );
	    }
	    my @a = split(/\t/, $_, 12);
	    next if ( defined($opt_Q) && $a[4] < $opt_Q ); # ignore low mapping quality reads
	    #next if ( $opt_F && ($a[1] & 0x200 || $a[1] & 0x100) ); # ignore "non-primary alignment", or quality failed reads.
	    #next if ( $opt_F && ($a[1] & 0x100) ); # ignore "non-primary alignment"
	    if ( $a[1] & 0x100 ) {
		next if ($opt_F); # -F 0 to turn it off
	    }
	    next if ( $a[9] eq "*" );
	    # filter duplicated reads if option -t is set
	    $totalreads++;
	    if ( $opt_t ) {
		%DUP = () if ( $a[3] ne $DUPP );
		if ( $a[7] =~ /^\d/o ) {
		    if ( $DUP{ "$a[3]-$a[6]-$a[7]" } ) {
			$dupreads++;
			next;
		    }
		    $DUP{ "$a[3]-$a[6]-$a[7]" } = 1;
		    $DUPP = $a[3];
		} elsif ( $a[1] & 0x8 ) {  # mate not mapped
		    if ( $DUP{ "$a[3]-$a[5]" } ) {
			$dupreads++;
			next;
		    }
		    $DUP{ "$a[3]-$a[5]" } = 1;
		    $DUPP = $a[3];
		}
	    }
	    my $nm = 0;
	    my @segid = $a[5] =~ /(\d+)[ID]/g; # For total indels
	    my $idlen = 0; $idlen += $_ foreach(@segid);
	    if ( $a[11] && $a[11] =~ /NM:i:(\d+)/io ) {  # number of mismatches.  Don't use NM since it includes gaps, which can be from indels
		$nm = $1 - $idlen;
		next if ( $opt_m && $nm > $opt_m ); # edit distance - indels is the # of mismatches
	    } else {
		print STDERR "No NM tag for mismatches. $a[0] $a[2] $a[3] $a[11]\n" if ( $opt_y && $a[5] ne "*" );
		#next;
	    }
	    my $n = 0; # keep track the read position, including softclipped
	    my $p = 0; # keep track the position in the alignment, excluding softclipped
	    my $dir = $a[1] & 0x10 ? -1 : 1;

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
	    # Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's algorithm to figure out indels
	    my $offset = 0;
	    if( $a[5] =~ s/^(\d+)D// ) { $a[3] += $1; }
	    $a[5] =~ s/\d+D$//;
	    $a[5] =~ s/^(\d+)I/$1S/;
	    $a[5] =~ s/(\d+)I$/$1S/;
	    while( $idlen > 0 && $opt_k ) {
		my $flag = 0;
		if ($a[5] =~ /^(\d+)S(\d+)([ID])/o) {
		    my $tslen = $1 + ($3 eq "I" ? $2 : 0);
		    $tslen .= "S";
		    $a[3] += $3 eq "D" ? $2 : 0;
		    $a[5] =~ s/^\d+S\d+[ID]/$tslen/;
		    $flag = 1;
		}
		if ($a[5] =~ /(\d+)([ID])(\d+)S$/o) {
		    my $tslen = $3 + ($2 eq "I" ? $1 : 0);
		    $tslen .= "S";
		    $a[5] =~ s/\d+[ID]\d+S$/$tslen/;
		    $flag = 1;
		}
		if ($a[5] =~ /^(\d+)S(\d+)M(\d+)([ID])/o) {
		    if ( $2 <= 10 ) {
			my $tslen = $1 + $2 + ($4 eq "I" ? $3 : 0);
			$tslen .= "S";
			$a[3] += $2 + ($4 eq "D" ? $3 : 0);
			$a[5] =~ s/^\d+S\d+M\d+[ID]/$tslen/;
			$flag = 1;
		    }
		}
		if ($a[5] =~ /(\d+)([ID])(\d+)M(\d+)S$/o) {
		    if ( $3 <= 10 ) {
			my $tslen = $4 + $3 + ($2 eq "I" ? $1 : 0);
			$tslen .= "S";
			$a[5] =~ s/\d+[ID]\d+M\d+S$/$tslen/;
			$flag = 1;
		    }
		}

		# The following two clauses to make indels at the end of reads as softly
		# clipped reads and let VarDict's algorithm identify indels
		if ($a[5] =~ /^(\d)M(\d+)([ID])(\d+)M/o) {
		    my $tslen = $1 + ($3 eq "I" ? $2 : 0);
		    my $mlen = $4;
		    $a[3] += $1 + ($3 eq "D" ? $2 : 0);
		    my $tn = 0;
		    $tn++ while( $tn < $mlen && $REF->{ $a[3]+$tn } && $REF->{ $a[3]+$tn } ne substr($a[9], $tslen+$tn, 1));
		    $tslen += $tn;
		    $mlen -= $tn;
		    $a[3] += $tn;
		    $a[5] =~ s/^\dM\d+[ID]\d+M/${tslen}S${mlen}M/;
		    $flag = 1;
		}
		if ($a[5] =~ /(\d+)([ID])(\d)M$/o) {
		    my $tslen = $3 + ($2 eq "I" ? $1 : 0);
		    $tslen .= "S";
		    $a[5] =~ s/\d+[ID]\d+M$/$tslen/;
		    $flag = 1;
		}

		# Combine two deletions and insertion into one complex if they are close
		if ($a[5] =~ /^(.*?)(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M/o) {
		    my $tslen = $4 + $5 + $6;
		    my $dlen = $3 + $4 + $6 + $7;
		    my $mid = $4 + $6;
		    my $ov5 = $1;
		    my $refoff = $a[3] + $2;
		    my $rdoff = $2;
		    my $RDOFF = $2;
		    my $rm = $8;
		    if ( $ov5 ) {
			my @rdp = $ov5 =~ /(\d+)[MIS]/g;  # read position
			my @rfp = $ov5 =~ /(\d+)[MND]/g;  # reference position
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
			$tslen = $dlen . "D" . $rm . "M";
		    } else {
			$tslen = "${dlen}D${tslen}I${rm}M";
		    }
		    if ( $mid <= 15 ) {
			#print STDERR "B: $rn $RDOFF $dlen $tslen $a[5]\n";
			#$a[5] =~ s/\d+M\d+D\d+M\d+I\d+M\d+D/${RDOFF}M${dlen}D${tslen}I/;
			$a[5] =~ s/\d+M\d+D\d+M\d+I\d+M\d+D\d+M/${RDOFF}M$tslen/;
			$flag = 1;
		    }
		} elsif ($a[5] =~ /^(.*?)(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M/o) { # three deletions
		    my $tslen = $4 + $6;
		    my $dlen = $3 + $4 + $5 + $6 + $7;
		    my $mid = $4 + $6;
		    my $ov5 = $1;
		    my $refoff = $a[3] + $2;
		    my $rdoff = $2;
		    my $RDOFF = $2;
		    my $rm = $8;
		    if ( $ov5 ) {
			my @rdp = $ov5 =~ /(\d+)[MIS]/g;  # read position
			my @rfp = $ov5 =~ /(\d+)[MND]/g;  # reference position
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
			$tslen = $dlen . "D" . $rm . "M";
		    } else {
			$tslen = "${dlen}D${tslen}I${rm}M";
		    }
		    if ( $mid <= 15 ) {
			#print STDERR "B: $rn $RDOFF $dlen $tslen $a[5]\n";
			$a[5] =~ s/\d+M\d+D\d+M\d+D\d+M\d+D\d+M/${RDOFF}M$tslen/;
			$flag = 1;
		    }
		} elsif ($a[5] =~ /^(.*?)(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M/o) { # three indels
		    my $tslen = $5 + $8;
		    $tslen += $3 if ( $4 eq "I" );
		    $tslen += $6 if ( $7 eq "I" );
		    $tslen += $9 if ( $10 eq "I" );
		    my $dlen = $5 + $8;
		    $dlen += $3 if ( $4 eq "D" );
		    $dlen += $6 if ( $7 eq "D" );
		    $dlen += $9 if ( $10 eq "D" );
		    my $mid = $5 + $8;
		    my $ov5 = $1;
		    my $refoff = $a[3] + $2;
		    my $rdoff = $2;
		    my $RDOFF = $2;
		    my $rm = $11;
		    if ( $ov5 ) {
			my @rdp = $ov5 =~ /(\d+)[MIS]/g;  # read position
			my @rfp = $ov5 =~ /(\d+)[MND]/g;  # reference position
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
			$tslen = $dlen . "D" . $rm . "M";
		    } else {
			$tslen = "${dlen}D${tslen}I${rm}M";
		    }
		    if ( $mid <= 15 ) {
			#print STDERR "B: $rn $RDOFF $dlen $tslen $a[5]\n";
			$a[5] =~ s/\d+M\d+[DI]\d+M\d+[DI]\d+M\d+[DI]\d+M/${RDOFF}M$tslen/;
			#print STDERR "A: $rn $RDOFF $dlen $tslen $a[5]\n";
			$flag = 1;
		    }
		}

		# Combine two close deletions (<10bp) into one
		if ($a[5] =~ /(\d+)D(\d+)M(\d+)([DI])(\d+I)?/o ) {
		    if ( $2 <= 15 ) {
			my $dlen = $1 + $2 + ($4 eq "D" ? $3 : 0);
			my $ilen = $2 + ($4 eq "I" ? $3 : 0);
			my $istr = $5;
			if ( $istr && $4 eq "D") {
			    $ilen += $1 if ( $istr =~ /(\d+)I/ );
			}
			$a[5] =~ s/\d+D\d+M\d+[DI](\d+I)?/${dlen}D${ilen}I/;
			$flag = 1;
		    }
		}

		# Combine two close indels (<10bp) into one
		if ($a[5] =~ /(\D)(\d+)I(\d+)M(\d+)([DI])(\d+I)?/o && $1 ne "D") {
		    if ( $3 <= 15 ) {
			my $dlen = $3 + ($5 eq "D" ? $4 : 0);
			my $ilen = $2 + $3 + ($5 eq "I" ? $4 : 0);
			my $istr = $6;
			if ( $istr && $5 eq "D" ) {
			    $ilen += $1 if ( $istr =~ /(\d+)I/ );
			}
			$a[5] =~ s/\d+I\d+M\d+[DI](\d+I)?/${dlen}D${ilen}I/;
			$flag = 1;
		    }
		}
		if ( $a[5] =~ /(\d+)D(\d+)D/ ) {
		    my $dlen = $1 + $2;
		    $a[5] =~ s/\d+D\d+D/${dlen}D/;
		    $flag = 1;
		}
		if ( $a[5] =~ /(\d+)I(\d+)I/ ) {
		    my $ilen = $1 + $2;
		    $a[5] =~ s/\d+I\d+I/${ilen}I/;
		    $flag = 1;
		}
		last unless ( $flag );
	    }
	    # The following two clauses to capture sometimes mis-softly clipped reads by aligner
	    if ($a[5] =~ /^(.*?)(\d+)M(\d+)S$/o) {
		my $ov5 = $1;
		my $mch = $2;
		my $soft = $3;
		my $refoff = $a[3] + $2;
		my $rdoff = $2;
		if ( $ov5 ) {
		    my @rdp = $ov5 =~ /(\d+)[MIS]/g;  # read position
		    my @rfp = $ov5 =~ /(\d+)[MND]/g;  # reference position
		    $refoff += $_ foreach(@rfp); 
		    $rdoff += $_ foreach(@rdp); 
		}
		my $rn = 0;
		my %RN = ();
		$rn++ while( $rn < $soft && $REF->{ $refoff + $rn } && $REF->{ $refoff + $rn } eq substr($a[9], $rdoff + $rn, 1) );
		if ( $rn > 0 ) {
		    $mch += $rn;
		    $soft -= $rn;
		    if ( $soft > 0 ) {
			$a[5] =~ s/\d+M\d+S$/${mch}M${soft}S/;
		    } else {
			$a[5] =~ s/\d+M\d+S$/${mch}M/;
		    }
		    $rn = 0;
		}
		if ( $soft > 0 ) {
		    while( $rn + 1 < $soft && $REF->{ $refoff + $rn + 1} && $REF->{ $refoff + $rn + 1} eq substr($a[9], $rdoff + $rn + 1, 1) && ord(substr($a[10], $rdoff + $rn + 1, 1))-33 > $LOWQUAL ) {
			$rn++;
			$RN{ $REF->{ $refoff + $rn + 1} }++;
		    }
		    #my $rn2 = $rn + 1;
		    #$rn2++ while( $rn2 + 1 < $soft && $REF->{ $refoff + $rn2 + 1} && $REF->{ $refoff + $rn2 + 1} eq substr($a[9], $rdoff + $rn2 + 1, 1) && ord(substr($a[10], $rdoff + $rn2 + 1, 1))-33 > $LOWQUAL );
		    #if ( $rn2 - $rn > 3 ) {
			#print STDERR "RN: $rn $rn2 $a[9] $a[0] $a[5] $dir $a[3]\n";
			#$rn = $rn2; 
		    #}
		    my @rn_nt = keys %RN; # Don't adjust for homopolymers
		    if ( ($rn > 4 && @rn_nt > 1) || ($REF->{ $refoff } && $REF->{ $refoff } eq substr($a[9], $rdoff, 1) )) {
			$mch += $rn + 1;
			$soft -= $rn + 1;
			if ( $soft > 0 ) {
			    $a[5] =~ s/\d+M\d+S$/${mch}M${soft}S/;
			} else {
			    $a[5] =~ s/\d+M\d+S$/${mch}M/;
			}
			#$rn++;
		    }
		    if ( $rn == 0 ) {
			my ($rrn, $rmch) = (0, 0);
			while( $rrn < $mch && $rn < $mch) {
			    last unless( $REF->{ $refoff - $rrn - 1 } );
			    if ( $REF->{ $refoff - $rrn - 1 } ne substr($a[9], $rdoff - $rrn - 1, 1) ) {
				$rn = $rrn+1;
				$rmch = 0;
			    } elsif ( $REF->{ $refoff - $rrn - 1 } eq substr($a[9], $rdoff - $rrn - 1, 1) ) {
				$rmch++;
			    }
			    $rrn++;
			    last if ( $rmch >= 3 ); # Stop at three consecute matches
			}

			#$rn++ while( $REF->{ $refoff - $rn - 1 } && $REF->{ $refoff - $rn - 1 } ne substr($a[9], $rdoff - $rn - 1, 1) );
			if ( $rn > 0 && $rn < $mch ) {
			    $soft += $rn;
			    $mch -= $rn;
			    $a[5] =~ s/\d+M\d+S$/${mch}M${soft}S/;
			}
		    }
		}
	    } elsif ( $a[5] =~ /^(.*?)(\d+)M$/o ) {  # Make >=3 mismatches in the end as soft clipping
		my $ov5 = $1;
		my $mch = $2;
		my $refoff = $a[3] + $2;
		my $rdoff = $2;
		if ( $ov5 ) {
		    my @rdp = $ov5 =~ /(\d+)[MIS]/g;  # read position
		    my @rfp = $ov5 =~ /(\d+)[MND]/g;  # reference position
		    $refoff += $_ foreach(@rfp); 
		    $rdoff += $_ foreach(@rdp); 
		}
		my $rn = 0;
		my ($rrn, $rmch) = (0, 0);
		while( $rrn < $mch && $rn < $mch ) {
		    last unless( $REF->{ $refoff - $rrn - 1 } );
		    if ( $REF->{ $refoff - $rrn - 1 } ne substr($a[9], $rdoff - $rrn - 1, 1) ) {
			$rn = $rrn+1;
			$rmch = 0;
		    } elsif ( $REF->{ $refoff - $rrn - 1 } eq substr($a[9], $rdoff - $rrn - 1, 1) ) {
			$rmch++;
		    }
		    $rrn++;
		    last if ( $rmch >= 3 );
		}
		$mch -= $rn;
		$a[5] =~ s/\d+M$/${mch}M${rn}S/ if ( $rn >= 3 );
	    }
	    if ( $a[5] =~ /^(\d+)S(\d+)M/o ) {
		my $mch = $2;
		my $soft = $1;
		my $rn = 0;
		my %RN = ();
		$rn++ while($rn < $soft && $REF->{ $a[3] - $rn - 1 } && $REF->{ $a[3] - $rn - 1 } eq substr($a[9], $soft - $rn - 1, 1));
		if ( $rn > 0 ) {
		    $mch += $rn;
		    $soft -= $rn;
		    if ( $soft > 0 ) {
			$a[5] =~ s/^\d+S\d+M/${soft}S${mch}M/;
		    } else {
			$a[5] =~ s/^\d+S\d+M/${mch}M/;
		    }
		    $a[3] -= $rn;
		    $rn = 0;
		}
		if ( $soft > 0 ) {
		    while( $rn + 1 < $soft && $REF->{ $a[3] - $rn - 2} && $REF->{ $a[3] - $rn - 2} eq substr($a[9], $soft - $rn - 2, 1) && ord(substr($a[10], $soft - $rn - 2, 1))-33 > $LOWQUAL ) {
			$rn++;
			$RN{ $REF->{ $a[3] - $rn - 2} }++;
		    }
		    my @rn_nt = keys %RN; # Don't adjust for homopolymers
		    if ( ($rn > 4 && @rn_nt > 1)  || ($REF->{ $a[3] - 1 } && $REF->{ $a[3] - 1 } eq substr($a[9], $soft - 1, 1))) {
			$mch += $rn + 1;
			$soft -= $rn + 1;
			if ( $soft > 0 ) {
			    $a[5] =~ s/^\d+S\d+M/${soft}S${mch}M/;
			} else {
			    $a[5] =~ s/^\d+S\d+M/${mch}M/;
			}
			$a[3] -= $rn + 1;
		    }
		    if ( $rn == 0 ) {
			my ($rrn, $rmch) = (0, 0);
			while( $rrn < $mch && $rn < $mch ) {
			    last unless( $REF->{ $a[3] + $rrn } );
			    if ( $REF->{ $a[3] + $rrn } ne substr($a[9], $soft + $rrn, 1) ) {
				$rn = $rrn+1;
				$rmch = 0;
			    } elsif ( $REF->{ $a[3] + $rrn } eq substr($a[9], $soft + $rrn, 1) ) {
				$rmch++;
			    }
			    $rrn++;
			    last if ( $rmch >= 3 ); # Stop at three consecute matches
			}
			#$rn++ while( $REF->{ $a[3] + $rn } && $REF->{ $a[3] + $rn } ne substr($a[9], $soft + $rn, 1) );
			if( $rn > 0 && $rn < $mch ) {
			    $soft += $rn;
			    $mch -= $rn;
			    $a[5] =~ s/^\d+S\d+M/${soft}S${mch}M/;
			    $a[3] += $rn;
			}
		    }
		}
	    } elsif ( $a[5] =~ /^(\d+)M/o ) { # Make >=3 mismatches in the end as soft clipping
		my $mch = $1;
		my $rn = 0;
		my ($rrn, $rmch) = (0, 0);
		while( $rrn < $mch && $rn < $mch ) {
		    last unless( $REF->{ $a[3] + $rrn } );
		    if( $REF->{ $a[3] + $rrn } ne substr($a[9], $rrn, 1) ) {
			$rn = $rrn+1;
			$rmch = 0;
		    } elsif ( $REF->{ $a[3] + $rrn } eq substr($a[9], $rrn, 1) ) {
			$rmch++;
		    }
		    $rrn++;
		    last if ( $rmch >= 3 );
		}
		if ( $rn >= 3 ) {
		    $mch -= $rn;
		    $a[5] =~ s/^\d+M/${rn}S${mch}M/;
		    $a[3] += $rn;
		}
	    }

	    next if ( $a[5] =~ /^\d\dS.*\d\dS$/ ); # Ignore reads that are softclipped at both ends and both greater than 10bp
	    my $start = $a[3];
	    my @cigar = $a[5] =~ /(\d+)([A-Z])/g;
	    my @segs = $a[5] =~ /(\d+)[MI]/g; # Only match and insertion counts toward read length
	    my @segs2 = $a[5] =~ /(\d+)[MIS]/g; # For total length, including soft-clipped bases
	    my $rlen = 0; $rlen += $_ foreach(@segs); #$stat->sum(\@segs); # The read length for matched bases
	    next if ( $opt_M && $rlen < $MINMATCH );
	    my $rlen2= 0; $rlen2 += $_ foreach(@segs2); #$stat->sum(\@segs2); # The total length, including soft-clipped bases
	    $RLEN = $rlen2 if ($rlen2 > $RLEN); # Determine the read length

	    next if ( $opt_F && $a[1] & 0x800 ); # Ignore the supplementary alignment so that it won't skew the coverage
	    
	    # Determine whether to filter a read in CRISPR mode
	    if ( $opt_J ) {
		my $rlen3= 0; $rlen3 += $1 while( $a[5] =~ /(\d+)[MD=X]/g ); # The total aligned length, excluding soft-clipped bases and insertions
		if ( $opt_j ) {
		    next unless ( $opt_J - $start > $opt_j && $start + $rlen3 - $opt_J > $opt_j );
		}
	    }

	    if ( $a[1] & 0x8 ) { # Mate unmapped, potential insertion
		# to be implemented
	    } elsif($a[4] > 10 && (!$opt_U)) { # Consider high mapping quality mates only
		my $mdir = $a[1] & 0x20 ? -1 : 1;
		my $mstart = $a[7];
		my $mend = $a[7] + $rlen2;
		my @msegs = $a[5] =~ /(\d+)[MND]/g;
		my $end = $a[3]; $end += $_ foreach( @msegs );
		my $soft5 = 0;
		if ( $a[5] =~ /^(\d+)S|^\d+H/ ) {
		    $soft5 = $start if ( $1 && ord(substr($a[10], $1 - 1, 1)) - 33 > $GOODQ );
		}
		my $soft3 = 0;
		if ( $a[5] =~ /(\d+)S$|H$/ ) {
		   $soft3 = $end if ( $1 && ord(substr($a[10], -$1, 1)) - 33 > $GOODQ );
		}
		my $MIN_D = 75;
		if ( $a[6] eq "=" ) {
		    my $mlen = $a[8];
		    if ( $dir * $mdir == -1 && $mlen * $dir > 0 ) { # deletion candidate
			$mlen = $mstart > $start ? $mend - $start : $end - $mstart;
			if( abs($mlen) > $INSSIZE + $INSSTDAMT * $INSSTD ) {
			    if ( $dir == 1 ) {
				push(@svfdel, { cnt => 0 } ) if ( @svfdel == 0 || $start - $svdelfend > $MINSVCDIST*$RLEN );
				addSV($svfdel[$#svfdel], $start, $end, $mstart, $mend, $dir, $rlen2, $mlen, $soft3, $RLEN/2, ord(substr($a[10], 15, 1))-33, $a[4], $nm);
				$svdelfend = $end;
			    } else {
				push(@svrdel, { cnt => 0 } ) if ( @svrdel == 0 || $start - $svdelrend > $MINSVCDIST*$RLEN );
				addSV($svrdel[$#svrdel], $start, $end, $mstart, $mend, $dir, $rlen2, $mlen, $soft5, $RLEN/2, ord(substr($a[10], 15, 1))-33, $a[4], $nm);
				$svdelrend = $end;
			    }
			    adddisccnt( $svfdel[$#svfdel] ) if ( @svfdel && abs($start - $svdelfend) <= $MINSVCDIST*$RLEN );
			    adddisccnt( $svrdel[$#svrdel] ) if ( @svrdel && abs($start - $svdelrend) <= $MINSVCDIST*$RLEN );
			    adddisccnt( $svfdup[$#svfdup] ) if ( @svfdup && abs($start - $svdupfend) <= $MIN_D );
			    adddisccnt( $svrdup[$#svrdup] ) if ( @svrdup && abs($start - $svduprend) <= $MIN_D );
			    adddisccnt( $svfinv5[$#svfinv5] ) if ( @svfinv5 && abs($start - $svinvfend5) <= $MIN_D );
			    adddisccnt( $svrinv5[$#svrinv5] ) if ( @svrinv5 && abs($start - $svinvrend5) <= $MIN_D );
			    adddisccnt( $svfinv3[$#svfinv3] ) if ( @svfinv3 && abs($start - $svinvfend3) <= $MIN_D );
			    adddisccnt( $svrinv3[$#svrinv3] ) if ( @svrinv3 && abs($start - $svinvrend3) <= $MIN_D );
			}
		    } elsif ( $dir * $mdir == -1 && $dir * $mlen < 0 ) { # duplication
			if ( $dir == 1 ) {
			    push(@svfdup, { cnt => 0 } ) if ( @svfdup == 0 || $start - $svdupfend > $MINSVCDIST*$RLEN );
			    addSV($svfdup[$#svfdup], $start, $end, $mstart, $mend, $dir, $rlen2, $mlen, $soft3, $RLEN/2, ord(substr($a[10], 15, 1))-33, $a[4], $nm);
			    $svdupfend = $end;
			} else {
			    push(@svrdup, { cnt => 0 } ) if ( @svrdup == 0 || $start - $svduprend > $MINSVCDIST*$RLEN );
			    addSV($svrdup[$#svrdup], $start, $end, $mstart, $mend, $dir, $rlen2, $mlen, $soft5, $RLEN/2, ord(substr($a[10], 15, 1))-33, $a[4], $nm);
			    $svduprend = $end;
			}
			$svfdup[$#svfdup]->{ disc }++ if ( @svfdup && abs($start - $svdupfend) <= $MINSVCDIST*$RLEN );
			$svrdup[$#svrdup]->{ disc }++ if ( @svrdup && abs($start - $svduprend) <= $MINSVCDIST*$RLEN );
			adddisccnt( $svfdel[$#svfdel] ) if ( @svfdel && abs($start - $svdelfend) <= $MIN_D );
			adddisccnt( $svrdel[$#svrdel] ) if ( @svrdel && abs($start - $svdelrend) <= $MIN_D );
			adddisccnt( $svfinv5[$#svfinv5] ) if ( @svfinv5 && abs($start - $svinvfend5) <= $MIN_D );
			adddisccnt( $svrinv5[$#svrinv5] ) if ( @svrinv5 && abs($start - $svinvrend5) <= $MIN_D );
			adddisccnt( $svfinv3[$#svfinv3] ) if ( @svfinv3 && abs($start - $svinvfend3) <= $MIN_D );
			adddisccnt( $svrinv3[$#svrinv3] ) if ( @svrinv3 && abs($start - $svinvrend3) <= $MIN_D );
		    } elsif ( $dir * $mdir == 1 ) { # Inversion
			if ( $dir == 1 && $mlen ) {
			    if ( $mlen < -3 * $RLEN ) {
				push(@svfinv3, { cnt => 0 } ) if ( @svfinv3 == 0 || $start - $svinvfend3 > $MINSVCDIST*$RLEN );
				addSV($svfinv3[$#svfinv3], $start, $end, $mstart, $mend, $dir, $rlen2, $mlen, $soft3, $RLEN/2, ord(substr($a[10], 15, 1))-33, $a[4], $nm);
				$svinvfend3 = $end;
				$svfinv3[$#svfinv3]->{ disc }++;
			    } elsif ( $mlen > 3 * $RLEN ) {
				push(@svfinv5, { cnt => 0 } ) if ( @svfinv5 == 0 || $start - $svinvfend5 > $MINSVCDIST*$RLEN );
				addSV($svfinv5[$#svfinv5], $start, $end, $mstart, $mend, $dir, $rlen2, $mlen, $soft3, $RLEN/2, ord(substr($a[10], 15, 1))-33, $a[4], $nm);
				$svinvfend5 = $end;
				$svfinv5[$#svfinv5]->{ disc }++;
			    }
			} elsif ( $mlen ) {
			    if ( $mlen < -3 * $RLEN ) {
				push(@svrinv3, { cnt => 0 } ) if ( @svrinv3 == 0 || $start - $svinvrend3 > $MINSVCDIST*$RLEN );
				addSV($svrinv3[$#svrinv3], $start, $end, $mstart, $mend, $dir, $rlen2, $mlen, $soft5, $RLEN/2, ord(substr($a[10], 15, 1))-33, $a[4], $nm);
				$svinvrend3 = $end;
				$svrinv3[$#svrinv3]->{ disc }++;
			    } elsif ( $mlen > 3 * $RLEN ) {
				push(@svrinv5, { cnt => 0 } ) if ( @svrinv5 == 0 || $start - $svinvrend5 > $MINSVCDIST*$RLEN );
				addSV($svrinv5[$#svrinv5], $start, $end, $mstart, $mend, $dir, $rlen2, $mlen, $soft5, $RLEN/2, ord(substr($a[10], 15, 1))-33, $a[4], $nm);
				$svinvrend5 = $end;
				$svrinv5[$#svrinv5]->{ disc }++;
			    }
			}
			if ( $mlen ) {
			    adddisccnt( $svfdel[$#svfdel] ) if ( @svfdel && ($start - $svdelfend) <= $MIN_D );
			    adddisccnt( $svrdel[$#svrdel] ) if ( @svrdel && ($start - $svdelrend) <= $MIN_D );
			    adddisccnt( $svfdup[$#svfdup] ) if ( @svfdup && ($start - $svdupfend) <= $MIN_D );
			    adddisccnt( $svrdup[$#svrdup] ) if ( @svrdup && ($start - $svduprend) <= $MIN_D );
			    #adddisccnt( $svfinv5[$#svfinv5] ) if ( @svfinv5 && abs($start - $svinvfend5) <= $MINSVCDIST*$RLEN );
			    #adddisccnt( $svrinv5[$#svrinv5] ) if ( @svrinv5 && abs($start - $svinvrend5) <= $MINSVCDIST*$RLEN );
			    #adddisccnt( $svfinv3[$#svfinv3] ) if ( @svfinv3 && abs($start - $svinvfend3) <= $MINSVCDIST*$RLEN );
			    #adddisccnt( $svrinv3[$#svrinv3] ) if ( @svrinv3 && abs($start - $svinvrend3) <= $MINSVCDIST*$RLEN );
			}
		    }
		} else { # Inter-chr translocation
		    # to be implemented
		    my $mchr = $a[6];
		    if ( $dir == 1 ) {
			push(@{ $svffus{ $mchr } }, { cnt => 0 } ) if ( (! $svffus{ $mchr } ) || $start - $svfusfend{ $mchr } > $MINSVCDIST*$RLEN );
			my $svn = @{ $svffus{ $mchr } } - 1;
			addSV($svffus{ $mchr }->[$svn], $start, $end, $mstart, $mend, $dir, $rlen2, 0, $soft3, $RLEN/2, ord(substr($a[10], 15, 1))-33, $a[4], $nm);
			$svfusfend{ $mchr } = $end;
			$svffus{ $mchr }->[$svn]->{ disc }++;
		    } else {
			push(@{ $svrfus{ $mchr } }, { cnt => 0 } ) if ( (! $svrfus{ $mchr } ) || $start - $svfusrend{ $mchr } > $MINSVCDIST*$RLEN );
			my $svn = @{ $svrfus{ $mchr } } - 1;
			addSV($svrfus{ $mchr }->[$svn], $start, $end, $mstart, $mend, $dir, $rlen2, 0, $soft5, $RLEN/2, ord(substr($a[10], 15, 1))-33, $a[4], $nm);
			$svfusrend{ $mchr } = $end;
			$svrfus{ $mchr }->[$svn]->{ disc }++;
		    }
		    adddisccnt( $svfdel[$#svfdel] ) if ( @svfdel && ($start - $svdelfend) <= 25 );
		    adddisccnt( $svrdel[$#svrdel] ) if ( @svrdel && ($start - $svdelrend) <= 25 );
		    adddisccnt( $svfdup[$#svfdup] ) if ( @svfdup && ($start - $svdupfend) <= 25 );
		    adddisccnt( $svrdup[$#svrdup] ) if ( @svrdup && ($start - $svduprend) <= 25 );
		    adddisccnt( $svfinv5[$#svfinv5] ) if ( @svfinv5 && ($start - $svinvfend5) <= 25 );
		    adddisccnt( $svrinv5[$#svrinv5] ) if ( @svrinv5 && ($start - $svinvrend5) <= 25 );
		    adddisccnt( $svfinv3[$#svfinv3] ) if ( @svfinv3 && ($start - $svinvfend3) <= 25 );
		    adddisccnt( $svrinv3[$#svrinv3] ) if ( @svrinv3 && ($start - $svinvrend3) <= 25 );
		}
	    }

	    for(my $ci = 0; $ci < @cigar; $ci += 2) {
		last if ($opt_u && $dir == 1 && $start >= $a[7]);
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
			# ignore large soft clip due to chimeric reads in library construction
			if ( (! $CHIMERIC) && $m >= 20 && $a[11] && $a[11] =~ /SA:Z:(\S+)/ ) {
			    my ($sachr, $sapos, $sadir, $sacig, $saq, $saunk) = split(/,/, $1);
			    $sadir = $sadir eq "+" ? 1 : -1;
			    if ( $dir * $sadir == -1 && $sachr eq $a[2] && abs($sapos - $a[3]) < 2*$RLEN && $sacig =~ /^\d\dS/ ) {
				$n += $m;
				$offset = 0;
				$start = $a[3];  # had to reset the start due to softclipping adjustment
				print STDERR "$a[0] @a[2..5] is ignored as chimeric with SA:$sapos,$sadir,$sacig\n" if ( $opt_y );
				next;
			    }
			}
			# align softclipped but matched sequences due to mis-softclipping
			while( $m-1 >= 0 && $start - 1 > 0 && $start - 1 <= $CHRS{ $chr } && $REF->{ $start-1 } && $REF->{ $start-1 } eq substr($a[9], $m-1, 1) && ord(substr($a[10],$m-1, 1))-33 > 10) {
			    if ( $start - 1 >= $START && $start - 1 <= $END ) {
				$hash->{ $start - 1 }->{ $REF->{ $start - 1 } }->{ cnt } = 0 unless( $hash->{ $start - 1 }->{ $REF->{ $start - 1 } }->{ cnt } );
				addCnt($hash->{ $start - 1 }->{ $REF->{ $start - 1 } }, $dir, $m, ord(substr($a[10],$m-1, 1))-33, $a[4], $nm);
				$cov->{ $start - 1 }++;
			    }
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
			    if ( $qn >= 1 && $qn > $lowqcnt && $start >= $START && $start <= $END ) {
				for(my $si = $m-1; $m - $si <= $qn; $si--) {
				    $sclip5->{ $start }->{ nt }->[$m-1-$si]->{ substr($a[9], $si, 1) }++;
				    $sclip5->{ $start }->{ seq }->[$m-1-$si]->{ substr($a[9], $si, 1) }->{ cnt } = 0 unless( $sclip5->{ $start }->{ seq }->[$m-1-$si]->{ substr($a[9], $si, 1) }->{ cnt });
				    addCnt($sclip5->{ $start }->{ seq }->[$m-1-$si]->{ substr($a[9], $si, 1) }, $dir, $si - ($m-$qn), ord(substr($a[10], $si, 1))-33, $a[4], $nm);
				}
				$sclip5->{ $start }->{ cnt } = 0 unless( $sclip5->{ $start }->{ cnt } );
				addCnt( $sclip5->{ $start }, $dir, $m, $q/$qn, $a[4], $nm);
			    }
			}
			$m = $cigar[$ci];
		    } elsif ( $ci == @cigar - 2 ) { # 3' soft clipped
			# ignore large soft clip due to chimeric reads in library construction
			if ( (! $CHIMERIC) && $m >= 20 && $a[11] && $a[11] =~ /SA:Z:(\S+)/ ) {
			    my ($sachr, $sapos, $sadir, $sacig, $saq, $saunk) = split(/,/, $1);
			    $sadir = $sadir eq "+" ? 1 : -1;
			    if ( $dir * $sadir == -1 && $sachr eq $a[2] && abs($sapos - $a[3]) < 2*$RLEN && $sacig =~ /\d\dS$/ ) {
				$n += $m;
				$offset = 0;
				$start = $a[3];  # had to reset the start due to softclipping adjustment
				print STDERR "$a[0] @a[2..5] is ignored as chimeric with SA:$sapos,$sadir,$sacig\n" if ( $opt_y );
				next;
			    }
			}
			while( $n < length($a[9]) && $REF->{ $start } && $REF->{ $start } eq substr($a[9], $n, 1) && ord(substr($a[10], $n, 1))-33 > 10) {
			    if ( $start >= $START && $start <= $END ) {
				$hash->{ $start }->{ $REF->{ $start } }->{ cnt } = 0 unless( $hash->{ $start }->{ $REF->{ $start } }->{ cnt } );
				addCnt($hash->{ $start }->{ $REF->{ $start } }, $dir, $rlen2-$p, ord(substr($a[10], $n, 1))-33, $a[4], $nm);
				$cov->{$start}++;
			    }
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
			    if ( $qn >= 1 && $qn > $lowqcnt && $start >= $START && $start <= $END ) {
				for(my $si = 0; $si < $qn; $si++) {
				    $sclip3->{ $start }->{ nt }->[$si]->{ substr($a[9], $n+$si, 1) }++;
				    $sclip3->{ $start }->{ seq }->[$si]->{ substr($a[9], $n+$si, 1) }->{ cnt } = 0 unless( $sclip3->{ $start }->{ seq }->[$si]->{ substr($a[9], $n+$si, 1) }->{ cnt } );
				    addCnt($sclip3->{ $start }->{ seq }->[$si]->{ substr($a[9], $n+$si, 1) }, $dir, $qn - $si, ord(substr($a[10], $n+$si, 1))-33, $a[4], $nm);
				}
				$sclip3->{ $start }->{ cnt } = 0 unless( $sclip3->{ $start }->{ cnt } );
				addCnt( $sclip3->{ $start }, $dir, $m, $q/$qn, $a[4], $nm);
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
		    if ( $opt_k &&  $cigar[$ci+2] && $cigar[$ci+2] <= $VEXT && $cigar[$ci+3] eq "M" && $cigar[$ci+5] && $cigar[$ci+5] =~ /[ID]/ ) {
			$s .= "#" . substr($a[9], $n+$m, $cigar[$ci+2]);
			$q .= substr($a[10], $n+$m, $cigar[$ci+2]); 
			$s .= $cigar[$ci+5] eq "I" ? ("^" . substr($a[9], $n+$m+$cigar[$ci+2], $cigar[$ci+4])) : ("^" . $cigar[$ci+4]);
			$q .= $cigar[$ci+5] eq "I" ? substr($a[10], $n+$m+$cigar[$ci+2], $cigar[$ci+4]) : substr($a[10], $n+$m+$cigar[$ci+2], 1);
			$multoffs += $cigar[$ci+2] + ($cigar[$ci+5] eq "D" ? $cigar[$ci+4] : 0);
			$multoffp += $cigar[$ci+2] + ($cigar[$ci+5] eq "I" ? $cigar[$ci+4] : 0);
			if ( $cigar[$ci+6] && $cigar[$ci+7] eq "M" ) {
			    my ($toffset, $tss, $tq) = findOffset($start+$multoffs, $n+$m+$multoffp, $cigar[$ci+6], \$a[9], \$a[10], $REF, $cov);
			    ($offset, $ss) = ($toffset, $tss);
			    $q .= $tq;
			}
			$ci += 4;
		    } else {
			if ( $opt_k && $cigar[$ci+3] && $cigar[$ci+3] eq "M" ) {
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
				    $cov->{ $start + $osi }++;
				}
			    }
			}
		    }
		    $s = "$s&$ss" if ( $offset > 0);
		    if ( $start - 1 >= $START && $start -1 <= $END && $s !~ /N/ ) {
			my $inspos = $start - 1;
			if( $s =~ /^[ATGC]+$/ ) {
			    ($inspos, $s) = adjInsPos($start-1, $s, $REF);
			}
			$ins{ $inspos }->{ "+$s" }++;
			$hash->{ $inspos }->{ I }->{ "+$s" }->{ $dir }++;
			my $hv = $hash->{ $inspos }->{ I }->{ "+$s" };
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
			#if ( $REF->{ $inspos } && $hash->{ $inspos }->{ $REF->{ $inspos } } && substr($a[9], $n-1-($start-$inspos), 1) eq $REF->{ $inspos } ) {
			    #subCnt($hash->{ $inspos }->{ $REF->{ $inspos } }, $dir, $tp, $tmpq, $a[4], $nm);
			    subCnt($hash->{ $inspos }->{ substr($a[9], $n-1-($start-1-$inspos), 1) }, $dir, $tp, ord(substr($a[10], $n-1-($start-1-$inspos), 1))-33, $a[4], $nm - $nmoff) if ( $inspos > $a[3] );
			#}
			# Adjust count if the insertion is at the edge so that the AF won't > 1
			if ( $ci == 2 && ($cigar[1] eq "S" || $cigar[1] eq "H") ) {
			    my $ttref = $hash->{ $inspos }->{ $REF->{ $inspos } };
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
			    #$cov->{ $inspos }->{ $REF->{ $inspos } }++;
			    $cov->{ $inspos }++;
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
		    if ( $opt_k && $cigar[$ci+2] && $cigar[$ci+2] <= $VEXT && $cigar[$ci+3] eq "M" && $cigar[$ci+5] && $cigar[$ci+5] =~ /[ID]/ ) {
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
			if ( $opt_k && $cigar[$ci+3] && $cigar[$ci+3] eq "M" ) {
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
			$hash->{ $start }->{ $s }->{ $dir }++;
			$dels5{ $start }->{ $s }++;
			my $hv = $hash->{ $start }->{ $s };
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
			    #$cov->{ $start+$i }->{ $s }++;
			    $cov->{ $start+$i }++;
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
			if ( $dir == 1 ) {
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
			    my ($toffset, $tss, $tq, $tnmoff) = findOffset($start+$ddlen+1, $n+1, $cigar[$ci+2], \$a[9], \$a[10], $REF, $cov);
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
			    $hash->{ $start - $qbases + 1 }->{ $s }->{ $dir }++;
			    $mnp{ $start - $qbases + 1 }->{ $s }++ if ( $s =~ /^[ATGC]&[ATGC]+$/ );
			    my $hv = $hash->{ $start - $qbases + 1 }->{ $s };
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
				$cov->{ $start - $qi + 1 }++;
			    }
			    if ( $s =~ /^-/ ) {
				$dels5{ $start - $qbases + 1 }->{ $s }++;
				for(my $qi = 1; $qi < $ddlen; $qi++) {
				    $cov->{ $start + $qi }++;
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
		    last if ($opt_u && $dir == 1 && $start >= $a[7]);
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
	print STDERR "TIME: Finish parsing SAM: ", time(), " $bami $chr:$START-$END\n" if ( $opt_y );
    }

    if ( $opt_i ) {
	while( my ($intron, $icnt) = each %SPLICE ) {
	    print join("\t", $sample, $chr, $intron, $icnt), "\n";
	}
	return;
    }
    return ($hash, $cov) if ( $svflag );
    #use Object; print STDERR Object::Perl(\@svfdel), Object::Perl(\@svrdel);# if ( $opt_y );
    #use Object; print STDERR Object::Perl(\@svfinv3), Object::Perl(\@svrinv3);
    #use Object; print STDERR Object::Perl(\@svfinv5), Object::Perl(\@svrinv5);
    print STDERR "TIME: Start realign: ", time(), "\n" if ( $opt_y );
    unless($opt_U) {
	filterSV(\@svfinv3);
	filterSV(\@svrinv3);
	filterSV(\@svfinv5);
	filterSV(\@svrinv5);
	filterSV(\@svfdel);
	filterSV(\@svrdel);
	filterSV(\@svfdup);
	filterSV(\@svrdup);
	while( my ($mchr, $svv) = each %svffus ) {
	    filterSV($svv);
	}
	while( my ($mchr, $svv) = each %svrfus ) {
	    filterSV($svv);
	}
	while(my ($psv, $psvv) = each %SOFTP2SV) {
	    $psvv = [sort {$b->{ cnt } <=> $a->{ cnt }} @$psvv];
	}
    }
    adjMNP($hash, \%mnp, $cov, $chr, $REF, $sclip5, $sclip3);
    if ( $opt_k ) {
	print STDERR "Start Realigndel\n" if ( $opt_y );
	realigndel($hash, \%dels5, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $START, $END);
	print STDERR "\n\nStart Realignins\n" if ( $opt_y );
	realignins($hash, \%ins, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $START, $END);
	print STDERR "\n\nStart Realignlgdel\n" if ( $opt_y );
	realignlgdel($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $START, $END, \@svfdel, \@svrdel);
	print STDERR "\n\nStart Realignlgins30\n" if ( $opt_y );
	realignlgins30($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $START, $END);
	print STDERR "\n\nStart Realignlgins\n" if ( $opt_y );
	realignlgins($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $START, $END, \@svfdup, \@svrdup);
    }
    unless($opt_U) {
	print STDERR "\n\nStart Structral Variants: DEL\n" if ( $opt_y );
	findDEL($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, \@svfdel, \@svrdel);
	print STDERR "\n\nStart Structral Variants: INV\n" if ( $opt_y );
	findINV($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, \@svfinv3, \@svrinv3, \@svfinv5, \@svrinv5);
	print STDERR "\n\nStart Structral Variants\n" if ( $opt_y );
	findsv($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, \@svfdel, \@svrdel, \@svfdup, \@svrdup, \@svfinv5, \@svrinv5, \@svfinv3, \@svrinv3);
	print STDERR "\n\nStart Structral Variants: DEL discordant pairs only\n" if ( $opt_y );
	findDELdisc($hash, $cov, $REF, $chr, $bams, \@svfdel, \@svrdel, $sclip5, $sclip3);
	print STDERR "\n\nStart Structral Variants: INV discordant pairs only\n" if ( $opt_y );
	#findINVdisc($hash, $cov, $REF, $chr, $bams, \@svfinv3, \@svrinv3, \@svfinv5, \@svrinv5, $sclip5, $sclip3);
	print STDERR "\n\nStart Structral Variants: DUP discordant pairs only\n" if ( $opt_y );
	findDUPdisc($hash, $cov, $REF, $chr, $bams, \@svfdup, \@svrdup, $sclip5, $sclip3);
    }
    outputClipping($sclip5, $sclip3) if ( $opt_y );
    print STDERR "TIME: Finish realign: ", time(), "\n" if ( $opt_y );
    return ($hash, $cov, $opt_t && $totalreads ? sprintf("%.3f", $dupreads/$totalreads) : 0);
}

sub toVars {
    my ($chr, $START, $END, $bam, $REF) = @_;
    my @bams = $bam =~ /^http/ ? ($bam) : split(/:/, $bam);
    my ($hash, $cov, $duprate) = parseSAM($chr, $START, $END, \@bams, $REF);
    my %vars; # the variant structure
    while( my ($p, $v) = each %$hash ) {
	next unless( %$v );
	unless( $v->{ SV } ) {
	    next unless( $p >= $START && $p <= $END); 
	    next unless ( $cov->{ $p } );
	}
	my @tmp = ();
	my $vcov = 0; #the variance coverage
	my @var = ();
	my @vn = keys %$v;
	if ( @vn == 1 && $vn[0] eq $REF->{ $p } ) {
	    next unless( $opt_p || $BAM2 || $opt_a ); # ignore if only reference were seen and no pileup to avoid computation
	}
	#my @v = values %{ $cov->{ $p } };
	#my $tcov = 0; $tcov += $_ foreach(@v); #$stat->sum(\@v);
	my $tcov = $cov->{ $p };
	my $hicov = 0;
	#$hicov += $_->{ hicnt } ? $_->{ hicnt } : 0 foreach( values %$v );
	print STDERR "Error tcov: $chr $p $START $END $v->{ SV }->{ type }\n" unless( $tcov );
	next if ( $tcov == 0 ); # ignore when there's no coverage
	while( my ($n, $cnt) = each %$v ) {
	    next if ( $n eq "SV" );
	    if ( $n eq "I" ) {
		while( my ($in, $icnt) = each %$cnt ) {
		    $hicov += $icnt->{ hicnt } ? $icnt->{ hicnt } : 0;
		}
	    } else {
		$hicov += $cnt->{ hicnt } ? $cnt->{ hicnt } : 0;
	    }
	}
	while( my ($n, $cnt) = each %$v ) {
	    if ( $n eq "SV" ) {
		$vars{ $p }->{ SV } = "$cnt->{ splits }-$cnt->{ pairs }-$cnt->{ clusters }";
		next;
	    }
	    unless( $n eq "I") {
		next unless( $cnt->{ cnt } );
		my $fwd = $cnt->{ 1 } ? $cnt->{ 1 } : 0;
		my $rev = $cnt->{ -1 } ? $cnt->{ -1 } : 0;
		my $bias = strandBias($fwd, $rev);
		my $vqual = sprintf("%.1f", $cnt->{ qmean }/$cnt->{ cnt }); # base quality
		my $MQ = sprintf("%.1f", $cnt->{ Qmean }/$cnt->{ cnt }); # mapping quality
		my ($hicnt, $locnt) = ($cnt->{ hicnt } ? $cnt->{ hicnt } : 0, $cnt->{ locnt } ? $cnt->{ locnt } : 0);
		my $ttcov = ( $cnt->{ cnt } > $tcov && $cnt->{ extracnt } && $cnt->{ cnt } - $tcov < $cnt->{ extracnt } ) ? $cnt->{ cnt } : $tcov;
		my $tvref = {n => $n, cov => $cnt->{ cnt }, fwd => $fwd, rev => $rev, bias => $bias, freq => $cnt->{ cnt }/$ttcov, pmean => sprintf("%.1f", $cnt->{ pmean }/$cnt->{ cnt } ), pstd => $cnt->{ pstd }, qual => $vqual, qstd => $cnt->{ qstd }, mapq => $MQ, qratio => sprintf("%.3f", $hicnt/($locnt ? $locnt : $locnt+0.5)), hifreq => ($hicov > 0 ? $hicnt/$hicov : 0), extrafreq => $cnt->{ extracnt } ? $cnt->{ extracnt }/$ttcov : 0, shift3 => 0, msi => 0, nm => sprintf("%.1f", $cnt->{ nm }/$cnt->{ cnt } ), hicnt => $hicnt, hicov => $hicov, duprate => $duprate };
		push(@var, $tvref);
		if ( $opt_D ) {
		    push( @tmp, "$n:" . ($fwd + $rev) . ":F-$fwd:R-$rev:" . sprintf("%.3f", $tvref->{freq}) . ":$tvref->{bias}:$tvref->{pmean}:$tvref->{pstd}:$vqual:$tvref->{qstd}:" . sprintf("%.3f", $tvref->{hifreq}) . ":$tvref->{mapq}:$tvref->{qratio}");
		}
	    }
	}
	if ( $v->{ I } ) {
	    while( my ($n, $cnt) = each %{ $v->{ I } } ) {
		my $fwd = $cnt->{ 1 } ? $cnt->{ 1 } : 0;
		my $rev = $cnt->{ -1 } ? $cnt->{ -1 } : 0;
		my $bias = strandBias($fwd, $rev);
		my $vqual = sprintf("%.1f", $cnt->{ qmean }/$cnt->{ cnt }); # base quality
		my $MQ = sprintf("%.1f", $cnt->{ Qmean }/$cnt->{ cnt }); # mapping quality
		my ($hicnt, $locnt) = ($cnt->{ hicnt } ? $cnt->{ hicnt } : 0, $cnt->{ locnt } ? $cnt->{ locnt } : 0);
		#$hicov += $hicnt ? $hicnt : 0;
		my $ttcov = ( $cnt->{ cnt } > $tcov && $cnt->{ extracnt } && $cnt->{ cnt } - $tcov < $cnt->{ extracnt } ) ? $cnt->{ cnt } : $tcov;
		if ( $ttcov < $cnt->{ cnt } ) {
		    $ttcov = $cnt->{ cnt };
		    if ( $cov->{ $p + 1 } && $ttcov < $cov->{ $p+1 } - $cnt->{ cnt } ) {
			$ttcov = $cov->{ $p + 1 };
			$hash->{ $p + 1 }->{ $REF->{ $p + 1 } }->{ 1 } -= $fwd; # Adjust the reference
			$hash->{ $p + 1 }->{ $REF->{ $p + 1 } }->{ -1 } -= $rev;
		    }
		    $tcov = $ttcov;
		}
		my $tvref = {n => $n, cov => $cnt->{ cnt }, fwd => $fwd, rev => $rev, bias => $bias, freq => $cnt->{ cnt }/$ttcov, pmean => sprintf("%.1f", $cnt->{ pmean }/$cnt->{ cnt } ), pstd => $cnt->{ pstd }, qual => $vqual, qstd => $cnt->{ qstd }, mapq => $MQ, qratio => sprintf("%.3f", $hicnt/($locnt ? $locnt : $locnt+0.5)), hifreq => ($hicov > 0 ? $hicnt/$hicov : 0), extrafreq => $cnt->{ extracnt } ? $cnt->{ extracnt }/$ttcov : 0, shift3 => 0, msi => 0, nm => sprintf("%.1f", $cnt->{ nm }/$cnt->{ cnt } ), hicnt => $hicnt, hicov => $hicov, duprate => $duprate };
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
	#use Object; print STDERR "G: '$genotype1' $p\n", Object::Perl($v) unless( $genotype1 );
	if ( $genotype1 =~ /^\+/ ) {
	    $genotype1 = $genotype1 =~ /<dup(\d+)/ ? ("+" . ($SVFLANK + $1)) : ("+" . (length($genotype1) - 1 ));
	}
	my $genotype2 = "";
	my $vn;

	if ( $vars{ $p }->{ REF } ) {
	    ($rfc, $rrc) = ($vars{ $p }->{ REF }->{ fwd }, $vars{ $p }->{ REF }->{ rev });
	}
	if ( $tcov > $cov->{ $p } && $hash->{ $p + 1 } && $hash->{ $p + 1 }->{ $REF->{ $p + 1 } } ) {
	    my $tpref = $hash->{ $p + 1 }->{ $REF->{ $p + 1 } };
	    ($rfc, $rrc) = ($tpref->{ 1 }, $tpref->{ -1 });
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
		    unless( $vn =~ /&/ || $vn =~ /#/ || $vn =~ /<dup/ ) {
			my $tseq1 = $vn;
			$tseq1 =~ s/^\+//;
			my $leftseq = join( "", (map { $REF->{ $_ }; } (($p-50 > 1 ? $p-50 : 1) .. $p)) ); # left 10 nt
			my $tseq2 = join("", (map { $REF->{ $_ }; } (($p+1) .. ($p+70 > $CHRS{ $chr } ? $CHRS{ $chr } : ($p+70)))));
			($msi, $shift3, $msint) = findMSI($tseq1, $tseq2, $leftseq);
			my ($tmsi, $tshift3, $tmsint) = findMSI($leftseq, $tseq2);
			($msi, $msint) = ($tmsi, $tmsint) if ( $msi < $tmsi ); # Don't change shift3
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
		    $varallele = "<DUP>" if ( length($varallele) > $SVMINLEN );
		    if ( $varallele =~ /<dup(\d+)/ ) {
			$ep = $sp + (2*$SVFLANK + $1) - 1;
			$genotype2 = "+" . (2*$SVFLANK + $1);
			$varallele = "<DUP>";
		    }
		} elsif ( $vn =~ /^-/ ) {
		    $varallele = $vn;
		    if ( $dellen < $SVMINLEN ) {
			$varallele =~ s/^-\d+//;
			my $leftseq = join( "", (map { $REF->{ $_ }; } (($p-70 > 1 ? $p - 70 : 1) .. ($p - 1))) ); # left 10 nt
			my $tseq = join("", (map { $REF->{ $_ }; } (($p) .. ($p+$dellen+70 > $CHRS{ $chr } ? $CHRS{ $chr } : ($p+$dellen+70)) )));
			($msi, $shift3, $msint) = findMSI(substr($tseq, 0, $dellen), substr($tseq, $dellen), $leftseq);
			my ($tmsi, $tshift3, $tmsint) = findMSI($leftseq, substr($tseq, $dellen));
			($msi, $msint) = ($tmsi, $tmsint) if ( $msi < $tmsi ); # Don't change shift3
			$msi = $shift3/$dellen unless( $msi > $shift3/$dellen ); # if ( $shift3%$dellen == 0 );
			if ( $vn =~ /<inv\d+/ ) {
			    $varallele = "<INV>";
			    $genotype2 = "<INV$dellen>";
			}
		    } elsif ( $vn =~ /^-\d+\^/ ) {
			$varallele = "<INV>";
			$genotype2 = "<INV$dellen>";
		    } else {
			$varallele = "<DEL>";
		    }
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
			$varallele = $REF->{ $p - 1 } unless( $varallele eq "<DEL>" );
			$refallele = $REF->{ $p - 1 };
			$sp--;
		    }
		    if ( $vn =~ /<(...)\d+>/ ) {
			#$refallele .= join("", map { $REF->{ $_ }; } ($p .. ($p+74))) . $1 . join("", map { $REF->{ $_ }; } (($p+$dellen-74) .. ($p+$dellen)));
			$refallele = $REF->{ $p };
			#$refallele = $REF->{ $p - 1 };
			#$sp--;
		    } elsif( $dellen < $SVMINLEN ) {
			for(my $i = 0; $i < $dellen; $i++) {
			    $refallele .= $REF->{ $p + $i };
			}
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
		    if ($varallele eq "<DEL>" && length($refallele) > 1) {
			$refallele = $REF->{ $sp };
			$tcov = $cov->{ $sp - 1 };
			$tcov = $vref->{ cov } if ( $vref->{ cov } > $tcov );
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
		if ( $opt_J && length($refallele) > 1 && length($varallele) > 1 ) { # fix 5' for complex in CRISPR mode
		    my $n = 0;
		    $n++ while ( length($refallele) > $n + 1 && length($varallele) > $n + 1 && substr($refallele, $n, 1) eq substr($varallele, $n, 1) );
		    if ( $n ) {
			$sp += $n;
			$refallele = substr($refallele, $n);
			$varallele = substr($varallele, $n);
		    }
		    # Let adjComplex to take care the 3'
		}

		if ( $opt_J && (length($refallele) != length($varallele) && substr($refallele, 0, 1) eq substr($varallele, 0, 1) ) ) {
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
		$vref->{ leftseq } = join( "", (map { $REF->{ $_ } ? $REF->{ $_ } : ""; } (($sp-20 < 1 ? 1 : $sp-20) .. ($sp - 1))) ); # left 20 nt
		$vref->{ rightseq } = join( "", (map { $REF->{ $_ } ? $REF->{ $_ } : ""; } (($ep+1) .. ($ep+20 > $CHRS{ $chr } ? $CHRS{ $chr } : $ep+20))) ); # right 20 nt
		$vref->{ extrafreq } = sprintf("%.4f", $vref->{ extrafreq } ) if ($vref->{ extrafreq });
		$vref->{ freq } = sprintf("%.4f", $vref->{ freq }) if ($vref->{ freq });
		$vref->{ hifreq } = sprintf("%.4f", $vref->{ hifreq }) if ($vref->{ hifreq });
		$vref->{ msi } = $msi ? sprintf("%.3f", $msi) : $msi;
		$vref->{ msint } = $msint ? length($msint) : 0;
		$vref->{ shift3 } = $shift3;
		$vref->{ sp } = $sp;
		#$vref->{ ep } = $genotype =~ /\+(\d+)/ ? $sp + $1 : $ep;
		$vref->{ ep } = $ep;
		$vref->{ refallele } = $refallele;
		$vref->{ varallele } = $varallele;
		$vref->{ genotype } = $genotype;
		$vref->{ tcov } = $tcov;
		$vref->{ rfc } = $rfc;
		$vref->{ rrc } = $rrc;
		$vref->{ bias } = $vars{$p}->{ REF } ? $vars{$p}->{ REF }->{ bias } . ";" . $vref->{ bias } : "0;" . $vref->{ bias };
		$vref->{ DEBUG } = join(" & ", @tmp) if ( $opt_D );

	    }
	} elsif ( $vars{$p}->{ REF } ) {
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
	    $vref->{ duprate } = $duprate;
	}
    }
    print STDERR "TIME: Finish preparing vars: ", time(), "\n" if ( $opt_y );
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
	    print STDERR "  AdjMnt: $p $vn $hash->{ $p }->{ $vn }->{ cnt }\n" if ( $opt_y );
	    my $mnt = $vn; $mnt =~ s/&//;
	    for(my $i = 0; $i < length($mnt)-1; $i++) {
		my $left = substr($mnt, 0, $i+1);
		substr($left, 1, 0) = "&" if ( length($left) > 1 );
		my $right = substr($mnt, -(length($mnt)-$i-1)); 
		substr($right, 1, 0) = "&" if ( length($right) > 1 );
		if ( $hash->{ $p }->{ $left } ) {
		    next unless( $hash->{ $p }->{ $left } && $hash->{ $p }->{ $left }->{ cnt } > 0 );
		    my $tref = $hash->{ $p }->{ $left };
		    if ($tref->{ cnt } < $hash->{ $p }->{ $vn }->{ cnt } && $tref->{ pmean }/$tref->{ cnt } <= $i + 1) {
			print STDERR "    AdjMnt Left: $p $vn Left: $left Cnt: $tref->{ cnt }\n" if ( $opt_y );
			adjCnt($hash->{ $p }->{ $vn }, $tref);
			delete $hash->{ $p }->{ $left };
		    }
		}
		if ( $hash->{ $p+$i+1 }->{ $right } ) {
		    next unless( $hash->{ $p+$i+1 }->{ $right }  && $hash->{ $p+$i+1 }->{ $right }->{ cnt } > 0);
		    my $tref = $hash->{ $p+$i+1 }->{ $right };
		    if ($tref->{ cnt } < $hash->{ $p }->{ $vn }->{ cnt } ) { #&& $tref->{ pmean }/$tref->{ cnt } <= length($mnt)-$i-1) {
			print STDERR "    AdjMnt Right: $p $vn Right: $right Cnt: $tref->{ cnt }\n" if ( $opt_y );
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
	my $curmsi = length($msimatch)/$nmsi;
	my $curshift3 = 0;
	if ( $tseq2 =~ /^(($msint)+)/ ) {
	    $curmsi += length($1)/$nmsi;
	}
	($maxmsi, $msicnt) = ($msint, $curmsi) if ( $curmsi > $msicnt );
	$nmsi++;
    }
    my $tseq = "$tseq1$tseq2";
    $shift3++ while(substr($tseq, $shift3, 1) eq substr($tseq2, $shift3, 1) && $shift3 < length($tseq2) );
    #print STDERR "$tseq1 $tseq2 $maxmsi $msicnt $shift3\n" if ( $tseq2 =~ /^GGCGCAGC/ );
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
    my $ntcnt = 0;
    my $a = $seq =~ s/A/A/g;
    $ntcnt++ if ( $a > 0 );
    return 1 if ( $a/$len > 0.75 );
    my $t = $seq =~ s/T/T/g;
    $ntcnt++ if ( $t > 0 );
    return 1 if ( $t/$len > 0.75 );
    my $g = $seq =~ s/G/G/g;
    $ntcnt++ if ( $g > 0 );
    return 1 if ( $g/$len > 0.75 );
    my $c = $seq =~ s/C/C/g;
    $ntcnt++ if ( $c > 0 );
    return 1 if ( $c/$len > 0.75 );
    return 1 if ( $ntcnt < 3 );
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
	next if ($sc5v->{ used });
	for( my $j = 0; $j < @tmp3; $j++ ) {
	    my ($p3, $sc3v, $cnt3) = @{ $tmp3[$j] };
	    last if ($sc5v->{ used });
	    next if ($sc3v->{ used });
	    next if ( $p5 - $p3 > $RLEN*2.5 );
	    next if ( $p3 - $p5 > $RLEN-10 );  # if they're too far away, don't even try
	    my $seq5 = findconseq($sc5v, 5);
	    my $seq3 = findconseq($sc3v, 3);
	    next unless(length($seq5) > 10 && length($seq3) > 10);
	    print STDERR "  Working lgins30: $p3 $p5 3: $seq3 $cnt3 5: ", scalar reverse($seq5), " $cnt5\n" if ( $opt_y );
	    next unless( $cnt5/$cnt3 >= 0.08 && $cnt5/$cnt3 <= 12 ); # Don't even try if there're extreme bias
	    my ($bp5, $bp3, $score) = find35match($seq5, $seq3, $p5, $p3, $REF);
	    next unless($score);
	    my $smscore = int($score/2); # to ensure higher quality of bases are used as read ends are usually poor quality
	    my $ins = $bp3 + $smscore > 1 ? substr($seq3, 0, -($bp3+$smscore) + 1) : $seq3;
	    $ins .= reverse(substr($seq5, 0, $bp5+$smscore)) if ( $bp5+$smscore > 0 );
	    if ( islowcomplexseq($ins) ) {
		print STDERR "  Discard low complexity insertion found $ins.\n" if ( $opt_y );
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
		    my $rpt = 2;
		    my $tnr = 3;
		    while( (($p3 - $p5 + length($ins))/$tnr)/length($ins) > 1 ) {
			$rpt++ if (($p3 - $p5 + length($ins))%$tnr == 0);
			$tnr++;
		    }
		    for(my $p = $p5; $p < $p5 + ($p3 - $p5 + length($ins))/$rpt - length($ins); $p++) {
			$tmp .= $REF->{ $p };
		    }
		    $ins = "+$tmp$ins";
		} else {
		    for(my $p = $p5; $p < $p3; $p++) {
			$tmp .= $REF->{ $p };
		    }
		    if( (length($ins)-length($tmp))%2 == 0 ) {
			my $tex = (length($ins)-length($tmp))/2;
			$ins = ($tmp . substr($ins, 0, $tex)) eq substr($ins, $tex) ? ("+" . substr($ins, $tex)) : "+$tmp$ins";
		    } else {
			$ins = "+$tmp$ins";
		    }
		}
		print STDERR "  Found lgins30: $p3 $p5 ", length($ins), " $tmp + $ins\n" if ( $opt_y );
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

sub findINV {
    my ($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $svfinv3, $svrinv3, $svfinv5, $svrinv5) = @_;
    findINVsub($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $svfinv5, 1, 5);
    findINVsub($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $svrinv5, -1, 5);
    findINVsub($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $svfinv3, 1, 3);
    findINVsub($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $svrinv3, -1, 3);
}

# Would only consider those with supports from both orientations
#  (svfinv5) --> | <-- (svrinv5) ..... (svfinv3) --> | <-- (svrinv3)
sub findINVsub {
    my ($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $svref, $dir, $side) = @_; # $dir == 1 means 3' soft clipping
    foreach my $inv (@$svref) {
	next if ( $inv->{ used } );
	my ($start, $end, $ms, $me, $cnt, $mlen) = map { $inv->{ $_ }; } qw(start end mstart mend cnt mlen);
	next unless( $cnt >= $MINR );
	my @soft = ();
	while( my ($sp, $v) = each %{ $inv->{ soft } } ) {
	    push(@soft, [$sp, $v]);
	}
	@soft = sort { $b->[1] <=> $a->[1]; } @soft;
	my $softp = @soft > 0 ? $soft[0]->[0] : 0;
	my $sclip = $dir == 1 ? $sclip3 : $sclip5;
	print STDERR "\n\nWorking INV $softp $dir $side pair_cnt: $cnt\n" if ( $opt_y );
	unless( isLoaded($chr, $ms, $me, $REF) ) {
	    getREF($chr, $ms, $me, $REF, 500);
	    parseSAM($chr, $ms-200, $me+200, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1);
	}
	my $bp = 0;
	my $scv = undef;
	my $seq = "";
	my $EXTRA = "";
	if ( $softp ) {
	    next unless( $sclip->{ $softp } );
	    $scv = $sclip->{ $softp };
	    next if ($scv->{ used });
	    $seq = findconseq($scv);
	    next unless( $seq );
	    ($bp, $EXTRA) = findMatchRev($seq, $REF, $softp, $dir);
	    ($bp, $EXTRA) = findMatchRev($seq, $REF, $softp, $dir, $SEED2, 0) unless( $bp );
	    next unless( $bp );
	} else { # Look within 100bp to see whether a soft cliping can be found but not associated with discordant pairs
	    my $sp = $dir == 1 ? $end : $start; # starting position
	    for(my $i = 1; $i <= 2 * $RLEN; $i++) {
		my $cp = $sp + $i * $dir;
		next unless($sclip->{ $cp });
		$scv = $sclip->{ $cp };
		next if ( $scv->{ used } );
		$seq = findconseq($scv);
		next unless( $seq );
		($bp, $EXTRA) = findMatchRev($seq, $REF, $cp, $dir);
		($bp, $EXTRA) = findMatchRev($seq, $REF, $cp, $dir, $SEED2, 0) unless( $bp );
		next unless( $bp );
		$softp = $cp;
		last if ( ($dir == 1 && abs($bp - $me) < $MINSVCDIST*$RLEN) || ($dir == -1 && abs($bp - $ms) < $MINSVCDIST*$RLEN) );
	    }
	    next unless( $bp );
	}
	print STDERR "    $softp $bp $dir $side $seq pair_cnt: $cnt soft_cnt: $scv->{cnt}\n" if ( $opt_y );
	if ( $side == 5 ) {
	    $bp-- if ( $dir == -1 );
	} else {
	    $dir == 1 ? ($bp++ && $softp--) : $softp--;
	}
	($bp, $softp) = ($softp, $bp) if ( $side == 3 );
	if ( ($dir == -1 && $side == 5) || ($dir == 1 && $side == 3) ) {
	    $softp++ && $bp-- while( $REF->{ $softp } && $REF->{ $softp } eq $REVCOMP{ $REF->{ $bp } } );
	}
	$softp-- && $bp++ while( $REF->{ $softp-1 } && $REF->{ $softp-1 } eq $REVCOMP{ $REF->{ $bp+1 } } );
	if ( $bp > $softp && $bp - $softp > 150 && ($bp - $softp)/abs($mlen) < 1.5 ) {
	    my $len = $bp - $softp + 1;
	    my $ins5 = reverse(join("", map { $REF->{ $_ } ? $REF->{ $_ } : "" } (($bp-$SVFLANK+1) .. $bp)));
	    my $ins3 = reverse(join("", map { $REF->{ $_ } ? $REF->{ $_ } : "" } ($softp .. ($softp+$SVFLANK-1))));
	    $ins3 =~ y/ATGC/TACG/;
	    $ins5 =~ y/ATGC/TACG/;
	    my $ins = "$ins5<inv" . ($len - 2*$SVFLANK) . ">$ins3";
	    if ( $dir == 1 && $EXTRA ) {
		$EXTRA =~ y/ATGC/TACG/;
		$EXTRA = reverse($EXTRA);
		$ins = $EXTRA . $ins;
	    } elsif ( $dir == -1 && $EXTRA ) {
		$ins .= $EXTRA;
	    }
	    my $gt = "-$len^$ins";
	    $hash->{ $softp }->{ $gt }->{ cnt } = 0 unless( $hash->{ $softp }->{ $gt }->{ cnt } );
	    my $vref = $hash->{ $softp }->{ $gt };
	    $inv->{ used } = 1;
	    $vref->{ pstd } = 1;
	    $vref->{ qstd } = 1;
	    $hash->{ $softp }->{ SV }->{ type } = "INV";
	    $hash->{ $softp }->{ SV }->{ splits } += $scv ? $scv->{ cnt } : 0;
	    $hash->{ $softp }->{ SV }->{ pairs } += $cnt;
	    $hash->{ $softp }->{ SV }->{ clusters }++;
	    my $ref = $dir == -1 ? ($hash->{ $softp }->{ $REF->{ $softp } } ? $hash->{ $softp }->{ $REF->{ $softp } } : undef) : undef;
	    adjCnt($vref, $scv, $ref);
	    #adjCnt($vref, $inv);
	    my %dels5 = ();
	    $dels5{ $softp }->{ $gt } = $cnt;
	    $cov->{ $softp } = $cov->{ $softp-1 } ? $cov->{ $softp-1 } : $cnt;
	    $scv->{ used } = 1;
	    realigndel($hash, \%dels5, $cov, $sclip5, $sclip3, $REF, $chr, $bams);
	    print STDERR "  Found INV SV: $seq $softp $gt BP: $bp cov: $cov->{ $softp } Cnt: $cnt EXTRA: $EXTRA $ms $me $mlen cnt: $scv->{cnt} ", ($bp - $softp)/abs($mlen), "\t DIR: $dir Side: $side\n" if ( $opt_y );
	    return $hash->{ $softp }->{ $gt };
	}
    }
    return undef;
}

sub adddisccnt {
    my $svref = shift;
    $svref->{ disc }++; # Discordant count
}

sub filterSV {
    my $sva = shift;
    foreach my $sv (@$sva) {
	my ($ms, $me, $cnt, $mlen, $s, $e, $pmean, $qmean, $Qmean, $nm) = checkCluster( $sv->{ mates } );
	if ( $ms ) {
	    $sv->{ mstart } = $ms;
	    $sv->{ mend } = $me;
	    $sv->{ cnt } = $cnt;
	    $sv->{ mlen } = $mlen;
	    $sv->{ start } = $s;
	    $sv->{ end } = $e;
	    $sv->{ pmean } = $pmean;
	    $sv->{ qmean } = $qmean;
	    $sv->{ Qmean } = $Qmean;
	    $sv->{ nm } = $nm;
	} else {
	    $sv->{ used } = 1;
	}
	if ( $sv->{ disc } && $sv->{ cnt }/$sv->{ disc } < 0.5 ) { # Too many unhappy mates are false positives
	    $sv->{ used } = 1 unless( $sv->{ cnt }/$sv->{ disc } >= 0.35 && $sv->{ cnt } >= 5 ); 
	}
	my @soft = ();
	while( my ($sp, $v) = each %{ $sv->{ soft } } ) {
	    push(@soft, [$sp, $v]);
	}
	@soft = sort { $b->[1] <=> $a->[1]; } @soft;
	$sv->{ softp } = @soft > 0 ? $soft[0]->[0] : 0;
	push(@{ $SOFTP2SV{ $sv->{ softp } } }, $sv) if ( $sv->{ softp } );
	print STDERR "SV cluster: $s $e $ms $me Cnt: $sv->{ cnt } Discordant Cnt: $sv->{ disc } Softp: $sv->{ softp } Used: ", $sv->{ used } ? 1 : 0, "\n" if ( $opt_y );
	#print STDERR "SV cluster: $s $e $ms $me Cnt: $sv->{ cnt } Discordant Cnt: $sv->{ disc } Softp: $sv->{ softp }\n";# if ( $opt_y );
    }
}

sub checkCluster {
    my $mates = shift;
    my @mates = sort { $a->[0] <=> $b->[0] } @$mates;
    my @cluster = ({ cnt => 0, ms => $mates[0]->[0], me => $mates[0]->[1], s => $mates[0]->[3], e => $mates[0]->[4] });
    my $cur = 0;
    foreach my $m (@mates) {
	if ( $m->[0] - $cluster[$cur]->{ me } > $MINSVCDIST*$RLEN ) {
	    $cur++; 
	    $cluster[$cur] = { cnt => 0, ms => $m->[0], me => $m->[1], s => $m->[3], e => $m->[4] };
	}
	$cluster[$cur]->{ cnt }++;
	$cluster[$cur]->{ mlen } += $m->[2];
	$cluster[$cur]->{ me } = $m->[1] if ( $m->[1] > $cluster[$cur]->{ me } );
	$cluster[$cur]->{ s } = $m->[3] if ( $m->[3] < $cluster[$cur]->{ s } );
	$cluster[$cur]->{ e } = $m->[4] if ( $m->[4] > $cluster[$cur]->{ e } );
	$cluster[$cur]->{ pmean } += $m->[5];
	$cluster[$cur]->{ qmean } += $m->[6];
	$cluster[$cur]->{ Qmean } += $m->[7];
	$cluster[$cur]->{ nm } += $m->[8];
    }
    @cluster = sort { $b->{ cnt } <=> $a->{ cnt } } @cluster;
    print STDERR join("; ", "Clusters", (map {($_->{ cnt }, $_->{s}, $_->{e}, $_->{ms}, $_->{me});} @cluster), "out of", @mates+0), "\n" if ( $opt_y );
    return $cluster[0]->{ cnt }/(@mates+0) >= 0.60 ? ($cluster[0]->{ ms }, $cluster[0]->{ me }, $cluster[0]->{ cnt }, $cluster[0]->{ mlen }/$cluster[0]->{ cnt}, $cluster[0]->{ s }, $cluster[0]->{ e }, $cluster[0]->{ pmean }, $cluster[0]->{ qmean }, $cluster[0]->{ Qmean }, $cluster[0]->{ nm }) : (0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

# Mark SV clusters as used
sub markSV {
    my ($s, $e, $sv) = @_;
    my ($cov, $pairs, $cnt) = (0, 0, 0);
    foreach my $sr (@$sv) {
	foreach my $r (@$sr) {
	    my ($rs, $re) = $r->{ start } < $r->{ mstart } ? ($r->{ end }, $r->{ mstart }) : ($r->{ mend }, $r->{ start });
	    print STDERR "   Marking SV $s $e $rs $re cnt: $r->{ cnt }\n" if ( $opt_y );
	    if ( isOverlap($s, $e, $rs, $re) ) {
		    print STDERR "       SV $s $e $rs $re cnt: $r->{ cnt } marked\n" if ( $opt_y );
		    $r->{ used } = 1;
		    $cnt++;
		    $pairs += $r->{ cnt };
		    $cov += int(($r->{ cnt } * $RLEN)/($r->{ end } - $r->{ start })) + 1;
	    }
	}
    }
    return ($cov, $cnt, $pairs);
}

# Mark DUP clusters as used
sub markDUPSV {
    my ($s, $e, $sv) = @_;
    my ($cov, $pairs, $cnt) = (0, 0, 0);
    foreach my $sr (@$sv) {
	foreach my $r (@$sr) {
	    my ($rs, $re) = $r->{ start } < $r->{ mstart } ? ($r->{ start }, $r->{ mend }) : ($r->{ mstart}, $r->{ end });
	    print STDERR "   Marking DUP SV $s $e $rs $re cnt: $r->{ cnt }\n" if ( $opt_y );
	    if ( isOverlap($s, $e, $rs, $re) ) {
		print STDERR "       DUP SV $s $e $rs $re cnt: $r->{ cnt } marked\n" if ( $opt_y );
		$r->{ used } = 1;
		$cnt++;
		$pairs += $r->{ cnt };
		$cov += int(($r->{ cnt } * $RLEN)/($r->{ end } - $r->{ start })) + 1;
	    }
	}
    }
    return ($cnt, $pairs);
}

sub findDEL {
    my ($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $svfdel, $svrdel) = @_;
    foreach my $del (@$svfdel) {
	next if ( $del->{ used } );
	my ($ms, $me, $cnt, $mlen, $start, $end, $pmean, $qmean, $Qmean, $nm) = map { $del->{ $_ }; } qw(mstart mend cnt mlen start end pmean qmean Qmean nm);
	next unless( $cnt >= $MINR );
	#use Object; print STDERR "Pass: ", Object::Perl($inv);
	my @soft = ();
	while( my ($sp, $v) = each %{ $del->{ soft } } ) {
	    push(@soft, [$sp, $v]);
	}
	@soft = sort { $b->[1] <=> $a->[1]; } @soft;
	my $softp = @soft > 0 ? $soft[0]->[0] : 0; 
	print STDERR "\n\nWorking DEL 5' $softp mate cluster cnt: $cnt\n" if ( $opt_y );
	if ( $softp ) {
	    print STDERR "\n\nWorking DEL 5' $softp mate cluster cnt: $cnt\n" if ( $opt_y );
	    next unless( $sclip3->{ $softp } );
	    my $scv = $sclip3->{ $softp };
	    next if ($scv->{ used });
	    my $seq = findconseq($scv);
	    next unless( $seq && length($seq) >= $SEED2 );
	    unless( isLoaded($chr, $ms, $me, $REF) ) {
		getREF($chr, $ms, $me, $REF, 500);
		parseSAM($chr, $ms-200, $me+200, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1);
	    }
	    my ($bp, $EXTRA) = findMatch($seq, $REF, $softp, 1);
	    next unless( $bp );
	    next unless( $bp - $softp > 30 && isOverlap($softp, $bp, $end, $ms) );
	    $bp--;
	    my $dellen = $bp - $softp + 1;
	    $bp-- && $softp-- while( $REF->{ $bp } eq $REF->{ $softp - 1 } );
	    $hash->{ $softp }->{ "-$dellen" }->{ cnt } = 0;
	    $hash->{ $softp }->{ SV }->{ type } = "DEL";
	    $hash->{ $softp }->{ SV }->{ pairs } += $cnt;
	    $hash->{ $softp }->{ SV }->{ splits } += $scv->{ cnt };
	    $hash->{ $softp }->{ SV }->{ clusters }++;
	    my $ref = $hash->{ $softp }->{ "-$dellen" };
	    $cov->{ $softp } = $cnt unless( $cov->{ $softp } && $cov->{ $softp } > $cnt );
	    $cov->{ $softp } = $cov->{ $bp } if ( $cov->{ $bp } && $cov->{ $softp } < $cov->{ $bp } );
	    adjCnt($ref, $scv, $hash->{ $softp }->{ $REF->{ $softp } });
	    #my $intcov = 0;
	    #map { $intcov += $cov->{ $_ } ? $cov->{ $_ } : 0; } ($start .. $end);
	    #my $mcnt = int($cov->{ $softp  - 1 } * $cnt * 2 * $RLEN/$intcov) + 1;
	    my $mcnt = $cnt;
	    adjCnt($ref, {cnt => $mcnt, hicnt => $mcnt, 1 => int($mcnt/2), -1 => $mcnt - int($mcnt/2), qmean => $qmean*$mcnt/$cnt, pmean => $pmean*$mcnt/$cnt, Qmean => $Qmean*$mcnt/$cnt, nm => $nm*$mcnt/$cnt}); 
	    $del->{ used } = 1;
	    my ($svcov, $clusters, $pairs) = markSV($softp, $bp, [$svrdel]);
	    print STDERR "    Found DEL SV from 5' softclip unhappy reads: $bp -$dellen Cnt: $cnt AdjCnt: $mcnt\n" if ( $opt_y );
	} else { # Look within a read length
	    getREF($chr, $ms, $me, $REF, 500) unless( isLoaded($chr, $ms, $me, $REF) );
	    print STDERR "\n\nWorking DEL 5' no softp mate cluster cnt: $cnt\n" if ( $opt_y );
	    #for(my $i = $end + 1; $i < $end + $RLEN; $i++) 
	    while( my ($i, $scv) = each %$sclip3 ) {
		#next unless( $sclip3->{ $i } );
		#my $scv = $sclip3->{ $i };
		next if ( $scv->{ used } );
		next unless( $i >= $end - 3 && $i - $end < 3 * $RLEN );
		my $seq = findconseq($scv);
		next unless( $seq && length($seq) >= $SEED2 );
		$softp = $i;
		my ($bp, $EXTRA) = findMatch($seq, $REF, $softp, 1);
		($bp, $EXTRA) = findMatch($seq, $REF, $softp, 1, $SEED2, 0) unless( $bp );
		next unless( $bp );
		next unless( $bp - $softp > 30 && isOverlap($softp, $bp, $end, $ms) );
		$bp--;
		my $dellen = $bp - $softp + 1;
		$hash->{ $softp }->{ "-$dellen" }->{ cnt } = 0;
		$hash->{ $softp }->{ SV }->{ type } = "DEL";
		$hash->{ $softp }->{ SV }->{ pairs } += $cnt;
		$hash->{ $softp }->{ SV }->{ splits } += $scv->{ cnt };
		$hash->{ $softp }->{ SV }->{ clusters }++;
		my $ref = $hash->{ $softp }->{ "-$dellen" };
		$cov->{ $softp } = $cnt unless( $cov->{ $softp } && $cov->{ $softp } > $cnt );
		$cov->{ $softp } = $cov->{ $bp } if ( $cov->{ $bp } && $cov->{ $softp } < $cov->{ $bp } );
		adjCnt($ref, $scv);
		#my $intcov = 0;
		#map { $intcov += $cov->{ $_ } ? $cov->{ $_ } : 0; } ($start .. $end);
		#my $mcnt = int(($cov->{ $softp  - 1 } ? $cov->{ $softp  - 1 } : $cov->{ $softp }) * $cnt * 2 * $RLEN/$intcov) + 1;
		my $mcnt = $cnt;
		adjCnt($ref, {cnt => $mcnt, hicnt => $mcnt, 1 => int($mcnt/2), -1 => $mcnt - int($mcnt/2), qmean => $qmean*$mcnt/$cnt, pmean => $pmean*$mcnt/$cnt, Qmean => $Qmean*$mcnt/$cnt, nm => $nm*$mcnt/$cnt}); 
		$del->{ used } = 1;
		my ($svcov, $clusters, $pairs) = markSV($softp, $bp, [$svrdel]);
		print STDERR "    Found DEL SV from 5' softclip happy reads: $bp -$dellen Cnt: $cnt AdjCnt: $mcnt\n" if ( $opt_y );
		last;
	    }
	}
    }
    foreach my $del (@$svrdel) {
	next if ( $del->{ used } );
	my ($ms, $me, $cnt, $mlen, $start, $end, $pmean, $qmean, $Qmean, $nm) = map { $del->{ $_ }; } qw(mstart mend cnt mlen start end pmean qmean Qmean nm);
	next unless( $cnt >= $MINR );
	#use Object; print STDERR "Pass: ", Object::Perl($del);
	my @soft = ();
	while( my ($sp, $v) = each %{ $del->{ soft } } ) {
	    push(@soft, [$sp, $v]);
	}
	@soft = sort { $b->[1] <=> $a->[1]; } @soft;
	my $softp = @soft > 0 ? $soft[0]->[0] : 0;
	if ( $softp ) {
	    print STDERR "\n\nWorking DEL 3' $softp mate cluster cnt: $cnt\n" if ( $opt_y );
	    next unless( $sclip5->{ $softp } );
	    my $scv = $sclip5->{ $softp };
	    next if ($scv->{ used });
	    my $seq = findconseq($scv);
	    next unless( $seq && length($seq) >= $SEED2 );
	    unless( isLoaded($chr, $ms, $me, $REF) ) {
		getREF($chr, $ms, $me, $REF, 500);
		parseSAM($chr, $ms-200, $me+200, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1);
	    }
	    my ($bp, $EXTRA) = findMatch($seq, $REF, $softp, -1);
	    ($bp, $EXTRA) = findMatch($seq, $REF, $softp, -1, $SEED2, 0) unless( $bp );
	    next unless( $bp );
	    next unless( $softp - $bp > 30 && isOverlap($bp, $softp, $me, $start) );
	    $bp++;
	    $softp--;
	    my $dellen = $softp - $bp + 1;
	    $hash->{ $bp }->{ "-$dellen" }->{ cnt } = 0;
	    $hash->{ $bp }->{ SV }->{ type } = "DEL";
	    $hash->{ $bp }->{ SV }->{ pairs } += $cnt;
	    $hash->{ $bp }->{ SV }->{ splits } += $scv->{ cnt };
	    $hash->{ $bp }->{ SV }->{ clusters }++;
	    my $ref = $hash->{ $bp }->{ "-$dellen" };
	    adjCnt($ref, $scv);
	    $cov->{ $bp } = $cnt unless( $cov->{ $bp } && $cov->{ $bp } > $cnt );
	    $cov->{ $bp } = $cov->{ $softp } if ( $cov->{ $softp } && $cov->{ $softp } > $cov->{ $bp } );
	    #my $intcov = 0;
	    #map { $intcov += $cov->{ $_ } ? $cov->{ $_ } : 0; } ($start .. $end);
	    #my $mcnt = int(($cov->{ $bp  - 1 } ? $cov->{ $bp  - 1 } : $cov->{ $bp }) * $cnt * 2 * $RLEN/$intcov) + 1;
	    my $mcnt = $cnt;
	    adjCnt($ref, {cnt => $mcnt, hicnt => $mcnt, 1 => int($mcnt/2), -1 => $mcnt - int($mcnt/2), qmean => $qmean*$mcnt/$cnt, pmean => $pmean*$mcnt/$cnt, Qmean => $Qmean*$mcnt/$cnt, nm => $nm*$mcnt/$cnt}); 
	    $del->{ used } = 1;
	    my ($svcov, $clusters, $pairs) = markSV($bp, $softp, [$svfdel]);
	    print STDERR "    Found DEL SV from 3' softclip unhappy reads: $bp -$dellen Cnt: $cnt AdjCnt: $mcnt\n" if ( $opt_y );
	} else {
	    print STDERR "\n\nWorking DEL 3' no softp mate cluster $chr $ms $me cnt: $cnt\n" if ( $opt_y );
	    getREF($chr, $ms, $me, $REF, 500) unless( isLoaded($chr, $ms, $me, $REF) );
	    #for(my $i = $start - 3*$RLEN; $i < $start; $i++)
	    while( my ($i, $scv) = each %$sclip5 ) {
		#next unless( $sclip5->{ $i } );
		#my $scv = $sclip5->{ $i };
		next if ( $scv->{ used } );
		next unless( $i <= $start + 3 && $start - $i < 3 * $RLEN );
		my $seq = findconseq($scv);
		next unless( $seq && length($seq) >= $SEED2 );
		$softp = $i;
		my ($bp, $EXTRA) = findMatch($seq, $REF, $softp, -1);
		($bp, $EXTRA) = findMatch($seq, $REF, $softp, -1, $SEED2, 0) unless( $bp );
		next unless( $bp );
		next unless( $softp - $bp > 30 && isOverlap($bp, $softp, $me, $start) );
		$bp++;
		$softp--;
		my $dellen = $softp - $bp + 1;
		$hash->{ $bp }->{ "-$dellen" }->{ cnt } = 0;
		$hash->{ $bp }->{ SV }->{ type } = "DEL";
		$hash->{ $bp }->{ SV }->{ pairs } += $cnt;
		$hash->{ $bp }->{ SV }->{ splits } += $scv->{ cnt };
		$hash->{ $bp }->{ SV }->{ clusters }++;
		my $ref = $hash->{ $bp }->{ "-$dellen" };
		adjCnt($ref, $scv);
		$cov->{ $bp } = $cnt unless( $cov->{ $bp } );
		$cov->{ $bp } = $cov->{ $softp } if ( $cov->{ $softp } && $cov->{ $softp } > $cov->{ $bp } );
		$cov->{ $bp } += $scv->{ cnt };
		#my $intcov = 0;
		#map { $intcov += $cov->{ $_ } ? $cov->{ $_ } : 0; } ($start .. $end);
		#my $mcnt = int(($cov->{ $bp  - 1 } ? $cov->{ $bp  - 1 } : $cov->{ $bp }) * $cnt * 2 * $RLEN/$intcov) + 1;
		my $mcnt = $cnt;
		adjCnt($ref, {cnt => $mcnt, hicnt => $mcnt, 1 => int($mcnt/2), -1 => $mcnt - int($mcnt/2), qmean => $qmean*$mcnt/$cnt, pmean => $pmean*$mcnt/$cnt, Qmean => $Qmean*$mcnt/$cnt, nm => $nm*$mcnt/$cnt}); 
		$del->{ used } = 1;
		my ($svcov, $clusters, $pairs) = markSV($bp, $softp, [$svfdel]);
		print STDERR "    Found DEL SV from 3' softclip happy reads: $bp -$dellen Cnt: $cnt AdjCnt: $mcnt\n" if ( $opt_y );
		last;
	    }
	}
    }
}

# With discordant only
# (svfdel) --> | ............ | <-- (svrdel)
sub findDELdisc { 
    my ($hash, $cov, $REF, $chr, $bams, $svfdel, $svrdel, $sclip5, $sclip3) = @_;
    my $MINDIST = 8*$RLEN; # the minimum distance between two clusters
    foreach my $del (@$svfdel) {
	next if ( $del->{ used } );
	my ($ms, $me, $cnt, $mlen, $start, $end, $pmean, $qmean, $Qmean, $nm) = map { $del->{ $_ }; } qw(mstart mend cnt mlen start end pmean qmean Qmean nm);
	next if ( %SPLICE && abs($mlen) < 250000 ); # more stringent for RNA-Seq
	next unless( $cnt >= $MINR + 5 );
	next unless( $ms > $end + $MINDIST );
	next unless( $Qmean/$cnt > $DISCPAIRQUAL );
	$mlen = $ms - $end - int($RLEN/($cnt+1));
	next unless( $mlen > 0 && $mlen > $MINDIST );
	my $bp = $end + int(($RLEN/($cnt+1))/2);
	$bp = $del->{ softp } if ( $del->{ softp } );
	$hash->{ $bp }->{ "-$mlen" }->{ cnt } = 0;
	$hash->{ $bp }->{ SV }->{ type } = "DEL";
	$hash->{ $bp }->{ SV }->{ splits } += $sclip3->{ $end+1 } ? $sclip3->{ $end+1 }->{ cnt } : 0;
	$hash->{ $bp }->{ SV }->{ splits } += $sclip5->{ $ms } ? $sclip5->{ $ms }->{ cnt } : 0;
	$hash->{ $bp }->{ SV }->{ pairs } += $cnt;
	$hash->{ $bp }->{ SV }->{ clusters }++;
	my $ref = $hash->{ $bp }->{ "-$mlen" };
	print STDERR "  Found DEL with discordant pairs only: cnt: $cnt BP: $bp Len: $mlen $start-$end<->$ms-$me\n" if ( $opt_y );
	adjCnt($ref, {cnt => 2*$cnt, hicnt => 2*$cnt, "1" => $cnt, "-1" => $cnt, qmean => $qmean*2, pmean => $pmean*2, Qmean => $Qmean*2, nm => $nm*2}); 
	$cov->{ $bp } = 2*$cnt unless( $cov->{ $bp } );
	$del->{ used } = 1;
	my ($svcov, $clusters, $pairs) = markSV( $end, $ms, [$svrdel] );
    }
    foreach my $del (@$svrdel) {
	next if ( $del->{ used } );
	my ($ms, $me, $cnt, $mlen, $start, $end, $pmean, $qmean, $Qmean, $nm) = map { $del->{ $_ }; } qw(mstart mend cnt mlen start end pmean qmean Qmean nm);
	next if ( %SPLICE && abs($mlen) < 250000 ); # more stringent for RNA-Seq
	next unless( $cnt >= $MINR + 5 );
	next unless( $start > $me + $MINDIST );
	next unless( $Qmean/$cnt > $DISCPAIRQUAL );
	$mlen = $start - $me - int($RLEN/($cnt+1));
	#use Object; print STDERR "3' $mlen ", Object::Perl($del), "\n";
	next unless( $mlen > 0 && $mlen > $MINDIST );
	my $bp = $me + int(($RLEN/($cnt+1))/2);
	$hash->{ $bp }->{ "-$mlen" }->{ cnt } = 0;
	$hash->{ $bp }->{ SV }->{ type } = "DEL";
	$hash->{ $bp }->{ SV }->{ splits } += $sclip3->{ $me+1 } ? $sclip3->{ $me+1 }->{ cnt } : 0;
	$hash->{ $bp }->{ SV }->{ splits } += $sclip5->{ $start } ? $sclip5->{ $start }->{ cnt } : 0;
	$hash->{ $bp }->{ SV }->{ pairs } += $cnt;
	$hash->{ $bp }->{ SV }->{ clusters }++;
	my $ref = $hash->{ $bp }->{ "-$mlen" };
	print STDERR "  Found DEL with discordant pairs only: cnt: $cnt BP: $bp Len: $mlen $start-$end<->$ms-$me\n" if ( $opt_y );
	$sclip5->{ $del->{ softp } }->{ used } = 1 if ( $del->{ softp } && $sclip5->{ $del->{ softp } } );
	adjCnt($ref, {cnt => 2*$cnt, hicnt => 2*$cnt, "1" => $cnt, "-1" => $cnt, qmean => $qmean*2, pmean => $pmean*2, Qmean => $Qmean*2, nm => $nm*2}); 
	$cov->{ $bp } = 2*$cnt unless( $cov->{ $bp } );
	$cov->{ $bp } = $cov->{ $start } if ( $cov->{ $start } && $cov->{ $bp } < $cov->{ $start } );
	$del->{ used } = 1;
	getREF($chr, $ms-100, $me+100, $REF, 200);
	my ($svcov, $clusters, $pairs) = markSV( $me, $start, [$svfdel] );
    }
}

# With discordant only
#  (svrdup) |<--.........-->| (svfdup)
sub findDUPdisc { 
    my ($hash, $cov, $REF, $chr, $bams, $svfdup, $svrdup, $sclip5, $sclip3) = @_;
    foreach my $dup (@$svfdup) {
	next if ( $dup->{ used } );
	my ($ms, $me, $cnt, $mlen, $start, $end, $pmean, $qmean, $Qmean, $nm) = map { $dup->{ $_ }; } qw(mstart mend cnt mlen start end pmean qmean Qmean nm);
	next unless( $cnt >= $MINR + 5 );
	next unless( $Qmean/$cnt > $DISCPAIRQUAL );
	$mlen = $end - $ms + int($RLEN/$cnt);
	my $bp = $ms - int(($RLEN/$cnt)/2);
	my $pe = $end; #$mlen + $bp - 1;
	unless( isLoaded($chr, $ms, $me, $REF) ) {
	    getREF($chr, $bp - 150, $bp + 150, $REF, 200) unless( $REF->{ $bp } );
	    parseSAM($chr, $ms-200, $me+200, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1);
	}
	my ($cntf, $cntr) = ($cnt, $cnt);
	my ($qmeanf, $qmeanr) = ($qmean, $qmean);
	my ($Qmeanf, $Qmeanr) = ($Qmean, $Qmean);
	my ($pmeanf, $pmeanr) = ($pmean, $pmean);
	my ($nmf, $nmr) = ($nm, $nm);
	if ( %{ $dup->{ soft } } ) {
	    my @soft;
	    while( my ($sp, $v) = each %{ $dup->{ soft } } ) {
		push(@soft, [$sp, $v]);
	    }
	    @soft = sort { $b->[1] <=> $a->[1]; } @soft;
	    $pe = $soft[0]->[0] if ( @soft > 0 );
	    next unless ( $sclip3->{ $pe } );
	    next if ( $sclip3->{ $pe }->{ used } );
	    $cntf = $sclip3->{ $pe }->{ cnt };
	    $qmeanf = $sclip3->{ $pe }->{ qmean };
	    $Qmeanf = $sclip3->{ $pe }->{ Qmean };
	    $pmeanf = $sclip3->{ $pe }->{ pmean };
	    $nmf = $sclip3->{ $pe }->{ nm };
	    my $seq = findconseq($sclip3->{ $pe });
	    my ($tbp, $EXTRA) = findMatch($seq, $REF, $bp, 1);
	    if ( $tbp && $tbp < $pe ) {
		$sclip3->{ $pe }->{ used } = 1;
		$tbp-- && $pe-- while( $REF->{ $tbp-1 } eq $REF->{ $pe - 1 } );
		$mlen = $pe - $tbp;
		$bp = $tbp;
		$pe--;
		$end = $pe;
		if ( $sclip5->{ $bp } ) {
		    $cntr = $sclip5->{ $bp }->{ cnt };
		    $qmeanr = $sclip5->{ $bp }->{ qmean };
		    $Qmeanr = $sclip5->{ $bp }->{ Qmean };
		    $pmeanr = $sclip5->{ $bp }->{ pmean };
		    $nmr = $sclip5->{ $bp }->{ nm };
		}
	    }
	}
	my $ins5 = join("", map { $REF->{ $_ } } ($bp .. ($bp+$SVFLANK-1)));
	my $ins3 = join("", map { $REF->{ $_ } } (($pe-$SVFLANK+1) .. $pe));
	my $ins = "$ins5<dup" . ($mlen - 2*$SVFLANK) . ">$ins3";
	$hash->{ $bp }->{ I }->{ "+$ins" }->{ cnt } = 0;
	$hash->{ $bp }->{ SV }->{ type } = "DUP";
	$hash->{ $bp }->{ SV }->{ pairs } += $cnt;
	$hash->{ $bp }->{ SV }->{ splits } += $dup->{ softp } && $sclip3->{ $dup->{ softp } } ? $sclip3->{ $dup->{ softp } }->{ cnt } : 0;
	$hash->{ $bp }->{ SV }->{ clusters }++;
	my $ref = $hash->{ $bp }->{ I }->{ "+$ins" };
	print STDERR "  Found DUP with discordant pairs only (forward): cnt: $cnt BP: $bp END: $pe $ins Len: $mlen $start-$end<->$ms-$me\n" if ( $opt_y );
	my $tcnt = $cntr + $cntf;
	#adjCnt($ref, {cnt => 2*$cnt, extracnt => 2*$cnt, hicnt => 2*$cnt, "1" => $cnt, "-1" => $cnt, qmean => $qmean*2, pmean => $pmean*2, Qmean => $Qmean*2, nm => $nm*2}); 
	adjCnt($ref, {cnt => $tcnt, extracnt => $tcnt, hicnt => $tcnt, "1" => $cntf, "-1" => $cntr, qmean => $qmeanf + $qmeanr, pmean => $pmeanf + $pmeanr, Qmean => $Qmeanf + $Qmeanr, nm => $nmf + $nmr}); 
	$dup->{ used } = 1;
	$cov->{ $bp } = $tcnt unless( $cov->{ $bp } );
	$cov->{ $bp } = $cov->{ $end } if ( $cov->{ $end } && $cov->{ $bp } < $cov->{ $end } );
	my ($clusters, $pairs) = markDUPSV( $bp, $pe, [$svrdup] );
	$hash->{ $bp }->{ SV }->{ clusters } += $clusters;
    }
    foreach my $dup (@$svrdup) {
	next if ( $dup->{ used } );
	my ($ms, $me, $cnt, $mlen, $start, $end, $pmean, $qmean, $Qmean, $nm) = map { $dup->{ $_ }; } qw(mstart mend cnt mlen start end pmean qmean Qmean nm);
	next unless( $cnt >= $MINR + 5 );
	next unless( $Qmean/$cnt > $DISCPAIRQUAL );
	$mlen = $me - $start + int($RLEN/$cnt);
	my $bp = $start - int(($RLEN/$cnt)/2);
	my $pe = $mlen + $bp - 1;
	my $tpe = $pe;
	unless( isLoaded($chr, $ms, $me, $REF) ) {
	    getREF($chr, $pe - 150, $pe + 150, $REF, 200) unless( $REF->{ $pe } );
	    parseSAM($chr, $ms-200, $me+200, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1);
	}
	my ($cntf, $cntr) = ($cnt, $cnt);
	my ($qmeanf, $qmeanr) = ($qmean, $qmean);
	my ($Qmeanf, $Qmeanr) = ($Qmean, $Qmean);
	my ($pmeanf, $pmeanr) = ($pmean, $pmean);
	my ($nmf, $nmr) = ($nm, $nm);
	if ( %{ $dup->{ soft } } ) {
	    my @soft;
	    while( my ($sp, $v) = each %{ $dup->{ soft } } ) {
		push(@soft, [$sp, $v]);
	    }
	    @soft = sort { $b->[1] <=> $a->[1]; } @soft;
	    $bp = $soft[0]->[0] if ( @soft > 0 );
	    next unless( $sclip5->{ $bp } );
	    next if ( $sclip5->{ $bp }->{ used } );
	    $cntr = $sclip5->{ $bp }->{ cnt };
	    $qmeanr = $sclip5->{ $bp }->{ qmean };
	    $Qmeanr = $sclip5->{ $bp }->{ Qmean };
	    $pmeanr = $sclip5->{ $bp }->{ pmean };
	    $nmr = $sclip5->{ $bp }->{ nm };
	    my $seq = findconseq($sclip5->{ $bp });
	    my ($tbp, $EXTRA) = findMatch($seq, $REF, $pe, -1);
	    if ( $tbp && $tbp > $bp ) {
		$sclip5->{ $bp }->{ used } = 1;
		$pe = $tbp;
		$mlen = $pe - $bp + 1;
		$tpe = $pe + 1;
		$tpe++ while( $REF->{ $bp + ($tpe-$pe-1) } eq $REF->{ $tpe } );
		if ( $sclip3->{ $tpe } ) {
		    $cntf = $sclip3->{ $tpe }->{ cnt };
		    $qmeanf = $sclip3->{ $tpe }->{ qmean };
		    $Qmeanf = $sclip3->{ $tpe }->{ Qmean };
		    $pmeanf = $sclip3->{ $tpe }->{ pmean };
		    $nmf = $sclip3->{ $tpe }->{ nm };
		}
	    }
	}
	my $ins5 = join("", map { $REF->{ $_ } } ($bp .. ($bp+$SVFLANK-1)));
	my $ins3 = join("", map { $REF->{ $_ } } (($pe-$SVFLANK+1) .. $pe));
	my $ins = "$ins5<dup" . ($mlen - 2*$SVFLANK) . ">$ins3";
	$hash->{ $bp }->{ I }->{ "+$ins" }->{ cnt } = 0;
	$hash->{ $bp }->{ SV }->{ type } = "DUP";
	$hash->{ $bp }->{ SV }->{ pairs } += $cnt;
	$hash->{ $bp }->{ SV }->{ splits } += $sclip5->{ $bp } ? $sclip5->{ $bp }->{ cnt } : 0;
	$hash->{ $bp }->{ SV }->{ splits } += $sclip3->{ $tpe } ? $sclip3->{ $tpe }->{ cnt } : 0;
	#$hash->{ $bp }->{ SV }->{ splits } += $sclip3->{ $me+1 } ? $sclip3->{ $me+1 }->{ cnt } : 0;
	$hash->{ $bp }->{ SV }->{ clusters }++;
	my $ref = $hash->{ $bp }->{ I }->{ "+$ins" };
	print STDERR "  Found DUP with discordant pairs only (reverse): cnt: $cnt BP: $bp Len: $mlen $start-$end<->$ms-$me\n" if ( $opt_y );
	#adjCnt($ref, {cnt => 2*$cnt, extracnt => 2*$cnt, hicnt => 2*$cnt, "1" => $cnt, "-1" => $cnt, qmean => $qmean*2, pmean => $pmean*2, Qmean => $Qmean*2, nm => $nm*2}); 
	my $tcnt = $cntr + $cntf;
	adjCnt($ref, {cnt => $tcnt, extracnt => $tcnt, hicnt => $tcnt, "1" => $cntf, "-1" => $cntr, qmean => $qmeanf + $qmeanr, pmean => $pmeanf + $pmeanr, Qmean => $Qmeanf + $Qmeanr, nm => $nmf + $nmr}); 
	$dup->{ used } = 1;
	$cov->{ $bp } = $tcnt unless( $cov->{ $bp } );
	$cov->{ $bp } = $cov->{ $me } if ( $cov->{ $me } && $cov->{ $bp } < $cov->{ $me } );
	my ($clusters, $pairs) = markDUPSV( $bp, $pe, [$svfdup] );
	$hash->{ $bp }->{ SV }->{ clusters } += $clusters;
    }
}

# Would only consider those with supports from both orientations
#  (svfinv5) --> | <-- (svrinv5) ..... (svfinv3) --> | <-- (svrinv3)
sub findINVdisc {
    my ($hash, $cov, $REF, $chr, $bams, $svfinv3, $svrinv3, $svfinv5, $svrinv5, $sclip5, $sclip3) = @_;
    foreach my $invf5 (@$svfinv5) {
	next if ( $invf5->{ used } );
	my ($ms, $me, $cnt, $mlen, $start, $end, $pmean, $qmean, $Qmean, $nm) = map { $invf5->{ $_ }; } qw(mstart mend cnt mlen start end pmean qmean Qmean nm);
	next unless( $Qmean/$cnt > $DISCPAIRQUAL );
	foreach my $invr5 (@$svrinv5) {
	    next if ( $invr5->{ used } );
	    my ($rms, $rme, $rcnt, $rmlen, $rstart, $rend, $rpmean, $rqmean, $rQmean, $rnm) = map { $invr5->{ $_ }; } qw(mstart mend cnt mlen start end pmean qmean Qmean nm);
	    next unless( $rQmean/$rcnt > $DISCPAIRQUAL );
	    next unless( $cnt + $rcnt > $MINR + 5 );
	    if (isOverlap($end, $me, $rstart, $rms)) {
		my $bp = int(abs($end+$rstart)/2);
		my $pe = int(abs($me+$rms)/2);
		getREF($chr, $pe - 150, $pe + 150, $REF, 100) unless( $REF->{ $pe } );
		my $len = $pe - $bp + 1;
		my $ins5 = reverse(join("", map { $REF->{ $_ } ? $REF->{ $_ } : "" } ($bp .. ($bp+$SVFLANK-1))));
		my $ins3 = reverse(join("", map { $REF->{ $_ } ? $REF->{ $_ } : "" } (($pe-$SVFLANK+1) .. $pe)));
		$ins3 =~ y/ATGC/TACG/;
		$ins5 =~ y/ATGC/TACG/;
		my $ins = "$ins3<inv" . ($len - 2*$SVFLANK) . ">$ins5";
		#print STDERR "R: $rms, $rme, $rcnt, $rmlen, $rstart, $rend, $rpmean, $rqmean, $rQmean, $rnm $bp $pe $len $ins\n";
		print STDERR "  Found INV with discordant pairs only: cnt: $cnt Len: $len $end-$rstart<->$me-$rms $ins\n" if ( $opt_y );
		$hash->{ $bp }->{ "-${len}^$ins" }->{ cnt } = 0 unless( $hash->{ $bp }->{ "-${len}^$ins" }->{ cnt } );
		my $vref = $hash->{ $bp }->{ "-${len}^$ins" };
		$invf5->{ used } = 1;
		$invr5->{ used } = 1;
		$vref->{ pstd } = 1;
		$vref->{ qstd } = 1;
		adjCnt($vref, {cnt => $cnt + $rcnt, hicnt => $cnt + $rcnt, "1" => $cnt, "-1" => $rcnt, qmean => $qmean + $rqmean, pmean => $pmean + $rpmean, Qmean => $Qmean + $rQmean, nm => $nm + $rnm});
		$hash->{ $bp }->{ SV }->{ type } = "INV";
		$hash->{ $bp }->{ SV }->{ pairs } += $cnt;
		$hash->{ $bp }->{ SV }->{ splits } += $sclip5->{ $start } ? $sclip5->{ $start }->{ cnt } : 0;
		$hash->{ $bp }->{ SV }->{ splits } += $sclip5->{ $ms } ? $sclip5->{ $ms }->{ cnt } : 0;
		$hash->{ $bp }->{ SV }->{ clusters }++;
		$cov->{ $bp } = 2*$cnt unless( $cov->{ $bp } );
		my ($svcov, $clusters, $pairs) = markSV( $bp, $pe, [$svfinv3, $svrinv3] );
	    }
	}
    }
    foreach my $invf3 (@$svfinv3) {
	next if ( $invf3->{ used } );
	my ($ms, $me, $cnt, $mlen, $start, $end, $pmean, $qmean, $Qmean, $nm) = map { $invf3->{ $_ }; } qw(mstart mend cnt mlen start end pmean qmean Qmean nm);
	foreach my $invr3 (@$svrinv3) {
	    next if ( $invr3->{ used } );
	    my ($rms, $rme, $rcnt, $rmlen, $rstart, $rend, $rpmean, $rqmean, $rQmean, $rnm) = map { $invr3->{ $_ }; } qw(mstart mend cnt mlen start end pmean qmean Qmean nm);
	    next unless( $rQmean/$rcnt > $DISCPAIRQUAL );
	    next unless( $cnt + $rcnt > $MINR + 5 );
	    if (isOverlap($me, $end, $rms, $rstart)) {
		my $pe = int(abs($end+$rstart)/2);
		my $bp = int(abs($me+$rms)/2);
		getREF($chr, $bp - 150, $bp + 150, $REF, 100) unless( $REF->{ $bp } );
		my $len = $pe - $bp + 1;
		my $ins5 = reverse(join("", map { $REF->{ $_ } ? $REF->{ $_ } : "" } ($bp .. ($bp+$SVFLANK-1))));
		my $ins3 = reverse(join("", map { $REF->{ $_ } ? $REF->{ $_ } : "" } (($pe-$SVFLANK+1) .. $pe)));
		$ins3 =~ y/ATGC/TACG/;
		$ins5 =~ y/ATGC/TACG/;
		my $ins = "$ins3<inv" . ($len - 2*$SVFLANK) . ">$ins5";
		print STDERR "  Found INV with discordant pairs only: cnt: $cnt Len: $len $me-$rms<->$end-$rstart $ins\n" if ( $opt_y );
		$hash->{ $bp }->{ "-${len}^$ins" }->{ cnt } = 0 unless( $hash->{ $bp }->{ "-${len}^$ins" }->{ cnt } );
		my $vref = $hash->{ $bp }->{ "-${len}^$ins" };
		$invf3->{ used } = 1;
		$invr3->{ used } = 1;
		$vref->{ pstd } = 1;
		$vref->{ qstd } = 1;
		adjCnt($vref, {cnt => $cnt + $rcnt, hicnt => $cnt + $rcnt, "1" => $cnt, "-1" => $rcnt, qmean => $qmean + $rqmean, pmean => $pmean + $rpmean, Qmean => $Qmean + $rQmean, nm => $nm + $rnm});
		$hash->{ $bp }->{ SV }->{ type } = "INV";
		$hash->{ $bp }->{ SV }->{ pairs } += $cnt;
		$hash->{ $bp }->{ SV }->{ splits } += $sclip3->{ $end+1 } ? $sclip3->{ $end+1 }->{ cnt } : 0;
		$hash->{ $bp }->{ SV }->{ splits } += $sclip3->{ $me+1 } ? $sclip3->{ $me+1 }->{ cnt } : 0;
		$hash->{ $bp }->{ SV }->{ clusters }++;
		$cov->{ $bp } = 2*$cnt unless( $cov->{ $bp } );
		my ($svcov, $clusters, $pairs) = markSV( $bp, $pe, [$svfinv5, $svrinv5] );
	    }
	}
    }
}

sub temp {
    my ($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $svfinv3, $svrinv3, $svfinv5, $svrinv5) = @_;
    foreach my $invf3 (@$svfinv3) {
	next if ( $invf3->{ used } );
	my ($ms, $me, $cnt, $mlen) = checkCluster($invf3->{ mates });
	next unless( $ms );
	next unless( $cnt >= $MINR );
	my @soft = ();
	while( my ($sp, $v) = each %{ $invf3->{ soft } } ) {
	    push(@soft, [$sp, $v]);
	}
	@soft = sort { $b->[1] <=> $a->[1]; } @soft;
	my $softp = @soft > 0 ? $soft[0]->[0] : 0;
	if ( $softp ) {
	    my $sc3v = $sclip3->{ $softp };
	    next if ($sc3v->{ used });
	    my $seq = findconseq($sc3v);
	    next unless( $seq );
	    my ($bp) = findMatchRev($seq, $REF, $softp, 1);
	    next unless( $bp );
	    $softp-- && $bp++ while( $REF->{ $softp - 1 } eq $REVCOMP{ $REF->{ $bp + 1 } } );
	    print STDERR "  Found INV SV 3': $seq $softp $bp $ms $me $mlen cnt: $sc3v->{cnt}\n" if ( $opt_y );
	    if ( $bp > $softp && ($bp - $softp)/$mlen < 1.5 ) {
		my $len = $bp - $softp + 1;
		my $ins = reverse(join("", map { $REF->{ $_ } } ($softp .. $bp)));
		$ins =~ y/ATGC/TACG/;
		$hash->{ $softp }->{ "-${len}^$ins" }->{ cnt } = 0;
		my $vref = $hash->{ $softp }->{ "-${len}^$ins" };
		$invf3->{ used } = 1;
		$vref->{ pstd } = 1;
		$vref->{ qstd } = 1;
		adjCnt($vref, $sc3v, $hash->{ $softp }->{ $REF->{ $softp } });
		my %dels5 = ();
		$dels5{ $softp }->{ "-${len}^$ins" } = $cnt;
		realigndel($hash, \%dels5, $cov, $sclip5, $sclip3, $REF, $chr, $bams);
	    }
	    $sc3v->{ used } = 1;
	} else { # Look within 100bp
	}
    }
}

sub findsv {
    my ($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $svfdel, $svrdel, $svfdup, $svrdup, $svfinv5, $svrinv5, $svfinv3, $svrinv3) = @_;
    my @tmp5;
    while(my($p, $sc5v) = each %$sclip5) {
	next if ( $sc5v->{ used } );
	push(@tmp5, [$p, $sc5v, $sc5v->{cnt}]);
    }
    @tmp5 = sort {$b->[2] <=> $a->[2];} @tmp5;
    my @tmp3;
    while(my($p, $sc3v) = each %$sclip3) {
	next if ( $sc3v->{ used } );
	push(@tmp3, [$p, $sc3v, $sc3v->{cnt}]);
    }
    @tmp3 = sort {$b->[2] <=> $a->[2];} @tmp3;
    for(my $i = 0; $i < @tmp5; $i++) {
	my ($p5, $sc5v, $cnt5) = @{ $tmp5[$i] };
	last if ($cnt5 < $MINR);
	next if ($sc5v->{ used });
	next if ($SOFTP2SV{ $p5 } && $SOFTP2SV{ $p5 }->[0]->{ used });
	my $seq = findconseq($sc5v);
	next unless( $seq && length($seq) >= $SEED2 );
	print STDERR "  Finding SV 5': $seq $p5 cnt: $cnt5\n" if ( $opt_y );
	my ($bp, $EXTRA) = findMatch($seq, $REF, $p5, -1);
	if ( $bp ) {
	    if ( $bp < $p5 ) { # candidate deletion
		my ($pairs, $pmean, $qmean, $Qmean, $nm) = checkPairs($chr, $bp, $p5, [$svfdel, $svrdel]);
		next unless( $pairs );
		$p5--; $bp++;
		my $dellen = $p5 - $bp + 1;
		$hash->{ $bp }->{ "-$dellen" }->{ cnt } = 0;
		$hash->{ $bp }->{ SV }->{ type } = "DEL";
		$hash->{ $bp }->{ SV }->{ pairs } += $pairs;
		$hash->{ $bp }->{ SV }->{ splits } += $cnt5;
		$hash->{ $bp }->{ SV }->{ clusters } += $pairs ? 1 : 0;
		my $ref = $hash->{ $bp }->{ "-$dellen" };
		$cov->{ $bp } = $pairs + $sc5v->{ cnt } unless( $cov->{ $bp } );
		$cov->{ $bp } = $cov->{ $p5+1 } if ( $cov->{ $bp } < $cov->{ $p5+1 } );
		adjCnt($ref, $sc5v);
		adjCnt($ref, {cnt => $pairs, hicnt => $pairs, 1 => int($pairs/2), -1 => $pairs - int($pairs/2), pmean => $pmean, qmean => $qmean, Qmean => $Qmean, nm => $nm});
	    } else { # candidate duplication
	    }
	} else { # candidate inversion
	    ($bp, $EXTRA) = findMatchRev($seq, $REF, $p5, -1);
	    #($bp, $EXTRA) = findMatchRev($seq, $REF, $p5, -1, $SEED2) unless($bp);
	    next unless( $bp );
	    next unless( abs($bp-$p5) > $SVFLANK );
	    if ( $bp > $p5 ) { # bp at 3' side
	    } else { # bp at 5' side
		($bp, $p5) = ($p5, $bp);
	    }
	    $bp--;
	    $p5-- && $bp++ while( $REF->{ $p5-1 } && $REF->{ $p5-1 } eq $REVCOMP{ $REF->{ $bp+1 } } );
	    my $ins5 = reverse(join("", map { $REF->{ $_ } } (($bp-$SVFLANK+1) .. $bp)));
	    my $ins3 = reverse(join("", map { $REF->{ $_ } ? $REF->{ $_ } : "" } ($p5 .. ($p5+$SVFLANK-1))));
	    $ins3 =~ y/ATGC/TACG/;
	    $ins5 =~ y/ATGC/TACG/;
	    my $mid = $bp - $p5 - length($ins5) - length($ins3) + 1;
	    my $vn = "-" . ($bp - $p5 + 1) . "^$ins5<inv$mid>$ins3$EXTRA";
	    $hash->{ $p5 }->{ $vn }->{ cnt } = 0 unless( $hash->{ $p5 }->{ $vn }->{ cnt } );
	    $hash->{ $p5 }->{ SV }->{ type } = "INV";
	    $hash->{ $p5 }->{ SV }->{ pairs } += 0;
	    $hash->{ $p5 }->{ SV }->{ splits } += $cnt5;
	    $hash->{ $p5 }->{ SV }->{ clusters } += 0;
	    my $ref = $hash->{ $p5 }->{ $vn };
	    adjCnt($ref, $sc5v);
	    $cov->{ $p5 } += $cnt5; # unless( $cov->{ $p5 } && $cov->{ $p5 } > $cnt5 );
	    $cov->{ $p5 } = $cov->{ $bp } if ( $cov->{ $bp } && $cov->{ $p5 } < $cov->{ $bp } );
	    print STDERR "    Found INV: $p5 $vn Cnt: $cnt5\n" if ( $opt_y );
	}
    }
    for(my $i = 0; $i < @tmp3; $i++) {
	my ($p3, $sc3v, $cnt3) = @{ $tmp3[$i] };
	last if ($cnt3 < $MINR);
	next if ($sc3v->{ used });
	next if ($SOFTP2SV{ $p3 } && $SOFTP2SV{ $p3 }->[0]->{ used });
	my $seq = findconseq($sc3v);
	next unless( $seq && length($seq) >= $SEED2 );
	print STDERR "  Finding SV 3': $seq $p3 cnt: $cnt3\n" if ( $opt_y );
	my ($bp, $EXTRA) = findMatch($seq, $REF, $p3, 1);
	if ( $bp ) {
	    if ( $bp > $p3 ) { # candidate deletion
		my ($pairs, $pmean, $qmean, $Qmean, $nm) = checkPairs($chr, $p3, $bp, [$svfdel, $svrdel]);
		next unless( $pairs );
		my $dellen = $bp - $p3;
		$bp--;
		$bp-- && $p3-- while( $REF->{ $p3 - 1 } eq $REF->{ $bp } );
		$hash->{ $p3 }->{ "-$dellen" }->{ cnt } = 0;
		$hash->{ $p3 }->{ SV }->{ type  } = "DEL";
		$hash->{ $p3 }->{ SV }->{ pairs } += $pairs;
		$hash->{ $p3 }->{ SV }->{ splits } += $cnt3;
		$hash->{ $p3 }->{ SV }->{ clusters } += $pairs ? 1 : 0;
		my $ref = $hash->{ $p3 }->{ "-$dellen" };
		adjCnt($ref, $sc3v);
		adjCnt($ref, {cnt => $pairs, hicnt => $pairs, 1 => int($pairs/2), -1 => $pairs - int($pairs/2), pmean => $pmean, qmean => $qmean, Qmean => $Qmean, nm => $nm});
	    } else { # candidate duplication
	    }
	} else { # candidate inversion
	    ($bp, $EXTRA) = findMatchRev($seq, $REF, $p3, 1);
	    #($bp, $EXTRA) = findMatchRev($seq, $REF, $p3, 1, $SEED2) unless($bp);
	    next unless( $bp );
	    next unless( abs($bp-$p3) > $SVFLANK );
	    if ( $bp < $p3 ) { # bp at 5' side 
		($bp, $p3) = ($p3, $bp);
		$p3++; $bp--;
	    } else { # bp at 3' side
	    }
	    $p3-- && $bp++ while( $REF->{ $p3-1 } && $REF->{ $p3-1 } eq $REVCOMP{ $REF->{ $bp+1 } } );
	    my $ins5 = reverse(join("", map { $REF->{ $_ } ? $REF->{ $_ } : "" } (($bp-$SVFLANK+1) .. $bp)));
	    my $ins3 = reverse(join("", map { $REF->{ $_ } ? $REF->{ $_ } : "" } ($p3 .. ($p3+$SVFLANK-1))));
	    $ins3 =~ y/ATGC/TACG/;
	    $ins5 =~ y/ATGC/TACG/;
	    my $mid = $bp - $p3 - 2*$SVFLANK + 1;
	    my $vn = "-" . ($bp - $p3 + 1) . "^$EXTRA$ins5<inv$mid>$ins3";
	    $hash->{ $p3 }->{ $vn }->{ cnt } = 0 unless( $hash->{ $p3 }->{ $vn }->{ cnt } );
	    $hash->{ $p3 }->{ SV }->{ type } = "INV";
	    $hash->{ $p3 }->{ SV }->{ pairs } += 0;
	    $hash->{ $p3 }->{ SV }->{ splits } += $cnt3;
	    $hash->{ $p3 }->{ SV }->{ clusters } += 0;
	    my $ref = $hash->{ $p3 }->{ $vn };
	    adjCnt($ref, $sc3v);
	    $cov->{ $p3 } += $cnt3; # unless( $cov->{ $p3 } && $cov->{ $p3 } > $cnt3 );
	    $cov->{ $p3 } = $cov->{ $bp } if( $cov->{ $bp } && $cov->{ $p3 } < $cov->{ $bp } );
	    print STDERR "    Found INV: $p3 BP: $bp Cov: $cov->{ $p3 } $cov->{ $bp } $vn EXTRA: $EXTRA Cnt: $cnt3\n" if ( $opt_y );
	}
    }
}

# Given a candidate SV identified by clipped reads, check whether there're mates to support it
sub checkPairs {
    my ($chr, $start, $end, $sv) = @_;
    my ($pairs, $pmean, $qmean, $Qmena, $nm) = (0, 0, 0, 0, 0);
	    #print STDERR "Pair ($pairs, $pmean, $qmean, $Qmena, $nm)\n";
    foreach my $svcluster (@$sv) {
	foreach my $svr (@$svcluster) {
	    next if ( $svr->{ used } );
	    my $s = ($svr->{ start } + $svr->{ end })/2;
	    my $e = ($svr->{ mstart } + $svr->{ mend })/2;
	    ($s, $e) = ($e, $s) if ( $s > $e );
	    next unless( isOverlap($start, $end, $s, $e) );
	    ($pairs, $pmean, $qmean, $Qmena, $nm) = ($svr->{ cnt }, $svr->{ pmean }, $svr->{ qmean }, $svr->{ Qmean }, $svr->{ nm }) if ( $svr->{ cnt } > $pairs );
	    $svr->{ used } = 1;
	    print STDERR "      Pair [$chr:$s-$e] overlapping [$chr:$start-$end] found and marked.\n" if ( $opt_y );
	}
    }
    return ($pairs, $pmean, $qmean, $Qmena, $nm);
}

# Determine 
sub isOverlap {
    my ($s1, $e1, $s2, $e2) = @_;
    return 0 unless( $s1 < $e2 && $s2 < $e1 );
    my @pos = sort {$a <=> $b} ($s1, $e1, $s2, $e2);
    my $ins = $pos[2] - $pos[1];
    return 1 if ( $ins/($e1 - $s1) > 0.75 && $ins/($e2 - $s2) > 0.75 );
    return 1 if ( $pos[1] - $pos[0] + $pos[3] - $pos[2] < 3 * $RLEN );
    return 0;
}

sub findMatch {
    my ($seq, $REF, $p, $dir, $SEED, $MM) = @_;
    $seq = reverse($seq) if ( $dir == -1 ); # dir==-1 means 5' clip
    $SEED = $SEED ? $SEED : $SEED1;
    print STDERR "    Working Match $p $seq $dir SEED: $SEED\n" if ( $opt_y );
    my $extra = "";
    for(my $i = length($seq) - $SEED; $i >= 0; $i--) {
	my $seed = substr($seq, $i, $SEED);
	if ( $REF->{ $seed } ) {
	    if ( @{ $REF->{ $seed } } == 1 ) {
		my $bp = $dir == 1 ? $REF->{ $seed }->[0] - $i : $REF->{ $seed }->[0] + length($seq) - $i - 1;
		if ( ismatchref($seq, $REF, $bp, $dir, $MM) ) {
		    my $mm = $dir == -1 ? -1 : 0;
		    while( $REF->{ $bp } ne substr($seq, $mm, 1) ) {
			$extra .= substr($seq, $mm, 1);
			$bp += $dir;
			$mm += $dir;
		    }
		    $extra = reverse($extra) if ( $extra && $dir == -1 );
		    print STDERR "      Found SV BP: $dir BP: $bp SEEDpos$REF->{ $seed }->[0] $p $seed $i $seq extra: $extra\n" if ( $opt_y ); 
		    return ($bp, $extra); 
		} else { # for complex indels, allowing some mismatches at the end up to 15bp or 20% length
		    my $sseq = $seq;
		    my $eqcnt = 0;
		    for(my $ii = 1; $ii <= 15; $ii++) {
			$bp += $dir;
			$sseq = $dir == 1 ? substr($sseq, 1) : substr($sseq, 0, -1);
			if( $dir == 1 ) {
			    next unless( $REF->{ $bp } && substr($sseq, 0, 1) eq $REF->{ $bp } ); $eqcnt++;
			    if ( substr($sseq, 1, 1) ne $REF->{ $bp + 1 } ) {
				next;
			    }
			    $extra = substr($seq, 0, $ii);
			} else {
			    next unless( $REF->{ $bp } && substr($sseq, -1, 1) eq $REF->{ $bp } ); $eqcnt++;
			    if ( substr($sseq, -2, 1) ne $REF->{ $bp - 1 } ) {
				next;
			    }
			    $extra = substr($seq, -$ii);
			}
			last if  ($eqcnt >= 3 && $eqcnt/$ii > 0.5);
			print STDERR "      FoundSEED SV BP: $dir BP: $bp SEEDpos$REF->{ $seed }->[0] $p $seed $i $seq EXTRA: $extra\n" if ( $opt_y ); 
			return ($bp, $extra) if ( ismatchref($sseq, $REF, $bp, $dir, 1) );
		    }
		}
	    }
	}
    }
    return (0, "");
}

sub findMatchRev {
    my ($seq, $REF, $p, $dir, $SEED, $MM) = @_;
    $seq = reverse($seq) if ( $dir == 1 ); # $dir = 1 means from 3' soft clip
    $seq =~ y/ATGC/TACG/;
    print STDERR "    Working MatchRev $p $seq $dir\n" if ( $opt_y );
    my $extra = "";
    $SEED = $SEED ? $SEED : $SEED1;
    for(my $i = length($seq) - $SEED; $i >= 0; $i--) {
	my $seed = substr($seq, $i, $SEED);
	if ( $REF->{ $seed } ) {
	    if ( @{ $REF->{ $seed } } == 1 ) {
		my $bp = $dir == 1 ? $REF->{ $seed }->[0] + length( $seq ) - $i -1 : $REF->{ $seed }->[0] - $i;
		if ( ismatchref($seq, $REF, $bp, -1 * $dir, $MM) ) {
		    print STDERR "      Found SV BP (reverse): $dir BP: $bp SEEDpos: $REF->{ $seed }->[0] $p $seed $i $seq\n" if ( $opt_y ); 
		    return ($bp, $extra); 
		} else { # for complex indels, allowing some mismatches at the end up to 15bp or 20% length
		    my $sseq = $seq;
		    my $eqcnt = 0;
		    for(my $ii = 1; $ii <= 15; $ii++) {
			$bp -= $dir;
			$sseq = $dir == -1 ? substr($sseq, 1) : substr($sseq, 0, -1);
			if( $dir == -1 ) {
			    next unless( $REF->{ $bp } && substr($sseq, 0, 1) eq $REF->{ $bp } ); $eqcnt++;
			    if ( substr($sseq, 1, 1) ne $REF->{ $bp + 1 } ) {
				next;
			    }
			    $extra = substr($seq, 0, $ii);
			} else {
			    next unless( $REF->{ $bp } && substr($sseq, -1, 1) eq $REF->{ $bp } ); $eqcnt++;
			    if ( substr($sseq, -2, 1) ne $REF->{ $bp - 1 } ) {
				next;
			    }
			    $extra = substr($seq, -$ii);
			}
			last if  ($eqcnt >= 3 && $eqcnt/$ii > 0.5);
			print STDERR "      FoundSEED SV BP (reverse): $dir BP: $bp SEEDpos$REF->{ $seed }->[0] $p $seed $i $seq EXTRA: $extra\n" if ( $opt_y ); 
			return ($bp, $extra) if ( ismatchref($sseq, $REF, $bp, -1*$dir, 1) );
		    }
		}
	    }
	}
    }
    return (0, "");
}

# test whether the two soft-clipped reads match
# Returns the breakpoints in 5 and 3 prime soft-clipped reads
sub find35match {
    my ($seq5, $seq3, $p5, $p3, $REF) = @_;
    my $LONGMM = 2;
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
		print STDERR "      Found 35 Match, $seq5 $seq3 ", join("\t", $n+$j, length($seq3), $i+$n, length($seq5), $n, $nm, $i, $j), "\n" if ( $opt_y );
		return ($B5, $B3, $max);
	    }
	}
    }
    return ($B5, $B3, $max);
}

# Realign large insertions that are not present in alignment
sub realignlgins {
    my ($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $START, $END, $svfdup, $svrdup) = @_;
    my @tmp;
    while(my($p, $sc5v) = each %$sclip5) {
	push(@tmp, [$p, $sc5v, $sc5v->{cnt}]);
    }
    @tmp = sort {$b->[2] <=> $a->[2];} @tmp;
    foreach my $t (@tmp) {
	my ($p, $sc5v, $cnt) = @$t;
	last if ( $cnt < $MINR );
	next if ($sc5v->{ used }); # already been used in 
	my $seq = findconseq($sc5v);
	next unless( $seq );
	print STDERR "  Working lgins: 5: $p $seq cnt: $cnt\n" if ( $opt_y );
	#next if ($seq =~ /^.AAAAAAAA/ || $seq =~ /^.TTTTTTTT/ );
	next if ( length($seq) < 12 );
	my ($bi, $ins) = findbi($seq, $p, $REF, -1, $chr);
	my $EXTRA = "";
	unless( $bi ) {
	    #next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
	    next if ( islowcomplexseq($seq) );
	    ($bi, $EXTRA) = findMatch($seq, $REF, $p, -1, $SEED1, 1);
	    next unless( $bi && $bi - $p > 15 && $bi - $p < $SVMAXLEN );
	    parseSAM($chr, $bi - $RLEN <= $END ? $END + 1 : $bi - $RLEN, $bi + $RLEN, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1) if( $bi > $END );
	    if ( $bi - $p > $SVMINLEN + 2*$SVFLANK ) {
		$ins = join("", (map { $REF->{ $_ }; } ($p .. ($p+$SVFLANK-1))));
		$ins .= "<dup" . ($bi-$p-2*$SVFLANK+1) .">";
		$ins .= join("", (map {$bi - $_ < length($seq) - length($EXTRA) ? substr($seq, $bi-$_+length($EXTRA), 1) : $REF->{ $_ };} (($bi-$SVFLANK+1) .. $bi)));
	    } else {
		$ins = join("", (map {$bi - $_ < length($seq) - length($EXTRA) ? substr($seq, $bi-$_+length($EXTRA), 1) : $REF->{ $_ };} ($p .. $bi)));
	    }
	    $ins .= $EXTRA;
	    my ($clusters, $pairs) = markDUPSV($p, $bi, [$svfdup, $svrdup]);
	    if ( (! $cov->{ $p-1 }) || ($cov->{ $bi } && $cov->{ $p-1 } < $cov->{ $bi }) ) {
		$cov->{ $p-1 } = $cov->{ $bi } ? $cov->{ $bi } : $sc5v->{ cnt };
	    } else {
		$cov->{ $p-1 } += $sc5v->{ cnt } if( $sc5v->{ cnt } > $cov->{ $p-1 } ) ;
	    }
	    $bi = $p - 1;
	    $hash->{ $bi }->{ SV }->{ type } = "DUP";
	    $hash->{ $bi }->{ SV }->{ pairs } += $pairs;
	    $hash->{ $bi }->{ SV }->{ splits } += $cnt;
	    $hash->{ $bi }->{ SV }->{ clusters } += $clusters;
	}
	print STDERR "  Found candidate lgins from 5: $bi +$ins $p $seq\n" if ($opt_y);
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
	$cov->{ $bi } += $sc5v->{ cnt } unless( $hash->{ $bi }->{ SV } );
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
	$sc5v->{ used } = $bi + $len;
	my %tins = ();
	$tins{ $bi }->{ "+$ins" } = $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt };
	my $newins = realignins($hash, \%tins, $cov, $sclip5, $sclip3, $REF, $chr);
	$newins = "+$ins" unless( $newins );
	$hash->{ $bi }->{ SV }->{ splits } += $hash->{ $bi }->{ I }->{ $newins }->{ cnt } - $tins{ $bi }->{ "+$ins" } if ( $hash->{ $bi }->{ SV } );
	adjCnt($hash->{ $bi }->{ I }->{ $newins }, $hash->{$bi}->{ $REF->{$bi} }, $hash->{$bi}->{ $REF->{$bi} }) if ($rpflag && $bams && length($ins) >= 5 && length($ins) < $RLEN-10 && $hash->{$bi}->{ $REF->{$bi} } && $hash->{$bi}->{ $REF->{$bi} }->{ cnt } && noPassingReads($chr, $bi, $bi+length($ins), $bams) && $hash->{ $bi }->{ I }->{ $newins }->{ cnt } > 2 * $hash->{$bi}->{ $REF->{$bi} }->{ cnt });
    }
    @tmp = ();
    while(my($p, $sc3v) = each %$sclip3) {
	push(@tmp, [$p, $sc3v, $sc3v->{cnt}]);
    }
    @tmp = sort {$b->[2] <=> $a->[2];} @tmp;
    foreach my $t (@tmp) {
	my ($p, $sc3v, $cnt) = @$t;
	last if ( $cnt < $MINR );
	next if ($sc3v->{ used }); # already been used in 
	my $seq = findconseq($sc3v);
	next unless( $seq );
	print STDERR "  Working lgins 3: $p $seq cnt: $cnt\n" if ( $opt_y );
	#next if ($seq =~ /^.AAAAAAA/ || $seq =~ /^.TTTTTTT/ );
	next if ( length($seq) < 12 );
	my ($bi, $ins, $be) = findbi($seq, $p, $REF, 1, $chr);
	#print STDERR "Here $seq $p $sc3v '$bi' $cnt\n";
	#next unless( $bi );
	my $EXTRA = "";
	unless( $bi ) {
	    #next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
	    next if ( islowcomplexseq($seq) );
	    ($bi, $EXTRA) = findMatch($seq, $REF, $p, 1, $SEED1, 1);
	    next unless( $bi && $p - $bi > 15 && $p - $bi < $SVMAXLEN );
	    parseSAM($chr, $bi - $RLEN, $bi + $RLEN >= $START ? $START - 1 : $bi + $RLEN, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1) if ($bi < $START);
	    my $shift5 = 0;
	    while( $REF->{ $p - 1 } eq $REF->{ $bi - 1 } ) {
		$p--; $bi--; $shift5++;
	    }
	    if ( $p - $bi > $SVMINLEN + 2*$SVFLANK ) {
		$ins = join("", (map {$_-$bi >= $shift5 && $_-$bi-$shift5 < length($seq) - length($EXTRA) ? substr($seq, $_-$bi-$shift5+length($EXTRA), 1) : $REF->{ $_ };} ($bi .. ($bi+$SVFLANK-1))));
		$ins .= "<dup" . ($p-$bi-2*$SVFLANK) .">";
		$ins .= join("", (map {$REF->{ $_ };} (($p-$SVFLANK) .. ($p-1))));
	    } else {
		$ins = join("", (map {$_-$bi >= $shift5 && $_-$bi-$shift5 < length($seq) - length($EXTRA) ? substr($seq, $_-$bi-$shift5+length($EXTRA), 1) : $REF->{ $_ };} ($bi .. ($p-1))));
	    }
	    $ins .= $EXTRA;
	    my ($clusters, $pairs) = markDUPSV($bi, $p-1, [$svfdup, $svrdup]);
	    $bi = $bi - 1;
	    $hash->{ $bi }->{ SV }->{ type } = "DUP";
	    $hash->{ $bi }->{ SV }->{ pairs } += $pairs;
	    $hash->{ $bi }->{ SV }->{ splits } += $cnt;
	    $hash->{ $bi }->{ SV }->{ clusters } += $clusters;
	    if ( (! $cov->{ $bi }) || ($cov->{ $p } && $cov->{ $bi } < $cov->{ $p }) ) {
		$cov->{ $bi } = $cov->{ $p } ? $cov->{ $p } : $sc3v->{ cnt };
	    } else {
		$cov->{ $bi } += $sc3v->{ cnt } if( $sc3v->{ cnt } > $cov->{ $bi });
	    }
	}
	print STDERR "  Found candidate lgins from 3: $bi +$ins $p $seq\n" if ( $opt_y );
	$hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } = 0 unless( $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } );
	$hash->{ $bi }->{ I }->{ "+$ins" }->{ pstd } = 1;
	$hash->{ $bi }->{ I }->{ "+$ins" }->{ qstd } = 1;
	    #print STDERR "3B: $bi $cov->{ $bi } $ins $sc3v->{cnt}\n";
	adjCnt($hash->{ $bi }->{ I }->{ "+$ins" }, $sc3v, $p-$bi > $sc3v->{pmean}/$cnt ? undef : $hash->{ $bi }->{ $REF->{ $bi } });
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
	$hash->{ $bi }->{ SV }->{ splits } += $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } - $tins{ $bi }->{ "+$ins" } if ( $hash->{ $bi }->{ SV } );
	adjCnt($hash->{ $bi }->{ I }->{ "+$ins" }, $hash->{$bi}->{ $REF->{$bi} }, $hash->{$bi}->{ $REF->{$bi} }) if ($rpflag && $bams && length($ins) >= 5 && length($ins) < $RLEN-10 && $hash->{$bi}->{ $REF->{$bi} } && $hash->{$bi}->{ $REF->{$bi} }->{ cnt } && noPassingReads($chr, $bi, $bi+length($ins), $bams) && $hash->{ $bi }->{ I }->{ "+$ins" }->{ cnt } > 2 * $hash->{$bi}->{ $REF->{$bi} }->{ cnt });
    }
}

# Realign large deletions that are not present in alignment
sub realignlgdel {
    my ($hash, $cov, $sclip5, $sclip3, $REF, $chr, $bams, $START, $END, $svfdel, $svrdel) = @_;
    my $LONGMM = 3;
    my @tmp = ();
    while(my($p, $sc5v) = each %$sclip5) {
	push(@tmp, [$p, $sc5v, $sc5v->{cnt}]);
    }
    @tmp = sort {$b->[2] <=> $a->[2];} @tmp;
    my ($svcov, $clusters, $pairs) = (0, 0, 0);
    foreach my $t (@tmp) {
	my ($p, $sc5v, $cnt) = @$t;
	last if ( $cnt < $MINR );
	next if ($sc5v->{ used }); # already been used in 
	my $seq = findconseq($sc5v, 5);
	next unless( $seq );
	#next if ($seq =~ /^.AAAAAAA/ || $seq =~ /^.TTTTTTT/ );
	next if ( length($seq) < 7 );
	print STDERR "  Working Realignlgdel: 5' $p '$seq' $cnt\n" if ( $opt_y );
	my $bp = findbp($seq, $p - 5, $REF, $INDELSIZE, -1, $chr);
	my ($extra, $EXTRA) = ("", "");
	unless( $bp ) {
	    #next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
	    next if ( islowcomplexseq($seq) );
	    ($bp, $EXTRA) = findMatch($seq, $REF, $p, -1, $SEED1, 1);
	    next unless( $bp && $p - $bp > 15 && $p - $bp < $SVMAXLEN );
	    $bp++;
	    ($svcov, $clusters, $pairs) = markSV($bp, $p, [$svfdel, $svrdel]);
	    unless( $svcov ) {
		next if ( $cnt <= $MINR );
	    }
	    $hash->{ $bp }->{ SV }->{ type } = "DEL";
	    $hash->{ $bp }->{ SV }->{ pairs } += $pairs;
	    $hash->{ $bp }->{ SV }->{ splits } += $cnt;
	    $hash->{ $bp }->{ SV }->{ clusters } += $clusters;
	    parseSAM($chr, $bp - $RLEN, $bp + $RLEN >= $START ? $START - 1 : $bp + $RLEN, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1) if( $bp < $START );
	}
	my $dellen = $p - $bp;
	my $en = 0;
	my $gt = -$dellen;
	unless( $EXTRA ) {
	    while( substr($seq, $en, 1) ne $REF->{ $bp - $en - 1 } && $en < length($seq) ) {
		$extra .= substr($seq, $en, 1);
		$en++;
	    }
	    if ( $extra ) {
		$extra = reverse($extra);
		$gt = "-$dellen&$extra";
		$bp -= length($extra);
	    }
	} else {
	    $dellen -= length($EXTRA);
	    $gt = "-$dellen&$EXTRA";
	}
	print STDERR "  Found Realignlgdel: $bp $gt 5' $p $seq $cnt\n" if ( $opt_y );

	# Work on the softclipped read at 3'
	my $n = 0;
	unless( $extra || $EXTRA ) {
	    $n++ while( $REF->{ $bp+$n } && $REF->{ $bp + $dellen + $n} && $REF->{ $bp+$n } eq $REF->{ $bp + $dellen + $n} );
	}
	my $sc3p = $bp + $n;
	my $str = "";
	my $mcnt = 0;
	while( $REF->{ $bp+$n } && $REF->{ $bp+$dellen+$n } && $REF->{ $bp+$n } ne $REF->{ $bp + $dellen + $n} && $mcnt <= $LONGMM) {
	    $str .= $REF->{ $bp + $dellen + $n};
	    $n++; $mcnt++;
	}
	if ( length($str) == 1) {
	    my $nm = 0;
	    $n++ && $nm++ while( $REF->{ $bp+$n } && $REF->{ $bp+$dellen+$n } && $REF->{ $bp+$n } eq $REF->{ $bp + $dellen + $n} );
	    $sc3p = $bp + $n if ( $nm >= 3 && (! $sclip3->{ $sc3p }) );
	}
	if ( $hash->{$bp}->{SV} && ! $sclip3->{ $sc3p } ) {  # likely a false positive
	    unless( $svcov || $cnt > $MINR ) {
		delete $hash->{$bp}->{ SV };
		next;
	    }
	}

	$hash->{$bp}->{ $gt }->{ cnt } = 0 unless( $hash->{$bp}->{ $gt }->{ cnt } );
	$hash->{$bp}->{ $gt }->{ qstd } = 1; # more accurate implementation later
	$hash->{$bp}->{ $gt }->{ pstd } = 1; # more accurate implementation later
	adjCnt( $hash->{$bp}->{ $gt }, $sc5v );
	$sc5v->{ used } = $bp;
	$cov->{ $bp } = $cov->{ $p } unless( $cov->{ $bp } );
	if ( $dellen < $INDELSIZE ) {
	    for(my $tp = $bp; $tp < $bp + $dellen; $tp++) {
		$cov->{$tp} += $sc5v->{ cnt };
	    }
	}
	if ( $sclip3->{ $sc3p } && (! $sclip3->{ $sc3p }->{ used }) ) {
	    $sc3p > $bp ? adjCnt($hash->{$bp}->{ $gt }, $sclip3->{ $sc3p }, $hash->{$bp}->{ $REF->{ $bp } }) : adjCnt($hash->{$bp}->{ $gt }, $sclip3->{ $sc3p });
	    if ( $sc3p == $bp ) {
		if ( $dellen < $INDELSIZE ) {
		    for(my $tp = $bp; $tp < $bp + $dellen; $tp++) {
			$cov->{$tp} += $sclip3->{ $sc3p }->{ cnt };
		    }
		}
	    }
	    for(my $ip = $bp+1; $ip < $sc3p; $ip++) {
		rmCnt($hash->{$ip}->{ $REF->{$dellen + $ip} }, $sclip3->{$sc3p});
		delete $hash->{$ip}->{ $REF->{$dellen + $ip} } if ( $hash->{$ip}->{ $REF->{$dellen + $ip} } && $hash->{$ip}->{ $REF->{$dellen + $ip} }->{ cnt } == 0);
		delete $hash->{$ip} if ( (keys %{$hash->{$ip}} ) == 0);
	    }
	    $sclip3->{ $sc3p }->{ used } = $bp;
	}
	my %dels5;
	$dels5{ $bp }->{$gt} = $hash->{$bp}->{ $gt }->{cnt};
	realigndel($hash, \%dels5, $cov, $sclip5, $sclip3, $REF, $chr, $bams);
	$hash->{$bp}->{ SV }->{ splits } += $hash->{$bp}->{ $gt }->{cnt} - $dels5{ $bp }->{$gt} if ( $hash->{$bp}->{ SV } );
	if ( $svcov > $hash->{$bp}->{ $gt }->{ cnt } ) {
	    addVarFactor($hash->{$bp}->{ $gt }, ($svcov - $hash->{$bp}->{ $gt }->{cnt})/$hash->{$bp}->{ $gt }->{cnt});
	}
	print STDERR "  Found lgdel done: $bp $gt $p 5' $seq $hash->{$bp}->{ $gt }->{cnt}\n\n" if ( $opt_y );
    }

    # Work on 3' clipped reads
    @tmp = ();
    ($svcov, $clusters, $pairs) = (0, 0, 0);
    while(my($p, $sc3v) = each %$sclip3) {
	push(@tmp, [$p, $sc3v, $sc3v->{cnt}]);
    }
    @tmp = sort {$b->[2] <=> $a->[2];} @tmp;
    foreach my $t (@tmp) {
	my ($p, $sc3v, $cnt) = @$t;
	last if ($cnt < $MINR);
	next if ($sc3v->{ used }); # already been used in 
	my $seq = findconseq($sc3v, 3);
	next unless( $seq );
	#next if ($seq =~ /^.AAAAAAA/ || $seq =~ /^.TTTTTTT/ );
	next if ( length($seq) < 7 );
	print STDERR "  Working Realignlgdel: 3' $p '$seq' $cnt\n" if ( $opt_y );
	my $bp = findbp($seq, $p + 5, $REF, $INDELSIZE, 1, $chr);
	my ($extra, $EXTRA) = ("", "");
	unless( $bp ) {
	    #next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
	    next if ( islowcomplexseq($seq) );
	    ($bp, $EXTRA) = findMatch($seq, $REF, $p, 1, $SEED1, 1);
	    next unless( $bp && $bp - $p > 15 && $p - $bp < $SVMAXLEN );
	    ($svcov, $clusters, $pairs) = markSV($p, $bp, [$svfdel, $svrdel]);
	    unless( $svcov ) {
		next if ( $cnt <= $MINR ); # a little more stringent
	    }
	    $hash->{ $p }->{ SV }->{ type } = "DEL";
	    $hash->{ $p }->{ SV }->{ pairs } += $pairs;
	    $hash->{ $p }->{ SV }->{ splits } += $cnt;
	    $hash->{ $p }->{ SV }->{ clusters } += $clusters;
	    parseSAM($chr, $bp - $RLEN <= $END ? $END + 1 : $bp - $RLEN, $bp + $RLEN, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1) if ( $bp > $END );
	}
	my $dellen = $bp - $p;
	my $en = 0;
	if ( $EXTRA ) {
	    $dellen -= length($EXTRA);
	} else {
	    while( substr($seq, $en, 1) ne $REF->{ $bp + $en } && $en < length($seq) ) {
		$extra .= substr($seq, $en, 1);
		$en++;
	    }
	}
	my $gt = -$dellen;
	my $sc5p = $bp;
	$bp = $p; # Set it to 5'
	if ( $extra ) {
	    $gt = "-$dellen&$extra";
	    $sc5p += length($extra);
	} elsif ($EXTRA ) {
	    $gt = "-$dellen&$EXTRA";
	} else { # 5' adjustment
	    $bp-- && $sc5p-- while($REF->{ $bp - 1 } eq $REF->{ $bp + $dellen - 1 });
	    if ( $bp != $p && $hash->{ $p }->{ SV } ) {
		$hash->{ $bp }->{ SV } = $hash->{ $p }->{ SV };
		delete $hash->{ $p }->{ SV };
	    }
	}
	if ( $hash->{$bp}->{ SV } && (!$sclip5->{ $sc5p }) ) {
	    unless( $svcov || $cnt > $MINR ) {
		delete $hash->{ $bp }->{ SV } if ( $hash->{ $bp }->{ SV } );
		next;
	    }
	}
	print STDERR "  Found Realignlgdel: bp: $bp $gt 3' $p 5'clip: $sc5p '$seq' $cnt\n" if ( $opt_y );
	$hash->{$bp}->{ $gt }->{ cnt } = 0 unless( $hash->{$bp}->{ $gt }->{ cnt } );
	$hash->{$bp}->{ $gt }->{ qstd } = 1; # more accurate implementation later
	$hash->{$bp}->{ $gt }->{ pstd } = 1; # more accurate implementation later
	if ( $dellen < $INDELSIZE ) {
	    for(my $tp = $bp; $tp < $bp + $dellen + length($extra) + length($EXTRA); $tp++) {
		$cov->{$tp} += $sc3v->{ cnt };
	    }
	}
	$cov->{$bp} = $cov->{ $p - 1 } ? $cov->{ $p - 1 } : $sc3v->{ cnt } unless( $cov->{ $bp } );
	$sc3v->{ pmean } += $dellen*$sc3v->{ cnt };
	adjCnt( $hash->{$bp}->{ $gt }, $sc3v );
	$sc3v->{ used } = $p + $dellen;

	my %dels5;
	$dels5{ $bp }->{$gt} = $hash->{$bp}->{ $gt }->{cnt};
	realigndel($hash, \%dels5, $cov, $sclip5, $sclip3, $REF, $chr, $bams);
	$hash->{$bp}->{ SV }->{ splits } += $hash->{$bp}->{ $gt }->{cnt} - $dels5{ $bp }->{$gt} if ( $hash->{$bp}->{ SV } );
	print STDERR "  Found lgdel: $bp $gt $p 3' '$seq' $hash->{$bp}->{ $gt }->{cnt}\n\n" if ( $opt_y );
	if ( $svcov > $hash->{$bp}->{ $gt }->{ cnt } ) {
	    addVarFactor($hash->{$bp}->{ $gt }, ($svcov - $hash->{$bp}->{ $gt }->{cnt})/$hash->{$bp}->{ $gt }->{cnt});
	}
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
    $vref->{ 1 } -= $tv->{ 1 } ? $tv->{ 1 } : 0;
    $vref->{ -1 } -= $tv->{ -1 } ? $tv->{ -1 } : 0;
    foreach my $k (qw(cnt hicnt locnt pmean qmean Qmean 1 -1)) {
	$vref->{ $k } = 0 if ( $vref->{ $k } && $vref->{ $k } < 0 );
    }
}

# Find the consensus sequence in soft-clipped reads.  Consensus is called if
# the matched nucleotides are >90% of all softly clipped nucleotides.
sub findconseq {
    my ($scv, $dir) = @_;
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
	    if ( $ncnt > $max || $scv->{ seq }->[$i]->{ $nt }->{ qmean } > $maxq) {
		$max = $ncnt;
		$mnt = $nt;
		$maxq = $scv->{ seq }->[$i]->{ $nt }->{ qmean };
	    }
	}
	last if ( $i == 3 && @{ $scv->{ nt } } >= 6 && $tt/$scv->{ cnt } < 0.2 && $tt <= 2);
	unless ( ($tt-$max <= 2 && $max > $tt - $max) || $max/$tt >= 0.8) {
	    last if ( $flag );
	    $flag = 1;
	}
	$total += $tt;
	$match += $max;
	$seq .= $mnt;
    }
    my $SEQ = ($total && $match/$total > 0.9 && length($seq)/1.5 > @{ $scv->{ nt } } - length($seq) && (length($seq)/(@{ $scv->{ nt } } +0) > 0.8 || @{ $scv->{ nt } } - length($seq) < 10 || length($seq) > 25)) ? $seq : "";
    if ( $SEQ && length($SEQ) > $SEED2) {
	$scv->{ used } = 1 if ($SEQ =~ /^.AAAAAAA/ || $SEQ =~ /^.TTTTTTT/ );
	$scv->{ used } = 1 if (islowcomplexseq($SEQ));
    }
    if ( $dir && $SEQ && length($SEQ) >= $ADSEED ) {
	if ( $dir == 3 ) { # 3'
	    if ( $adaptor{ substr($SEQ, 0, $ADSEED) } ) {
		$SEQ = "";
	    }
	} elsif ( $dir == 5 ) { # 5'
	    if ( $adaptor_rev{ reverse(substr($SEQ, 0, $ADSEED)) } ) {
		$SEQ = "";
	    }
	}
    }
    $scv->{ SEQ } = $SEQ;
    print STDERR "  Candidate consensus: $seq Reads: $scv->{ cnt } M: $match T: $total Final: $SEQ\n" if ( $opt_y );
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
	my $extrains = "";
	my $extra = "";
	my ($inv5, $inv3) = ("", "");
	if ( $vn =~ /^-\d+\^([ATGNC]+)<...\d+>([ATGNC]+)$/o ) {
	    $inv5 = $1;
	    $inv3 = $2;
	} elsif ($vn =~ /^-\d+(.*)/) {
	    $extra = $1;
	    $extra =~ s/\^|&|#//g;
	    $extrains = $1 if ($vn =~ /\^([ATGNC]+)/o);
	}
	#my $wustart = ($p - $dellen - 100) > 1 ? ($p - $dellen - 100) : 1;
	my $wustart = ($p - 200) > 1 ? ($p - 200) : 1;
	my $wupseq = join( "", (map { $REF->{ $_ } ? $REF->{ $_ } : ""; } ($wustart .. ($p-1)))) . $extra; # . $extra; # 5' flanking seq
	$wupseq = $inv3 if ( $inv3 );
	#my $sanend = ($p + 2*$dellen + 100) > $CHRS{$chr} ? $CHRS{ $chr } : ($p + 2*$dellen + 100);
	my $sanend = ($p + 200) > $CHRS{$chr} ? $CHRS{ $chr } : ($p + 200);
	my $sanpseq = $extra . join( "", (map { $REF->{ $_ }; } (($p + $dellen + length($extra) - length($extrains)) .. $sanend))); # 3' flanking seq
	$sanpseq = $inv5 if ( $inv5 );
	my ($mm3, $sc3p, $nm3, $misp3, $misnt3) = findMM3($REF, $p, $sanpseq, $sclip3); # mismatches, mismatch positions, 5 or 3 ends
	my ($mm5, $sc5p, $nm5, $misp5, $misnt5) = findMM5($REF, $p+$dellen+length($extra)-length($extrains)-1, $wupseq, $sclip5);
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
	delete $hash->{ $misp3 }->{ $misnt3 } if ( $misp3 && $hash->{ $misp3 }->{ $misnt3 } && @$mm3 == 1 && $hash->{ $misp3 }->{ $misnt3 }->{ cnt } < $dcnt );
	delete $hash->{ $misp5 }->{ $misnt5 } if ( $misp5 && $hash->{ $misp5 }->{ $misnt5 } && @$mm5 == 1 && $hash->{ $misp5 }->{ $misnt5 }->{ cnt } < $dcnt );

	#next if ( $dellen < 2 ); # soft-clipping only happens when deletion is at least 2 bp
	foreach my $sc5pp (@$sc5p) {
	    if ( $sclip5->{ $sc5pp } && (! $sclip5->{ $sc5pp }->{ used }) ) {
		my $tv = $sclip5->{ $sc5pp };
		my $seq = findconseq( $tv );
		print STDERR "  Realigndel 5: $p $sc5pp seq: '$seq' Wuseq: ", scalar reverse($wupseq), " cnt: $tv->{ cnt } $dcnt $vn $p cov: $cov->{ $p }\n" if ( $opt_y );
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
		print STDERR "  Realigndel 3: $p $sc3pp seq '$seq' Sanseq: $sanpseq cnt: $tv->{ cnt } $dcnt $vn $p $dellen ", substr($sanpseq, $sc3pp-$p), "\n" if ( $opt_y );
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
    $vref->{ 1 } += $tv->{ 1 } ? $tv->{ 1 } : 0;
    $vref->{ -1 } += $tv->{ -1 } ? $tv->{ -1 } : 0;
    print STDERR "    AdjCnt: '+' $vref->{ cnt } $tv->{ cnt } ", $vref->{ 1 }, " ", $tv->{ 1 } ? $tv->{ 1 } : 0, " Ref: ", $ref && $ref->{ cnt } ? $ref->{ cnt } : "NA", "\n" if ( $opt_y );
    print STDERR "    AdjCnt: '-' $vref->{ cnt } $tv->{ cnt } ", $vref->{ -1 }, " ", $tv->{ -1 } ? $tv->{ -1 } : 0, " Ref: ", $ref && $ref->{ cnt } ? $ref->{ cnt } : "NA", "\n" if ( $opt_y );
    return unless($ref);
    $ref->{ cnt } -= $tv->{ cnt };
    $ref->{ hicnt } -= $tv->{ hicnt } ? $tv->{ hicnt } : 0;
    $ref->{ locnt } -= $tv->{ locnt } ? $tv->{ locnt } : 0;
    $ref->{ pmean } -= $tv->{ pmean };
    $ref->{ qmean } -= $tv->{ qmean };
    $ref->{ Qmean } -= $tv->{ Qmean };
    $ref->{ nm } -= $tv->{ nm };
    $ref->{ 1 } -= $tv->{ 1 } ? $tv->{ 1 } : 0;
    $ref->{ -1 } -= $tv->{ -1 } ? $tv->{ -1 } : 0;
    foreach my $k (qw(cnt hicnt locnt pmean qmean Qmean 1 -1)) {
	$ref->{ $k } = 0 if ( $ref->{ $k } && $ref->{ $k } < 0 );
    }
}

# Adjust the reference count.
sub adjRefCnt {
    my ($tv, $ref, $len) = @_;
    return unless($ref);
    print STDERR "    AdjRefCnt: '+' $ref->{ cnt } $tv->{ cnt } ", $ref->{ 1 }, " ", $tv->{ 1 } ? $tv->{ 1 } : 0, " Ref: ", $ref && $ref->{ cnt } ? $ref->{ cnt } : "NA", "\n" if ( $opt_y );
    print STDERR "    AdjRefCnt: '-' $ref->{ cnt } $tv->{ cnt } ", $ref->{ -1 }, " ", $tv->{ -1 } ? $tv->{ -1 } : 0, " Ref: ", $ref && $ref->{ cnt } ? $ref->{ cnt } : "NA", "\n" if ( $opt_y );
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
    $ref->{ 1 } -= $tv->{ 1 } ? int($f * $tv->{ 1 }) : 0;
    $ref->{ -1 } -= $tv->{ -1 } ? int($f * $tv->{ -1 }) : 0;
    foreach my $k (qw(cnt hicnt locnt pmean qmean Qmean 1 -1)) {
	$ref->{ $k } = 0 if ( $ref->{ $k } && $ref->{ $k } < 0 );
    }
}

# Adjust the reference by factor
sub adjRefFactor {
    my ($ref, $f) = @_;
    return unless($ref);
    $f = 1 if ( $f > 1 );
    return if ( $f < -1 );
    print STDERR "    AdjRefFactor: $ref->{ cnt } $f\n" if ( $opt_y );
    $ref->{ cnt } -= int($f*$ref->{ cnt });
    $ref->{ hicnt } -= $ref->{ hicnt } ? int($f*$ref->{ hicnt }) : 0;
    $ref->{ locnt } -= $ref->{ locnt } ? int($f*$ref->{ locnt }) : 0;
    $ref->{ pmean } -= $f * $ref->{ pmean };
    $ref->{ qmean } -= $f * $ref->{ qmean };
    $ref->{ Qmean } -= $f * $ref->{ Qmean };
    $ref->{ nm } -= $f * $ref->{ nm };
    $ref->{ 1 } -= $ref->{ 1 } ? int($f * $ref->{ 1 }) : 0;
    $ref->{ -1 } -= $ref->{ -1 } ? int($f * $ref->{ -1 }) : 0;
    foreach my $k (qw(cnt hicnt locnt pmean qmean Qmean 1 -1)) {
	$ref->{ $k } = 0 if ( $ref->{ $k } && $ref->{ $k } < 0 );
    }
}

# Add variant by factor
sub addVarFactor {
    my ($ref, $f) = @_;
    return unless($ref);
    return if ( $f < -1 );
    $ref->{ cnt } += int($f*$ref->{ cnt });
    $ref->{ hicnt } += $ref->{ hicnt } ? int($f*$ref->{ hicnt }) : 0;
    $ref->{ locnt } += $ref->{ locnt } ? int($f*$ref->{ locnt }) : 0;
    $ref->{ pmean } += $f * $ref->{ pmean };
    $ref->{ qmean } += $f * $ref->{ qmean };
    $ref->{ Qmean } += $f * $ref->{ Qmean };
    $ref->{ nm } += $f * $ref->{ nm };
    $ref->{ 1 } += $ref->{ 1 } ? int($f * $ref->{ 1 }) : 0;
    $ref->{ -1 } += $ref->{ -1 } ? int($f * $ref->{ -1 }) : 0;
}

# Given a variant sequence, find the mismatches and potential softclipping positions
sub findMM3 {
    my ($REF, $p, $seq, $sclip3) = @_;
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
		$sclip3->{ $p+$n+$n2 }->{ used } = 1 if ( $sclip3->{ $p+$n+$n2 } );
		$mn += $n2;
	    } else {
		push(@sc3p, $p + $n);
		$sclip3->{ $p+$n }->{ used } = 1 if ( $sclip3->{ $p+$n } );
	    }
	}
    }
    return (\@mm, \@sc3p, $mn, $misp, $misnt);
}

sub findMM5 {
    my ($REF, $p, $seq, $sclip5) = @_;
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
		$sclip5->{ $p-$n-$n2 }->{ used } = 1 if ( $sclip5->{ $p-$n-$n2 } );
		$mn += $n2;
	    } else {
		push(@sc5p, $p - $n);
		$sclip5->{ $p-$n }->{ used } = 1 if ( $sclip5->{ $p-$n } );
	    }
	}
    }
    return (\@mm, \@sc5p, $mn, $misp, $misnt);
}

sub realignins {
    my ($hash, $ins, $cov, $sclip5, $sclip3, $REF, $chr) = @_;
    my $NEWINS = "";
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
	my $ins3 = "";
	my $inslen = length($ins);
	if ( $vn =~ /<dup(\d+)>([ATGC]+)$/ ) {
	    $ins3 = $2;
	    $inslen += $1 + length($2);
	}
	my $extra = ($vn =~ /&([ATGC]+)/) ? $1 : "";
	my $compm = ($vn =~ /#([ATGC]+)/) ? $1 : "";  # the match part for a complex variant
	my $newins = ($vn =~ /\^([ATGC]+)$/) ? $1 : "";  # the adjacent insertion
	my $newdel = ($vn =~ /\^(\d+)$/) ? $1 : 0;  # the adjacent deletion
	my $tn = $vn;  $tn =~ s/^\+//; $tn =~ s/&//; $tn =~ s/#//; $tn =~ s/\^\d+$//; $tn =~ s/\^//;
	my $wustart = ($p - 150) > 1 ? ($p - 150) : 1;
	my $wupseq = join( "", (map { $REF->{ $_ } ? $REF->{ $_ } : ""; } ($wustart .. $p))) . $tn; # 5prime flanking seq
	my $sanend = $CHRS{ $chr } < ($p + length($vn) + 100) ? $CHRS{ $chr } : ($p + length($vn) + 100);
	my $sanpseq = "";
	my ($mm3, $sc3p, $nm3, $misp3, $misnt3);
	if ( $ins3 ) {
	    my $p3 = $p + $inslen - length($ins3) + $SVFLANK;
	    $sanpseq = substr($ins3, $SVFLANK - length($ins3)) if ( length($ins3) > $SVFLANK );
	    $sanpseq .= join("", (map { $REF->{ $_ }; } (($p+1) .. ($p + 101))));
	    ($mm3, $sc3p, $nm3, $misp3, $misnt3) = findMM3($REF, $p3+1, $sanpseq, $sclip3); # mismatches, mismatch positions, 5 or 3 ends
	} else {
	    $sanpseq = $tn . join( "", (map { $REF->{ $_ }; } (($p+length($extra)+1+length($compm)+$newdel) .. $sanend))); # 3prime flanking seq
	    ($mm3, $sc3p, $nm3, $misp3, $misnt3) = findMM3($REF, $p+1, $sanpseq, $sclip3); # mismatches, mismatch positions, 5 or 3 ends
	}
	my ($mm5, $sc5p, $nm5, $misp5, $misnt5) = findMM5($REF, $p+length($extra)+length($compm)+$newdel, $wupseq, $sclip5);
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
	delete $hash->{ $misp3 }->{ $misnt3 } if ( $misp3 && $hash->{ $misp3 }->{ $misnt3 } && @$mm3 == 1 && $hash->{ $misp3 }->{ $misnt3 }->{ cnt } < $icnt );
	delete $hash->{ $misp5 }->{ $misnt5 } if ( $misp5 && $hash->{ $misp5 }->{ $misnt5 } && @$mm5 == 1 && $hash->{ $misp5 }->{ $misnt5 }->{ cnt } < $icnt );

	foreach my $sc5pp (@{ $sc5p }) {
	    print STDERR "    55: $p $sc5pp VN: '$vn'  5' seq: ^$wupseq^\n" if ( $opt_y );
	    if ( $sclip5->{ $sc5pp } && (! $sclip5->{ $sc5pp }->{ used }) ) {
		my $tv = $sclip5->{ $sc5pp };
		my $seq = findconseq( $tv );
		print STDERR "    ins5: $p $sc5pp $seq $wupseq VN: $vn iCnt: $icnt vCnt: $tv->{ cnt }\n" if ( $opt_y );
		if ( $seq && ismatch($seq, $wupseq, -1) ) {
		    print STDERR "      ins5: $p $sc5pp $seq $wupseq VN: $vn iCnt: $icnt vCnt: $tv->{ cnt } used\n" if ( $opt_y );
		    $cov->{ $p } += $tv->{ cnt } if ( $sc5pp > $p );
		    adjCnt($vref, $tv);
		    $sclip5->{ $sc5pp }->{ used } = 1;
		    if ( length($ins) + 1 == length($vn) && $sc5pp <= $p ) {
			#print STDERR "$sc5pp $p\n";
			# To find a case and implement later
		    }
		}
	    }
	}

	foreach my $sc3pp (@$sc3p) {
	    print STDERR "    33: $p $sc3pp VN: $vn'  3' seq: ^$sanpseq^\n" if ( $opt_y );
	    if ( $sclip3->{ $sc3pp } && (! $sclip3->{ $sc3pp }->{ used }) ) {
		my $tv = $sclip3->{ $sc3pp };
		my $seq = findconseq( $tv );
		print STDERR "    ins3: $p $sc3pp $seq $sanpseq VN: $vn iCnt: $icnt vCnt: $tv->{cnt}\n" if ($opt_y);
		my $mseq = $ins3 ? $sanpseq : substr($sanpseq, $sc3pp - $p - 1);
		if ( $seq && ismatch($seq, $mseq, 1) ) {
		    print STDERR "      ins3: $p $sc3pp $seq VN: $vn iCnt: $icnt vCnt: $tv->{cnt} used\n" if ( $opt_y );
		    $cov->{ $p } += $tv->{ cnt } if ( $sc3pp <= $p || length($ins) > $tv->{ pmean }/$tv->{ cnt } );
		    my $ref = $sc3pp <= $p ? "" : $hash->{$p}->{ $REF->{$p} };
		    adjCnt($vref, $tv, length($ins) <= $tv->{ pmean }/$tv->{ cnt } ? $ref : undef);
		    $sclip3->{ $sc3pp }->{ used } = 1;
		    if ( length($ins) + 1 == length($vn) && length($ins) > $RLEN && $sc3pp >= $p + 1 + length($ins) ) {
			my $flag = 0;
			my $offset = ($sc3pp-$p-1) % length($ins);
			my $tvn = $vn;
			for(my $seqi = 0; $seqi < length($seq) && $seqi < length($ins); $seqi++ ) {
			    if ( substr($seq, $seqi, 1) ne substr($ins, $seqi + $offset, 1) ) {
				$flag++;
				substr($tvn, $seqi + $offset + 1, 1) = substr($seq, $seqi, 1);
			    }
			}
			if ( $flag ) {
			    $hash->{ $p }->{ I }->{ $tvn } = $hash->{ $p }->{ I }->{ $vn };
			    delete $hash->{ $p }->{ I }->{ $vn };
			    $NEWINS = $tvn;
			}
		    }
		}
	    }
	}
	if ( @$sc3p && @$sc5p && $sc3p->[0] > $sc5p->[0] + 3 && $sc3p->[0] - $sc5p->[0] < $RLEN * 0.75) {
	    adjRefFactor($hash->{$p}->{ $REF->{$p} }, ($sc3p->[0] - $sc5p->[0] - 1)/$RLEN);
	    adjRefFactor($vref, -($sc3p->[0] - $sc5p->[0] - 1)/$RLEN);
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
    return $NEWINS;
}

sub ismatch {
    my ($seq1, $seq2, $dir, $MM) = @_;
    $MM = 3 unless( defined($MM) );
    print STDERR "      Matching two seqs $seq1 $seq2 $dir MM: $MM\n" if ( $opt_y );
    $seq2 =~ s/#|\^//g;
    my ($mm, $n) = (0, 0);
    for(my $n = 0; $n < length($seq1) && $n < length($seq2); $n++) {
	$mm++ if ( substr($seq1, $n, 1) ne substr($seq2, $dir*$n - ($dir == -1 ? 1 : 0), 1) );
    }
    return ($mm <= $MM && $mm/length($seq1) < 0.15) ? 1 : 0;
}

sub ismatchref {
    my ($seq, $REF, $p, $dir, $MM) = @_;
    $MM = 3 unless( defined($MM) );
    print STDERR "      Matching REF $seq $p $dir MM: $MM\n" if ( $opt_y );
    my ($mm, $n) = (0, 0);
    for(my $n = 0; $n < length($seq); $n++) {
	return 0 unless( $REF->{ $p + $dir*$n } );
	$mm++ if ( substr($seq, $dir == 1 ? $n : $dir * $n - 1, 1) ne $REF->{ $p + $dir*$n } );
    }
    return ($mm <= $MM && $mm/length($seq) < 0.15) ? 1 : 0;
}

sub strandBias {
    my ($fwd, $rev) = @_;
    if ( $fwd + $rev <= 12 ) {  # using p=0.01, because prop.test(1,12) = 0.01
	return $fwd*$rev > 0 ? 2 : 0;
    }
    return ($fwd/($fwd+$rev) >= $BIAS && $rev/($fwd+$rev) >= $BIAS && $fwd >= $MINB && $rev >= $MINB) ? 2 : 1;
}

# Output remaining soft-clipped reads that haven't been used
sub outputClipping {
    my ($sc5, $sc3) = @_;
    print STDERR "5' Remaining clipping reads\n";
    while( my ($p, $sc) = each %$sc5 ) {
	next if ( $sc->{ used } );
	next if ( $sc->{ cnt } < $MINR );
	my $seq = findconseq( $sc );
	if ( $seq && length( $seq ) > $SEED2 ) {
	    $seq = reverse( $seq );
	    print STDERR "  P: $p Cnt: $sc->{ cnt } Seq: $seq\n";
	}
    }
    print STDERR "3' Remaining clipping reads\n";
    while( my ($p, $sc) = each %$sc3 ) {
	next if ( $sc->{ used } );
	next if ( $sc->{ cnt } < $MINR );
	my $seq = findconseq( $sc );
	if ( $seq && length( $seq ) > $SEED2 ) {
	    print STDERR "  P: $p Cnt: $sc->{ cnt } Seq: $seq\n";
	}
    }
}

#FCB02N4ACXX:3:2206:20108:2526#GATGGTTC  163     chr3    38181981	50      79M188N11M      =       38182275	667     TGAAGTTGTGTGTGTCTGACCGCGATGTCCTGCCTGGCACCTGTGTCTGGTCTATTGCTAGTGAGCTCATCGTAAAGAGGTGCCGCCGGG \YY`c`\ZQPJ`e`b]e_Sbabc[^Ybfaega_^cafhR[U^ee[ec][R\Z\__ZZbZ\_\`Z`d^`Zb]bBBBBBBBBBBBBBBBBBB AS:i:-8 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:72A16A0    YT:Z:UU XS:A:+  NH:i:1  RG:Z:15
sub USAGE {
    print STDERR <<USAGE;
    $0 [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] [-x #_nu] [-g gene] [-f freq] [-r #_reads] [-B #_reads] region_info

    VarDict is a variant calling program for SNV, MNV, indels (<120 bp default, but can be set using -I option), and complex variants.  It accepts any BAM format, either
    from DNA-seq or RNA-seq.  There are several distinct features over other variant callers.  First, it can perform local
    realignment over indels on the fly for more accurate allele frequencies of indels.  Second, it rescues softly clipped reads
    to identify indels not present in the alignments or support existing indels.  Third, when given the PCR amplicon information,
    it will perform amplicon-based variant calling and filter out variants that show amplicon bias, a common false positive in PCR
    based targeted deep sequencing.  Forth, it has very efficient memory management and memory usage is linear to the region of
    interest, not the depth.  Five, it can handle ultra-deep sequencing and the performance is only linear to the depth.  It has
    been tested on depth over 2M reads.  Finally, it has a build-in capability to perform paired sample analysis, intended for
    somatic mutation identification, comparing DNA-seq and RNA-seq, or resistant vs sensitive in cancer research.  By default,
    the region_info is an entry of refGene.txt from IGV, but can be any region or bed files.

    -H|-? Print this help page

    -h|--header
       Print a header row decribing columns

    -v VCF format output

    -i|--splice 
       Output splicing read counts

    -p Do pileup regarless the frequency

    -C Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2 (deprecated)

    -D|--debug
       Debug mode.  Will print some error messages and append full genotype at the end.

    -y|verbose
       Verbose mode.  Will output variant calling process.

    -M Similar to -D, but will append individual quality and position data instead of mean

    -t|--dedup
       Indicate to remove duplicated reads.  Only one pair with same start positions will be kept

    -3 Indicate to move indels to 3-prime if alternative alignment can be achieved.

    -U|--nosv Turn off structural variant calling

    -u Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using forward read only.

    -F bit
       The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and duplicates).  Use -F 0 to turn it off.

    -z 0/1
       Indicate wehther is zero-based cooridates, as IGV does.  Default: 1 for BED file or amplicon BED file.  Use 0 to turn it off.
       When use -R option, it is set to 0
       
    -a|--amplicon int:float
       Indicate it is amplicon based calling.  Reads do not map to the amplicon will be skipped.  A read pair is considered belonging
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
       The threshold for allele frequency, default: 0.01 or 1%
       
    -r minimum reads
       The minimum # of variance reads, default 2
       
    -B INT
       The minimum # of reads to determine strand bias, default 2
       
    -Q INT
       If set, reads with mapping quality less than INT will be filtered and ignored
       
    -q INT
       The phred score for a base to be considered a good call.  Default: 22.5 (for Illumina)
       For PGM, set it to ~15, as PGM tends to under estimate base quality.
       
    -m INT
       If set, reads with mismatches more than INT will be filtered and ignored.  Gaps are not counted as mismatches.  
       Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem is calculated as NM - Indels.  Default: 8,
       or reads with more than 8 mismatches will not be used.
   
    -T|--trim INT
       Trim bases after [INT] bases in the reads
       
    -X INT
       Extension of bp to look for mismatches after insersion or deletion.  Default to 3 bp, or only calls when they are within 3 bp.
       
    -P number
       The read position filter.  If the mean variants position is less that specified, it is considered false positive.  Default: 5
       
    -Z|--downsample double
       For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with caution.  The
       downsampling will be random and non-reproducible.
       
    -o Qratio
       The Qratio of (good_quality_reads)/(bad_quality_reads+0.5).  The quality is defined by -q option.  Default: 1.5
       
    -O MapQ
       The reads should have at least mean MapQ to be considered a valid variant.  Default: no filtering
       
    -V freq
       The lowest frequency in normal sample allowed for a putative somatic mutations.  Default to 0.05
       
    -I INT
       The indel size.  Default: 50bp
       
    -M INT
       The minimum matches for a read to be considered.  If, after soft-clipping, the matched bp is less than INT, then the 
       read is discarded.  It is meant for PCR based targeted sequencing where there is no insert and the matching is only the primers.
       Default: 25, or reads with matches less than 25bp will be filtered
       
    -L INT
       The minimum structural variant length to be presented using <DEL> <DUP> <INV> <INS>, etc.  Default: 500.  Any indel, complex
       variants less than this will be spelled out with exact nucleotides
       
    -w|--insert-size INSERT_SIZE
       The insert size.  Used for SV calling.  Default: 300
       
    -W|--insert-std INSERT_STD
       The insert size STD.  Used for SV calling.  Default: 100
       
    -A INSERT_STD_AMT
       The number of STD.  A pair will be considered for DEL if INSERT > INSERT_SIZE + INSERT_STD_AMT * INSERT_STD.  Default: 4
       
    -J|--crispr CRISPR_cutting_site
       The genomic position that CRISPR/Cas9 suppose to cut, typically 3bp from the PAM NGG site and within the guide.  For
       CRISPR mode only.  It will adjust the variants (mostly In-Del) start and end sites to as close to this location as possible,
       if there are alternatives. The option should only be used for CRISPR mode.
       
    -j CRISPR_filtering_bp
       In CRISPR mode, the minimum amount in bp that a read needs to overlap with cutting site.  If a read does not meet the criteria,
       it will not be used for variant calling, since it is likely just a partially amplified PCR.  Default: not set, or no filtering 

    --adaptor adaptor_seq
       Filter adaptor sequences so that they are not used in realignment.  Multiple adaptors can be supplied by multiple of this option.

    --chimeric
       Indicate to turn off chimeric reads filtering.  Chimeric reads are artifacts from library construction, where a read can be split
       into two segments, each will be aligned within 1-2 read length distance, but in opposite direction.

AUTHOR
       Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
       Report bugs to zhongwu\@yahoo.com

COPYRIGHT
       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

USAGE
   exit(0);
}
