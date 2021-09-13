#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin";
use my_bio;

(@ARGV != 2) and die "usage: $0 seq.fa around_length\n";

$around_len = $ARGV[1];

open($in, $ARGV[0]);
while (($name, $seq) = fasta_nonwhite_get($in)) {
	$seq = uc($seq);
	while ($seq =~ /[ACGT]+/g) {
		$contig_st = length($`);
		$contig_len = length($&);
		if ($contig_len < $around_len) {
#			substr($seq, $contig_st, $contig_len, lc(substr($seq, $contig_st, $contig_len)));
			push(@pos_buf, [$contig_st, $contig_len]);
		}
	}

	for $pos ($pos_buf) {
		substr($seq, $pos->[0], $pos->[1], lc(substr($seq, $pos->[0], $pos->[1])));
	}

	while ($seq =~ /[Nnacgt]+/g) {
		$gap_st = length($`);
		$gap_len = length($&);

		if ($gap_st - $around_len < 0 or
			$gap_st + $gap_len + $around_len > length($seq) or
			substr($seq, $gap_st - $around_len, $around_len) =~ /[Nnacgt]/ or
			substr($seq, $gap_st + $gap_len, $around_len) =~ /[Nnacgt]/) {

			next
		}

		printf(">%s;%d;%d;%s;%s\n", $name, $gap_st, $gap_st + $gap_len, $around_len, 'L');
		print(substr($seq, $gap_st - $around_len, $around_len), "\n");

		printf(">%s;%d;%d;%s;%s\n", $name, $gap_st, $gap_st + $gap_len, $around_len, 'R');
		print(substr($seq, $gap_st + $gap_len, $around_len), "\n");
	}
}
