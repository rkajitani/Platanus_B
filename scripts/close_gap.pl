#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin";
use my_bio;

(@ARGV != 2) and die "usage: $0 base.fa gap_filled_info_sorted.tsv\n";

open($in, $ARGV[0]);
$i = 0;
while (($name, $seq) = fasta_nonwhite_get($in)) {
	$base_seq{$name} = $seq;
	push(@base_names, $name);
}
close $in;

open($in, "tac $ARGV[1] |");
while (chomp($l = <$in>)) {
	($name, $start, $end, $new_seq) = (split(/\t/, $l))[0, 4, 5, 12];
	if ($end <= length($base_seq{$name}) and substr($base_seq{$name}, $start, $end - $start) !~ /[a-z]/) {
		substr($base_seq{$name}, $start, $end - $start, lc($new_seq));
	}
}
close $in;

for $name (@base_names) {
	$seq = uc($base_seq{$name});
	$len = length($seq);

	$name =~ s/len\d+/len$len/;
	print ">$name\n";
	print_seq($seq, 80);
}
