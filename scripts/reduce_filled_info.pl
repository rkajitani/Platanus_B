#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin";
use my_bio;

(@ARGV == 0) and die "usage: $0 base.fa filled_info1.tsv filled_info2.tsv [filled_info3.tsv ...]\n";

open($in, $ARGV[0]);
$i = 0;
while (($name, $seq) = fasta_nonwhite_get($in)) {
	$flag{$name} = '0' x length($seq);
}
close $in;

for $file (@ARGV[1..($#ARGV)]) {
	open($in, $file);
	while (chomp($l = <$in>)) {
		($name, $start, $end) =  (split(/\t/, $l))[0, 1, 2];
		if (substr($flag{$name}, $start, $end - $start) !~ /1/) {
			substr($flag{$name}, $start, $end - $start, '1' x ($end - $start));
			push(@out_line, $l);
		}
	}
	close $in;
}

for (@out_line) {
	print($_, "\n");
}
