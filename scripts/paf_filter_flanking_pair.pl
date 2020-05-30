#!/usr/bin/perl

use List::Util qw/min max/;

(@ARGV != 2) and die "usage: $0 in.paf max_distance(bp)\n";

$max_dist = $ARGV[1];

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($q_name, $q_len, $q_start, $q_end, $strand, $t_name, $t_start, $t_end, $aln_len) = (split(/\t/, $l))[0, 1, 2, 3, 4, 5, 7, 8, 9, 10];

	@aln = (
		q_start =>  $q_start,
		q_end =>  $q_end,
		t_name => $t_name,
		q_len =>  $q_len,
		t_start => $t_start,
		t_end => $t_end,
		strand => $strand,
		line => $l,
	);

	if ($q_name =~ /^(\S+);([LR])$/) {
		$gap_id = $1;
		if ($2 eq 'R') {
			$right_flag = 1;
		}
		else {
			$right_flag = 0;
		}
	}
	else {
		next;
	}

	$flanking_aln{$gap_id}[$right_flag] = {@aln};
}
close $in;

for $aln (values %flanking_aln) {
	if (not(defined($aln->[0]) and defined($aln->[1])) or 
		$aln->[0]{t_name} ne $aln->[1]{t_name} or
		$aln->[0]{strand} ne $aln->[1]{strand}) {

		next;
	}

	if ($aln->[0]{strand} eq '+') {
		$dist =  ($aln->[1]{t_start} - $aln->[1]{q_start}) - ($aln->[0]{t_end} + ($aln->[0]{q_len} - $aln->[0]{q_end}));
	}
	else {
		$dist =  ($aln->[0]{t_start} - ($aln->[0]{q_len} - $aln->[0]{q_end})) - ($aln->[1]{t_end} + $aln->[1]{q_start});
	}

#	if ($dist < 0 or $dist > $max_dist) {
	if (abs($dist) > $max_dist) {
		next;
	}

	for $i (0, 1) {
		print($aln->[$i]{line}, "\n");
	}
}
