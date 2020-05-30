#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin";
use my_bio;
use List::Util qw/min max/;

(@ARGV != 2) and die "usage: $0 in.paf target.fa\n";

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
	$t_name_flag{$t_name} = 1;
}
close $in;


open($in, $ARGV[1]);
while (($name, $seq) = fasta_nonwhite_get($in)) {
	if ($t_name_flag{$name}) {
		$t_seq{$name} = $seq;
	}
}
close $in;


while (($gap_id, $aln) = each %flanking_aln) {
	($gap_seq_name, $gap_start, $gap_end) = (split(/;/, $gap_id))[0, 1, 2];

	$inner_aln_start = $gap_start - ($aln->[0]{q_len} - $aln->[0]{q_start});
	$inner_aln_end = $gap_end + $aln->[1]{q_end};

	if ($aln->[0]{strand} eq '+') {
		$left_inner_aln_edge = $aln->[0]{t_start};
		$right_inner_aln_edge = $aln->[1]{t_end};
		$filling_seq = substr($t_seq{$aln->[0]{t_name}}, $left_inner_aln_edge, $right_inner_aln_edge - $left_inner_aln_edge);
	}
	else {
		$left_inner_aln_edge = $aln->[1]{t_start};
		$right_inner_aln_edge = $aln->[0]{t_end};
		$filling_seq = rev_comp(substr($t_seq{$aln->[0]{t_name}}, $left_inner_aln_edge, $right_inner_aln_edge - $left_inner_aln_edge));
	}
	
	if ($right_inner_aln_edge - $left_inner_aln_edge < 0) {
		next;
	}

	print(join("\t", (
		$gap_seq_name,
		$gap_start,
		$gap_end,
		$gap_end - $gap_start,
		$inner_aln_start,
		$inner_aln_end,
		$inner_aln_end - $inner_aln_start,
		$aln->[0]{strand},
		$aln->[0]{t_name},
		$left_inner_aln_edge,
		$right_inner_aln_edge,
		$right_inner_aln_edge - $left_inner_aln_edge,
		$filling_seq,
	)), "\n");
}
