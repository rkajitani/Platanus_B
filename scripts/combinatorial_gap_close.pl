#!/usr/bin/perl

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";

$sub_prog_dir = "$FindBin::Bin";
$minimap2_bin = "$sub_prog_dir/minimap2";

$n_threads = 1;
$max_dist = 100000;
$min_idt = 0.90;
$min_qcov = 0.25;
$out_prefix = 'out';
$flank_len_str = '500,1000,5000,10000,20000,40000,80000,160000';
$n_iteration = 2;

GetOptions(
	'b=s' => \$base_fasta,
	'c=s' => \$comp_fasta,
	'o=s' => \$out_prefix,
	't=i' => \$n_threads,
	'i=f' => \$min_idt,
	'q=f' => \$min_qcov,
	'd=i' => \$max_dist,
	'f=s' => \$flank_len_str,
	'n=i' => \$n_iteration,
	'outer' => \$outer_flag
);

if (not(defined($base_fasta) and defined($comp_fasta))) {
	print_usage();
	exit;
}

@flank_len_list = split(',', $flank_len_str);

$int_dir_root = sprintf("%s_gap_closing_intermediates", $out_prefix);
mkdir($int_dir_root);

for $n (1..($n_iteration)) {
	$int_dir = sprintf("%s/iteration%d", $int_dir_root, $n);
	mkdir($int_dir);

	if ($n == 1) {
		$int_base_fasta = $base_fasta;
	}
	else {
		$int_base_fasta = sprintf("%s/iteration%d/final_closed.fa", $int_dir_root, $n - 1);
	}

	for $flank_len (@flank_len_list) {
		system("$sub_prog_dir/fasta_around_gap.pl $int_base_fasta ${flank_len} > $int_dir/flank${flank_len}.fa");
#		system("$minimap2_bin -c -t $n_threads -p 1 -k 19 $comp_fasta $int_dir/flank${flank_len}.fa > $int_dir/flank${flank_len}.paf 2> $int_dir/minimap2.log");
		system("$minimap2_bin -c -t $n_threads -p 1 -x asm10 $comp_fasta $int_dir/flank${flank_len}.fa > $int_dir/flank${flank_len}.paf 2> $int_dir/minimap2.log");
		system("$sub_prog_dir/paf_max_match_unique.pl $int_dir/flank${flank_len}.paf > $int_dir/flank${flank_len}_top.paf");
		system("$sub_prog_dir/paf_filter_qcov.pl $int_dir/flank${flank_len}_top.paf $min_idt $min_qcov > $int_dir/flank${flank_len}_top_filt.paf");
		system("$sub_prog_dir/paf_filter_flanking_pair.pl $int_dir/flank${flank_len}_top_filt.paf $max_dist > $int_dir/flank${flank_len}_pair.paf");

		if (not $outer_flag) {
			system("$sub_prog_dir/get_flanked_region_info.pl $int_dir/flank${flank_len}_pair.paf $comp_fasta | sort -snk2,2 | sort -sk1,1 > $int_dir/flank${flank_len}_info.tsv");
		}
		else {
			system("$sub_prog_dir/get_flanked_region_info_outer.pl $int_dir/flank${flank_len}_pair.paf $comp_fasta | sort -snk2,2 | sort -sk1,1 > $int_dir/flank${flank_len}_info.tsv");
		}

		system("perl -ane \'print if (\$F[12] !~ /[Nn]/)\' $int_dir/flank${flank_len}_info.tsv > $int_dir/flank${flank_len}_filled_info.tsv");
		system("perl -ane \'print if (\$F[12] =~ /[Nn]/)\' $int_dir/flank${flank_len}_info.tsv > $int_dir/flank${flank_len}_unfilled_info.tsv");
	}

	$cmd = "$sub_prog_dir/reduce_filled_info.pl $int_base_fasta ";
	for $flank_len (sort{$b <=> $a} @flank_len_list) {
		$cmd .= "$int_dir/flank${flank_len}_filled_info.tsv ";
	}
	$cmd .= " | sort -snk2,2 | sort -sk1,1 > $int_dir/final_filled_info.tsv";
	system($cmd);

	system("$sub_prog_dir/close_gap.pl $int_base_fasta $int_dir/final_filled_info.tsv >$int_dir/final_closed.fa");
}

$int_dir = sprintf("%s/iteration%d", $int_dir_root, $n_iteration);
system("cp $int_dir/final_closed.fa ${out_prefix}_closed.fa");



sub print_usage
{
	print

<<'USAGE';
Usage:
    combinatorial_gap_close.pl [OPTION] -b base.fa -c complement.fa

Options:
	-b FILE        fasta file of sequences to be gap-closed (mandatory)
	-c FILE        fasta file of sequences used for gap-closing (mandatory)
	-o STR         prefix of output files (default, out)
	-t INT         number of threads (default, 1)
	-i FLOAT[0-1]  minnimum sequence identity for alignment of gap-flanking regions (default, 0.90)
	-q FLOAT[0-1]  minnimum query-coverage for alignment of gap-flanking regions (default, 0.25)
	-d INT         maximum distance between aligned positions of gap-flanking region pair (default, 50000)
	-f STR         list of lengths of gap-flonking regions; specified as comma-separated string (default, '500,1000,5000,10000,20000,40000,80000,160000')
	-n STR         number of iteration (default, 2)
	-outer         flag to fill gaps using outer alignments (dfault, false)

USAGE
}
