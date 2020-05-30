#!/usr/bin/perl

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";

$sub_prog_dir = "$FindBin::Bin";
$minimap2_bin = "$sub_prog_dir/minimap2";

$n_threads = 1;
$out_str = 'rmred';

GetOptions(
	'o=s' => \$out_str,
	't=i' => \$n_threads,
);

if (@ARGV < 1) {
	print_usage();
	exit;
}

$input_list_str = join(" ", @ARGV);

$int_dir = sprintf("%s_redundancy_removal_intermediates", $out_str);
mkdir($int_dir);

system("cat $input_list_str > $int_dir/all.fa");
system("$minimap2_bin -c -t $n_threads -p 1 -k 19 $int_dir/all.fa $int_dir/all.fa > $int_dir/raw.paf 2> $int_dir/minimap2.log");
system("$sub_prog_dir/paf_filter_qcov.pl $int_dir/raw.paf 1 1 | perl -ane \'print if (\$F[0] ne \$F[5])\' > $int_dir/exact_nonself.paf");
system("$sub_prog_dir/paf_contained_short_seq_list_bubble_aware.pl $int_dir/exact_nonself.paf > $int_dir/contained.list");

for $file (@ARGV) {
	system("$sub_prog_dir/fasta_grep.pl -v -f $int_dir/contained.list $file > $file.rmred");
}


sub print_usage
{
	print

<<'USAGE';
Usage:
    remove_redundant_seq.pl [OPTION] seq1.fa [seq2.fa ...]

Options:
	-o STR  string of prefix of intermediate directory (default, rmred)
	-t INT  number of threads (default, 1)

USAGE
}
