#!/usr/bin/perl

(@ARGV != 1) and die "usage: $0 exact_match_without_self.paf\n";

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($q_name, $q_len, $t_name, $t_len) = (split(/\t/, $l))[0, 1, 5, 6];

	if ($q_len > $t_len or ($q_len == $t_len and $q_name lt $t_name)) {
		next;
	}
		
	if ($q_name =~ /primary_bubble(\d+)/) {
		$p_bubble_flag{$1} = 1;
	}
	elsif ($q_name =~ /secondary_bubble(\d+)/) {
		$s_bubble_flag{$1} = 1;
	}
}
close $in;

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($q_name, $q_len, $t_name, $t_len) = (split(/\t/, $l))[0, 1, 5, 6];

	if ($q_len > $t_len or ($q_len == $t_len and $q_name lt $t_name) or $print_flag{$q_name}) {
		next;
	}

	if ($q_name =~ /primary_bubble(\d+)/) {
		if ($s_bubble_flag{$1}) {
			print "$q_name\n";
			$print_flag{$q_name} = 1;
		}
	}
	elsif ($q_name =~ /secondary_bubble(\d+)/) {
		if ($p_bubble_flag{$1}) {
			print "$q_name\n";
			$print_flag{$q_name} = 1;
		}
	}
	else {
		print "$q_name\n";
		$print_flag{$q_name} = 1;
	}
}
close $in;
