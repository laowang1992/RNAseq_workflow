#!/bin/env perl
# only suitable for paired-end read mapped by hisat2 log

use strict;
use warnings;

my $logfile = join(".", $ARGV[0], "log");
#print "$logfile\n";

open LOG, "$logfile";

my $total_read;
my $alg_conc_1;
my $unp_alg_1;
my $unmapped;
while(<LOG>){
	chomp;
	if(/Total pairs: (\d+)/){
		$total_read = $1 * 2;
	}
	if(/Aligned concordantly 1 time: (\d+)/){
		$alg_conc_1 = $1 * 2;
	}
	if(/Aligned 1 time: (\d+)/){
		$unp_alg_1 = $1;
	}
	if(/Aligned 0 time: (\d+)/){
		$unmapped = $1;
	}
}

my $mapped = $total_read - $unmapped;
my $unique = $alg_conc_1 + $unp_alg_1;
my $mapped_rate = $mapped / $total_read;
my $unique_rate = $unique / $total_read;

print "$ARGV[0]\t$total_read\t$mapped\t$mapped_rate\t$unique\t$unique_rate\n";
