#!/bin/env perl
use strict;
use warnings;

my $sampleInfo = $ARGV[0];
my $alignPath = $ARGV[1];
my $outfile = $ARGV[2];

if(!$alignPath){
	$alignPath = "./"
}

open SAMPLES, "$sampleInfo";
open OUT, ">$outfile";

print OUT "Sample\tNumber of input reads\tUniquely mapped reads number\tUniquely mapped reads %\tNumber of reads mapped to multiple loci\t% of reads mapped to multiple loci\tNumber of reads mapped to too many loci\t% of reads mapped to too many loci\tNumber of reads unmapped: too many mismatches\t% of reads unmapped: too many mismatches\tNumber of reads unmapped: too short\t% of reads unmapped: too short\tNumber of reads unmapped: other\t% of reads unmapped: other\tNumber of chimeric reads\t% of chimeric reads\n";
while(<SAMPLES>){
	chomp;
	my @a = split(/\t/, $_);
	my $sampleName = $a[1];
	my $sampleFile = join("", $alignPath, "/", $sampleName, "/", $sampleName, "Log.final.out");
	
	my ($inputReads, $uniqMapReads, $uniqMapRate, $multiMapReads, $multiMapRate, $manyLociReads, $manyLociRate, $unmapMismatchReads, $unmapMismatchRate, $unmapShortReads, $unmapShortRate, $unmapOtherReads, $unmapOtherRate, $chimericReads, $chimericRate) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	if(-e $sampleFile){
		open ALIGNLOG, $sampleFile;
		while(<ALIGNLOG>){
			chomp;
			$inputReads         = $1 if /Number of input reads \|\t(\d+)$/;
			$uniqMapReads       = $1 if /Uniquely mapped reads number \|\t(\d+)$/;
			$uniqMapRate        = $1 if /Uniquely mapped reads % \|\t(.+%)$/;
			$multiMapReads      = $1 if /Number of reads mapped to multiple loci \|\t(\d+)$/;
			$multiMapRate       = $1 if /% of reads mapped to multiple loci \|\t(.+%)$/;
			$manyLociReads      = $1 if /Number of reads mapped to too many loci \|\t(\d+)$/;
			$manyLociRate       = $1 if /% of reads mapped to too many loci \|\t(.+%)$/;
			$unmapMismatchReads = $1 if /Number of reads unmapped: too many mismatches \|\t(\d+)$/;
			$unmapMismatchRate  = $1 if /% of reads unmapped: too many mismatches \|\t(.+%)$/;
			$unmapShortReads    = $1 if /Number of reads unmapped: too short \|\t(\d+)$/;
			$unmapShortRate     = $1 if /% of reads unmapped: too short \|\t(.+%)$/;
			$unmapOtherReads    = $1 if /Number of reads unmapped: other \|\t(\d+)$/;
			$unmapOtherRate     = $1 if /% of reads unmapped: other \|\t(.+%)$/;
			$chimericReads      = $1 if /Number of chimeric reads \|\t(\d+)$/;
			$chimericRate       = $1 if /% of chimeric reads \|\t(.+%)$/;
		}
		close ALIGNLOG;
	}else{
		print "$sampleName not find ...\n";
	}
	print OUT "$sampleName\t$inputReads\t$uniqMapReads\t$uniqMapRate\t$multiMapReads\t$multiMapRate\t$manyLociReads\t$manyLociRate\t$unmapMismatchReads\t$unmapMismatchRate\t$unmapShortReads\t$unmapShortRate\t$unmapOtherReads\t$unmapOtherRate\t$chimericReads\t$chimericRate\n";
}

close SAMPLES;
close OUT;
