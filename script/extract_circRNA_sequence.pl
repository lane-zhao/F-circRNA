#!/usr/bin/perl -w
use strict;
my $filename=$ARGV[0];
my $result_name=$ARGV[1];
my $reference=$ARGV[2];
open (FILE,$filename) or die $!;
my (@arr,$line);
while (<FILE>){
	@arr=split(/\:|\||\t/,$_);
	$line.="$arr[0]\t$arr[1]\t$arr[2]\t$arr[0]\:$arr[1]\|$arr[2]\t1\t+\n";
	
	

}

open OUT, ">$result_name.fusion_circRNA.bed";
print OUT $line;
close OUT;
system "sed -i '1d' $result_name.fusion_circRNA.bed";
system "bedtools getfasta -fi $reference -bed $result_name.fusion_circRNA.bed -s -name -fullHeader -fo $result_name.fusion_circRNA_sequence.fasta";

system "rm -rf $result_name.fusion_circRNA.bed";
