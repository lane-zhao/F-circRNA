#!/usr/bin/perl -w
#this script is  used to progress the star-fusion result to an uniq format ####
use strict;
my $fusion_result=$ARGV[0];
my $result_file=$ARGV[1];
open (FILE,$fusion_result) or die $!;
readline FILE;
my (@arr,$result);
while (<FILE>){
@arr=split(/\t|--|:/,$_);
$result.="$arr[0]\t$arr[6]\t$arr[7]\t$arr[8]\t$arr[1]\t$arr[10]\t$arr[11]\t$arr[12]\n";
}

open OUT, ">$result_file.star-fusion_result.txt";
print OUT $result;
close OUT;
