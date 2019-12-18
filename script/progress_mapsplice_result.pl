#!/use/bin/perl -w
use strict;
my $filename=$ARGV[0];
my $result_file=$ARGV[1];
open (FILE,$filename) or die $!;
my (@arr,@chr,$strand,@gene1,@gene2,$result,$strand1,$strand2);
while (<FILE>){
@arr=split(/\t/,$_);
@chr=split(/~/,$arr[0]);
@gene1=split(/,/,$arr[60]);
@gene2=split(/,/,$arr[61]);
$strand=$arr[5];
$strand1=substr($strand,0,1);
$strand2=substr($strand,1,2);
$result.="$gene1[0]\t$chr[0]\t$arr[1]\t$strand1\t$gene2[0]\t$chr[1]\t$arr[2]\t$strand2\n";
}

open OUT, ">$result_file.fusion.result.txt";
print OUT $result;
close OUT;
