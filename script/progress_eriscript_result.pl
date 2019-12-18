#!/usr/bin/perl -w
my $filename=$ARGV[0];
my $result_file=$ARGV[1];
open (FILE,$filename) or die $!;
my (@arr,$result);
readline FILE;
while (<FILE>){
	@arr=split(/\t/,$_);
	if ($arr[3]=~/^\d/ && $arr[6]=~/^\d+/){
	$result.="$arr[0]\tchr$arr[2]\t$arr[3]\t$arr[4]\t$arr[1]\tchr$arr[5]\t$arr[6]\t$arr[7]\n";
	}
	else{next
	}
}

open OUT, ">$result_file.eriscript_result.txt";
print OUT $result;
close OUT;