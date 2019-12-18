#!/usr/bin/perl -w
use strict;
my $circRNA=$ARGV[0];
my $location=$ARGV[1];
my $result=$ARGV[2];
my (@cir,@site,@c,@s,$f_circRNA,@arr,$re,%hash);
open (CIR,$circRNA) or die $!;
readline CIR;
while (<CIR>){
	
	chomp($_);
	
	push (@cir,$_)

}

open (L,$location) or die $!;

while (<L>){
	chomp($_);
	push (@site,$_)

}

foreach my $line(@cir){
	
	@c=split(/\:|\||\t/,$line);
	foreach my $line2(@site){
		@s=split (/\t/,$line2);
		if ($c[0]."_1" eq $s[0]){
				#print $c[0];
				#print $s[0];
			if ($c[1]>=$s[3] && $c[1] <=$s[4]){
				push (@arr,$line);
			}
		}
		if ($c[0]."_2" eq $s[0]){
				#print $c[0];
				#print $s[0];
			if ($c[2]>=$s[3] && $c[2] <=$s[4]){
				push (@arr,$line)
				
			}
		}
					
	}

}
#print @arr;


@arr = grep { ++$hash{$_} == 2 } @arr;

if (@arr){
	foreach my $line(@arr){
	$re.=$line."\n";
	}
}
else{
$re="There has no f-circRNA in this sample!\n";
}
my $head="f-circRNA_ID\thost_gene\tf-circRNA_start\tf-circRNA_end\tjunction_reads\tSM_MS_SMS\tnon_junction_reads\tjunction_reads_ratio\tf-circRNA_type\tgene_id\tstrand\tjunction_reads_ID\n";
open OUT,">$result";
print OUT $head.$re;

