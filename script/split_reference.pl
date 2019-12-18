#!/usr/bin/env perl
use strict;
use warnings;
my $reference=$ARGV[0];
system("mkdir -p reference_dir");
my %hash;
my $id = "";
open FIN,$reference;
while(my $line = <FIN>){
    if($line =~ /^>(.+?)\s/){
        $id = $1;
    }else{
        $hash{$id} .= $line;
    }
}
close FOUT;

foreach my $id(sort keys %hash){
    open FOUT,">reference_dir/$id.fa";
    print FOUT ">$id\n";
    print FOUT $hash{$id};
    close FOUT;
}

