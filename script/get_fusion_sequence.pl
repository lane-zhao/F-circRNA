#!/usr/bin/perl -w
use strict;
use File::Path;
use Getopt::Long;
my ($fusion_file,$gene_location,$reference,$result_file,$software_name,$help,@locate,@fusion,@L,@f,$f,$L,$bed,$f_bed,@arr,@location,$line,$fusion_fasta);
my (@filter,$filter_name,@filter_arr,@filter_samples,$filter_result);
Getopt::Long::GetOptions(

	"help!"        => \$help,

	"fusion_file|F:s"      => \$fusion_file,

	"gtf_annotation_file|gtf:s"     => \$gene_location,
	"reference|R:s"				=>\$reference,
	"out_fuison_location|L:s"    => \$result_file,

	"output_prefix|P:s"      => \$software_name,
);

if(!defined($fusion_file) and !defined($help) and !defined($gene_location) and !defined($result_file) and !defined($software_name)and !defined($reference)){
    
    print "Please use the --help or -H option to get usage information.\n";

}elsif(defined($help)){
					print '
Version  : 1.0
Aurthor  : Lei Zhao
Contact  : zhaolei_lane@163.com

Options  :
    -F,  --fusion_file                 <string> fusion location file of genes,extract from the fusion software directly. [ required ]
    -gtf,--gtf_annotation_file         <string> gene annotation file, usually download from Genecode. [ required ]
    -R,--reference					<string>	reference file, usually download from Genecode. Reference index file build by samtools "faidx" command also needed. [ required ]
	-L, --out_fuison_location     <string> the bed file we will generate,this file will be used in the script "get_fusion_new_gtf.pl" . [ required]
    -O, --output_prefix                <string> the prefix of fasta output. [ required ]
    -H ,--help                         output help information to screen

###################
#<fusion_file> contain 8 colums
gene1	chr1	100	+	gene2	chr12	123	-
gene3	chr3	134	-	gene4	chr7	111	+

';
}else{

open (LO,$gene_location) or die $!;


###progress the gencode-gtf files ##
while (<LO>) {
 	if ($_=~/^#/){
		next
	}else{								
			@arr=split(/\t/,$_);
			if ($arr[2]=~/gene/){
			$_=~s/\sgene_name\s"(\S+)"/$1/;
			@location=split(/\t|;/,$_);
			}else{
			next
			}
	}
	$line="$location[0]\t$location[3]\t$location[4]\t$location[6]\t$location[10]\n"; #####bed format####
	chomp($line);
	push(@locate,$line); 
} 
close(LO);

###extract  the  location of fusion genes###
open (FUSION,$fusion_file) or die $!;
while (my $line=<FUSION>){
	chomp($line);
	push(@fusion,$line);
}
close FUSION;

foreach my $line(@locate){
	@L=split(/\t/,$line);

	foreach my $line2(@fusion){
		@f=split(/\t/,$line2);
		if ($f[0]=~/^$L[4]$/){
			if($L[3]=~/\+/){
				if ($L[1]>$f[2]){
			
				next;
				}
				else{
				$bed="$L[0]\t$L[1]\t$f[2]\t$f[0]_$f[4]_$f[2]_$f[6]_1\t1\t+\n";
				}
			}
			else{
					if ($f[2]>$L[2]){   ####eliminate the event that fusion location out of the gene location caused by some little change of gtf release####
					next;
					}
					else{
					$bed="$L[0]\t$f[2]\t$L[2]\t$f[0]_$f[4]_$f[2]_$f[6]_1\t1\t-\n";
					}
					
				
			}
			
		$f_bed.=$bed
		}
		
			if ($f[4]=~/^$L[4]$/){
				if($L[3]=~/\+/){
					if ($f[6]>$L[2]){
					next
					}
					else{
			
					$bed="$L[0]\t$f[6]\t$L[2]\t$f[0]_$f[4]_$f[2]_$f[6]_2\t1\t+\n";
					}
				}
				else{
					if ($L[1]>$f[6]){   ####eliminate the event that fusion location out of the gene location caused by some little change of gtf release####
					next;
					}
					else{
					$bed="$L[0]\t$L[1]\t$f[6]\t$f[0]_$f[4]_$f[2]_$f[6]_2\t1\t-\n";
					}
				}
			$f_bed.= $bed;
			}
	}

}

open OUT ,">$result_file.tmp"; ######fusion bed format of  every gene####

print OUT $f_bed;

close OUT;

######remove the  gene name not match the genesymbol of gtf file caused by version variation###### 
open (FILTER,$result_file.".tmp");
my %filter_hash;
while (<FILTER>){
@filter=split(/\t/,$_);
$filter_hash{$filter[3]}=$_;

}

foreach my $keys (keys%filter_hash){
	@filter_arr=split(/_/,$keys);
	$filter_name=$filter_arr[0]."_".$filter_arr[1]."_".$filter_arr[2]."_".$filter_arr[3];
	my %filter_uniq;
	push (@filter_samples,$filter_name);
	@filter_samples=grep ++$filter_uniq{$_}<2,@filter_samples;
}
foreach my $filter_name(@filter_samples){
	if (exists $filter_hash{$filter_name."_1"} && $filter_hash{$filter_name."_2"} ){
		$filter_result.=$filter_hash{$filter_name."_1"}.$filter_hash{$filter_name."_2"};
	}
	else{
		next;
	}
}
open FILTER ,">$result_file"; ######fusion bed format of  every gene####

print FILTER $filter_result;
close FILTER;

system "rm -rf $result_file.tmp" ;

####extrcat sequence of fusion genes####
system "bedtools getfasta -fi $reference -bed $result_file -s -name -fullHeader -fo $result_file.tmp.fasta";

#######merge sequnece of two fusion genes##############
my $filename=$result_file.".tmp.fasta";
open (FASTA,$filename) or die $!;
my (%hash,@line,$name,$temp,$key,@samples);
while (<FASTA>){
	chomp($_);
	if ($_=~/^>/){
	$key=$_;
	}else{
		 $hash{$key}=$_;
		 
	}

	
}


foreach my $keys(keys%hash){
	
	@line=split(/_/,$keys);
	$name=$line[0]."_".$line[1]."_".$line[2]."_".$line[3];
	my %uniq;
	push (@samples,$name);
	@samples=grep ++$uniq{$_}<2,@samples;
		
}

foreach my $name(@samples){
	if (exists $hash{$name."_1"} &&$hash{$name."_2"} ){
		$hash{$name."_2"}=lc($hash{$name."_2"});
		$fusion_fasta.=$name."\n".$hash{$name."_1"}.$hash{$name."_2"}."\n";
	}
	else{
		next
	}
}
open FA, ">$software_name.fasta";
print FA $fusion_fasta;
close FA;

}
