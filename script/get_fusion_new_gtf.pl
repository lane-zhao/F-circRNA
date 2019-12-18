#!/usr/bin/perl -w
use strict;

use File::Path;
use Getopt::Long;
my ($help,$gene_location,$new_gtf_file,@fusion1,@locate_exon,@fusion,@L,@f,@transform,$fusion_file,$f,$L,$bed,$f_bed,$gene,$length,$site,@location,$line2,@locate_gene);
my (@arr,$line,@total_gtf,$gene_gtf);
Getopt::Long::GetOptions(

	"help!"        => \$help,

	"fusion_file|F:s"      => \$fusion_file,

	"gtf_annotation_file|gtf:s"     => \$gene_location,

	"output_prefix|P:s"      => \$new_gtf_file,
);

if(!defined($fusion_file) and !defined($help) and !defined($new_gtf_file)){
    
    print "Please use the --help or -H option to get usage information.\n";

}elsif(defined($help)){
					print '
Version  : 1.0
Aurthor  : Lei Zhao
Contact  : zhaolei_lane@163.com
Description: This script is used to generate the new gtf file according the fusion file.

Options  :
    -F,  --fusion_file                 <string> Fusion file generate by get_fusion_sequence.pl [ required ]
    -gtf,--gtf_annotation_file         <string> gene annotation file, usually download from Genecode. [ required ]
    -O, --output_prefix                <string> the prefix of new gtf file output. [ required ]
    -H ,--help                         output help information to screen

###################
#<fusion_file> contain 6 colums
chr1	100	200	gene1_gene2_130_350	1	+
chr1	400	900	gene1_gene2_130_350	2	-

';
}else{

open (LOCATE,$gene_location) or die $!;
while (<LOCATE>) {

	chomp($_);
	if ($_=~/^##/){
		next
		
	}else{   #####eliminate the format influence###
		$_=~s/start_codon/startcodon/g;			
		$_=~s/stop_codon/stopcodon/g;	
		$_=~s/five_prime_utr/fiveprimeutr/g;	
		$_=~s/three_prime_utr/threeprimeutr/g;	
		$_=~s/ensembl_havana/ensemblhavana/g;	
		push(@locate_exon,$_);			
	}
	
}
close(LOCATE);


###extract  the  location of fusion genes###
open (FUSION,$fusion_file) or die $!;
while (my $line=<FUSION>){
	chomp($line);
	push(@fusion,$line);
}
close FUSION;


foreach my $line(@fusion){

	@f=split(/\t|_/,$line);
	
	foreach my $line2(@locate_exon){
		#print $line2;
		@L=split(/\t/,$line2);
		if ($L[2]=~/gene/){
			next;
		}else{
		
###generate gene1 gtf file ######
			if ($L[8]=~/\"$f[3]\"/ && $f[7]=~/1/){
					if ($L[4]<$f[1]){#####fusion 位点大于该exon###
						next;
					}elsif($L[4]>$f[1]&&$L[3]<=$f[1] && $f[2]>=$L[4]){###fusion位点在exon中间并跨越该exon####
						$L[4]=$L[4]-$f[1];
						$L[3]=1;
					}elsif($L[4]>=$f[1]&&$L[3]<=$f[1] && $f[2]<=$L[4]){ ###fusion位点在exon中间并被该exon包含####
						$L[4]=$f[2]-$f[1];
						$L[3]=1;
					}elsif($L[3]>=$f[1] &&$L[4]<=$f[2]){######fusion位点在exon前面并跨越该exon####
						$L[3]=$L[3]-$f[1];
						$L[4]=$L[4]-$f[1];
						if ($L[3]=~/^0$/){$L[3]=1};
					}elsif($L[3]>=$f[1] &&$L[3]<=$f[2] && $L[4]>=$f[2]){######fusion位点在exon前面并被该exon包含####	
						$L[3]=$L[3]-$f[1];
						$L[4]=$f[2]-$f[1];
						if ($L[3]=~/^0$/){$L[3]=1};
					}else{
						next;
					}
						$length=$f[2]-$f[1];
						$gene="$f[3]_$f[4]_$f[5]_$f[6]_$f[7]\tNA\tgene\t1\t$length\n";
						$site="$f[3]_$f[4]_$f[5]_$f[6]_$f[7]\tNA\t$L[2]\t$L[3]\t$L[4]\n";
						push (@total_gtf,$site);
			}
			
###generate gene2 gtf file ###					
			if ($L[8]=~/\"$f[4]\"/ && $f[7]=~/2/){
					if ($L[4]<$f[1]){#####fusion 位点大于该exon###
						next;
					}elsif($L[4]>$f[1]&&$L[3]<=$f[1] && $f[2]>=$L[4]){###fusion位点在exon中间并跨越该exon####
						$L[4]=$L[4]-$f[1];
						$L[3]=1;
					}elsif($L[4]>=$f[1]&&$L[3]<=$f[1] && $f[2]<=$L[4]){ ###fusion位点在exon中间并被该exon包含####
						$L[4]=$f[2]-$f[1];
						$L[3]=1;
					}elsif($L[3]>=$f[1] &&$L[4]<=$f[2]){######fusion位点在exon前面并跨越该exon####
						$L[3]=$L[3]-$f[1];
						$L[4]=$L[4]-$f[1];
						if ($L[3]=~/^0$/){$L[3]=1};
					}elsif($L[3]>=$f[1] &&$L[3]<=$f[2] && $L[4]>=$f[2]){######fusion位点在exon前面并被该exon包含####	
						$L[3]=$L[3]-$f[1];
						$L[4]=$f[2]-$f[1];
						if ($L[3]=~/^0$/){$L[3]=1};
					}else{
						next;
					}
						$length=$f[2]-$f[1];
						$gene="$f[3]_$f[4]_$f[5]_$f[6]_$f[7]\tNA\tgene\t1\t$length\n";
						$site="$f[3]_$f[4]_$f[5]_$f[6]_$f[7]\tNA\t$L[2]\t$L[3]\t$L[4]";	
						push (@total_gtf,$site);
			
			}
		}
					
	}#print $gene;
#$gene_gtf.=$gene;
push (@total_gtf,$gene);
}

#####merge gtf files of gene1 and gene2 #####
#print @total_gtf;
#open (FILE,$total_gtf) or die $!;
my (@gene1,@gene2,$name,%hash,@g2,@exon,@arr,$line2,$gtf1,$gtf2);
my ($temp_gtf1,$temp_gtf2,@temp_gene1,$last_file,$temp_g2);
foreach my $line (@total_gtf){
chomp($line);
	
	if ($line=~/^(\S+)_(\S+)_(\S+)_(\S+)_1/){
		if ($line=~/^(\S+)_(\S+)_(\S+)_(\S+)_1\t(\S+)\tgene/){
		@gene1=split(/\t/,$line);
		$hash{$gene1[0]}=$gene1[4];
		push (@temp_gene1,$line);
		}else{
			push (@exon,$line); #### extract gene1  exon rows ####
		}
	}else{
		push (@gene2,$line);#####extract gene2 rows #####
	}
}


foreach my $line(@gene2){
	
	@g2=split(/_|\t/,$line);
	$name="$g2[0]_$g2[1]_$g2[2]_$g2[3]_1";
#### gene2 location add gene1 length #######
	if ($g2[6]=~/^gene$/){
		$g2[7]=$g2[7];
		$g2[8]=$g2[8]+$hash{$name};
		$temp_g2=$g2[7]+$hash{$name};
		$temp_gtf1.="$g2[0]_$g2[1]_$g2[2]_$g2[3]_2\t$g2[5]\t$g2[6]\t$temp_g2\t$g2[8]\t.\t+\t.\tgene_id \"$g2[0]_$g2[1]_$g2[2]_$g2[3].1\"\;\n";
	}else{
	
		$g2[7]=$g2[7]+$hash{$name};
		$g2[8]=$g2[8]+$hash{$name};   
	}
		$gtf1.="$g2[0]_$g2[1]_$g2[2]_$g2[3]\t$g2[5]\t$g2[6]\t$g2[7]\t$g2[8]\t.\t+\t.\tgene_id \"$g2[0]_$g2[1]_$g2[2]_$g2[3].1\"\;\n";
		
}

foreach my $line(@exon){
	chomp($line);
	@arr=split(/_|\t/,$line);
	$line2="$arr[0]_$arr[1]_$arr[2]_$arr[3]\t$arr[5]\t$arr[6]\t$arr[7]\t$arr[8]\t.\t+\t.\tgene_id \"$arr[0]_$arr[1]_$arr[2]_$arr[3].1\"\;\n";
	$gtf2.=$line2;
}
foreach my $line(@temp_gene1){
	chomp($line);
	@arr=split(/_|\t/,$line);
	$line2="$arr[0]_$arr[1]_$arr[2]_$arr[3]_1\t$arr[5]\t$arr[6]\t$arr[7]\t$arr[8]\t.\t+\t.\tgene_id \"$arr[0]_$arr[1]_$arr[2]_$arr[3].1\"\;\n";
	$temp_gtf2.=$line2;
}

open OUT ,">$new_gtf_file.gtf";
print OUT "$gtf1$gtf2";
close OUT;
###generate gene1 and gene2 location in order to identifile the cirRNA breatpoint  site in the last step of f_cirRNA software### 
#$last_file= "temp."."$new_gtf_file";
$last_file="$new_gtf_file.tmp";
open LAST , ">$last_file";
print LAST "$temp_gtf1$temp_gtf2";  
close LAST;
}