#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Config::IniFiles;
use Cwd qw/getcwd abs_path/;
use File::Path;
use FindBin qw($Bin $Script $RealBin);
use File::Basename qw(basename dirname);
use Parallel::ForkManager;

my $BEGIN_TIME=time();
my $script   = "$Bin/script";
my $database = "$Bin/database"; 
my $src = "$Bin/src"; 

my $workdir  = getcwd;

my ($help, $input, $outdir, $species, $show, $group_file,$fusion_method,$cir_method,$qc_method);
my ($gtf_annotation_file,$reference,$STAR_lib,$mapsplice2_reference_dir,$mapsplice2_bowtie1_index,$ericscript_db_lib);
my (%samples, %groups);

my $sample_num;
my $group_num;

&GetOptions(
	"help!"        			=> \$help,
	"input=s"     			=> \$input,
	"reference=s"			=>\$reference,
	"gtf_annotation_file=s"		=> \$gtf_annotation_file,
	"STAR_lib=s"			=> \$STAR_lib,
	"mapsplice2_reference_dir=s" 	=> \$mapsplice2_reference_dir,
	"mapsplice2_bowtie1_index=s" 	=> \$mapsplice2_bowtie1_index,
	"ericscript_db_lib=s" 		=> \$ericscript_db_lib,
	"outdir=s"     			=> \$outdir,
	"qc_method=s" 			=> \$qc_method,
	"fusion_method=s"		=> \$fusion_method,
	"cir_method=s"			=> \$cir_method,
	
);
&parsing_parameters();

&load_sample_info;
&fastq_qc();
&find_fusions();
&find_fusion_cirRNA();
print "all jobs done!\n";


###########################################################################
sub add_cmd{
	my ($cmd_array,$cmd)=@_;
	push (@$cmd_array,$cmd);
}
sub run_task {
	my ($cmd_array, $task_id, $task_name)=@_;
	#$name ||= "task"."-".$$."-".localtime();
	open W ,">run.bash/$task_id.sh";
	foreach my $cmd(@$cmd_array){
		print W $cmd,"\n";
	}
	close W;
	print STDERR scalar(localtime())." : ". "run  $task_name [start]..."."\n";
	system("sh run.bash/$task_id.sh 1>run.bash/$task_id.o 2>run.bash/$task_id.log");
	print STDERR scalar(localtime())." : ". "run  $task_name [finish]..."."\n";
}




sub parsing_parameters{

	&usage if ($help);
	
	unless(defined( $input )){
		print STDERR "Error: -input must be specified!\n";
		&usage;
	}
	
	$qc_method     ||= "trim_galore";
	$fusion_method ||= "STAR-fusion";
	$cir_method    ||= "CIRI2";
	$outdir        ||= "output";

	if($qc_method ne "trim_galore" && $qc_method ne "trimmomatic"){
		print STDERR "Error: -qc_method must be trim_galore or trimmomatic!\n";
		&usage;
	}
	
	if($fusion_method ne "STAR-fusion" && $fusion_method ne "mapsplice2" && $fusion_method ne "ericscript" ){
		print STDERR "Error: -fusion_method must be STAR-fusion, mapsplice2 or ericscript!\n";
		&usage;
	}
	if($cir_method ne "CIRI2" && $cir_method ne "find_circ"  ){
		print STDERR "Error: -cir_method must be CIRI2 or find_circ!\n";
		&usage;
	}
	unless($gtf_annotation_file && -e $gtf_annotation_file){
		print STDERR "Error: -gtf_annotation_file must be specified!\n";
		&usage;
	}
	
	$outdir = abs_path( $outdir );
	system("mkdir -p $outdir");
	system("mkdir -p run.bash");

}

sub load_sample_info{
	open  FIN,$input or die "can not open raw fq list file!";
	while(<FIN>){
		chomp;
		next if /^#/;
		my ($sn,@fqs) = split /\t/,$_;
		if(scalar @fqs == 2){
			
		}
		if(exists $samples{$sn}){
			die "sample_name [$sn] is duplicate! please check input file[$input]!";
		}
		$samples{$sn}{R1} = "$outdir/fastq_qc/$sn/$sn.R1_trimmed.fq";
		$samples{$sn}{R2} = "$outdir/fastq_qc/$sn/$sn.R2_trimmed.fq";
	}
	close FIN;
}


sub fastq_qc{
	my $pm = Parallel::ForkManager->new(10);
	open FIN,$input or die "can not open $input!";
	while(<FIN>){
		$pm->start and next;
		chomp;
		my ($sn,@fqs) = split /\t/,$_;
		
		my @cmds;
		my $task_id = "qc";
		my $task_name = "[$sn] raw data qc";
		
		if(scalar @fqs == 2){
		
		
			add_cmd(\@cmds,"mkdir -p $outdir/fastq_qc/$sn");
			#add_cmd(\@cmds,"echo -e '$sn\\t$outdir/fastq_qc/$sn/$sn.R1.fq\\t$outdir/fastq_qc/$sn/$sn.R2.fq\n' > $outdir/fastq_qc/$sn/$sn.rawdata.list");
			add_cmd(\@cmds,"echo -e '$sn\\t$outdir/fastq_qc/$sn/$sn.R1_trimmed.fq\\t$outdir/fastq_qc/$sn/$sn.R2_trimmed.fq\n' > $outdir/fastq_qc/$sn/$sn.cleandata.list");
			my @fq1 = split /,/,$fqs[0];
			my @fq1_tmp;
			foreach my $x(@fq1){
				my $x_name = basename($x);
				if($x =~ /.gz$/){
					add_cmd(\@cmds,"zcat $x > $outdir/fastq_qc/$sn/tmp.$x_name.R1.fq");
				}else{
					add_cmd(\@cmds,"cat  $x > $outdir/fastq_qc/$sn/tmp.$x_name.R1.fq");
				}
				push @fq1_tmp,"$outdir/fastq_qc/$sn/tmp.$x_name.R1.fq";
			}
			add_cmd(\@cmds, "cat @fq1_tmp > $outdir/fastq_qc/$sn/$sn.R1.fq");
			
			my @fq2 = split /,/,$fqs[1];
			my @fq2_tmp;
			foreach my $x(@fq2){
				my $x_name = basename($x);
				if($x =~ /.gz$/){
					add_cmd(\@cmds,"zcat $x > $outdir/fastq_qc/$sn/tmp.$x_name.R2.fq");
				}else{
					add_cmd(\@cmds,"cat  $x > $outdir/fastq_qc/$sn/tmp.$x_name.R2.fq");
				}
				push @fq2_tmp,"$outdir/fastq_qc/$sn/tmp.$x_name.R2.fq";
			}
			add_cmd(\@cmds,"cat @fq2_tmp > $outdir/fastq_qc/$sn/$sn.R2.fq");
			
			add_cmd(\@cmds,"rm @fq1_tmp @fq2_tmp");
			if($qc_method eq "trimmomatic"){
				add_cmd(\@cmds,"java -Xms20g -Xmx20g -jar $src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 20 -phred33 $outdir/fastq_qc/$sn/$sn.R1.fq $outdir/fastq_qc/$sn/$sn.R2.fq $outdir/fastq_qc/$sn/$sn.R1_trimmed.fq $outdir/fastq_qc/$sn/$sn.R1_unpaired.fq $outdir/fastq_qc/$sn/$sn.R2_trimmed.fq $outdir/fastq_qc/$sn/$sn.R2_unpaired.fq ILLUMINACLIP:/home/pub/software/Trimmomatic/0.38/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:75");
				add_cmd(\@cmds,"java -jar $script/FastqStat.jar -i $outdir/fastq_qc/$sn/$sn.rawdata.list > $outdir/fastq_qc/$sn/$sn.rawdata.stat.main.xls");
				add_cmd(\@cmds,"java -jar $script/FastqStat.jar -i $outdir/fastq_qc/$sn/$sn.cleandata.list > $outdir/fastq_qc/$sn/$sn.cleandata.stat.main.xls");
				#add_cmd(\@cmds,"fastqc -t 10 -o $outdir/fastq_qc/$sn $outdir/fastq_qc/$sn/$sn.R1.fq &");
				#add_cmd(\@cmds,"fastqc -t 10 -o $outdir/fastq_qc/$sn $outdir/fastq_qc/$sn/$sn.R2.fq &");
				add_cmd(\@cmds,"fastqc -t 10 -o $outdir/fastq_qc/$sn $outdir/fastq_qc/$sn/$sn.R1_trimmed.fq &");
				add_cmd(\@cmds,"fastqc -t 10 -o $outdir/fastq_qc/$sn $outdir/fastq_qc/$sn/$sn.R2_trimmed.fq &");
				add_cmd(\@cmds,"wait");
				add_cmd(\@cmds,"rm $outdir/fastq_qc/$sn/$sn.R1.fq $outdir/fastq_qc/$sn/$sn.R2.fq ");
			}elsif($qc_method eq "trim_galore"){
				
				#add_cmd("cd $outdir/fastq_qc/$sn");
			
				add_cmd(\@cmds,"$src/trim_galore -q 20 --phred33 --fastqc --stringency 1 -e 0.1 --dont_gzip -o $outdir/fastq_qc/$sn/ --paired $outdir/fastq_qc/$sn/$sn.R1.fq $outdir/fastq_qc/$sn/$sn.R2.fq");
				
				#add_cmd("cd -");
				
				#chdir "";
				add_cmd(\@cmds,"wait");
				
				add_cmd(\@cmds,"mv $outdir/fastq_qc/$sn/$sn.R1_val_1.fq $outdir/fastq_qc/$sn/$sn.R1_trimmed.fq");
				add_cmd(\@cmds,"mv $outdir/fastq_qc/$sn/$sn.R1_val_1_fastqc.zip  $outdir/fastq_qc/$sn/$sn.R1_trimmed_fastqc.zip");
				add_cmd(\@cmds,"mv $outdir/fastq_qc/$sn/$sn.R1_val_1_fastqc.html $outdir/fastq_qc/$sn/$sn.R1_trimmed_fastqc.html");
				add_cmd(\@cmds,"mv $outdir/fastq_qc/$sn/$sn.R2_val_2.fq $outdir/fastq_qc/$sn/$sn.R2_trimmed.fq");
				add_cmd(\@cmds,"mv $outdir/fastq_qc/$sn/$sn.R2_val_2_fastqc.zip  $outdir/fastq_qc/$sn/$sn.R2_trimmed_fastqc.zip");
				add_cmd(\@cmds,"mv $outdir/fastq_qc/$sn/$sn.R2_val_2_fastqc.html $outdir/fastq_qc/$sn/$sn.R2_trimmed_fastqc.html");
				#add_cmd(\@cmds,"fastqc -t 10 -o $outdir/fastq_qc/$sn $outdir/fastq_qc/$sn/$sn.R1.fq &");
				#add_cmd(\@cmds,"fastqc -t 10 -o $outdir/fastq_qc/$sn $outdir/fastq_qc/$sn/$sn.R2.fq &");
				add_cmd(\@cmds,"rm $outdir/fastq_qc/$sn/$sn.R1.fq $outdir/fastq_qc/$sn/$sn.R2.fq ");
				
			}else{
				die "only support trim_galore and trimmomatic now !";
			}
			
			
		}
		run_task(\@cmds,$sn.".".$task_id,$task_name);
		$pm->finish;
	}
	$pm->wait_all_children;
	
}


sub find_fusions{
		my @cmds;
		my $pm = Parallel::ForkManager->new(10);
		
	foreach my $sn(sort keys %samples){
		my $task_id = "find_fusions";
		my $task_name = "[$sn] find fusions";
		$pm->start and next;
		
		if($fusion_method eq "STAR-fusion"){
				add_cmd(\@cmds,"mkdir -p $outdir/fusion_results/$sn");
				add_cmd(\@cmds,"STAR-Fusion  --no_annotation_filter  --min_FFPM 0 --genome_lib_dir $STAR_lib --left_fq $outdir/fastq_qc/$sn/$sn.R1_trimmed.fq --right_fq $outdir/fastq_qc/$sn/$sn.R2_trimmed.fq --output_dir $outdir/fusion_results/$sn");
				add_cmd(\@cmds,"cd $outdir/fusion_results/$sn");
				add_cmd(\@cmds,"perl $script/progress_star_fusion_result.pl star-fusion.fusion_predictions.tsv $sn");
				add_cmd(\@cmds,"perl $script/get_fusion_sequence.pl --fusion_file $sn.star-fusion_result.txt --reference $reference --gtf_annotation_file $gtf_annotation_file --out_fuison_location $sn.location.txt --output_prefix $sn.fusion");
				add_cmd(\@cmds,"perl $script/get_fusion_new_gtf.pl --fusion_file $sn.location.txt --gtf_annotation_file $gtf_annotation_file --output_prefix $sn.fusion");
				add_cmd(\@cmds,"cd $outdir");
			}elsif($fusion_method eq "mapsplice2"){
				add_cmd(\@cmds,"mkdir -p $outdir/fusion_results/$sn");
				add_cmd(\@cmds,"python $src/MapSplice-v2.2.1/mapsplice.py --fusion-non-canonical -c $mapsplice2_reference_dir -x $mapsplice2_bowtie1_index -1 $outdir/fastq_qc/$sn/$sn.R1_trimmed.fq -2 $outdir/fastq_qc/$sn/$sn.R2_trimmed.fq -p 10 --gene-gtf  $gtf_annotation_file -o $outdir/fusion_results/$sn");
				add_cmd(\@cmds,"cd $outdir/fusion_results/$sn");
				add_cmd(\@cmds,"perl $script/progress_mapsplice_result.pl fusions_well_annotated.txt $sn");
				add_cmd(\@cmds,"perl $script/get_fusion_sequence.pl --fusion_file $sn.fusion.result.txt --reference $reference --gtf_annotation_file $gtf_annotation_file --out_fuison_location $sn.location.txt --output_prefix $sn.fusion");
				add_cmd(\@cmds,"perl $script/get_fusion_new_gtf.pl --fusion_file $sn.location.txt --gtf_annotation_file $gtf_annotation_file --output_prefix $sn.fusion");
				add_cmd(\@cmds,"cd $outdir");
			}elsif($fusion_method eq "ericscript"){
				add_cmd(\@cmds,"mkdir -p $outdir/fusion_results");
				add_cmd(\@cmds,"perl $src/ericscript-0.5.5/ericscript.pl  -p 10 -name $sn --dbfolder $ericscript_db_lib -o $outdir/fusion_results/$sn --readlength 150 $outdir/fastq_qc/$sn/$sn.R1_trimmed.fq $outdir/fastq_qc/$sn/$sn.R2_trimmed.fq");
				add_cmd(\@cmds,"cd $outdir/fusion_results/$sn");
				add_cmd(\@cmds,"perl $script/progress_eriscript_result.pl $sn.results.filtered.tsv $sn");
				add_cmd(\@cmds,"perl $script/get_fusion_sequence.pl --fusion_file $sn.eriscript_result.txt --reference $reference --gtf_annotation_file $gtf_annotation_file --out_fuison_location $sn.location.txt --output_prefix $sn.fusion");
				add_cmd(\@cmds,"perl $script/get_fusion_new_gtf.pl --fusion_file $sn.location.txt --gtf_annotation_file $gtf_annotation_file --output_prefix $sn.fusion");
				add_cmd(\@cmds,"cd $outdir");
			}else{
				die "only support STAR-Fusion,mapsplice2 and ericscript now !";
			}
		
		
		run_task(\@cmds,$sn.".".$task_id,$task_name);
		$pm->finish;
	}
	$pm->wait_all_children;
}

sub find_fusion_cirRNA{
		my $pm = Parallel::ForkManager->new(10);
		my @cmds;
		
		
		foreach my $sn(sort keys %samples){
		my $task_id = "find_fusions_circRNA";
		my $task_name = "[$sn] find fusion circRNA";
		$pm->start and next;
		add_cmd(\@cmds,"mkdir -p $outdir/cirRNA_results");
			if($cir_method eq "CIRI2"){
				add_cmd(\@cmds,"mkdir -p $outdir/new_reference_alignment/");
				add_cmd(\@cmds,"cd $outdir/fusion_results/$sn");
				add_cmd(\@cmds,"bwa index $sn.fusion.fasta");
				add_cmd(\@cmds,"bwa mem -M -t 20 $sn.fusion.fasta $outdir/fastq_qc/$sn/$sn.R1_trimmed.fq $outdir/fastq_qc/$sn/$sn.R2_trimmed.fq -o $outdir/new_reference_alignment/$sn.sam 1>$outdir/new_reference_alignment/$sn.stdout 2> $outdir/new_reference_alignment/$sn.stderr");
				add_cmd(\@cmds,"cd $outdir/cirRNA_results/");
				add_cmd(\@cmds,"perl $src/CIRI2.pl -I $outdir/new_reference_alignment/$sn.sam -O $outdir/cirRNA_results/$sn.ciri.txt -0 -F $outdir/fusion_results/$sn/$sn.fusion.fasta -T 10 -A $outdir/fusion_results/$sn/$sn.fusion.gtf -G $outdir/cirRNA_results/$sn.ciri.log");
				add_cmd(\@cmds,"perl $script/identify_f_circRNA_in_circRNA_CIRI2.pl $outdir/cirRNA_results/$sn.ciri.txt $outdir/fusion_results/$sn/$sn.fusion.tmp  $outdir/cirRNA_results/$sn.fusion_circRNA.txt");
				add_cmd(\@cmds,"perl $script/extract_circRNA_sequence.pl $outdir/cirRNA_results/$sn.fusion_circRNA.txt $sn  $outdir/fusion_results/$sn/$sn.fusion.fasta");
			}elsif($cir_method eq "find_circ"){
				add_cmd(\@cmds,"mkdir -p $outdir/new_reference_alignment/");
				add_cmd(\@cmds,"cd $outdir/fusion_results/$sn");
				add_cmd(\@cmds,"bowtie2-build --threads 10 $outdir/fusion_results/$sn/$sn.fusion.fasta $sn.fusion.fasta.bowtie2 ");
				add_cmd(\@cmds,"bowtie2 -p10 --very-sensitive --mm  --score-min=C,-15,0 -x $sn.fusion.fasta.bowtie2  -1 $outdir/fastq_qc/$sn/$sn.R1_trimmed.fq -2 $outdir/fastq_qc/$sn/$sn.R2_trimmed.fq -S $outdir/new_reference_alignment/$sn.sam 2> 1.log ");
				add_cmd(\@cmds,"cd $outdir/");
				
				add_cmd(\@cmds,"$src/samtools-1.6/bin/samtools view $outdir/new_reference_alignment/$sn.sam -hbuS | $src/samtools-1.6/bin/samtools sort -o $outdir/new_reference_alignment/$sn.bam");
				add_cmd(\@cmds,"$src/samtools-1.6/bin/samtools view -hf 4 $outdir/new_reference_alignment/$sn.bam | $src/samtools-1.6/bin/samtools view -Sb -o  $outdir/new_reference_alignment/$sn.unmapped.bam");
				add_cmd(\@cmds,"$src/unmapped2anchors.py $outdir/new_reference_alignment/$sn.unmapped.bam > $outdir/new_reference_alignment/$sn.qfa");
				
		
				add_cmd(\@cmds,"bowtie2 -p 10 --reorder --mm  --score-min=C,-15,0 -q -x $outdir/fusion_results/$sn/$sn.fusion.fasta.bowtie2 -U $outdir/new_reference_alignment/$sn.qfa -S  $outdir/new_reference_alignment/$sn.align.sam 2> bt2_secondpass.log");
	
				add_cmd(\@cmds,"cat $outdir/new_reference_alignment/$sn.align.sam | $src/find_circ.py  -G $outdir/fusion_results/$sn/$sn.fusion.fasta -p f_cir. -s $outdir/cirRNA_results/sites.log > $outdir/cirRNA_results/sites.bed 2> $outdir/cirRNA_results/sites.reads");
				add_cmd(\@cmds,"cd $outdir/cirRNA_results");
				add_cmd(\@cmds,"grep CIRCULAR sites.bed  > $sn.circ_candidates.bed");
				add_cmd(\@cmds,"perl $script/identify_f_circRNA_in_circRNA_find_circ.pl $outdir/cirRNA_results/$sn.circ_candidates.bed $outdir/fusion_results/$sn/$sn.fusion.tmp  $outdir/cirRNA_results/$sn.fusion_circRNA.txt");
				add_cmd(\@cmds,"perl $script/extract_circRNA_sequence.pl $outdir/cirRNA_results/$sn.fusion_circRNA.txt $sn  $outdir/fusion_results/$sn/$sn.fusion.fasta");
				
				
			}else{
				die "only support CIRI2 and find_circ now !";
			}
			run_task(\@cmds,$sn.".".$task_id,$task_name);
			$pm->finish;
		}
		$pm->wait_all_children;
}
#foreach my $line (@cmd){print $line}
sub usage{
	die(qq/
Usage:	$0 <arguments>
Arguments:
  Input Options:
	-input          		<string>	input fastq list, users need to prepare a file to tell f-circRNA the name and path of fastq, which separate with tab, required
							example: 	Test1	test1_R1.fastq	test2_R2.fastq
									Test2	test2_R1.fastq	test2_R2.fastq

	-gtf_annotation_file 		<string>	gene annotation file, required

  Output Options:
        -outdir         		<string>	output dir, default is output.
	
  Analysis Options:
	-qc_method			<string>	qc method, only support "trimmomatic" or "trim_galore", default is "trim_galore"
	-fusion_method			<string>	fusion method，only support "STAR-fusion", "mapsplice2" or "ericscript", default is "STAR-fusion"
	-cir_method			<string>	cirRNA identify method, only support "CIRI2" or "find_circ", default is "CIRI2"
	-STAR_lib			<string>	STAR-fusion lib file，if fusion method is "STAR-fusion", this file is  required.
	-mapsplice2_reference_dir	<string>	Reference file, users must prepare it as mapsplice2 required format,if fusion method is "mapslice2", this file is required
	-mapsplice2_bowtie1_index	<string>	Reference index file, bowtie1 is required to budild this index, if fusion method is "mapslice2", this file is required
	-ericscript_db_lib		<string>	ericscript lib file, if fusion method is "ericscript", this file is required
\n/);
}
