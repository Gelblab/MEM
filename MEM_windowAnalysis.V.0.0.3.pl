#!/usr/bin/perl -w

#####################################################################
##
## File Name:
##   MEM_windowAnalysis.pl
##
## Description:
##   Runs sliding window analysis to detect mendelian error clusters.
##
## Usage:
##   Default Usage: perl MEM_windowAnalysis.pl -c cases_list
##
## Author:
##   Nihir Patel (nihir.patel@mssm.edu)
##
## Date:
##   Fri Jul 28 16:00:00 EST 2017
## 
##	Changes
##	Wed Aug 16 14:52:48 EDT 2017 : Fixed Merging sub routine 
##	Wed Aug 16 17:30:25 EDT 2017 : Resolve issue where two distinct events were merged due to same values in Min,max and Mean
##  Implemented single filter for both rounds
##	Implemented Separate of Sample and create individuale VCF 
##	Added bedtools merging for split calls.
#####################################################################

use strict;
use warnings;
use IO::Zlib;
use Getopt::Long qw(:config no_ignore_case);
use Scalar::Util qw(looks_like_number);
use Carp qw/croak/;
use File::Basename;
 
my $case_input_list = '';  # list_of_cases_case_input_lists.txt
my $bground_input_list = 'none';  # list_of_bground_case_input_lists.txt
my $name = 'MEM';
my $window = 2000000;
my $slide = 100000;
my $filter_bground_ME = 2; 
my $filter_cases_ME = 1;
my $filter_bground_sample = 3; 
my $filter_rescue_cases_ME = 5; # max me in cases
my $run_from = 'Round1';
my $step_check = 'incomplete';
my $tool = 'bedtools'; 
my $genome = 'b37.genome.bed';
my $analysis = 'both';


help() if (@ARGV < 1);
GetOptions(
	'analysis|a=s' => \$analysis,
	'bground_input_list|b=s' => \$bground_input_list,
	'case_input_list|c=s' => \$case_input_list,
	'filter_bground_ME|X=i' => \$filter_bground_ME,
	'filter_cases_ME|Y=i' => \$filter_cases_ME,
	'filter_bground_sample|Z=i' => \$filter_bground_sample,
	'filter_rescue_cases_ME|M=i' => \$filter_rescue_cases_ME,
	'genome|g=s' => \$genome,
	'name|n=s' => \$name,
	'run_from|r=s' => \$run_from,
	'slide|s=i' => \$slide,
	'tool|t=s' => \$tool,
	'window|w=i' => \$window,
   	'help|h' => \my $help,
   	'verbose|v' => \my $verbose);

if($bground_input_list eq 'none'){
	$analysis = 'cases';
}


checkNexit("$genome");
my $overlap = $window * (-1);
my $windowFile_R1="$name.b37.Round1.$window.window.$slide.slide.bed";
my $windowFile_R2="$name.b37.Round2.$window.window.$slide.slide.bed";

help() if(defined $help);
verbose() if(defined $verbose or !defined($case_input_list)); ## If input is not provided

#General Running argument
sub help{
	print STDERR "Default Usage:\n	MEM_windowAnalysis.pl -c cases

Arguments:
	-a, --analysis			Analysis type [both or cases]
	-b, --bground			List of Control Files, ME files must be in bed format [Require if -a is both]
	-c, --cases			List of Cases Files [Require if run_from is not provided]
	-g, --genome			Genome file containing start and end of autosomal chromosome in b37 version [required, provided along with tool]
	-n, --name			Name to attched to intermediate File
	-r, --run_from			Round1 analysis for intermediate steps, available options are Round1, Merge1, Round2, Merge2, Split, split
	-s, --slide 			Size of slide [Required][Default 100000 ]
	-w, --window			Size of window [Required][Default 2000000 ]
	-t, --tool			Full path to bedtools [required] [Default : bedtools]
	-M, --rescue_cases_ME 		Min ME count to Rescue the window  after Step1[Default 5 ]
	-X, --filter_bground_ME 	Max ME count in bground to filter window after Step1 [Default 2 ]	
	-Y, --filter_cases_ME 		Min ME count required cases to consider as cluster  [Default 1 ]
	-Z, --filter_bground_sample 	Max Sample count in bground to filter window after Step1 [Default 3 ]	
	-h, --help			Usage summary
	-v, --verbose			Detailed Usage Information\n";
exit();
}

# Verbose information for running tool
sub verbose{
	print STDERR "Default Usage:\n	MEM_windowAnalysis.pl  -c cases -b bground -r run_from -w 2000000 -s 100000 -X 5 -Z 3 -Y 1 -M 5

Arguments:
	-a, --analysis			Analysis type [both or cases]
	-b, --bground			List of Control Files, ME files must be in bed format [Require if -a is both]
	-c, --cases			List of Cases Files [Require if run_from is not provided]
	-g, --genome			Genome file containing start and end of autosomal chromosome in b37 version [required, provided along with tool]
	-n, --name			Name to attched to intermediate File
	-r, --run_from			Round1 analysis for intermediate steps, available options are Round1, Merge1, Round2, Merge2
	-s, --slide 			Size of slide [Required][Default 100000 ]
	-w, --window			Size of window [Required][Default 2000000 ]
	-t, --tool			Full path to bedtools [required] [Default : bedtools]
	-M, --rescue_cases_ME 		Min ME count to Rescue the window[Default 5 ]
	-X, --filter_bground_ME 	Max ME count in bground to filter window after Step1 [Default 5 ]	
	-Y, --filter_min_cases_ME 		Min ME count in cases to rescue window after Step1 [Default 1 ]
	-Z, --filter_bground_sample 	Max Sample count in bground to filter window after Step1 [Default 3 ]	
	-h, --help			Usage summary
	-v, --verbose			Detailed Usage Information\n";
exit();
}

## Sub for checking files and Exist if file is missing
sub checkNexit{
	my $file = $_[0];
	unless (-e $file and !(-z $file)){print STDERR "File $file Dosen't exists or is Empty\n"; exit ();}
	unless (!(-d $file)){print STDERR "File $file is a Directory\n"; exit ();}
}

sub checkFile{
	my $file = $_[0];
	my $fname = $_[1]; 
	my $check = 'false';
	if(-z $file or (!(-e $file)) or -d $file){print STDERR "File $fname do not contains information in Autosomes  or it is empty. The file will be ignored\n"; $check = "true";}
	return $check; 
}

## Check if bedtools is available 
my $BEDTOOLS = ''; 
for my $path ( split /:/, $ENV{PATH} ){
	if (-f "$path/$tool" && -x _){
		$BEDTOOLS="$path/$tool";
	    last;
	}
	elsif  (-f "$tool" && -x _){
		$BEDTOOLS="$tool";
	    last;
	} 
}
die "bedtool is not available, please install or provide full path\n" unless ( -f $BEDTOOLS );


if($run_from ne 'round1' and $run_from ne 'Round1' and $run_from ne 'Merge1' and $run_from ne 'Round2' and $run_from ne 'Round2' and $run_from ne 'merge2' and $run_from ne 'Merge2' and $run_from ne 'split' and $run_from ne 'SPLIT' ){
	print STDERR "###########################\nError Occured !!!!!!!! \n\n"; 
	print STDERR "Please choose correct Start Point : Valid options are Round1, Merge1, Round2 and Merge2\n###########################\n\n"; 
	verbose();
}

print STDERR "######################################################################################################\nOptions USED: 
analysis(a) => $analysis
bground_input_list(b) => $bground_input_list
case_input_list(c) => $case_input_list
window(w) => $window
slide(s) => $slide
output_id(s) => $name\n

Parameters for filtering after Round 1 are as below 
bground max ME Filter(X) => $filter_bground_ME
Control max Sample filter (Z) => $filter_bground_sample
Cases min ME filter (Y) => $filter_cases_ME
Max ME in cases to rescue the window(M) => $filter_rescue_cases_ME
######################################################################################################\n";

if ($analysis eq 'both'){
	system "printf '#CHROM\tStart\tEND\tSample\tNum_Sample\tMin_ME\tMax_ME\tMean_ME\t#CHROM\tStart\tEnd\tSample_bground\tNum_Sample_bground\tMin_ME_bground\tMax_ME_bground\tMean_ME_bground\n' > .merged.header.txt";
}else{
	system "printf '#CHROM\tStart\tEND\tSample\tNum_Sample\tMin_ME\tMax_ME\tMean_ME\n' > .merged.header.txt";
}
if(lc($run_from) eq 'round1' or $run_from eq 'Round1'){

	checkNexit($case_input_list); 
	checkNexit($bground_input_list) if ($analysis eq 'both');
	
	#### Runwindow analysis ## case_input_list containing Filtered ME bed case_input_lists (ME script output VCF)  
	print STDERR "\nRuning window Analysis for samples in $case_input_list file with window of $window and slide of $slide\n";
	print STDERR "######################################################################################################\n\n";
	
	system "$BEDTOOLS makewindows -b $genome -w $window -s $slide > $windowFile_R1";
	
	checkNexit("$windowFile_R1"); 
	
	countMEPerWindow ($case_input_list,$windowFile_R1,"$name.Round1");
	checkNexit("$case_input_list.$name.Round1.data.mat");
	
	#system "sh RunwindowAnalysis.sh $case_input_list $window $slide cases $BEDTOOLS";
	
	print STDERR "Completed window Analysis for samples in $case_input_list file\n";
	print STDERR "######################################################################################################\n\n";
	
	if ($analysis eq 'both'){
		print STDERR "Running window Analysis for samples in $bground_input_list file with window of $window and slide of $slide\n";
			
		#system "sh RunwindowAnalysis.sh $bground_input_list $window $slide bground";
		print STDERR "Completed window Analysis for samples in $bground_input_list file\n";
		print STDERR "######################################################################################################\n\n";
		
		# countMEPerWindow(input_sample_list,windowFil,mid_string_for_round recognition)
		countMEPerWindow($bground_input_list,$windowFile_R1,"$name.Round1");
		checkNexit("$bground_input_list.$name.Round1.data.mat");
	}
	print STDERR "Now Checking Size of the $case_input_list.$name.Round1.data.mat \n";
	
	### Check if case_input_lists are okey 
	my $f2 = 'Ok';
	my $f1 = 'Ok';
	if ($analysis eq 'both'){$f2 = get_dims("$bground_input_list.$name.Round1.data.mat");}
	$f1 = get_dims("$case_input_list.$name.Round1.data.mat");
	

	#print STDERR "\n\n$f2\n\n"; 

	if($f1 eq 'Exit'){
		print STDERR "File $case_input_list.$name.Round1.data.mat have incorrect number of dimension please check your cases samples\n";
		print STDERR "After correcting both $case_input_list.$name.Round1.data.mat file and $bground_input_list.$name.Round1.data.mat file you can  re-run the analysis from the next step\n";
		print STDERR "perl Run_window_Analysis.pl -r mid -w $window -s $slide -p $name\n"; 
		exit;
	}
	
	elsif($f2 eq 'Exit' ){
		print STDERR "File $bground_input_list.$name.Round1.data.mat have incorrect number of dimension please check your cases samples\n"; 
		print STDERR "After correcting both $case_input_list.$name.Round1.data.mat file and $bground_input_list.$name.Round1.data.mat file you can  re-run the analysis from the next step\n";
		print STDERR "perl Run_window_Analysis.pl -r mid -w $window -s $slide -p $name\n"; 
		exit;
	}
	else{
		$step_check = 'complete'; 
	}
}	
##################################################################

if(lc($run_from) eq 'merge1' or $run_from eq 'Merge1' or $step_check eq 'complete'){
	$step_check = 'incomplete';
	#Step2 : Seprate windows for each sample (Convert to Row) BOTH for cases and bground 

	checkNexit("$case_input_list.$name.Round1.data.mat");	
	checkNexit("$bground_input_list.$name.Round1.data.mat") if ($analysis eq 'both');
	
	print STDERR "Converting matrix to row file for Cases and bground\n"; 
	convertTOrow("$case_input_list.$name.Round1.data.mat", "$case_input_list.$name.$window.$slide.row.txt");
	convertTOrow("$bground_input_list.$name.Round1.data.mat", "$bground_input_list.$name.$window.$slide.row.txt") if ($analysis eq 'both');
	
	checkNexit("$case_input_list.$name.$window.$slide.row.txt");
	checkNexit("$bground_input_list.$name.$window.$slide.row.txt") if ($analysis eq 'both');
	
	print STDERR "Done Conversion\n"; 
	print STDERR "############################################################################################\n\n";
	print STDERR "Combining files\nOverlap used is $overlap BPs\n";

	system "$BEDTOOLS merge -d $overlap -i $case_input_list.$name.$window.$slide.row.txt -c 5,5,4,4,4 -o distinct,count_distinct,min,max,mean -delim '|' > $case_input_list.$window.$slide.uniq.windows.merged.bed";
	
	if ($analysis eq 'both'){
		system "$BEDTOOLS merge -d $overlap -i $bground_input_list.$name.$window.$slide.row.txt -c 5,5,4,4,4 -o distinct,count_distinct,min,max,mean -delim '|' > $bground_input_list.$window.$slide.uniq.windows.merged.bed";
		system "$BEDTOOLS intersect -a $case_input_list.$window.$slide.uniq.windows.merged.bed -b $bground_input_list.$window.$slide.uniq.windows.merged.bed -loj -wa -wb -f 1 -r  > $name.Round1_merged.bed";
	
	}else{
		system "mv -v $case_input_list.$window.$slide.uniq.windows.merged.bed $name.Round1_merged.bed";
	}

	checkNexit("$name.Round1_merged.bed");
	checkNexit("$case_input_list.$name.$window.$slide.row.txt");
	print STDERR "Done Merging Files, output of First Round is $name.Round1_merged.bed\n"; 
	
	print STDERR "############################################################################################\n\n";
	
	if ($analysis eq 'both'){	
		## Step 4: Filter based on bground
		# If Max ME in cases <= 5 and bground has 3 or more samples with avg ME > 2 - filter out window
		# Filter out windows with Max ME = 1 in cases		
		
		print STDERR "Filtering the merged file using the bground data and provided filter\n"; 
		PostMergeCasesbgroundFilter ("$name.Round1_merged.bed", $filter_bground_ME, $filter_cases_ME, $filter_bground_sample, $filter_rescue_cases_ME, "$name.Round1_merged_filtered.bed");		
	
	}else{
		print STDERR "Removing empty rows filter\n"; 	
		PostMergeCasesOnlyFilter ("$name.Round1_merged.bed","$name.Round1_merged_filtered.bed");	
	}	
	checkNexit("$name.Round1_merged_filtered.bed");

	print STDERR "Re-merge windows .......\n"; 
	print STDERR "IF you have same number of Samples in columns 3, same number of min, max and mean ME, that means it's single window and such window will be merged\n";
	mergeFinalFile("$name.Round1_merged_filtered.bed","$name.Round1_merged_filtered_collapse.bed");

	checkNexit("$name.Round1_merged_filtered_collapse.bed");
	$step_check = 'complete'; 
}

################################################################################
######## Round 2 ############
################################################################################
	
if(lc($run_from) eq 'round2' or $run_from eq 'Round2' or $step_check eq 'complete'){	
	checkNexit("$name.Round1_merged_filtered_collapse.bed");
	$step_check = 'incomplete'; 

	checkNexit("$case_input_list.$name.updated.list");
	checkNexit("$bground_input_list.$name.updated.list") if ($analysis eq 'both');

	print STDERR "Round 1 Analysis completed, initializing Round 2 analysis \n";
	print STDERR "############################################################################################\n\n"; 

	print STDERR "Getting ME counts for merge windows\n";
	
	system "cut -f1-3 $name.Round1_merged_filtered_collapse.bed > $windowFile_R2";
	checkNexit("$windowFile_R2");
		
	print STDERR "Second Part has been sucessfully completed, continuing analysis with the Mering Script\n";  
	print STDERR "############################################################################################\n\n"; 
	print STDERR  "";

	#####system "sh Round2MergedCombined.sh $name";
	countMEPerWindow("$case_input_list.$name.updated.list","$windowFile_R2","Round2");
	countMEPerWindow("$bground_input_list.$name.updated.list","$windowFile_R2","Round2") if ($analysis eq 'both');
	

	my $f1_R2 = 'Ok'; my $f2_R2 = 'OK';
	
	$f1_R2 = get_dims("$case_input_list.$name.updated.list.Round2.data.mat");
	$f2_R2 = get_dims("$bground_input_list.$name.updated.list.Round2.data.mat") if ($analysis eq 'both');
		
	checkNexit("$case_input_list.$name.updated.list.Round2.data.mat");
	checkNexit("$bground_input_list.$name.updated.list.Round2.data.mat") if ($analysis eq 'both');
	
	if($f1_R2 eq 'Exit'){
		print STDERR "File $case_input_list.$name.updated.list.Round2.data.mat have incorrect number of dimension please check your cases samples\n";
		print STDERR "After correcting both $case_input_list.$name.updated.list.Round2.data.mat file and $bground_input_list.$name.updated.list.Round2.data.mat file you can  re-run the analysis from the next step\n";
		print STDERR "perl Run_window_Analysis.pl -r Round2 -w $window -s $slide -i $name\n"; 
		exit;
	}
	elsif($f2_R2 eq 'Exit' and $analysis eq 'both'){
		print STDERR "File $bground_input_list.$name.updated.list.Round2.data.mat have incorrect number of dimension please check your cases samples\n"; 
		print STDERR "After correcting both $case_input_list.$name.updated.list.Round2.data.mat file and $bground_input_list.$name.updated.list.Round2.data.mat file you can  re-run the analysis from the next step\n";
		print STDERR "perl Run_window_Analysis.pl -r Round2 -w $window -s $slide -i $name\n"; 
		exit;
	}
	else{
			$step_check = 'complete'; 
	}
}

if(lc($run_from) eq 'merge2' or $run_from eq 'Merge2' or $step_check = 'complete'){
	
	checkNexit("$case_input_list.$name.updated.list.Round2.data.mat");
	checkNexit("$bground_input_list.$name.updated.list.Round2.data.mat") if ($analysis eq 'both');

	convertTOrow ("$case_input_list.$name.updated.list.Round2.data.mat","$case_input_list.$name.Round2.$window.$slide.row.txt");
	convertTOrow ("$bground_input_list.$name.updated.list.Round2.data.mat","$bground_input_list.$name.Round2.$window.$slide.row.txt") if ($analysis eq 'both');

	checkNexit("$case_input_list.$name.Round2.$window.$slide.row.txt");
	checkNexit("$bground_input_list.$name.Round2.$window.$slide.row.txt") if ($analysis eq 'both');

	
	print STDERR "Done Conversion\n"; 
	print STDERR "#####################################################################\n\n";
	
	print STDERR "Combined files\nOverlap used is $overlap BPs\n";

	system "$BEDTOOLS merge -d $overlap -i $case_input_list.$name.Round2.$window.$slide.row.txt -c 5,5,4,4,4 -o distinct,count_distinct,min,max,mean -delim '|' > $case_input_list.Round2.$window.$slide.uniq.windows.merged.bed";
	 
	if ($analysis eq 'both') {
		system "$BEDTOOLS merge -d $overlap -i $bground_input_list.$name.Round2.$window.$slide.row.txt -c 5,5,4,4,4 -o distinct,count_distinct,min,max,mean -delim '|' > $bground_input_list.Round2.$window.$slide.uniq.windows.merged.bed";
		system "$BEDTOOLS intersect -a $case_input_list.Round2.$window.$slide.uniq.windows.merged.bed -b $bground_input_list.Round2.$window.$slide.uniq.windows.merged.bed -loj -wa -wb -f 1 -r  > $name.Round2_merged.bed";

	}else{
		system "mv -v $case_input_list.Round2.$window.$slide.uniq.windows.merged.bed $name.Round2_merged.bed";
	}
	checkNexit("$name.Round2_merged.bed");

	if ($analysis eq 'both'){	
		## Step 4: Filter based on bground
		# If Max ME in cases <= 5 and bground has 3 or more samples with avg ME > 2 - filter out window
		# Filter out windows with Max ME = 1 in cases		
		
		print STDERR "Filtering the merged file using the bground data and provided filter\n"; 
		PostMergeCasesbgroundFilter ("$name.Round2_merged.bed", $filter_bground_ME, $filter_cases_ME, $filter_bground_sample, $filter_rescue_cases_ME, "$name.Round2_merged_filtered.bed");		
	
	}else{
		print STDERR "Removing empty rows filter\n"; 	
		PostMergeCasesOnlyFilter ("$name.Round2_merged.bed","$name.Round2_merged_filtered.bed");	
	}	
	checkNexit("$name.Round2_merged_filtered.bed");
	
	system "cat .merged.header.txt $name.Round1_merged.bed > temp && mv temp $name.Round1_merged.bed"; 	
	system "cat .merged.header.txt $name.Round2_merged.bed > temp && mv temp $name.Round2_merged.bed"; 	
	print STDERR "Done Merging Files For Round2, output of Second Round is $name.Round2_merged.bed\nContinue further analysis and Filtering as needed\n"; 
	print STDERR "#######################################################################\n\n";
	$step_check = 'complete';

}	


if(lc($run_from) eq 'split' or $run_from eq 'SPLIT' or $step_check = 'complete'){

	SeperateSample("$name.Round2_merged_filtered.bed","$name.Round2_merged_filtered.SeperateSample.txt");
	checkNexit("$name.Round2_merged_filtered.SeperateSample.txt");
	SplitBySample("$name.Round2_merged_filtered.SeperateSample.txt",4);
	checkNexit("$name.ME.SepSample.list"); 
	
	open(MEM,"$name.ME.SepSample.list");
	while(<MEM>){
		chomp $_;
		my @s = split(/\t/,$_); 
		my $original_file = $s[0];
		my $me_file = $s[1];
		$me_file =~ s{\.[^.]+$}{}; 
		system "printf \"#CHROM\tStart\tEnd\tSample\tME_Per_Merged_Window\tNum_Window\tMin_ME\tMax_ME\tMean_ME\n\" > $me_file.merged.bed";
		system "bedtools intersect -a $me_file.bed -b $original_file -c -nonamecheck | bedtools merge -i - -c 4,9,9,9,9,9 -o distinct,collapse,count,min,max,mean -delim \"|\" >> $me_file.merged.bed";
		# FilterRound3(input,output,number of me to filer)
		FilterRound3("$me_file.merged.bed","$me_file.merged.filter.bed",2);
	}
	close (MEM); 
}

## Filter last step 

cleanFile("$case_input_list.Round2.$window.$slide.uniq.windows.merged.bed");
cleanFile("$bground_input_list.Round2.$window.$slide.uniq.windows.merged.bed");
cleanFile("$case_input_list.$name.updated.list.$name.updated.list");
cleanFile("$bground_input_list.$name.updated.list.$name.updated.list");

if ($analysis eq 'both'){
	cleanFile ("$bground_input_list.$window.$slide.uniq.windows.merged.bed");
	cleanFile ("$case_input_list.$window.$slide.uniq.windows.merged.bed");
	cleanFile ("$bground_input_list.$name.$window.$slide.row.txt");
}


sub convertTOrow {
	my $in = $_[0]; 
	my $out = $_[1]; 
	open(OUT, ">$out");
	open(IN,$in) or die "File $in can not be opened\n"; 
	my $j = 0; my @data = ();
	while(<IN>){
		chomp $_; 
		my @s = split(/\t/,$_);
		if($j ==0){
			for(my $i = 3 ; $i < scalar@s; $i++){
				$data[$i] = $s[$i];
			}
			$j++;
		} else {
			for(my $i = 3 ; $i < scalar@s; $i++){ 
	        	if($s[$i] ne '.' and $s[$i] > 0){print OUT "$s[0]\t$s[1]\t$s[2]\t$s[$i]\t$data[$i]\n";}
	        }
		}
	}
	close(IN);
	close(OUT);
}

sub mergeFinalFile {
	my $in=$_[0];
	my $out=$_[1];
	my $win=$_[2];
	
	open(IN2,$in); my @data = <IN2> ; chomp @data;
	open(OUT,">$out");
	my $str = ""; 
	my $i = 0 ;
	my $end = ""; 
	my @s2 = ();
	for(my $i = 0; $i < scalar@data; $i++){
		my $l = $data[$i];
		chomp $l; 
		my  @s = split(/\t/,$l); 
        my $ll = $data[$i-1];
       	@s2  = split(/\t/,$ll);		
		my $check = "$s[0]\t$s[3]\t$s[4]\t$s[5]\t$s[6]"; #Chrom,SampleID,Min,Max,Mean
	
		if($str ne $check or ($s2[2] < $s[1]) and $s2[0] == $s[0]){
			$str = $check;
		 	if($i == 0){
		 		print OUT "$s[0]\t$s[1]\t";
		 	}## First Line
			else{
				print OUT "$s2[2]\t$s2[3]\t$s2[4]\t$s2[5]\t$s2[6]\t$s2[7]\n";
				print OUT "$s[0]\t$s[1]\t";
			}
		}
		$end = "$s[2]\t$s[3]\t$s[4]\t$s[5]\t$s[6]\t$s[7]\n";
	}
	print OUT "$end";
	close(IN2);
	close(OUT);
}	

sub PostMergeCasesbgroundFilter {

	my $in = $_[0];
	my $filter_bground_ME = $_[1] // 2; 
	my $filter_cases_ME = $_[2] // 1;
	my $filter_bground_sample = $_[3] // 3; 
	my $filter_rescue_cases_ME = $_[4] // 5;
	my $out = $_[5] // "$name.Round1_merged_filtered.bed";
	
	open(IN,$in); 	
	open(OUT,"> $out"); 	
	while(<IN>){
		chomp $_; 
		my @s = split(/\t/,$_);
		if($s[7] > $filter_cases_ME){
			if($s[12] eq '.'){
				print OUT "$_\n";
			} 
			elsif(($s[12] >= $filter_bground_sample and $s[15] > $filter_bground_ME)){
				if($s[6] >= $filter_rescue_cases_ME){print OUT "$_\n";}
			}
		}
	}
	close(IN);
	close(OUT);

}

sub PostMergeCasesOnlyFilter {

	my $in = $_[0];
	my $out = $_[1] // "$name.Round1_merged_filtered.bed";
	
	open(IN,$in); 	
	open(OUT,"> $out");
	while(<IN>){
		chomp $_; 
		my @s = split(/\t/,$_); 
		if($s[7] > $filter_cases_ME){
			print OUT "$_\n";
		}
	}
	close(IN);
	close(OUT);
}
### Fix Window file option such that it can be used by both scripts
sub countMEPerWindow {
	my $in = $_[0];
	my $windowFile_R1 = $_[1]; 
	my $out = $_[2];
	my @outdata = (); 
	my @outdataHeader = (); # Header	 
	my $row = 0; 
	my $col = 0; 
	my %hash_pos = (); 

	open(LIST,$in);		
	open(OUTPUT, ">$in.$out.data.mat"); 
	open(OUTPUT_list, ">$in.$name.updated.list"); 
	
	checkNexit($windowFile_R1);
	open(IN,$windowFile_R1);

	my @data =(<IN>);  chomp @data; 

	while (<LIST>){
		chomp $_;
		#print STDERR "$_\n"; 
		my $file = $_;  
		
		my $check1 = checkFile($file,$file);
		## Remove other redundant chromsomes by intersecting with genome file 
		if($check1 ne 'true'){system "sed 's/chr//g' $file | $BEDTOOLS intersect -a - -b $windowFile_R1 -u >$file.temp.bed";} 

		my $check = checkFile("$file.temp.bed","$file"); # If file is empty exclude from analysis
		
		if($check ne 'true'){
			$row= 0; 		
			print OUTPUT_list "$file\n";
			print STDERR "Processing $file\n";
			$outdataHeader[$col] = "$file";	# Add header info			
			
			#system "$BEDTOOLS intersect -a $windowFile_R1 -b $file.temp.bed -wa | $BEDTOOLS groupby -c 3 -o count > $file.count";
			system "$BEDTOOLS intersect -a $windowFile_R1 -b $file.temp.bed -wa -c > $file.count";
			#cleanFile ("$file.temp.bed");
			
			my %hash = (); 
			open(IN1,"$file.count"); 
			while (<IN1>){
				chomp $_; 
				my $ll = $_; 
				my @ss = split(/\t/,$ll);
				my $pos = "$ss[0]\t$ss[1]\t$ss[2]";
				#print "$pos\n";
				$hash{$pos} = $ss[3]; ## Counts
			}
			close(IN1);
			for(my $i = 0 ;$i < scalar@data;$i++){
				my $ll = $data[$i];
				my @ss = split(/\t/,$ll); 
				my $pos = "$ss[0]\t$ss[1]\t$ss[2]";
				if($col == 0){
					$hash_pos{$i} = $pos;
				}				
				if(exists($hash{$pos})){
					$outdata[$row][$col] = $hash{$pos}; # Counts
				}else{
					$outdata[$row][$col] = '.'; # Zero MEs in window
				}
				$row++; # Number of rows 
			}
			$col++; # New File
			#cleanFile("$file.count"); 
		}else{
			#cleanFile("$file.temp.bed");
		}
	}
	
	close(LIST);
	print OUTPUT "#CHROM\tStat\tEnd\t"; 
	for(my $j = 0 ; $j < $col ;$j++){
		print OUTPUT "$outdataHeader[$j]"; 
		
		if($j < $col-1){
			print OUTPUT "\t";
		}
	}
	print OUTPUT "\n"; 
	
	for(my $i = 0 ; $i < $row ;$i++){
		for(my $j = 0; $j < $col; $j++){
		
			if ($j == 0) {
				print OUTPUT "$hash_pos{$i}\t";
			}
			
			print OUTPUT "$outdata[$i][$j]";  
			if($j < $col-1){
				print OUTPUT "\t";
			}
		}
		print OUTPUT "\n";  
	}
	close(OUTPUT);
	close(IN);
	close(OUTPUT_list);
}
## Get Dimension of resulting file to make sure everything Ran Smoothly 
sub get_dims{
	my $file = $_[0];
	my @s; 
	my $dim; 
	open(INN,$file); my $j=0; my $first=0;
	my $check = 'OK';
	while(<INN>){
		chomp $_;
  		if($first == 0) { @s=split(/\t/,$_); $dim = scalar@s; $first = 1; }
  		else{
  			@s=split(/\t/,$_); 
  			if (scalar@s != $dim) { if($check eq 'OK'){$check = 'Exit';}} 
  		}
  		$j++;
	} 
	return $check;
	close INN;
}

sub cleanFile{
	my $file = $_[0];
	if(-e $file){
		system "rm -f $file"; 
	}
}	

sub SeperateSample{
	my $file = $_[0];
	my $out = $_[1] // "$file.SeperateSample.txt";	
	open(SEP,$file);
	open(OUT,"> $out");
	 while(<SEP>) {
		chomp $_; 
		my @s = split(/\t/,$_);
		my @s1 = split(/\|/,$s[3]);
		for(my $i = 0 ; $i < scalar@s1; $i++){
			print OUT "$s[0]\t$s[1]\t$s[2]\t$s1[$i]\t$s[4]\t$s[5]\t$s[6]\t$s[7]\n";
			#if($i < scalar@s1-1){print OUT "\n"; }
		}
	}
	close(SEP);
	close(OUT); 
}

sub SplitBySample{
	my $file =$_[0]; 
	my $col  =$_[1]; 
	system "mkdir -p $name.work";

	$col = $col-1;
	open(IN,$file); 
	my %hash = (); 
	
	system "rm -f $name.work/*.MEM.sample.bed"; 
	system "rm -f $name.ME.SepSample.list"; 
	open(OUTPUT2, ">> $name.ME.SepSample.list"); 
	while(<IN>){
		chomp $_; 
		my $l = $_; 
		my @s  = split(/\t/,$l); 
		my $outFile = basename("$s[$col]");
		#print STDERR "$outFile\n";
		open(OUTPUT, ">> $name.work/$outFile.MEM.sample.bed"); 		
		print OUTPUT "$l\n"; 
		if(!(exists($hash{$outFile}))){
			print OUTPUT2 "$s[$col]\t$name.work/$outFile.MEM.sample.bed\n"; 
			$hash{$outFile} = $outFile; 
		}
	}
	close(IN);
	close(OUTPUT); 
	close(OUTPUT2); 
}

sub FilterRound3{
	my $file3 =$_[0]; 
	my $out3  =$_[1]; 
	my $filter  =$_[2]; 
	open(INR3,$file3); 
	system "rm -f $out3"; 
	open(OUTPUT3, ">> $out3"); 
	while(<INR3>){
		chomp $_; 
		my $l = $_; 
		if($l !~ /#CHROM/){
			my @s  = split(/\t/,$l); 
			if($s[8] > $filter){
				print OUTPUT3 "$l\n"; 
			}
		}
		else{
			print OUTPUT3 "$l\n";
		}			
	}	
	close(INR3);
	close(OUTPUT3); 
}
