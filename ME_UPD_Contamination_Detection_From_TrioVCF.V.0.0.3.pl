#!/usr/bin/perl 
#########################################################################
#  Script     : ME_UPD_Contamination_Detection_From_TrioVCF.V2.pl
#  Author     : Nihir Patel
#  E-mail     : nihir.patel@mssm.edu
#  Date       : 08/23/2016
#  Last Edited: 06/20/2017, Nihir Patel
#  Description: Extract mendalian error given VCF input and pedigree input
#	
##########################################################################

use strict;
use Getopt::Long qw(:config posix_default bundling no_ignore_case);;
use Scalar::Util qw(looks_like_number);

my $input = '';
my $maxDP = 1000;
my $minDP = 5;
my $altBAF = 0.9;
my $refBAF = 0.1;
my $hetBAF_max = 0.8; 
my $hetBAF_min = 0.3;
my $minGQ = 30;
my $minAAD = 3;
my $type = 'SNP'; # Hidden feature , for INDEL inclusion, not recommended and not tested 
my $vqsr = 'PASS';
my $out_prefix = 'output';
my $pedigree = "none";
my $SNP_all = 'No';  # Hidden feature ,for output all SNPs
my %CHROMOSOME = (chr1 => 1, chr2 => 2, chr3 => 3, chr4 => 4, chr5 => 5, chr6 => 6, chr7 => 7, chr8 => 8, chr9 => 9, chr10 => 10, chr11 => 11, chr12 => 12, chr13 => 13, chr14 => 14, chr15 => 15, chr16 => 16, chr17 => 17, chr18 => 18, chr19 => 19, chr20 => 20, chr21 => 21, chr22 => 22);

my @data = ();

help() if (@ARGV < 1);
GetOptions('input|i=s' => \$input,
	'out_prefix|o=s' => \$out_prefix,
	'pedigree|p=s' => \$pedigree,
	'VQSRfilter|V=s' => \$vqsr,
	'minDP|d=i' => \$minDP,
	'maxDP|D=i' => \$maxDP,
	'minAAD|c=i' => \$minAAD,
	'altBAF|a=f' => \$altBAF,
	'refBAF|R=f' => \$refBAF,
	'hetBAF_min|m=f' => \$hetBAF_min,
	'hetBAF_max|H=f' => \$hetBAF_max,
	'minGQ|G=i' => \$minGQ,
	'SNP_all|S=s' => \$SNP_all,
 	'help|h' => \my $help,
 	'verbose|v' => \my $verbose);

help() if ( defined $help );
verbose() if(defined $verbose or !defined($input) or $pedigree eq 'none'); ## If input is not provided

print STDERR "Options USED: 
input(i) => $input
out_prefix(o) => $out_prefix
pedigree(p) => $pedigree
VQSRfilter(V) => $vqsr
minDP(d) => $minDP
maxDP(D) => $maxDP
min_AAD(c) => $minAAD
altBAF(a) => $altBAF
refBAF(R) => $refBAF
hetBAF_min(m) => $hetBAF_min
hetBAF_max(H) => $hetBAF_max
minGQ(G) => $minGQ
\n";

my %sample = ();
my %sample1 = ();
my %hash = ();
my %hash1 = ();
my @ped_sample = (); 
my $error = 0;
my $snp = 0;
my $chr = 1;

open(DATA, $pedigree) or die "cannot open pedigree input $pedigree\n"; 
my @data = <DATA>; chomp @data;

if ($input =~ /\.gz$/){
	open(IN, "gzip -dc $input |") || die "can't open $input\n";
}
else{
	open(IN, $input) or die "cannot open vcf input $input\n"; #vcf input
}

my $vcf ="$out_prefix".".ME.snps.txt";
open(ERR, ">$vcf"); #

if($SNP_all ne 'No'){
	my $vcf2 ="$out_prefix".".all.snps.txt";
	open(ALL, ">$vcf2"); #
}
while(<IN>){
	my $l = $_; chomp $l; 
	my @s = split(/\t/,$l);
	
	if($l =~ /#CHROM/){
		for(my $i = 9 ; $i < scalar@s ; $i++){ # Add sample info to hash
			$sample{$i} = $s[$i];
			$sample1{$s[$i]} = $i;
		}

		#### Validate Pedigree input

		my $j = 0; # It is important 
		print STDERR "#######################################################################\n";
		print STDERR "Make sure Pedigree order is correct\nFamily\tProband\tFather\tMother\n$data[0]\n";
		print STDERR "#######################################################################\n";

		for(my $i = 0 ; $i < scalar@data; $i++){
			my @D = split(/\t/,$data[$i]);
			if(exists($sample1{$D[1]}) and exists($sample1{$D[2]}) and exists($sample1{$D[3]}) and $data[$i] !~ /#/){
				@ped_sample[$j] = $data[$i];$j++;
			}
			elsif($data[$i] !~ /#/){
				print STDERR "Family $D[0] do not have all Trio sample within $input, The $D[0] will be ignored\n"; 
			}
		}
		print ERR "#CHROM\tSTART\tEND\tREF\tALT\tFilter\tBlindID\tME?\tchild_genotype\tchild_ref_count\tchild_alt_count\tchild_read_depth\tchild_GQ\tchild_BAF\tmom_genotype\tmom_ref_count\tmom_alt_count\tmom_read_depth\tmom_GQ\tmom_BAF\tdad_genotype\tdad_ref_count\tdad_alt_count\tdad_read_depth\tdad_GQ\tdad_BAF\n";
		print ALL "#CHROM\tSTART\tEND\tREF\tALT\tFilter\tBlindID\tME?\tchild_genotype\tchild_ref_count\tchild_alt_count\tchild_read_depth\tchild_GQ\tchild_BAF\tmom_genotype\tmom_ref_count\tmom_alt_count\tmom_read_depth\tmom_GQ\tmom_BAF\tdad_genotype\tdad_ref_count\tdad_alt_count\tdad_read_depth\tdad_GQ\tdad_BAF\n" if($SNP_all ne 'No');

	}
	my	$chr = "chr"."$s[0]";
	if(exists($CHROMOSOME{$s[0]}) or exists($CHROMOSOME{$chr})){
		### Initialization 
		my $idx_GT = 0;
		my $idx_AD = '';
		my $idx_DP = 0;
		my $idx_GQ = 0;
		my $idx_PL = '';
		my $child_BAF = -1;
		my $mom_BAF  = -1;
		my $dad_BAF = -1;

		##### Get Format Fields of interest
		my @FORMAT = split(/:/,$s[8]);
		for(my $j = 0 ; $j < scalar@FORMAT; $j++){
			if($FORMAT[$j] eq 'GT'){$idx_GT = $j;}	#FORMAT[0]
			if($FORMAT[$j] eq 'AD'){$idx_AD = $j;}	#FORMAT[1]
			if($FORMAT[$j] eq 'DP'){$idx_DP = $j;}	#FORMAT[2]
			if($FORMAT[$j] eq 'GQ'){$idx_GQ = $j;}	#FORMAT[3]
			if($FORMAT[$j] eq 'PL'){$idx_PL = $j;}	#FORMAT[3]
		}

		if($s[6] eq 'PASS' or lc($vqsr) ne 'yes'){
			for(my $i = 0 ; $i < scalar@ped_sample; $i++){

				my @D = split(/\t/,$ped_sample[$i]);
				my @child_FORMAT = split(/:/,$s[$sample1{$D[1]}]);
				my @mom_FORMAT = split(/:/,$s[$sample1{$D[3]}]);
				my @dad_FORMAT = split(/:/,$s[$sample1{$D[2]}]);

				## Separate Genotype inputd into 1. Genotype Quality 2.Ref RedDepth 3.Alt ReadDepth 4.Genotype Quality

				#1. Genotype
				my $child_GT = $child_FORMAT[$idx_GT];
				my $mom_GT = $mom_FORMAT[$idx_GT];
				my $dad_GT = $dad_FORMAT[$idx_GT];
				#2.Allele Depth
				my @child_AD = split(/,/,$child_FORMAT[$idx_AD]);
				my @mom_AD = split(/,/,$mom_FORMAT[$idx_AD]);
				my @dad_AD = split(/,/,$dad_FORMAT[$idx_AD]);

				### Check for multi allelic SNPs/INDEL
				my @REF = split(/,/,$s[4]);
				my $isMulti= scalar@REF;

				## IF it's not Multi Allelic check if it is an INDEL
				if(!($child_GT eq '0/0' and $dad_GT eq '0/0' and $mom_GT eq '0/0') and !($child_GT eq '0/0' and $dad_GT eq './.' and $mom_GT eq '0/0') and !($child_GT eq '0/0' and $dad_GT eq '0/0' and $mom_GT eq './.') and !($child_GT eq '0/0' and $dad_GT eq './.' and $mom_GT eq './.')){
					if((length($s[4]) < 2 and length($s[3]) < 2 ) or $type ne 'SNP'){	#if4 #### Look of simple Denovo not complex ones also exclude Indels

						## Ratio if the coverage is non Zero 
						if(($child_AD[0]+ $child_AD[1]) > 0 and ($mom_AD[0] + $mom_AD[1]) > 0 and ($dad_AD[0] + $dad_AD[1]) > 0){

							## Get BAF = B Allele Ratio
							$child_BAF =$child_AD[1]/($child_AD[0] + $child_AD[1]);
							$mom_BAF   =$mom_AD[1]/($mom_AD[0] + $mom_AD[1]) ;
							$dad_BAF   =$dad_AD[1]/($dad_AD[0] + $dad_AD[1]) ;

##### Filter based on Read Depth and Genotype Quality
							if($child_FORMAT[$idx_DP] <= $maxDP and $child_FORMAT[$idx_DP] >= $minDP and $child_FORMAT[$idx_GQ] >= $minGQ and
 							$mom_FORMAT[$idx_DP] <= $maxDP and $mom_FORMAT[$idx_DP] >= $minDP and $mom_FORMAT[$idx_GQ] >= $minGQ and
 							$dad_FORMAT[$idx_DP] <= $maxDP and $dad_FORMAT[$idx_DP] >= $minDP and $dad_FORMAT[$idx_GQ] >= $minGQ){
 
 								my $PrintStart = "$s[0]\t$s[1]\t$s[1]\t$s[3]\t$s[4]\t$s[6]\t$D[0]";
 								my $PrintEnd = "$child_FORMAT[$idx_GT]\t$child_AD[0]\t$child_AD[1]\t$child_FORMAT[$idx_DP]\t$child_FORMAT[$idx_GQ]\t$child_BAF\t$mom_FORMAT[$idx_GT]\t$mom_AD[0]\t$mom_AD[1]\t$mom_FORMAT[$idx_DP]\t$mom_FORMAT[$idx_GQ]\t$mom_BAF\t$dad_FORMAT[$idx_GT]\t$dad_AD[0]\t$dad_AD[1]\t$dad_FORMAT[$idx_DP]\t$dad_FORMAT[$idx_GQ]\t$dad_BAF"; 
 								
  								# Check If the SNP is mendel error
  								
 								## Case 1 de_novo								
								if($mom_GT eq '0/0' and $dad_GT eq '0/0' and $child_GT ne '0/0' and $child_GT ne './.' and $mom_BAF <= $refBAF and $dad_BAF <= $refBAF){
									if($child_GT eq '0/1' and $child_BAF <= $hetBAF_max and $child_BAF >= $hetBAF_min and $child_AD[1] >= $minAAD){
										print ALL "$PrintStart\tde_novo\t$PrintEnd\n" if($SNP_all ne 'No');
									}
									if($child_GT eq '1/1' and $child_BAF >= $altBAF){
										print ALL "$PrintStart\tde_novo\t$PrintEnd\n" if($SNP_all ne 'No');
									}
								}
								## Case 2 Re-checked
								elsif($mom_GT eq '0/1' and $dad_GT eq '0/0' and $child_GT eq '1/1' and $mom_BAF <= $hetBAF_max and $mom_BAF >= $hetBAF_min and $mom_AD[1] >= $minAAD and $child_BAF >= $altBAF and $dad_BAF <= $refBAF){
 										print ERR "$PrintStart\tYES\t$PrintEnd\n";
										print ALL "$PrintStart\tYES\t$PrintEnd\n" if($SNP_all ne 'No');
								}

								## Case 3 Re-checked
								elsif($mom_GT eq '0/1' and $dad_GT eq '1/1' and $child_GT eq '0/0' and $mom_BAF <= $hetBAF_max and $mom_BAF >= $hetBAF_min and $mom_AD[1] >= $minAAD and $child_BAF <= $refBAF and $dad_BAF >= $altBAF){
 										print ERR "$PrintStart\tYES\t$PrintEnd\n";
										print ALL "$PrintStart\tYES\t$PrintEnd\n" if($SNP_all ne 'No');
								}

								## Case 4 Re-checked
								elsif($mom_GT eq '1/1' and $dad_GT eq '0/1' and $child_GT eq '0/0' and $mom_BAF >= $altBAF and $dad_BAF <= $hetBAF_max and $dad_BAF >= $hetBAF_min and $dad_AD[1] >= $minAAD and $child_BAF <= $refBAF){
 										print ERR "$PrintStart\tYES\t$PrintEnd\n";
										print ALL "$PrintStart\tYES\t$PrintEnd\n" if($SNP_all ne 'No');
								}

								## Case 5 Re-checked
								elsif($mom_GT eq '0/0' and $dad_GT eq '0/1' and $child_GT eq '1/1' and $mom_BAF <= $refBAF and $dad_BAF <= $hetBAF_max and $dad_BAF >= $hetBAF_min and $dad_AD[1] >= $minAAD and $child_BAF >= $altBAF){
 										print ERR "$PrintStart\tYES\t$PrintEnd\n";
										print ALL "$PrintStart\tYES\t$PrintEnd\n" if($SNP_all ne 'No');
								}

								## Case 6 Re-checked
								elsif($mom_GT eq '0/0'and $dad_GT eq '1/1' and $child_GT ne '0/1' and $mom_BAF <= $refBAF and $dad_BAF >= $altBAF){
									if($child_GT eq '0/0' and $child_BAF <= $refBAF){
 										print ERR "$PrintStart\tYES\t$PrintEnd\n";
										print ALL "$PrintStart\tYES\t$PrintEnd\n" if($SNP_all ne 'No');
									}
									if($child_GT eq '1/1' and $child_BAF >= $altBAF){
 										print ERR "$PrintStart\tYES\t$PrintEnd\n";
										print ALL "$PrintStart\tYES\t$PrintEnd\n" if($SNP_all ne 'No');
									}
								}

								## Case 7 Re-checked
								elsif($mom_GT eq '1/1' and $dad_GT eq '0/0' and $child_GT ne '0/1' and $mom_BAF >= $altBAF and $dad_BAF <= $refBAF){
									if($child_GT eq '0/0' and $child_BAF <= $refBAF){
 										print ERR "$PrintStart\tYES\t$PrintEnd\n";
										print ALL "$PrintStart\tYES\t$PrintEnd\n" if($SNP_all ne 'No');
									}
									if($child_GT eq '1/1' and $child_BAF >= $altBAF){
 										print ERR "$PrintStart\tYES\t$PrintEnd\n";
										print ALL "$PrintStart\tYES\t$PrintEnd\n" if($SNP_all ne 'No');
									}
								}

								## Case 8 Re-checked
								elsif($mom_GT eq '1/1'	and $dad_GT eq '1/1' and $mom_BAF >= $altBAF and $dad_BAF >= $altBAF){
									if($child_GT eq '0/0' and $child_BAF <= $refBAF ){
										print ALL "$PrintStart\tContamination\t$PrintEnd\n" if($SNP_all ne 'No');
									}
									if($child_GT eq '0/1' and $child_BAF <= $hetBAF_max and $child_BAF >= $hetBAF_min and $child_AD[1] >= $minAAD){
										print ALL "$PrintStart\tContamination\t$PrintEnd\n" if($SNP_all ne 'No');
									}
									if($child_GT eq '1/1' and $child_BAF >= $altBAF){
 										print ALL "$PrintStart\tNo\t$PrintEnd\n" if($SNP_all ne 'No');
									}
								}
								# Normal SNPs
								else{
									if($child_GT eq '0/0' and $child_BAF <= $refBAF ){
	 									print ALL "$PrintStart\tNo\t$PrintEnd\n" if($SNP_all ne 'No');
									}
									if($child_GT eq '0/1' and $child_BAF <= $hetBAF_max and $child_BAF >= $hetBAF_min and $child_AD[1] >= $minAAD){
										print ALL "$PrintStart\tNo\t$PrintEnd\n" if($SNP_all ne 'No');
									}
									if($child_GT eq '1/1' and $child_BAF >= $altBAF){
										print ALL "$PrintStart\tNo\t$PrintEnd\n" if($SNP_all ne 'No');
									}
								}
							} #end if3
						} #end if2
					}
				}
			} #end for
		}
	}
}

close(IN);

sub help{
"Default Usage:\n	perl ME_UPD_Contamination_Detection_From_TrioVCF.pl -i input -p pedigreeFile -o output -V no -d 5 -D 1000 -c 3 -a 0.9 -R 0.1 -m 0.2 -H 0.8 -G 30
												or
	perl ME_UPD_Contamination_Detection_From_TrioVCF.pl --input input -p pedigreeFile --out_prefix output --VQSRfilter no --minDP 5 --maxDP 1000 --minAAD 3 --altBAF 0.9 --refBAF 0.1 --hetBAF_min 0.2 --hetBAF_max 0.8 --minGQ 30
exit();
}
 
sub verbose{
print STDERR 

"Default Usage:\n	perl ME_UPD_Contamination_Detection_From_TrioVCF.pl -i input -p pedigreeFile -o output -V no -d 5 -D 1000 -c 3 -a 0.9 -R 0.1 -m 0.2 -H 0.8 -G 30
												or
	perl ME_UPD_Contamination_Detection_From_TrioVCF.pl --input input -p pedigreeFile --out_prefix output --VQSRfilter no --minDP 5 --maxDP 1000 --minAAD 3 --altBAF 0.9 --refBAF 0.1 --hetBAF_min 0.2 --hetBAF_max 0.8 --minGQ 30

	Arguments:
		-i, --input		Input inputname [requried]
		-o, --out_prefix	Output input name prefix [default 'output']
		-p, --pedigree		pedigree file [requried]
		-V, --VQSRfilter	VQSR filter to apply [defalut 'yes']
		-d, --minDP		Minimum read depth required across a Trio[5]
		-D, --maxDP		Maximum read depth allowed across a Trio[1000]
		-c, --minAAD		Minimum alternative allele depth [3]
		-a, --altBAF		Minimum B-allele frequency threshold for alternate site i.e. 1/1 [0.9]
		-R, --refBAF		Maximum B-allele frequency threshold for refference site i.e. 0/0 [0.1]
		-m, --hetBAF_min	Minimum B-allele frequency threshold for heterozygous sites i.e. 0/1 [0.2]
		-H, --hetBAF_max	Maximum B-allele frequency threshold for heterozygous sites i.e. 0/1 [0.8]
		-G, --minGQ		Minimum genotype quality cutoff[30]
		-h, --help		Usage summary
		-v, --verbose		Detailed Usage Information\n";
exit();
}
