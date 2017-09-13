# MEM
### MEM : Mendelian Error Method to rapidly detect deletions in whole exome and genome trio sequence data

#### Step 1: Call SNVs
	- Variant calling must be performed with one of the following methods
		1. Joint calling for the trio
		2. Joint genotyping (utilizing gVCFs) followed by VCF merge
		
#### Step 2: Extract Mendelian Errors (MEs) 
	- Script provided: ME_UPD_Contamination_Detection_From_TrioVCF.V.0.0.3.pl
	- Default Usage: perl ME_UPD_Contamination_Detection_From_TrioVCF.V.0.0.3.pl -i VCF -p pedigreeFile -o output_name
	- Additional Arguments:
		-i, --input		Path to input VCF [required]
		-o, --out_prefix	Output input name prefix [Default 'output']
		-p, --pedigree		pedigree file [required]
		-V, --VQSRfilter	VQSR filter to apply [Default 'yes']
		-d, --minDP		Minimum read depth required across a Trio [Default 5]
		-D, --maxDP		Maximum read depth allowed across a Trio [Default 1000]
		-c, --minAAD		Minimum alternative allele depth [Default 3]
		-a, --altBAF		Minimum B-allele frequency threshold for homozygous alternate sites i.e. 1/1 [Default 0.9]
		-R, --refBAF		Maximum B-allele frequency threshold for reference sites i.e. 0/0 [Default 0.1]
		-m, --hetBAF_min	Minimum B-allele frequency threshold for heterozygous sites i.e. 0/1 [Default 0.2]
		-H, --hetBAF_max	Maximum B-allele frequency threshold for heterozygous sites i.e. 0/1 [Default 0.8]
		-G, --minGQ		Minimum genotype quality [Default 30]
		-h, --help		Usage summary
		-v, --verbose		Detailed Usage Information

#### Step 3: Filter MEs
	- Read depth, Genotype quality and B-allele frequency filters are built into extraction script
	- Other recommended filters:
		- Exclude segmental duplication regions
		- Exclude common CNV regions
		- For WGS data
			- Mappability = 1
			- Exclude simple repeat regions
			- Remove sites that fail Hardy-Weinberg Equilibrium
			
#### Step 4: Run MEM pipeline

Note : This step assumes that bedtols V2.26.0 or higher (https://github.com/arq5x/bedtools2) is installed. It can be in the path or full path can be provided using '--tool' argument.

	- Script provided: MEM_windowAnalysis.pl
	- Default Usage: perl MEM_windowAnalysis.pl -a cases -c caseSample.list_of_files -b backgroundSample.list_of_files -w 2000000 -s 100000
	- Additional arguments:
		-a, --analysis		Analysis type [both or cases]
		-b, --bground		List of Control Files, ME files must be in bed format [Required if -a is both]
		-c, --cases		List of Cases Files [Required if run_from is not provided]
		-g, --genome		Genome file containing start and end of autosomal chromosome in b37 version [Required, provided]
		-n, --name		Name to attched to intermediate File
		-r, --run_from		Begin from intermediate steps, available options are Round1, Merge1, Round2, Merge2
		-s, --slide 		Size of slide [Required][Default 100000 ]
		-w, --window		Size of window [Required][Default 2000000 ]
		-t, --tool		Full path to bedtools [Required] [Default : bedtools]
		-M, --rescue_cases_ME 	Min ME count to rescue the window in cases after Step1 [Default 5 ]
		-X, --filter_bground_ME Max ME count in bground to filter window after Step1 [Default 2 ]	
		-Y, --filter_cases_ME 	Min ME count required in cases to consider as a ME cluster  [Default 1 ]
		-Z, --filter_bground_sample Max Sample count in bground to filter window after Step1 [Default 3 ]	
		-h, --help		Usage summary
		-v, --verbose		Detailed Usage Information
		
#### Output is a list of regions with ME clusters per case sample
	- We suggest using the first ME and the last ME in the region with the ME cluster as the minimum coordinates for the deletions 
	
	
	
	
	
	
