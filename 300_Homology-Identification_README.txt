Homology Identification README file

Author: Andrew Mason
Date:	18th Feb 2016

All steps in 100_Genome-Processing_README.txt should have been completed before analysis begins.


Three Homology approaches will be used.


1) RepeatMasker (RM)

RepeatMasker is the standard tool for identifying repetitive elements. For this type of annotation the search effort can be reduced by ignoring low complexity regions (-nolow flag) and by specifying a taxon area, even if it is still quite vague (-species vertebrates for instance). A good standard search would be:

	(i) repeatmasker -species vertebrates -nolow species_name_processed.fa &> /dev/null &

Positions can be extracted using a simple shell script, which also extracts the sequence:

	(ii) ./301_extract_RM_positions.sh species_name_processed.fa.out species_name_processed.fa

If you have a custom library of known elements not in RepBase then this can also be used. Be careful with this, as RM does not like long header names (must be less than 50 characters including ">") and this must include a repeat class definition. For us, "#LTR" is enough. So, a custom library must be a fasta file of sequences, with the names formatted to be: ">seq_name#LTR". An example for the avian lineage based off my work on the chicken Galgal4 assembly is provided in 302_Gg4_custom_RM_lib.fas. 
The analysis can be run as so:

	(iii) repeatmasker -lib custom.fa -nolow species_name_processed.fa &> /dev/null &
	
Following both of these analyses you only need to keep the files ending .out and .tbl. Even the .tbl is not necessary, but it shows nice summary information about general repeat content etc.

Be careful with custom libraries, as any non-LTR elements in the custom library will find all of those. In mammalian genomes, the addition of a LINE would hugely increase annotated repeat content, but they are mostly going to be false positives in LTR retrotransposon identification. To ensure no non-LTR elements slip through, it is often prudent to run a second RM (using the -species vertebrate flag in this case) to check for non-LTR elements and remove them. Positions and sequences were extracted and then validated using the same script:

	(iv) ./302_extract_RM_custom_lib_positions.sh species_name_processed.fa.out species_name_processed.fa



2) RefBLASTsearch






3) ReDoSt

ReDoSt detects DIRS1 elements within genomes using a database of existing elements contained within the installation. The program is an easy install, but needs Python 2.7 with NumPY and BioPython modules, and a standalone BLAST program. The genome to be studied should be put in the folder "genomes" in its own directory (use the same name, a softlink from the species_name_processed.fa to species_name/species_name.fasta is good - must use fasta end term). The program is then run (see below) and produces a directory in the results directory with results files. The file "annotCoGenom.tmp" will be created if there are any hits. Sometimes there will not be, there has been significant lineage dependent sorting of these elements in Eukaryotes. 

	() ReDoSt folder_where_genome_is > log.txt 

Any hits will have positions, sequences, and containing 10kb genomic fragments, for the DIRS element pol gene, but not the full element. These must be gained using....[COMPLETE]. Positions can then be extracted to a GTF style file, exterior positions file and element sequences using the script below.

	() 30#_extract_ReDoSt_positions.py ####
