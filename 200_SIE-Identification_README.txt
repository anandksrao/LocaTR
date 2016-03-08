SIE Identification README file

Author: Andrew Mason
Date:	18th Feb 2016

All steps in 100_Genome-Processing_README.txt should have been completed before analysis begins.

Four structural identification programs will be used.


1) LTR Harvest (LH)

Part of the Genome Tools package, this is powerful and very quick following construction of the suffix array in Genome Processing. The three options below for -minlenltr, -maxlenltr and -similar are based on my work. These do not have to be followed exactly, and the values may have very different effects dependent on the genome. Please refer to the latest LTR Harvest manual for other options. Program should be run, then the gff3 sorted.

	(i) gt ltrharvest -index species_name_processed.fa -minlenltr 80 -maxlenltr 2000 -similar 75 -gff3 species_name_LH_out.gff3
	(ii) gt gff3 -sort species_name_LH_out.gff3 > species_name_LH_out.sorted.gff3

The gff3 file keeps information on features (mostly just the LTR pair), but for validation we need the sequences and positions using the external positions. These data can be extracted using "201_extract_LH_positions.py"

	(iii) python 201_extract_LH_positions.py species_name_LH_out.sorted.gff3 species_name_processed.fa

Positions are unstranded at this point, so strandedness will be added later. All extracted sequences come from the +ve strand.
LTR Harvest creates a lot of false positives. Positions should be thoroughly validated before being used (see 400_Validation_README.txt)


2) LTR_STRUC (LS)

It is quite a basic piece of software and has issues with memory allocation
and will not reverse complement sequences to search in both the forward and 
reverse orientations. It will also not perform batch processing. This pack,
created by Andrew Mason (andrewmason1991@gmail.com) in May 2014, automates
this process to a certain degree and produces a set of putative element
positions which can then be studied further.

This automation falls into three main parts:
1. Sequence pre-processing.
2. Executing LTR_STRUC
3. Creating a positions list

"Extras" you will need:
- win32com <Python-version> Python module 
- Windows Python - available at https://www.python.org/download/

1. Sequence pre-processing
LTR_STRUC has memory allocation problems when files are too large. A simple 
solution is to split your multi-sequence FASTA file up into chromosomes or 
smaller contigs. Creation of reverse complement sequences is also required 
for LTR_STRUC. 

	(iv) python 202_LS_seq_formatter.py [-h] [-rc REV_COMP] input_seq.fasta

As standard this splits a multi fasta into individual sequences. If there is
only one sequence it will automatically create a reverse complement. If you
have a multi fasta and need to create reverse comp sequences of all files then
you need to specify the -rc flag (just the flag). You need to do this for
LTR_STRUC to fully analyse your genome. 

All generated sequences should be moved into a directory called "sequences"
and copied to your Windows desktop ready for use in stage 2.

2. Executing LTR_STRUC
Copy the "sequences" directory from above into the LTR_STRUC directory. Copy the 
"ltrstruc_batch.py" file into the LTR_STRUC directory. Then, move this whole 
directory onto the C:\ drive so that it executes from there. The LTR_STRUC directory 
must be at the position: C:\LTR_STRUC
Navigate into this directory and run ltrstruc_batch.py. If Python 2.6 (or 
above) is on your machine double clicking this will start the batch processing. The
script will immediately ask for user input (just once, at the start) asking you to 
choose run sensitivity. This can either be a number in the range 1-10, where 1 is the
most sensitive (also slowest) and 10 is the quickest, or by entering the string "d" 
or "D" to represent "default". This is roughly equivalent to selecting 5 on the 
sensitivity scale. Press enter after you have typed your choice. If you choose an
invalid character the program will display an error message and exit.
Your machine must remain on at all times during this stage. There is no 
ability to put this process in the background.
This will run for some time (large genomes may take 3-4 days depending on your 
machine) and will produce a directory called "results". Extract this and 
return it to the LINUX server (preferably) for the final step.

WARNINGS!
i) This process runs in the foreground. If you lock your PC the current sequence
analysis will complete, but the next will not initiate. 
ii) As each new execution of LTR_STRUC_1_1.exe begins, the new window is in the
foreground. As a result you can accidentally type in the window causing the analysis
for that particular sequence to crash. This is annoying I know, but you don't have 
to restart the whole analysis. Just check your output and see which sequence has 
been skipped and rerun the python script with just that sequence as input.

3. Creating a positions list
The results file contains single sequences predicted by LTR_STRUC as putative
LTR retrotransposons. However, LTR_STRUC does not collect positional information
on these elements. The final step "ls_pos_extract.py" uses BLASTn searches to
extract positional information for each of the hits. From the LTR_STRUC analysis
we know the chromosome and the length. This data is used to validate the BLASTn
results. 
For this you need a BLAST database of the sequences you used as inputs for the
LTR_STRUC analysis. If this does not exist it is created.
The analysis is run using:

	(v) python 204_LS_ltrstruc_batch_pos_extract.py [-h] ref_genome

Following position extraction the positions are ordered and merged (in case there 
are duplicates following detection of both orientations of the sequence(s)). This 
final positions file is a tab delimited file with columns: 
chromosome-start-end-strand-source
In this file the source is "LS" for LTR_STRUC. This is helpful if you later use 
further identification methods and merge all datasets. 

An error file is created (with helpful naming) if the top BLAST hit is not congruent
with the LTR_STRUC origin and length classification. This should not happen, but 
could arise if an element has a particularly high copy number and is a recent
insertion to the genome. 


3) MGEScan_LTR (MGS)

Running MGEScan_LTR......

	(vi)
	
Positions can be extracted from the ltr.out file using a simple positions script. This also extracts the sequences.

	(vii) ./205_extract_MGS_positions.sh MGS_ltr.out species_name_processed.fa


4) RetroTector (ReTe/RT)




















