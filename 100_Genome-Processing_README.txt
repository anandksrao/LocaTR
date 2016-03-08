Genome-Processing README file

Author: Andrew Mason
Date:	18th Feb 2016

Long sequence headers, and those including spaces, full stops (periods), underscores within names etc., are often problematic for analysis software. All genome files should therefore be preprocessed before any analysis is done.
This README file documents how to do this, and then perform other preparatory processes on the genome sequence.

1). Genome Processing

If you already have the genome file on your system, ignore the next line. 
Genomes can be obtained using the FTP sites of NCBI or Ensembl adn downloaded in zipped format (.gz) using 'wget' and the file depository address. Once the file is downloaded you are ready to process it. 

The script '101_format_genome_file.py' takes the genome file as it's single argument, in unzipped (i) or zipped (ii) format:

	(i) python 101_format_genome_file.py genome.fa
	(ii) python 101_format_genome_file.py genome.fa.gz

You can also pass the download directly to the script, but make sure to change the name from the long download name to a more usable one, such as species:

	(iii) wget fttp://link/address/file.gz ; mv file.fa.gz species_name.fa.gz; python 101_format_genome_file.py species_name.fa.gz

This will remove long sequence headers and replace with a system based on contig number. All contig numbers are filled with zeros based on total sequence number. The first sequence in a file of 100 sequences would be "seq001", in a file of 1 million sequences it would be "seq0000001". All original information is stored in a file called "species_name_information.txt", and the processed genome is called "species_name_processed.fa"


2) BLAST database creation

A usual requirement for any homology searching you perform. Some identification programs may require it too.

	(iv) makeblastdb -dbtype nucl -in species_name_processed.fa -out species_name_processed
	
This generates index files ending in ".nhr", ".nin" and ".nsq"


3) Genome Tools suffix array

Any analysis with LTR Harvest will require a suffix array for analysis. This can be created with the Genome Tools tool gt suffixerator.

	(v) gt suffixerator -db species_name_processed.fa -tis -suf -lcp -des -ssp -sds -dna


4) tRNA gene identification

Some programs use tRNA genes to annotate the PPT of LTR retrotransposons. This can be done with the program tRNAscan-SE:

	(vi) tRNAscan-SE -o species_name_tRNA_output.txt species_name_processed.fa

The output file created can then be used to extract the tRNA sequences, create alignments, and build the pHMMs required by LTR Digest and other validation scripts.

	(vii) python 102_extract_tRNA_seq.py --threads 2 species_name_processed.fa species_name_tRNA_output.txt

