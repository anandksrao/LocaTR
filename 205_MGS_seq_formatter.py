### python 205_MGS_seq_formatter.py [-h] input_seq.fasta

import argparse
parser = argparse.ArgumentParser(
	description="seq_formatter is a multiple function script for splitting multi fasta files. It also produces reverse complement sequences (user option). If you have a single sequence fasta file it will detect this and automatically do a reverse complement.",
	epilog="Author: Andrew Mason; Release: 04/11/14; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("input_seq", help="Fasta file to be analysed")
usr_args = parser.parse_args()

import sys
import re
import subprocess
sys.path.append("/groups2/avian_genomes/chicken_erv/04_avian_lineage_LTR_retrotransposon_full_analysis/00_pipeline")
import static_functions

fasta_format_checker = int(subprocess.check_output("grep \"^>\" " + usr_args.input_seq + " | wc -l | awk \'{print $1}\'", shell=True))
if (fasta_format_checker == 0):
        sys.exit("Input file not in fasta format. Ensure sequences have headers starting with \">\".")
        
seq_file = open(usr_args.input_seq).read().rstrip("\n")
headers = static_functions.header_extractor(seq_file)
seq_list = static_functions.seq_only_extractor(seq_file)
static_functions.fast_fasta_splitter(headers, seq_list, "Y")
