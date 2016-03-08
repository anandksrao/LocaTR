## ./extract_MGS_positions.sh <MGS_ltr.out> <ref_genome.fa>
## process the MGEScan_LTR output file to a more usable format

#!/bin/sh

outname=`echo $2 | cut -d. -f1`

grep -v "^$" $1 \
| grep -v '^[0-9][0-9]*--' \
| sort -n \
| awk '{OFS="\t"}$2<1{$2="1"}1' \
| awk '{OFS="\t"}{end=index($1,".fa")}{$1=substr($1,0,end-1)}1' \
| sed s/_1/.1/g \
| awk '{print $1 "\t" $2 "\t" $5 "\t" $6 "\tMGS"}' > "${outname}_MGS_positions.txt"


## Notes
#
# Script reads file and excludes blank lines and 'clade lines'
# If 'start pos' values are less than 1 these are replaced with '1'
# Chromosome names are formatted to not have their file extension (will deal with .fa*)
# _1 is replaced with .1 for proper referencing
# File output in standard positions ordering
#
##

python /groups2/avian_genomes/chicken_erv/04_avian_lineage_LTR_retrotransposon_full_analysis/00_pipeline/001_seq_extract.py -u --prefix "${outname}_MGS" "${outname}_MGS_positions.txt" $2
