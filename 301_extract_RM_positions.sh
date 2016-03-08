#!/bin/sh

## scripts takes a RepeatMasker .out file and extracts the LTR annotated positions
## ./301_extract_RM_positions.sh path/to/file.out ref_genome

outname=`echo $1 | cut -d. -f1`

tail -n +4 $1 \
| grep "LTR" \
| awk '{print $5 "\t" $6 "\t" $7 "\t" $9 "\tRM"}' \
| awk '{OFS="\t"}$4=="C"{$4="-"}1' \
| sort -k1,1 -k2,2n \
| uniq > "${outname}_RM_positions.txt"

python /groups2/avian_genomes/chicken_erv/04_avian_lineage_LTR_retrotransposon_full_analysis/00_pipeline/002_pos_merger.py --prefix "${outname}_RM" "${outname}_RM_positions.txt"
rm "${outname}_RM_positions.txt"
python /groups2/avian_genomes/chicken_erv/04_avian_lineage_LTR_retrotransposon_full_analysis/00_pipeline/001_seq_extract.py -u --prefix "${outname}_RM" "${outname}_RM_merged_positions.txt" $2



