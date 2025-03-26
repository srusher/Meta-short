#!/bin/bash
R_CONTAINER="/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/R/ggplot_dplyr/ggplot_dplyr.sif"

prefix=$1
alignment_report=$2
rscript=$3

input_tsv="$prefix-rscript-input.csv"

echo -e "Species.name,num_reads,percent_classified_reads,algorithm" > $input_tsv

awk -F'\t' '{print $1 "," $2 "," $3 "," "alignment"}' $alignment_report | sort -k3,3nr >> $input_tsv


singularity exec --bind /scicomp $R_CONTAINER Rscript $rscript $input_tsv "$(pwd)"

mv plot.png $prefix-kraken-alignment-comparison-plot.png