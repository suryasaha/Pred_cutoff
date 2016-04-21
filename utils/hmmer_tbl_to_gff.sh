#!/bin/sh

# Surya Saha
# Solgenomics@BTI
# Purpose: Create GFF from hmmer-3.1b2 TBL files

usage(){
	echo "usage:
	$0 <input .TBL file> <output .GFF file>"
	exit 1
}

grep -v '^#' PtoDC3000-HrpL.hmm.tbl|  awk -v OFS='\t' '{print $1,"nhmmer","match_region",$9,$10,$13,$12,".","Score="$14}' > PS_14_PtoDC3000-HrpL.hmm.gff

if [ "$#" -ne 2 ]
then
	usage
fi

printf "Input .TBL file : %s \n" "$1"
printf "Input .TBL file : %s \n" "$2"

printf "##gff3\n##Original lines from hmmer TBL file in comments\n" > "$2"

while read LINE
do
	echo $LINE| grep '^#'
	echo $LINE| awk '$12=="+"' | awk -v OFS='\t' '{print $1,"nhmmer","match_region",$9,$10,$13,$12,".","Score="$14}'
	echo $LINE| awk '$12=="-"' | awk -v OFS='\t' '{print $1,"nhmmer","match_region",$10,$9,$13,$12,".","Score="$14}'

done < "$1" >> "$2"
