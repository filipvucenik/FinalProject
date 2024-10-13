#!/bin/bash

# Step 1: get input files
# 1 first bed
# 2 second bed that will be substracted from the first bed
# 3 fasta reads file
# 4 bam file for igv track

a=$1
b=$2
reads=$3
bam_track=$4

if [ -z "$a" ] || [ -z "$b" ] || [ -z "$reads" ] || [ -z "$bam_track" ]; then
    echo "Please provide the following arguments: 1) first bed file, 2) second bed file, 3) fasta reads file, 4) bam file for igv track (if you want to add more tracks, write them separated by space in the same string)"
    exit 1
fi

echo "Current Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Step 2: run bedtools to subtract the second bed from the first bed
IFS='/' read -r -a arr1 <<< "$a"
IFS='/' read -r -a arr2 <<< "$b"
rmpath1=${arr1[-1]}
rmpath2=${arr2[-1]}
IFS='.' read -r -a ar1 <<< "$rmpath1"
IFS='.' read -r -a ar2 <<< "$rmpath2"
sub_bed="${ar1[0]}_subtracted_${ar2[0]}.bed"
bedtools subtract -a $a -b $b > $sub_bed

wc_out=$(wc -l $sub_bed)
line_count=$(echo $wc_out | awk '{print $1}')

if [ $line_count -gt 300 ]; then
    echo "The number of lines in the output bed file is $line_count, it will be truncated to 300 lines"
    head -n 300 $sub_bed > temp.bed
    mv temp.bed $sub_bed
fi

# Step 3: make report
output="${ar1[0]}_subtracted_${ar2[0]}_report.html"
create_report $sub_bed \
--fasta $reads \
--flanking 100 \
--tracks $sub_bed $a $b $bam_track \
--output $output

if [ ! -d "reports" ]; then
    mkdir reports
fi

if [ ! -d "reports/bed_subtraction" ]; then
    mkdir reports/bed_subtraction
fi

echo "moving files to reports directory"
mv $sub_bed reports/bed_subtraction
mv $output reports

