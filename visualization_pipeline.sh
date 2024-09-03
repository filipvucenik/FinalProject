#!/bin/bash

# Step 1: get input file and name of output name
paf_file=$1
reads=$2
output_name=$3

# Step 2: run the python script
minimap_output="${output_name}_mp2.paf"
python Scripts/tominimapPaf.py $paf_file $minimap_output
echo "converted to minimap2 paf ${minimap_output}"

# Step 3: run rustybam to convert paf to sam
conda activate snp
echo "Current Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"
sam_output="${output_name}_mp2.sam"
rustybam paf-to-sam -f $reads $minimap_output > $sam_output
echo "converted to sam ${sam_output}"

# Step 4: run samtools to convert sam to bam
bam_output="${output_name}_mp2.bam"
samtools view -bho $bam_output $sam_output
echo "converted to bam ${bam_output}"

# Step 5: run samtools to sort the bam file
sorted_bam_output="${output_name}_mp2_sorted.bam"
samtools sort -o $sorted_bam_output $bam_output
echo "sorted bam file ${sorted_bam_output}"

# Step 6: run samtools to index the bam file
samtools index $sorted_bam_output
echo "indexed bam file ${sorted_bam_output}"

# Step 7: index the reads
seqkit faidx $reads
echo "indexed reads"

# Step 8: initialize the output directory
if [ ! -d "visualization" ]; then
    mkdir visualization
fi
if [ ! -d "visualization/minimap2Paf" ]; then
    mkdir visualization/minimap2Paf
fi
if [ ! -d "visualization/SAM_files" ]; then
    mkdir visualization/SAM_files
fi
if [ ! -d "visualization/BAM_files" ]; then
    mkdir visualization/BAM_files
fi
if [ ! -d "visualization/indexed_reads" ]; then
    mkdir visualization/indexed_reads
fi
echo "created output directories"

# Step 9: move the files to the output directory
mv $minimap_output visualization/minimap2Paf
mv $sam_output visualization/SAM_files
mv $bam_output visualization/BAM_files
mv $sorted_bam_output visualization/BAM_files
mv ${sorted_bam_output}.bai visualization/BAM_files
mv ${reads}.fai visualization/indexed_reads
echo "done"