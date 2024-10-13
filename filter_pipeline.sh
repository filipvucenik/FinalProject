#!/bin/bash

# Initialize variables with default values
reads=""
overlaps=""
number=10
start=0
file=""
name=""
annos=""

# Function to display usage
usage() {
    echo "Usage: $0 -r <reads> -o <overlaps> -a <annotations> --name <name> [-n <number>] [-s <start>] [-f <file>]"
    exit 1
}

# Parse command-line arguments
while getopts ":r:o:n:s:f:a:-:" opt; do
    case $opt in
        r) reads="$OPTARG" ;;
        o) overlaps="$OPTARG" ;;
        n) number="$OPTARG" ;;
        s) start="$OPTARG" ;;
        f) file="$OPTARG" ;;
        a) annos="$OPTARG" ;;
        -)
            case "${OPTARG}" in
                name)
                    val="${!OPTIND}"; OPTIND=$((OPTIND + 1))
                    name="$val"
                    ;;
                *)
                    echo "Unknown option --${OPTARG}"
                    usage
                    ;;
            esac
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

# Check required arguments
if [ -z "$reads" ] || [ -z "$overlaps" ] || [ -z "$annos" ] || [ -z $name ]; then
    usage
fi

# Run filterReadsOverlaps.py
echo "Current Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"
 
python Scripts/filterReadsOverlaps.py -r "$reads" -o "$overlaps" -n "$number" -s "$start" -f "$file" --output "$name"


output_reads="$name.fasta"
output_overlaps="$name.paf"
output_read_names="$name.txt"

# Run visualization script
./visualization_pipeline.sh $output_overlaps $output_reads $name

# Run BED file creation script

python Scripts/createBed.py $output_read_names $annos $name

if [ ! -d "filtered" ]; then
    mkdir filtered
fi

if [ ! -d "filtered/bed_files" ]; then
    mkdir filtered/bed_files
fi

if [ ! -d "filtered/reads" ]; then
    mkdir filtered/reads
fi

if [ ! -d "filtered/overlaps" ]; then
    mkdir filtered/overlaps
fi

if [ ! -d  "filtered/bam_files" ]; then
    mkdir filtered/bam_files
fi

echo "moving files to filtered directory"
mv $output_reads filtered/reads
mv "${output_reads}.fai" filtered/reads
mv $output_overlaps filtered/overlaps
mv $output_read_names filtered/reads
mv $name.bed filtered/bed_files
bam_file="${name}_mp2_sorted.bam"
bai_file="${name}_mp2_sorted.bam.bai"
mv visualization/BAM_files/$bam_file filtered/bam_files
mv visualization/BAM_files/$bai_file filtered/bam_files
