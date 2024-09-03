import pysam
import sys

def paf_to_sam(paf_file, sam_file):
    # Open the PAF file
    with open(paf_file, 'r') as paf:
        # Create a SAM file writer
        header = {
            'HD': {'VN': '1.0'},
            'SQ': []
        }
        with pysam.AlignmentFile(sam_file, 'w', header=header) as sam:
            for line in paf:
                fields = line.strip().split('\t')
                
                if(len(fields) < 13):
                    print("Weird line: ", line)
                    continue
                qname = fields[0]
                flag = 0  # You may need to set this based on the alignment
                rname = fields[5]
                pos = int(fields[7])
                mapq = 255  # Placeholder value
                cigar = fields[11]
                pnext = 0


                # Write the SAM record
                sam.write(a)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python P2S.py <input.paf> <output.sam>")
        sys.exit(1)
    
    paf_file = sys.argv[1]
    sam_file = sys.argv[2]
    paf_to_sam(paf_file, sam_file)