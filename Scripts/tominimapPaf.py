import sys 


def convertCigarToMinimap(cigar):
    minimap_cigar = ""
    star_pos = 0
    prev_match = False
    current_num = 0
    total_matches = 0
    total_alignment_length = 0
    for i in range(len(cigar)):
        if cigar[i].isdigit():
            continue
       
        if cigar[i] == '=' or cigar[i] == 'X':
            total_alignment_length += int(cigar[star_pos:i])
            total_matches += int(cigar[star_pos:i])
            if prev_match:
                current_num += int(cigar[star_pos:i])
            else:
                prev_match = True
                current_num = int(cigar[star_pos:i])
        elif cigar[i] == 'I':
            total_alignment_length += int(cigar[star_pos:i])
            if prev_match:
                minimap_cigar += f"{current_num}M"
                current_num = 0
            minimap_cigar += f"{int(cigar[star_pos:i])}I"
            prev_match = False
        elif cigar[i] == 'D':
            total_alignment_length += int(cigar[star_pos:i])
            if prev_match:
                minimap_cigar += f"{current_num}M"
                current_num = 0
            minimap_cigar += f"{int(cigar[star_pos:i])}D"
            prev_match = False
        star_pos = i + 1
            
    if prev_match:
        minimap_cigar += f"{current_num}M"
        
    return minimap_cigar, total_matches, total_alignment_length
            
            
        
        
    


def main():
    if len(sys.argv) < 3:
        print("Usage: python tominimapPaf.py <paf_file> <output_file>")
        exit(1)
        
        
    paf_file = sys.argv[1]
    out_file = sys.argv[2]
    
    with open(paf_file, 'r') as file:
        with open(out_file, "w") as out:
            for line in file:
                split = line.split("\t")
                read1_name = split[0]
                read1_len = split[1]
                read1_start = split[2]
                read1_end = split[3]
                strand = split[4]
                read2_name = split[5]
                read2_len = split[6]
                read2_start = split[7]
                read2_end = split[8]
                cigar  = split[11]
                minimap_cigar, total_matches, total_alignment_length = convertCigarToMinimap(cigar)
                out.write(f"{read1_name}\t{read1_len}\t{read1_start}\t{read1_end}\t{strand}\t{read2_name}\t{read2_len}\t{read2_start}\t{read2_end}\t{total_matches}\t{total_alignment_length}\t{0}\tNM:i:0\tms:i:0\tAS:i:0\tnn:i:0\ttp:A:S\tcm:i:0\ts1:i:0\tde:f:0\trl:i:0\tcg:Z:{minimap_cigar}\n")
                



if __name__ == '__main__':
    main()

