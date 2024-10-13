import sys


def main():
    if(len(sys.argv)) < 4:
        print("Usage: python chr19_p_GT_constructor.py <vcf_file> <pref> <parantage> <optional indel flag I>")
        exit(1)
    input_file = sys.argv[1]
    pref = sys.argv[2]
    parantage = sys.argv[3]
    out_file = parantage + "_" +  pref + "_p_GT.txt"
    indel = False
    if len(sys.argv) == 5:
        if sys.argv[4] == "I":
            out_file = parantage + "_" + pref + "_p_GT_indel.txt"  
            indel = True         
    with open (input_file, 'r') as file:
        with open(out_file, 'w') as out:
            out.write(f"#Chromosome {pref} {parantage}\n")
            for line in file:
                if line[0] == '#' or not line.startswith(pref):
                    continue
                split = line.split()
                if(len(split[3]) == 1 and len(split[4]) == 1 and split[3] != split[4]):
                    out.write(f"{split[1]}\n")
                elif len(split[3]) > len(split[4]) and indel: # Deltion
                    s_pos = int(split[1])
                    e_pos = s_pos + len(split[3])
                    for i in range(s_pos, e_pos):
                        out.write(f"{i}\n")
                elif len(split[3]) < len(split[4]) and indel: # Insertion
                    s_pos = int(split[1])
                    e_pos = s_pos + len(split[4])
                    for i in range(s_pos, e_pos):
                        out.write(f"{i}\n")
                    
                
        

if __name__ == "__main__":
    main()
    
