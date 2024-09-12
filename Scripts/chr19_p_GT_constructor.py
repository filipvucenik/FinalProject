import sys


def main():
    if(len(sys.argv)) < 4:
        print("Usage: python chr19_p_GT_constructor.py <vcf_file> <pref> <parantage>")
        exit(1)
    input_file = sys.argv[1]
    pref = sys.argv[2]
    parantage = sys.argv[3]
    out_file = parantage + "_" +  pref + "_p_GT.txt"
    with open (input_file, 'r') as file:
        with open(out_file, 'w') as out:
            out.write(f"#Chromosome {pref} {parantage}\n")
            for line in file:
                if line[0] == '#' or not line.startswith(pref):
                    continue
                split = line.split()
                if(len(split[3]) == 1 and len(split[4]) == 1 and split[3] != split[4]):
                    out.write(f"{split[1]}\n")
                
        

if __name__ == "__main__":
    main()
    
