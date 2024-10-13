import sys



def main():
    if len(sys.argv) != 2:
        print("Usage: python herroToSNPAnno.py <herro_file>")
        sys.exit(1)
    herro_file = sys.argv[1]
    output_file = "chr19_e_snp_GT_no_indel.anno"
    with open(herro_file, 'r') as f:
        with open(output_file, 'w') as o:
            i = 0
            for line in f:
                split = line.split("\t")[1:]
                o.write(f"{i} ")
                o.write(" ".join(split))     
                i += 1

if __name__ == "__main__":
    main()