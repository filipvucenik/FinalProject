import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python separateChr19Ref.py <reference>")
        exit(1)
        
    referencefile = sys.argv[1]
    mat = True
    mat_file = "maternal_chr19.fa"
    pat_file = "paternal_chr19.fa"
    
    with open(referencefile, 'r') as ref:
        for line in ref:
            if line[0] == '>':
                pref = line
                continue
            reference = line.rstrip()
            if mat:
                mat = False
                with open(mat_file, 'w') as mat_file:
                    mat_file.write(f"{pref}{reference}")
            else:
                with open(pat_file, 'w') as pat_file:
                    pat_file.write(f"{pref}{reference}")
                    

if __name__ == "__main__":
    main()