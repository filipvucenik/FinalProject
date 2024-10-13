import argparse

def main():
    parser = argparse.ArgumentParser(description='takes n read names from BED files (default=10)')
    parser.add_argument('-b','--bed', type=str, help='bed file', required=True )
    parser.add_argument('-n','--number', type=int, help='read names', required=False, default=10)

    args = parser.parse_args()
    file = args.bed
    number = args.number
    
    reads = set()
    
    with open (file, 'r') as f:
        for line in f:
            if number == 0:
                break
            split = line.split("\t")
            if split[0] not in reads:
                reads.add(split[0])
                number-=1
                
    out_file = file.replace(".bed", ".txt")
    if out_file == file:
        out_file+=".txt"
        
    with open(out_file, "w") as f:
        f.writelines([line + "\n" for line in reads])

if __name__ == "__main__":
    main()
