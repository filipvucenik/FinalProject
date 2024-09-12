import argparse


def main():
    parser = argparse.ArgumentParser(description='Filter reads and overlaps for easier visualization')
    parser.add_argument('-r','--reads', type=str, help='reads file', required=True )
    parser.add_argument('-n', '--number', type=int, help='number of reads to filter', required=False, default=10)
    parser.add_argument('-s', '--start', type=int, help='start filtering from this read', required=False, default=0)

    args = parser.parse_args()
    reads = args.reads
    number = args.number
    start = args.start
    
    reads_suffix = '.'+reads.split(".")[-1]
    out_reads = reads.replace(reads_suffix, "_filtered.txt")
    
    if out_reads == reads:
        out_reads = reads + "_filtered.txt"
    
    with open(reads, 'r') as f:
        reads_lines = f.readlines()
        filtered_reads = reads_lines[start:start+2*number]
        reads_set = set(x[1:] for x in filtered_reads if x[0] == '>')
        
        with open(out_reads, 'w') as f2:
            f2.writelines(reads_set)
            



if __name__ == "__main__":
    main()