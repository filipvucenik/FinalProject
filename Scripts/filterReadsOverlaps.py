import argparse


def main():
    parser = argparse.ArgumentParser(description='Filter reads and overlaps for easier visualization')
    parser.add_argument('-r','--reads', type=str, help='reads file', required=True )
    parser.add_argument('-o','--overlaps', type=str, help='overlaps file', required=True )
    parser.add_argument('-n', '--number', type=int, help='number of reads to filter', required=False, default=10)
    parser.add_argument('-s', '--start', type=int, help='start filtering from this read', required=False, default=0)
    parser.add_argument('-f', '--file', type=str, help='read names file', required=False)
    
     
    args = parser.parse_args()
    reads = args.reads
    overlaps = args.overlaps
    number = args.number
    start = args.start
    
    reads_suffix = '.'+reads.split(".")[-1]
    out_reads = reads.replace(reads_suffix, "_filtered"+reads_suffix)
    if out_reads == reads:
        out_reads = reads+"_filtered.fasta"
    overlaps_suffix = '.'+overlaps.split(".")[-1]
    out_overlaps = overlaps.replace(overlaps_suffix, "_filtered"+overlaps_suffix)
    if out_overlaps == overlaps:
        out_overlaps = overlaps+"_filtered.paf"
    out_reads_notation = reads.replace(reads_suffix, "_filtered.txt")
    if out_reads_notation == reads:
        out_reads_notation = reads+"_filtered.txt"
    
    reads_set = set()
    if args.file:
        with open(args.file, 'r') as f:
            reads_set = set(x[:-1] for x in f.readlines() if x[0] != '#')
        with open(reads, 'r') as f2:
            with open(out_reads, 'w') as f3:
                reads_lines = f2.readlines()
                for i in range(len(reads_lines)):
                    if reads_lines[i][1:-1] in reads_set:           
                        f3.write(reads_lines[i])
                        f3.write(reads_lines[i+1])
    else:
        with open(reads, 'r') as f:
            reads_lines = f.readlines()
            filtered_reads = reads_lines[start:start+2*number]
            reads_set = set(x[1:-1] for x in filtered_reads if x[0] == '>')
            with open(out_reads, 'w') as f2:
                f2.writelines(filtered_reads)
            
    added_reads = set(reads_set)
    with open(overlaps, 'r') as f:
        with open(out_overlaps, 'w') as f2:
                with open(out_reads, 'a') as ro:
                    for line in f:
                        split_line = line.split("\t")
                        if split_line[0] in reads_set:
                            f2.write(line)
                            for i in range(len(reads_lines)):
                                if reads_lines[i][1:-1] == split_line[5] and split_line[5] not in added_reads:
                                    ro.write(reads_lines[i])
                                    ro.write(reads_lines[i+1])
                                    added_reads.add(split_line[5])
                                    break
                                
                        if split_line[5] in reads_set:
                            f2.write(line)
                            for i in range(len(reads_lines)):
                                if reads_lines[i][1:-1] == split_line[0] and split_line[0] not in added_reads:
                                    ro.write(reads_lines[i])
                                    ro.write(reads_lines[i+1])
                                    added_reads.add(split_line[0])
                                    break
                                

    with open(out_reads_notation, 'w') as f:
        f.write("#Filtered reads:\n")
        for read in reads_set:
            f.write(read+'\n')
                            
                    
        
          
if __name__ == "__main__":
    main()
    
    
    