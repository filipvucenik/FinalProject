import sys 

def main ():
    if len(sys.argv) < 3:
        print("Usage: readsNotInAnno <anno_file> <read_index>")
        exit(1)
        
    anno_file = sys.argv[1]
    read_index_file = sys.argv[2]
    
    outFile = "not_in_anno.fasta.fai"
    
    with open(anno_file, 'r') as f:
        anno_lines = [x.split(" ")[0] for x in f.readlines()]
    
    with open(read_index_file, 'r') as f:
        reads_lines = [x.split("\t")[0] for x in f.readlines()]
    
    i = 0
    j = 0
    notInAnno = []
    while(i < len(reads_lines)):
        if anno_lines[j] != reads_lines[i]:
            notInAnno.append(reads_lines[i])
            i+=1
            continue
        i+=1
        j+=1
        
    with open(outFile, 'w') as f:
        f.write("\n".join(notInAnno))

if __name__ == "__main__":
    main()

