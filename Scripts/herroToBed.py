import sys

def writeToBed (output, annos, read_name):
    for anno in annos:
        output.write(f"{read_name}\t{int(anno)}\t{int(anno)+1}\t1\n")

        
def main():
    if len(sys.argv) < 4:
        print("Usage: python createBed.py <read_names> <annotations> <output_file_name>")
        sys.exit(1)
    reads = sys.argv[1]
    annos = sys.argv[2]
    output_file = sys.argv[3] + ".bed"
    outputFile = open(output_file, "w")
    
    annosDict = {}
    
    with open(annos, "r") as annosFile:
        for line in annosFile:
            if line == "\n":
                continue
            split = line.split("\t")
            ind = split[0]
            annosDict[ind] = list(map(str, split[1:]))  
    
    with open(reads, "r") as readsFile:
        for line in readsFile:
            if line[0] == "#":
                continue
            ind = line[:-1]
            if ind not in annosDict:
                print(f"Warning: Index {ind} not found in annotations.")
                continue
            writeToBed(outputFile, annosDict[ind], line[:-1])

            
           
    

if __name__ == "__main__":
    main()