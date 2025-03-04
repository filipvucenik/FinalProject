import sys

def writeToBed (output, annos, read_name):
    error_file = "error"
    es = open(error_file, "w")
    for anno in annos:
        split = anno.split(":")
        if len(split) == 1:
            output.write(f"{read_name}\t{int(anno)}\t{int(anno)+1}\t1\n")
        else:
            if split[0].startswith("mq"):
                continue
            pos = int(split[0])
            marker = split[1]
            if marker != 'N':
                output.write(f"{read_name}\t{pos}\t{pos+1}\t{marker}\n")
            else:
                es.write(f"{read_name}\t{pos}\t{pos+1}\t{marker}\n")
                
        
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
            split = line.split(" ")
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

            
           
    
    
    name = sys.argv[3]
    

if __name__ == "__main__":
    main()