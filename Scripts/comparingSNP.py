import sys 


def main ():
    if len(sys.argv) < 4:
        print("Usage: python comparingSNP.py <file1> <file2> <output_file>")
        sys.exit(1)
        
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    d1 = {}
    d2 = {}
    
    with open(file1, "r") as f1:
        for line in f1:
            split = line.split()
            read_ind = int(split[0])
            d1[read_ind] = set([int(x) for x in split[1:]])
            
    with open(file2, "r") as f2:
        for line in f2:
            split = line.split()
            read_ind = int(split[0])
            d2[read_ind]  = set([int(x) for x in split[1:]])

    out_file = sys.argv[3]
    with open(out_file, "w") as out:
        out.write(f"#good read wrong snp calls 1:{file1} 2:{file2}\n")
        for key in d1:
            if key in d2:
                if d1[key] != d2[key]:
                    l1 = list(d1[key])
                    l2 = list(d2[key])
                    l1.sort()
                    l2.sort()
                    out.write(f"{key}\n{' '.join(map(str, l1))}\n{' '.join(map(str, l2))}\n")
        out.write(f"#reads present in {file1} but not in {file2}\n")
        for key in d1:
            if key not in d2:
                l = list(d1[key])
                l.sort()
                out.write(f"{key}\n{' '.join(map(str, l))}\n")
        out.write(f"#reads present in {file2} but not in {file1}\n")
        for key in d2:
            if key not in d1:
                l = list(d2[key])
                l.sort()
                out.write(f"{key}\n{' '.join(map(str, l))}\n")
                
if __name__ == "__main__":
    main()
    