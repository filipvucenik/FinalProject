import sys

def main ():
    if len(sys.argv) != 2:
        print("Usage: python changeCmaForDsh.py <reads>")
        sys.exit(1)
        
    reads = sys.argv[1]
    suffix = "." + reads.split(".")[-1]
    out = reads.replace(suffix, "_dsh"+suffix)
    if out == reads:
        print("Error: reads file must have a suffix")
        sys.exit(1)
    with open(reads, "r") as f:
        with open(out, "w") as f2:
            for line in f:
                f2.write(line.replace(",", "_"))

                    
if __name__ == "__main__":
    main()