import sys 

def check(reference, vcf):
    cout_same = 0
    cout_diff = 0
    for line in vcf:
        if line[0] == '#' or line[0] == ' ' or not line.startswith("chr19") :
            continue
        split = line.split()
        if len(split[3]) == 1:
            if split[3] == reference[int(split[1])-1]:
                cout_same += 1
            else:
                cout_diff += 1
                    
    print(f"Same: {cout_same} Different: {cout_diff}")

def main():
    if len(sys.argv) < 4:
        print("Usage: python vcfCheck.py <vcfPat> <vcfMat> <reference>")
        exit(1)
    vcfpat = sys.argv[1]
    vcfmat = sys.argv[2]
    referencefile = sys.argv[3]
    
    mat = True
    with open(referencefile, 'r') as ref:
        for line in ref:
            if line[0] == '>':
                continue
            reference = line.rstrip()
            if mat:
                mat = False
                with open(vcfmat, 'r') as vcf:
                    print("Checking maternal")
                    check(reference, vcf)
            else:
                with open(vcfpat, 'r') as vcf:
                    print("Checking paternal")
                    check(reference, vcf)

    
if __name__ == "__main__":
    main()