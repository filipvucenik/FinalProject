import sys




def main():
    if len(sys.argv) < 3:
        print("Usage: python HifiasmRunAnalysis.py <directory with hifiasm output> <sh1>")
        sys.exit(1)


    h0 = sys.argv[1]+"/hifiasm.asm.0.ovlp.paf"
    h1 = sys.argv[1]+"/hifiasm.asm.1.ovlp.paf"
    sh1 = sys.argv[2]
    correct0 = 0
    allovl0 = 0
    with open(h0) as f:
        for line in f:
            line = line.strip().split("\t")
            r1P = line[0].__contains__(sh1)
            r2P = line[5].__contains__(sh1)
            if r1P == r2P:
                correct0 += 1
            allovl0 += 1

    print("Same haplotype")
    print(correct0/allovl0)

    correct1 = 0
    allovl1 = 0
    with open(h1) as f:
        for line in f:
            line = line.strip().split("\t")
            r1P = line[0].__contains__(sh1)
            r2P = line[5].__contains__(sh1)
            if r1P != r2P:
                correct1 += 1
            allovl1 += 1

    print("Different haplotype")
    print(correct1/allovl1)

    print("Total accuracy")
    print((correct0+correct1)/(allovl0+allovl1))

if __name__ == "__main__":
    main()