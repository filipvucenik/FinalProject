import sys


def main():
    if len(sys.argv) < 3:
        print("Usage: python snpBenchmark.py <snp_annotations> <ground_truth>")
        exit(1)
        
    snp_annotations = sys.argv[1]
    ground_truth = sys.argv[2]
    
    reads = set()
    anno = {}
    with open(snp_annotations, "r") as f:
        for line in f:
            split = line.split(" ")
            anno[int(split[0])] = set([int(x) for x in split[1:]])
            reads.add(int(split[0]))
    print("Loaded annotations")
    gt = {}
    
    with open(ground_truth, "r") as f:
        for line in f:
            split = line.split(" ")
            gt[int(split[0])] = set([int(x) for x in split[1:]])
            reads.add(int(split[0]))
            
    print("Loaded ground truth")
    tp = 0
    fp = 0
    fn = 0
    
    for key in reads:
        if key in gt and key in anno:
            s1 = anno[key]
            s2 = gt[key]
            i = s1.intersection(s2)
            tp += len(i)
            fp += len(s1) - len(i)
            fn += len(s2) - len(i)
        elif key in gt:
            fn += len(gt[key])
        elif key in anno:
            fp += len(anno[key])
    print("Benchmark for all positions")
    print("True Positives: ", tp)
    print("False Positives: ", fp)
    print("False Negatives: ", fn)
    print("True Negatives: X")
    
    
    
            
if __name__ == "__main__":
    main()