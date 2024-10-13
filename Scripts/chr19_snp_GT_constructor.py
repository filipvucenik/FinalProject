import sys

badread = True

def binary_search_closest(arr, target):
    left = 0
    right = len(arr) - 1
    while left <= right:
        mid = (left + right) // 2
        if arr[mid] == target:
            return mid, mid
        elif arr[mid] > target:
            right = mid - 1
        else:
            left = mid + 1
            
    return max(0, left), min(len(arr) - 1, right)
        
def main():
    if(len(sys.argv) < 4):
        print("Usage: python chr19_snp_GT_constructor.py <reads_file> <GT_PATERNAL> <GT_MATERNAL> <OPTIONAL: vcf flag>")
        sys.exit(1)
    reads_file = sys.argv[1]
    GT_PATERNAL = sys.argv[2]
    GT_MATERNAL = sys.argv[3]
    off = 0
    if len(sys.argv) == 5:
        vcf = sys.argv[4]
        if vcf == "vcf":
            off = 1
    
    snp_out = "chr19_p_snp_GT.anno"
    if badread:
        snp_out = "badread_50X_GT.anno"
    with open(GT_PATERNAL, "r") as f:
        parental_gt = list(map(int,f.readlines()[1:]))
    with open(GT_MATERNAL, "r") as f:
        maternal_gt = list(map(int,f.readlines()[1:]))
    with open(reads_file, "r") as f:
        with open(snp_out, "w") as out:
            lines = f.readlines()
            for i in range(0, len(lines), 2):
                if lines[i] == "\n":
                    continue
                
                if badread:
                    split = lines[i][1:-1].split("|")
                    pos = split[2].split(";")[0].split("-")
                    start =int(pos[0])
                    end = int(pos[1])
                    
                    reversed = split[1] == "-strand"
                    
                    if split[0].endswith("PATERNAL"):
                        gt = parental_gt
                    else:
                        gt = maternal_gt
                    
                else:
                    split = lines[i][1:].split("_")
                    pos = split[2].split("=")[1].split("-")
                    start = int(pos[0])
                    end = int(pos[1])
                    
                    reversed = split[1] == "reverse"
                                
                    if split[4].endswith("PATERNAL"):
                        gt = parental_gt
                    else:
                        gt = maternal_gt
                start_idx, _ = binary_search_closest(gt, start)
                _, end_idx = binary_search_closest(gt, end)
                if start_idx > end_idx:
                    continue
                gt_segment = gt[start_idx:end_idx]
               
                if len(gt_segment) == 0:
                    continue
                
                if reversed:
                    reduced_segment = [(end - x - 1 + off) for x in gt_segment]
                else:
                    reduced_segment = [(x - start - off) for x in gt_segment]
                out.write(f"{lines[i][1:-1]}")
                for x in reduced_segment:
                    if x < len(lines[i+1]) and x >= 0:
                        if lines[i+1][x] != "\n":
                            out.write(f" {x}:{lines[i+1][x]}")
                    else:
                        out.write(f" {x}:N")
                out.write(f"\n")
            

        
if __name__ == "__main__":
    main()