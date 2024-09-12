import sys


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
    if len(sys.argv) < 3:
        print("Usage: python chr19_maternal_GT_constructor.py <reads> <maternal_GT_ref>")
        exit(1)
    reads_file = sys.argv[1]
    GT_MATERNAL = sys.argv[2]
    
    with open(GT_MATERNAL, "r") as f:
        gt = list(map(int, f.readlines()[1:]))
        
    snp_out = "chr19_m_snp_GT.anno"
    
    with open(reads_file, "r") as f:
        with open(snp_out, "w") as out:
            for line in f:
                if line[0] != ">":
                    continue
                split = line[1:].split("_")
                rn = int(split[0].split("=")[1]) - 1
                reverse = split[1] == "reverse"
                pos = split[2].split("=")[1].split("-")
                start = int(pos[0])
                end = int(pos[1])    
                
                start_idx, _ = binary_search_closest(gt, start)
                _, end_idx = binary_search_closest(gt, end)
                gt_segment = gt[start_idx:end_idx]
                
                if len(gt_segment) == 0:
                    continue
                
                if reverse:
                    reduced_segment = [end - x - 1 for x in gt_segment]
                    reduced_segment.reverse()
                else:
                    reduced_segment = [x - start for x in gt_segment]
                out.write(f"{rn} {' '.join(map(str, reduced_segment))}\n")
    
    
if __name__ == "__main__":
    main()
