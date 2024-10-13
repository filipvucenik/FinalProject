import sys

def check(regions, begin, end):
    for r in regions:
        if r[0] >= begin and r[1] <= end:
            return True
    return False

def check_read(r_fw,r_origin, r_begin, r_end, r_ovlp_begin, r_ovlp_end, maternal_regions, paternal_regions):
    if r_fw:
        begin = r_begin + r_ovlp_begin
        end = r_begin + r_ovlp_end
    else:
        begin = r_end - r_ovlp_end
        end = r_end - r_ovlp_begin
    
    if (r_origin.__contains__("MATERNAL")):
            if check(maternal_regions, begin, end):
                return True
    else:
        if check(paternal_regions, begin, end):
            return True

def badread(line, maternal_regions, paternal_regions):
    split = line.split("\t")
    r1 = split[0]
    r1_split = r1.split("|")
    r1_origin = r1_split[0].split(";")[1]
    r1_fw = r1_split[1] == "+strand"
    r1_begin = int(r1_split[2].split(";")[0].split("-")[0]) 
    r1_end = int(r1_split[2].split(";")[0].split("-")[1])
    r1_ovlp_begin = int(split[2])
    r1_ovlp_end = int(split[3])
    
    r1_res = check_read(r1_fw, r1_origin, r1_begin, r1_end, r1_ovlp_begin, r1_ovlp_end, maternal_regions, paternal_regions)

    if r1_res:
        return True
    
    r2 = split[5]
    r2_split = r2.split("|")
    r2_origin = r2_split[0].split(";")[1]
    r2_fw = r2_split[1] == "+strand"
    r2_begin = int(r2_split[2].split(";")[0].split("-")[0]) 
    r2_end = int(r2_split[2].split(";")[0].split("-")[1])
    r2_ovlp_begin = int(split[7])
    r2_ovlp_end = int(split[8])
    
    r2_res = check_read(r2_fw, r2_origin, r2_begin, r2_end, r2_ovlp_begin, r2_ovlp_end, maternal_regions, paternal_regions)

    if r2_res:
        return True

    return False  
    
def seqrequester(line, maternal_regions, paternal_regions):
    split = line.split("\t")
    
    r1 = split[0]
    r1_split = r1.split("_")
    r1_origin = r1_split[-1]
    r1_fw = r1_split[1] == "forward"
    r1_begin = int(r1_split[2].split("=")[1].split("-")[0])
    r1_end = int(r1_split[2].split("=")[1].split("-")[1])
    r1_ovlp_begin = int(split[2])
    r1_ovlp_end = int(split[3])
    
    r1_res = check_read(r1_fw, r1_origin, r1_begin, r1_end, r1_ovlp_begin, r1_ovlp_end, maternal_regions, paternal_regions)

    if r1_res:
        return True
    
    r2 = split[0]
    r2_split = r2.split("_")
    r2_origin = r2_split[-1]
    r2_fw = r2_split[1] == "forward"
    r2_begin = int(r2_split[2].split("=")[1].split("-")[0])
    r2_end = int(r2_split[2].split("=")[1].split("-")[1])
    r2_ovlp_begin = int(split[2])
    r2_ovlp_end = int(split[3])
    
    r2_res = check_read(r2_fw, r2_origin, r2_begin, r2_end, r2_ovlp_begin, r2_ovlp_end, maternal_regions, paternal_regions)
    
    if r2_res:
        return True

    return False  

def main():
    if len(sys.argv) < 3:
        print("Usage: homozygousPercent.py <bed_file> <paf_file>")
        exit(1)
    
    bdr = True
    
    bed_file = sys.argv[1]
    paf_file = sys.argv[2]
    
    maternal_regions = []
    paternal_regions = []
    
    with open(bed_file, 'r') as f:
        for line in f:
            split = line[:-1].split("\t")
            if(split[0].__contains__("MATERNAL")):
                maternal_regions.append((int(split[1]), int(split[2])))
            else:
                paternal_regions.append((int(split[1]), int(split[2])))
    
    homozygous = 0
    not_homozygous = 0
    with open(paf_file, 'r') as f:
        for line in f:
            if bdr:
                hmz = badread(line[:-1], maternal_regions, paternal_regions)
            else:
                hmz = seqrequester(line[:-1], maternal_regions, paternal_regions)    
            
            if hmz:
                homozygous+=1
            else:
                not_homozygous+=1
    
    print(f"homozygous={homozygous} non_homozygous={not_homozygous} percentage={homozygous/(homozygous + not_homozygous)}")
    
if __name__ == "__main__":
    main()