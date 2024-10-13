import argparse

param = 0

def checkSameHaplotype(ovlp):
    split = ovlp.split("\t")
    lhs_origin = split[0].split("|")[0].split(";")[1]
    rhs_origin = split[5].split("|")[0].split(";")[1]
    
    return lhs_origin == rhs_origin
    
    
def checkSameRegion(ovlp, same_hpl=False):
    extend=0
    if not same_hpl:
        extend = param
    
    split = ovlp.split("\t")
    lhs_ovlp_begin =  int(split[0].split("|")[2].split(";")[0].split("-")[0]) - extend
    lhs_ovlp_end =   int(split[0].split("|")[2].split(";")[0].split("-")[1]) + extend

    rhs_ovlp_begin =  int(split[5].split("|")[2].split(";")[0].split("-")[0]) - extend
    rhs_ovlp_end =  int(split[5].split("|")[2].split(";")[0].split("-")[1]) + extend
    
    if rhs_ovlp_begin >= lhs_ovlp_begin and rhs_ovlp_begin <= lhs_ovlp_end or rhs_ovlp_end >= lhs_ovlp_begin and rhs_ovlp_end <= lhs_ovlp_end:
        return True
    return False
    


def main():
    parser = argparse.ArgumentParser(description="Checks how good badread overlaps are classified")
    parser.add_argument('-s', '--same_haplotype', type=str, help="path to paf file with same haplotype overlaps", required=True)
    parser.add_argument('-d', '--different_haplotype', type=str, help="path to paf file with different haplotype overlaps", required=True)
    parser.add_argument('-w', '--wrong', type=str, help="path to paf file with wrong overlaps", required=False)
    parser.add_argument('-n', '--no_snps', type=str, help="path to paf file with no snps positions", required=False)
    
    
    args = parser.parse_args()
    same_file = args.same_haplotype
    different_file = args.different_haplotype
    if args.wrong:
        wrong_file = args.wrong
    if args.no_snps:
        no_snps_file = args.no_snps

    counter_same_true = 0
    counter_same_false = 0    
    counter_good_ovlp_st = 0
    counter_wrong_ovlp_st = 0
    counter_good_ovlp_sf = 0
    counter_wrong_ovlp_sf = 0
    
    wrong_same_hpl = open("wrong_same", 'w')
    wrong_different = open("wrong_different", "w")
    good_region_not_called = open("godd_region_not_called", "w")
    
    counter1 = 10
    counter1 = 10    
    with open(same_file, "r") as f:
        for line in f:
            same = checkSameHaplotype(line)
            if same:
                counter_same_true+=1
            else:
                counter_same_false+=1
            good_ovlp = checkSameRegion(line, same)
            
            if same and good_ovlp:
                counter_good_ovlp_st+=1
            elif same and not good_ovlp:
                counter_wrong_ovlp_st+=1
            elif not same and good_ovlp:
                counter_good_ovlp_sf+=1
            elif not same and not good_ovlp:
                counter_wrong_ovlp_sf+=1
        
    print("Same haplotype:")
    print(f"same haplotype and overlaping region {counter_good_ovlp_st}")
    print(f"different haplotype and overlaping region {counter_good_ovlp_sf}")
    print(f"same haplotype and non overlaping region{counter_wrong_ovlp_st}")
    print(f"different haplotype and non overlaping region {counter_wrong_ovlp_sf}")

    counter_same_true = 0
    counter_same_false = 0    
    counter_good_ovlp_st = 0
    counter_wrong_ovlp_st = 0
    counter_good_ovlp_sf = 0
    counter_wrong_ovlp_sf = 0
    
    with open(different_file, 'r') as f:
        for line in f:
            same = checkSameHaplotype(line)
            if same:
                counter_same_true+=1
            else:
                counter_same_false+=1
            good_ovlp = checkSameRegion(line, same)
            if same and good_ovlp:
                counter_good_ovlp_st+=1
            elif same and not good_ovlp:
                counter_wrong_ovlp_st+=1
            elif not same and good_ovlp:
                counter_good_ovlp_sf+=1
            elif not same and not good_ovlp:
                counter_wrong_ovlp_sf+=1
    
    print("Different haplotype:")
    print(f"same haplotype and overlaping region {counter_good_ovlp_st}")
    print(f"different haplotype and overlaping region {counter_good_ovlp_sf}")
    print(f"same haplotype and non overlaping region{counter_wrong_ovlp_st}")
    print(f"different haplotype and non overlaping region {counter_wrong_ovlp_sf}")
    
       
    if args.wrong: 
        counter_same_true = 0
        counter_same_false = 0    
        counter_good_ovlp_st = 0
        counter_wrong_ovlp_st = 0
        counter_good_ovlp_sf = 0
        counter_wrong_ovlp_sf = 0
        with open(wrong_file, 'r') as f:
            for line in f:
                same = checkSameHaplotype(line)
                if same:
                    counter_same_true+=1
                else:
                    counter_same_false+=1
                good_ovlp = checkSameRegion(line, same)
                if same and good_ovlp:
                    counter_good_ovlp_st+=1
                elif same and not good_ovlp:
                    counter_wrong_ovlp_st+=1
                elif not same and good_ovlp:
                    counter_good_ovlp_sf+=1
                elif not same and not good_ovlp:
                    counter_wrong_ovlp_sf+=1
        
        print("wrong overlaps:")
        print(f"same haplotype and overlaping region {counter_good_ovlp_st}")
        print(f"different haplotype and overlaping region {counter_good_ovlp_sf}")
        print(f"same haplotype and non overlaping region{counter_wrong_ovlp_st}")
        print(f"different haplotype and non overlaping region {counter_wrong_ovlp_sf}")
    
    if args.no_snps:
        counter_same_true = 0
        counter_same_false = 0    
        counter_good_ovlp_st = 0
        counter_wrong_ovlp_st = 0
        counter_good_ovlp_sf = 0
        counter_wrong_ovlp_sf = 0
        
        with open(no_snps_file, 'r') as f:
            for line in f:
                same = checkSameHaplotype(line)
                if same:
                    counter_same_true+=1
                else:
                    counter_same_false+=1
                good_ovlp = checkSameRegion(line, same)
                if same and good_ovlp:
                    counter_good_ovlp_st+=1
                elif same and not good_ovlp:
                    counter_wrong_ovlp_st+=1
                elif not same and good_ovlp:
                    counter_good_ovlp_sf+=1
                elif not same and not good_ovlp:
                    counter_wrong_ovlp_sf+=1
        
        print("no snps overlaps:")
        print(f"same haplotype and overlaping region {counter_good_ovlp_st}")
        print(f"different haplotype and overlaping region {counter_good_ovlp_sf}")
        print(f"same haplotype and non overlaping region{counter_wrong_ovlp_st}")
        print(f"different haplotype and non overlaping region {counter_wrong_ovlp_sf}")
            
if __name__ == "__main__":
    main()