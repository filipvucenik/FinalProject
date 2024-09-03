import sys



def main():
    if len(sys.argv) < 5:
        print(sys.argv)
        print("Usage: python filterDifferentCigars.py <first_paf_file> <first_alignment_tool_name> <second_paf_file> <second_alignment_tool_name> ")
        exit(1)
    first_paf_file = sys.argv[1]
    fatn = sys.argv[2]
    second_paf_file = sys.argv[3]
    satn = sys.argv[4]
    
    pref =  "_comp_" + fatn + "_" + satn 
    
    filtered_first_paf = first_paf_file.replace(".paf", pref + "_filtered.paf")
    filtered_second_paf =  second_paf_file.replace(".paf", pref + "_filtered.paf")
    cigars_different = "different_cigars"+pref+".txt"
    
    with open(first_paf_file, 'r') as first_file:
        first_lines = first_file.readlines()
    with open(second_paf_file, 'r') as second_file:
        second_lines = second_file.readlines()
        
    with open(filtered_first_paf, 'w') as first_out:
        with open(filtered_second_paf, 'w') as second_out:
            with open(cigars_different, 'w') as diff_out:
                for i in range(len(first_lines)):
                    first_line = first_lines[i].split()
                    second_line = second_lines[i].split()
                    if first_line[len(first_line)-1] != second_line[len(second_line)-1]:
                        first_out.write("\t".join(first_line) + "\n")
                        second_out.write("\t".join(second_line) + "\n")
                        diff_out.write(f"{first_line[len(first_line)-1]}\n{second_line[len(second_line)-1]}\n\n")
    
    
if __name__ == "__main__":
    main()