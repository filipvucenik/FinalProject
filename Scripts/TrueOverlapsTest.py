import sys


def count_matches(s):
    count = 0
    pos = 0
    start_pos = 0

    if not s:
        return count

    for i in range(len(s)):
        if s[i].isdigit():
            if pos == 0:
                start_pos = i
            pos += 1
        else:
            if pos == 0:
                continue
            total_num = int(s[start_pos:start_pos + pos])
            pos = 0
            if s[i] == '=':
                count += total_num
            elif s[i] == 'X':
                count += total_num
    return count


def read_reference_genome(file_path):
    genome = ""

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            genome += line[:len(line) - 1]

    return genome


def main():
    if len(sys.argv) < 4:
        print("Usage: python TrueOverlapsTest.py <reference_genome> <mutated_reference_genome> <paf_file_PA>")
        sys.exit(1)
    reference_genome = read_reference_genome(sys.argv[1])
    mutated_reference_genome = read_reference_genome(sys.argv[2])
    paf_file = sys.argv[3]
    errors = 0
    all = 0
    out_file = paf_file + ".true"
    out_csv = paf_file + ".csv"
    with open(paf_file, 'r') as file:
        with open(out_file, "w") as out:
            for line in file:
                split = line.split("\t")
                read1 = split[0]
                lhs_pos = int(read1.split(",")[2].split("=")[1].split("-")[0])
                lhs_begin = int(split[2])
                lhs_end = int(split[3])
                read2 = split[5]
                rhs_pos = int(read2.split(",")[2].split("=")[1].split("-")[0])
                rhs_begin = int(split[7])
                rhs_end = int(split[8])
                lhs_rf = reference_genome[lhs_pos + lhs_begin:lhs_pos + lhs_end]
                lhs_rf_mut = mutated_reference_genome[lhs_pos + lhs_begin:lhs_pos + lhs_end]
                rhs_rf = reference_genome[rhs_pos + rhs_begin:rhs_pos + rhs_end]
                rhs_rf_mut = mutated_reference_genome[rhs_pos + rhs_begin:rhs_pos + rhs_end]

                match = count_matches(split[len(split)-1][:-1])
                len1 = split[1]
                len2 = split[6]

                if lhs_rf == rhs_rf or lhs_rf == rhs_rf_mut or lhs_rf_mut == rhs_rf or lhs_rf_mut == rhs_rf_mut:
                    out.write(read1 + " " + read2 + " "+str(match)+" "+len1+" "+len2+" 1\n")
                else:
                    out.write(read1 + " " + read2 + " "+str(match)+" "+len1+" "+len2+" 0\n")
                    errors += 1
                all += 1
    with open(out_csv, 'w') as csv:
        csv.write("matches,len1,len2,match\n")
        with open(out_file, 'r') as file:
            for line in file:
                split = line.split(" ")
                csv.write(split[2] + "," + split[3] + "," + split[4] + "," + split[5] + "\n")

    print("Done")
    print("Error rate: ", errors/all)


if __name__ == '__main__':
    main()
