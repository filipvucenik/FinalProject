import sys


def load_true_annotations(path):
    ret = []
    try:
        with open(path, 'r') as file:
            for line in file:
                line = line.replace("\n", "")
                split = line.split(" ")
                if split[-1] == "":
                    split = split[:-1]
                ret.append([int(num) for num in split])
    except FileNotFoundError:
        print(f"File '{path}' not found.")
    except Exception as e:
        print("An error occurred:", e)

    return ret


def load_snp_annotations(path):
    ret = {}
    try:
        with open(path, 'r') as file:
            for line in file:
                split = line.split(" ")
                split = [int(num) for num in split]
                ret[int(split[0])] = set(split[1:])
    except FileNotFoundError:
        print(f"File '{path}' not found.")
    except Exception as e:
        print("An error occurred:", e)

    return ret


def load_read_lengths(path):
    ret = []
    try:
        with open(path, 'r') as file:
            for line in file:
                if line.startswith(">"):
                    continue
                ret.append(len(line))
    except FileNotFoundError:
        print(f"File '{path}' not found.")
    except Exception as e:
        print("An error occurred:", e)

    return ret


def main():
    if len(sys.argv) != 4:
        print("Usage: python annotation_comparator <reads> <true_annotations> <snps_annotations>")
        sys.exit(1)

    reads_file = sys.argv[1]
    true_annotations_file = sys.argv[2]
    snps_annotations_file = sys.argv[3]

    lens = load_read_lengths(reads_file)
    true_annotations = load_true_annotations(true_annotations_file)
    snps_annotations = load_snp_annotations(snps_annotations_file)

    TT = 0
    FF = 0
    TF = 0
    FT = 0

    print(len(lens))
    print(len(true_annotations))
    print(len(snps_annotations))

    for i in range(0, len(lens)):
        for j in range(0, lens[i]):
            if j in true_annotations[i] and i in snps_annotations and j in snps_annotations[i]:
                TT += 1
            if j in true_annotations[i] and i in snps_annotations and j not in snps_annotations[i]:
                TF += 1
            if j not in true_annotations[i] and i in snps_annotations and j in snps_annotations[i]:
                FT += 1
            if j not in true_annotations[i] and i in snps_annotations and j not in snps_annotations[i]:
                FF += 1

    print(f"{str(TT)} {str(TF)}")
    print(f"{str(FT)} {str(FF)}")


if __name__ == "__main__":
    main()
