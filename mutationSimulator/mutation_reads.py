import sys


def read_mutations(file_path):
    mutations = []
    with open(file_path, 'r') as file:
        for line in file:
            mutations.append(int(line[:-1]))
    mutations.sort()
    return mutations


def read_reads(file_path):
    reads = []
    with open(file_path, 'r') as file:
        for line in file:
            reads.append(line[:-1])
    return reads


def equal_or_higher(arr, x):
    l = 0
    r = len(arr) - 1
    closest_ind = -1
    while l <= r:
        mid = l + (r - l) // 2
        if arr[mid] >= x:
            closest_ind = mid
            r = mid - 1
        else:
            l = mid + 1

    return closest_ind


def equal_or_lower(arr, x):
    l = 0
    r = len(arr) - 1
    closest_ind = -1
    while l <= r:
        mid = l + (r - l) // 2
        if arr[mid] <= x:
            closest_ind = mid
            l = mid + 1
        else:
            r = mid - 1

    return closest_ind


def main():
    if len(sys.argv) < 3:
        print("Usage: python mutation_reads.py <mutations_file> <reads_file>")
        sys.exit(1)

    mutations_file = sys.argv[1]
    reads_file = sys.argv[2]

    mutations = read_mutations(mutations_file)
    reads = read_reads(reads_file)

    reads_split = reads_file.split(".")
    output_file = ".".join(reads_split[:-1]) + ".out"

    output = []

    for i in range(0, len(reads), 2):
        read = reads[i]
        parts = read.split(",")
        interval = parts[2].split("=")[1]

        lower_bound = int(interval.split("-")[0])
        upper_bound = int(interval.split("-")[1])

        lower_bound_ind = equal_or_higher(mutations, lower_bound)
        upper_bound_ind = equal_or_lower(mutations, upper_bound)
        output.append(read)
        positions = []
        for j in range(lower_bound_ind, upper_bound_ind + 1):
            positions.append(str(mutations[j] - lower_bound))
        output.append("Positions=" + ",".join(positions))

    with open(output_file, 'w') as file:
        for line in output:
            file.write(line + "\n")


if __name__ == "__main__":
    main()
