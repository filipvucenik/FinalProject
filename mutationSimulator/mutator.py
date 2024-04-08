import sys
import argparse
import random


def read_reference_genome(file_path):
    genome = ""
    name = ""

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                parts = line.split(" ")
                parts[0] = parts[0] + "_mutated"
                name = " ".join(parts)
            else:
                genome += line

    return name, genome


parser = argparse.ArgumentParser(description='Mutate a reference genome'
                                             'Usage: python mutate.py <reference_genome> r <rate_of_mutation> o '
                                             '<output_file>')
parser.add_argument('reference_genome_file', type=str, help='Reference genome')
parser.add_argument('-r', '--rate_of_mutation', type=str, help='Rate of mutation (%)')
parser.add_argument('-o', '--output_file', type=str, help='Output file')


def main():
    args = parser.parse_args()
    if not args.reference_genome_file:
        print("Error: Reference genome is required")
        sys.exit(1)

    print(args.reference_genome_file)

    rate_of_mutation = args.rate_of_mutation
    output_file = args.output_file
    if rate_of_mutation is None:
        rate_of_mutation = 0.01
    if output_file is None:
        ref_len_d = len(args.reference_genome_file.split("."))
        suff = args.reference_genome_file.split(".")[ref_len_d - 1]
        name = args.reference_genome_file.split(".")[ref_len_d - 2]

        output_file = name + "_mutated." + suff
    result = read_reference_genome(args.reference_genome_file)
    name = result[0]
    genome = result[1]

    positions = []

    mutations = round(len(genome) * float(rate_of_mutation))
    print("Reference genome has", len(genome), "nucleotides")
    print("Mutating", mutations, "positions")
    while mutations > 0:
        position = random.randint(0, len(genome) - 1)
        positions.append(position)
        genome = genome[:position] + random.choice('ACGT') + genome[position + 1:]
        mutations -= 1

    with open(output_file, 'w') as file:
        file.write(name)
        file.write(genome)
        print("Mutated genome has been written to", output_file)

    sorted(positions)

    with open("mutations.txt", 'w') as file:
        for position in positions:
            file.write(str(position) + "\n")
        print("Mutations have been written to mutations.txt")


if __name__ == "__main__":
    main()
