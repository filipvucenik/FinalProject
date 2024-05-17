import sys


def main():
    if len(sys.argv) < 3:
        print("Usage: python bas_file_generator.py <maternal_genome> <paternal_genome>")
        sys.exit(1)

    maternal_genome_input = sys.argv[1]
    paternal_genome_input = sys.argv[2]

    maternal_name = ""
    paternal_name = ""
    maternal_genome = ""
    paternal_genome = ""
    homozygous_regions = []

    with open(maternal_genome_input, 'r') as maternal_genome_file:
        maternal_genome_lines = maternal_genome_file.readlines()
        for line in maternal_genome_lines:
            if line.startswith('>'):
                maternal_name = line.split(" ")[0][1:] + "_maternal"
                continue
            maternal_genome += line.strip()

    with open(paternal_genome_input, 'r') as paternal_genome_file:
        paternal_genome_lines = paternal_genome_file.readlines()
        for line in paternal_genome_lines:
            if line.startswith('>'):
                paternal_name = line.split(" ")[0][1:] + "_paternal"
                continue
            paternal_genome += line.strip()

    start_ind = 0
    counter = 0
    mutation_flag = False

    for j in range(0, 10):
        if maternal_genome[j] != paternal_genome_lines[j]:
            mutation_flag = True
        if mutation_flag:
            counter += 1

    for j in range(10, len(maternal_genome)):
        if mutation_flag and counter < 10:
            counter += 1
        elif counter == 10:
            counter = 0
            start_ind = j
            mutation_flag = False

        if maternal_genome[j] != paternal_genome[j]:
            if mutation_flag:
                counter = 0
            else:
                homozygous_regions.append((start_ind, j - 11))
                mutation_flag = True

    with open("homozygosity.bed", "w") as output_file:
        for region in homozygous_regions:
            output_file.write(
                maternal_name + " " + str(region[0]) + " " + str(region[1]) +
                " + " + paternal_name + " " + str(region[0]) + " " + str(region[1]) + "\n")


if __name__ == "__main__":
    main()
