import sys
from sklearn.metrics import confusion_matrix
from matplotlib import pyplot as plt

def main():
    if len(sys.argv) != 3:
        print("Usage: python ModifiedPafParsers.py <modified PAF file> <h1>")
        sys.exit(1)

    data_file = sys.argv[1]
    informative_positions = []
    informative_positions_matches = []
    true_predictions = []
    h1 = sys.argv[2]
    paf_predicted = []

    correct_counter = 0
    all_counter = 0

    counter0 = 0
    with open(data_file) as f:
        for line in f:
            line = line.strip().split("\t")
            ip_split = line[13].split("/")
            ipm = int(ip_split[0])
            ip = int(ip_split[1])
            if ipm > 3000:
                print(line)
            informative_positions.append(ip)
            informative_positions_matches.append(ipm)
            r1P = line[0].__contains__(h1)
            r2P = line[5].__contains__(h1)
            paf_predicted.append(int(line[12]))
            all_counter += 1
            if r1P == r2P:
                true_predictions.append(0)
                counter0 +=1
                if int(line[12]) == 0:
                    correct_counter += 1
            else:
                true_predictions.append(1)
                if int(line[12]) == 1:
                    correct_counter += 1

    output_file = data_file + ".csv"
    with open(output_file, "w") as f:
        f.write("Informative_Positions,Informative_Positions_Matches,True_Predictions\n")
        for i in range(len(informative_positions)):
            if informative_positions[i] == 0:
                continue
            f.write(str(informative_positions[i]) + "," + str(informative_positions_matches[i]) + "," + str(
                true_predictions[i]) + "\n")

    print("Paf predicted confusion matrix")
    print(confusion_matrix(true_predictions, paf_predicted))

    plt.matshow(confusion_matrix(true_predictions, paf_predicted))
    plt.show()

    print("Accuracy: " + str(correct_counter / all_counter))
    print("Counter 0:", str(counter0 / all_counter))
    print("Number counter 0: " + str(counter0))


if __name__ == "__main__":
    main()
