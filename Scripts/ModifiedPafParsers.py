from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
import pandas as pd
import sys
import matplotlib.pyplot as plt



def main():
    if len(sys.argv) != 3:
        print("Usage: python InformativePositionsAnalysis.py <modified PAF file> <h1>")
        sys.exit(1)


    data_file = sys.argv[1]
    informative_positions = []
    informative_positions_matches = []
    true_predictions = []
    h1 = sys.argv[2]
    with open(data_file) as f:
        for line in f:
            line = line.strip().split("\t")
            ip_split = line[13].split("/")
            ipm = int(ip_split[0])
            ip = int(ip_split[1])
            informative_positions.append(ip)
            informative_positions_matches.append(ipm)
            r1P = line[0].__contains__(h1)
            r2P = line[5].__contains__(h1)
            if r1P and r2P:
                true_predictions.append(0)
            else:
                true_predictions.append(1)





if __name__ == "__main__":
    main()
