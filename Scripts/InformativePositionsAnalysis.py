from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import pandas as pd
import sys
from sklearn.tree import DecisionTreeClassifier
import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import load_iris
from sklearn.inspection import DecisionBoundaryDisplay


def main():
    if len(sys.argv) < 2:
        print("Usage: python InformativePositionAnalysis.py <data_file>")
        sys.exit(1)

    data = pd.read_csv(sys.argv[1])
    X = data[['Informative_Positions', 'Informative_Positions_Matches']]
    print(X.head(4))
    y = data['True_Predictions']

    X_train, X_test, y_train, y_test = train_test_split(*(X, y), test_size=0.2, random_state=0)
    model = LogisticRegression()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_pred_train = model.predict(X_train)
    x0 = X[y == 0]
    print (len(x0))
    x1 = X[y == 1]
    print(len(x1))


    plt.scatter(x0['Informative_Positions'], x0['Informative_Positions_Matches'], color='red', alpha=0.5)
    #plt.scatter(x1['Informative_Positions'], x1['Informative_Positions_Matches'], color='blue', alpha=0.5)
    plt.xlabel('Informativne Pozicije')
    plt.ylabel('Informativne Pozicije Poklapanja')
    plt.legend()
    plt.show()

    print("Logistic Regression")
    count = 0
    print(classification_report(y_test, y_pred))
    print(confusion_matrix(y_test, y_pred))
    print(accuracy_score(y_test, y_pred))

    print(classification_report(y_train, y_pred_train))
    print(confusion_matrix(y_train, y_pred_train))

    print(data)

    classifier = LogisticRegression().fit(X_train, y_train)
    disp = DecisionBoundaryDisplay.from_estimator(
        classifier, X_train, response_method="predict",
    )
    disp.ax_.scatter(X_train[:, 0], X_train[:, 1], c=y_train, edgecolor="k")
    plt.show()


    print("Coefficients: ", model.coef_)
    print("Intercept: ", model.intercept_)
    print("SVM linear kernel")
    model = SVC(class_weight='balanced', kernel='rbf', max_iter=1000, C=0.01)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_pred_train = model.predict(X_train)


    print(classification_report(y_test, y_pred))
    print(confusion_matrix(y_test, y_pred))

    print(classification_report(y_train, y_pred_train))
    print(confusion_matrix(y_train, y_pred_train))

    print(accuracy_score(y_test, y_pred))
    #print("Coefficients: ", model.coef_)
    print("Intercept: ", model.intercept_)


    print("Decision Tree")
    model = DecisionTreeClassifier()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_pred_train = model.predict(X_train)

    print(classification_report(y_test, y_pred))
    print(confusion_matrix(y_test, y_pred))

    print(classification_report(y_train, y_pred_train))
    print(confusion_matrix(y_train, y_pred_train))

    informative_sum_1 =0
    informative_sum_0 = 0

    #for i in range(len(data)):
    #    if data['True_Predictions'][i] == 1:
    #        informative_sum_1 += data['Informative_Positions_Matches'] / data['Informative_Positions'][i]
    #    else:
    #        informative_sum_0 += data['Informative_Positions_Matches'] / data['Informative_Positions'][i]

    #print("Informative sum 1: ", informative_sum_1 / len(data[data['True_Predictions'] == 1]))
    #print("Informative sum 0: ", informative_sum_0 / len(data[data['True_Predictions'] == 0]))


if __name__ == "__main__":
    main()
