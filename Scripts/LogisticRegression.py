from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
import pandas as pd
import sys
import matplotlib.pyplot as plt


def main():
    if len(sys.argv) < 2:
        print("Usage: python LogisticRegression.py <data_file>")
        sys.exit(1)
    data = pd.read_csv(sys.argv[1])
    X = data[['len1', 'len2', 'matches']]
    X['diff'] = abs(X['len1'] - X['len2'])
    print(X.head(4))
    y = data['match']
    X1 = pd.DataFrame()
    X1['matches'] = X['matches']
    X1['len2'] = X['len2']
    X_train, X_test, y_train, y_test = train_test_split(*(X1, y), test_size=0.2, random_state=0)
    model = LogisticRegression()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    x0 = X[y == 0]
    x1 = X[y == 1]

    plt.scatter(x0['matches'], x0['diff'], color='red')
    plt.scatter(x1['matches'], x1['diff'], color='blue')
    plt.xlabel('matches')
    plt.ylabel('diff')
    plt.legend()
    plt.show()

    count = 0
    for i, prediction in enumerate(y_pred):
        if prediction == 1:
            count += 1

    print("Actual 1: ", sum(y_test))
    print("Predicted 1: ", count)

    print(accuracy_score(y_test, y_pred))
    print("Coefficients: ", model.coef_)
    print("Intercept: ", model.intercept_)


if __name__ == "__main__":
    main()
