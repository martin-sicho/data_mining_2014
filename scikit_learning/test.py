from sklearn.naive_bayes import GaussianNB
from sklearn import cross_validation
from sklearn import datasets

iris = datasets.load_iris()
x = iris.data
y = iris.target

x_train, x_test, y_train, y_test = cross_validation.train_test_split(x, y, test_size=0.3)

classifier = GaussianNB()
classifier.fit(x_train, y_train)
y_predictions = classifier.predict(x_test)

print y_predictions, y_test

