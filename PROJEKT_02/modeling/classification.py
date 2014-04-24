import numpy
from params import *
from sklearn.naive_bayes import GaussianNB
from sklearn import cross_validation
from sklearn.metrics import confusion_matrix

def naiveBayesClassification(compounds_all):
    print "Building naive Bayes classifier (" + str(FOLDS) + "-fold cross-validation)..."
    print
    # get the data
    keys = compounds_all.keys()
    fingerprint_data = [compounds_all[cmpnd_id]['fingerprint'] for cmpnd_id in keys]
    fingerprint_data = numpy.asarray(fingerprint_data)
    activity_data = [compounds_all[cmpnd_id]['active'] for cmpnd_id in keys]
    activity_data = numpy.asarray(activity_data)

    # perform K-fold cross-validation
    kfold_xv_strat = cross_validation.StratifiedKFold(activity_data, FOLDS, indices=False)
    confusion_matrices = []
    probabilities = []
    scores = []
    models = []
    true_activities = []
    for train, test in kfold_xv_strat:
        fingerprint_data_train = fingerprint_data[train]
        fingerprint_data_test = fingerprint_data[test]
        activity_data_train = activity_data[train]
        activity_data_test = activity_data[test]

        # model building
        classifier = GaussianNB()
        classifier.fit(fingerprint_data_train, activity_data_train)

        # testing
        activity_data_predictions = classifier.predict(fingerprint_data_test)
        models.append(classifier)
        probabilities.append(classifier.predict_proba(fingerprint_data_test))
        scores.append(classifier.score(fingerprint_data_test, activity_data_test))
        activity_confusion_matrix = confusion_matrix(activity_data_test, activity_data_predictions)
        confusion_matrices.append(activity_confusion_matrix)
        true_activities.append(activity_data_test)
    print "Done."
    return {
        'confusion_matrices' : confusion_matrices
        , 'probabilities' : probabilities
        , 'scores' : scores
        , 'models' : models
        , 'true_activity_data' : true_activities
    }
