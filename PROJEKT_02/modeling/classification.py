import numpy
from params import *
from sklearn.naive_bayes import GaussianNB
from sklearn import cross_validation
from sklearn.metrics import confusion_matrix, roc_curve, auc

def naiveBayesClassifierTraining(compounds_all):
    print "Building naive Bayes classifier (" + str(NB_FOLDS) + "-fold cross-validation)..."
    # get the data
    keys = compounds_all.keys()
    fingerprint_data = [compounds_all[cmpnd_id]['fingerprint'] for cmpnd_id in keys]
    fingerprint_data = numpy.asarray(fingerprint_data)
    activity_data = [compounds_all[cmpnd_id]['active'] for cmpnd_id in keys]
    activity_data = numpy.asarray(activity_data)

    # perform K-fold cross-validation
    classifier = GaussianNB()
    kfold_xv_strat = cross_validation.StratifiedKFold(activity_data, NB_FOLDS, indices=False)
    confusion_matrices = []
    probabilities = []
    scores = []
    models = []
    true_activities = []
    aucs = []
    for train, test in kfold_xv_strat:
        fingerprint_data_train = fingerprint_data[train]
        fingerprint_data_test = fingerprint_data[test]
        activity_data_train = activity_data[train]
        activity_data_test = activity_data[test]

        # model building
        classifier.fit(fingerprint_data_train, activity_data_train)

        # testing
        activity_data_predictions = classifier.predict(fingerprint_data_test)
        models.append(classifier)

        probability_estimates = classifier.predict_proba(fingerprint_data_test)
        probabilities.append(probability_estimates)

        scores.append(classifier.score(fingerprint_data_test, activity_data_test))

        activity_confusion_matrix = confusion_matrix(activity_data_test, activity_data_predictions)
        confusion_matrices.append(activity_confusion_matrix)

        true_activities.append(activity_data_test)

        # ROC curves
        fpr, tpr, thresholds = roc_curve(activity_data_test, probability_estimates[:, 1])
        aucs.append(auc(fpr, tpr))
    classifier.fit(fingerprint_data, activity_data)
    print "Done."
    return {
        'confusion_matrices' : confusion_matrices
        , 'probabilities' : probabilities
        , 'scores' : scores
        , 'models' : models
        , 'true_activity_data' : true_activities
        , 'AUCs' : aucs
        , 'fingerprint_data' : fingerprint_data
        , 'activity_data' : activity_data
        , 'final_model' : classifier
    }

def playWithResults(classification_results):
    best_model_idx = classification_results['scores'].index(max(classification_results['scores']))
    print "Average AUC: " + str(numpy.mean(classification_results['AUCs']))
    print "BEST MODEL DETAILS:"
    print "Confusion matrix: "
    print classification_results['confusion_matrices'][best_model_idx]
    print "AUC: " + str(classification_results['AUCs'][best_model_idx])
    # for idx, probabilities in enumerate(classification_results['probabilities'][best_model_idx]):
    #     print "Active probability: " + str(probabilities[1])
    #     print "True Activity: " + str(classification_results['true_activity_data'][best_model_idx][idx])
