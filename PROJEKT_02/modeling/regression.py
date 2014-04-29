import numpy
from sklearn import svm, cross_validation
from params import *
import matplotlib.pyplot as plt

def supportVectorRegression(actives):
    print "Performing Support Vector Regression..."
    # get the data
    keys = actives.keys()
    fingerprint_data = [actives[cmpnd_id]['fingerprint'] for cmpnd_id in keys]
    fingerprint_data = numpy.asarray(fingerprint_data)
    activity_data = [-1.0*numpy.log10(actives[cmpnd_id]['ic50'] / 10E9) for cmpnd_id in keys]
    activity_data = numpy.asarray(activity_data)

    # cross-validate multiple regression models
    fingerprint_data, fingerprint_data_validation_set, activity_data, activity_data_validation_set = cross_validation.train_test_split(fingerprint_data, activity_data, test_size=VALIDATION_SET_SIZE)
    kfold_xv_strat = cross_validation.KFold(len(activity_data), n_folds=SVR_FOLDS, indices=False)
    clf = svm.SVR(
                        C=C
                      , epsilon=EPSILON
                      , kernel=KERNEL
                      , probability=PROBABILITY
                      , verbose=VERBOSE
                      , coef0=COEF0
                      , degree=DEGREE
                      , gamma=GAMMA
                )
    predicted_values = []
    true_values = []
    scores = []
    models = []
    for train, test in kfold_xv_strat:
        fingerprint_data_train = fingerprint_data[train]
        fingerprint_data_test = fingerprint_data[test]
        activity_data_train = activity_data[train]
        activity_data_test = activity_data[test]

        clf.fit(fingerprint_data_train, activity_data_train)
        models.append(clf)
        predicted_values.append(clf.predict(fingerprint_data_test))
        true_values.append(activity_data_test)
        scores.append(clf.score(fingerprint_data_test, activity_data_test))
    return {
        'predicted_values' : predicted_values
        , 'true_values' : true_values
        , 'scores' : scores
        , 'models' : models
        , 'fingerprint_data' : fingerprint_data
        , 'activity_data' : activity_data
        , 'fingerprint_data_validation_set' : fingerprint_data_validation_set
        , 'activity_data_validation_set' : activity_data_validation_set
        , 'final_model' : clf.fit(fingerprint_data, activity_data)
    }

def playWithResults(results):
    best_model_idx = results['scores'].index(max(results['scores']))
    print "Best score: " + str(results['scores'][best_model_idx])
    print "Average score: " + str(numpy.mean(results['scores']))
    predicted_best = results['predicted_values'][best_model_idx]
    true_best = results['true_values'][best_model_idx]
    final_model = results['final_model']
    predictions_all = final_model.predict(results['fingerprint_data_validation_set'])
    print "Score of final model on the validation set: " + str(final_model.score(results['fingerprint_data_validation_set'], results['activity_data_validation_set']))
    xspan = (8,11)
    plt.plot((xspan[0],xspan[1]), (xspan[0],xspan[1]), linestyle='--')
    plt.plot(results['activity_data_validation_set'], predictions_all, marker='o', linestyle='None', label="Validation set performance")
    plt.plot(true_best, predicted_best, marker='+', linestyle='None', label="Predictions of the best model in the particular X-validation step")
    plt.xlabel('True values')
    plt.ylabel('Predicted values')
    plt.ylim(xspan)
    plt.xlim(xspan)
    plt.legend()
    plt.show()