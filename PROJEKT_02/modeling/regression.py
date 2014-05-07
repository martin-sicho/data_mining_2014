import numpy
from sklearn import svm, cross_validation
from params import *
import matplotlib.pyplot as plt
from datageneration import utilities

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
    kfold_xv_strat = cross_validation.KFold(len(activity_data), SVR_FOLDS, indices=False)
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

def estimateDistanceThreshold(mols):
    fps = utilities.getFingerprintList(mols)[0]
    dists = utilities.generateDistMatrix(fps)[0]
    return numpy.median(dists)

def compareDistances(actives, decoys):
    actives_fps = utilities.getFingerprintList(actives)[0]
    decoys_fps = utilities.getFingerprintList(decoys)[0]
    distances_min = []
    distances_max = []
    for decoy in decoys_fps:
        dists = utilities.getMolDistFromSet(decoy, actives_fps)
        distances_min.append(dists[0])
        distances_max.append(dists[1])
    return numpy.mean(distances_min), numpy.mean(distances_max)

def playWithResults(results, decoys, actives_test_set):
    actives_test_set_fps, keys = utilities.getFingerprintList(actives_test_set)
    actives_test_set_pic50 = [-1.0*numpy.log10(actives_test_set[chmblid]['ic50'] / 10E9) for chmblid in keys]
    actives_test_set_fps = numpy.asarray(actives_test_set_fps)
    actives_test_set_pic50 = numpy.asarray(actives_test_set_pic50)

    #keys = decoys.keys()
    #decoys_fingerprint_data = [decoys[cmpnd_id]['fingerprint'] for cmpnd_id in keys]
    #decoys_fingerprint_data = numpy.asarray(decoys_fingerprint_data)
    #zeros = [10.75 for x in keys]

    # best model from cross-validation
    best_model_idx = results['scores'].index(max(results['scores']))
    print "Best score: " + str(results['scores'][best_model_idx])
    print "Average score: " + str(numpy.mean(results['scores']))
    #predicted_best = results['predicted_values'][best_model_idx]
    #true_best = results['true_values'][best_model_idx]

    # final model on training set
    final_model = results['final_model']
    predicted_train = final_model.predict(results['fingerprint_data'])
    print "Score of final model on the molecules from the training set: " + str(final_model.score(results['fingerprint_data'], results['activity_data']))

    #predictions_all = final_model.predict(results['fingerprint_data_validation_set'])
    #predictions_decoys = final_model.predict(decoys_fingerprint_data)
    predictions_test_set = final_model.predict(actives_test_set_fps)
    #print "Score of final model on the validation set: " + str(final_model.score(results['fingerprint_data_validation_set'], results['activity_data_validation_set']))
    print "Score of final model on the molecules filtered out during clustering: " + str(final_model.score(actives_test_set_fps, actives_test_set_pic50))

    span = (8, 11)
    plt.plot((span[0],span[1]), (span[0],span[1]), linestyle='--')
    #plt.plot(results['activity_data_validation_set'], predictions_all, marker='o', linestyle='None', label="Validation set performance")
    #plt.plot(true_best, predicted_best, marker='+', linestyle='None', label="Performance of the best model in the particular X-validation step")
    #plt.plot(zeros, predictions_decoys, marker='o', linestyle='None', label="decoys")
    plt.plot(actives_test_set_pic50, predictions_test_set, marker='o', linestyle='None', label="Performance on the filtered out molecules during clustering")
    plt.plot(results['activity_data'], predicted_train, marker='o', linestyle='None', label="Performance on the training set")
    plt.xlabel('True values')
    plt.ylabel('Predicted values')
    plt.ylim(span)
    plt.xlim(span)
    plt.legend()
    plt.show()