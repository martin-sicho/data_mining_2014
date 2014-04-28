import sys
import os
import pickle

from datageneration import chembl, dud, fingerprinter
from datageneration.params import *
from modeling import classification


def main(args):
    """
        1. load ChEMBL data for molecules with IC50 (in nM) less than or equal to IC_50_THRESHOLD and convert them to RDKit molecules
            - all molecules are saved to tha data/ directory
            - if the pickle file already exists, no action is taken unless RELOAD_CHEMBL_DATA = True
        2. do the same with decoys from DUD
        3. merge both datasets
        4. compute fingerprints
        5. train Naive Bayes Classifier
        6. use SVM to predict IC50 values of active molecules
    """

    # load actives from ChEMBL
    actives = {}
    if not os.path.exists(DATA_FOLDER_PATH):
        os.mkdir(DATA_FOLDER_PATH)
    actives_file = [x for x in os.listdir(DATA_FOLDER_PATH) if x.startswith('actives_chembl') and x.endswith('.p')]
    if not actives_file or RELOAD_DATA:
        actives = chembl.loadChEMBLData(ACCESSION, IC_50_THRESHOLD, DATA_FOLDER_PATH)
        chembl.computeConsensualIC50(actives, DATA_FOLDER_PATH)
        chembl.appendRDKitMols(actives, DATA_FOLDER_PATH)
    else:
        actives = pickle.load(open(DATA_FOLDER_PATH + actives_file[0], 'rb'))

    # load decoys downloaded from DUD
    decoys = {}
    if os.path.exists(DECOYS_SDF_FILE_PATH[:-4] + ".p") and not RELOAD_DATA:
        decoys = pickle.load(open(DECOYS_SDF_FILE_PATH[:-4] + ".p", 'rb'))
    else:
        decoys = dud.getDecoys(DECOYS_SDF_FILE_PATH)

    # merge both data sets
    actives.update(decoys)
    compounds_all = actives

    # compute Morgan fingerprints
    if os.path.exists(MERGED_DATASET_PATH) and not RELOAD_DATA:
        print "Loading previously created dataset..."
        compounds_all = pickle.load(open(MERGED_DATASET_PATH, 'rb'))
    else:
        fingerprinter.appendMorganFingerprints(compounds_all, MORGAN_RADIUS)

    # train and cross-validate multiple Naive Bayes Classifiers
    classification_results = dict()
    if not os.path.exists(RESULTS_SAVE_FILE_PATH):
        classification_results = classification.naiveBayesClassification(compounds_all)
        print "Saving results..."
        pickle.dump(classification_results, open(RESULTS_SAVE_FILE_PATH, 'wb'))
        print "Finished analysis."
    else:
        print "Loading previous results..."
        classification_results = pickle.load(open(RESULTS_SAVE_FILE_PATH, 'rb'))

    # have fun with the classification results
    classification.playWithResults(classification_results)

    # SVM part
    

if __name__ == '__main__':
    main(sys.argv)
