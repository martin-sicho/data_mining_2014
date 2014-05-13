import sys
import os
import pickle
from rdkit import Chem
from rdkit.Chem import PropertyMol

from datamanipulation import chembl, dud, fingerprinter, utilities
from datamanipulation.params import *
from modeling import classification, regression

def trainModels():
    """
        Uses datamanipulation and modeling modules to train and validate models for a specific target.
        Use datamanipulation.params and modeling.params to set training and validation options.

        This method roughly proceeds as follows:
            1. load ChEMBL data for molecules with IC50 (in nM) less than or equal to IC_50_THRESHOLD and convert them to RDKit molecules
                - all molecules are saved to tha data/ directory
                - if the pickle file already exists, no action is taken unless RELOAD_CHEMBL_DATA = True
            2. do the same with decoys from DUD
            3. merge both datasets
            4. compute fingerprints
            5. train Naive Bayes Classifier
            6. cluster active molecules -> pick a representative compound of each cluster for training
            7. use SVM to predict IC50 values of active molecules (SV Regression)
            8. test regression model on the remaining molecules from the clusters
    """

    # load actives from ChEMBL
    actives = {}
    if not os.path.exists(DATA_FOLDER_PATH):
        os.mkdir(DATA_FOLDER_PATH)
    actives_file = [x for x in os.listdir(DATA_FOLDER_PATH) if x.startswith('actives_chembl') and x.endswith('.p')]
    if not actives_file or RELOAD_DATA and not USE_DOWNLOADED_STRUCTS:
        actives = chembl.loadChEMBLData(ACCESSION, IC_50_THRESHOLD, DATA_FOLDER_PATH)
    else:
        actives = pickle.load(open(DATA_FOLDER_PATH + actives_file[0], 'rb'))

    if not actives_file or RELOAD_DATA and not USE_DOWNLOADED_STRUCTS:
        chembl.computeConsensualIC50(actives, DATA_FOLDER_PATH)
        chembl.appendRDKitMols(actives, DATA_FOLDER_PATH)

    # load decoys downloaded from DUD
    decoys = {}
    if os.path.exists(DECOYS_SDF_FILE_PATH[:-4] + ".p"):
        decoys = pickle.load(open(DECOYS_SDF_FILE_PATH[:-4] + ".p", 'rb'))
    else:
        if os.path.exists(DECOYS_SDF_FILE_PATH):
            decoys = dud.getDecoys(DECOYS_SDF_FILE_PATH)
        else:
            print "Decoys not found in: " + DECOYS_SDF_FILE_PATH
            print "Make sure you set the right path."
            exit()

    # merge both data sets
    compounds_all = {}
    compounds_all.update(actives)
    compounds_all.update(decoys)

    # compute Morgan fingerprints
    if os.path.exists(MERGED_DATASET_PATH) and not RELOAD_DATA:
        print "Loading previously created dataset..."
        compounds_all = pickle.load(open(MERGED_DATASET_PATH, 'rb'))
    else:
        fingerprinter.appendMorganFingerprints(compounds_all)

    actives = { cmpndid : compounds_all[cmpndid] for cmpndid in compounds_all.keys() if compounds_all[cmpndid]['active']}
    pickle.dump(actives, open(ACTIVES_DUMP, 'wb'))
    decoys = { cmpndid : compounds_all[cmpndid] for cmpndid in compounds_all.keys() if not compounds_all[cmpndid]['active']}

    # train and cross-validate multiple Naive Bayes Classifiers
    classification_results = dict()
    if not os.path.exists(CLASS_RESULTS_SAVE_FILE_PATH) or RELOAD_DATA:
        classification_results = classification.naiveBayesClassifierTraining(compounds_all)
        print "Saving results..."
        pickle.dump(classification_results, open(CLASS_RESULTS_SAVE_FILE_PATH, 'wb'))
        print "Finished analysis."
    else:
        print "Loading previous results..."
        classification_results = pickle.load(open(CLASS_RESULTS_SAVE_FILE_PATH, 'rb'))

    # have fun with the classification results
    print "# CLASSIFICATION STATISTICS #"
    classification.playWithResults(classification_results)

    # cluster actives according to their similarity and keep only the diverse molecules
    actives_testset = dict()
    if CLUSTER:
        clusters = utilities.clusterMols(actives)
        actives_kept = dict()
        for cluster in clusters:
            actives_kept[cluster[0]] = actives[cluster[0]]
            remains = cluster[1:]
            actives_filtered_out = {chmblid : actives[chmblid] for chmblid in remains}
            actives_testset.update(actives_filtered_out)
        actives = actives_kept

    # estimate maximum distances between active molecules to set threshold for the application domain
    # distance_actives = regression.estimateDistanceThreshold(actives) # median of distances between two actives
    # min_distance_decoys, max_distance_decoys = regression.compareDistances(actives, decoys) # average min/max distance of closest/farthest decoy from any of the actives
    # print "median of distances between two actives: " + str(distance_actives)
    # print "average min/max distance of closest/farthest decoy from any of the actives: " + str(min_distance_decoys) + "/" + str(max_distance_decoys)

    # Support vector regression
    regression_results = dict()
    if not os.path.exists(REGRESS_RESULTS_SAVE_FILE_PATH) or RELOAD_DATA:
        regression_results = regression.supportVectorRegression(actives)
        pickle.dump(regression_results, open(REGRESS_RESULTS_SAVE_FILE_PATH, 'wb'))
    else:
        regression_results = pickle.load(open(REGRESS_RESULTS_SAVE_FILE_PATH, 'rb'))


    # do something with the regression results
    print "# REGRESSION STATISTICS #"
    regression.playWithResults(regression_results, decoys, actives_testset)

    return classification_results['final_model'], regression_results['final_model']

def predict(classmodel, regressmodel, molfile_path):
    print "Starting predictions for: " + molfile_path
    suppl = Chem.SDMolSupplier(molfile_path)
    mols = dict()
    for mol in suppl:
        pmol = PropertyMol.PropertyMol(mol)
        mols[pmol.GetProp("_Name")] = {"RDKit" : pmol}
    fingerprinter.appendMorganFingerprints(mols, dump=None)
    actives = pickle.load(open(ACTIVES_DUMP, 'rb'))

    found_sth = False
    for mol in mols:
        prediction = classmodel.predict(mols[mol]['fingerprint'])
        fingerprints_actives = utilities.getFingerprintList(actives)[0]
        min_distance = utilities.getMolDistFromSet(mols[mol]['fingerprint'], fingerprints_actives)[0]
        if min_distance <= APPLICABILITY_DOMAIN_DISTANCE_THRESHOLD and prediction[0]:
            print  mol + " is active"
            print "Predicted pIC50: " + str(regressmodel.predict(mols[mol]['fingerprint'])[0])
            found_sth = True
    if not found_sth:
        print "None of the molecules within the specified set were found to be active."


def main(args):
    # build the models
    classmodel, regressmodel = trainModels()

    # test the activity prediction on the set of downloaded decoys
    molfile = DECOYS_SDF_FILE_PATH
    predict(classmodel, regressmodel, molfile)

if __name__ == '__main__':
    main(sys.argv)
