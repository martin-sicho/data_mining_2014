import sys, os, pickle
from datageneration import chembl, dud, fingerprinter
from params.params import *

def main(args):
    """
        1. load ChEMBL data for molecules with IC50 (in nM) less than or equal to IC_50_THRESHOLD and convert them to RDKit molecules
            - all molecules are saved to tha data/ directory
            - if the pickle file already exists, no action is taken unless RELOAD_CHEMBL_DATA = True
        2. do the same with decoys from DUD
        3. merge both datasets
        4. compute fingerprints
    """

    # load actives from ChEMBL
    actives = {}
    actives_file = [x for x in os.listdir(DATA_FOLDER) if x.startswith('actives_chembl') and x.endswith('.p')]
    if not actives_file or RELOAD_DATA:
        actives = chembl.loadChEMBLData(ACCESSION, IC_50_THRESHOLD, DATA_FOLDER)
        chembl.computeConsensualIC50(actives, DATA_FOLDER)
        chembl.appendRDKitMols(actives, DATA_FOLDER)
    else:
        actives = pickle.load(open(DATA_FOLDER + actives_file[0], 'rb'))

    # load decoys downloaded from DUD
    decoys = {}
    if os.path.exists(DECOYS_FILE_PATH[:-4] + ".p") and not RELOAD_DATA:
        decoys = pickle.load(open(DECOYS_FILE_PATH[:-4] + ".p", 'rb'))
    else:
        decoys = dud.getDecoys(DECOYS_FILE_PATH)

    # merge both data sets
    actives.update(decoys)
    compounds_all = actives

    # compute Morgan fingerprints
    if os.path.exists(PICKLE_PATH_ALL) and not RELOAD_DATA:
        compounds_all = pickle.load(open(PICKLE_PATH_ALL, 'rb'))
    else:
        fingerprinter.appendMorganFingerprints(compounds_all, MORGAN_RADIUS)


if __name__ == '__main__':
    main(sys.argv)
