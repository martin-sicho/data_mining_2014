import sys, os, pickle
from datamanipulation import chembl, dud

ACCESSION = 'P03372'
IC_50_THRESHOLD = 30
DATA_FOLDER = "data/"
RELOAD_DATA = False
DECOYS_FILE_PATH = DATA_FOLDER + "decoys_DUD_structures.sdf"

def main(args):
    """
        1. load ChEMBL data for molecules with IC50 (in nM) less than or equal to IC_50_THRESHOLD
            - all molecules are saved to tha data/ directory
            - if the pickle file already exists, no action is taken unless RELOAD_CHEMBL_DATA = True
    """
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
        pickle.load(open(DECOYS_FILE_PATH[:-4] + ".p", 'rb'))
    else:
        decoys = dud.getDecoys(DECOYS_FILE_PATH)


if __name__ == '__main__':
    main(sys.argv)
