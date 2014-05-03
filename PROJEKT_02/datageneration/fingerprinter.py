import pickle
from params import *

from rdkit.Chem import AllChem


def appendMorganFingerprints(compounds):
    print "Computing Morgan Fingerprints..."
    for cmnd_id in compounds.keys():
        mol = compounds[cmnd_id]['RDKit']
        info = {}
        #fp = AllChem.GetMorganFingerprint(mol, radius, bitInfo=info)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_LENGTH, bitInfo=info)
        compounds[cmnd_id]['fingerprint'] = fp #list(fp.ToBitString())
        #compounds[cmnd_id]['fingerprint_info'] = info
    pickle.dump(compounds, open(MERGED_DATASET_PATH, "wb"))
    print "Done."