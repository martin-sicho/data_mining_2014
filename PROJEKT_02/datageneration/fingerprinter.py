import pickle
from rdkit.Chem import AllChem
from params.params import *

def appendMorganFingerprints(compounds, radius):
    print "Computing Morgan Fingerprints..."
    for cmnd_id in compounds.keys():
        mol = compounds[cmnd_id]['RDKit']
        info = {}
        #fp = AllChem.GetMorganFingerprint(mol, radius, bitInfo=info)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=1024, bitInfo=info)
        compounds[cmnd_id]['fingerprint'] = fp #list(fp.ToBitString())
        #compounds[cmnd_id]['fingerprint_info'] = info
    pickle.dump(compounds, open(MERGED_DATASET_PATH, "wb"))
    print "Done."