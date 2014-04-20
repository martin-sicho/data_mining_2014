import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from params.params import *

def appendMorganFingerprints(compounds, radius):
    print "Computing Morgan Fingerprints..."
    for cmnd_id in compounds.keys():
        mol = compounds[cmnd_id]['RDKit']
        info = {}
        #fp = AllChem.GetMorganFingerprint(mol, radius, bitInfo=info)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=1024, bitInfo=info)
        compounds[cmnd_id]['fingerprint'] = fp
        compounds[cmnd_id]['fingerprint_info'] = info
    pickle.dump(compounds, open(PICKLE_PATH_ALL, "wb"))
    print "Done."