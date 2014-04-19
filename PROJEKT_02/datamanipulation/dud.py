import pickle
from rdkit import Chem
from rdkit.Chem import PropertyMol

def getDecoys(decoys_path):
    print "Converting decoys..."
    decoys = {}
    suppl = Chem.SDMolSupplier(decoys_path)
    for mol in suppl:
        pmol = PropertyMol.PropertyMol(mol)
        decoys[pmol.GetProp("_Name")] = {"RDKit" : pmol, 'active': False}
    pickle.dump(decoys, open(decoys_path[:-4] + ".p", "wb"))
    print "Saved data for " + str(len(decoys)) + " decoys."
    return decoys
