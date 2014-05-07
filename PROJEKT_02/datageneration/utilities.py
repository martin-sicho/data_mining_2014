from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from params import *

def isFloat(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

def getFingerprintList(mols):
    keys = mols.keys()
    return [mols[cmpnd_id]['fingerprint'] for cmpnd_id in keys], keys

def generateDistMatrix(fps):
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])
    return dists, nfps

def getMolDistFromSet(mol, molset):
    dists = [1 - x for x in DataStructs.BulkTanimotoSimilarity(mol, molset)]
    return min(dists), max(dists)

def clusterMols(mols,cutoff=CLUSTERING_TANIMOTO_CUTOFF):
    fps, keys = getFingerprintList(mols)

    # first generate the distance matrix:
    dists, nfps = generateDistMatrix(fps)

    # now cluster the data:
    cs = list(Butina.ClusterData(dists,nfps,cutoff,isDistData=True))

    # transform indices to dictionary keys
    for i, cluster in enumerate(cs):
        cluster = list(cluster)
        for idx, mol_idx in enumerate(cluster):
            cluster[idx] = keys[mol_idx]
        cs[i] = tuple(cluster)
    return tuple(cs)