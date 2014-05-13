# main parameters
ACCESSION = 'P03372'
DATA_FOLDER_PATH = "data/"
DECOYS_SDF_FILE_PATH = DATA_FOLDER_PATH + "decoys_DUD_structures.sdf"
CLASS_RESULTS_SAVE_FILE_PATH = DATA_FOLDER_PATH + "results.p"
REGRESS_RESULTS_SAVE_FILE_PATH = DATA_FOLDER_PATH + "results_regr.p"
ACTIVES_DUMP = DATA_FOLDER_PATH + "actives.p"
IC_50_THRESHOLD = 30
RELOAD_DATA = False
MORGAN_RADIUS = 3
MORGAN_LENGTH = 512
USE_DOWNLOADED_STRUCTS = True
CLUSTER = True
CLUSTERING_TANIMOTO_CUTOFF = 0.15
APPLICABILITY_DOMAIN_DISTANCE_THRESHOLD = CLUSTERING_TANIMOTO_CUTOFF

# datamanipulation parameters
XML_FILENAME = "actives_chembl.xml"
PICKLE_FILENAME = XML_FILENAME[:-4] + ".p"
MERGED_DATASET_PATH = DATA_FOLDER_PATH + "complete_dataset.p"