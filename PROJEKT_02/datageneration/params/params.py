# main parameters
ACCESSION = 'P03372'
IC_50_THRESHOLD = 30
DATA_FOLDER_PATH = "data/"
RELOAD_DATA = False
DECOYS_SDF_FILE_PATH = DATA_FOLDER_PATH + "decoys_DUD_structures.sdf"
MORGAN_RADIUS = 7

# datageneration parameters
XML_FILENAME = "actives_chembl.xml"
PICKLE_FILENAME = XML_FILENAME[:-4] + ".p"
MERGED_DATASET_PATH = DATA_FOLDER_PATH + "complete_dataset.p"

# results save file
RESULTS_SAVE_FILE_PATH = DATA_FOLDER_PATH + "results.p"
