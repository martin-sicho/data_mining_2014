# main method parameters
ACCESSION = 'P03372'
DATA_FOLDER_PATH = "data/"
DECOYS_SDF_FILE_PATH = DATA_FOLDER_PATH + "decoys_DUD_structures.sdf"
RESULTS_SAVE_FILE_PATH = DATA_FOLDER_PATH + "results.p"
IC_50_THRESHOLD = 30
RELOAD_DATA = True
MORGAN_RADIUS = 3
USE_DOWNLOADED_STRUCTS = False

# datageneration parameters
XML_FILENAME = "actives_chembl.xml"
PICKLE_FILENAME = XML_FILENAME[:-4] + ".p"
MERGED_DATASET_PATH = DATA_FOLDER_PATH + "complete_dataset.p"
