# main parameters
ACCESSION = 'P03372'
IC_50_THRESHOLD = 30
DATA_FOLDER = "data/"
RELOAD_DATA = False
DECOYS_FILE_PATH = DATA_FOLDER + "decoys_DUD_structures.sdf"
MORGAN_RADIUS = 7

# datageneration parameters
XML_FILENAME = "actives_chembl.xml"
PICKLE_FILENAME = XML_FILENAME[:-4] + ".p"
PICKLE_PATH_ALL = DATA_FOLDER + "complete_dataset.p"
