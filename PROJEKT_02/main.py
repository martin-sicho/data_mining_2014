import sys
import datamanipulation

ACCESSION = 'P03372'
IC_50_THRESHOLD = 30
DATA_FOLDER = "data/"

def main(args):
    datamanipulation.loadChEMBLData(ACCESSION, IC_50_THRESHOLD, DATA_FOLDER)

if __name__ == '__main__':
    main(sys.argv)
