import xml.etree.ElementTree as etree
from xml.dom import minidom
import requests
import json
import re
import os
import pickle
import logging
import sys
import numpy

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PropertyMol
from standardise import standardise

from utilities import *
from params import *


def loadChEMBLData(uniprot_accesion, IC_50_threshold, data_folder):
    filename = XML_FILENAME
    xml_path = data_folder + filename
    if os.path.exists(xml_path) and os.path.exists(PICKLE_FILENAME):
        return pickle.load(open(PICKLE_FILENAME, "rb"))
    if not(os.path.exists(data_folder)):
        os.makedirs(data_folder)
    print "Getting bioactivity data from ChEMBL database..."
    # 1. Use UniProt accession to get target details
    target_data = json.loads(requests.get("https://www.ebi.ac.uk/chemblws/targets/uniprot/%s.json" % uniprot_accesion).content)

    # 2. Get all bioactivties for target CHEMBL_ID
    bioactivity_data_all = json.loads(requests.get("https://www.ebi.ac.uk/chemblws/targets/%s/bioactivities.json" % target_data['target']['chemblId']).content)

    # 3. Get compounds with high binding affinity (IC50 < threshold)
    root = etree.Element("compounds")
    compounds = {}
    for bioactivity in [record for record in bioactivity_data_all['bioactivities'] if re.search('IC50', record['bioactivity_type']) and record['units'] == "nM" and isFloat(record['value']) and float(record['value']) < IC_50_threshold]:
        chemblid = bioactivity['ingredient_cmpd_chemblid']
        print "Compound CHEMBLID: %s" % chemblid
        cmpd_data = etree.fromstring(requests.get("https://www.ebi.ac.uk/chemblws/compounds/%s" % chemblid).content)

        # add IC50 values
        ic50 = etree.Element("ic50")
        ic50.text = bioactivity['value']
        cmpd_data.append(ic50)

        # add activity stamp
        active = etree.Element("active")
        active.text = "true"
        cmpd_data.append(active)

        # add the node to the XML tree
        root.append(cmpd_data)

        # get SMILES for each molecule
        smiles = cmpd_data.findall("./smiles")[0].text
        compounds[chemblid] = {'smiles' : smiles}
        print "Smiles: %s" % smiles
    #tree = etree.ElementTree(root)
    #tree.write(folder + filename, encoding="utf8")
    with open(xml_path, "w") as outfile:
        outfile.write(minidom.parseString(etree.tostring(root)).toprettyxml(encoding="UTF-8").encode("utf8"))
    pickle.dump(compounds, open(data_folder + PICKLE_FILENAME, "wb"))
    print "Done. Downloaded data for " + str(len(compounds)) + " active compounds."
    return compounds

def computeConsensualIC50(compounds, data_folder):
    print "Computing consensual IC50 values..."
    tree = etree.parse(data_folder + XML_FILENAME)
    data = tree.getroot()
    for chemblid in compounds.keys():
        ic50_data = data.findall(".//compound[chemblId='%s']/ic50" % chemblid)
        ic50_values = []
        for ic50 in ic50_data:
            ic50_values.append(float(ic50.text))
        median = numpy.median(ic50_values)
        # mean = numpy.mean(ic50_values)
        # minimum = min(ic50_values)
        # maximum = max(ic50_values)
        compounds[chemblid]['ic50'] = median
    pickle.dump(compounds, open(data_folder + PICKLE_FILENAME, "wb"))
    return compounds

def appendRDKitMols(compounds, data_folder):
    print "Getting rdkit molecules..."
    failed_counter = 0
    for chemblid in compounds.keys():
        mol = Chem.MolFromSmiles(compounds[chemblid]['smiles'])
        if mol == None:
            print "RDKit molecule generation failed for " + chemblid
            del compounds[chemblid]
            failed_counter+=1
            continue
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        if not('win' in sys.platform):
            standardised_mol = None
            try:
                standardised_mol = standardise.apply(mol)
            except standardise.StandardiseException as e:
                logging.warn(e.message)
            pmol = PropertyMol.PropertyMol(standardised_mol)
            pmol.SetProp("_Name", chemblid)
            compounds[chemblid]['RDKit'] = pmol
        else:
            pmol = PropertyMol.PropertyMol(mol)
            pmol.SetProp("_Name", chemblid)
            compounds[chemblid]['RDKit'] = pmol
        compounds[chemblid]['active'] = True
        mol_path = data_folder + "actives_chembl_structures/" + chemblid + ".mol"
        if not(os.path.exists(mol_path[:mol_path.rindex('/')])):
            os.makedirs(mol_path[:mol_path.rindex('/')])
        with open(mol_path, "w") as outfile:
            outfile.write(Chem.MolToMolBlock(pmol))
    pickle.dump(compounds, open(data_folder + PICKLE_FILENAME, "wb"))
    print "Done. Saved data for " + str(len(compounds)) + " active compounds (" + str(failed_counter) + " molecules failed)."