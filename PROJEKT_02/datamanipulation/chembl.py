from xml.dom import minidom
import requests, json, re, os, pickle
import xml.etree.ElementTree as etree

def isFloat(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

def loadChEMBLData(uniprot_accesion, IC_50_threshold, folder):
    filename = "actives_chembl.xml"
    if os.path.exists(folder + filename):
        return
    if not(os.path.exists(folder)):
        os.makedirs(folder)
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
        compounds[chemblid] = smiles
        print "Smiles: %s" % smiles
    #tree = etree.ElementTree(root)
    #tree.write(folder + filename, encoding="utf8")
    with open(folder + filename, "w") as outfile:
        outfile.write(minidom.parseString(etree.tostring(root)).toprettyxml(encoding="UTF-8").encode("utf8"))
    pickle.dump(compounds, open(folder + filename[:-4] + ".p", "w"))



# # 4. Get assay details foe Ki actvity types
#
# print """
#
# # =========================================================
# # 4. Get assay details foe Ki actvity types
# # =========================================================
# """
#
# for bioactivity in [record for record in bioactivity_data['bioactivities'] if re.search('Ki', record['bioactivity_type'], re.IGNORECASE)]:
#
#     print "Assay CHEMBLID: %s" % bioactivity['assay_chemblid']
#
#     assay_data = json.loads(requests.get("https://www.ebi.ac.uk/chemblws/assays/%s.json" % bioactivity['assay_chemblid']).content)
#
#     print "  %s" % assay_data['assay']['assayDescription']
