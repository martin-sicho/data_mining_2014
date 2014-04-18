import requests
import json
import re

######

def isFloat(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

########################################################################

# 1. Use UniProt accession to get target details

print """
# =========================================================
# 1. Use UniProt accession to get target details
# =========================================================
"""

accession = 'Q00534'

target_data = json.loads(requests.get("https://www.ebi.ac.uk/chemblws/targets/uniprot/%s.json" % accession).content)

print "Target Description: %s" % target_data['target']['description']
print "Target CHEMBLID:    %s" % target_data['target']['chemblId']

# 2. Get all bioactivties for target CHEMBL_ID

print """

# =========================================================
# 2. Get all bioactivties for target CHEMBL_ID
# =========================================================
"""

bioactivity_data = json.loads(requests.get("https://www.ebi.ac.uk/chemblws/targets/%s/bioactivities.json" % target_data['target']['chemblId']).content)

print "Bioactivity Count:           %d" % len(bioactivity_data['bioactivities'])
print "Bioactivity Count (IC50's):  %d" % len([record for record in bioactivity_data['bioactivities'] if record['bioactivity_type'] == 'IC50'])

# 3. Get compounds with high binding affinity (IC50 < 100)

print """

# =========================================================
# 3. Get compounds with high binding affinity (IC50 < 100)
# =========================================================
"""

for bioactivity in [record for record in bioactivity_data['bioactivities'] if re.search('IC50', record['bioactivity_type']) and isFloat(record['value']) and float(record['value']) < 100]:

    print "Compound CHEMBLID: %s" % bioactivity['ingredient_cmpd_chemblid']

    cmpd_data = json.loads(requests.get("https://www.ebi.ac.uk/chemblws/compounds/%s.json" % bioactivity['ingredient_cmpd_chemblid']).content)

    print "  %s" % cmpd_data['compound']['smiles']

# 4. Get assay details foe Ki actvity types

print """

# =========================================================
# 4. Get assay details foe Ki actvity types
# =========================================================
"""

for bioactivity in [record for record in bioactivity_data['bioactivities'] if re.search('Ki', record['bioactivity_type'], re.IGNORECASE)]:

    print "Assay CHEMBLID: %s" % bioactivity['assay_chemblid']

    assay_data = json.loads(requests.get("https://www.ebi.ac.uk/chemblws/assays/%s.json" % bioactivity['assay_chemblid']).content)

    print "  %s" % assay_data['assay']['assayDescription']
