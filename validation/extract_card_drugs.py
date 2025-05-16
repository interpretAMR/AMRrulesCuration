import csv
import obonet

def parse_obo_file(file_path, term_id):
    children = []
    with open(file_path, 'r') as file:
        current_term = None
        for line in file:
            line = line.strip()
            if line.startswith("[Term]"):
                current_term = None
            elif line.startswith("id: "):
                current_term = line.split("id: ")[1]
            elif line.startswith("is_a: ") and current_term:
                parent = line.split("is_a: ")[1].split(" ! ")[0]
                if parent == term_id:
                    children.append(current_term)
    return children

card_ontology = obonet.read_obo('validation/card-ontology/aro.obo')
card_categories = csv.DictReader(open('validation/card-data/aro_categories.tsv'), delimiter='\t')

# get the drug names from the categories
drug_names = {} # key is the name, value is the aro accession
for row in card_categories:
    if row['ARO Category'] == 'Drug Class':
        drug_names[row['ARO Name']] = row['ARO Accession']

drug_names

id_to_name = {id_: data.get("name") for id_, data in card_ontology.nodes(data=True)}

output_rows = []

for drug_class, aro_accession in drug_names.items():
    try:
        children = parse_obo_file('validation/card-ontology/aro.obo', aro_accession)
        for child in children:
            child_name = id_to_name[child]
            output_rows.append([child, child_name, drug_class])
    except:
        print(f"Error processing {aro_accession} for {drug_class}")
        continue

# we then want to specifically pull out all the drug names for the different beta-lactam classes
# and also the polymyxin classes 
# additionally extracting 'antibiotic mixture', ARO:3000707 - this has betalactam + inhibitors
betalac_aros = ['ARO:3009105', 'ARO:3009106', 'ARO:3009107', 'ARO:3009108', 
             'ARO:3009109', 'ARO:3009123', 'ARO:3009124', 'ARO:3009125',
             'ARO:3000035', 'ARO:3007783', 'ARO:0000022', 'ARO:3007629',
             'ARO:3000707']

for aro_accession in betalac_aros:
    children = parse_obo_file('validation/card-ontology/aro.obo', aro_accession)
    drug_class = id_to_name[aro_accession]
    for child in children:
        child_name = id_to_name[child]
        output_rows.append([child, child_name, drug_class])


output_file = 'validation/card_drug_names.tsv'
with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(['ARO Accession', 'Drug Name', 'Drug Class'])
    for row in output_rows:
        writer.writerow(row)

