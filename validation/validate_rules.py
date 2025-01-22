import argparse as ArgumentParser
import csv
import re
# required for reading in the OBO onotology files
import obonet

# Validation checks (automatic)
# - check that columns are present as per current version of spec
# –	are all gene rule IDs unique
# –	are all gene rule IDs using the same prefix 
# –	does organism always have a value
# –	does gene always have a value
# - if the gene value specifies a combo rule, are all those rule IDs present in the file?
# –	check that one of nodeID/refseq accession/GenBank accession/HMM accession are entered
# –	check ARO accession is either ‘-‘ or a value starting with ARO:, if starts with ARO:, checks it's a valid accession using latest CARD ontology
# –	check mutation always has a value
# –	check variation type always has a value, and it is one of the allowable variation types
# –	check context always has a value and is either core/acquired
# –	check that one of drug/drug class always has a value that isn’t ‘-‘
# –	check phenotype is either wildtype/nonwildtype
# –	check clinical category is S/I/R
# –	check there is a value in breakpoint, breakpoint standard, PMID and that these values aren’t ‘NA’ or ‘-‘
# –	check evidence code has a value, and it is one or more of the suggested values. If it's a new ECO code, print that to console as we may want to add it to the list of suggested values/check it's a valid code.
# –	if evidence limitations has a value, check it is one of the allowable values

##TODO:
# - check accessions actually exist (aro, refseq, genbank, hmm)
# - check variation type and mutation in combination, eg if its a gene presence detected, mutation should be '-', otherwise if it's dna it should be c. etc
# drug or drug class must be approved value
# list the unique values of breakpoint standards
# logic for checking the breakpoints... MIC < or <= for sensitive, MIC > for resistant, disk zone > for sensitive, disk zone < for resistant
# check that certain columns have only a single value, not multiple (eg mutation)
# is it possible to download allowable values directly from a google sheet in google drive?
# summarise at the end which checks have failures
# Kara has some example rules that I can test on

def parse_args():
    parser = ArgumentParser.ArgumentParser(description="Perform auto validation of AMRrules file")
    
    parser.add_argument("rules", help="Path to the rules file to check, must be a tab separated file")
    
    return parser.parse_args()

def get_column(column_name, draftrules):
        return [row[column_name] for row in draftrules if column_name in row]

def check_if_allowed_value(value_list, col_name, allowable_values):

    print("\nChecking  " + col_name + " column...")

    invalid_indices = [index for index, value in enumerate(value_list) if value == '' or value in ['NA', '-'] or value not in allowable_values]

    if not invalid_indices:
        print("✅ All " + col_name + " values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print(col_name + " column must contain one of the following values:\n" + ", ".join(allowable_values))
        for index in invalid_indices:
            print(f"Row {index + 1}: {value_list[index]}")

def check_if_not_missing(value_list, col_name):

    print("\nChecking  " + col_name + " column...")

    invalid_indices = [index for index, value in enumerate(value_list) if value == '' or value in ['NA', '-']]

    if not invalid_indices:
        print("✅ All " + col_name + " values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print(col_name + " column must contain a value that is not NA or '-'")
        for index in invalid_indices:
            print(f"Row {index + 1}: {value_list[index]}")

def check_ruleIDs(id_list):
    
    print("\nChecking ruleID column...")
    
    invalid_indices = []
    
    if len(id_list) != len(set(id_list)):
        print(f"❌ Rule IDs are not unique")
        invalid_indices = [index for index, value in enumerate(id_list) if id_list.count(value) > 1]
    else:
        print("All rule IDs have passed auto validation")

    prefix = id_list[0][:3]
    prefix_options = [prefix]
    for rule_id in id_list:
        if rule_id[:3] != prefix:
            invalid_indices.append(id_list.index(rule_id))
            # add the different prefix to the list
            prefix_options.append(rule_id[:3])
    
    if not invalid_indices and len(prefix_options) == 1:
        print("✅ All values are valid")
        print(f"Rule prefix: {prefix}")
    else:
        print(f"{len(invalid_indices)} rows have failed the check. Rule IDs must be unique and have the same prefix.")
        for index in invalid_indices:
            print(f"Row {index + 1}: {id_list[index]}")
        print(f"Multiple rule ID prefixes: " + ", ".join(prefix_options))

    # return this for later   
    return(prefix_options)

def check_organism(organism_list):
    
    print("\nChecking organism column...")
    
    # want to check if any values are missing, NA, or '-', and return which index number in the list where that's the case
    invalid_indices = [index for index, value in enumerate(organism_list) if value in ['NA', '-'] or value == '']
    
    for index, value in enumerate(organism_list):
        if index not in invalid_indices and not value.startswith("s__"):
            invalid_indices.append(index)
            #print(f"Organism must start with s__ to denote species. Invalid name at {index + 1}: {value}")
    # now check that all values start with s__, skip any invalid indices
    for index, value in enumerate(organism_list):
        if index not in invalid_indices:
            if not value.startswith("s__"):
                invalid_indices.append(index)
    
    if not invalid_indices:
        print("✅ All organism names passed auto validation")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Organism names must be present, not 'NA' or '-', and start with 's__'")
        for index in invalid_indices:
            print(f"Row {index + 1}: {organism_list[index]}")
    
    unique_organisms = set(organism_list)
    unique_organisms_str = ', '.join(map(str, unique_organisms))
    print(f"\nUnique organism names: {unique_organisms_str}")

def check_gene(gene_list, rule_list):
    
    print("\nChecking gene column...")
    
    # want to check if any values are missing, NA, or '-', and return which index number in the list where that's the case
    invalid_indices = [index for index, value in enumerate(gene_list) if value == '' or value in ['NA', '-']]

    if not invalid_indices:
        print("✅ All gene values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Gene names must be present, not 'NA' or '-'")
        for index in invalid_indices:
            print(f"Row {index + 1}: {gene_list[index]}")

    # now we want to check for gene names that are actually combo rules - if there are any, we want to check that any rule IDs mentioned here are present in the file already
    # if there is a value in gene list that follows the format of three capital letters followed by a string of four numbers, this is one to compare against rule ids
    print("\nNow checking for combinatorial rules in gene column...")

    invalid_indices_dict = {}

    pattern = re.compile(r'[A-Z]{3}\d{4}')
    for index, gene in enumerate(gene_list):
        if index not in invalid_indices:
            matches = pattern.findall(gene)
            for match in matches:
                if match not in rule_list:
                    invalid_indices_dict[index + 1] = match
                    break
    
    if not invalid_indices_dict:
        print("✅ All gene combinatorial rule IDs are valid")
    
    else:
        print(f"❌ {len(invalid_indices_dict)} rows have failed the check")
        for index, rule_id in invalid_indices_dict.items():
            print(f"Row {index}: ruleID {rule_id} is not present in the list of rules")

def check_id_accessions(nodeID_list, refseq_list, genbank_list, hmm_list, refseq_file, refseq_nodes, hmm_file):
    
    print("\nChecking nodeID, refseq accession, GenBank accession and HMM accession columns...")

    # parse the refseq file - get the refseq nucl and prot accessions, genbank accessions and hmm accessions to check against
    # Combine the nucl and prot accessions together
    refseq_accessions = [] 
    genbank_accessions = []
    refseq = csv.DictReader(open(refseq_file, newline=''), delimiter='\t')
    for row in refseq:
        if "RefSeq protein" in row:
            refseq_accessions.append(row["RefSeq protein"])
        if "RefSeq nucleotide" in row:
            refseq_accessions.append(row["RefSeq nucleotide"])
        if "GenBank protein" in row:
            genbank_accessions.append(row["GenBank protein"])
        if "GenBank nucleotide" in row:
            genbank_accessions.append(row["GenBank nucleotide"])
    # remove any empty strings
    refseq_accessions = [value for value in refseq_accessions if value != ""]
    genbank_accessions = [value for value in genbank_accessions if value != ""]

    refseq_node_ids = []
    refseq_hierarchy = csv.DictReader(open(refseq_nodes, newline=''), delimiter='\t')
    for row in refseq_hierarchy:
        if "parent_node_id" in row:
            refseq_node_ids.append(row["parent_node_id"])
    # remove any duplicates and empty strings
    refseq_node_ids = set(refseq_node_ids)
    refseq_node_ids = [value for value in refseq_node_ids if value != ""]

    hmm_accessions = []
    hmm_table = csv.DictReader(open(hmm_file, newline=''), delimiter='\t')
    for row in hmm_table:
        if "#Accession" in row:
            hmm_accessions.append(row["#Accession"])
    # remove any empty strings
    hmm_accessions  = [value for value in hmm_accessions if value != ""]

    # Check that in combination, at least one of these columns has a value
    invalid_indices = []
    for index, values in enumerate(zip(nodeID_list, refseq_list, genbank_list, hmm_list)):
        if all(value == '' or value in ['NA', '-'] for value in values):
            invalid_indices.append(index)
    
    # now we need to assess each list individually against the relevant accession lists
    # fine in the value is '-' as this is allowed in individual columns, just not all in combo
    # store the output so we can easily print what rows and columns are the issues
    invalid_node = []
    for index, value in enumerate(nodeID_list):
        if index not in invalid_indices and value not in refseq_node_ids and value != '-':
            invalid_node.append(index)
    invalid_refseq = []
    for index, value in enumerate(refseq_list):
        if index not in invalid_indices and value not in refseq_accessions and value != '-':
            invalid_refseq.append(index)
    invalid_gb = []
    for index, value in enumerate(genbank_list):
        if index not in invalid_indices and value not in genbank_accessions and value != '-':
            invalid_gb.append(index)
    invalid_hmm = []
    for index, value in enumerate(hmm_list):
        if index not in invalid_indices and value not in hmm_accessions and value != '-':
            invalid_hmm.append(index)
    
    if not invalid_indices:
        print("✅ All rows contain at least one value in one of these columns.")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check, and must contain at least one value in either nodeID, refseq accession, GenBank accession and HMM accession.")
        for index in invalid_indices:
            print(f"Row {index + 1}")
    if invalid_node or invalid_refseq or invalid_gb or invalid_hmm:
        print(f"❌ One or more accessions aren't present in either the NCBI Reference Gene Catalog (for nodeID, refseq accession and genbank accession) or the NCBI Reference HMM Catalog (for HMM accession).")
        print("\nInvalid nodeID accessions:")
        for index in invalid_node:
            print(f"Row {index + 1}: {nodeID_list[index]}")
        print("\nInvalid refseq accessions:")
        for index in invalid_refseq:
            print(f"Row {index + 1}: {refseq_list[index]}")
        print("\nInvalid genbank accessions:")
        for index in invalid_gb:
            print(f"Row {index + 1}: {genbank_list[index]}")
        print("\nInvalid HMM accessions:")
        for index in invalid_hmm:
            print(f"Row {index + 1}: {hmm_list[index]}")
    else:
        print("✅ All accessions are present in the relevant catalogues.")

def check_aro(aro_list, aro_obo_file):
    
    print("\nChecking ARO accession column...")

    # grab all the ARO accessions from the ARO ontology
    # Load the OBO file
    aro_obo = obonet.read_obo(aro_obo_file)

    # Extract all ARO terms
    aro_terms = [node for node in aro_obo.nodes if node.startswith('ARO:')]

    # value is invalid if doesn't start with ARO or isn't '-'
    # need to expand this as well to check if the ARO accessions in the list are actually ones that exist in the CARD ontology.
    invalid_indices = []
    for index, value in enumerate(aro_list):
        # a dash is fine
        if value == '-':
            continue
        # if the value is an accession, check it
        if value.startswith("ARO:"):
            if value not in aro_terms:
                invalid_indices.append(index)
        # the cell can't be empty, if it's empty it must be '-' (checked above)
        elif value == '' or value == 'NA':
            invalid_indices.append(index)
        # if it's anything else then it's invalid
        else:
            invalid_indices.append(index)

    if not invalid_indices:
        print("✅ All ARO accession values are valid and exist in the CARD ontology")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("ARO accession column must contain a valid ARO accession, and cannot be empty. The following rows contain invalid or empty accessions:")
        for index in invalid_indices:
            print(f"Row {index + 1}: {aro_list[index]}")

def check_mutation(mutation_list):
    
    print("\nChecking mutation column...")
    
    # check that there is either a value or '-' in this column
    invalid_indices = [index for index, value in enumerate(mutation_list) if value == '']

    if not invalid_indices:
        print("✅ All mutation values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Mutation column must contain either a value or '-' if no mutation required.")
        for index in invalid_indices:
            print(f"Row {index + 1}: {mutation_list[index]}")

def check_drug_drugclass(drug_list, drug_class_list):
   
    print("\nChecking drug and drug class columns...")
    
    # one of these columns must have a value in it
    invalid_indices = [index for index, values in enumerate(zip(drug_list, drug_class_list)) if all(value == '' or value in ['NA', '-'] for value in values)]
    
    if not invalid_indices:
        print("✅ All drug and drug class values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("One of drug or drug class must contain a value that is not NA or '-'")
        for index in invalid_indices:
            print(f"Row {index + 1}")

def check_evidence_code(evidence_code_list):
    
    print("\nChecking evidence code column...")

    allowable_values = ["ECO:0001091 knockout phenotypic evidence", "ECO:0000012 functional complementation evidence", "ECO:0001113 point mutation phenotypic evidence", "ECO:0000024 protein-binding evidence", "ECO:0001034 crystallography evidence", "ECO:0000005 enzymatic activity assay evidence", "ECO:0000042 gain-of-function mutant phenotypic evidence", "ECO:0007000 high throughput mutant phenotypic evidence", "ECO:0001103 natural variation mutant evidence", "ECO:0005027 genetic transformation evidence", "ECO:0000020 protein inhibition evidence"]

    # can be more than one of those values in this column, so need to split on the , separating them
    invalid_indices = []
    invalid_codes = []
    for index, value in enumerate(evidence_code_list):
        if value == '' or value in ['NA', '-']:
            invalid_indices.append(index)
            continue
        codes = [code.strip() for code in value.split(',')]
        for code in codes:
            if code not in allowable_values:
                invalid_indices.append(index)
                # only add the code to the list if it's got an eco code prefix
                if code.startswith("ECO:"):
                    invalid_codes.append(code)

    if not invalid_indices:
        print("✅ All evidence codes are valid")
    # otherwise we want to check if the ECO code still starts with "ECO:" but is a code that we haven't got in our suggested list
    # users can have other ECO codes but we want to flag these for potential inclusion in the allowable values
    else:
        if invalid_codes:
            unique_codes = set(invalid_codes)
            print("The following evidence codes are new and not currently in the list of suggested values:")
            print(f"{', '.join(unique_codes)}")
        print(f"\n❌ {len(invalid_indices)} rows have failed the check. Each rule must have an evidence code and not be empty. Please check all evidence codes are ECO codes.")
        for index in invalid_indices:
            print(f"Row {index + 1}: {evidence_code_list[index]}")

def check_evidence_grade_limitations(evidence_grade_list, evidence_limitations_list):

    print("\nChecking evidence grade and limitations columns...")

    allowable_grades = ["strong", "moderate", "weak"]
    allowable_limitations = ["lacks evidence for this species", "lacks evidence for this genus", "lacks evidence for this allele", "lacks evidence of the degree to which MIC is affected", "low clinical relevance", "unknown clinical relevance", "statistical geno/pheno evidence but no experimental evidence"]
    
    invalid_indices_grades = []
    invalid_indices_limitations = []

    for index, (grade, limitations) in enumerate(zip(evidence_grade_list, evidence_limitations_list)):
        if grade not in allowable_grades:
            invalid_indices_grades.append(index)
            #print(f"Invalid evidence grade at row {index + 1}: {grade}")
        elif grade in ["moderate", "weak"] and (limitations == '' or limitations == '-' or limitations not in allowable_limitations):
            invalid_indices_limitations.append(index)
            #print(f"Evidence grade {grade} at row {index + 1} requires valid limitations: {limitations}")

    if not invalid_indices_grades:
        print("✅ All evidence grades are valid")
    if not invalid_indices_limitations:
        print("✅ All evidence limitations are valid")
    
    if invalid_indices_grades or invalid_indices_limitations:
        print(f"❌ {len(invalid_indices_grades) + len(invalid_indices_limitations)} rows have failed the check")
        if invalid_indices_grades:
            print("Evidence grade column must contain one of the following values:\n" + ", ".join(allowable_grades))
            for index in invalid_indices_grades:
                print(f"Row {index + 1}: {evidence_grade_list[index]}")
        if invalid_indices_limitations:
            print("If evidence grade is 'moderate' or 'weak', evidence limitations column must contain one of the following values: " + ", ".join(allowable_limitations))
            for index in invalid_indices_limitations:
                print(f"Row {index + 1}: {evidence_limitations_list[index]}")

def main():
    args = parse_args()

    # read in the draft rules
    print(f"\nValidating rules file: {args.rules}")
    with open(args.rules, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        columns = reader.fieldnames
        draftrules = list(reader)
    
    expected_columns = ["ruleID", "organism", "gene", "nodeID", "refseq accession", "GenBank accession", "HMM accession", "ARO accession", "mutation", "variation type", "context", "drug", "drug class", "phenotype", "clinical category", "breakpoint", "breakpoint standard", "PMID", "evidence code", "evidence grade", "evidence limitations"]
    print("\nChecking that all required columns for spec v0.5 are present...")
    for column in expected_columns:
        if column not in columns:
            print(f"❌ {column} column not found in file.")
    print("\nContinuting to validate values in each column...")

    # grab the rule IDs and check
    if "ruleID" in columns:
        rule_ids = get_column("ruleID", draftrules)
        check_ruleIDs(rule_ids)
    else:
        print("\n❌ No ruleID column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")
        rule_ids = None

    # check that the value in organism is always present (can't be blank, '-' or NA, must start with s__)
    if "organism" in columns:
        check_organism(get_column("organism", draftrules))
    else:
        print("\n❌ No organism column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    # check that gene always has a value, and that this value is not '-' or NA or missing, as it should always be some kind of string
    if "gene" in columns and "ruleID" in columns:
        check_gene(get_column("gene", draftrules), rule_ids)
    if "gene" in columns and not "ruleID" in columns:
        print("\n❌ No ruleID column found in file. Spec v0.5 requires this column to be present, and cannot validate gene without it. Continuing to validate other columns...")
    else:
        print("\n❌ No gene column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    # check that for columns nodeID, refseq accession, GenBank accession, and HMM accession, at least one of these columns has a value
    # TODO: Automate downloads of relevant files from AMRFP ftp server?
    if "nodeID" in columns and "refseq accession" in columns and "GenBank accession" in columns and "HMM accession" in columns:
        refseq_file = 'refgenes_2024-12-18.1.tsv'
        hmm_file = 'hmms_amrfp_2024-12-18.1.tsv'
        refseq_nodes_file = 'ReferenceGeneHierarchy_2024-12-18.1.txt'
        check_id_accessions(get_column("nodeID", draftrules), get_column("refseq accession", draftrules), get_column("GenBank accession", draftrules), get_column("HMM accession", draftrules), refseq_file, refseq_nodes_file, hmm_file)
    else:
        for column in ["nodeID", "refseq accession", "GenBank accession", "HMM accession"]:
            if column not in columns:
                print(f"\n❌ {column} column not found in file.")
        print("\n❌ Spec v0.5 requires all of nodeID, refseq accession, GenBank accession, and HMM accession columns to be present in order to validate. Continuing to validate other columns...")

    if "ARO accession" in columns:
        # TODO: Replace this with URL to obo file
        aro_obo_file = 'card-ontology/aro.obo'
        check_aro(get_column("ARO accession", draftrules), aro_obo_file)
    else:
        print("\n❌ No ARO accession column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    if "mutation" in columns:
        check_mutation(get_column("mutation", draftrules))
    else:
        print("\n❌ No mutation column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")
    
    if "variation type" in columns:
        variation_allowed_types = ["Gene presence detected", "Protein variant detected", "Nucleotide variant detected", "Promoter variant detected", "Inactivating mutation detected", "Gene copy number variant detected", "Nucleotide variant detected in multi-copy gene", "Low frequency variant detected", "Combination"]
        check_if_allowed_value(get_column("variation type", draftrules), "variation type", variation_allowed_types)
    else:
        print("\n❌ No variation type column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    if "context" in columns:
        check_if_allowed_value(get_column("context", draftrules), "context", ["core", "acquired"])
    else:
        print("\n❌ No context column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    # NOTE: we want child terms under ARO:1000003 (antibiotic molecule) to get our drugs and drug classes
    # We will also need child terms under ARO:3000042 (beta-lactamase inhibitor) to get our inhibitor classes and names
    if "drug" in columns and "drug class" in columns:
        check_drug_drugclass(get_column("drug", draftrules), get_column("drug class", draftrules))
    else:
        for column in ["drug", "drug class"]:
            if column not in columns:
                print(f"\n❌ {column} column not found in file.")
        print("\n❌ Spec v0.5 requires at least both drug and drug class columns to be present in order to validate. Continuing to validate other columns...")

    if "phenotype" in columns:
        check_if_allowed_value(get_column("phenotype", draftrules), "phenotype", ["wildtype", "nonwildtype"])
    else:
        print("\n❌ No phenotype column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    if "clinical category" in columns:
        check_if_allowed_value(get_column("clinical category", draftrules), "clinical category", ["S", "I", "R"])
    else:
        print("\n❌ No clinical category column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    if "breakpoint" in columns:
        check_if_not_missing(get_column("breakpoint", draftrules), "breakpoint")
    else:
        print("\n❌ No breakpoint column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    if "breakpoint standard" in columns:
        check_if_not_missing(get_column("breakpoint standard", draftrules), "breakpoint standard")
    else:
        print("\n❌ No breakpoint standard column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    if "PMID" in columns:
        check_if_not_missing(get_column("PMID", draftrules), "PMID")
    else:
        print("\n❌ No PMID column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    if "evidence code" in columns:
        check_evidence_code(get_column("evidence code", draftrules))
    else:
        print("\n❌ No evidence code column found in file. Spec v0.5 requires this column to be present. Continuing to validate other columns...")

    if "evidence grade" in columns and "evidence limitations" in columns:
        check_evidence_grade_limitations(get_column("evidence grade", draftrules), get_column("evidence limitations", draftrules))
    else:
        for column in ["evidence grade", "evidence limitations"]:
            if column not in columns:
                print(f"\n❌ {column} column not found in file.")
        print("\n❌ Both evidence grade and limitations columns required for spec v0.5.")

if __name__ == "__main__":
    main()
