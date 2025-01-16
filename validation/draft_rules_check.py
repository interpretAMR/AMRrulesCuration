import argparse as ArgumentParser
import pandas as pd
import re

# Validation checks (automatic)
# –	are all gene rule IDs unique
# –	are all gene rule IDs using the same prefix 
# –	does organism always have a value
# –	does gene always have a value
# - if the gene value specifies a combo rule, are all those rule IDs present in the file?
# –	check that one of nodeID/refseq accession/GenBank accession/HMM accession are entered
# –	check ARO accession is either ‘-‘ or a value starting with ARO:
# –	check mutation always has a value
# –	check variation type always has a value, and it is one of the allowable variation types
# –	check context always has a value and is either core/acquired
# –	check that one of drug/drug class always has a value that isn’t ‘-‘
# –	check phenotype is either wildtype/nonwildtype
# –	check clinical category is S/I/R
# –	check there is a value in breakpoint, breakpoint standard, PMID and that these values aren’t ‘NA’ or ‘-‘
# –	check evidence code has a value, and it is one or more of the allowable values
# –	if evidence limitations has a value, check it is one of the allowable values

def parse_args():
    parser = ArgumentParser.ArgumentParser(description="Check the draft rules file for any errors")
    
    parser.add_argument("rules", help="Path to the draft rules file to check, must be a tab separated file")
    
    return parser.parse_args()

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
    invalid_indices = [index for index, value in enumerate(organism_list) if pd.isna(value) or value in ['NA', '-']]
    
    for index, value in enumerate(organism_list):
        if index not in invalid_indices and not value.startswith("s__"):
            invalid_indices.append(index)
            print(f"Organism must start with s__ to denote species. Invalid name at {index + 1}: {value}")
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
    
    unique_organisms = pd.Series(organism_list).unique()
    unique_organisms_str = ', '.join(map(str, unique_organisms))
    print(f"\nUnique organism names: {unique_organisms_str}")

def check_gene(gene_list, rule_list):
    
    print("\nChecking gene column...")
    
    # want to check if any values are missing, NA, or '-', and return which index number in the list where that's the case
    invalid_indices = [index for index, value in enumerate(gene_list) if pd.isna(value) or value in ['NA', '-']]

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

def check_id_accessions(nodeID_list, refseq_list, genbank_list, hmm_list):
    
    print("\nChecking nodeID, refseq accession, GenBank accession and HMM accession columns...")

    invalid_indices = []
    for index, values in enumerate(zip(nodeID_list, refseq_list, genbank_list, hmm_list)):
        if all(pd.isna(value) or value in ['NA', '-'] for value in values):
            invalid_indices.append(index)
    
    if not invalid_indices:
        print("✅ All rows contain at least one value in one of these columns")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check,and must contain at least one value in either nodeID, refseq accession, GenBank accession and HMM accession")
        for index in invalid_indices:
            print(f"Row {index + 1}")

def check_aro(aro_list):
    
    print("\nChecking ARO accession column...")

    # value is invalid if doesn't start with ARO or isn't '-'
    invalid_indices = [index for index, value in enumerate(aro_list) if pd.isna(value) or value != '-' and not value.startswith("ARO")]

    if not invalid_indices:
        print("✅ All ARO accession values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("ARO accession column must contain either '-' (if no ARO accession) or start with 'ARO:'")
        for index in invalid_indices:
            print(f"Row {index + 1}: {aro_list[index]}")

def check_mutation(mutation_list):
    
    print("\nChecking mutation column...")
    
    # check that there is either a value or '-' in this column
    invalid_indices = [index for index, value in enumerate(mutation_list) if pd.isna(value)]

    if not invalid_indices:
        print("✅ All mutation values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Mutation column must contain either a value or '-' if no mutation required.")
        for index in invalid_indices:
            print(f"Row {index + 1}: {mutation_list[index]}")

def check_variation_type(variation_list):
    
    print("\nChecking variation type column...")

    # the allowable values for this column are:
    allowable_values = ["Gene presence detected", "Protein variant detected", "Nucleotide variant detected", "Promoter variant detected", "Inactivating mutation detected", "Gene copy number variant detected", "Nucleotide variant detected in multi-copy gene", "Low frequency variant detected", "Combination"]

    invalid_indices = []

    for index, value in enumerate(variation_list):
        if pd.isna(value) or value not in allowable_values:
            invalid_indices.append(index)

    if not invalid_indices:
        print("✅ All variation type values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Variation type column must contain one of the following values: " + ", ".join(allowable_values))
        for index in invalid_indices:
            print(f"Row {index + 1}: {variation_list[index]}")

def check_context(context_list):
    
    print("\nChecking context column...")

    allowable_values = ["core", "acquired"]
    # check that the value in this column is one of these
    invalid_indices = [index for index, value in enumerate(context_list) if value not in allowable_values]

    if not invalid_indices:
        print("✅ All context values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Context column must contain either 'core' or 'acquired'")
        for index in invalid_indices:
            print(f"Row {index + 1}: {context_list[index]}")

def check_drug_drugclass(drug_list, drug_class_list):
   
    print("\nChecking drug and drug class columns...")
    
    # one of these columns must have a value in it
    invalid_indices = [index for index, values in enumerate(zip(drug_list, drug_class_list)) if all(pd.isna(value) or value in ['NA', '-'] for value in values)]
    
    if not invalid_indices:
        print("✅ All drug and drug class values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("One of drug or drug class must contain a value that is not NA or '-'")
        for index in invalid_indices:
            print(f"Row {index + 1}")

def check_phenotype(phenotype_list):
    
    print("\nChecking phenotype column...")

    allowable_values = ["wildtype", "nonwildtype"]
    # check that the value in this column is one of these
    invalid_indices = [index for index, value in enumerate(phenotype_list) if value not in allowable_values]

    if not invalid_indices:
        print("✅ All phenotype values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Phenotype column must contain either 'wildtype' or 'nonwildtype'")
        for index in invalid_indices:
            print(f"Row {index + 1}: {phenotype_list[index]}")

def check_sir(sir_list):
    
    print("\nChecking clinical category column...")

    allowable_values = ["S", "I", "R"]
    # check that the value in this column is one of these
    invalid_indices = [index for index, value in enumerate(sir_list) if value not in allowable_values]

    if not invalid_indices:
        print("✅ All clinical category values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Clincal category column must contain either 'S', 'I' or 'R'")
        for index in invalid_indices:
            print(f"Row {index + 1}: {sir_list[index]}")

def check_breakpoint(breakpoint_list):
    
    print("\nChecking breakpoint column...")

    # check that there is a value in this column that isn't NA or '-'
    invalid_indices = [index for index, value in enumerate(breakpoint_list) if pd.isna(value) or value in ['NA', '-']]

    if not invalid_indices:
        print("✅ All breakpoint values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Breakpoint column must contain a value that is not NA or '-'")
        for index in invalid_indices:
            print(f"Row {index + 1}: {breakpoint_list[index]}")

def check_breakpoint_standard(breakpoint_std_list):

    print("\nChecking breakpoint standard column...")

    invalid_indices = [index for index, value in enumerate(breakpoint_std_list) if pd.isna(value) or value in ['NA', '-']]

    if not invalid_indices:
        print("✅ All breakpoint standard values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Breakpoint standard column must contain a value that is not NA or '-'")
        for index in invalid_indices:
            print(f"Row {index + 1}: {breakpoint_std_list[index]}")

def check_pmid(pmid_list):

    print("\nChecking PMID column...")

    invalid_indices = [index for index, value in enumerate(pmid_list) if pd.isna(value) or value in ['NA', '-']]

    if not invalid_indices:
        print("✅ All PMID values are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("PMID column must contain a value that is not NA or '-'")
        for index in invalid_indices:
            print(f"Row {index + 1}: {pmid_list[index]}")

def check_evidence_code(evidence_code_list):
    
    print("\nChecking evidence code column...")

    allowable_values = ["ECO:0001091 knockout phenotypic evidence", "ECO:0000012 functional complementation evidence", "ECO:0001113 point mutation phenotypic evidence", "ECO:0000024 protein-binding evidence", "ECO:0001034 crystallography evidence", "ECO:0000005 enzymatic activity assay evidence", "ECO:0000042 gain-of-function mutant phenotypic evidence", "ECO:0007000 high throughput mutant phenotypic evidence", "ECO:0001103 natural variation mutant evidence", "ECO:0005027 genetic transformation evidence", "ECO:0000020 protein inhibition evidence"]

    # can be more than one of those values in this column, so need to split on the , separating them
    invalid_indices = []
    for index, value in enumerate(evidence_code_list):
        if pd.isna(value) or value in ['NA', '-']:
            invalid_indices.append(index)
            break
        codes = [code.strip() for code in value.split(',')]
        for code in codes:
            if code not in allowable_values:
                invalid_indices.append(index)
                break
    
    if not invalid_indices:
        print("✅ All evidence codes are valid")
    else:
        print(f"❌ {len(invalid_indices)} rows have failed the check")
        print("Evidence code column must contain one or more of the following values:\n" + ", ".join(allowable_values))
        print("Please ensure all evidence codes match one of the allowable codes")
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
        elif grade in ["moderate", "weak"] and (pd.isna(limitations) or limitations == '' or limitations == '-' or limitations not in allowable_limitations):
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
    print(f"\nValidating draft rules file: {args.rules}")
    draftrules = pd.read_csv(args.rules, sep="\t")

    # grab the rule IDs and check
    rule_ids = draftrules["ruleID"].tolist()
    check_ruleIDs(rule_ids)

    # check that the value in organism is always present (can't be blank, '-' or NA, must start with s__)
    check_organism(draftrules["organism"].tolist())

    # check that gene always has a value, and that this value is not '-' or NA or missing, as it should always be some kind of string
    check_gene(draftrules["gene"].tolist(), rule_ids)

    # check that for columns nodeID, refseq accession, GenBank accession, and HMM accession, at least one of these columns has a value
    check_id_accessions(draftrules["nodeID"].tolist(), draftrules["refseq accession"].tolist(), draftrules["GenBank accession"].tolist(), draftrules["HMM accession"].tolist())

    check_aro(draftrules["ARO accession"].tolist())

    check_mutation(draftrules["mutation"].tolist())

    check_variation_type(draftrules["variation type"].tolist())

    check_context(draftrules["context"].tolist())

    check_drug_drugclass(draftrules["drug"].tolist(), draftrules["drug class"].tolist())

    check_phenotype(draftrules["phenotype"].tolist())

    check_sir(draftrules["clinical category"].tolist())

    check_breakpoint(draftrules["breakpoint"].tolist())

    check_breakpoint_standard(draftrules["breakpoint standard"].tolist())

    check_pmid(draftrules["PMID"].tolist())

    check_evidence_code(draftrules["evidence code"].tolist())

    check_evidence_grade_limitations(draftrules["evidence grade"].tolist(), draftrules["evidence limitations"].tolist())

if __name__ == "__main__":
    main()
