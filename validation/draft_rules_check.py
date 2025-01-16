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

# OUTPUT OF THE CHECK
# output rows where a check has failed, and state the value that appears to be incorrect
# make some text that helps with the manual check:
# –	print the rule prefix (to check is correct)
# –	print all unique organism names

def parse_args():
    parser = ArgumentParser.ArgumentParser(description="Check the draft rules file for any errors")
    
    parser.add_argument("rules", help="Path to the draft rules file to check, must be a tab separated file")
    
    return parser.parse_args()

def check_ruleIDs(id_list):
    if len(id_list) != len(set(id_list)):
        print("Rule IDs are not unique")
    else:
        print("All rule IDs are unique")

    prefix = id_list[0][:3]
    prefix_options = [prefix]
    for rule_id in id_list:
        if rule_id[:3] != prefix:
            # add the different prefix to the list
            prefix_options.append(rule_id[:3])
    if len(prefix_options) == 1:
        print("All rule IDs have the same prefix: " + prefix)
    else:
        print("Multiple rule ID prefixes: " + ", ".join(prefix_options))

    # return this for later   
    return(prefix_options)

def check_organism(organism_list):
    # want to check if any values are missing, NA, or '-', and return which index number in the list where that's the case
    invalid_indices = [index for index, value in enumerate(organism_list) if pd.isna(value) or value in ['NA', '-']]
    if invalid_indices:
        print(f"Invalid organism names found at rows: {[index + 1 for index in invalid_indices]}")
    # now check that all values start with s__, skip any invalid indices
    for index, value in enumerate(organism_list):
        if index not in invalid_indices:
            if not value.startswith("s__"):
                print(f"Organism must start with s__ to denote species. Invalid name at {index + 1}: {value}")
                invalid_indices.append(index)
    
    if not invalid_indices:
        print("All organism names are valid.")

def check_gene(gene_list, rule_list):
    # want to check if any values are missing, NA, or '-', and return which index number in the list where that's the case
    invalid_indices = [index for index, value in enumerate(gene_list) if pd.isna(value) or value in ['NA', '-']]
    if invalid_indices:
        print(f"Invalid gene names found at rows: {[index + 1 for index in invalid_indices]}")
    # now we want to check for gene names that are actually combo rules - if there are any, we want to check that any rule IDs mentioned here are present in the file already
    # if there is a value in gene list that follows the format of three capital letters followed by a string of four numbers, this is one to compare against rule ids
    pattern = re.compile(r'[A-Z]{3}\d{4}')
    for index, gene in enumerate(gene_list):
        if index not in invalid_indices:
            matches = pattern.findall(gene)
            for match in matches:
                if match not in rule_list:
                    print(f"Gene value {match} at row {index + 1} is not a valid rule ID")
                    invalid_indices.append(index)
                    break
    if not invalid_indices:
        print("All gene names are valid")

def check_id_accessions(nodeID_list, refseq_list, genbank_list, hmm_list):
    invalid_indices = []
    for index, values in enumerate(zip(nodeID_list, refseq_list, genbank_list, hmm_list)):
        if all(pd.isna(value) or value in ['NA', '-'] for value in values):
            invalid_indices.append(index)
    if invalid_indices:
        print(f"Invalid ID/accession values found at rows: {[index + 1 for index in invalid_indices]}")
    else:
        print("All ID/accession values are valid")

def check_aro(aro_list):
    # value is invalid if doesn't start with ARO or isn't '-'
    invalid_indices = [index for index, value in enumerate(aro_list) if pd.isna(value) or value != '-' and not value.startswith("ARO")]
    if invalid_indices:
        print(f"Invalid ARO accession values found at rows: {[index + 1 for index in invalid_indices]}")
        print("ARO accession column must contain either '-' (if no ARO accession) or start with ARO")
    else:
        print("All ARO accession values are valid")

def check_mutation(mutation_list):
    # check that there is either a value or '-' in this column
    invalid_indices = [index for index, value in enumerate(mutation_list) if pd.isna(value)]
    if invalid_indices:
        print(f"Invalid mutation values found at rows: {[index + 1 for index in invalid_indices]}")
        print("Mutation column must contain either a value or '-'")
    else:
        print("All mutation values are valid")

def check_variation_type(variation_list):
    # the allowable values for this column are:
    allowable_values = ["Gene presence detected", "Protein variant detected", "Nucleotide variant detected", "Promoter variant detected", "Inactivating mutation detected", "Gene copy number variant detected", "Nucleotide variant detected in multi-copy gene", "Low frequency variant detected", "Combination"]
    # check that the value in this column is one of these
    invalid_indices = [index for index, value in enumerate(variation_list) if value not in allowable_values]
    if invalid_indices:
        print(f"Invalid variation type values found at rows: {[index + 1 for index in invalid_indices]}")
        print("Variation type column must contain one of the following values: " + ", ".join(allowable_values))
    else:
        print("All variation type values are valid")

def check_context(context_list):
    allowable_values = ["core", "acquired"]
    # check that the value in this column is one of these
    invalid_indices = [index for index, value in enumerate(context_list) if value not in allowable_values]
    if invalid_indices:
        print(f"Invalid context values found at rows: {[index + 1 for index in invalid_indices]}")
        print("Context column must contain either 'core' or 'acquired'")
    else:
        print("All context values are valid")

def check_drug_drugclass(drug_list, drug_class_list):
    # one of these columns must have a value in it
    invalid_indices = [index for index, values in enumerate(zip(drug_list, drug_class_list)) if all(pd.isna(value) or value in ['NA', '-'] for value in values)]
    if invalid_indices:
        print(f"Invalid drug/drug class values found at rows: {[index + 1 for index in invalid_indices]}")
    else:
        print("All drug and drug class values are valid")

def check_phenotype(phenotype_list):
    allowable_values = ["wildtype", "nonwildtype"]
    # check that the value in this column is one of these
    invalid_indices = [index for index, value in enumerate(phenotype_list) if value not in allowable_values]
    if invalid_indices:
        print(f"Invalid phenotype values found at rows: {[index + 1 for index in invalid_indices]}")
        print("Phenotype column must contain either 'wildtype' or 'nonwildtype'")
    else:
        print("All phenotype values are valid")

def check_sir(sir_list):
    allowable_values = ["S", "I", "R"]
    # check that the value in this column is one of these
    invalid_indices = [index for index, value in enumerate(sir_list) if value not in allowable_values]
    if invalid_indices:
        print(f"Invalid clinical category values found at rows: {[index + 1 for index in invalid_indices]}")
        print("Clincal category column must contain either 'S', 'I' or 'R'")
    else:
        print("All clinical category values are valid")

def check_breakpoint(breakpoint_list):
    # check that there is a value in this column that isn't NA or '-'
    invalid_indices = [index for index, value in enumerate(breakpoint_list) if pd.isna(value) or value in ['NA', '-']]
    if invalid_indices:
        print(f"Invalid breakpoint values found at rows: {[index + 1 for index in invalid_indices]}")
        print("Breakpoint column must contain a value that is not NA or '-'")
    else:
        print("All breakpoint values are valid")

def check_breakpoint_standard(breakpoint_std_list):
    invalid_indices = [index for index, value in enumerate(breakpoint_std_list) if pd.isna(value) or value in ['NA', '-']]
    if invalid_indices:
        print(f"Invalid breakpoint standard values found at rows: {[index + 1 for index in invalid_indices]}")
        print("Breakpoint standard column must contain a value that is not NA or '-'")
    else:
        print("All breakpoint standard values are valid")

def check_pmid(pmid_list):
    invalid_indices = [index for index, value in enumerate(pmid_list) if pd.isna(value) or value in ['NA', '-']]
    if invalid_indices:
        print(f"Invalid PMID values found at rows: {[index + 1 for index in invalid_indices]}")
        print("PMID column must contain a value that is not NA or '-'")
    else:
        print("All PMID values are valid")

def check_evidence_code(evidence_code_list):
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

    if invalid_indices:
        print(f"Invalid evidence codes found at rows: {[index + 1 for index in invalid_indices]}")
    else:
        print("All evidence codes are valid")

def check_evidence_grade_limitations(evidence_grade_list, evidence_limitations_list):
    allowable_grades = ["strong", "moderate", "weak"]
    allowable_limitations = ["lacks evidence for this species", "lacks evidence for this genus", "lacks evidence for this allele", "lacks evidence of the degree to which MIC is affected", "low clinical relevance", "unknown clinical relevance", "statistical geno/pheno evidence but no experimental evidence"]
    
    invalid_indices = []
    for index, (grade, limitations) in enumerate(zip(evidence_grade_list, evidence_limitations_list)):
        if grade not in allowable_grades:
            invalid_indices.append(index)
            print(f"Invalid evidence grade at row {index + 1}: {grade}")
        elif grade in ["moderate", "weak"] and (pd.isna(limitations) or limitations == '' or limitations == '-' or limitations not in allowable_limitations):
            invalid_indices.append(index)
            print(f"Evidence grade {grade} at row {index + 1} requires valid limitations: {limitations}")

    if not invalid_indices:
        print("All evidence grades and limitations are valid")

def main():
    args = parse_args()

    # read in the draft rules
    draftrules = pd.read_csv(args.rules, sep="\t")

    # grab the rule IDs and check
    rule_ids = draftrules["ruleID"].tolist()
    rule_prefix = check_ruleIDs(rule_ids)

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

    """
    How do we want this output to look...

    Okay we want a subheading, that lists what column/s we're checking

    and then we want a new line - if every row has passed the check, then we want the tick symbol, and say "All values are valid"

    For the rule ID and organism columns, also list what the rule ID prefix is (so we can check it's an allowed one),
    and for the organims list all unique organisms (so we can make sure they're all accepted species or whatever)

    if not every row has passed the check, we want a summary line that says "X rows have failed the check"
    Second summary line should state what we're checking for (either give the list of allowed values, or the fact that it has be be '-', or whatever the combination rule for the set of columns are)
    Then list, line by line, each row index that has failed, and why it has failed

    """

    # print the rule prefixes and all unique organism names
    print(f"Rule prefixes: {', '.join(rule_prefix)}")
    unique_organisms = draftrules['organism'].dropna().unique()
    print(f"Unique organism names: {', '.join(unique_organisms)}")

if __name__ == "__main__":
    main()
