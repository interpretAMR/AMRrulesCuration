import argparse
import csv
import urllib.request
import io, re, os

#TODO: probably put this sort of thing in a utils file that we import
aa_conversion = {'G': 'Gly', 'A': 'Ala', 'S': 'Ser', 'P': 'Pro', 'T': 'Thr', 'C': 'Cys', 'V': 'Val', 'L': 'Leu', 'I': 'Ile', 
                 'M': 'Met', 'N': 'Asn', 'Q': 'Gln', 'K': 'Lys', 'R': 'Arg', 'H': 'His', 'D': 'Asp', 'E': 'Glu', 'W': 'Trp', 
                 'Y': 'Tyr', 'F': 'Phe', '*': 'STOP'}

def parse_args():
    parser = argparse.ArgumentParser(description="Interpreter for AMRrules.")

    parser.add_argument('--input', type=str, required=True, help='Path to the input file.')
    parser.add_argument('--output_prefix', type=str, required=True, help='Prefix name for the output files.')
    parser.add_argument('--output_dir', type=str, default=os.getcwd(), help='Directory to write the output files to. Default is current directory.')
    parser.add_argument('--rules', '-r', type=str, required=True, help='Path to the rules file, in tab-delimited format.')
    #TODO: implement card and resfinder options, currently only amrfp is supported
    parser.add_argument('--amr_tool', '-t', type=str, default='amrfp', help='AMR tool used to detect genotypes: options are amrfp, card, resfinder. Currently only amrfp is supported.')
    parser.add_argument('--hamronized', '-H', action='store_true', help='Input file has been hamronized')
    parser.add_argument('--amrfp_db_version', type=str, default='latest', help='Version of the AMRFP database used. Default is latest. NOTE STILL TO BE IMPLEMENTED')
    parser.add_argument('--annot_opts', '-a', type=str, default='minimal', help='Annotation options: minimal (context, drug, phenotype, category, evidence grade), full (everything including breakpoints, standards, etc)')
    
    return parser.parse_args()

def download_and_parse_reference_gene_hierarchy(url):
    print(f"Downloading Reference Gene Hierarchy from {url}...")
    
    # Use urllib to download the file
    with urllib.request.urlopen(url) as response:
        content = response.read().decode('utf-8')  # Decode the response as UTF-8
    
    # Parse the downloaded content as a TSV file
    content_io = io.StringIO(content)
    reader = csv.DictReader(content_io, delimiter='\t')
    
    amrfp_nodes = {}
    for row in reader:
        node_id = row.get('node_id')
        parent_node = row.get('parent_node_id')
        amrfp_nodes[node_id] = parent_node

    print(f"Downloaded and parsed {len(amrfp_nodes)} rows from the Reference Gene Hierarchy.")
    return amrfp_nodes

def parse_rules_file(rules_file, tool):
    rules = []
    #if tool == 'amrfp':
    #    key_value = 'nodeID'
    with open(rules_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rules.append(row)
    return rules

def extract_mutation(row):
    # deal with either header option from amrfp uggghhhh
    gene_or_element_symbol = row.get('Gene symbol') or row.get('Element symbol')
    # if we've got point mutations, we need to extract the actual mutation
    # and convert it to the AMRrules syntax so we can identify the correct rule
    gene_symbol, mutation = gene_or_element_symbol.rsplit("_", 1)

    _, ref, pos, alt, _ = re.split(r"(\D+)(\d+)(\D+)", mutation)
    # this means it is a protein mutation
    if row.get('Method') in ["POINTX", "POINTP"]:
        # convert the single letter AA code to the 3 letter code
        ref = aa_conversion.get(ref)
        alt = aa_conversion.get(alt)
        return(f"p.{ref}{pos}{alt}", "Protein variant detected")
    elif row.get('Method') == "POINTN":
        # e.g., 23S_G2032T -> c.2032G>T
        return(f"c.{pos}{ref}>{alt}", "Nucleotide variant detected")

def check_rules(row, rules, amrfp_nodes):

    # if the row is a point mutation, we need to extract that info and look only for those rules
    if row.get('Element subtype') == "POINT" or row.get('Subtype') == "POINT":
        amrrules_mutation, type = extract_mutation(row)
    elif row.get('Element subtype') == "AMR" or row.get('Subtype') == "AMR":
        amrrules_mutation = None
        type = 'Gene presence detected'

    # select rules that match our variation type
    rules_to_check = []
    for rule in rules:
            if rule['variation type'] == type:
                rules_to_check.append(rule)

    # we need to set some logic here about what accessions we're going to be looking for
    # if there's a nodeID, then that's where we want to start
    # otherwise we're going to be checking refseq, or HMM accessions
    #TODO Genbank accessions - not relevant for amrfp as these won't be reported in the output file
    hierarchy_id = row.get('Hierarchy node')
    seq_acc = row.get('Accession of closest sequence')

    # only grab HMM if it's actually listed
    if row.get('HMM id') != 'NA':
        HMM_acc = row.get('HMM id')
    else:
        HMM_acc = None

    # First we're going to check for the hierarchy ID, and if we have one or matches, we we return that
    matching_rules = [rule for rule in rules_to_check if rule.get('nodeID') == hierarchy_id]
    if len(matching_rules) > 0:
        if type == 'Gene presence detected':
            # just return these rules
            return matching_rules
        elif type == 'Protein variant detected' or type == 'Nucleotide variant detected':
            # only return rules with the appropriate mutation
            final_matching_rules = []
            # now we need to check the mutation, extracting any matching rules
            for rule in matching_rules:
                if rule['mutation'] == amrrules_mutation:
                    final_matching_rules.append(rule)
            return final_matching_rules
    
    # Okay so nothing matched directly to the nodeID, or we would've returned out of the function. 
    # So now we need to check if there's a parent node that matches
    # Note we don't need to check for variation type here, because we're not going to need to go up the hierarchy
    # for that type of rule
    parent_node = amrfp_nodes.get(hierarchy_id)
    while parent_node is not None and parent_node != 'AMR':
        matching_rules = [rule for rule in rules_to_check if rule.get('nodeID') == parent_node]
        if len(matching_rules) > 0:
            return(matching_rules)
        parent_node = amrfp_nodes.get(parent_node)

    # Okay so using the nodeID didn't work, so now we need to check the sequence accession
    matching_rules = [rule for rule in rules_to_check if rule.get('refseq accession') == seq_acc]
    if len(matching_rules) > 0:
        if type == 'Gene presence detected':
            # just return these rules
            return matching_rules
        elif type == 'Protein variant detected' or type == 'Nucleotide variant detected':
            # only return rules with the appropriate mutation
            final_matching_rules = []
            # now we need to check the mutation, extracting any matching rules
            for rule in matching_rules:
                if rule['mutation'] == amrrules_mutation:
                    final_matching_rules.append(rule)
            return final_matching_rules

    #TODO: HMM accession check

    # if nothing matched, then we reurn None
    return None

def annotate_rule(row, rules, annot_opts):
    minimal_columns = ['ruleID', 'context', 'drug', 'drug class', 'phenotype', 'clinical category', 'evidence grade']
    full_columns = ['breakpoint', 'breakpoint standard', 'evidence code', 'evidence limitations', 'PMID', 'rule curation note']

    if rules is None:
        # if we didn't find a matching rule, then we need to add new columns for each of the options but using '-' as the value
        if annot_opts == 'minimal':
            for col in minimal_columns:
                row[col] = '-'
        elif annot_opts == 'full':
            for col in full_columns:
                row[col] = '-'
        return [row]
    # if we found multiple rules, we want to return each rule as its own row in the output file
    if len(rules) > 1:
        # we need to create a new row for each rule
        output_rows = []
        for rule in rules:
            new_row = row.copy()
            if annot_opts == 'minimal':
                for col in minimal_columns:
                    new_row[col] = rule.get(col)
            elif annot_opts == 'full':
                for col in full_columns:
                    new_row[col] = rule.get(col)
            output_rows.append(new_row)
        return output_rows
    # if we found a single rule, we want to annotate the row with the rule info
    else:
        if annot_opts == 'minimal':
            for col in minimal_columns:
                row[col] = rules[0].get(col)
        elif annot_opts == 'full':
            for col in full_columns:
                row[col] = rules[0].get(col)
        return [row]

def write_summary(output_rows):

    # We now want to write a sumamary file that groups hits based on drug class or drug
    # In each drug/drug class, we want to list the highest category and phenotype for that drug/drug class
    # then all the markers that are associated with that drug/drug class, then all the singleton rule IDs
    #TODO: implement combo rules

    drugs_and_classes = [] # list of unique drugs or drug classes to parse
    for row in output_rows:
        drug = row.get('drug')
        drug_class = row.get('drug class')
        drugs_and_classes.append(drug)
        drugs_and_classes.append(drug_class)
    # remove duplicates
    drugs_and_classes = list(set(drugs_and_classes))
    # remove any '-' or '' values
    drugs_and_classes = [x for x in drugs_and_classes if x not in ['-', '']]

    summary_rows = []
    for drug_or_class in drugs_and_classes:
        summarised = {'drug_or_class': drug_or_class} # initalise the drug or class
        # extract all rows that match this drug or class
        matching_rows = [row for row in output_rows if row.get('drug') == drug_or_class or row.get('drug class') == drug_or_class]
        # if there's only one row, then just add the relevant info to the summary
        if len(matching_rows) == 1:
            summarised['category'] = matching_rows[0].get('clinical category')
            summarised['phenotype'] = matching_rows[0].get('phenotype')
            summarised['evidence grade'] = matching_rows[0].get('evidence grade')
            summarised['markers'] = matching_rows[0].get('Gene symbol') or matching_rows[0].get('Element symbol')
            summarised['ruleIDs'] = matching_rows[0].get('ruleID')
            summary_rows.append(summarised)
        # if there are multiple rows, we need to combine some of the info
        elif len(matching_rows) > 1:
            # we need to get the highest category and phenotype
            categories = [row.get('clinical category') for row in matching_rows]
            phenotypes = [row.get('phenotype') for row in matching_rows]
            # get the highest category
            highest_category = max(categories, key=lambda x: ['-', 'S', 'I', 'R'].index(x))
            summarised['category'] = highest_category
            # get the highest phenotype
            highest_phenotype = max(phenotypes, key=lambda x: ['-', 'wildtype', 'nonwildtype'].index(x))
            summarised['phenotype'] = highest_phenotype
            # get the highest evidence grade
            evidence_grades = [row.get('evidence grade') for row in matching_rows]
            highest_evidence_grade = max(evidence_grades, key=lambda x: ['-', 'weak', 'moderate', 'strong'].index(x))
            summarised['evidence grade'] = highest_evidence_grade
            # combine all the markers into a single string
            markers = []
            for row in matching_rows:
                gene_symbol = row.get('Gene symbol') or row.get('Element symbol')
                if gene_symbol not in markers:
                    markers.append(gene_symbol)
            summarised['markers'] = ';'.join(markers)
            # combine all the rule IDs into a single string
            ruleIDs = []
            for row in matching_rows:
                ruleID = row.get('ruleID')
                if ruleID not in ruleIDs:
                    ruleIDs.append(ruleID)
            summarised['ruleIDs'] = ';'.join(ruleIDs)
            summary_rows.append(summarised)
    
    return summary_rows

if __name__ == "__main__":
    args = parse_args()

    if args.amr_tool != 'amrfp':
        raise NotImplementedError("Currently only amrfp is supported. Please use amrfp as the AMR tool.")

    if args.amr_tool == 'amrfp':
        # then we need to grab the refgene heirarchy direct from the ncbi website (get latest for now)
        #TODO: user specifies version of amrfp database they used, or we extract this from hamronized file
        reference_gene_hierarchy_url = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneHierarchy.txt"
        amrfp_nodes = download_and_parse_reference_gene_hierarchy(reference_gene_hierarchy_url)

    # first we need to parse the rules file. For this first draft we just want to extract the gene, nodeID, context, drug, drug class, phenotype, clinical category options
    rules = parse_rules_file(args.rules, args.amr_tool)
    
    #TODO: check if we have multiple genomes in a single input file, if so, process them separately but combine together at the end
    # now we want to read in the amrfinder plus file, and for each row, look for a matching rule
    matched_hits = {} # key is the hierarchy ID, value is set of rules
    unmatched_hits = []
    output_rows = []
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        # Check that the input file has the heirarchy column otherwise return an error and ask user to re-run
        if 'Hierarchy node' not in reader.fieldnames:
            raise ValueError("Input file does not contain 'Hierarchy node' column. Please re-run AMRFinderPlus with the --print_node option to ensure this column is in the output file.")
        
        row_count = 1
        for row in reader:
            matched_rules = check_rules(row, rules, amrfp_nodes)
            if matched_rules is not None:
                # annotate the row with the rule info, based on whether we're using minimal or full annotation
                new_rows = annotate_rule(row, matched_rules, args.annot_opts)
                matched_hits[row_count] = matched_rules
            else:
                new_rows = annotate_rule(row, None, args.annot_opts)
                unmatched_hits.append(row.get('Hierarchy node'))
            row_count += 1
            # add the new rows to the output row list
            output_rows.extend(new_rows)
    
    # write the output files
    interpreted_output_file = os.path.join(args.output_dir, args.output_prefix + '_interpreted.tsv')
    summary_output_file = os.path.join(args.output_dir, args.output_prefix + '_summary.tsv')
    with open(interpreted_output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames + ['ruleID', 'context', 'drug', 'drug class', 'phenotype', 'clinical category', 'evidence grade'], delimiter='\t')
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"{len(matched_hits)} hits matched a rule and {len(unmatched_hits)} hits did not match a rule.")
    print(f"Output written to {interpreted_output_file}.")

    summary_output = write_summary(output_rows)

    with open(summary_output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['drug_or_class', 'category', 'phenotype', 'evidence grade', 'markers', 'ruleIDs'], delimiter='\t')
        writer.writeheader()
        writer.writerows(summary_output)
    print(f"Summary output written to {summary_output_file}.")



    
        