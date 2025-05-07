import argparse
import csv
import urllib.request
import io

def parse_args():
    parser = argparse.ArgumentParser(description="Interpreter for AMRrules.")

    parser.add_argument('--input', type=str, required=True, help='Path to the input file.')
    parser.add_argument('--output', type=str, required=True, help='Path to the output file.')
    parser.add_argument('--rules', '-r', type=str, required=True, help='Path to the rules file.')
    #TODO: implement card and resfinder options, currently only amrfp is supported
    parser.add_argument('--amr_tool', '-t', type=str, default='amrfp', help='AMR tool used to detect genotypes: options are amrfp, card, resfinder. Currently only amrfp is supported.')
    parser.add_argument('--hamronized', '-H', action='store_true', help='Input file has been hamronized')
    parser.add_argument('--amrfp_db_version', type=str, default='latest', help='Version of the AMRFP database used. Default is latest.')
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
    rules = {} # key: relevant ID for geno method, value: values for this rule
    if tool == 'amrfp':
        key_value = 'nodeID'
    with open(rules_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rules[row.get(key_value)] = row
    return rules

def check_rules(hierarchy_id, rules, amrfp_nodes):
    # check if the hierarchy ID is in the rules
    if hierarchy_id in rules.keys():
        # then we have a match, so we can extract the rule
        rule = rules[hierarchy_id]
        return rule   
    else:
        # no match, so we need to check the parent nodes
        parent_node = amrfp_nodes.get(hierarchy_id)
        while parent_node is not None and parent_node != 'AMR':
            if parent_node in rules:
                rule = rules[parent_node]
                break
            parent_node = amrfp_nodes.get(parent_node)
        else:
            # no match found in parent nodes
            return None
        return rule

def annotate_rule(row, rule, annot_opts):
    minimal_columns = ['ruleID', 'context', 'drug', 'drug class', 'phenotype', 'clinical category', 'evidence grade']
    full_columns = ['breakpoint', 'breakpoint standard', 'evidence code', 'evidence limitations', 'PMID', 'rule curation note']

    if rule is None:
        # if we didn't find a matching rule, then we need to add new columns for each of the options but using '-' as the value
        if annot_opts == 'minimal':
            for col in minimal_columns:
                row[col] = '-'
        elif annot_opts == 'full':
            for col in full_columns:
                row[col] = '-'
        return row
    else:
        if annot_opts == 'minimal':
            for col in minimal_columns:
                row[col] = rule.get(col)
        elif annot_opts == 'full':
            for col in full_columns:
                row[col] = rule.get(col)
        return row

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
        
    # now we want to read in the amrfinder plus file, and for each row, look for a matching rule
    # logic for finding a rule match:
    # 1. hierarchyID exatch match to nodeID = extract rule and apply
    # 2. no exact match to hierarchyID, lookup parent ID (all the way to 'AMR' root) and check if there is an exact match with a parent nodeID
    # 3. otherwise no match so do not apply a rule
    matched_hits = []
    unmatched_hits = []
    output_rows = []
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        # Check that the input file has the heirarchy column otherwise return an error and ask user to re-run
        if 'Hierarchy node' not in reader.fieldnames:
            raise ValueError("Input file does not contain 'Hierarchy node' column. Please re-run AMRFinderPlus with the --print_node option to ensure this column is in the output file.")
        for row in reader:
            hierarchy_id = row.get('Hierarchy node')
            matched_rule = check_rules(hierarchy_id, rules, amrfp_nodes)
            if matched_rule is not None:
                # annotate the row with the rule info, based on whether we're using minimal or full annotation
                new_row = annotate_rule(row, matched_rule, args.annot_opts)
                matched_hits.append(hierarchy_id)
            else:
                new_row = annotate_rule(row, None, args.annot_opts)
                unmatched_hits.append(hierarchy_id)
            output_rows.append(new_row)
    # write the output file
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames + ['ruleID', 'context', 'drug', 'drug class', 'phenotype', 'clinical category', 'evidence grade'], delimiter='\t')
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Matched {len(matched_hits)} hits and {len(unmatched_hits)} unmatched hits.")
    print(f"Output written to {args.output}.")
    
        