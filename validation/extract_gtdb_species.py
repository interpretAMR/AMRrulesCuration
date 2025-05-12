# read in the GTDB taxonomy file and extract all unique species names
gtdb_taxonomy_file = 'bac120_taxonomy_r220.tsv'
gtdb_species = []

with open(gtdb_taxonomy_file, 'r') as f:
    for line in f:
        cols = line.strip().split('\t')
        taxonomy_parts = cols[1].split(';')
        # Extract the species name from the taxonomy
        species_name = taxonomy_parts[-1].strip()
        # Check if the species name is already in the list
        if species_name not in gtdb_species:
            gtdb_species.append(species_name)
# Write the unique species names to a file
with open('gtdb_species_r220.txt', 'w') as f:
    for species in gtdb_species:
        f.write(f"{species}\n")