### Automatic validation of rules

`validate_rules.py`, found in this folder, is a simple python script that will perform a number of automatic checks on draft rulesets, to quickly identify typos or invalid entries.

In order to run this code, your rules file needs to be in **tab-delimited format**. Ensure that the first row in the file lists the column names, as per the [latest version of the specification](https://docs.google.com/spreadsheets/d/1F-J-_8Kyo3W0Oh6eDYyd0N8ahqVwiddM2112-Fg1gKc/edit?usp=sharing). The validation script will notify you if it is unable to find the columns it's expecting, and skip those columns when performing the automatic checks.

To use, download this GitHub repository. Navigate to the `validation` folder, and run:
```
./validate_rules.py /path/to/your/rules_file.txt
```

**Note that the script MUST be run from within its home directory (AMRrulesCuration/validation) as it depends on other files present in this location.**

Results will be printed to stdout. The script will list each each and provide an explanation for why particular rows have failed the check. The final summary lists the checks that have been passed or failed.

* *ruleID*: checks all values are a unique ruleID with contain the same prefix. Will print all detected prefixes to stdout
* *organism*: checks that all values start with `s__`. Will print all unique organism names to stdout out
* *gene*: checks all rows have a value that is not 'NA' or '-'
* *combinatorial rules*: if the gene column contains a combinatorial rule, checks that all ruleIDs are already defined in the rules file
* *nodeID, refseq accession, GenBank accession, HMM accession*: checks that each row contains a value in at least one of these columns. If accessions are present, compares to the current versions of those catalogues to ensure the accessions exist. Note that for `refseq accession` and `GenBank accession`, the validation script is only comparing to accessions listed in the AMRFinderPlus ReferenceGeneHierarchy. Therefore it will flag accessions that aren't present in the AMRFinderPlus table, but do in fact exist in the NCBI databases.
* *ARO accession*: checks the value starts with `ARO:` and that the accession is a valid ARO accession according to the latest CARD ontology
* *mutation*: checks that the row isn't empty or NA for this column (must be '-' if no mutation needs to be specified)
* *variation type*: checks the value is one of the allowable values according to the spec
* *mutation and variation type*: checks that the values in mutation and varaition type are compatiable: eg, if variation type is `Gene presence detected`, mutation should be `-`. If variation type is `Nucleotide variant detected`, the value in mutation should start with `c.`
* *context*: checks that the value is either `core` or `acquired`, cannot be empty, NA, or `-`
* *drug and drug class*: checks that at least one of these columns contains a value, cannot both be empty, NA, or `-`
* *phenotype*: checks that the value is either `wildtype` or `nonwildtype`, cannot be empty, NA, or `-`
* *clinical category*: checks that the value is either `S`, `I`, or `R`
* *breakpoint*: checks that there is a value that isn't empty, NA, or `-`
* *clinical category and breakpoint*: checks if clinical category is `S`, then breakpoint must start with `MIC <`, `MIC <=` or `disk >`. Reverse is applicable for clinical category of `R`. `not applicable` is also a permitted value in breakpoint, if the rule is associated with an expected resistance.
* *breakpoint standard*: checks the value isn't empty, NA, or `-`. Will print all unique values in this column to stdout
* *PMID*: checks the value isn't empty, NA, or `-`, must contain a PMID
* *evidence code*: checks that the value is one of the current suggested ECO codes. If it isn't, will print a list of new ECO codes. Value cannot be NA, `-` or empty. If multiple evidence codes are present, they must be separated by a `,`. All evidence codes must start with `ECO:`
* *evidence grade and limitations*: checks that evidence grade has a value, cannot be NA, `-` or empty. Evidence grade must be either `strong`, `moderate`, or `weak`. If `moderate` or `weak`, then evidence limitations must be filled in with one of the allowable values as per the specification
