import requests
import csv
import re

# Define input file and output file paths
input_file = 'intact_negative.txt'
output_file = 'concatenated_protein_sequences.fasta'

# UniProt API URL for bulk sequence retrieval
uniprot_url = "https://www.uniprot.org/uniprot/{}.fasta"

# Function to fetch protein sequence from UniProt
def fetch_sequence(uniprot_id):
    try:
        response = requests.get(uniprot_url.format(uniprot_id))
        if response.status_code == 200:
            # Extract sequence from FASTA format
            lines = response.text.split('\n')
            return ''.join(lines[1:])  # Skip the header line
        else:
            print(f"Error fetching {uniprot_id}: {response.status_code}")
            return None
    except Exception as e:
        print(f"Exception for {uniprot_id}: {e}")
        return None

# Function to extract UniProt IDs from interactor columns
def extract_uniprot_id(column_value):
    match = re.match(r"uniprotkb:([\w\-]+)", column_value)
    if match:
        return match.group(1)
    return None

# Process the input file
protein_pairs = set()  # To avoid duplicate pairs
with open(input_file, 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    for row in reader:
        interactor_a = extract_uniprot_id(row['#ID(s) interactor A'])
        interactor_b = extract_uniprot_id(row['ID(s) interactor B'])
        if interactor_a and interactor_b:
            protein_pairs.add((interactor_a, interactor_b))

# Fetch sequences and write concatenated sequences to output file
with open(output_file, 'w') as fasta_file:
    for interactor_a, interactor_b in protein_pairs:
        print(f"Fetching sequences for pair: {interactor_a}, {interactor_b}")
        
        # Fetch sequences for both interactors
        sequence_a = fetch_sequence(interactor_a)
        sequence_b = fetch_sequence(interactor_b)
        
        if sequence_a and sequence_b:
            # Concatenate sequences with an underscore
            concatenated_sequence = sequence_a + "_" + sequence_b
            
            # Write to file in FASTA format
            fasta_file.write(f">Pair_{interactor_a}_{interactor_b}\n")
            fasta_file.write(f"{concatenated_sequence}\n")

print(f"Concatenated sequences saved to {output_file}")
