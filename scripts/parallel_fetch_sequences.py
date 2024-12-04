# Fetches protein sequences from UniProt according to accession IDs found in IntAct
# Also concatenates sequences with an underscore (_) to create dataset

import requests
import csv
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

input_file = 'data/negatome.txt'
output_file = 'data/intact_negative.fasta'

# UniProt API URL
uniprot_url = "https://www.uniprot.org/uniprot/{}.fasta"

# Failure categories for logging purposes
failure_reasons = {
    "missing_sequence_a": 0,
    "missing_sequence_b": 0,
    "missing_both_sequences": 0,
    "network_error": 0,
    "other_errors": 0
}

# Function to fetch protein sequence from UniProt
def fetch_sequence(uniprot_id):
    try:
        response = requests.get(uniprot_url.format(uniprot_id))
        if response.status_code == 200:
            # Extract sequence from FASTA format
            lines = response.text.split('\n')
            return ''.join(lines[1:])
        elif response.status_code == 404:
            print(f"Missing sequence for {uniprot_id}: Status 404")
            return None
        else:
            print(f"Error fetching {uniprot_id}: Status {response.status_code}")
            return None
    except Exception as e:
        print(f"Network error fetching {uniprot_id}: {e}")
        return "network_error"

# Function to process a protein pair and return concatenated sequences
def process_protein_pair(pair):
    interactor_a, interactor_b = pair
    sequence_a = fetch_sequence(interactor_a)
    sequence_b = fetch_sequence(interactor_b)

    if sequence_a == "network_error" or sequence_b == "network_error":
        failure_reasons["network_error"] += 1
        print(f"Failed: Pair {interactor_a}, {interactor_b} (Network error)")
        return None

    if not sequence_a and not sequence_b:
        failure_reasons["missing_both_sequences"] += 1
        print(f"Failed: Pair {interactor_a}, {interactor_b} (Both sequences missing)")
        return None
    elif not sequence_a:
        failure_reasons["missing_sequence_a"] += 1
        print(f"Failed: Pair {interactor_a}, {interactor_b} (Missing sequence A)")
        return None
    elif not sequence_b:
        failure_reasons["missing_sequence_b"] += 1
        print(f"Failed: Pair {interactor_a}, {interactor_b} (Missing sequence B)")
        return None

    # Concatenate sequences if both are valid
    concatenated_sequence = sequence_a + "_" + sequence_b
    print(f"Completed: Pair {interactor_a}, {interactor_b}")
    return f">{interactor_a}_{interactor_b}\n{concatenated_sequence}\n"

# Extract protein pairs from the input file
def extract_protein_pairs(file_path):
    protein_pairs = set()  # To avoid duplicate pairs
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            interactor_a = extract_uniprot_id(row['#ID(s) interactor A'])
            interactor_b = extract_uniprot_id(row['ID(s) interactor B'])
            if interactor_a and interactor_b:
                protein_pairs.add((interactor_a, interactor_b))
    return protein_pairs

# Function to extract UniProt IDs from interactor columns
def extract_uniprot_id(column_value):
    match = re.match(r"uniprotkb:([\w\-]+)", column_value)
    if match:
        return match.group(1)
    return None

# Main execution
if __name__ == "__main__":
    protein_pairs = extract_protein_pairs(input_file)
    total_pairs = len(protein_pairs)
    print(f"Initial protein pairs count: {total_pairs}")

    # Counters for successes
    success_count = 0

    # Parallelize for faster processing
    with ThreadPoolExecutor(max_workers=60) as executor:
        futures = {executor.submit(process_protein_pair, pair): pair for pair in protein_pairs}
        with open(output_file, 'w') as fasta_file:
            for future in as_completed(futures):
                try:
                    result = future.result()
                    if result:
                        fasta_file.write(result)
                        success_count += 1
                except Exception as e:
                    pair = futures[future]
                    print(f"Unexpected error with pair {pair}: {e}")
                    failure_reasons["other_errors"] += 1

    # Print summary of results
    print(f"\nProcessing completed!")
    print(f"Total pairs: {total_pairs}")
    print(f"Successful pairs: {success_count}")
    print(f"Failed pairs: {total_pairs - success_count}")
    print("\nFailure Breakdown:")
    for reason, count in failure_reasons.items():
        print(f"{reason.replace('_', ' ').capitalize()}: {count}")

    print(f"\nConcatenated sequences saved to {output_file}")
