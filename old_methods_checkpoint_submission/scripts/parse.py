import os
import sys
import pandas as pd


def read_list_of_prots(file_path):
    """
    Reads list_of_prots.txt and returns a dictionary mapping (pdb_id, chain_id) to uniprot_id
    and a set of unique (pdb_id, chain_id) tuples.
    """
    pdb_chain_to_uniprot = {}
    unique_pdb_chains = set()
    try:
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    uniprot_id, pdb_id, chain = parts[0], parts[1], parts[2]
                    pdb_chain_to_uniprot[(pdb_id, chain)] = uniprot_id
                    unique_pdb_chains.add((pdb_id, chain))
    except FileNotFoundError:
        print(f"Error: {file_path} not found.", file=sys.stderr)
        sys.exit(1)
    return pdb_chain_to_uniprot, unique_pdb_chains
    

def read_interactions(file_path):
    """
    Reads interactions_data.txt and returns a list of tuples (pdb1, pdb2, interaction_type).
    """
    interactions = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    pdb1, pdb2, interaction_type = parts[0], parts[1], parts[2]
                    interactions.append((pdb1, pdb2, interaction_type))
    except FileNotFoundError:
        print(f"Error: {file_path} not found.", file=sys.stderr)
        sys.exit(1)
    return interactions


def get_sequence_from_pdb(pdb_id: str, chain_id: str, pdb_dir: str) -> str:
    """
    Extract amino acid sequence from a PDB file for a specified chain using SEQRES records.
    
    Args:
        pdb_id (str): PDB identifier
        chain_id (str): Chain identifier (e.g., 'A', 'B', etc.)
        pdb_dir (str): Directory containing PDB files
        
    Returns:
        str: Amino acid sequence in one-letter code
    """
    # Dictionary for converting three-letter amino acid codes to one-letter codes
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        # DNA/RNA nucleotides
        'DA': 'A', 'DT': 'T', 'DG': 'G', 'DC': 'C',
        'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C',
        'U': 'U'
    }
    
    try:
        with open(f"{pdb_dir}/{pdb_id}.pdb", 'r') as f:
            # Check for HTML error page
            first_line = f.readline().strip()
            if first_line.startswith('<!DOCTYPE') or first_line.startswith('<html'):
                print(f"Warning: PDB file {pdb_id}.pdb contains HTML error message", file=sys.stderr)
                return ''
                
            # Go back to start of file
            f.seek(0)
            
            sequence = []
            in_seqres_block = False
            
            for line in f:
                if line.startswith('SEQRES'):
                    current_chain = line[11].strip()
                    if current_chain == chain_id:
                        in_seqres_block = True
                        # Extract residues from the SEQRES line (positions 19-70)
                        residues = line[19:70].split()
                        for res in residues:
                            if res in aa_dict:
                                sequence.append(aa_dict[res])
                            else:
                                sequence.append('X')  # Unknown residue
                elif line.startswith('ATOM') and not in_seqres_block:
                    # If no SEQRES records found, fall back to ATOM records
                    current_chain = line[21].strip()
                    if current_chain == chain_id:
                        residue = line[17:20].strip()
                        residue_num = line[22:26].strip()
                        if residue in aa_dict and residue_num != prev_residue_num:
                            sequence.append(aa_dict[residue])
                            prev_residue_num = residue_num
            
            if sequence:
                return ''.join(sequence)
            else:
                print(f"Warning: No sequence found for chain {chain_id} in {pdb_id}.pdb", file=sys.stderr)
                return ''
                
    except FileNotFoundError:
        print(f"Warning: PDB file not found: {pdb_id}.pdb", file=sys.stderr)
        return ''
    except Exception as e:
        print(f"Warning: Error processing {pdb_id}.pdb: {str(e)}", file=sys.stderr)
        return ''



def cleanup_concatenated_sequences(input_file, output_file):
    """
    Clean up concatenated sequences file by removing zero-length sequences.
    Outputs only the concatenated sequence and interaction type to a .csv file.
    """
    df = pd.read_csv(input_file)
    cleaned_df = df[df['Total_Length'] > 0]  # Keep only rows with positive sequence length
    removed_count = len(df) - len(cleaned_df)
    cleaned_df[['Concatenated_Sequence', 'Interaction_Type']].to_csv(output_file, index=False)
    return removed_count



def main():
    # File paths
    list_of_prots_file = 'list_of_prots.txt'
    interactions_file = 'interactions_data.txt'
    pdb_dir = 'pdb_files'
    output_file = 'concatenated_sequences.csv'
    individual_sequences_file = 'individual_sequences.csv'
    cleaned_output_file = 'concatenated_sequences_cleaned.csv'

    # Step 1: Read list_of_prots.txt
    print("Reading list_of_prots.txt...")
    pdb_chain_to_uniprot, unique_pdb_chains = read_list_of_prots(list_of_prots_file)
    print(f"Found {len(unique_pdb_chains)} unique PDB-chain pairs.")

    # Step 2: Read interactions_data.txt
    print("Reading interactions_data.txt...")
    interactions = read_interactions(interactions_file)
    print(f"Found {len(interactions)} interactions.")

    # Step 3: Extract sequences for each unique PDB-chain pair
    sequences = []
    for pdb_id, chain_id in unique_pdb_chains:
        uniprot_id = pdb_chain_to_uniprot[(pdb_id, chain_id)]
        print(f"Extracting sequence for PDB ID: {pdb_id}, Chain: {chain_id}...")
        seq = get_sequence_from_pdb(pdb_id, chain_id, pdb_dir)
        sequences.append({'PDB_ID': pdb_id, 'Chain_ID': chain_id, 'Uniprot_ID': uniprot_id, 'Sequence': seq, 'Length': len(seq)})

    # Save individual sequences to CSV
    pd.DataFrame(sequences).to_csv(individual_sequences_file, index=False)
    print(f"Saved individual sequences to {individual_sequences_file}.")

    # Step 4: Concatenate sequences based on interactions
    concatenated_sequences = []
    for pdb1, pdb2, interaction_type in interactions:
        seq1 = next((s['Sequence'] for s in sequences if s['PDB_ID'] == pdb1), '')
        seq2 = next((s['Sequence'] for s in sequences if s['PDB_ID'] == pdb2), '')
        concatenated_sequences.append({
            'Interaction_Pair': f"{pdb1}_{pdb2}",
            'PDB1': pdb1, 'PDB2': pdb2,
            'PDB1_Length': len(seq1), 'PDB2_Length': len(seq2),
            'Concatenated_Sequence': seq1 + seq2,
            'Total_Length': len(seq1) + len(seq2),
            'Interaction_Type': interaction_type
        })

    # Save concatenated sequences to CSV
    pd.DataFrame(concatenated_sequences).to_csv(output_file, index=False)
    print(f"Saved concatenated sequences to {output_file}.")

    # Step 5: Clean up concatenated sequences
    removed_count = cleanup_concatenated_sequences(output_file, cleaned_output_file)
    print(f"Removed {removed_count} sequences with zero length.")
    print(f"Saved cleaned concatenated sequences to {cleaned_output_file}.")

if __name__ == '__main__':
    main()