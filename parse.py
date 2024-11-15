import os
import sys

def read_list_of_prots(file_path):
    """
    Reads list_of_prots.txt and returns a dictionary mapping (pdb_id, chain_id) to uniprot_id
    and a set of unique (pdb_id, chain_id) tuples.
    
    Returns:
        tuple: (dict of PDB mappings, set of unique (pdb_id, chain_id) pairs)
    """
    pdb_chain_to_uniprot = {}
    unique_pdb_chains = set()  # Store (pdb_id, chain_id) pairs
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
    Clean up the concatenated sequences file by removing entries with zero-length sequences.
    Outputs only the concatenated sequence and interaction type.
    
    Args:
        input_file (str): Path to original concatenated sequences file
        output_file (str): Path to cleaned output file
    
    Returns:
        int: Number of sequences removed
    """
    removed_count = 0
    valid_sequences = []
    
    # Read the original file
    with open(input_file, 'r') as f:
        f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 10:  # Make sure we have all fields
                total_length = int(parts[8]) if parts[8].isdigit() else 0
                if total_length > 0:
                    sequence = parts[7]  # Concatenated sequence
                    interaction_type = parts[9]  # Interaction type
                    valid_sequences.append((sequence, interaction_type))
                else:
                    removed_count += 1
    
    # Write the cleaned file with only sequence and class
    with open(output_file, 'w') as f:
        f.write("Sequence\tClass\n")  # New simplified header
        for sequence, interaction_type in valid_sequences:
            f.write(f"{sequence}\t{interaction_type}\n")
            
    return removed_count


def main():
    # File paths
    list_of_prots_file = 'list_of_prots.txt'
    interactions_file = 'interactions_data.txt'
    pdb_dir = 'pdb_files'
    output_file = 'concatenated_sequences.txt'
    individual_sequences_file = 'individual_sequences.txt'
    cleaned_output_file = 'concatenated_sequences_cleaned.txt'


    # Step 1: Read list_of_prots.txt and get unique PDB-chain pairs
    print("Reading list_of_prots.txt...")
    pdb_chain_to_uniprot, unique_pdb_chains = read_list_of_prots(list_of_prots_file)
    print(f"Found {len(unique_pdb_chains)} unique PDB-chain pairs in list_of_prots.txt")

    # Step 2: Read interactions_data.txt
    print("Reading interactions_data.txt...")
    interactions = read_interactions(interactions_file)
    print(f"Found {len(interactions)} interactions")

    # Step 3: Extract sequences for each unique PDB-chain pair
    sequences = {}  # Will store sequences keyed by (pdb_id, chain_id)
    with open(individual_sequences_file, 'w') as seq_f:
        seq_f.write("PDB_ID\tChain_ID\tUniprot_ID\tSequence\tLength\n")  # Header
        
        for pdb_id, chain_id in unique_pdb_chains:
            uniprot_id = pdb_chain_to_uniprot[(pdb_id, chain_id)]
            print(f"Extracting sequence for PDB ID: {pdb_id}, Chain: {chain_id}...")
            seq = get_sequence_from_pdb(pdb_id, chain_id, pdb_dir)
            
            # Store sequence with (pdb_id, chain_id) as key
            sequences[(pdb_id, chain_id)] = seq
            
            # Write to individual sequences file
            seq_length = len(seq) if seq else 0
            seq_f.write(f"{pdb_id}\t{chain_id}\t{uniprot_id}\t{seq}\t{seq_length}\n")

# Step 4: Concatenate sequences based on interactions
    print(f"Concatenating sequences and writing to {output_file}...")
    with open(output_file, 'w') as out_f:
        # Modified header to include interaction type
        out_f.write("Interaction_Pair\tPDB1\tChain1\tPDB2\tChain2\tPDB1_Length\tPDB2_Length\tConcatenated_Sequence\tTotal_Length\tInteraction_Type\n")
        
        for interaction in interactions:
            pdb1, pdb2, interaction_type = interaction
            
            # Find the corresponding chain IDs from our mapping
            pdb1_chains = [chain for (pid, chain) in sequences.keys() if pid == pdb1]
            pdb2_chains = [chain for (pid, chain) in sequences.keys() if pid == pdb2]
            
            if not pdb1_chains or not pdb2_chains:
                missing = []
                if not pdb1_chains:
                    missing.append(pdb1)
                if not pdb2_chains:
                    missing.append(pdb2)
                print(f"Warning: PDB(s) {', '.join(missing)} from interaction not found in list_of_prots.txt", file=sys.stderr)
                continue

            # Use the first chain found for each PDB
            chain1 = pdb1_chains[0]
            chain2 = pdb2_chains[0]
            
            seq1 = sequences.get((pdb1, chain1), '')
            seq2 = sequences.get((pdb2, chain2), '')
            len1 = len(seq1)
            len2 = len(seq2)
            
            if seq1 and seq2:
                concatenated_seq = seq1 + " " + seq2
                total_length = len1 + len2
                out_f.write(f"{pdb1}_{pdb2}\t{pdb1}\t{chain1}\t{pdb2}\t{chain2}\t{len1}\t{len2}\t{concatenated_seq}\t{total_length}\t{interaction_type}\n")
            else:
                print(f"Warning: Missing sequence for pair {pdb1}_{pdb2}.", file=sys.stderr)
                out_f.write(f"{pdb1}_{pdb2}\t{pdb1}\t{chain1}\t{pdb2}\t{chain2}\t{len1}\t{len2}\t\t0\t{interaction_type}\n")


    print(f"Processing complete. Results saved in:")
    print(f"1. Individual sequences: {individual_sequences_file}")
    print(f"2. Concatenated sequences: {output_file}")

    print("\nCleaning up concatenated sequences...")
    removed_count = cleanup_concatenated_sequences(output_file, cleaned_output_file)
    print(f"Removed {removed_count} sequences with length 0")
    print(f"3. Cleaned concatenated sequences: {cleaned_output_file}")


if __name__ == '__main__':
    main()