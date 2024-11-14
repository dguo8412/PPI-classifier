from Bio import PDB
import os

def extract_sequence_from_pdb(pdb_file):
    """
    Extract amino acid sequence from a PDB file.
    
    Args:
        pdb_file (str): Path to PDB file
    
    Returns:
        str: Amino acid sequence
    """
    # Dictionary to convert three letter amino acid codes to one letter
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_file)
        
        # Get the first model and chain
        model = structure[0]
        chain = next(iter(model))
        
        # Extract sequence
        sequence = ''
        for residue in chain:
            res_name = residue.get_resname()
            if res_name in aa_dict:
                sequence += aa_dict[res_name]
                
        return sequence
    except Exception as e:
        print(f"Error processing {pdb_file}: {str(e)}")
        return None

def process_interaction_data(interaction_file, pdb_folder):
    """
    Process interaction data and create sequences dataset.
    
    Args:
        interaction_file (str): Path to interaction data file
        pdb_folder (str): Path to folder containing PDB files
    
    Returns:
        tuple: Lists of concatenated sequences and labels
    """
    sequences = []
    labels = []
    
    with open(interaction_file, 'r') as f:
        for line in f:
            # Split the line into components
            prot1_id, prot2_id, label = line.strip().split()
            
            # Construct file paths
            pdb1_path = os.path.join(pdb_folder, f"{prot1_id}.pdb")
            pdb2_path = os.path.join(pdb_folder, f"{prot2_id}.pdb")
            
            # Check if both files exist
            if not (os.path.exists(pdb1_path) and os.path.exists(pdb2_path)):
                print(f"Missing PDB file for {prot1_id} or {prot2_id}")
                continue
            
            # Extract sequences
            seq1 = extract_sequence_from_pdb(pdb1_path)
            seq2 = extract_sequence_from_pdb(pdb2_path)
            
            if seq1 and seq2:
                # Concatenate sequences with space separator
                concat_seq = f"{seq1} {seq2}"
                sequences.append(concat_seq)
                labels.append(int(label))
            else:
                print(f"Failed to extract sequence for {prot1_id} or {prot2_id}")
    
    return sequences, labels

def save_processed_data(sequences, labels, output_file):
    """
    Save processed sequences and labels to a file.
    
    Args:
        sequences (list): List of concatenated sequences
        labels (list): List of interaction labels
        output_file (str): Path to output file
    """
    with open(output_file, 'w') as f:
        for seq, label in zip(sequences, labels):
            f.write(f"{seq}\t{label}\n")

# Example usage
if __name__ == "__main__":
    interaction_file = "interaction_data.txt"
    pdb_folder = "pdb_files"
    output_file = "processed_sequences.txt"
    
    # Process the data
    sequences, labels = process_interaction_data(interaction_file, pdb_folder)
    
    # Save to file
    save_processed_data(sequences, labels, output_file)
    
    # Print summary
    print(f"Processed {len(sequences)} interaction pairs")
    print(f"Positive interactions: {sum(labels)}")
    print(f"Negative interactions: {len(labels) - sum(labels)}")