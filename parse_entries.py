from Bio import SeqIO
import os
import warnings
warnings.filterwarnings('ignore')

def extract_sequence_pdb(pdb_file, chain='0'):
    """
    Extract amino acid sequence using SeqIO
    """
    try:
        # Try using SeqIO first
        for record in SeqIO.parse(pdb_file, "pdb-seqres"):
            return str(record.seq)
    except Exception as e:
        print(f"Error processing {pdb_file}: {str(e)}")
        return None

def process_interaction_data(interaction_file, pdb_folder, missing_file="missing_structures.txt"):
    """
    Process interaction data and create sequences dataset.
    
    Args:
        interaction_file (str): Path to file containing interaction pairs and labels
        pdb_folder (str): Path to folder containing PDB files
        missing_file (str): Path to output file for logging missing/problematic structures
    
    Returns:
        tuple: (sequences, labels) where sequences is a list of concatenated protein sequences
               and labels is a list of corresponding interaction labels (0 or 1)
    """
    sequences = []
    labels = []
    missing_pairs = []
    processed_count = 0
    total_pairs = sum(1 for line in open(interaction_file))
    
    with open(interaction_file, 'r') as f:
        for line in f:
            processed_count += 1
            if processed_count % 100 == 0:
                print(f"Processing pair {processed_count}/{total_pairs}")
                
            prot1_id, prot2_id, label = line.strip().split()
            
            pdb1_path = os.path.join(pdb_folder, f"{prot1_id}.pdb")
            pdb2_path = os.path.join(pdb_folder, f"{prot2_id}.pdb")
            
            # Check file existence
            if not os.path.exists(pdb1_path) or not os.path.exists(pdb2_path):
                missing = []
                if not os.path.exists(pdb1_path):
                    missing.append(prot1_id)
                if not os.path.exists(pdb2_path):
                    missing.append(prot2_id)
                missing_pairs.append((prot1_id, prot2_id, missing))
                continue
            
            try:
                # Extract sequences for both proteins
                seq1 = extract_sequence_pdb(pdb1_path)
                seq2 = extract_sequence_pdb(pdb2_path)
                
                if seq1 and seq2:
                    concat_seq = f"{seq1} {seq2}"
                    sequences.append(concat_seq)
                    labels.append(int(label))
                else:
                    missing = []
                    if not seq1:
                        missing.append(prot1_id)
                    if not seq2:
                        missing.append(prot2_id)
                    missing_pairs.append((prot1_id, prot2_id, missing))
                    
            except Exception as e:
                print(f"Error processing {prot1_id} or {prot2_id}: {str(e)}")
                missing_pairs.append((prot1_id, prot2_id, [prot1_id, prot2_id]))
    
    # Log missing structures
    if missing_pairs:
        with open(missing_file, 'w') as f:
            f.write("The following structures had issues:\n")
            for prot1, prot2, missing in missing_pairs:
                f.write(f"{prot1}\t{prot2}\tMissing: {', '.join(missing)}\n")
    
    return sequences, labels

def save_processed_data(sequences, labels, output_file):
    """
    Save processed sequences and labels to a file.
    
    Args:
        sequences (list): List of concatenated protein sequences
        labels (list): List of interaction labels
        output_file (str): Path to output file
    """
    with open(output_file, 'w') as f:
        for seq, label in zip(sequences, labels):
            f.write(f"{seq}\t{label}\n")

def main():
    """
    Main function to run the sequence extraction and processing pipeline.
    """
    # Define input/output paths
    interaction_file = "interaction_data.txt"
    pdb_folder = "pdb_files"
    output_file = "processed_sequences.txt"
    missing_file = "missing_structures.txt"
    
    # Process the data
    print("Starting sequence extraction...")
    sequences, labels = process_interaction_data(interaction_file, pdb_folder, missing_file)
    
    # Save to file
    print("\nSaving processed sequences...")
    save_processed_data(sequences, labels, output_file)
    
    # Print summary
    print("\nProcessing Summary:")
    print(f"Successfully processed {len(sequences)} interaction pairs")
    print(f"Positive interactions: {sum(labels)}")
    print(f"Negative interactions: {len(labels) - sum(labels)}")
    
    # Calculate success rate
    total_processed = len(sequences)
    total_attempted = total_processed + sum(1 for _ in open(missing_file)) - 1  # -1 for header
    if total_attempted > 0:
        success_rate = (total_processed / total_attempted) * 100
        print(f"Success rate: {success_rate:.2f}%")
    
    print(f"\nResults saved to: {output_file}")
    print(f"Missing structures logged to: {missing_file}")

if __name__ == "__main__":
    main()