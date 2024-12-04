# Generates ESM embeddings given input fasta files containing protein sequences


import os
import pandas as pd
import torch
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, SamplingConfig
from huggingface_hub import login
import numpy as np
from tqdm import tqdm
from esm.utils.constants.models import ESM3_OPEN_SMALL

def read_fasta_pairs(file_path):
    pairs = {}
    current_id = None
    current_seq = ''
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    pairs[current_id] = current_seq
                current_id = line[1:] 
                current_seq = ''
            else:
                current_seq += line
        if current_id:
            pairs[current_id] = current_seq
            
    # Convert to DataFrame with split sequences
    data = []
    for pair_id, seq in pairs.items():
        sequences = seq.rsplit('_', 1)
        if len(sequences) != 2:
            print(f"Warning: Skipping malformed entry {pair_id}")
            continue
            
        seq1, seq2 = sequences
        
        data.append({
            'pair_id': pair_id,
            'seq1': seq1,
            'seq2': seq2
        })
    return pd.DataFrame(data)


DEVICE = "cuda"
login()
model: ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to(DEVICE)
OUTPUT_FOLDER = "positive_protein_embeddings/"
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# Read paired sequences
df = read_fasta_pairs('../data/positive_pairs.fasta')

# Process each sequence pair
with torch.no_grad():
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing protein pairs"):
        pair_id = row['pair_id']
        output_path = os.path.join(OUTPUT_FOLDER, f"{pair_id}.pt")
        
        if os.path.exists(output_path):
            print(f"Skipping {pair_id}: embedding already exists")
            continue
            
        try:
            # Process first sequence
            protein1 = ESMProtein(sequence=row['seq1'])
            protein_tensor1 = model.encode(protein1).to(DEVICE)
            output1 = model.forward_and_sample(
                protein_tensor1,
                SamplingConfig(return_per_residue_embeddings=True)
            )
            embedding1 = output1.per_residue_embedding.mean(dim=0)
            
            # Process second sequence
            protein2 = ESMProtein(sequence=row['seq2'])
            protein_tensor2 = model.encode(protein2).to(DEVICE)
            output2 = model.forward_and_sample(
                protein_tensor2,
                SamplingConfig(return_per_residue_embeddings=True)
            )
            embedding2 = output2.per_residue_embedding.mean(dim=0)
            
            # Concatenate embeddings
            combined_embedding = torch.cat([embedding1, embedding2], dim=0)
            
            # Save combined embedding
            torch.save(combined_embedding, output_path)
            
        except Exception as e:
            print(f"\nFailed to process {pair_id}: {str(e)}")
            continue

print("\nProcessing complete!")
print(f"Embeddings saved to: {OUTPUT_FOLDER}/")

# Print example of how to load embeddings
print("\nTo load embeddings for a specific protein pair:")
print("embedding = torch.load(f'{OUTPUT_FOLDER}/PAIR_ID.pt')")

# Optional: Print embedding tensor shape for verification
if os.listdir(OUTPUT_FOLDER):
    example_file = os.path.join(OUTPUT_FOLDER, os.listdir(OUTPUT_FOLDER)[0])
    example_embedding = torch.load(example_file)
    print(f"\nCombined embedding shape: {example_embedding.shape}")
    print(f"Individual sequence embedding shape: {example_embedding.shape[0]//2}")