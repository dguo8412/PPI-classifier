# client = esm.sdk.client(MODEL_NAME, url="https://forge.evolutionaryscale.ai", token=FORGE_TOKEN)
import os
import pandas as pd
import torch
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, SamplingConfig
import numpy as np

# Set the model parameters and device
MODEL_NAME = "esm3-large-2024-03"  # Replace with the specific model name if available
DEVICE = "cuda"  # Use 'cpu' if GPU is not available
PDB_FOLDER = "pdb_files/"
OUTPUT_FILE = "protein_embeddings.h5"
FORGE_TOKEN = ""

# Initialize the ESM3 client with the proprietary model
client = ESM3.from_pretrained(MODEL_NAME, device=DEVICE)

# Directory where your .pdb files are stored and where embeddings will be saved
PDB_FOLDER = "pdb_files/"
OUTPUT_FOLDER = "protein_embeddings/"

# Ensure the output folder exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# Read sequences from CSV
df = pd.read_csv(INPUT_CSV)

# Keep track of processed proteins
processed_proteins = set(
    os.path.splitext(f)[0] 
    for f in os.listdir(OUTPUT_FOLDER) 
    if f.endswith('.pt')
)

# Process each sequence
with torch.no_grad():
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing proteins"):
        pdb_id = row['pdb_id']
        chain_id = row['chain_id']
        sequence = row['sequence']
        
        # Create filename
        protein_id = f"{pdb_id}_{chain_id}"
        output_path = os.path.join(OUTPUT_FOLDER, f"{protein_id}.pt")
        
        # Skip if already processed
        if protein_id in processed_proteins:
            continue
        
        try:
            # Create ESM protein object and generate embeddings
            protein = ESMProtein(sequence=sequence)
            protein_tensor = client.encode(protein).to(DEVICE)
            
            output = client.forward_and_sample(
                protein_tensor, 
                SamplingConfig(return_per_residue_embeddings=True)
            )
            
            # Mean pool the embeddings
            embedding = output.per_residue_embedding.mean(dim=0)
            
            # Save individual embedding
            torch.save(embedding, output_path)
            processed_proteins.add(protein_id)
            
        except Exception as e:
            print(f"\nFailed to process {protein_id}: {str(e)}")
            continue
        
        # Optional: Clear CUDA cache periodically
        if torch.cuda.is_available() and len(processed_proteins) % 50 == 0:
            torch.cuda.empty_cache()

print("\nProcessing complete!")
print(f"Embeddings saved to: {OUTPUT_FOLDER}/")

# Print example of how to load embeddings
print("\nTo load embeddings for a specific protein:")
print("embedding = torch.load(f'{OUTPUT_FOLDER}/PDB_ID_CHAIN.pt')")

# Optional: Print embedding tensor shape for verification
if processed_proteins:
    example_file = os.path.join(OUTPUT_FOLDER, f"{next(iter(processed_proteins))}.pt")
    example_embedding = torch.load(example_file)
    print(f"\nEmbedding shape: {example_embedding.shape}")