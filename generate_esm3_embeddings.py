import os
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

# Loop through each file in the pdb_files folder
for pdb_file in os.listdir(PDB_FOLDER):
    if pdb_file.endswith(".pdb"):
        file_path = os.path.join(PDB_FOLDER, pdb_file)
        
        # Load the protein from the .pdb file
        protein = ESMProtein.from_pdb(file_path)
        
        # Generate embeddings
        try:
            protein_tensor = client.encode(protein)
            output = client.forward_and_sample(
                protein_tensor, SamplingConfig(return_per_residue_embeddings=True)
            )
            
            # Extract per-residue embeddings and mean-pool to get a fixed-size embedding
            per_residue_embeddings = np.array(output.per_residue_embedding)
            mean_embedding = np.mean(per_residue_embeddings, axis=0)
            
            # Define the output path for this embedding
            embedding_file = os.path.join(OUTPUT_FOLDER, f"{pdb_file}.npy")
            
            # Save the mean-pooled embedding as a .npy file
            np.save(embedding_file, mean_embedding)
            print(f"Processed {pdb_file} successfully. Saved to {embedding_file}")

        except Exception as e:
            print(f"Failed to process {pdb_file}: {e}")
