import pandas as pd
from sklearn.model_selection import train_test_split
import os

def split_data(input_file, output_dir, train_size=0.8, val_size=0.1, test_size=0.1, random_state=42):
    """
    Split the cleaned sequence data into train, dev, and test sets.
    Adjusted for a dataset of ~10,000 samples.
    
    Args:
        input_file (str): Path to cleaned concatenated sequences file
        output_dir (str): Directory to save the split files
        train_size (float): Proportion of data for training (default: 0.8)
        val_size (float): Proportion of data for validation (default: 0.1)
        test_size (float): Proportion of data for testing (default: 0.1)
        random_state (int): Random seed for reproducibility
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read the data
    print(f"Reading data from {input_file}...")
    df = pd.read_csv(input_file, sep='\t')
    total_samples = len(df)
    
    print(f"\nTotal dataset size: {total_samples} samples")
    print(f"Expected split sizes:")
    print(f"Training: ~{int(total_samples * train_size)} samples")
    print(f"Validation: ~{int(total_samples * val_size)} samples")
    print(f"Test: ~{int(total_samples * test_size)} samples")
    
    # First split: separate test set
    train_val, test = train_test_split(
        df, 
        test_size=test_size, 
        random_state=random_state,
        stratify=df['Class']
    )
    
    # Second split: separate train and validation from remaining data
    relative_val_size = val_size / (train_size + val_size)
    train, val = train_test_split(
        train_val,
        test_size=relative_val_size,
        random_state=random_state,
        stratify=train_val['Class']
    )
    
    # Save the splits
    print("\nSaving split datasets...")
    train.to_csv(os.path.join(output_dir, 'train.tsv'), sep='\t', index=False)
    val.to_csv(os.path.join(output_dir, 'dev.tsv'), sep='\t', index=False)
    test.to_csv(os.path.join(output_dir, 'test.tsv'), sep='\t', index=False)
    
    # Print detailed statistics
    print(f"\nActual split sizes:")
    print(f"Training set: {len(train)} sequences ({len(train)/total_samples*100:.1f}%)")
    print(f"Validation set: {len(val)} sequences ({len(val)/total_samples*100:.1f}%)")
    print(f"Test set: {len(test)} sequences ({len(test)/total_samples*100:.1f}%)")
    
    print("\nClass distribution:")
    print("\nTraining set:")
    train_dist = train['Class'].value_counts()
    print(train_dist)
    print("\nValidation set:")
    val_dist = val['Class'].value_counts()
    print(val_dist)
    print("\nTest set:")
    test_dist = test['Class'].value_counts()
    print(test_dist)
    
    # Check if classes are balanced
    print("\nClass balance analysis:")
    train_ratio = train_dist.max() / train_dist.min()
    print(f"Training set imbalance ratio (max/min): {train_ratio:.2f}")
    if train_ratio > 3:
        print("Warning: Training set shows significant class imbalance. Consider using class weights or balancing techniques.")

if __name__ == '__main__':
    # Configuration
    input_file = 'concatenated_sequences_cleaned.txt'
    output_dir = 'split_data'
    
    # Split the data
    split_data(input_file, output_dir)