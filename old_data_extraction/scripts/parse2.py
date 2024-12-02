import pandas as pd
import requests
import csv
from tqdm import tqdm
import time
import asyncio
import aiohttp
from typing import List, Dict
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import warnings
warnings.filterwarnings('ignore')

class AsyncRateLimiter:
    """Asynchronous rate limiter for API requests"""
    def __init__(self, max_requests_per_second: int = 10):
        self.delay = 1.0 / max_requests_per_second
        self.last_call = time.time()
        self._lock = asyncio.Lock()
    
    async def acquire(self):
        async with self._lock:
            now = time.time()
            time_passed = now - self.last_call
            if time_passed < self.delay:
                await asyncio.sleep(self.delay - time_passed)
            self.last_call = time.time()

async def fetch_sequence(session: aiohttp.ClientSession, 
                        pdb_id: str, 
                        chain_id: str, 
                        rate_limiter: AsyncRateLimiter,
                        pbar: tqdm) -> Dict:
    """Fetch sequence for a single PDB entry asynchronously"""
    try:
        await rate_limiter.acquire()
        
        url = f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id.lower()}'
        async with session.get(url) as response:
            if response.status != 200:
                pbar.set_description(f"✗ {pdb_id}_{chain_id} (HTTP {response.status})")
                return None
                
            data = await response.json()
            
            # Process response
            molecules = data.get(pdb_id.lower(), [])
            for molecule in molecules:
                if (molecule.get('molecule_type') == 'polypeptide(L)' and 
                    'in_chains' in molecule and 
                    chain_id in molecule.get('in_chains', [])):
                    
                    sequence = molecule.get('sequence')
                    if sequence:
                        pbar.set_description(f"✓ {pdb_id}_{chain_id}")
                        return {
                            'pdb_id': pdb_id,
                            'chain_id': chain_id,
                            'sequence': sequence
                        }
            
            pbar.set_description(f"✗ {pdb_id}_{chain_id} (No sequence)")
            return None
            
    except Exception as e:
        pbar.set_description(f"✗ {pdb_id}_{chain_id} (Error)")
        return None

async def process_batch(entries: List[Dict], 
                       batch_size: int = 50,
                       max_concurrent: int = 50) -> List[Dict]:
    """Process a batch of entries concurrently"""
    
    rate_limiter = AsyncRateLimiter(max_requests_per_second=20)
    results = []
    
    # Configure connection pooling
    conn = aiohttp.TCPConnector(limit=max_concurrent, ttl_dns_cache=300)
    timeout = aiohttp.ClientTimeout(total=30)
    
    async with aiohttp.ClientSession(connector=conn, timeout=timeout) as session:
        with tqdm(total=len(entries), desc="Fetching sequences", unit="entries") as pbar:
            # Process in smaller chunks to maintain progress bar responsiveness
            for i in range(0, len(entries), batch_size):
                chunk = entries[i:i + batch_size]
                tasks = [
                    fetch_sequence(session, entry['pdb_id'], entry['chain_id'], 
                                 rate_limiter, pbar) 
                    for entry in chunk
                ]
                
                # Wait for all tasks in chunk to complete
                chunk_results = await asyncio.gather(*tasks)
                results.extend([r for r in chunk_results if r is not None])
                
                # Update progress
                pbar.update(len(chunk))
    
    return results

def save_results(output_file: str, results: List[Dict]):
    """Save results to CSV with proper encoding"""
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        if results:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)

async def main(input_file: str, output_file: str, batch_size: int = 50, max_concurrent: int = 50):
    """Main async function to process PDB entries"""
    
    # Read input file
    df = pd.read_csv(input_file, sep='\t', names=['uniprot_id', 'pdb_id', 'chain_id'])
    
    # Convert DataFrame to list of dictionaries
    entries = df[['pdb_id', 'chain_id']].to_dict('records')
    
    print(f"Processing {len(entries)} entries...")
    start_time = time.time()
    
    # Process entries
    results = await process_batch(entries, batch_size, max_concurrent)
    
    # Save results
    save_results(output_file, results)
    
    # Print summary
    elapsed_time = time.time() - start_time
    print("\n=== Summary ===")
    print(f"Successfully retrieved: {len(results)} sequences")
    print(f"Failed: {len(entries) - len(results)} entries")
    print(f"Success rate: {(len(results) / len(entries) * 100):.1f}%")
    print(f"Total time: {elapsed_time:.2f} seconds")
    print(f"Average time per entry: {(elapsed_time / len(entries)):.3f} seconds")
    print(f"\nResults saved to: {output_file}")

if __name__ == "__main__":
    input_file = "list_of_prots.txt"
    output_file = "pdb_sequences.csv"
    
    # Run with asyncio
    asyncio.run(main(
        input_file=input_file,
        output_file=output_file,
        batch_size=50,          # Adjust based on memory
        max_concurrent=50       # Adjust based on API limits
    ))