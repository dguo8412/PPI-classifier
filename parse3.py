import pandas as pd
import requests
import csv
from tqdm import tqdm
import time
import asyncio
import aiohttp
from typing import List, Dict, Tuple
import warnings
warnings.filterwarnings('ignore')

class AsyncRateLimiter:
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
                         pbar: tqdm) -> Tuple[Dict, Dict]:
    """
    Fetch sequence for a single PDB entry asynchronously
    Returns: Tuple of (success_result, failure_result)
    """
    try:
        await rate_limiter.acquire()
        
        # Convert chain_id to string and uppercase
        chain_id = str(chain_id).upper()
        
        url = f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id.lower()}'
        async with session.get(url) as response:
            if response.status != 200:
                error_info = {
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'error': f"HTTP {response.status}",
                    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
                }
                pbar.set_description(f"✗ {pdb_id}_{chain_id} (HTTP {response.status})")
                return None, error_info
                
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
                        }, None
            
            error_info = {
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'error': "No sequence found",
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
            }
            pbar.set_description(f"✗ {pdb_id}_{chain_id} (No sequence)")
            return None, error_info
            
    except Exception as e:
        error_info = {
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'error': str(e),
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        pbar.set_description(f"✗ {pdb_id}_{chain_id} (Error)")
        return None, error_info

async def process_batch(entries: List[Dict], 
                       batch_size: int = 50,
                       max_concurrent: int = 50) -> Tuple[List[Dict], List[Dict]]:
    """Process a batch of entries concurrently"""
    
    rate_limiter = AsyncRateLimiter(max_requests_per_second=20)
    successes = []
    failures = []
    
    conn = aiohttp.TCPConnector(limit=max_concurrent, ttl_dns_cache=300)
    timeout = aiohttp.ClientTimeout(total=30)
    
    async with aiohttp.ClientSession(connector=conn, timeout=timeout) as session:
        with tqdm(total=len(entries), desc="Fetching sequences", unit="entries") as pbar:
            for i in range(0, len(entries), batch_size):
                chunk = entries[i:i + batch_size]
                tasks = [
                    fetch_sequence(session, entry['pdb_id'], entry['chain_id'], 
                                 rate_limiter, pbar) 
                    for entry in chunk
                ]
                
                chunk_results = await asyncio.gather(*tasks)
                
                # Separate successes and failures
                for success, failure in chunk_results:
                    if success:
                        successes.append(success)
                    if failure:
                        failures.append(failure)
                
                pbar.update(len(chunk))
    
    return successes, failures

def save_results(successes: List[Dict], failures: List[Dict], 
                success_file: str, failure_file: str):
    """Save results to separate success and failure CSV files"""
    
    # Save successful results
    if successes:
        with open(success_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=successes[0].keys())
            writer.writeheader()
            writer.writerows(successes)
    
    # Save failed results
    if failures:
        with open(failure_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=failures[0].keys())
            writer.writeheader()
            writer.writerows(failures)

async def main(input_file: str, success_file: str, failure_file: str, 
               batch_size: int = 50, max_concurrent: int = 50):
    """Main async function to process PDB entries"""
    
    # Read input file
    df = pd.read_csv(input_file, sep='\t', names=['uniprot_id', 'pdb_id', 'chain_id'])
    entries = df[['pdb_id', 'chain_id']].to_dict('records')
    
    print(f"Processing {len(entries)} entries...")
    start_time = time.time()
    
    # Process entries
    successes, failures = await process_batch(entries, batch_size, max_concurrent)
    
    # Save results
    save_results(successes, failures, success_file, failure_file)
    
    # Print summary
    elapsed_time = time.time() - start_time
    print("\n=== Summary ===")
    print(f"Total entries processed: {len(entries)}")
    print(f"Successful retrievals: {len(successes)}")
    print(f"Failed retrievals: {len(failures)}")
    print(f"Success rate: {(len(successes) / len(entries) * 100):.1f}%")
    print(f"Total time: {elapsed_time:.2f} seconds")
    print(f"Average time per entry: {(elapsed_time / len(entries)):.3f} seconds")
    print(f"\nResults saved to:")
    print(f"  Successes: {success_file}")
    print(f"  Failures: {failure_file}")

if __name__ == "__main__":
    input_file = "list_of_prots.txt"
    success_file = "pdb_sequences_success.csv"
    failure_file = "pdb_sequences_failed.csv"
    
    # Run with asyncio
    asyncio.run(main(
        input_file=input_file,
        success_file=success_file,
        failure_file=failure_file,
        batch_size=100,
        max_concurrent=100
    ))