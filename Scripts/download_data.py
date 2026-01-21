#!/usr/bin/env python
"""
Download and prepare M. tuberculosis genomes for comparative genomics.
Run this once to set up your project data.
"""

import os
import gzip
import shutil
import requests
from pathlib import Path
from tqdm import tqdm  # pip install tqdm

def download_file(url, output_path):
    """Download a file with progress bar."""
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    
    with open(output_path, 'wb') as f:
        with tqdm(total=total_size, unit='B', unit_scale=True, desc=output_path.name) as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                pbar.update(len(chunk))

def decompress_file(compressed_path, output_path):
    """Decompress .gz file."""
    with gzip.open(compressed_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Decompressed: {output_path}")

def main():
    # Project paths
    project_root = Path(r"C:\Users\dalto\Desktop\Education Resources\Personal\Projects\Mycobacterium_tuberculosis_Comparative_Genomics")
    raw_dir = project_root / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    
    # URLs for the genomes (updated for direct download)
    genome_data = {
        "H37Rv": {
            "fasta": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz",
            "gff": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz"
        },
        "CDC1551": {
            "fasta": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/585/GCF_000008585.1_ASM858v1/GCF_000008585.1_ASM858v1_genomic.fna.gz",
            "gff": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/585/GCF_000008585.1_ASM858v1/GCF_000008585.1_ASM858v1_genomic.gff.gz"
        }
    }
    
    print("Downloading and processing genome data...")
    print("=" * 60)
    
    for strain, files in genome_data.items():
        print(f"\nProcessing {strain}:")
        print("-" * 40)
        
        for file_type, url in files.items():
            # Create file names
            compressed_file = raw_dir / f"{strain}_{file_type}.gz"
            final_file = raw_dir / f"{strain}_{file_type}.{file_type}"
            
            # Download compressed file
            print(f"Downloading {file_type}...")
            download_file(url, compressed_file)
            
            # Decompress
            print(f"Decompressing...")
            decompress_file(compressed_file, final_file)
            
            # Remove compressed file
            compressed_file.unlink()
            print(f"Removed compressed file: {compressed_file.name}")
            
            # Verify file
            if final_file.exists():
                size_mb = final_file.stat().st_size / (1024 * 1024)
                print(f"✓ Created: {final_file.name} ({size_mb:.1f} MB)")
    
    print("\n" + "=" * 60)
    print("Data setup complete!")
    print(f"Files available in: {raw_dir}")
    
    # List created files
    print("\nCreated files:")
    for f in sorted(raw_dir.glob("*")):
        if not f.name.endswith('.gz'):  # Only show uncompressed files
            size_mb = f.stat().st_size / (1024 * 1024)
            print(f"  - {f.name} ({size_mb:.1f} MB)")

if __name__ == "__main__":
    main()