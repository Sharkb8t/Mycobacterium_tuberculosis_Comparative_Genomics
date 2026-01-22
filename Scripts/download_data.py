#!/usr/bin/env python
"""
Author: Dalton A. Schmidt
GitHub: https://github.com/Sharkb8t

Download and prepare Mycobacterium tuberculosis genome data for comparative genomics.
This script downloads FASTA and GFF files for H37Rv and CDC1551 strains from NCBI.
"""

import os
import gzip
import shutil
import requests
import argparse
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqUtils import GC
import warnings
warnings.filterwarnings('ignore')

class MTBDataDownloader:
    def __init__(self, project_root=None):
        if project_root is None:
            # Get the project root (parent directory of Scripts)
            self.project_root = Path(__file__).parent.parent
        else:
            self.project_root = Path(project_root)
        
        # Define paths
        self.raw_data_dir = self.project_root / "Data" / "raw"
        self.processed_data_dir = self.project_root / "Data" / "processed"
        
        # Create directories
        self.raw_data_dir.mkdir(parents=True, exist_ok=True)
        self.processed_data_dir.mkdir(parents=True, exist_ok=True)
        
        # Define genome URLs
        self.genome_urls = {
            "H37Rv": {
                "fasta": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz",
                "gff": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz",
                "genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gbff.gz"
            },
            "CDC1551": {
                "fasta": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/585/GCF_000008585.1_ASM858v1/GCF_000008585.1_ASM858v1_genomic.fna.gz",
                "gff": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/585/GCF_000008585.1_ASM858v1/GCF_000008585.1_ASM858v1_genomic.gff.gz",
                "genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/585/GCF_000008585.1_ASM858v1/GCF_000008585.1_ASM858v1_genomic.gbff.gz"
            }
        }
    
    def download_file(self, url, output_path):
        """Download a file with progress bar."""
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        
        with open(output_path, 'wb') as f:
            with tqdm(total=total_size, unit='B', unit_scale=True, 
                    desc=f"Downloading {output_path.name}") as pbar:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))
        
        return output_path
    
    def decompress_file(self, compressed_path, output_path):
        """Decompress .gz file."""
        with gzip.open(compressed_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"  Decompressed to: {output_path.name}")
        
        # Remove the compressed file
        compressed_path.unlink()
        return output_path
    
    def verify_genome_file(self, file_path):
        """Verify genome file integrity and get basic stats."""
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
            if len(records) > 0:
                total_length = sum(len(rec.seq) for rec in records)
                gc_content = sum(GC(rec.seq) * len(rec.seq) for rec in records) / total_length
                return {
                    "status": "valid",
                    "records": len(records),
                    "total_length": total_length,
                    "gc_content": round(gc_content, 2)
                }
            else:
                return {"status": "empty"}
        except Exception as e:
            return {"status": f"error: {str(e)}"}
    
    def download_all_genomes(self):
        """Download and process all genome files."""
        print("=" * 70)
        print("Mycobacterium tuberculosis Comparative Genomics - Data Download")
        print("=" * 70)
        
        results = {}
        
        for strain, urls in self.genome_urls.items():
            print(f"\n{'='*40}")
            print(f"Processing: {strain}")
            print(f"{'='*40}")
            
            strain_results = {}
            
            for file_type, url in urls.items():
                # Define file paths
                compressed_name = f"{strain}_{file_type}.gz"
                final_name = f"{strain}_{file_type}"
                
                compressed_path = self.raw_data_dir / compressed_name
                final_path = self.raw_data_dir / final_name
                
                # Skip if already exists
                if final_path.exists():
                    print(f"✓ {file_type.upper()} already exists: {final_name}")
                    strain_results[file_type] = "already_exists"
                    continue
                
                try:
                    # Download
                    print(f"  Downloading {file_type}...")
                    self.download_file(url, compressed_path)
                    
                    # Decompress
                    self.decompress_file(compressed_path, final_path)
                    
                    # Verify FASTA files
                    if file_type == "fasta":
                        stats = self.verify_genome_file(final_path)
                        print(f"  Genome stats: {stats['records']} contig(s), "
                            f"{stats['total_length']:,} bp, GC: {stats['gc_content']}%")
                    
                    strain_results[file_type] = "success"
                    
                except Exception as e:
                    print(f"  ✗ Error downloading {file_type}: {e}")
                    strain_results[file_type] = f"error: {e}"
            
            results[strain] = strain_results
        
        # Print summary
        self.print_summary(results)
        
        # Generate a README for the data
        self.generate_data_readme()
        
        return results
    
    def print_summary(self, results):
        """Print download summary."""
        print("\n" + "=" * 70)
        print("DOWNLOAD SUMMARY")
        print("=" * 70)
        
        for strain, files in results.items():
            print(f"\n{strain}:")
            for file_type, status in files.items():
                if "success" in str(status) or "already_exists" in str(status):
                    print(f"  ✓ {file_type}: OK")
                else:
                    print(f"  ✗ {file_type}: {status}")
    
    def generate_data_readme(self):
        """Generate a README file for the data directory."""
        readme_content = """# Mycobacterium tuberculosis Genome Data

## Files Description

### H37Rv (Reference Strain)
- **H37Rv_fasta**: Complete genome sequence in FASTA format
- **H37Rv_gff**: Gene annotations in GFF3 format
- **H37Rv_genbank**: Annotated genome in GenBank format

### CDC1551 (Clinical Strain)
- **CDC1551_fasta**: Complete genome sequence in FASTA format
- **CDC1551_gff**: Gene annotations in GFF3 format
- **CDC1551_genbank**: Annotated genome in GenBank format

## Data Sources

All data downloaded from NCBI RefSeq:
- H37Rv: GCF_000195955.2 (ASM19595v2)
- CDC1551: GCF_000008585.1 (ASM858v1)

## File Formats

1. **FASTA (.fna)**: Contains nucleotide sequences
2. **GFF3 (.gff)**: Contains gene, CDS, and feature annotations
3. **GenBank (.gbff)**: Comprehensive annotation format

## Usage in Analysis

These files serve as the basis for comparative genomics analysis including:
- Genome alignment
- Gene content comparison
- Variant calling
- Functional annotation comparison

---
*Generated by download_data.py*
"""
        
        readme_path = self.raw_data_dir / "README.md"
        with open(readme_path, 'w') as f:
            f.write(readme_content)
        
        print(f"\n✓ Data README generated: {readme_path}")

def main():
    parser = argparse.ArgumentParser(description="Download M. tuberculosis genome data")
    parser.add_argument('--project-root', type=str, default=None,
                    help='Path to project root directory')
    
    args = parser.parse_args()
    
    downloader = MTBDataDownloader(args.project_root)
    downloader.download_all_genomes()

if __name__ == "__main__":
    main()