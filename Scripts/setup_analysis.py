#!/usr/bin/env python
"""
Setup script for M. tuberculosis comparative genomics
Author: [Your Name]
GitHub: https://github.com/Sharkb8t
"""

import os
import subprocess
import requests
import gzip
import shutil
from pathlib import Path

class MTBComparativeGenomics:
    def __init__(self, project_dir):
        self.project_dir = Path(project_dir)
        self.setup_directories()
        
    def setup_directories(self):
        """Create project directory structure"""
        dirs = [
            'data/raw',
            'data/processed',
            'scripts',
            'notebooks',
            'results/alignments',
            'results/annotations',
            'results/variants',
            'results/plots',
            'docs'
        ]
        
        for d in dirs:
            (self.project_dir / d).mkdir(parents=True, exist_ok=True)
        print("Directory structure created.")
    
    def download_genome(self, url, output_name):
        """Download and extract genome files"""
        output_path = self.project_dir / 'data/raw' / output_name
        
        print(f"Downloading {url}...")
        response = requests.get(url, stream=True)
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        # Extract if compressed
        if output_path.suffix == '.gz':
            with gzip.open(output_path, 'rb') as f_in:
                with open(output_path.with_suffix(''), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            output_path.unlink()  # Remove gz file
        
        print(f"Downloaded: {output_name}")
        return output_path.with_suffix('') if output_path.suffix == '.gz' else output_path

def main():
    project = MTBComparativeGenomics(r"C:\Users\dalto\Desktop\Education Resources\Personal\Projects\Mycobacterium_tuberculosis_Comparative_Genomics")
    
    # Download genomes (uncomment to use)
    '''
    urls = {
        'H37Rv': {
            'genome': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz',
            'annotation': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz'
        },
        'CDC1551': {
            'genome': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/585/GCF_000008585.1_ASM858v1/GCF_000008585.1_ASM858v1_genomic.fna.gz',
            'annotation': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/585/GCF_000008585.1_ASM858v1/GCF_000008585.1_ASM858v1_genomic.gff.gz'
        }
    }
    
    for strain, files in urls.items():
        for file_type, url in files.items():
            output_name = f"{strain}_{file_type}{Path(url).suffix}"
            project.download_genome(url, output_name)
    '''
    
    print("\nSetup complete!")
    print("Next steps:")
    print("1. Download genome files using the URLs provided")
    print("2. Place them in data/raw/ directory")
    print("3. Begin analysis in Jupyter Notebook")

if __name__ == "__main__":
    main()