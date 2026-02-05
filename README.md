# Mycobacterium tuberculosis Comparative Genomics

[![Python 3.14+](https://img.shields.io/badge/python-3.14+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue)](https://github.com/Sharkb8t/Mycobacterium_tuberculosis_Comparative_Genomics)

A comprehensive bioinformatics pipeline for comparative genomics analysis of Mycobacterium tuberculosis strains H37Rv and CDC1551. This project integrates genome sequence analysis, gene annotation parsing, and variant calling to study genomic diversity and evolution in tuberculosis pathogens.

## 📋 Features

- **Complete Data Pipeline**: Download, process, and analyze Mtb genome data from NCBI
- **Multi-format Support**: Handles FASTA, GBFF, and GFF3 file formats
- **Integrated Analysis**: Combines genome statistics, gene annotations, and comparative genomics
- **Variant Analysis**: Simulates and analyzes SNPs, insertions, deletions, and MNPs
- **Gene Annotation**: Parses GFF3 and GBFF files for functional insights
- **Visualization**: Generates comprehensive plots and reports
- **Reproducible**: Conda environment and requirements specification

## 🏗️ Project Structure
---
Mycobacterium_tuberculosis_Comparative_Genomics
│
├── README.md # This file
├── LICENSE
├── .gitignore
│
├── Data/
│ ├── Raw/ # Original genome files (FASTA, GBFF, & GFF3 formats)
│ │ ├── CDC1551_fasta.fna
│ │ ├── CDC1551_gbff.gbff
│ │ ├── CDC1551_gff.gff
│ │ ├── H37Rv_fasta.fna
│ │ ├── H37Rv_gbff.gbff
│ │ └── H37Rv_gff.gff
│ │
│ └── Processed/ # Cleaned/processed data
│
├── Notebooks/ # Jupyter notebooks
│ ├── 01_Data_Acquisition.ipynb
│ ├── 02_Genome_Statistics.ipynb
│ └── 03_Comparative_Analysis.ipynb
│
├── Scripts/ # Python scripts
│ ├── download_data.py # Download genome data from NCBI
│ ├── genome_analyzer.py # Genome statistics and analysis
│ ├── annotation_parser.py # Parse GFF3 and GBFF annotations
│ ├── alignment_analysis.py # Comparative alignment and variant analysis
│ └── integrated_analysis.py # Complete analysis pipeline
│
├── Results/ # Analysis outputs (generated during analysis)
│ ├── Alignments/ # Alignment coordinates and statistics
│ ├── Annotations/ # Parsed gene annotations
│ ├── Variants/ # Variant calls and statistics
│ ├── Plots/ # Visualization files
│ └── Reports/ # Comprehensive analysis reports
│
├── Documentation/
│ ├── methods.md # Detailed methodology
│ └── references.bib # Bibliography
│
├── environment.yml # Conda environment specification
└── requirements.txt # Python package requirements
---

## 🚀 Quick Start

### 1. Prerequisites

- Python 3.8 or higher
- Conda (recommended) or pip

### 2. Installation

**Using Conda (recommended):**
```bash
# Clone the repository
git clone https://github.com/Sharkb8t/Mycobacterium_tuberculosis_Comparative_Genomics.git
cd Mycobacterium_tuberculosis_Comparative_Genomics

# Create and activate conda environment
conda env create -f environment.yml
conda activate mtb_comp_genomics
```

**Using pip**
```bash
# Clone the repository
git clone https://github.com/Sharkb8t/Mycobacterium_tuberculosis_Comparative_Genomics.git
cd Mycobacterium_tuberculosis_Comparative_Genomics

# Install dependencies
pip install -r requirements.txt
```

### 3. Download Data
```bash
# Run from the project root directory
python Scripts/download_data.py
```
This will download the following files if not downloaded from the project directory:

- H37Rv and CDC1551 FASTA files
- GFF3 annotation files
- GBFF (GenBank) annotation files

### 4. Run Complete Analysis Pipeline
```bash
# Run integrated analysis (recommended for first-time users)
python Scripts/integrated_analysis.py
```

## 📊 Analysis Pipeline

### Phase 1: Data Acquisition & Genome Statistics
```bash
# Download genome data
python Scripts/download_data.py

# Analyze genome characteristics
python Scripts/genome_analyzer.py --strains H37Rv CDC1551 --save-csv
```
**Outputs:**

- Genome length, GC content, contig statistics
- N50, N75, and L50 calculations
- Comparative genome statistics

### Phase 2: Gene Annotation Analysis
```bash
# Parse and analyze gene annotations
python Scripts/annotation_parser.py --strains H37Rv CDC1551 --parse-all --compare
```
**Outputs**

- Parsed GFF3 and GBFF annotations
- Gene statistics (count, length, density)
- Gene content comparison between strains 

### Phase 3: Comparative Genomics
```bash
# Run comparative alignment and variant analysis
python Scripts/alignment_analysis.py
```
**Outputs**

- Simulated whole genome alignments
- Variant calls (SNPs, indels, MNPs)
- Variant-gene associations
- Comprehensive visualizations

### Phase 4: Integrated Analysis
```bash
# Run complete pipeline with single command
python Scripts/integrated_analysis.py --strains H37Rv CDC1551
```
**Outputs**

- Master analysis report (`Results/integrated_analysis_report.md`)
- All intermediate results
- Publication-ready figures

## 🔬 Key Analyses Performed

### 1. Genome Statistics

- Total genome length and GC content
- Contig/assembly statistics (N50, N75, L50)
- Genome size comparison between strains
- GC skew and nucleotide composition analysis

### 2. Gene Annotation

- Gene count and density calculations
- Functional annotation parsing
- Pseudogene identification
- Gene content conservation analysis

### 3. Comparative Genomics

- Whole genome alignment simulation
- Variant calling and characterization
- Ti/Tv (transition/transversion) ratio calculation
- Variant impact prediction (synonymous/missense/nonsense)

### 4. Integration

- Variant annotation with gene context
- Identification of affected biological pathways
- Strain-specific gene analysis

## 📈 Sample Results

### Genome Comparison Statistics


### Variant Analysis

- **Total Variants**: ~1,500-2,000 between strains
- **SNPs**: 75% of total variants
- **Ti/Tv Ratio**: ~2.0-3.0 (typical for bacterial genomes)
- **Coding Variants**: ~90% in gene regions
- **Synonymous**: ~75% of coding variants

## 📖 Jupyter Notebooks

The project includes three interactive notebooks for step-by-step analysis:
1. `01_Data_Acquisition.ipynb` - Download and validate genome data
2. `02_Genome_Statistics.ipynb` - Calculate and visualize genome statistics
3. `03_Comparative_Analysis.ipynb` - Perform comparative genomics analysis

To use the notebooks, be sure to launch Jupyter Notebook through the installed desktop application or via bash:
```bash
# Launch Jupyter
jupyter notebook
# Navigate to the Notebooks/ directory
```

## 🛠️ Script Documentation

### download_data.py
Downloads M. tuberculosis genome data for the CDC1551 & H37RV strains, from NCBI RefSeq.
```bash
python Scripts/download_data.py [--project-root PATH]
```

### genome_analyzer.py
Analyzes genomes sequences for bacteria strains in the `Raw/` directory.
```bash
python Scripts/genome_analyzer.py [--project-root PATH] [--strains H37Rv CDC1551] [--save-csv]
```

### annotation_parser.py
Parses GFF3 and GBFF annotation files for data analysis.
```bash
python Scripts/annotation_parser.py [--project-root PATH] [--strains H37Rv CDC1551] [--parse-all] [--compare]
```

### alignment_analysis.py
Performs comparative alignment and variant analysis on parsed genome data.
```bash
python Scripts/alignment_analysis.py [--project-dir PATH]
```

### integrated_analysis.py
This file will completely run the analysis pipeline in one command.
```bash
python Scripts/integrated_analysis.py [--project-dir PATH] [--strains H37Rv CDC1551]
```

## 🧪 Dependencies

Core dependencies (see `requirements.txt` for complete list):
- **Biopython** (≥1.86): Biological sequence analysis
- **NumPy** (≥2.4.1): Numerical computing
- **pandas** (≥2.3.3): Data manipulation
- **matplotlib** (≥3.10.8) & seaborn (≥0.13.2): Visualization
- **scikit-learn** (≥1.8.0): Machine learning utilities
- **Jupyter** (≥1.1.1): Interactive notebooks
- **gffutils** (≥0.13.0): GFF3 file parsing

## 📝 Output Files
Mycobacterium_tuberculosis_Comparative_Genomics/
│
└── Results/              # Analysis outputs
   ├── Alignments/
   │   ├── alignment_statistics.csv
   │   └── simulated_alignment_coordinates.csv
   │
   ├── Annotations/
   │   ├── all_strains_basic_stats.csv
   │   ├── CDC1551_basic_stats
   │   ├── CDC1551_gbff_annotations.csv
   │   ├── CDC1551_gff_annotations.csv
   │   ├── gene_comparison_H37Rv_vs_CDC1551.csv
   │   ├── H37Rv_basic_stats.csv
   │   ├── H37Rv_gbff_annotations.csv
   │   ├── H37Rv_gff_annotations.csv
   │   └── genome_analysis_report.txt
   │
   ├── Variants/
   │   ├── annotated_variants.csv
   │   ├── annotation_statistics.csv
   │   ├── simulated_variants.csv
   │   └── variant_statistics.csv
   │
   ├── Plots/
   │   ├── comparative_analysis_plots.png
   │   ├── gene_variant_analysis.png
   │   ├── genome_statistics_comparison.png
   │   ├── basic_genome_comparison.png
   │   └── gc_profiles.png
   │
   ├── enhanced_comparative_genomics_report.md
   ├── integrated_analysis_report.md
   ├── data_ascuisition_summary.csv
   ├── cdc1551_sliding_window.csv
   ├── contig_details.csv
   ├── genome_comparison_stats.csv
   ├── h37rv_sliding_window.csv
   └── custom_analysis_summary.csv

## 🎯 Biological Significance

This project addresses key questions in M. tuberculosis genomics:

- Genomic conservation between reference and clinical strains
- Functional implications of strain-specific variants
- Evolutionary patterns in tuberculosis pathogens
- Potential targets for strain-specific interventions

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Development Guidelines
- Follow PEP 8 style guide for Python code
- Add docstrings to all functions and classes
- Include tests for new functionality
- Update documentation as needed

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

- **NCBI** for providing M. tuberculosis genome data
- **Biopython** community for excellent bioinformatics tools
- **Open source community** for numerous scientific Python packages

## 📧 Contact
