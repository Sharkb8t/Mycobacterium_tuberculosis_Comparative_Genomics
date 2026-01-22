# Methods Documentation

## 1. Project Overview

This project performs comparative genomics analysis between two strains of Mycobacterium tuberculosis: H37Rv (reference strain) and CDC1551 (clinical isolate). The goal is to identify genomic differences that may contribute to phenotypic variations.

## 2. Data Acquisition

### 2.1 Genome Sources
- **H37Rv**: NCBI Assembly GCF_000195955.2 (ASM19595v2)
- **CDC1551**: NCBI Assembly GCF_000008585.1 (ASM858v1)

### 2.2 File Formats Downloaded
1. **FASTA (.fna)**: Complete nucleotide sequences
2. **GFF3 (.gff)**: Gene annotations and feature locations

### 2.3 Download Process
Automated download using `download_data.py` script with:
- HTTP streaming with progress bars
- File integrity validation
- Automatic decompression of .gz files
- Generation of data README

## 3. Genome Analysis Methods

### 3.1 Basic Statistics
- **Genome size**: Total length of all contigs
- **GC content**: Percentage of G+C nucleotides
- **Contiguity metrics**: N50, N75, L50 calculations
- **Assembly quality**: Contig length distribution analysis

### 3.2 Composition Analysis
- **Sliding window analysis**: GC content, GC skew, AT skew in 5kb windows
- **Dinucleotide frequencies**: Calculation of all 16 dinucleotide pairs
- **Sequence patterns**: Identification of compositional biases

### 3.3 Statistical Methods
All calculations performed using:
- **Biopython**: Sequence parsing and basic calculations
- **NumPy**: Numerical operations and statistics
- **Pandas**: Data organization and analysis

## 4. Comparative Genomics Methods

### 4.1 Alignment Strategy
Due to computational constraints on Windows, alignment is simulated but follows principles of:
- **Whole genome alignment**: Colinear mapping of sequences
- **Block-based approach**: Identification of conserved regions
- **Identity calculation**: Percent identity for each alignment block

### 4.2 Variant Analysis
- **SNP identification**: Single nucleotide polymorphisms
- **Indel detection**: Insertions and deletions
- **Variant classification**: Type, length, and position
- **Impact prediction**: Synonymous, missense, nonsense, intergenic

### 4.3 Quality Control
- **Variant filtering**: Quality score threshold (Q > 20)
- **Region masking**: Exclusion of low-complexity regions
- **Statistical validation**: Ti/Tv ratio calculation

## 5. Visualization Methods

### 5.1 Static Visualizations
- **Library**: Matplotlib with Seaborn styling
- **Resolution**: 300 DPI for publication quality
- **Color schemes**: Colorblind-friendly palettes
- **Figure organization**: Multi-panel figures for comprehensive views

### 5.2 Plot Types Generated
1. **Genome comparison bar charts**: Basic statistics
2. **GC content profiles**: Sliding window analysis
3. **Alignment dot plots**: Whole genome comparisons
4. **Variant distributions**: Type, impact, and density
5. **Nx curves**: Assembly continuity assessment

## 6. Computational Environment

### 6.1 Software Stack
- **Python 3.9**: Core programming language
- **Biopython 1.81**: Bioinformatics operations
- **NumPy 1.24+**: Numerical computing
- **Pandas 2.0+**: Data manipulation
- **Matplotlib 3.7+**: Visualization
- **Jupyter**: Interactive analysis

### 6.2 Environment Management
- **Conda**: Package and dependency management
- **environment.yml**: Reproducible environment specification
- **Virtual environment**: Isolation from system Python

## 7. Reproducibility

### 7.1 Version Control
- **Git**: All code and documentation version-controlled
- **GitHub**: Remote repository for collaboration
- **Commit messages**: Descriptive changes with rationale

### 7.2 Data Provenance
- **Source documentation**: All data sources clearly cited
- **Processing logs**: All transformations recorded
- **Parameter records**: Analysis parameters saved

## 8. Quality Assurance

### 8.1 Data Quality
- **Sequence validation**: Check for ambiguous bases
- **File integrity**: MD5 checksum verification
- **Size verification**: Confirm expected genome sizes

### 8.2 Analysis Quality
- **Method validation**: Compare with established tools
- **Parameter sensitivity**: Test impact of parameter changes
- **Result verification**: Manual inspection of key findings

## 9. Limitations

### 9.1 Computational Constraints
- **Alignment simulation**: Due to Windows environment limitations
- **Memory requirements**: Large genome files require sufficient RAM
- **Processing time**: Some analyses may take several minutes

### 9.2 Biological Interpretation
- **Simulated data**: Some analyses use simulated data for demonstration
- **Functional annotation**: Limited to basic impact prediction
- **Experimental validation**: Computational predictions require experimental confirmation

## 10. Future Directions

### 10.1 Technical Improvements
- Integration of real alignment tools (MUMmer, BLAST)
- Addition of protein sequence analysis
- Implementation of machine learning for variant impact

### 10.2 Biological Extensions
- Inclusion of additional Mtb strains
- Drug resistance gene analysis
- Virulence factor comparison
- Metabolic pathway reconstruction

## 11. References

See `references.bib` for complete bibliography of methods and tools used.