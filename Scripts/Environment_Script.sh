# Create environment with bioinformatics tools
conda create -n mtb_comparative python=3.9 biopython jupyter notebook pandas numpy matplotlib seaborn scipy scikit-learn

# Activate environment
conda activate mtb_comparative

# Install bioinformatics packages
conda install -c bioconda samtools bcftools bedtools
conda install -c bioconda prokka mummer blast
conda install -c bioconda roary fasttree
conda install -c bioconda snpeff snpsift
conda install -c bioconda pygenomeviz