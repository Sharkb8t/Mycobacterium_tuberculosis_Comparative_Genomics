# Install datasets tool
conda install -c conda-forge ncbi-datasets-cli

# Download both genomes
datasets download genome taxon "Mycobacterium tuberculosis" --reference --include genome,gff3 --filename mtb_genomes.zip

# Extract
unzip mtb_genomes.zip -d data/raw/