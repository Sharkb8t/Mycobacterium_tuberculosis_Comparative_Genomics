# First annotate with Prokka
prokka --outdir annotations/H37Rv --prefix H37Rv H37Rv.fna
prokka --outdir annotations/CDC1551 --prefix CDC1551 CDC1551.fna

# Run Roary for pan-genome
roary -f results/pan_genome -e -n -v annotations/*/*.gff