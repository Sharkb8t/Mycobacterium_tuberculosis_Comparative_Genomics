# Create BWA index
bwa index H37Rv.fna
# Align CDC1551 to H37Rv
bwa mem H37Rv.fna CDC1551.fna > alignment.sam
# Call variants
bcftools mpileup -f H37Rv.fna alignment.bam | bcftools call -mv -Ov -o variants.vcf