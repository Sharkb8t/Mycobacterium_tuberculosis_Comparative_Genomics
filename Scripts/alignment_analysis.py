"""
Author: Dalton A. Schmidt
GitHub: https://github.com/Sharkb8t

Comparative genomics alignment analysis for Mycobacterium tuberculosis strains.
Performs whole genome alignment simulation and variant calling between H37Rv and CDC1551.
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction as GC
import warnings
warnings.filterwarnings('ignore')

class ComparativeAlignment:
    def __init__(self, project_root=None):
        if project_root is None:
            # Try to auto-detect project root
            current_dir = Path.cwd()
            if current_dir.name == "Scripts":
                self.project_root = current_dir.parent
            else:
                self.project_root = current_dir
        else:
            self.project_root = Path(project_root)
    
        # Resolve to absolute path
        self.project_root = self.project_root.resolve()
        
        # Define paths
        self.raw_dir = self.project_root / "Data" / "Raw"
        self.results_dir = self.project_root / "Results"
        self.alignments_dir = self.results_dir / "Alignments"
        self.variants_dir = self.results_dir / "Variants"
        self.plots_dir = self.results_dir / "Plots"
        self.annotations_dir = self.results_dir / "Annotations"
        
        # Create directories
        self.alignments_dir.mkdir(parents=True, exist_ok=True)
        self.variants_dir.mkdir(parents=True, exist_ok=True)
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.annotations_dir.mkdir(parents=True, exist_ok=True)
        
        # Load genome lengths for simulation
        self.genome_lengths = self.get_genome_lengths()
    
    def get_genome_lengths(self):
        """Get actual genome lengths from FASTA files."""
        lengths = {}
        strains = ["H37Rv", "CDC1551"]
        
        for strain in strains:
            fasta_file = self.raw_dir / f"{strain}_fasta.fna"
            if not fasta_file.exists():
                # Try alternative naming
                fasta_file = self.raw_dir / f"{strain}_fasta"
            
            if fasta_file.exists():
                records = list(SeqIO.parse(fasta_file, "fasta"))
                lengths[strain] = sum(len(rec.seq) for rec in records)
            else:
                print(f"Warning: FASTA file not found for {strain}")
                # Use approximate lengths if files not found
                lengths = {"H37Rv": 4411532, "CDC1551": 4403837}
                break
        
        return lengths
    
    def simulate_alignment(self):
        """Simulate whole genome alignment for demonstration."""
        print("Simulating whole genome alignment...")
        print("Note: In a real analysis, you would use MUMmer4 or similar tools")
        
        np.random.seed(42)  # For reproducibility
        
        # Use actual genome lengths
        h37rv_length = self.genome_lengths.get("H37Rv", 4411532)
        
        # Simulate alignment blocks
        num_blocks = 150
        data = []
        
        # Start positions
        ref_pos = 1
        query_pos = 1
        
        for i in range(num_blocks):
            # Generate block with realistic properties
            block_length = int(np.random.exponential(30000))
            block_length = max(5000, min(block_length, 100000))
            
            # Most blocks are colinear (95%), some have rearrangements (5%)
            if np.random.random() < 0.95:
                # Colinear block with small variations
                query_start = query_pos + np.random.randint(-500, 500)
            else:
                # Rearranged block
                query_start = np.random.randint(1, h37rv_length - block_length)
            
            # Calculate identity (mostly high with some variation)
            identity = np.random.normal(99.7, 0.15)
            identity = max(98.0, min(100.0, identity))
            
            # Create block data
            block_data = {
                'block_id': i + 1,
                'ref_start': ref_pos,
                'ref_end': ref_pos + block_length,
                'query_start': query_start,
                'query_end': query_start + block_length,
                'ref_length': block_length,
                'query_length': block_length,
                'identity': identity,
                'ref_id': 'H37Rv',
                'query_id': 'CDC1551'
            }
            
            data.append(block_data)
            
            # Update positions for next block
            ref_pos += block_length + np.random.randint(0, 2000)
            query_pos = query_start + block_length + np.random.randint(0, 2000)
            
            # Stop if we exceed genome length
            if ref_pos > h37rv_length * 0.95:
                break
        
        df = pd.DataFrame(data)
        
        # Calculate coverage percentages
        total_aligned = df['ref_length'].sum()
        df['ref_coverage'] = (df['ref_length'] / h37rv_length) * 100
        df['query_coverage'] = (df['ref_length'] / self.genome_lengths.get("CDC1551", 4403837)) * 100
        
        # Save to file
        output_file = self.alignments_dir / "simulated_alignment_coordinates.csv"
        df.to_csv(output_file, index=False)
        
        print(f"✓ Simulated alignment saved to: {output_file}")
        print(f"  Total aligned bases: {total_aligned:,}")
        print(f"  Coverage: {(total_aligned / h37rv_length * 100):.2f}%")
        print(f"  Mean identity: {df['identity'].mean():.2f}%")
        
        return df
    
    def analyze_alignment(self, alignment_df):
        """Analyze alignment data and generate statistics."""
        print("\nAnalyzing alignment statistics...")
        
        if alignment_df.empty:
            print("  No alignment data to analyze")
            return {}
        
        h37rv_length = self.genome_lengths.get("H37Rv", 4411532)
        cdc1551_length = self.genome_lengths.get("CDC1551", 4403837)
        
        total_aligned = alignment_df['ref_length'].sum()
        
        stats = {
            'total_aligned_bp': total_aligned,
            'h37rv_coverage': (total_aligned / h37rv_length) * 100,
            'cdc1551_coverage': (total_aligned / cdc1551_length) * 100,
            'num_blocks': len(alignment_df),
            'mean_identity': alignment_df['identity'].mean(),
            'median_identity': alignment_df['identity'].median(),
            'mean_block_size': alignment_df['ref_length'].mean(),
            'max_block_size': alignment_df['ref_length'].max(),
            'min_block_size': alignment_df['ref_length'].min()
        }
        
        # Print statistics
        print(f"  Alignment blocks: {stats['num_blocks']}")
        print(f"  H37Rv coverage: {stats['h37rv_coverage']:.2f}%")
        print(f"  CDC1551 coverage: {stats['cdc1551_coverage']:.2f}%")
        print(f"  Mean identity: {stats['mean_identity']:.2f}%")
        print(f"  Mean block size: {stats['mean_block_size']:.0f} bp")
        
        # Save statistics
        stats_df = pd.DataFrame([stats])
        stats_df.to_csv(self.alignments_dir / "alignment_statistics.csv", index=False)
        
        return stats
    
    def simulate_variants(self):
        """Simulate genomic variants between strains."""
        print("\nSimulating genomic variants...")
        
        np.random.seed(42)
        
        # Use actual genome length
        genome_length = self.genome_lengths.get("H37Rv", 4411532)
        
        # Simulate realistic number of variants for Mtb strains
        # Typically ~1000-2000 SNPs between strains
        num_variants = np.random.randint(1200, 1800)
        
        # Generate random positions
        positions = np.random.choice(range(1, genome_length), num_variants, replace=False)
        positions.sort()
        
        data = []
        
        for i, pos in enumerate(positions):
            # Determine variant type
            variant_types = ['SNP', 'Insertion', 'Deletion', 'MNP']
            vtype = np.random.choice(variant_types, p=[0.75, 0.15, 0.09, 0.01])
            
            # Define nucleotide bases
            bases = ['A', 'C', 'G', 'T']
            
            # Generate variant details
            if vtype == 'SNP':
                ref_base = np.random.choice(bases)
                alt_base = np.random.choice([b for b in bases if b != ref_base])
                length = 1
                
            elif vtype == 'Insertion':
                ref_base = '-'
                ins_length = np.random.randint(1, 5)
                alt_base = ''.join(np.random.choice(bases, ins_length))
                length = ins_length
                
            elif vtype == 'Deletion':
                del_length = np.random.randint(1, 5)
                ref_base = ''.join(np.random.choice(bases, del_length))
                alt_base = '-'
                length = del_length
                
            else:  # MNP
                mnp_length = np.random.randint(2, 4)
                ref_base = ''.join(np.random.choice(bases, mnp_length))
                alt_base = ''.join([np.random.choice([b for b in bases if b != r]) for r in ref_base])
                length = mnp_length
            
            # Determine if in coding region (simplified)
            # Mtb genomes are ~90% coding
            in_gene = np.random.random() < 0.9
            
            # Determine impact
            if in_gene:
                impact = np.random.choice(['Synonymous', 'Missense', 'Nonsense'], 
                                        p=[0.75, 0.23, 0.02])
            else:
                impact = 'Intergenic'
            
            # Generate variant ID
            variant_id = f"VAR_{i+1:05d}"
            
            data.append({
                'position': pos,
                'variant_id': variant_id,
                'type': vtype,
                'ref': ref_base,
                'alt': alt_base,
                'length': length,
                'in_gene': in_gene,
                'impact': impact,
                'quality': np.random.uniform(30, 100),
                'filter': 'PASS'
            })
        
        df = pd.DataFrame(data)
        
        # Save to file
        output_file = self.variants_dir / "simulated_variants.csv"
        df.to_csv(output_file, index=False)
        
        print(f"✓ Simulated variants saved to: {output_file}")
        print(f"  Total variants: {len(df):,}")
        print(f"  SNPs: {(df['type'] == 'SNP').sum():,}")
        print(f"  Insertions: {(df['type'] == 'Insertion').sum():,}")
        print(f"  Deletions: {(df['type'] == 'Deletion').sum():,}")
        
        return df
    
    def analyze_variants(self, variants_df):
        """Analyze variant data and generate statistics."""
        print("\nAnalyzing variant statistics...")
        
        if variants_df.empty:
            print("  No variant data to analyze")
            return {}
        
        # Calculate basic statistics
        stats = {
            'total_variants': len(variants_df),
            'snps': (variants_df['type'] == 'SNP').sum(),
            'insertions': (variants_df['type'] == 'Insertion').sum(),
            'deletions': (variants_df['type'] == 'Deletion').sum(),
            'mnps': (variants_df['type'] == 'MNP').sum(),
            'in_genes': variants_df['in_gene'].sum(),
            'synonymous': (variants_df['impact'] == 'Synonymous').sum(),
            'missense': (variants_df['impact'] == 'Missense').sum(),
            'nonsense': (variants_df['impact'] == 'Nonsense').sum(),
            'intergenic': (variants_df['impact'] == 'Intergenic').sum()
        }
        
        # Calculate percentages
        stats['snps_percent'] = (stats['snps'] / stats['total_variants']) * 100
        stats['in_genes_percent'] = (stats['in_genes'] / stats['total_variants']) * 100
        
        # Calculate Ti/Tv ratio (simplified)
        snps = variants_df[variants_df['type'] == 'SNP']
        if len(snps) > 0:
            transitions = ['AG', 'GA', 'CT', 'TC']
            transversions = ['AC', 'CA', 'AT', 'TA', 'CG', 'GC', 'GT', 'TG']
            
            ti_count = 0
            tv_count = 0
            
            for _, snp in snps.iterrows():
                change = snp['ref'] + snp['alt']
                if change in transitions:
                    ti_count += 1
                elif change in transversions:
                    tv_count += 1
            
            if tv_count > 0:
                stats['ti_tv_ratio'] = ti_count / tv_count
            else:
                stats['ti_tv_ratio'] = float('inf')
        
        # Print statistics
        print(f"  Total variants: {stats['total_variants']:,}")
        print(f"  SNPs: {stats['snps']:,} ({stats['snps_percent']:.1f}%)")
        print(f"  In coding regions: {stats['in_genes']:,} ({stats['in_genes_percent']:.1f}%)")
        print(f"  Synonymous: {stats['synonymous']:,}")
        print(f"  Missense: {stats['missense']:,}")
        print(f"  Nonsense: {stats['nonsense']:,}")
        if 'ti_tv_ratio' in stats:
            print(f"  Ti/Tv ratio: {stats['ti_tv_ratio']:.2f}")
        
        # Save statistics
        stats_df = pd.DataFrame([stats])
        stats_df.to_csv(self.variants_dir / "variant_statistics.csv", index=False)
        
        return stats
    
    def create_visualizations(self, alignment_df, variants_df):
        """Create visualizations for alignment and variant analysis."""
        print("\nCreating visualizations...")
        
        # Create plots directory
        plots_dir = self.results_dir / "Plots"
        plots_dir.mkdir(exist_ok=True)
        
        # Set plotting style
        plt.style.use('seaborn-v0_8-whitegrid')
        sns.set_palette("husl")
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # Plot 1: Alignment dot plot
        if not alignment_df.empty:
            scatter = axes[0, 0].scatter(
                alignment_df['ref_start'] / 1e6,
                alignment_df['query_start'] / 1e6,
                c=alignment_df['identity'],
                cmap='viridis',
                s=15,
                alpha=0.6
            )
            axes[0, 0].set_xlabel('H37Rv Position (Mbp)')
            axes[0, 0].set_ylabel('CDC1551 Position (Mbp)')
            axes[0, 0].set_title('Whole Genome Alignment (Simulated)')
            axes[0, 0].grid(True, alpha=0.3)
            plt.colorbar(scatter, ax=axes[0, 0], label='Identity (%)')
        
        # Plot 2: Identity distribution
        if not alignment_df.empty:
            axes[0, 1].hist(alignment_df['identity'], bins=30,
                        color='skyblue', edgecolor='black')
            axes[0, 1].axvline(alignment_df['identity'].mean(),
                            color='red', linestyle='--',
                            label=f'Mean: {alignment_df["identity"].mean():.2f}%')
            axes[0, 1].set_xlabel('Alignment Identity (%)')
            axes[0, 1].set_ylabel('Frequency')
            axes[0, 1].set_title('Alignment Identity Distribution')
            axes[0, 1].legend()
            axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Block size distribution
        if not alignment_df.empty:
            axes[0, 2].hist(alignment_df['ref_length'] / 1000, bins=30,
                        color='lightgreen', edgecolor='black')
            axes[0, 2].axvline(alignment_df['ref_length'].mean() / 1000,
                            color='red', linestyle='--',
                            label=f'Mean: {alignment_df["ref_length"].mean()/1000:.1f} kbp')
            axes[0, 2].set_xlabel('Alignment Block Size (kbp)')
            axes[0, 2].set_ylabel('Frequency')
            axes[0, 2].set_title('Block Size Distribution')
            axes[0, 2].legend()
            axes[0, 2].grid(True, alpha=0.3)
        
        # Plot 4: Variant type distribution
        if not variants_df.empty:
            type_counts = variants_df['type'].value_counts()
            axes[1, 0].bar(type_counts.index, type_counts.values,
                        color=['lightblue', 'lightgreen', 'lightcoral', 'gold'])
            axes[1, 0].set_xlabel('Variant Type')
            axes[1, 0].set_ylabel('Count')
            axes[1, 0].set_title('Variant Type Distribution')
            axes[1, 0].grid(True, alpha=0.3, axis='y')
            
            # Add count labels
            for i, v in enumerate(type_counts.values):
                axes[1, 0].text(i, v, str(v), ha='center', va='bottom')
        
        # Plot 5: Variant impact distribution
        if not variants_df.empty:
            impact_counts = variants_df['impact'].value_counts()
            axes[1, 1].bar(impact_counts.index, impact_counts.values,
                        color=['lightgreen', 'gold', 'salmon', 'lightgray'])
            axes[1, 1].set_xlabel('Impact Type')
            axes[1, 1].set_ylabel('Count')
            axes[1, 1].set_title('Variant Impact Distribution')
            axes[1, 1].grid(True, alpha=0.3, axis='y')
            
            # Add count labels
            for i, v in enumerate(impact_counts.values):
                axes[1, 1].text(i, v, str(v), ha='center', va='bottom')
        
        # Plot 6: Variant density along genome
        if not variants_df.empty:
            window_size = 100000
            genome_length = self.genome_lengths.get("H37Rv", 4411532)
            bins = np.arange(0, genome_length + window_size, window_size)
            variant_counts, _ = np.histogram(variants_df['position'], bins=bins)
            
            axes[1, 2].bar(bins[:-1] / 1e6, variant_counts,
                        width=window_size/1e6,
                        color='steelblue',
                        edgecolor='black',
                        alpha=0.7)
            axes[1, 2].set_xlabel('Genome Position (Mbp)')
            axes[1, 2].set_ylabel('Variant Count per 100kb')
            axes[1, 2].set_title('Variant Density Along Genome')
            axes[1, 2].grid(True, alpha=0.3)
        
        plt.suptitle('M. tuberculosis H37Rv vs CDC1551 Comparative Genomics Analysis',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        # Save figure
        output_file = plots_dir / "comparative_analysis_plots.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"✓ Visualizations saved to: {output_file}")
    
    def integrate_annotations(self, variants_df):
        """Integrate gene annotations with variant analysis."""
        print("\nIntegrating gene annotations...")
        
        try:
            # Import AnnotationParser from the new module
            from annotation_parser import AnnotationParser
        except ImportError:
            print("Warning: annotation_parser module not found. Skipping annotation integration.")
            return variants_df
        
        # Initialize annotation parser
        annotator = AnnotationParser(self.project_root)
        
        # Annotate variants for both strains
        print("Annotating variants with gene information...")
        
        # Annotate from H37Rv perspective
        annotated_variants_h37rv = annotator.annotate_variants(variants_df, "H37Rv")
        
        # Annotate from CDC1551 perspective
        annotated_variants_cdc = annotator.annotate_variants(variants_df, "CDC1551")
        
        # Merge annotations
        annotated_variants = annotated_variants_h37rv.copy()
        annotated_variants['gene_annotation_cdc'] = annotated_variants_cdc['gene_annotation']
        annotated_variants['gene_product_cdc'] = annotated_variants_cdc['gene_product']
        
        # Save annotated variants
        output_file = self.variants_dir / "annotated_variants.csv"
        annotated_variants.to_csv(output_file, index=False)
        
        print(f"✓ Annotated variants saved to: {output_file}")
        
        # Analyze annotated variants
        self.analyze_annotated_variants(annotated_variants)
        
        return annotated_variants
    
    def analyze_annotated_variants(self, annotated_variants_df):
        """Analyze variants with gene annotations."""
        print("\nAnalyzing annotated variants...")
        
        if annotated_variants_df.empty:
            print("  No annotated variants to analyze")
            return {}
        
        # Calculate annotation statistics
        total_variants = len(annotated_variants_df)
        in_gene_variants = annotated_variants_df['in_gene_region'].sum()
        intergenic_variants = total_variants - in_gene_variants
        
        # Gene-specific statistics
        gene_variants = annotated_variants_df[annotated_variants_df['in_gene_region']]
        
        stats = {
            'total_variants': total_variants,
            'variants_in_genes': in_gene_variants,
            'variants_intergenic': intergenic_variants,
            'percent_in_genes': (in_gene_variants / total_variants * 100) if total_variants > 0 else 0,
            'unique_genes_affected': gene_variants['gene_annotation'].nunique() if not gene_variants.empty else 0,
            'average_variants_per_gene': (in_gene_variants / gene_variants['gene_annotation'].nunique()) 
                                        if not gene_variants.empty and gene_variants['gene_annotation'].nunique() > 0 else 0
        }
        
        # Print statistics
        print(f"  Variants in gene regions: {stats['variants_in_genes']:,} ({stats['percent_in_genes']:.1f}%)")
        print(f"  Intergenic variants: {stats['variants_intergenic']:,}")
        print(f"  Unique genes affected: {stats['unique_genes_affected']:,}")
        print(f"  Average variants per gene: {stats['average_variants_per_gene']:.2f}")
        
        # Find genes with most variants
        if not gene_variants.empty:
            gene_variant_counts = gene_variants['gene_annotation'].value_counts().head(10)
            print("\n  Top 10 genes with most variants:")
            for gene, count in gene_variant_counts.items():
                if gene != 'Intergenic':
                    print(f"    {gene}: {count} variants")
        
        # Save statistics
        stats_df = pd.DataFrame([stats])
        stats_df.to_csv(self.variants_dir / "annotation_statistics.csv", index=False)
        
        # Generate gene variant plot
        self.create_gene_variant_plot(gene_variants)
        
        return stats
    
    def create_gene_variant_plot(self, gene_variants_df):
        """Create visualization for gene variant analysis."""
        if gene_variants_df.empty:
            return
        
        # Set style
        plt.style.use('seaborn-v0_8-whitegrid')
        
        # Create figure
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Plot 1: Variant distribution by gene type
        if 'gene_type' in gene_variants_df.columns:
            type_counts = gene_variants_df['gene_type'].value_counts()
            axes[0].bar(type_counts.index, type_counts.values, color=['lightblue', 'lightgreen', 'salmon'])
            axes[0].set_xlabel('Gene Type')
            axes[0].set_ylabel('Variant Count')
            axes[0].set_title('Variants by Gene Type')
            axes[0].tick_params(axis='x', rotation=45)
        
        # Plot 2: Top 10 genes with most variants
        top_genes = gene_variants_df['gene_annotation'].value_counts().head(10)
        if len(top_genes) > 0:
            axes[1].barh(range(len(top_genes)), top_genes.values)
            axes[1].set_yticks(range(len(top_genes)))
            axes[1].set_yticklabels(top_genes.index)
            axes[1].set_xlabel('Variant Count')
            axes[1].set_title('Top 10 Genes with Most Variants')
            axes[1].invert_yaxis()
        
        plt.suptitle('Gene Annotation Analysis of Variants', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        # Save figure
        output_file = self.plots_dir / "gene_variant_analysis.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"✓ Gene variant analysis plot saved to: {output_file}")

    def generate_enhanced_report(self, alignment_stats, variant_stats, annotation_stats=None):
        """Generate comprehensive analysis report with annotation data."""
        from datetime import datetime
        
        # Format Ti/Tv ratio properly
        ti_tv_ratio = variant_stats.get('ti_tv_ratio', 'N/A')
        if isinstance(ti_tv_ratio, (int, float)):
            ti_tv_str = f"{ti_tv_ratio:.2f}"
        else:
            ti_tv_str = str(ti_tv_ratio)

        # Add annotation section if available
        annotation_section = ""
        if annotation_stats:
            annotation_section = f"""
### Gene Annotation Analysis
- **Variants in gene regions:** {annotation_stats.get('variants_in_genes', 0):,} ({annotation_stats.get('percent_in_genes', 0):.1f}%)
- **Intergenic variants:** {annotation_stats.get('variants_intergenic', 0):,}
- **Unique genes affected:** {annotation_stats.get('unique_genes_affected', 0):,}
- **Average variants per gene:** {annotation_stats.get('average_variants_per_gene', 0):.2f}

### Annotation Insights

The integration of GFF3 and GBFF annotations provides biological context to the variant analysis:
- Identified specific genes most affected by variants
- Provided functional annotations for variant locations
- Enabled analysis of variant distribution across different gene types
"""
        else:
            annotation_section = """
### Gene Annotation Analysis
*Annotation analysis was not performed in this run. Run with annotation integration enabled.*
"""

        report_content = f"""# M. tuberculosis Comparative Genomics Analysis Report

## Overview
Comparative genomics analysis between Mycobacterium tuberculosis strains H37Rv and CDC1551.

**Analysis Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
**Project Directory:** {self.project_root}

## Methods

### Data Sources
- H37Rv: FASTA, GBFF, and GFF3 files from NCBI
- CDC1551: FASTA, GBFF, and GFF3 files from NCBI

### Analysis Pipeline
1. Genome data acquisition and validation
2. Basic genome statistics calculation
3. Whole genome alignment simulation
4. Variant identification and characterization
5. Gene annotation parsing (GFF3/GBFF)
6. Functional impact prediction with gene context

## Results

### Genome Alignment (Simulated)
- **Total aligned bases:** {alignment_stats.get('total_aligned_bp', 0):,} bp
- **H37Rv coverage:** {alignment_stats.get('h37rv_coverage', 0):.2f}%
- **CDC1551 coverage:** {alignment_stats.get('cdc1551_coverage', 0):.2f}%
- **Alignment blocks:** {alignment_stats.get('num_blocks', 0)}
- **Mean identity:** {alignment_stats.get('mean_identity', 0):.2f}%

### Genomic Variants (Simulated)
- **Total variants:** {variant_stats.get('total_variants', 0):,}
- **SNPs:** {variant_stats.get('snps', 0):,} ({variant_stats.get('snps', 0)/variant_stats.get('total_variants', 1)*100:.1f}%)
- **Insertions:** {variant_stats.get('insertions', 0):,}
- **Deletions:** {variant_stats.get('deletions', 0):,}
- **Ti/Tv ratio:** {ti_tv_str}
{annotation_section}
### Variant Impact (Gene Context)
- **Synonymous:** {variant_stats.get('synonymous', 0):,}
- **Missense:** {variant_stats.get('missense', 0):,}
- **Nonsense:** {variant_stats.get('nonsense', 0):,}
- **Intergenic:** {variant_stats.get('intergenic', 0):,}

## Key Findings

1. **High Genomic Similarity:** The two strains show >99% identity in aligned regions.
2. **Limited Structural Variation:** Few large rearrangements detected.
3. **Variant Distribution:** Most variants are SNPs, with Ti/Tv ratio consistent with purifying selection.
4. **Functional Impact:** Majority of coding variants are synonymous, suggesting conserved protein function.

## Files Generated

### Alignment Data
- `simulated_alignment_coordinates.csv`: Alignment block coordinates
- `alignment_statistics.csv`: Alignment statistics summary

### Variant Data
- `simulated_variants.csv`: Variant calls with annotations
- `annotated_variants.csv`: Variants with gene annotations
- `variant_statistics.csv`: Variant statistics summary
- `annotation_statistics.csv`: Gene annotation statistics

### Annotation Data
- `H37Rv_gff_annotations.csv`: Parsed GFF3 annotations for H37Rv
- `H37Rv_gbff_annotations.csv`: Parsed GBFF annotations for H37Rv
- `CDC1551_gff_annotations.csv`: Parsed GFF3 annotations for CDC1551
- `CDC1551_gbff_annotations.csv`: Parsed GBFF annotations for CDC1551

### Visualizations
- `comparative_analysis_plots.png`: 6-panel visualization of results
- `gene_variant_analysis.png`: Gene-specific variant analysis

## Next Steps
1. Validate findings with experimental data
2. Perform functional enrichment analysis on affected genes
3. Correlate genomic differences with phenotypic variations
4. Expand analysis to include additional Mtb strains

---
*Report generated automatically by enhanced alignment_analysis.py*
"""
        
        # Save report
        report_file = self.results_dir / "enhanced_comparative_genomics_report.md"
        with open(report_file, 'w') as f:
            f.write(report_content)
        
        print(f"✓ Enhanced analysis report saved to: {report_file}")
    
    def run_complete_analysis(self):
        """Run complete comparative genomics analysis."""
        print("=" * 80)
        print("M. TUBERCULOSIS COMPARATIVE GENOMICS ANALYSIS")
        print("=" * 80)
        
        # Step 1: Simulate alignment
        alignment_df = self.simulate_alignment()
        
        # Step 2: Analyze alignment
        alignment_stats = self.analyze_alignment(alignment_df)
        
        # Step 3: Simulate variants
        variants_df = self.simulate_variants()
        
        # Step 4: Analyze variants
        variant_stats = self.analyze_variants(variants_df)
        
        # Step 5: Integrate annotations (NEW STEP)
        try:
            annotated_variants = self.integrate_annotations(variants_df)
            annotation_stats = self.analyze_annotated_variants(annotated_variants)
        except Exception as e:
            print(f"\nNote: Annotation integration failed: {e}")
            print("Continuing without annotation analysis...")
            annotation_stats = None
        
        # Step 6: Create visualizations
        self.create_visualizations(alignment_df, variants_df)
        
        # Step 7: Generate report
        if annotation_stats:
            self.generate_enhanced_report(alignment_stats, variant_stats, annotation_stats)
        else:
            # Fall back to original report
            self.generate_report(alignment_stats, variant_stats)
        
        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE")
        print("=" * 80)
        print(f"\nResults available in:")
        print(f"  • Alignments: {self.alignments_dir.relative_to(self.project_root)}")
        print(f"  • Variants: {self.variants_dir.relative_to(self.project_root)}")
        print(f"  • Plots: {self.plots_dir.relative_to(self.project_root)}")

def main():
    """Command-line interface for comparative analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Comparative genomics analysis of M. tuberculosis strains"
    )
    parser.add_argument(
        '--project-dir',
        type=str,
        default=None,
        help='Path to project root directory'
    )
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = ComparativeAlignment(args.project_dir)
    
    # Run complete analysis
    analyzer.run_complete_analysis()

if __name__ == "__main__":
    main()