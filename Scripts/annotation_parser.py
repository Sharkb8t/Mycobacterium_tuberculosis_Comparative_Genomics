"""
Author: Dalton A. Schmidt
GitHub: https://github.com/Sharkb8t

Annotation file parsers for M. tuberculosis genomic data.
Handles GFF3 and GBFF file formats to extract gene annotations.
"""
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')

class AnnotationParser:
    """Parse and process genome annotation files (GFF3, GBFF)."""
    
    def __init__(self, project_root: Optional[Path] = None):
        """
        Initialize the annotation parser.
        
        Args:
            project_root: Path to project root directory
        """
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
        self.annotations_dir = self.results_dir / "Annotations"
        
        # Create results directories if they don't exist
        self.annotations_dir.mkdir(parents=True, exist_ok=True)
    
    def parse_gff3(self, strain_name: str) -> pd.DataFrame:
        """
        Parse GFF3 annotation file.
        
        Args:
            strain_name: Name of the strain
        
        Returns:
            DataFrame with gene annotations
        """
        gff_file = self.raw_dir / f"{strain_name}_gff.gff"
        
        if not gff_file.exists():
            raise FileNotFoundError(f"GFF3 file not found for {strain_name}")
        
        annotations = []
        
        with open(gff_file, 'r') as f:
            for line in f:
                # Skip comment lines and empty lines
                if line.startswith('#') or line.strip() == '':
                    continue
                
                parts = line.strip().split('\t')
                
                if len(parts) < 9:
                    continue
                
                # Parse basic fields
                seqid = parts[0]
                source = parts[1]
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                score = parts[5] if parts[5] != '.' else None
                strand = parts[6]
                phase = parts[7] if parts[7] != '.' else None
                attributes = parts[8]
                
                # Parse attributes into dictionary
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key.strip()] = value.strip()
                
                # Only process gene/CDS features
                if feature_type in ['gene', 'CDS', 'rRNA', 'tRNA']:
                    annotation = {
                        'strain': strain_name,
                        'seqid': seqid,
                        'source': source,
                        'type': feature_type,
                        'start': start,
                        'end': end,
                        'length': end - start + 1,
                        'strand': strand,
                        'score': score,
                        'phase': phase,
                        'gene_id': attr_dict.get('ID', attr_dict.get('gene_id', '')),
                        'name': attr_dict.get('Name', attr_dict.get('gene', '')),
                        'product': attr_dict.get('product', attr_dict.get('Note', '')),
                        'locus_tag': attr_dict.get('locus_tag', ''),
                        'gene_biotype': attr_dict.get('gene_biotype', 'protein_coding'),
                        'attributes': attributes
                    }
                    annotations.append(annotation)
        
        df = pd.DataFrame(annotations)
        
        if not df.empty:
            # Save parsed annotations
            output_file = self.annotations_dir / f"{strain_name}_gff_annotations.csv"
            df.to_csv(output_file, index=False)
            print(f"✓ Parsed {len(df)} features from GFF3 for {strain_name}")
        
        return df
    
    def parse_gbff(self, strain_name: str) -> pd.DataFrame:
        """
        Parse GBFF (GenBank Flat File) annotation file.
        
        Args:
            strain_name: Name of the strain
        
        Returns:
            DataFrame with gene annotations
        """
        gbff_file = self.raw_dir / f"{strain_name}_gbff.gbff"
        
        if not gbff_file.exists():
            raise FileNotFoundError(f"GBFF file not found for {strain_name}")
        
        annotations = []
        record = SeqIO.read(gbff_file, "genbank")
        
        for feature in record.features:
            if feature.type in ['gene', 'CDS', 'rRNA', 'tRNA']:
                # Extract basic information
                annotation = {
                    'strain': strain_name,
                    'seqid': record.id,
                    'source': 'GenBank',
                    'type': feature.type,
                    'start': int(feature.location.start) + 1,  # Convert to 1-based
                    'end': int(feature.location.end),
                    'length': int(feature.location.end) - int(feature.location.start),
                    'strand': feature.location.strand,
                    'gene_id': feature.qualifiers.get('gene_id', [''])[0],
                    'name': feature.qualifiers.get('gene', [''])[0],
                    'locus_tag': feature.qualifiers.get('locus_tag', [''])[0],
                    'product': feature.qualifiers.get('product', [''])[0],
                    'protein_id': feature.qualifiers.get('protein_id', [''])[0],
                    'note': '; '.join(feature.qualifiers.get('note', [])),
                    'function': feature.qualifiers.get('function', [''])[0],
                    'db_xref': '; '.join(feature.qualifiers.get('db_xref', [])),
                    'transl_table': feature.qualifiers.get('transl_table', [''])[0],
                    'translation': feature.qualifiers.get('translation', [''])[0] if feature.type == 'CDS' else ''
                }
                
                # Handle pseudogenes
                if 'pseudo' in feature.qualifiers:
                    annotation['pseudogene'] = True
                    annotation['gene_biotype'] = 'pseudogene'
                else:
                    annotation['pseudogene'] = False
                    annotation['gene_biotype'] = 'protein_coding' if feature.type == 'CDS' else feature.type
                
                annotations.append(annotation)
        
        df = pd.DataFrame(annotations)
        
        if not df.empty:
            # Save parsed annotations
            output_file = self.annotations_dir / f"{strain_name}_gbff_annotations.csv"
            df.to_csv(output_file, index=False)
            print(f"✓ Parsed {len(df)} features from GBFF for {strain_name}")
        
        return df
    
    def get_gene_annotations(self, strain_name: str) -> pd.DataFrame:
        """
        Get gene annotations from either GFF3 or GBFF file.
        Prioritizes GFF3 if available, falls back to GBFF.
        
        Args:
            strain_name: Name of the strain
        
        Returns:
            DataFrame with gene annotations
        """
        try:
            # Try GFF3 first
            gff_file = self.raw_dir / f"{strain_name}_gff.gff"
            if gff_file.exists():
                return self.parse_gff3(strain_name)
        except Exception as e:
            print(f"Warning: Could not parse GFF3 for {strain_name}: {e}")
        
        try:
            # Fall back to GBFF
            return self.parse_gbff(strain_name)
        except Exception as e:
            print(f"Error: Could not parse GBFF for {strain_name}: {e}")
            return pd.DataFrame()
    
    def calculate_gene_statistics(self, annotations_df: pd.DataFrame) -> Dict:
        """
        Calculate statistics from gene annotations.
        
        Args:
            annotations_df: DataFrame with gene annotations
        
        Returns:
            Dictionary with gene statistics
        """
        if annotations_df.empty:
            return {}
        
        # Filter for genes and CDS features
        gene_features = annotations_df[annotations_df['type'].isin(['gene', 'CDS'])]
        
        if gene_features.empty:
            return {}
        
        stats = {
            'total_genes': len(gene_features),
            'mean_gene_length': gene_features['length'].mean(),
            'median_gene_length': gene_features['length'].median(),
            'total_coding_length': gene_features['length'].sum(),
            'min_gene_length': gene_features['length'].min(),
            'max_gene_length': gene_features['length'].max(),
            'forward_strand_genes': (gene_features['strand'] == 1).sum() if 'strand' in gene_features.columns else 0,
            'reverse_strand_genes': (gene_features['strand'] == -1).sum() if 'strand' in gene_features.columns else 0,
            'unique_gene_types': gene_features['type'].nunique(),
            'genes_with_product': gene_features['product'].notna().sum() if 'product' in gene_features.columns else 0,
            'pseudogenes': gene_features.get('pseudogene', pd.Series([False])).sum() if 'pseudogene' in gene_features.columns else 0
        }
        
        # Calculate gene density
        if 'start' in gene_features.columns and 'end' in gene_features.columns:
            genome_length = gene_features['end'].max() if not gene_features.empty else 0
            if genome_length > 0:
                stats['gene_density'] = stats['total_genes'] / (genome_length / 1000)  # genes per kb
                stats['coding_percentage'] = (stats['total_coding_length'] / genome_length) * 100
        
        # Get gene type distribution
        type_counts = gene_features['type'].value_counts().to_dict()
        stats['gene_type_distribution'] = type_counts
        
        return stats
    
    def compare_gene_annotations(self, strain1: str, strain2: str) -> Dict:
        """
        Compare gene annotations between two strains.
        
        Args:
            strain1: Name of first strain
            strain2: Name of second strain
        
        Returns:
            Dictionary with comparison results
        """
        # Get annotations for both strains
        annotations1 = self.get_gene_annotations(strain1)
        annotations2 = self.get_gene_annotations(strain2)
        
        if annotations1.empty or annotations2.empty:
            print(f"Warning: Missing annotations for {strain1} or {strain2}")
            return {}
        
        # Filter for CDS features (protein-coding genes)
        cds1 = annotations1[annotations1['type'] == 'CDS'].copy()
        cds2 = annotations2[annotations2['type'] == 'CDS'].copy()
        
        # Get unique gene identifiers
        # Try locus_tag first, then gene_id, then name
        if 'locus_tag' in cds1.columns and 'locus_tag' in cds2.columns:
            cds1['gene_key'] = cds1['locus_tag']
            cds2['gene_key'] = cds2['locus_tag']
        elif 'gene_id' in cds1.columns and 'gene_id' in cds2.columns:
            cds1['gene_key'] = cds1['gene_id']
            cds2['gene_key'] = cds2['gene_id']
        else:
            cds1['gene_key'] = cds1['name']
            cds2['gene_key'] = cds2['name']
        
        # Find common and unique genes
        genes1 = set(cds1['gene_key'].dropna().unique())
        genes2 = set(cds2['gene_key'].dropna().unique())
        
        common_genes = genes1.intersection(genes2)
        unique_to_strain1 = genes1 - genes2
        unique_to_strain2 = genes2 - genes1
        
        comparison = {
            'strain1': strain1,
            'strain2': strain2,
            'total_genes_strain1': len(genes1),
            'total_genes_strain2': len(genes2),
            'common_genes': len(common_genes),
            'unique_to_strain1': len(unique_to_strain1),
            'unique_to_strain2': len(unique_to_strain2),
            'common_genes_percentage': (len(common_genes) / len(genes1) * 100) if len(genes1) > 0 else 0,
            'unique_genes_list_strain1': list(unique_to_strain1),
            'unique_genes_list_strain2': list(unique_to_strain2),
            'gene_similarity': len(common_genes) / max(len(genes1), len(genes2)) if max(len(genes1), len(genes2)) > 0 else 0
        }
        
        # Save comparison results
        output_file = self.annotations_dir / f"gene_comparison_{strain1}_vs_{strain2}.csv"
        comparison_df = pd.DataFrame([comparison])
        comparison_df.to_csv(output_file, index=False)
        
        return comparison
    
    def annotate_variants(self, variants_df: pd.DataFrame, strain: str) -> pd.DataFrame:
        """
        Annotate variants with gene information.
        
        Args:
            variants_df: DataFrame with variants
            strain: Strain name for annotation reference
        
        Returns:
            DataFrame with annotated variants
        """
        if variants_df.empty:
            return variants_df
        
        # Get gene annotations for the strain
        annotations = self.get_gene_annotations(strain)
        
        if annotations.empty:
            print(f"Warning: No annotations found for {strain}")
            variants_df['gene_annotation'] = 'Unknown'
            variants_df['gene_product'] = ''
            variants_df['in_gene_region'] = False
            return variants_df
        
        # Filter for gene regions
        gene_regions = annotations[annotations['type'].isin(['gene', 'CDS'])]
        
        if gene_regions.empty:
            variants_df['gene_annotation'] = 'Unknown'
            variants_df['gene_product'] = ''
            variants_df['in_gene_region'] = False
            return variants_df
        
        # Annotate each variant
        annotated_variants = []
        
        for _, variant in variants_df.iterrows():
            pos = variant['position']
            
            # Find overlapping genes
            overlapping_genes = gene_regions[
                (gene_regions['start'] <= pos) & (gene_regions['end'] >= pos)
            ]
            
            if not overlapping_genes.empty:
                # Get the first overlapping gene
                gene = overlapping_genes.iloc[0]
                
                variant_annotated = variant.copy()
                variant_annotated['gene_annotation'] = gene.get('name', gene.get('locus_tag', 'Unknown'))
                variant_annotated['gene_product'] = gene.get('product', '')
                variant_annotated['in_gene_region'] = True
                variant_annotated['gene_start'] = gene['start']
                variant_annotated['gene_end'] = gene['end']
                variant_annotated['gene_strand'] = gene.get('strand', '')
                variant_annotated['gene_type'] = gene['type']
                
                # Determine if variant affects coding sequence
                if gene['type'] == 'CDS':
                    # Calculate position within CDS
                    cds_pos = pos - gene['start'] + 1
                    if gene.get('strand') == -1:
                        cds_pos = gene['end'] - pos + 1
                    variant_annotated['cds_position'] = cds_pos
                else:
                    variant_annotated['cds_position'] = None
            else:
                variant_annotated = variant.copy()
                variant_annotated['gene_annotation'] = 'Intergenic'
                variant_annotated['gene_product'] = ''
                variant_annotated['in_gene_region'] = False
                variant_annotated['gene_start'] = None
                variant_annotated['gene_end'] = None
                variant_annotated['gene_strand'] = ''
                variant_annotated['gene_type'] = ''
                variant_annotated['cds_position'] = None
            
            annotated_variants.append(variant_annotated)
        
        return pd.DataFrame(annotated_variants)


def main():
    """Command-line interface for annotation parsing."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Parse and analyze genome annotation files"
    )
    parser.add_argument(
        '--project-dir',
        type=str,
        default=None,
        help='Path to project root directory'
    )
    parser.add_argument(
        '--strains',
        nargs='+',
        default=['H37Rv', 'CDC1551'],
        help='Strain names to analyze'
    )
    parser.add_argument(
        '--parse-all',
        action='store_true',
        help='Parse both GFF3 and GBFF files for each strain'
    )
    parser.add_argument(
        '--compare',
        action='store_true',
        help='Compare gene annotations between strains'
    )
    
    args = parser.parse_args()
    
    # Initialize parser
    parser_obj = AnnotationParser(args.project_dir)
    
    if args.parse_all:
        for strain in args.strains:
            print(f"\nParsing annotations for {strain}:")
            try:
                gff_result = parser_obj.parse_gff3(strain)
                print(f"  GFF3: {len(gff_result)} features")
            except Exception as e:
                print(f"  GFF3: Error - {e}")
            
            try:
                gbff_result = parser_obj.parse_gbff(strain)
                print(f"  GBFF: {len(gbff_result)} features")
            except Exception as e:
                print(f"  GBFF: Error - {e}")
            
            # Calculate statistics
            annotations = parser_obj.get_gene_annotations(strain)
            if not annotations.empty:
                stats = parser_obj.calculate_gene_statistics(annotations)
                print(f"  Statistics: {stats.get('total_genes', 0)} genes, "
                    f"{stats.get('mean_gene_length', 0):.0f} bp average length")
    
    if args.compare and len(args.strains) >= 2:
        print(f"\nComparing gene annotations:")
        comparison = parser_obj.compare_gene_annotations(args.strains[0], args.strains[1])
        if comparison:
            print(f"  Common genes: {comparison['common_genes']}")
            print(f"  Unique to {args.strains[0]}: {comparison['unique_to_strain1']}")
            print(f"  Unique to {args.strains[1]}: {comparison['unique_to_strain2']}")
            print(f"  Gene similarity: {comparison['gene_similarity']:.2%}")


if __name__ == "__main__":
    main()