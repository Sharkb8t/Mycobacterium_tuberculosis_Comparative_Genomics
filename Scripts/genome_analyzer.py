"""
Author: Dalton A. Schmidt
GitHub: https://github.com/Sharkb8t

Genome analysis utilities for M. tuberculosis comparative genomics.
Provides functions for genome statistics, GC analysis, and basic comparisons.
"""
from Bio import SeqIO
from Bio.SeqUtils import GC, GC_skew, molecular_weight
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

class GenomeAnalyzer:
    """Analyze genome sequences and generate comparative statistics."""
    
    def __init__(self, project_root: Optional[Path] = None):
        """
        Initialize the genome analyzer.
        
        Args:
            project_root: Path to project root directory
        """
        if project_root is None:
            self.project_root = Path.cwd()
        else:
            self.project_root = Path(project_root)
        
        # Define paths
        self.raw_dir = self.project_root / "Data" / "raw"
        self.results_dir = self.project_root / "Results"
        
        # Create results directories if they don't exist
        self.results_dir.mkdir(exist_ok=True)
    
    def load_genome(self, strain_name: str) -> List:
        """
        Load genome sequences for a strain.
        
        Args:
            strain_name: Name of the strain (e.g., "H37Rv")
        
        Returns:
            List of Bio.SeqRecord objects
        """
        fasta_file = self.raw_dir / f"{strain_name}_fasta"
        
        if not fasta_file.exists():
            raise FileNotFoundError(f"FASTA file not found for {strain_name}: {fasta_file}")
        
        records = list(SeqIO.parse(fasta_file, "fasta"))
        
        if not records:
            raise ValueError(f"No sequences found in {fasta_file}")
        
        return records
    
    def calculate_basic_stats(self, records: List) -> Dict:
        """
        Calculate basic genome statistics.
        
        Args:
            records: List of Bio.SeqRecord objects
        
        Returns:
            Dictionary with genome statistics
        """
        if not records:
            return {}
        
        # Calculate lengths
        lengths = [len(rec.seq) for rec in records]
        total_length = sum(lengths)
        
        # Calculate GC content
        gc_total = 0
        for rec in records:
            gc_total += GC(rec.seq) * len(rec.seq)
        gc_content = gc_total / total_length
        
        # Sort lengths for N50 calculation
        sorted_lengths = sorted(lengths, reverse=True)
        
        # Calculate N50, N75, L50
        cumulative_sum = 0
        n50 = 0
        l50 = 0
        n75 = 0
        
        for i, length in enumerate(sorted_lengths):
            cumulative_sum += length
            if cumulative_sum >= total_length * 0.5 and n50 == 0:
                n50 = length
                l50 = i + 1
            if cumulative_sum >= total_length * 0.75 and n75 == 0:
                n75 = length
        
        stats = {
            'num_contigs': len(records),
            'total_length': total_length,
            'average_length': np.mean(lengths),
            'median_length': np.median(lengths),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'gc_content': gc_content,
            'n50': n50,
            'n75': n75,
            'l50': l50,
            'contig_names': [rec.id for rec in records],
            'contig_lengths': lengths
        }
        
        return stats
    
    def calculate_sliding_window_stats(self, sequence: str, window_size: int = 1000) -> pd.DataFrame:
        """
        Calculate sliding window statistics for a genome sequence.
        
        Args:
            sequence: DNA sequence as string
            window_size: Size of sliding window
        
        Returns:
            DataFrame with window statistics
        """
        sequence = sequence.upper()
        data = []
        
        for i in range(0, len(sequence) - window_size + 1, window_size):
            window = sequence[i:i + window_size]
            
            if len(window) < window_size:
                continue
            
            # Calculate statistics for this window
            a_count = window.count('A')
            c_count = window.count('C')
            g_count = window.count('G')
            t_count = window.count('T')
            n_count = window.count('N')
            
            total = a_count + c_count + g_count + t_count + n_count
            
            if total > 0:
                gc_content = (g_count + c_count) / total * 100
                
                # Calculate GC skew
                if g_count + c_count > 0:
                    gc_skew_val = (g_count - c_count) / (g_count + c_count)
                else:
                    gc_skew_val = 0
                
                # Calculate AT skew
                if a_count + t_count > 0:
                    at_skew_val = (a_count - t_count) / (a_count + t_count)
                else:
                    at_skew_val = 0
            else:
                gc_content = 0
                gc_skew_val = 0
                at_skew_val = 0
            
            data.append({
                'position': i,
                'gc_content': gc_content,
                'gc_skew': gc_skew_val,
                'at_skew': at_skew_val,
                'a_content': a_count / len(window) * 100 if total > 0 else 0,
                'c_content': c_count / len(window) * 100 if total > 0 else 0,
                'g_content': g_count / len(window) * 100 if total > 0 else 0,
                't_content': t_count / len(window) * 100 if total > 0 else 0,
                'n_content': n_count / len(window) * 100
            })
        
        return pd.DataFrame(data)
    
    def calculate_dinucleotide_frequencies(self, sequence: str) -> pd.DataFrame:
        """
        Calculate dinucleotide frequencies in a genome.
        
        Args:
            sequence: DNA sequence as string
        
        Returns:
            DataFrame with dinucleotide frequencies
        """
        sequence = sequence.upper()
        nucleotides = ['A', 'C', 'G', 'T']
        dinucleotides = [a + b for a in nucleotides for b in nucleotides]
        
        # Initialize counts
        counts = {dinuc: 0 for dinuc in dinucleotides}
        total_dinucleotides = 0
        
        # Count dinucleotides
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            if dinuc in counts:
                counts[dinuc] += 1
                total_dinucleotides += 1
        
        # Calculate frequencies
        frequencies = {}
        if total_dinucleotides > 0:
            frequencies = {dinuc: count / total_dinucleotides for dinuc, count in counts.items()}
        
        # Create DataFrame
        df = pd.DataFrame({
            'dinucleotide': list(frequencies.keys()),
            'frequency': list(frequencies.values()),
            'count': [counts[d] for d in frequencies.keys()]
        })
        
        # Add expected frequency based on mononucleotide frequencies
        mono_counts = {nuc: sequence.count(nuc) for nuc in nucleotides}
        total_mono = sum(mono_counts.values())
        
        if total_mono > 0:
            mono_freq = {nuc: count / total_mono for nuc, count in mono_counts.items()}
            df['expected_frequency'] = df['dinucleotide'].apply(
                lambda x: mono_freq[x[0]] * mono_freq[x[1]] if total_mono > 0 else 0
            )
            df['odds_ratio'] = df['frequency'] / df['expected_frequency']
        
        return df
    
    def compare_genomes(self, strain1: str, strain2: str) -> Dict:
        """
        Compare two genomes and calculate differences.
        
        Args:
            strain1: Name of first strain
            strain2: Name of second strain
        
        Returns:
            Dictionary with comparison results
        """
        # Load genomes
        records1 = self.load_genome(strain1)
        records2 = self.load_genome(strain2)
        
        # Calculate statistics
        stats1 = self.calculate_basic_stats(records1)
        stats2 = self.calculate_basic_stats(records2)
        
        # Calculate differences
        comparison = {
            'strain1': strain1,
            'strain2': strain2,
            'length_difference': abs(stats1['total_length'] - stats2['total_length']),
            'length_ratio': stats1['total_length'] / stats2['total_length'] if stats2['total_length'] > 0 else 0,
            'gc_difference': abs(stats1['gc_content'] - stats2['gc_content']),
            'contig_difference': abs(stats1['num_contigs'] - stats2['num_contigs']),
            'stats1': stats1,
            'stats2': stats2
        }
        
        return comparison
    
    def generate_report(self, strain_names: List[str]) -> str:
        """
        Generate a text report of genome analysis.
        
        Args:
            strain_names: List of strain names to analyze
        
        Returns:
            Formatted report string
        """
        report_lines = []
        report_lines.append("=" * 70)
        report_lines.append("MYCOBACTERIUM TUBERCULOSIS GENOME ANALYSIS REPORT")
        report_lines.append("=" * 70)
        report_lines.append(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append(f"Project Directory: {self.project_root}")
        report_lines.append("")
        
        # Analyze each strain
        all_stats = {}
        for strain in strain_names:
            try:
                records = self.load_genome(strain)
                stats = self.calculate_basic_stats(records)
                all_stats[strain] = stats
                
                report_lines.append(f"STRAIN: {strain}")
                report_lines.append("-" * 40)
                report_lines.append(f"Number of contigs: {stats['num_contigs']}")
                report_lines.append(f"Total length: {stats['total_length']:,} bp")
                report_lines.append(f"Average contig length: {stats['average_length']:,.0f} bp")
                report_lines.append(f"GC content: {stats['gc_content']:.2f}%")
                report_lines.append(f"N50: {stats['n50']:,} bp")
                report_lines.append(f"N75: {stats['n75']:,} bp")
                report_lines.append(f"L50: {stats['l50']}")
                report_lines.append("")
                
            except Exception as e:
                report_lines.append(f"ERROR analyzing {strain}: {str(e)}")
                report_lines.append("")
        
        # Compare strains if we have at least two
        if len(strain_names) >= 2 and len(all_stats) >= 2:
            report_lines.append("COMPARISON BETWEEN STRAINS")
            report_lines.append("-" * 40)
            
            strains = list(all_stats.keys())
            for i in range(len(strains)):
                for j in range(i + 1, len(strains)):
                    strain1 = strains[i]
                    strain2 = strains[j]
                    
                    comp = self.compare_genomes(strain1, strain2)
                    
                    report_lines.append(f"{strain1} vs {strain2}:")
                    report_lines.append(f"  Length difference: {comp['length_difference']:,} bp")
                    report_lines.append(f"  GC difference: {comp['gc_difference']:.3f}%")
                    report_lines.append(f"  Contig count difference: {comp['contig_difference']}")
                    report_lines.append("")
        
        report_lines.append("=" * 70)
        report_lines.append("END OF REPORT")
        report_lines.append("=" * 70)
        
        return "\n".join(report_lines)
    
    def save_stats_to_csv(self, strain_names: List[str], output_dir: Optional[Path] = None):
        """
        Save genome statistics to CSV files.
        
        Args:
            strain_names: List of strain names to analyze
            output_dir: Directory to save CSV files (default: Results directory)
        """
        if output_dir is None:
            output_dir = self.results_dir
        
        output_dir.mkdir(exist_ok=True)
        
        all_data = []
        
        for strain in strain_names:
            try:
                records = self.load_genome(strain)
                stats = self.calculate_basic_stats(records)
                
                # Save individual strain stats
                strain_df = pd.DataFrame([stats])
                strain_df.to_csv(output_dir / f"{strain}_basic_stats.csv", index=False)
                
                # Collect for combined file
                stats['strain'] = strain
                all_data.append(stats)
                
                print(f"✓ Saved {strain} statistics to {output_dir / f'{strain}_basic_stats.csv'}")
                
            except Exception as e:
                print(f"✗ Error processing {strain}: {e}")
        
        # Save combined statistics
        if all_data:
            combined_df = pd.DataFrame(all_data)
            combined_df.to_csv(output_dir / "all_strains_basic_stats.csv", index=False)
            print(f"✓ Saved combined statistics to {output_dir / 'all_strains_basic_stats.csv'}")

def main():
    """Command-line interface for genome analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Analyze M. tuberculosis genome sequences"
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
        help='Strain names to analyze (default: H37Rv CDC1551)'
    )
    parser.add_argument(
        '--output-report',
        type=str,
        default='genome_analysis_report.txt',
        help='Output file for analysis report'
    )
    parser.add_argument(
        '--save-csv',
        action='store_true',
        help='Save statistics to CSV files'
    )
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = GenomeAnalyzer(args.project_dir)
    
    # Generate report
    report = analyzer.generate_report(args.strains)
    
    # Print report to console
    print(report)
    
    # Save report to file
    report_path = Path(args.output_report)
    with open(report_path, 'w') as f:
        f.write(report)
    
    print(f"\n✓ Report saved to: {report_path}")
    
    # Save CSV files if requested
    if args.save_csv:
        analyzer.save_stats_to_csv(args.strains)

if __name__ == "__main__":
    main()