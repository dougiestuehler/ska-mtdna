#!/usr/bin/env python3

import argparse
import subprocess
import dendropy
import pandas as pd
from pathlib import Path
import glob
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import shutil
import sys

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ASCII art for default message
ASCII_ART = """⢠⠀⠀⠀⠀⣀⣠⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠼⡯⢛⠟⢋⠕⡡⢊⢷⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⢻⣀⠔⢁⠔⢡⠂⠚⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠣⣠⠊⡰⠁⡔⡘⣇⣀⡠⢤⠤⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠉⠑⠚⠉⠉⣯⠊⡰⢁⠌⠜⣱⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⢣⠎⡠⢃⠔⠔⢁⢽⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠶⣠⣊⣊⠴⠷⢺⠉⠏⠩⢑⠦⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡞⠠⢁⢁⠔⠙⣆⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣧⡡⡡⣫⣴⣽⣾⡲⠤⣄⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⠋⠉⠀⠈⢺⡿⢀⠈⠉
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⣧⠃⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⡃⠀⠀⠀"""

def parse_args():
    """Parse command-line arguments."""
    description = """Optimize SKA2 parameters for mtDNA haplotyping with mapping.
    Modes:
    --network (default): Uses weights wb=0.87, wh=0.06, wn=0.07, ws=0.0, wg=0.0
    --phylo: Uses weights wb=0.87, wh=0.06, wn=0.07, ws=0.0, wg=-0.1; uses ska align only.
    """
    parser = argparse.ArgumentParser(description=description)
    
    # Mode selection
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument('--network', action='store_true', help='Network mode (default)')
    mode_group.add_argument('--phylo', action='store_true', help='Phylogenetic mode')
    
    parser.add_argument('-k', type=str, default='11,13,15,17,19,21,23,25,27,29', help='Comma-separated k-mer sizes: default=11,13,15,17,19,21,23,25,27,29')
    parser.add_argument('-m', type=str, default='0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1', help='Comma-separated m values: 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1')
    parser.add_argument('--fasta-dir', type=str, required=True, help='Directory containing *.fa or *.fasta files (required)')
    parser.add_argument('--output-dir', type=str, default='output', help='Directory for output files')
    parser.add_argument('--num_top_params', type=int, default=10, help='Number of top-scoring k,m parameter sets to examine with ska map (default: 10, network mode only)')
    parser.add_argument('--wb', type=float, default=None, help='Weight for bootstrap score (network: 0.87, phylo: 0.87)')
    parser.add_argument('--wh', type=float, default=None, help='Weight for haplotypes (network: 0.06, phylo: 0.06)')
    parser.add_argument('--wn', type=float, default=None, help='Weight for valid nucleotide sites (network: 0.07, phylo: 0.07)')
    parser.add_argument('--ws', type=float, default=None, help='Weight for segregating sites (network: 0.0, phylo: 0.0)')
    parser.add_argument('--wg', type=float, default=None, help='Weight for gap proportion (network: 0.0, phylo: -0.1)')
    
    args = parser.parse_args()
    
    # Set default mode to network if neither is specified
    if not args.phylo:
        args.network = True
    
    # Validate num_top_params
    if args.num_top_params < 1:
        parser.error("--num_top_params must be a positive integer")
    
    # Convert k values to integers
    try:
        args.k = [int(k) for k in args.k.split(',')]
    except ValueError:
        parser.error("k-mer sizes must be integers")
    
    # Set default weights based on mode
    if args.network:
        args.wb = args.wb if args.wb is not None else 0.87
        args.wh = args.wh if args.wh is not None else 0.06
        args.wn = args.wn if args.wn is not None else 0.07
        args.ws = args.ws if args.ws is not None else 0.0
        args.wg = args.wg if args.wg is not None else 0.0
    else:  # phylo mode
        args.wb = args.wb if args.wb is not None else 0.87
        args.wh = args.wh if args.wh is not None else 0.06
        args.wn = args.wn if args.wn is not None else 0.07
        args.ws = args.ws if args.ws is not None else 0.0
        args.wg = args.wg if args.wg is not None else -0.1
    
    return args

def update_fasta_headers(fasta_dir):
    """Update FASTA headers to match filename (without .fa or .fasta extension)."""
    fasta_files = glob.glob(str(Path(fasta_dir) / "*.fa*"))
    if not fasta_files:
        raise FileNotFoundError(f"No *.fa or *.fasta files found in {fasta_dir}")
    
    for fasta_file in fasta_files:
        fasta_path = Path(fasta_file)
        sample_name = fasta_path.stem
        temp_file = fasta_path.with_suffix('.temp.fa')
        
        with open(fasta_path, 'r') as infile, open(temp_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    outfile.write(f'>{sample_name}\n')
                else:
                    outfile.write(line)
        
        shutil.move(temp_file, fasta_path)
        logging.info(f"Updated headers in {fasta_path} to {sample_name}")
    
    return [Path(f) for f in fasta_files]

def validate_alignment(aln_file):
    """Check if alignment file is valid and compute metrics, counting only valid sites."""
    if not aln_file.is_file() or aln_file.stat().st_size == 0:
        return False, "Alignment file is empty", 0, 0.0, 0.0, 0, 0
    
    sequences = []
    headers = []
    with open(aln_file) as f:
        seq = ""
        header = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    headers.append(header)
                header = line.strip()[1:]
                seq = ""
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
            headers.append(header)
    
    if not sequences:
        return False, "No sequences in alignment", 0, 0.0, 0.0, 0, 0
    
    seq_count = len(sequences)
    if seq_count < 2:
        return False, "Too few sequences", seq_count, 0.0, 0.0, 0, 0
    
    seq_len = len(sequences[0])
    if not all(len(s) == seq_len for s in sequences):
        return False, "Sequences have unequal lengths", seq_count, 0.0, 0.0, 0, 0
    
    segregating_sites = 0
    gap_count = 0
    n_count = 0
    valid_nucleotide = 0
    total_bases = seq_count * seq_len
    valid_bases = {'A', 'T', 'C', 'G'}
    
    for i in range(seq_len):
        column = [s[i].upper() for s in sequences]
        if '-' in column or 'N' in column:
            gap_count += column.count("-")
            n_count += column.count("N")
            continue
        valid_nucleotide += 1
        bases = set(column)
        if len(bases) > 1:
            segregating_sites += 1
    
    gap_prop = gap_count / total_bases if total_bases > 0 else 0.0
    n_prop = n_count / total_bases if total_bases > 0 else 0.0
    
    if segregating_sites == 0:
        return False, "No segregating sites", seq_count, gap_prop, n_prop, valid_nucleotide, segregating_sites
    
    return True, f"{seq_count} sequences, {segregating_sites} segregating sites", seq_count, gap_prop, n_prop, valid_nucleotide, segregating_sites

def compute_haplotype_count(aln_file):
    """Compute the number of haplotypes, ignoring positions with gaps or Ns."""
    sequences = []
    with open(aln_file) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                seq = ""
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
    
    if not sequences:
        return 0
    
    seq_len = len(sequences[0])
    filtered_seqs = [list(seq.upper()) for seq in sequences]
    
    # Identify positions with gaps ('-') or Ns ('N') to ignore
    positions_to_ignore = set()
    for i in range(seq_len):
        column = [seq[i] for seq in sequences]
        if '-' in column or 'N' in column:
            positions_to_ignore.add(i)
    
    # Create haplotype strings using only valid positions
    haplotype_strings = [''.join(seq[i] for i in range(seq_len) if i not in positions_to_ignore) for seq in filtered_seqs]
    return len(set(haplotype_strings))

def run_ska_build(k, fasta_files, outdir):
    """Run ska build for given k-mer length."""
    skf_file = outdir / f"Combined-k{k}.skf"
    if skf_file.is_file() and skf_file.stat().st_size > 0:
        logging.info(f"Using existing {skf_file}")
        return skf_file
    cmd = ["ska", "build", "-k", str(k), "-o", str(skf_file)] + [str(f) for f in fasta_files]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.debug(f"ska build output for k={k}: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logging.error(f"ska build failed for k={k}: {e.stderr}")
        raise
    return skf_file

def run_ska_align(skf_file, m, outdir):
    """Run ska align for given min sample fraction."""
    m_str = str(m).replace(".", "p")
    aln_file = outdir / f"Combined-k{skf_file.stem.split('-k')[1]}-m{m_str}.aln"
    if aln_file.is_file() and aln_file.stat().st_size > 0:
        logging.info(f"Using existing {aln_file}")
        return aln_file
    cmd = ["ska", "align", str(skf_file), "-m", str(m), "--filter", "no-filter", "--ambig-mask", "-o", str(aln_file)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"ska align failed for k={skf_file.stem.split('-k')[1]}, m={m}: {e.stderr}")
        raise
    return aln_file

def run_phylogenetic_analysis(aln_file, k, m, outdir):
    """Run FastTree for phylogenetic analysis."""
    m_str = str(m).replace(".", "p")
    tree_file = outdir / f"Combined-k{k}-m{m_str}.tree"
    if tree_file.is_file() and tree_file.stat().st_size > 0:
        logging.info(f"Using existing {tree_file}")
    else:
        cmd = ["FastTree", "-gtr", "-nt", "-quote", "-quiet", "-boot", "100", str(aln_file)]
        with open(tree_file, "w") as f:
            subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
    
    try:
        tree = dendropy.Tree.get(path=tree_file, schema="newick")
        bootstrap_values = [float(node.label) for node in tree.internal_nodes() if node.label and node.label.replace(".", "").isdigit()]
        bootstrap_score = sum(bootstrap_values) / len(bootstrap_values) if bootstrap_values else 0
    except Exception as e:
        logging.error(f"Failed to parse {tree_file}: {str(e)}")
        return tree_file, 0
    return tree_file, bootstrap_score

def normalize_metrics(df):
    """Normalize metrics with 5% tolerance for segregating sites."""
    if df.empty:
        logging.warning("Empty DataFrame; returning empty")
        return df
    
    df_norm = df.copy()
    metrics = ['bootstrap_score', 'haplotype_count', 'valid_nucleotide', 'segregating', 'gap_prop']
    metrics = [m for m in metrics if m in df.columns]  # Only normalize available metrics
    
    for metric in metrics:
        min_val = df[metric].min() if not df.empty else 0
        max_val = df[metric].max() if not df.empty else 1
        if max_val == min_val:
            df_norm[f'norm_{metric}'] = 1.0
        else:
            if metric == 'segregating':
                norm = np.where(
                    df[metric] >= max_val * 0.95,
                    1.0,
                    (df[metric] - min_val) / (max_val - min_val)
                )
                df_norm[f'norm_{metric}'] = norm
            else:
                df_norm[f'norm_{metric}'] = (df[metric] - min_val) / (max_val - min_val)
            df_norm[f'norm_{metric}'] = df_norm[f'norm_{metric}'].clip(0.0, 1.0)
    
    return df_norm

def compute_composite_score(df, wb, wh, wn, ws, wg):
    """Compute composite score using normalized metrics, ensuring non-negative values."""
    if df.empty:
        logging.warning("Empty DataFrame; cannot compute composite score")
        return df
    
    df['temp_score'] = (
        wb * df.get('norm_bootstrap_score', 0.0) +
        wh * df['norm_haplotype_count'] +
        wn * df['norm_valid_nucleotide'] +
        ws * df['norm_segregating'] +
        wg * df['norm_gap_prop']
    )
    df['temp_score'] = np.clip(df['temp_score'], 0.0, None)
    return df

def plot_heatmap(df, outdir, metric='temp_score', title="Composite Score by k-mer Length and Min Sample Fraction", filename="composite_score_heatmap.png"):
    """Plot heatmap of specified metric."""
    if df.empty:
        logging.warning(f"No data to plot heatmap for {metric}")
        return
    
    all_k = sorted(df["k"].unique())
    all_m = sorted(df["m"].unique())
    pivot_score = df.pivot(index="m", columns="k", values=metric).reindex(index=all_m, columns=all_k, fill_value=0)
    pivot_haplo = df.pivot(index="m", columns="k", values="haplotype_count").reindex(index=all_m, columns=all_k, fill_value=0)
    
    plt.figure(figsize=(14, 10))
    plt.imshow(pivot_score, cmap="YlOrRd", aspect="auto", interpolation="nearest", norm=mcolors.Normalize(vmin=0, vmax=1.0))
    plt.colorbar(label=metric.replace('_', ' ').title(), ticks=np.linspace(0, 1.0, 11))
    plt.xticks(ticks=range(len(all_k)), labels=all_k, rotation=45, ha='right')
    plt.yticks(ticks=range(len(all_m)), labels=[f"{m:.1f}" for m in all_m])
    plt.xlabel("k-mer Length")
    plt.ylabel("Min Sample Fraction")
    plt.title(title)
    
    for i in range(len(all_m)):
        for j in range(len(all_k)):
            haplo_count = pivot_haplo.iloc[i, j]
            score = pivot_score.iloc[i, j]
            if not np.isnan(haplo_count):
                plt.text(j, i, f"{int(haplo_count)}\n{score:.2f}", ha="center", va="center", color="white" if score > 0.5 else "black", fontsize=10)
    
    plt.tight_layout()
    heatmap_path = outdir / filename
    plt.savefig(heatmap_path, dpi=300)
    plt.close()
    logging.info(f"Heatmap for {metric} saved as {heatmap_path}")

def plot_bootstrap_vs_haplotypes(final_df, df, outdir):
    """Plot bootstrap scores (from ska align) vs. haplotypes (from ska map) for top k,m pairs."""
    if final_df.empty or df.empty:
        logging.warning("Empty final_df or df; cannot plot bootstrap vs haplotypes")
        return
    
    # Merge final_df (ska map haplotypes) with df (ska align bootstrap) on k and m
    plot_df = final_df[['k', 'm', 'haplotype_count']].merge(
        df[['k', 'm', 'bootstrap_score']],
        on=['k', 'm'],
        how='left'
    )
    
    if plot_df.empty:
        logging.warning("No matching k,m pairs for plotting bootstrap vs haplotypes")
        return
    
    plt.figure(figsize=(12, 8))
    
    # Define markers for different m values
    markers = ['o', 's', '^', 'v', 'D', '*', 'p', 'h', '<', '>']
    unique_m = sorted(plot_df['m'].unique())
    marker_dict = {m: markers[i % len(markers)] for i, m in enumerate(unique_m)}
    
    # Plot scatter points with different markers for each m
    for m in unique_m:
        subset = plot_df[plot_df['m'] == m]
        plt.scatter(
            subset['haplotype_count'],
            subset['bootstrap_score'] * 100,
            c=subset['k'],
            cmap='viridis',
            s=100,
            alpha=0.6,
            marker=marker_dict[m],
            label=f'm={m:.1f}'
        )
    
    plt.colorbar(label='k-mer Length')
    plt.legend(title='Min Sample Fraction', bbox_to_anchor=(1.2, 1), loc='upper left')
    plt.xlabel('Number of Haplotypes (ska map)')
    plt.ylabel('Average Bootstrap Score (ska align)')
    plt.title('Bootstrap Score (ska align) vs. Haplotypes (ska map)')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.ylim(0, 100)
    
    plt.tight_layout()
    
    plot_path = outdir / 'bootstrap_vs_haplotypes.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"Scatter plot saved as {plot_path}")

def run_ska_weed(skf_file, m, outdir):
    """Run ska weed to filter k-mers."""
    m_str = str(m).replace(".", "p")
    filtered_skf = outdir / f"Filtered-k{skf_file.stem.split('-k')[1]}-m{m_str}.skf"
    if filtered_skf.is_file() and filtered_skf.stat().st_size > 0:
        logging.info(f"Using existing {filtered_skf}")
        return filtered_skf
    cmd = ["ska", "weed", str(skf_file), "-m", str(m), "-o", str(filtered_skf)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"ska weed failed for {skf_file}, m={m}: {e.stderr}")
        raise
    return filtered_skf

def run_ska_distance(skf_file, outdir):
    """Run ska distance and return output file."""
    distance_file = outdir / f"{skf_file.stem}.distance"
    if distance_file.is_file() and distance_file.stat().st_size > 0:
        logging.info(f"Using existing {distance_file}")
        return distance_file
    cmd = ["ska", "distance", str(skf_file), "-o", str(distance_file)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"ska distance failed for {skf_file}: {e.stderr}")
        raise
    return distance_file

def calculate_averages(distance_file):
    """Calculate average distances per sample."""
    value_sums = {}
    value_counts = {}
    
    with open(distance_file, 'r') as file:
        for line in file:
            if line.startswith("Sample1\tSample2\tDistance\tMismatches"):
                continue
            row = line.split()
            for i in [0, 1]:
                key = row[i]
                value = float(row[2])
                if key in value_sums:
                    value_sums[key] += value
                    value_counts[key] += 1
                else:
                    value_sums[key] = value
                    value_counts[key] = 1
    
    averages = {key: value_sums[key] / value_counts[key] for key in value_sums}
    return averages

def get_most_divergent_sample(distance_file):
    """Identify the sample with the greatest average distance."""
    averages = calculate_averages(distance_file)
    if not averages:
        logging.warning(f"No distances found in {distance_file}")
        return None
    most_divergent = max(averages, key=averages.get)
    logging.debug(f"Most divergent sample: {most_divergent} with average distance {averages[most_divergent]:.4f}")
    return most_divergent

def run_ska_map(skf_file, reference_fasta, m, outdir):
    """Run ska map with the most divergent sample as reference, including --ambig-mask."""
    m_str = str(m).replace(".", "p")
    k = skf_file.stem.split('-k')[1]
    aln_file = outdir / f"{skf_file.stem}-ska-map-{reference_fasta.stem}.aln"
    
    if aln_file.is_file() and aln_file.stat().st_size > 0:
        logging.info(f"Using existing {aln_file}")
        return aln_file
    
    # Validate input files
    if not skf_file.is_file() or skf_file.stat().st_size == 0:
        logging.error(f"SKF file {skf_file} is missing or empty")
        raise FileNotFoundError(f"SKF file {skf_file} is invalid")
    
    if not reference_fasta.is_file() or reference_fasta.stat().st_size == 0:
        logging.error(f"Reference FASTA {reference_fasta} is missing or empty")
        raise FileNotFoundError(f"Reference FASTA {reference_fasta} is invalid")
    
    cmd = ["ska", "map", "-o", str(aln_file), str(reference_fasta), str(skf_file), "--ambig-mask", "--repeat-mask"]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"ska map failed for {skf_file}, m={m}: {e.stderr}")
        raise
    return aln_file

def convert_to_nexus(aln_file, outdir):
    """Convert alignment to NEXUS format using seqret."""
    nexus_file = outdir / f"{aln_file.stem}.nex"
    if nexus_file.is_file() and nexus_file.stat().st_size > 0:
        logging.info(f"Using existing {nexus_file}")
        return nexus_file
    cmd = ["seqret", "-sequence", str(aln_file), "-outseq", str(nexus_file), "-osformat", "nexus"]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"seqret failed for {aln_file}: {e.stderr}")
        raise
    return nexus_file

def write_parameters_log(outdir, mode, k_values, m_values, wb, wh, wn, ws, wg, num_top_params):
    """Write parameters used to a log file."""
    log_path = outdir / "parameters.log"
    with open(log_path, 'w') as f:
        f.write(f"Mode: {mode}\n")
        f.write(f"k-mer sizes: {', '.join(map(str, k_values))}\n")
        f.write(f"m values: {', '.join(map(str, m_values))}\n")
        f.write(f"Weights: wb={wb}, wh={wh}, wn={wn}, ws={ws}, wg={wg}\n")
        f.write(f"Number of top k,m parameter sets: {num_top_params}\n")
    logging.info(f"Parameters logged to {log_path}")

def main():
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        print("Use --help to print program info.")
        print(ASCII_ART)
        sys.exit(0)
    
    args = parse_args()
    
    k_values = sorted(args.k)
    m_values = sorted([float(m) for m in args.m.split(',')])
    wb, wh, wn, ws, wg = args.wb, args.wh, args.wn, args.ws, args.wg
    
    mode = 'network' if args.network else 'phylo'
    logging.info(f"Running in {mode} mode")
    logging.info(f"Parameters: k={k_values}, m={m_values}, weights=[wb={wb}, wh={wh}, wn={wn}, ws={ws}, wg={wg}], num_top_params={args.num_top_params}")
    
    outdir = Path(args.output_dir)
    outdir.mkdir(exist_ok=True)
    
    # Write parameters to log file
    write_parameters_log(outdir, mode, k_values, m_values, wb, wh, wn, ws, wg, args.num_top_params)
    
    # Update FASTA headers
    try:
        fasta_files = update_fasta_headers(args.fasta_dir)
        logging.info(f"Processed {len(fasta_files)} FASTA files")
    except FileNotFoundError as e:
        logging.error(str(e))
        return
    
    # Run initial SKA analysis
    results = []
    for k in k_values:
        try:
            skf_file = run_ska_build(k, fasta_files, outdir)
            for m in m_values:
                try:
                    aln_file = run_ska_align(skf_file, m, outdir)
                    is_valid, message, seq_count, gap_prop, n_prop, valid_nucleotide, segregating = validate_alignment(aln_file)
                    logging.info(f"Alignment for k={k}, m={m}: {message}")
                    if not is_valid:
                        logging.warning(f"Skipping k={k}, m={m}: {message}")
                        continue
                    
                    haplotype_count = compute_haplotype_count(aln_file)
                    
                    # Run phylogenetic analysis in both modes
                    bootstrap_score = 0.0
                    tree_file = ""
                    try:
                        tree_file, bootstrap_score = run_phylogenetic_analysis(aln_file, k, m, outdir)
                        logging.info(f"Bootstrap score for k={k}, m={m}: {bootstrap_score:.4f}")
                    except Exception as e:
                        logging.error(f"Phylogenetic analysis failed for k={k}, m={m}: {str(e)}")
                        bootstrap_score = 0.0
                        tree_file = ""
                    
                    results.append({
                        'k': int(k),  # Ensure k is integer
                        'm': m,
                        'haplotype_count': haplotype_count,
                        'valid_nucleotide': valid_nucleotide,
                        'segregating': segregating,
                        'gap_prop': gap_prop,
                        'bootstrap_score': bootstrap_score,
                        'aln_file': str(aln_file),
                        'skf_file': str(skf_file),
                        'tree_file': str(tree_file)
                    })
                    logging.info(f"Processed k={k}, m={m}: Segregating={segregating}, Bootstrap={bootstrap_score:.4f}")
                
                except Exception as e:
                    logging.error(f"Error for k={k}, m={m}: {str(e)}")
                    continue
        except Exception as e:
            logging.error(f"Error for k={k}: {str(e)}")
            continue
    
    if not results:
        logging.error("No valid results obtained")
        return
    
    df = pd.DataFrame(results)
    df['k'] = df['k'].astype(int)  # Ensure k is integer in DataFrame
    df = normalize_metrics(df)
    df = compute_composite_score(df, wb, wh, wn, ws, wg)
    df = df.sort_values('temp_score', ascending=False)
    
    df.to_csv(outdir / 'initial_results.csv', index=False)
    logging.info("Initial results saved to initial_results.csv")
    
    # Plot heatmaps
    plot_heatmap(df, outdir, metric="temp_score", title="Composite Score by k-mer Length and Min Sample Fraction", filename="composite_score_heatmap.png")
    plot_heatmap(df, outdir, metric="bootstrap_score", title="Bootstrap Score by k-mer Length and Min Sample Fraction", filename="bootstrap_score_heatmap.png")
    
    # In phylo mode, print top 10 parameter sets and exit
    if args.phylo:
        logging.info("Phylogenetic mode: Top 10 parameter sets (sorted by composite score):")
        top_df = df.nlargest(10, 'temp_score')[['k', 'm', 'haplotype_count', 'segregating', 'bootstrap_score', 'valid_nucleotide', 'gap_prop', 'temp_score']]
        if not top_df.empty:
            print("\nTop 10 parameter sets for phylogenetic mode (sorted by composite score):")
            print(top_df.to_string(index=False, float_format=lambda x: f"{x:g}"))
        else:
            print("\nNo valid parameter sets available for phylogenetic mode.")
        return
    
    # Network mode: Process top k,m pairs with weed, distance, and map
    top_params = df.nlargest(args.num_top_params, 'temp_score')[['k', 'm', 'skf_file']]
    if len(top_params) < args.num_top_params:
        logging.warning(f"Requested {args.num_top_params} top k,m pairs, but only {len(top_params)} available; proceeding with available pairs")
        if top_params.empty and not df.empty:
            top_params = df[['k', 'm', 'skf_file']].iloc[:1]  # Fallback to highest-scoring pair
    logging.info(f"Top {len(top_params)} k,m pairs for ska map analysis: {[(row['k'], row['m']) for _, row in top_params.iterrows()]}")
    
    # Process top k,m pairs
    final_results = []
    for _, row in top_params.iterrows():
        best_k = int(row['k'])  # Ensure k is integer
        m = row['m']
        skf_file = Path(row['skf_file'])
        try:
            if not skf_file.is_file() or skf_file.stat().st_size == 0:
                logging.error(f"Unfiltered SKF file {skf_file} for k={best_k} is missing or empty")
                continue
            
            # Run ska weed
            filtered_skf = run_ska_weed(skf_file, m, outdir)
            if not filtered_skf.is_file() or filtered_skf.stat().st_size == 0:
                logging.warning(f"Filtered SKF {filtered_skf} for k={best_k}, m={m} is empty or missing")
                continue
            
            # Run ska distance
            distance_file = run_ska_distance(filtered_skf, outdir)
            
            # Calculate average distances
            most_divergent = get_most_divergent_sample(distance_file)
            if not most_divergent:
                logging.warning(f"No divergent sample found for k={best_k}, m={m}")
                continue
            logging.info(f"Most divergent sample for k={best_k}, m={m}: {most_divergent}")
            
            # Find reference FASTA
            reference_fasta = Path(args.fasta_dir) / f"{most_divergent}.fa"
            if not reference_fasta.is_file():
                logging.error(f"Reference FASTA {reference_fasta} for k={best_k}, m={m} not found")
                continue
            logging.debug(f"Using reference FASTA: {reference_fasta}")
            
            # Run ska map
            map_aln_file = run_ska_map(filtered_skf, reference_fasta, m, outdir)
            
            # Convert to NEXUS
            nexus_file = convert_to_nexus(map_aln_file, outdir)
            
            # Calculate metrics
            is_valid, message, _, gap_prop, _, valid_nucleotide, segregating = validate_alignment(map_aln_file)
            if not is_valid:
                logging.warning(f"Invalid mapped alignment for k={best_k}, m={m}: {message}")
                continue
            
            haplotype_count = compute_haplotype_count(map_aln_file)
            
            # Get bootstrap score from initial results
            bootstrap_score = df[(df['k'] == best_k) & (df['m'] == m)]['bootstrap_score'].iloc[0] if not df[(df['k'] == best_k) & (df['m'] == m)].empty else 0.0
            tree_file = df[(df['k'] == best_k) & (df['m'] == m)]['tree_file'].iloc[0] if not df[(df['k'] == best_k) & (df['m'] == m)].empty else ""
            
            final_results.append({
                'k': best_k,
                'm': m,
                'bootstrap_score': bootstrap_score,
                'most_divergent': most_divergent,
                'haplotype_count': haplotype_count,
                'valid_nucleotide': valid_nucleotide,
                'segregating': segregating,
                'gap_prop': gap_prop,
                'aln_file': str(map_aln_file),
                'tree_file': str(tree_file),
                'nexus_file': str(nexus_file)
            })
            logging.info(f"Mapped k={best_k}, m={m}: Bootstrap={bootstrap_score:.2f}, Hap={haplotype_count}, Valid={valid_nucleotide}, Seg={segregating}")
            
        except Exception as e:
            logging.error(f"Error processing k={best_k}, m={m}: {str(e)}")
            continue
    
    # Save and log final results
    final_df = pd.DataFrame(final_results)
    if not final_df.empty:
        final_df = normalize_metrics(final_df)
        final_df = compute_composite_score(final_df, wb, wh, wn, ws, wg)
        # Sort by haplotype_count, segregating (SNPs), and bootstrap_score (descending)
        final_df = final_df.sort_values(by=['haplotype_count', 'segregating', 'bootstrap_score'], ascending=[False, False, False])
        
        final_df.to_csv(outdir / 'final_results.csv', index=False)
        logging.info("Final results saved to final_results.csv, ranked by haplotype count (ska map), SNPs (segregating sites, ska map), and bootstrap score (ska align)")
        logging.info("Final results for top k,m pairs (sorted by haplotype count, SNPs, then bootstrap score):")
        logging.info(final_df[['k', 'm', 'haplotype_count', 'segregating', 'bootstrap_score', 'most_divergent', 'valid_nucleotide']].to_string(index=False))
        
        # Plot bootstrap vs. haplotypes using map haplotypes and align bootstrap
        plot_bootstrap_vs_haplotypes(final_df, df, outdir)
    else:
        logging.warning("No valid final results generated")

if __name__ == "__main__":
    main()