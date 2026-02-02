#!/usr/bin/env python3
"""
Generate G1: Per-Dataset Accuracy Heatmap
Shows validation accuracy (%) for 10 system configurations across all datasets
Outputs: Heatmap (systems × datasets) with color-coded accuracy
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
import os

def parse_system_name(system):
    """Convert system name to human-readable label"""
    if system == 'no_agent':
        return 'No-Agent\nBaseline'
    elif system == 'single_gpt':
        return 'Single-LLM\nGPT-5.2'
    elif system == 'single_claude':
        return 'Single-LLM\nClaude'
    elif system == 'single_gemini':
        return 'Single-LLM\nGemini'
    elif 'parallel_gp_st_cl_pl_gm_bl' in system:
        return 'PolyLLM\nGPT=Stats\nClaude=Pipeline\nGemini=Biology'
    elif 'parallel_gp_st_cl_bl_gm_pl' in system:
        return 'PolyLLM\nGPT=Stats\nClaude=Biology\nGemini=Pipeline'
    elif 'parallel_gp_pl_cl_st_gm_bl' in system:
        return 'PolyLLM\nGPT=Pipeline\nClaude=Stats\nGemini=Biology'
    elif 'parallel_gp_pl_cl_bl_gm_st' in system:
        return 'PolyLLM\nGPT=Pipeline\nClaude=Biology\nGemini=Stats'
    elif 'parallel_gp_bl_cl_st_gm_pl' in system:
        return 'PolyLLM\nGPT=Biology\nClaude=Stats\nGemini=Pipeline'
    elif 'parallel_gp_bl_cl_pl_gm_st' in system:
        return 'PolyLLM\nGPT=Biology\nClaude=Pipeline\nGemini=Stats'
    elif 'sequential_gp_st_cl_pl_gm_bl' in system:
        return 'PolyLLM-Seq\nGPT=Stats\nClaude=Pipeline\nGemini=Biology'
    elif 'sequential_gp_st_cl_bl_gm_pl' in system:
        return 'PolyLLM-Seq\nGPT=Stats\nClaude=Biology\nGemini=Pipeline'
    elif 'sequential_gp_pl_cl_st_gm_bl' in system:
        return 'PolyLLM-Seq\nGPT=Pipeline\nClaude=Stats\nGemini=Biology'
    elif 'sequential_gp_pl_cl_bl_gm_st' in system:
        return 'PolyLLM-Seq\nGPT=Pipeline\nClaude=Biology\nGemini=Stats'
    elif 'sequential_gp_bl_cl_st_gm_pl' in system:
        return 'PolyLLM-Seq\nGPT=Biology\nClaude=Stats\nGemini=Pipeline'
    elif 'sequential_gp_bl_cl_pl_gm_st' in system:
        return 'PolyLLM-Seq\nGPT=Biology\nClaude=Pipeline\nGemini=Stats'
    else:
        return system.replace('_', '\n').title()

def clean_dataset_name(dataset):
    """Convert dataset name to human-readable label"""
    # Remove numeric prefixes and simplify
    dataset = dataset.replace('1_', '').replace('2_', '').replace('3_', '').replace('4_', '').replace('5_', '').replace('6_', '').replace('7_', '')
    dataset = dataset.replace('_pe_clean', '').replace('_se_clean', '')
    dataset = dataset.replace('_batch', '\n(Batch)')
    dataset = dataset.replace('_contam', '\n(Contam)')
    dataset = dataset.replace('10-k-brain-cells_healthy_mouse', '10k Mouse Brain')
    dataset = dataset.replace('Hu_REH_Parental', 'REH')
    dataset = dataset.replace('Hu_SUP_B15_Parental', 'SUP-B15')
    dataset = dataset.replace('pbmc_healthy_human', 'PBMC Human')
    return dataset

def generate_heatmap(csv_file, output_dir, modality):
    """Generate per-dataset accuracy heatmap for 10 system configurations"""

    print(f"Reading {csv_file}...")
    df = pd.read_csv(csv_file)

    # Define the 10 systems we want (in order)
    # 6 role permutations (parallel) + 3 single-LLM + 1 no-agent
    system_order = [
        'parallel_gp_st_cl_pl_gm_bl',  # 1. GPT=stats, Claude=pipeline, Gemini=biology
        'parallel_gp_st_cl_bl_gm_pl',  # 2. GPT=stats, Claude=biology, Gemini=pipeline
        'parallel_gp_pl_cl_st_gm_bl',  # 3. GPT=pipeline, Claude=stats, Gemini=biology
        'parallel_gp_pl_cl_bl_gm_st',  # 4. GPT=pipeline, Claude=biology, Gemini=stats
        'parallel_gp_bl_cl_st_gm_pl',  # 5. GPT=biology, Claude=stats, Gemini=pipeline
        'parallel_gp_bl_cl_pl_gm_st',  # 6. GPT=biology, Claude=pipeline, Gemini=stats
        'single_gpt',                  # 7. Single-LLM GPT-5.2
        'single_claude',               # 8. Single-LLM Claude
        'single_gemini',               # 9. Single-LLM Gemini
        'no_agent'                     # 10. No-agent baseline
    ]

    # Filter to only include these systems
    df = df[df['System'].isin(system_order)]

    # Create pivot table: rows=systems, columns=datasets, values=accuracy
    pivot_data = df.pivot_table(
        index='System',
        columns='Dataset',
        values='Accuracy_%',
        aggfunc='mean'  # In case of duplicates, take mean
    )

    # Reorder rows according to system_order
    pivot_data = pivot_data.reindex([s for s in system_order if s in pivot_data.index])

    # Clean dataset names for column labels
    pivot_data.columns = [clean_dataset_name(col) for col in pivot_data.columns]

    # Clean system names for row labels
    pivot_data.index = [parse_system_name(idx) for idx in pivot_data.index]

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10))

    # Create heatmap with red-yellow-green colormap
    sns.heatmap(
        pivot_data,
        annot=True,
        fmt='.1f',
        cmap='RdYlGn',
        vmin=0,
        vmax=100,
        cbar_kws={'label': 'Accuracy (%)', 'shrink': 0.8},
        linewidths=0.5,
        linecolor='gray',
        ax=ax
    )

    # Formatting
    title = f'Per-Dataset Accuracy Heatmap ({modality})\n10 System Configurations Across Datasets'
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Dataset', fontsize=12, fontweight='bold')
    ax.set_ylabel('System Configuration', fontsize=12, fontweight='bold')

    # Rotate labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)

    plt.tight_layout()

    # Save figure
    output_prefix = f'G1_per_dataset_heatmap_{modality.lower().replace(" ", "_")}'
    png_path = os.path.join(output_dir, f'{output_prefix}.png')
    pdf_path = os.path.join(output_dir, f'{output_prefix}.pdf')

    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")

    plt.close()

    # Print summary statistics
    print(f"\n{modality} Summary:")
    print(f"  Total systems: {len(pivot_data)}")
    print(f"  Total datasets: {len(pivot_data.columns)}")
    print(f"  Mean accuracy: {pivot_data.values.flatten().mean():.1f}%")
    print(f"  PolyLLM configs (rows 1-6) mean: {pivot_data.iloc[:6].values.flatten().mean():.1f}%")
    print(f"  Single-LLM (rows 7-9) mean: {pivot_data.iloc[6:9].values.flatten().mean():.1f}%")
    if len(pivot_data) >= 10:
        print(f"  No-agent (row 10) mean: {pivot_data.iloc[9].values.mean():.1f}%")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python generate_G1_per_dataset_heatmap.py <csv_file>")
        print("\nExample:")
        print("  python bin/generate_G1_per_dataset_heatmap.py experiments/bulk_rna_csv_figures/bulk_evaluation_per_experiment.csv")
        print("  python bin/generate_G1_per_dataset_heatmap.py experiments/scrna_csv_figures/scrna_evaluation_per_experiment.csv")
        sys.exit(1)

    csv_file = sys.argv[1]

    if not os.path.exists(csv_file):
        print(f"Error: File not found: {csv_file}")
        sys.exit(1)

    # Determine output directory and modality from input path
    if 'bulk' in csv_file.lower():
        output_dir = 'experiments/bulk_rna_csv_figures'
        modality = 'Bulk RNA-seq'
    elif 'scrna' in csv_file.lower():
        output_dir = 'experiments/scrna_csv_figures'
        modality = 'scRNA-seq'
    else:
        output_dir = os.path.dirname(csv_file)
        modality = 'RNA-seq'

    os.makedirs(output_dir, exist_ok=True)

    generate_heatmap(csv_file, output_dir, modality)
    print("\nDone!")
