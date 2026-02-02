#!/usr/bin/env python3
"""
Generate Per-Dataset Performance Plots for scRNA-seq

Shows how each system performs on individual datasets.
Perfect for showing system robustness across different scRNA data types:
- Cell lines (REH, SUP-B15)
- PBMC (healthy human)
- Mouse brain (10k cells)
- Cell cycle datasets (GSE75748, GSE146773)
- H1/FUCCI normalized (GSE64016)

Usage:
    python3 bin/generate_per_dataset_plots_scrna.py scrna_evaluation.csv [output_prefix]
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ============================================================================
# Configuration
# ============================================================================

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 9

# Colorblind-friendly palette
COLORS = {
    'parallel': '#2166AC',
    'sequential': '#4393C3',
    'single_claude': '#D6604D',
    'single_gpt': '#F4A582',
    'single_gemini': '#FDB863',
    'baseline': '#999999'
}

# Dataset display names (short for x-axis)
DATASET_LABELS = {
    '1_GD428_21136_Hu_REH_Parental': 'REH\n(Parental)',
    '2_GD444_21136_Hu_SUP_Parental': 'SUP-B15\n(Parental)',
    '2_GD444_21136_Hu_Sup_Parental': 'SUP-B15\n(Parental)',
    'pbmc_healthy_human': 'PBMC\n(Human)',
    '10-k-brain-cells_healthy_mouse': 'Mouse Brain\n(10k cells)',
    'GSE75748': 'GSE75748\n(hESC)',
    'GSE146773': 'GSE146773\n(U-2 OS)',
    '3_GSE64016_H1andFUCCI_normalized_EC_original': 'GSE64016\n(H1/FUCCI)'
}

def parse_system_config(system_code):
    """Parse system configuration with full role extraction"""
    if 'no_agent' in system_code or system_code == 'no_agent':
        return {'type': 'no_agent', 'arch': None, 'roles': None}
    elif 'single_claude' in system_code or system_code == 'single_claude':
        return {'type': 'single_claude', 'arch': None, 'roles': None}
    elif 'single_gpt' in system_code or system_code == 'single_gpt':
        return {'type': 'single_gpt', 'arch': None, 'roles': None}
    elif 'single_gemini' in system_code or system_code == 'single_gemini':
        return {'type': 'single_gemini', 'arch': None, 'roles': None}

    arch = 'Parallel' if 'parallel' in system_code.lower() else 'Sequential'

    # Extract full role configuration
    # Default configuration: GPT-5.2=stats, Claude=pipeline, Gemini=biology
    if '_default' in system_code or 'default' == system_code:
        gpt_role = 'stats'
        claude_role = 'pipe'
        gemini_role = 'bio'
    else:
        gpt_role = None
        claude_role = None
        gemini_role = None

        if 'gp_st' in system_code:
            gpt_role = 'stats'
        elif 'gp_pl' in system_code:
            gpt_role = 'pipe'
        elif 'gp_bl' in system_code:
            gpt_role = 'bio'

        if 'cl_st' in system_code:
            claude_role = 'stats'
        elif 'cl_pl' in system_code:
            claude_role = 'pipe'
        elif 'cl_bl' in system_code:
            claude_role = 'bio'

        if 'gm_st' in system_code:
            gemini_role = 'stats'
        elif 'gm_pl' in system_code:
            gemini_role = 'pipe'
        elif 'gm_bl' in system_code:
            gemini_role = 'bio'

    roles_str = f"GPT={gpt_role or '?'}, Cl={claude_role or '?'}, Gem={gemini_role or '?'}"

    return {'type': 'multi_agent', 'arch': arch, 'roles': roles_str,
            'gpt': gpt_role, 'claude': claude_role, 'gemini': gemini_role}

def get_system_label(system_code, short=False):
    """Get display label for system"""
    parsed = parse_system_config(system_code)

    if parsed['type'] == 'no_agent':
        return 'No-Agent'
    elif parsed['type'] == 'single_claude':
        return 'Single-Claude' if short else 'Single-LLM (Claude)'
    elif parsed['type'] == 'single_gpt':
        return 'Single-GPT' if short else 'Single-LLM (GPT-5.2)'
    elif parsed['type'] == 'single_gemini':
        return 'Single-Gemini' if short else 'Single-LLM (Gemini)'
    elif parsed['type'] == 'multi_agent':
        arch_abbr = 'Par.' if parsed['arch'] == 'Parallel' else 'Seq.'
        if short:
            return f"PolyLLM ({arch_abbr})"
        return f"PolyLLM ({arch_abbr}) {parsed['roles']}"
    return system_code

def assign_color(system_code):
    """Assign color to system"""
    parsed = parse_system_config(system_code)

    if parsed['type'] == 'no_agent':
        return COLORS['baseline']
    elif parsed['type'] == 'single_claude':
        return COLORS['single_claude']
    elif parsed['type'] == 'single_gpt':
        return COLORS['single_gpt']
    elif parsed['type'] == 'single_gemini':
        return COLORS['single_gemini']
    elif parsed['type'] == 'multi_agent':
        return COLORS['sequential'] if parsed['arch'] == 'Sequential' else COLORS['parallel']
    return COLORS['baseline']

def select_top_systems(df, for_heatmap=False):
    """Select top 1 parallel + top 1 sequential PolyLLM + single-LLM baselines"""
    parallel_systems = []
    sequential_systems = []

    for system in df['System'].unique():
        parsed = parse_system_config(system)
        if parsed['type'] == 'multi_agent':
            if parsed['arch'] == 'Parallel':
                parallel_systems.append(system)
            else:
                sequential_systems.append(system)

    # Calculate mean accuracy for parallel systems
    parallel_accuracy = []
    for system in parallel_systems:
        system_df = df[df['System'] == system]
        mean_acc = system_df['Accuracy_%'].mean()
        parallel_accuracy.append({'System': system, 'Mean_Accuracy': mean_acc})

    # Calculate mean accuracy for sequential systems
    sequential_accuracy = []
    for system in sequential_systems:
        system_df = df[df['System'] == system]
        mean_acc = system_df['Accuracy_%'].mean()
        sequential_accuracy.append({'System': system, 'Mean_Accuracy': mean_acc})

    # Select top 1 from each
    top_polyllm = []
    if parallel_accuracy:
        parallel_df = pd.DataFrame(parallel_accuracy).sort_values('Mean_Accuracy', ascending=False)
        top_polyllm.append(parallel_df.head(1)['System'].values[0])

    if sequential_accuracy:
        sequential_df = pd.DataFrame(sequential_accuracy).sort_values('Mean_Accuracy', ascending=False)
        top_polyllm.append(sequential_df.head(1)['System'].values[0])

    # Add single-LLM baselines
    fixed_systems = []
    for system in df['System'].unique():
        parsed = parse_system_config(system)
        if parsed['type'] in ['single_gpt', 'single_claude', 'single_gemini']:
            fixed_systems.append(system)

    return top_polyllm + fixed_systems

# ============================================================================
# Plotting Functions
# ============================================================================

def plot_accuracy_per_dataset(df, systems, output_prefix):
    """Grouped bar chart: Accuracy per dataset"""
    datasets = sorted(df['Dataset'].unique())

    # Compact size for single-column layout
    fig, ax = plt.subplots(figsize=(7, 4))

    x = np.arange(len(datasets))
    width = 0.10  # Narrower bars

    for i, system in enumerate(systems):
        system_data = []
        for dataset in datasets:
            subset = df[(df['System'] == system) & (df['Dataset'] == dataset)]
            if len(subset) > 0:
                system_data.append(subset['Accuracy_%'].values[0])
            else:
                system_data.append(0)

        color = assign_color(system)
        label = get_system_label(system, short=True)
        offset = (i - len(systems)/2) * width

        bars = ax.bar(x + offset, system_data, width, label=label,
                     color=color, edgecolor='black', linewidth=0.8)

    ax.set_xlabel('Dataset', fontsize=14, fontweight='bold')
    ax.set_ylabel('Decision Accuracy (%)', fontsize=14, fontweight='bold')
    ax.set_title('Per-Dataset Decision Accuracy (scRNA-seq)', fontsize=15, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([DATASET_LABELS.get(d, d) for d in datasets], fontsize=9)
    ax.set_ylim(0, 105)
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), ncol=1, fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_per_dataset_accuracy.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_per_dataset_accuracy.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_per_dataset_accuracy.pdf / .png")

def plot_completion_per_dataset(df, systems, output_prefix):
    """Grouped bar chart: Completion rate per dataset"""
    datasets = sorted(df['Dataset'].unique())

    # Compact size for single-column layout
    fig, ax = plt.subplots(figsize=(7, 4))

    x = np.arange(len(datasets))
    width = 0.10  # Narrower bars

    for i, system in enumerate(systems):
        system_data = []
        for dataset in datasets:
            subset = df[(df['System'] == system) & (df['Dataset'] == dataset)]
            if len(subset) > 0:
                # Completion = reached stage 5
                completed = 1 if subset['Completed'].values[0] == 1 else 0
                system_data.append(completed * 100)
            else:
                system_data.append(0)

        color = assign_color(system)
        label = get_system_label(system, short=True)
        offset = (i - len(systems)/2) * width

        bars = ax.bar(x + offset, system_data, width, label=label,
                     color=color, edgecolor='black', linewidth=0.8)

    ax.set_xlabel('Dataset', fontsize=14, fontweight='bold')
    ax.set_ylabel('Completion Rate (%)', fontsize=14, fontweight='bold')
    ax.set_title('Per-Dataset Pipeline Completion Rate (scRNA-seq)', fontsize=15, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([DATASET_LABELS.get(d, d) for d in datasets], fontsize=9)
    ax.set_ylim(0, 105)
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), ncol=1, fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_per_dataset_completion.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_per_dataset_completion.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_per_dataset_completion.pdf / .png")

def plot_heatmap_accuracy(df, systems, output_prefix):
    """Heatmap: Systems × Datasets accuracy (excluding systems with multiple 0% datasets)"""
    datasets = sorted(df['Dataset'].unique())

    # Build matrix and filter out systems with multiple 0% datasets
    matrix = []
    row_labels = []
    filtered_systems = []

    for system in systems:
        row = []
        zero_count = 0
        for dataset in datasets:
            subset = df[(df['System'] == system) & (df['Dataset'] == dataset)]
            if len(subset) > 0:
                acc = subset['Accuracy_%'].values[0]
                row.append(acc)
                if acc == 0:
                    zero_count += 1
            else:
                row.append(np.nan)

        # Only filter systems with more than 1 dataset at 0%
        if zero_count <= 1:
            matrix.append(row)
            row_labels.append(get_system_label(system, short=False))
            filtered_systems.append(system)

    matrix = np.array(matrix)

    if len(filtered_systems) > 0:
        print(f"   → Filtered out {len(systems) - len(filtered_systems)} systems with >1 dataset at 0%")
        print(f"   → Showing {len(filtered_systems)} systems in heatmap")
    else:
        print(f"   ⚠️  All systems filtered out. Showing all systems instead.")
        # If all filtered, show all systems
        matrix = []
        row_labels = []
        filtered_systems = systems
        for system in systems:
            row = []
            for dataset in datasets:
                subset = df[(df['System'] == system) & (df['Dataset'] == dataset)]
                if len(subset) > 0:
                    row.append(subset['Accuracy_%'].values[0])
                else:
                    row.append(np.nan)
            matrix.append(row)
            row_labels.append(get_system_label(system, short=False))
        matrix = np.array(matrix)

    fig, ax = plt.subplots(figsize=(14, 6))

    im = ax.imshow(matrix, cmap='viridis', aspect='auto', vmin=0, vmax=100)

    ax.set_xticks(np.arange(len(datasets)))
    ax.set_yticks(np.arange(len(filtered_systems)))
    ax.set_xticklabels([DATASET_LABELS.get(d, d) for d in datasets], fontsize=9)
    ax.set_yticklabels(row_labels, fontsize=10)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Add text annotations
    for i in range(len(filtered_systems)):
        for j in range(len(datasets)):
            if not np.isnan(matrix[i, j]):
                text = ax.text(j, i, f'{matrix[i, j]:.0f}',
                             ha="center", va="center", color="black", fontsize=9, fontweight='bold')

    ax.set_title('Decision Accuracy Heatmap: Systems × Datasets (scRNA-seq)', fontsize=14, fontweight='bold', pad=20)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Accuracy (%)', rotation=270, labelpad=20, fontsize=12)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_heatmap_accuracy.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_heatmap_accuracy.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_heatmap_accuracy.pdf / .png")

# ============================================================================
# Main
# ============================================================================

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 bin/generate_per_dataset_plots_scrna.py scrna_evaluation.csv [output_prefix]")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else 'experiments/scrna_csv_figures/scrna_per_dataset_figure'

    # Ensure output directory exists
    output_dir = os.path.dirname(output_prefix) if '/' in output_prefix else '.'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"📁 Created output directory: {output_dir}\n")

    print("\n" + "="*80)
    print("GENERATING PER-DATASET PERFORMANCE PLOTS (scRNA-seq)")
    print("="*80 + "\n")

    print("📖 Loading evaluation data...")
    df = pd.read_csv(csv_file)

    # Rename columns to match expected format
    if 'Total_Stages' in df.columns:
        df = df.rename(columns={'Total_Stages': 'Total_Decisions'})

    print(f"   ✓ Loaded {len(df)} experiments\n")

    # Select systems
    systems = select_top_systems(df)

    print(f"📊 Plotting {len(systems)} systems across {len(df['Dataset'].unique())} datasets\n")

    print("📊 Generating Figure 1: Accuracy per Dataset (Grouped Bar)...")
    plot_accuracy_per_dataset(df, systems, output_prefix)

    print("📊 Generating Figure 2: Completion per Dataset (Grouped Bar)...")
    plot_completion_per_dataset(df, systems, output_prefix)

    print("📊 Generating Figure 3: Accuracy Heatmap (Systems × Datasets)...")
    plot_heatmap_accuracy(df, systems, output_prefix)

    print("\n" + "="*80)
    print("✨ All per-dataset plots generated successfully!")
    print("="*80)
    print(f"\nOutput files:")
    print(f"  • {output_prefix}_per_dataset_accuracy.pdf / .png")
    print(f"  • {output_prefix}_per_dataset_completion.pdf / .png")
    print(f"  • {output_prefix}_heatmap_accuracy.pdf / .png")
    print("\nThese plots show how each system performs on INDIVIDUAL datasets")
    print("Perfect for demonstrating robustness across different scRNA-seq data types!\n")

if __name__ == '__main__':
    main()
