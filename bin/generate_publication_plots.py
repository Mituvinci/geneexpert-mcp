#!/usr/bin/env python3
"""
Generate Publication-Quality Plots for ICML 2026 Paper

Creates 4 publication-ready figures:
1. Decision Accuracy Comparison
2. Pipeline Completion Rate
3. Total Cost Per Experiment
4. User Intervention Rate

Style: ICML conference format (Times font, 14-18pt, colorblind-safe)

Usage:
    python3 bin/generate_publication_plots.py bulk_evaluation_per_experiment.csv experiments/results
"""

import sys
import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy import stats

# ============================================================================
# ICML Publication Style Configuration
# ============================================================================

# Font configuration (ICML standard)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14

# Colorblind-friendly palette
COLORS = {
    'multi_foundational_1': '#2166AC',  # Dark blue
    'multi_foundational_2': '#4393C3',  # Light blue
    'single_llm_1': '#D6604D',          # Dark orange
    'single_llm_2': '#F4A582',          # Medium orange
    'single_llm_3': '#FDB863',          # Light orange
    'baseline': '#999999'               # Gray
}

# System display names (shortened for publication)
SYSTEM_LABELS = {
    'role_perm_3': 'Multi-Agent (GPT-pipe)',
    'default_sequential': 'Multi-Agent (Sequential)',
    'single_agent_gpt': 'Single-LLM (GPT-5.2)',
    'single_agent_gemini': 'Single-LLM (Gemini)',
    'single_agent_claude': 'Single-LLM (Claude)',
    'no_agent': 'No-Agent Baseline'
}

# ============================================================================
# Data Loading Functions
# ============================================================================

def load_evaluation_data(csv_file):
    """Load per-experiment evaluation CSV"""
    df = pd.read_csv(csv_file)
    return df

def extract_costs_from_experiments(results_dir):
    """Extract total costs from all experiment folders"""
    costs = {}

    if not os.path.exists(results_dir):
        print(f"⚠ Warning: Results directory not found: {results_dir}")
        return costs

    for folder in os.listdir(results_dir):
        folder_path = os.path.join(results_dir, folder)
        if not os.path.isdir(folder_path):
            continue

        jsonl_file = os.path.join(folder_path, 'staged_analysis_agent_decisions.jsonl')

        if not os.path.exists(jsonl_file):
            continue

        total_cost = 0.0
        try:
            with open(jsonl_file, 'r') as f:
                for line in f:
                    if not line.strip():
                        continue
                    entry = json.loads(line)
                    if 'costs' in entry and 'total_usd' in entry['costs']:
                        total_cost += entry['costs']['total_usd']
        except Exception as e:
            print(f"⚠ Warning: Error reading {jsonl_file}: {e}")
            continue

        costs[folder] = total_cost

    return costs

def prepare_plot_data(df, costs_dict, top_systems):
    """Prepare aggregated data for plotting"""

    plot_data = []

    for system in top_systems:
        system_df = df[df['System'] == system]

        if len(system_df) == 0:
            continue

        # Accuracy
        accuracy_vals = system_df['Accuracy_%'].values
        accuracy_mean = np.mean(accuracy_vals)
        accuracy_sem = stats.sem(accuracy_vals) if len(accuracy_vals) > 1 else 0

        # Completion rate
        completion_vals = system_df['Completed'].values
        completion_rate = (np.sum(completion_vals) / len(completion_vals)) * 100

        # User interventions
        intervention_vals = system_df['User_Interventions'].values
        intervention_rate = (np.sum(intervention_vals) / len(system_df)) * 100

        # Cost (extract from folders)
        experiment_names = system_df['Experiment'].values
        experiment_costs = [costs_dict.get(exp, 0.0) for exp in experiment_names]
        cost_mean = np.mean(experiment_costs) if experiment_costs else 0.0
        cost_sem = stats.sem(experiment_costs) if len(experiment_costs) > 1 else 0.0

        plot_data.append({
            'system': system,
            'label': SYSTEM_LABELS.get(system, system),
            'accuracy_mean': accuracy_mean,
            'accuracy_sem': accuracy_sem,
            'completion_rate': completion_rate,
            'intervention_rate': intervention_rate,
            'cost_mean': cost_mean,
            'cost_sem': cost_sem
        })

    return pd.DataFrame(plot_data)

# ============================================================================
# Statistical Significance Testing
# ============================================================================

def compute_significance(df, baseline_system, metric_column):
    """
    Compute p-values comparing each system to baseline.
    Returns dict of {system: p_value}
    """
    significance = {}

    baseline_vals = df[df['System'] == baseline_system][metric_column].values

    for system in df['System'].unique():
        if system == baseline_system:
            continue

        system_vals = df[df['System'] == system][metric_column].values

        if len(system_vals) > 1 and len(baseline_vals) > 1:
            _, p_value = stats.ttest_ind(system_vals, baseline_vals)
            significance[system] = p_value
        else:
            significance[system] = 1.0  # Not significant

    return significance

def get_significance_marker(p_value):
    """Convert p-value to marker (*, **, ***)"""
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''

# ============================================================================
# Plotting Functions
# ============================================================================

def assign_colors(systems):
    """Assign colors to systems based on type"""
    color_map = {}

    multi_idx = 0
    single_idx = 0

    for system in systems:
        if system == 'no_agent':
            color_map[system] = COLORS['baseline']
        elif 'single' in system:
            single_idx += 1
            color_map[system] = COLORS[f'single_llm_{single_idx}']
        else:
            multi_idx += 1
            color_map[system] = COLORS[f'multi_foundational_{multi_idx}']

    return color_map

def plot_accuracy(plot_df, df_raw, baseline_system, output_prefix):
    """Figure 1: Decision Accuracy Comparison"""

    fig, ax = plt.subplots(figsize=(8, 5))

    # Compute significance
    significance = compute_significance(df_raw, baseline_system, 'Accuracy_%')

    # Assign colors
    colors = assign_colors(plot_df['system'].values)
    bar_colors = [colors[s] for s in plot_df['system'].values]

    # Horizontal bar chart
    y_pos = np.arange(len(plot_df))
    bars = ax.barh(y_pos, plot_df['accuracy_mean'], xerr=plot_df['accuracy_sem'],
                   color=bar_colors, edgecolor='black', linewidth=1.5,
                   error_kw={'linewidth': 2, 'ecolor': 'black'})

    # Add significance markers
    for i, (system, mean_val) in enumerate(zip(plot_df['system'], plot_df['accuracy_mean'])):
        if system in significance:
            marker = get_significance_marker(significance[system])
            if marker:
                ax.text(mean_val + 2, i, marker, va='center', fontsize=16, fontweight='bold')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'])
    ax.set_xlabel('Decision Accuracy (%)', fontsize=18)
    ax.set_xlim(0, 100)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_accuracy.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_accuracy.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_accuracy.pdf / .png")

def plot_completion(plot_df, output_prefix):
    """Figure 2: Pipeline Completion Rate"""

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = assign_colors(plot_df['system'].values)
    bar_colors = [colors[s] for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    ax.barh(y_pos, plot_df['completion_rate'],
            color=bar_colors, edgecolor='black', linewidth=1.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'])
    ax.set_xlabel('Completion Rate (%)', fontsize=18)
    ax.set_xlim(0, 100)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_completion.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_completion.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_completion.pdf / .png")

def plot_cost(plot_df, df_raw, baseline_system, output_prefix):
    """Figure 3: Total Cost Per Experiment"""

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = assign_colors(plot_df['system'].values)
    bar_colors = [colors[s] for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    ax.barh(y_pos, plot_df['cost_mean'], xerr=plot_df['cost_sem'],
            color=bar_colors, edgecolor='black', linewidth=1.5,
            error_kw={'linewidth': 2, 'ecolor': 'black'})

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'])
    ax.set_xlabel('Total Cost per Experiment (USD)', fontsize=18)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_cost.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_cost.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_cost.pdf / .png")

def plot_user_intervention(plot_df, output_prefix):
    """Figure 4: User Intervention Rate"""

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = assign_colors(plot_df['system'].values)
    bar_colors = [colors[s] for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    ax.barh(y_pos, plot_df['intervention_rate'],
            color=bar_colors, edgecolor='black', linewidth=1.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'])
    ax.set_xlabel('User Intervention Rate (%)', fontsize=18)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_user_intervention.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_user_intervention.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_user_intervention.pdf / .png")

# ============================================================================
# Main Execution
# ============================================================================

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 bin/generate_publication_plots.py bulk_evaluation_per_experiment.csv experiments/results [output_prefix]")
        sys.exit(1)

    csv_file = sys.argv[1]
    results_dir = sys.argv[2]
    output_prefix = sys.argv[3] if len(sys.argv) > 3 else 'figure'

    print("\n" + "="*80)
    print("GENERATING PUBLICATION-QUALITY PLOTS (ICML 2026 Format)")
    print("="*80 + "\n")

    # Load data
    print("📖 Loading evaluation data...")
    df = load_evaluation_data(csv_file)
    print(f"   ✓ Loaded {len(df)} experiments\n")

    print("💰 Extracting cost data from experiment folders...")
    costs = extract_costs_from_experiments(results_dir)
    print(f"   ✓ Extracted costs from {len(costs)} experiments\n")

    # Select top 5 systems + baseline
    print("🎯 Selecting top 5 architectures + baseline...")
    top_systems = [
        'role_perm_3',          # Multi-foundational (best)
        'default_sequential',   # Multi-foundational (sequential)
        'single_agent_gpt',     # Single-LLM (best)
        'single_agent_gemini',  # Single-LLM
        'single_agent_claude',  # Single-LLM
        'no_agent'              # Baseline
    ]

    # Prepare aggregated data
    plot_df = prepare_plot_data(df, costs, top_systems)

    print("\nSelected Systems:")
    for _, row in plot_df.iterrows():
        print(f"  • {row['label']}")
    print()

    # Generate plots
    baseline_system = 'no_agent'

    print("📊 Generating Figure 1: Decision Accuracy...")
    plot_accuracy(plot_df, df, baseline_system, output_prefix)

    print("📊 Generating Figure 2: Completion Rate...")
    plot_completion(plot_df, output_prefix)

    print("📊 Generating Figure 3: Total Cost...")
    plot_cost(plot_df, df, baseline_system, output_prefix)

    print("📊 Generating Figure 4: User Intervention Rate...")
    plot_user_intervention(plot_df, output_prefix)

    print("\n" + "="*80)
    print("✨ All plots generated successfully!")
    print("="*80)
    print("\nOutput files:")
    print(f"  • {output_prefix}_accuracy.pdf / .png")
    print(f"  • {output_prefix}_completion.pdf / .png")
    print(f"  • {output_prefix}_cost.pdf / .png")
    print(f"  • {output_prefix}_user_intervention.pdf / .png")
    print("\nStatistical significance markers:")
    print("  * p < 0.05, ** p < 0.01, *** p < 0.001")
    print()

if __name__ == '__main__':
    main()
