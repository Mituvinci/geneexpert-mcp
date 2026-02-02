#!/usr/bin/env python3
"""
Generate Publication-Quality Plots from CSV Files (ICML 2026)

Uses pre-generated CSV files for faster processing.

Usage:
    python3 bin/generate_plots_from_csv.py bulk_evaluation_v2_per_experiment.csv [output_prefix]
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# ============================================================================
# ICML Publication Style Configuration
# ============================================================================

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
    'multi_foundational_1': '#2166AC',
    'multi_foundational_2': '#4393C3',
    'single_llm_1': '#D6604D',
    'single_llm_2': '#F4A582',
    'single_llm_3': '#FDB863',
    'baseline': '#999999'
}

# System display names
SYSTEM_LABELS = {
    'role_perm_3': 'Multi-Agent (GPT-pipe)',
    'default_sequential': 'Multi-Agent (Sequential)',
    'single_agent_gpt': 'Single-LLM (GPT-5.2)',
    'single_agent_gemini': 'Single-LLM (Gemini)',
    'single_agent_claude': 'Single-LLM (Claude)',
    'no_agent': 'No-Agent Baseline'
}

# ============================================================================
# Helper Functions
# ============================================================================

def prepare_plot_data(df, top_systems):
    """Prepare aggregated data for plotting"""
    plot_data = []

    for system in top_systems:
        system_df = df[df['System'] == system]

        if len(system_df) == 0:
            continue

        accuracy_vals = system_df['Accuracy_%'].values
        accuracy_mean = np.mean(accuracy_vals)
        accuracy_sem = stats.sem(accuracy_vals) if len(accuracy_vals) > 1 else 0

        completion_vals = system_df['Completed'].values
        completion_rate = (np.sum(completion_vals) / len(completion_vals)) * 100

        intervention_vals = system_df['User_Interventions'].values
        intervention_rate = (np.sum(intervention_vals > 0) / len(system_df)) * 100

        cost_vals = system_df['Total_Cost_USD'].values
        cost_mean = np.mean(cost_vals) if len(cost_vals) > 0 else 0.0
        cost_sem = stats.sem(cost_vals) if len(cost_vals) > 1 else 0.0

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

def compute_significance(df, baseline_system, metric_column):
    """Compute p-values comparing each system to baseline"""
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
            significance[system] = 1.0

    return significance

def get_significance_marker(p_value):
    """Convert p-value to marker"""
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''

def assign_colors(systems):
    """Assign colors to systems"""
    color_map = {}
    multi_idx = 0
    single_idx = 0

    for system in systems:
        if 'no_agent' in system:
            color_map[system] = COLORS['baseline']
        elif 'single' in system:
            single_idx += 1
            color_map[system] = COLORS[f'single_llm_{single_idx}']
        else:
            multi_idx += 1
            color_map[system] = COLORS[f'multi_foundational_{multi_idx}']

    return color_map

# ============================================================================
# Plotting Functions
# ============================================================================

def plot_accuracy(plot_df, df_raw, baseline_system, output_prefix):
    """Figure 1: Decision Accuracy Comparison"""
    fig, ax = plt.subplots(figsize=(8, 5))

    significance = compute_significance(df_raw, baseline_system, 'Accuracy_%')
    colors = assign_colors(plot_df['system'].values)
    bar_colors = [colors[s] for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    ax.barh(y_pos, plot_df['accuracy_mean'], xerr=plot_df['accuracy_sem'],
            color=bar_colors, edgecolor='black', linewidth=1.5,
            error_kw={'linewidth': 2, 'ecolor': 'black'})

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

def plot_cost(plot_df, output_prefix):
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
# Main
# ============================================================================

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 bin/generate_plots_from_csv.py bulk_evaluation_v2_per_experiment.csv [output_prefix]")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else 'figure'

    print("\n" + "="*80)
    print("GENERATING PUBLICATION-QUALITY PLOTS (ICML 2026 Format)")
    print("="*80 + "\n")

    print("📖 Loading evaluation data...")
    df = pd.read_csv(csv_file)
    print(f"   ✓ Loaded {len(df)} experiments\n")

    # Select top 5 systems + baseline
    top_systems = [
        'role_perm_3',
        'default_sequential',
        'single_agent_gpt',
        'single_agent_gemini',
        'single_agent_claude'
    ]

    # Filter to only systems that exist in data
    available_systems = df['System'].unique()
    top_systems = [s for s in top_systems if s in available_systems]

    print("🎯 Selected Systems:")
    for system_name in top_systems:
        print(f"  • {SYSTEM_LABELS.get(system_name, system_name)}")
    print()

    # Prepare data
    plot_df = prepare_plot_data(df, top_systems)
    baseline_system = 'no_agent' if 'no_agent' in available_systems else top_systems[-1]

    # Generate plots
    print("📊 Generating Figure 1: Decision Accuracy...")
    plot_accuracy(plot_df, df, baseline_system, output_prefix)

    print("📊 Generating Figure 2: Completion Rate...")
    plot_completion(plot_df, output_prefix)

    print("📊 Generating Figure 3: Total Cost...")
    plot_cost(plot_df, output_prefix)

    print("📊 Generating Figure 4: User Intervention Rate...")
    plot_user_intervention(plot_df, output_prefix)

    print("\n" + "="*80)
    print("✨ All plots generated successfully!")
    print("="*80)
    print(f"\nOutput files: {output_prefix}_*.pdf / .png")
    print("Statistical significance: * p<0.05, ** p<0.01, *** p<0.001\n")

if __name__ == '__main__':
    main()
