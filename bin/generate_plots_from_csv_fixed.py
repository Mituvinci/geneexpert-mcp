#!/usr/bin/env python3
"""
Generate Publication-Quality Plots from CSV Files (ICML 2026) - FIXED VERSION

Includes:
- All parallel and sequential systems (not just top 5)
- Descriptive labels showing model roles
- T-test output printed to console and file

Usage:
    python3 bin/generate_plots_from_csv_fixed.py bulk_evaluation_final_per_experiment.csv [output_prefix]
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
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['legend.fontsize'] = 11

# Colorblind-friendly palette
COLORS = {
    'parallel': '#2166AC',
    'sequential': '#4393C3',
    'single_claude': '#D6604D',
    'single_gpt': '#F4A582',
    'single_gemini': '#FDB863',
    'baseline': '#999999'
}

# Role permutation mappings
ROLE_DESCRIPTIONS = {
    'perm_1': 'GPT=stats, Claude=biology, Gemini=pipeline',
    'perm_2': 'GPT=stats, Claude=pipeline, Gemini=biology',
    'perm_3': 'GPT=pipeline, Claude=stats, Gemini=biology',
    'perm_4': 'GPT=pipeline, Claude=biology, Gemini=stats',
    'perm_5': 'GPT=biology, Claude=stats, Gemini=pipeline',
    'perm_6': 'GPT=biology, Claude=pipeline, Gemini=stats'
}

def get_system_label(system_code):
    """Convert system code to descriptive label"""
    if system_code == 'no_agent':
        return 'No-Agent Baseline'
    elif system_code == 'single_agent_claude':
        return 'Single-LLM (Claude)'
    elif system_code == 'single_agent_gpt':
        return 'Single-LLM (GPT-5.2)'
    elif system_code == 'single_agent_gemini':
        return 'Single-LLM (Gemini)'
    elif system_code == 'default_parallel':
        return 'PolyLLM Multi-Agent (Parallel)'
    elif system_code == 'default_sequential':
        return 'PolyLLM Multi-Agent (Sequential)'
    elif system_code.startswith('parallel_perm_'):
        return 'PolyLLM Multi-Agent (Parallel)'
    elif system_code.startswith('sequential_perm_'):
        return 'PolyLLM Multi-Agent (Sequential)'
    else:
        return system_code

def assign_color(system_code):
    """Assign color based on system type"""
    if 'no_agent' in system_code:
        return COLORS['baseline']
    elif 'single_claude' in system_code:
        return COLORS['single_claude']
    elif 'single_gpt' in system_code:
        return COLORS['single_gpt']
    elif 'single_gemini' in system_code:
        return COLORS['single_gemini']
    elif 'sequential' in system_code:
        return COLORS['sequential']
    elif 'parallel' in system_code:
        return COLORS['parallel']
    else:
        return COLORS['baseline']

# ============================================================================
# Helper Functions
# ============================================================================

def prepare_plot_data(df, systems_to_plot):
    """Prepare aggregated data for plotting"""
    plot_data = []

    for system in systems_to_plot:
        system_df = df[df['System'] == system]

        if len(system_df) == 0:
            continue

        accuracy_vals = system_df['Accuracy_%'].values
        accuracy_mean = np.mean(accuracy_vals)
        accuracy_sem = stats.sem(accuracy_vals) if len(accuracy_vals) > 1 else 0

        completion_vals = system_df['Completed'].values
        completion_rate = (np.sum(completion_vals) / len(completion_vals)) * 100 if len(completion_vals) > 0 else 0

        intervention_vals = system_df['User_Interventions'].values
        intervention_rate = (np.sum(intervention_vals > 0) / len(system_df)) * 100 if len(system_df) > 0 else 0

        cost_vals = system_df['Total_Cost_USD'].values
        cost_mean = np.mean(cost_vals) if len(cost_vals) > 0 else 0.0
        cost_sem = stats.sem(cost_vals) if len(cost_vals) > 1 else 0.0

        plot_data.append({
            'system': system,
            'label': get_system_label(system),
            'accuracy_mean': accuracy_mean,
            'accuracy_sem': accuracy_sem,
            'completion_rate': completion_rate,
            'intervention_rate': intervention_rate,
            'cost_mean': cost_mean,
            'cost_sem': cost_sem,
            'n_experiments': len(system_df)
        })

    return pd.DataFrame(plot_data)

def compute_all_ttests(df, baseline_system='single_agent_gpt'):
    """Compute all pairwise t-tests"""
    results = []

    systems = df['System'].unique()
    baseline_acc = df[df['System'] == baseline_system]['Accuracy_%'].values

    print("\n" + "="*80)
    print("T-TEST RESULTS: Decision Accuracy vs Single-LLM (GPT-5.2)")
    print("="*80)
    print(f"{'System':<40} {'Mean Acc':<12} {'p-value':<12} {'Significance':<12}")
    print("-"*80)

    for system in systems:
        if system == baseline_system:
            continue

        system_acc = df[df['System'] == system]['Accuracy_%'].values

        if len(system_acc) > 1 and len(baseline_acc) > 1:
            t_stat, p_value = stats.ttest_ind(system_acc, baseline_acc)
            sig = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
        else:
            t_stat, p_value, sig = 0, 1.0, 'N/A'

        mean_acc = np.mean(system_acc)
        label = get_system_label(system)

        print(f"{label:<40} {mean_acc:>10.1f}% {p_value:>11.4f} {sig:>12}")

        results.append({
            'System': system,
            'Label': label,
            'Mean_Accuracy': mean_acc,
            't_statistic': t_stat,
            'p_value': p_value,
            'significance': sig
        })

    print("="*80)
    print("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant")
    print("="*80 + "\n")

    return pd.DataFrame(results)

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

# ============================================================================
# Plotting Functions
# ============================================================================

def plot_accuracy(plot_df, ttest_results, output_prefix):
    """Figure 1: Decision Accuracy Comparison"""
    fig, ax = plt.subplots(figsize=(10, 8))

    colors = [assign_color(s) for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    ax.barh(y_pos, plot_df['accuracy_mean'], xerr=plot_df['accuracy_sem'],
            color=colors, edgecolor='black', linewidth=1.5,
            error_kw={'linewidth': 2, 'ecolor': 'black'})

    # Add significance markers
    for i, system in enumerate(plot_df['system']):
        mean_val = plot_df.iloc[i]['accuracy_mean']
        ttest_row = ttest_results[ttest_results['System'] == system]
        if len(ttest_row) > 0:
            marker = get_significance_marker(ttest_row.iloc[0]['p_value'])
            if marker:
                ax.text(mean_val + 2, i, marker, va='center', fontsize=14, fontweight='bold')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'], fontsize=10)
    ax.set_xlabel('Decision Accuracy (%)', fontsize=14)
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
    fig, ax = plt.subplots(figsize=(10, 8))

    colors = [assign_color(s) for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    ax.barh(y_pos, plot_df['completion_rate'],
            color=colors, edgecolor='black', linewidth=1.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'], fontsize=10)
    ax.set_xlabel('Completion Rate (%)', fontsize=14)
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
    fig, ax = plt.subplots(figsize=(10, 8))

    colors = [assign_color(s) for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    ax.barh(y_pos, plot_df['cost_mean'], xerr=plot_df['cost_sem'],
            color=colors, edgecolor='black', linewidth=1.5,
            error_kw={'linewidth': 2, 'ecolor': 'black'})

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'], fontsize=10)
    ax.set_xlabel('Total Cost per Experiment (USD)', fontsize=14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_cost.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_cost.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_cost.pdf / .png")

def plot_user_intervention(plot_df, output_prefix):
    """Figure 4: User Intervention Rate"""
    fig, ax = plt.subplots(figsize=(10, 8))

    colors = [assign_color(s) for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    ax.barh(y_pos, plot_df['intervention_rate'],
            color=colors, edgecolor='black', linewidth=1.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'], fontsize=10)
    ax.set_xlabel('User Intervention Rate (%)', fontsize=14)
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
        print("Usage: python3 bin/generate_plots_from_csv_fixed.py bulk_evaluation_final_per_experiment.csv [output_prefix]")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else 'figure'

    print("\n" + "="*80)
    print("GENERATING PUBLICATION-QUALITY PLOTS (ICML 2026 Format) - FIXED")
    print("="*80 + "\n")

    print("📖 Loading evaluation data...")
    df = pd.read_csv(csv_file)
    print(f"   ✓ Loaded {len(df)} experiments\n")

    # Use ALL systems (not just top 5)
    all_systems = sorted(df['System'].unique())

    print(f"🎯 Plotting {len(all_systems)} systems:")
    for system_name in all_systems:
        print(f"  • {get_system_label(system_name)}")
    print()

    # Compute t-tests
    ttest_df = compute_all_ttests(df, baseline_system='single_agent_gpt')
    ttest_file = output_prefix + '_ttests.csv'
    ttest_df.to_csv(ttest_file, index=False)
    print(f"✓ Saved t-test results to: {ttest_file}\n")

    # Prepare data
    plot_df = prepare_plot_data(df, all_systems)

    # Generate plots
    print("📊 Generating Figure 1: Decision Accuracy...")
    plot_accuracy(plot_df, ttest_df, output_prefix)

    print("📊 Generating Figure 2: Completion Rate...")
    plot_completion(plot_df, output_prefix)

    print("📊 Generating Figure 3: Total Cost...")
    plot_cost(plot_df, output_prefix)

    print("📊 Generating Figure 4: User Intervention Rate...")
    plot_user_intervention(plot_df, output_prefix)

    print("\n" + "="*80)
    print("✨ All plots generated successfully!")
    print("="*80)
    print(f"\nOutput files:")
    print(f"  • {output_prefix}_accuracy.pdf / .png")
    print(f"  • {output_prefix}_completion.pdf / .png")
    print(f"  • {output_prefix}_cost.pdf / .png")
    print(f"  • {output_prefix}_user_intervention.pdf / .png")
    print(f"  • {output_prefix}_ttests.csv")
    print("\nStatistical significance: * p<0.05, ** p<0.01, *** p<0.001\n")

if __name__ == '__main__':
    main()
