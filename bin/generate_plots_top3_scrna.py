#!/usr/bin/env python3
"""
Generate Publication-Quality Plots for scRNA-seq with TOP 3 PolyLLM + Baseline + 3 Single-LLM

Plots 7 systems total:
- Top 3 best-performing PolyLLM multi-agent systems (by accuracy)
- No-agent baseline
- Single-LLM (Claude, GPT, Gemini)

Usage:
    python3 bin/generate_plots_top3_scrna.py scrna_evaluation.csv [output_prefix]
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
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13
plt.rcParams['legend.fontsize'] = 13

# Colorblind-friendly palette
COLORS = {
    'parallel': '#2166AC',
    'sequential': '#4393C3',
    'single_claude': '#D6604D',
    'single_gpt': '#F4A582',
    'single_gemini': '#FDB863',
    'baseline': '#999999'
}

def parse_system_config(system_code):
    """Parse full system configuration name to extract architecture and roles"""
    if 'no_agent' in system_code or system_code == 'no_agent':
        return {'type': 'no_agent', 'arch': None, 'roles': None}
    elif 'single_claude' in system_code or system_code == 'single_claude':
        return {'type': 'single_claude', 'arch': None, 'roles': None}
    elif 'single_gpt' in system_code or system_code == 'single_gpt':
        return {'type': 'single_gpt', 'arch': None, 'roles': None}
    elif 'single_gemini' in system_code or system_code == 'single_gemini':
        return {'type': 'single_gemini', 'arch': None, 'roles': None}

    # Parse parallel/sequential multi-agent systems
    arch = 'Parallel' if 'parallel' in system_code.lower() else 'Sequential'

    # Extract role configuration from system name
    if 'gp_st_cl_bl_gm_pl' in system_code:
        roles = 'GPT=stats, Cl=bio, Gem=pipe'
    elif 'gp_st_cl_pl_gm_bl' in system_code:
        roles = 'GPT=stats, Cl=pipe, Gem=bio'
    elif 'gp_pl_cl_st_gm_bl' in system_code:
        roles = 'GPT=pipe, Cl=stats, Gem=bio'
    elif 'gp_pl_cl_bl_gm_st' in system_code:
        roles = 'GPT=pipe, Cl=bio, Gem=stats'
    elif 'gp_bl_cl_st_gm_pl' in system_code:
        roles = 'GPT=bio, Cl=stats, Gem=pipe'
    elif 'gp_bl_cl_pl_gm_st' in system_code:
        roles = 'GPT=bio, Cl=pipe, Gem=stats'
    elif '_default' in system_code or system_code in ['parallel_default', 'sequential_default']:
        roles = 'GPT=stats, Cl=pipe, Gem=bio'
    else:
        roles = 'default'

    return {'type': 'multi_agent', 'arch': arch, 'roles': roles}

def get_system_label(system_code):
    """Convert system code to descriptive label with role info"""
    parsed = parse_system_config(system_code)

    if parsed['type'] == 'no_agent':
        return 'No-Agent Baseline'
    elif parsed['type'] == 'single_claude':
        return 'Single-LLM (Claude)'
    elif parsed['type'] == 'single_gpt':
        return 'Single-LLM (GPT-5.2)'
    elif parsed['type'] == 'single_gemini':
        return 'Single-LLM (Gemini)'
    elif parsed['type'] == 'multi_agent':
        return f"PolyLLM ({parsed['arch']}): {parsed['roles']}"
    else:
        return system_code

def assign_color(system_code):
    """Assign color based on system type"""
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
        if parsed['arch'] == 'Sequential':
            return COLORS['sequential']
        else:
            return COLORS['parallel']
    else:
        return COLORS['baseline']

# ============================================================================
# Helper Functions
# ============================================================================

def select_top_systems(df):
    """Select top 3 PolyLLM systems + baseline + 3 single-LLM"""

    # Separate PolyLLM multi-agent systems
    polyllm_systems = []
    for system in df['System'].unique():
        parsed = parse_system_config(system)
        if parsed['type'] == 'multi_agent':
            polyllm_systems.append(system)

    # Get accuracy from aggregated data (Decision_Accuracy_% column)
    polyllm_accuracy = []
    for system in polyllm_systems:
        system_row = df[df['System'] == system]
        if len(system_row) > 0:
            # Use Decision_Accuracy_% from aggregated CSV (overall accuracy)
            acc = system_row['Decision_Accuracy_%'].values[0]
            polyllm_accuracy.append({'System': system, 'Mean_Accuracy': acc})

    polyllm_df = pd.DataFrame(polyllm_accuracy).sort_values('Mean_Accuracy', ascending=False)

    # Get top 3 PolyLLM systems
    top3_polyllm = polyllm_df.head(3)['System'].tolist()

    print("🏆 Top 3 PolyLLM Multi-Agent Systems (by accuracy):")
    for i, system in enumerate(top3_polyllm, 1):
        acc = polyllm_df[polyllm_df['System'] == system]['Mean_Accuracy'].values[0]
        label = get_system_label(system)
        print(f"   {i}. {label}: {acc:.1f}%")
    print()

    # Fixed systems: baseline + 3 single-LLM (find actual names in CSV)
    fixed_systems = []
    for system in df['System'].unique():
        parsed = parse_system_config(system)
        if parsed['type'] in ['single_gpt', 'single_claude', 'single_gemini', 'no_agent']:
            fixed_systems.append(system)

    # Combine
    selected_systems = top3_polyllm + fixed_systems

    return selected_systems

def prepare_plot_data(df_agg, df_per_exp, systems_to_plot):
    """Prepare aggregated data for plotting"""
    plot_data = []

    for system in systems_to_plot:
        system_row = df_agg[df_agg['System'] == system]

        if len(system_row) == 0:
            continue

        # Use aggregated metrics (matches evaluation summary)
        accuracy_mean = system_row['Decision_Accuracy_%'].values[0]
        completion_rate = system_row['Completion_Rate_%'].values[0]
        n_experiments = system_row['Total_Experiments'].values[0]

        # User intervention rate
        user_inputs = system_row['User_Interventions'].values[0]
        intervention_rate = (user_inputs / system_row['Total_Decisions'].values[0] * 100) if system_row['Total_Decisions'].values[0] > 0 else 0

        # Calculate cost from per-experiment data
        per_exp_system = df_per_exp[df_per_exp['System'] == system]
        cost_vals = per_exp_system['Total_Cost_USD'].values
        parsed = parse_system_config(system)
        if parsed['type'] in ['single_claude', 'single_gpt', 'single_gemini']:
            cost_vals = cost_vals * 3
        cost_total = np.sum(cost_vals) if len(cost_vals) > 0 else 0.0

        plot_data.append({
            'system': system,
            'label': get_system_label(system),
            'accuracy_mean': accuracy_mean,
            'accuracy_sem': 0,  # No SEM for aggregated data
            'completion_rate': completion_rate,
            'intervention_rate': intervention_rate,
            'cost_total': cost_total,
            'n_experiments': n_experiments
        })

    return pd.DataFrame(plot_data)

def compute_all_ttests(df_agg, df_per_exp, systems_to_plot, baseline_system='single_agent_gpt'):
    """Compute all pairwise t-tests"""
    results = []

    baseline_acc = df_per_exp[df_per_exp['System'] == baseline_system]['Accuracy_%'].values

    print("\n" + "="*80)
    print("T-TEST RESULTS: Decision Accuracy vs Single-LLM (GPT-5.2)")
    print("="*80)
    print(f"{'System':<45} {'Mean Acc':<12} {'n':<6} {'p-value':<12} {'Sig':<6}")
    print("-"*80)

    # Add baseline system first
    # Use aggregated accuracy for display
    baseline_row = df_agg[df_agg['System'] == baseline_system]
    baseline_mean_agg = baseline_row['Decision_Accuracy_%'].values[0] if len(baseline_row) > 0 else np.mean(baseline_acc)
    baseline_label = get_system_label(baseline_system)
    print(f"{baseline_label:<45} {baseline_mean_agg:>10.1f}% {len(baseline_acc):>5} {'---':>11} {'baseline':>6}")
    results.append({
        'System': baseline_system,
        'Label': baseline_label,
        'Mean_Accuracy': baseline_mean_agg,
        'n_experiments': len(baseline_acc),
        't_statistic': 0,
        'p_value': 1.0,
        'significance': 'baseline'
    })

    # Compare other systems to baseline
    for system in systems_to_plot:
        if system == baseline_system:
            continue

        system_acc = df_per_exp[df_per_exp['System'] == system]['Accuracy_%'].values

        if len(system_acc) > 1 and len(baseline_acc) > 1:
            t_stat, p_value = stats.ttest_ind(system_acc, baseline_acc)
            sig = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
        else:
            t_stat, p_value, sig = 0, 1.0, 'N/A'

        # Use aggregated accuracy for display
        system_row = df_agg[df_agg['System'] == system]
        mean_acc_agg = system_row['Decision_Accuracy_%'].values[0] if len(system_row) > 0 else np.mean(system_acc)
        label = get_system_label(system)

        print(f"{label:<45} {mean_acc_agg:>10.1f}% {len(system_acc):>5} {p_value:>11.4f} {sig:>6}")

        results.append({
            'System': system,
            'Label': label,
            'Mean_Accuracy': mean_acc_agg,
            'n_experiments': len(system_acc),
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
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = [assign_color(s) for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    # NO error bars on accuracy plot
    bars = ax.barh(y_pos, plot_df['accuracy_mean'],
                   color=colors, edgecolor='black', linewidth=1.5)

    # Add percentage labels AND significance markers
    for i, system in enumerate(plot_df['system']):
        mean_val = plot_df.iloc[i]['accuracy_mean']

        # Add percentage value
        ax.text(mean_val + 1, i, f'{mean_val:.1f}%', va='center', fontsize=12, fontweight='bold')

        # Add significance marker
        ttest_row = ttest_results[ttest_results['System'] == system]
        if len(ttest_row) > 0:
            marker = get_significance_marker(ttest_row.iloc[0]['p_value'])
            if marker:
                ax.text(mean_val + 8, i, marker, va='center', fontsize=16, fontweight='bold')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'])
    ax.set_xlabel('Decision Accuracy (%)', fontsize=16)
    ax.set_xlim(0, 105)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_accuracy.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_accuracy.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_accuracy.pdf / .png")

def plot_completion(plot_df, output_prefix):
    """Figure 2: Pipeline Completion Rate"""
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = [assign_color(s) for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    bars = ax.barh(y_pos, plot_df['completion_rate'],
                   color=colors, edgecolor='black', linewidth=1.5)

    # Add value labels on bars
    for i, rate in enumerate(plot_df['completion_rate']):
        x_pos = max(rate + 2, 5)  # Position text at least 5% from left edge
        ax.text(x_pos, i, f'{rate:.1f}%', va='center', fontsize=12, fontweight='bold')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'])
    ax.set_xlabel('Completion Rate (%)', fontsize=16)
    ax.set_xlim(0, 105)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_completion.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_completion.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_completion.pdf / .png")

def plot_cost(plot_df, output_prefix):
    """Figure 3: Total Cost Across All Datasets"""
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = [assign_color(s) for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    # No error bars, no text labels (cleaner visualization)
    bars = ax.barh(y_pos, plot_df['cost_total'],
                   color=colors, edgecolor='black', linewidth=1.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'])
    ax.set_xlabel('Total Cost Across All Datasets (USD)', fontsize=16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_cost.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_cost.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_cost.pdf / .png")

def plot_user_intervention(plot_df, output_prefix):
    """Figure 4: User Intervention Rate"""
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = [assign_color(s) for s in plot_df['system'].values]

    y_pos = np.arange(len(plot_df))
    bars = ax.barh(y_pos, plot_df['intervention_rate'],
                   color=colors, edgecolor='black', linewidth=1.5)

    # Add percentage labels on bars
    for i, rate in enumerate(plot_df['intervention_rate']):
        x_pos = max(rate + 1, 3)  # Position text at least 3% from left edge
        ax.text(x_pos, i, f'{rate:.1f}%', va='center', fontsize=12, fontweight='bold')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['label'])
    ax.set_xlabel('User Intervention Rate (%)', fontsize=16)
    ax.set_xlim(0, 105)
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
        print("Usage: python3 bin/generate_plots_top3_scrna.py scrna_evaluation.csv [output_prefix]")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else 'scrna_figure'

    print("\n" + "="*80)
    print("GENERATING scRNA-seq PLOTS: TOP 3 PolyLLM + Baseline + 3 Single-LLM")
    print("="*80 + "\n")

    print("📖 Loading evaluation data...")

    # Determine if this is aggregated or per-experiment CSV
    df = pd.read_csv(csv_file)

    if 'Decision_Accuracy_%' in df.columns:
        # Aggregated CSV - need to also load per-experiment for cost
        print("   Using aggregated CSV (matches evaluation summary)")
        per_exp_file = csv_file.replace('.csv', '_per_experiment.csv')
        df_per_exp = pd.read_csv(per_exp_file)
        print(f"   ✓ Loaded {len(df)} systems from aggregated CSV")
        print(f"   ✓ Loaded {len(df_per_exp)} experiments from per-experiment CSV\n")
    else:
        # Per-experiment CSV - need to aggregate
        print("   Per-experiment CSV detected - using aggregated CSV instead")
        agg_file = csv_file.replace('_per_experiment.csv', '.csv')
        df = pd.read_csv(agg_file)
        df_per_exp = pd.read_csv(csv_file)
        print(f"   ✓ Loaded {len(df)} systems from aggregated CSV")
        print(f"   ✓ Loaded {len(df_per_exp)} experiments from per-experiment CSV\n")

    # Select top systems
    selected_systems = select_top_systems(df)

    print("📊 Selected 7 systems for plotting:")
    for sys_code in selected_systems:
        system_row = df[df['System'] == sys_code]
        n_exp = system_row['Total_Experiments'].values[0] if len(system_row) > 0 else 0
        print(f"   • {get_system_label(sys_code)}: {int(n_exp)} experiments")
    print()

    # Find actual single_gpt system name for t-test baseline
    baseline_system = None
    for system in selected_systems:
        if parse_system_config(system)['type'] == 'single_gpt':
            baseline_system = system
            break

    if not baseline_system:
        baseline_system = selected_systems[-1]  # Fallback to last system

    # Compute t-tests
    ttest_df = compute_all_ttests(df, df_per_exp, selected_systems, baseline_system=baseline_system)
    ttest_file = output_prefix + '_ttests.csv'
    ttest_df.to_csv(ttest_file, index=False)
    print(f"✓ Saved t-test results to: {ttest_file}\n")

    # Prepare data
    plot_df = prepare_plot_data(df, df_per_exp, selected_systems)

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
    print("✨ All scRNA-seq plots generated successfully!")
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
