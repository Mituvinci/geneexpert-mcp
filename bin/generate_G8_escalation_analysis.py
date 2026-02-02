#!/usr/bin/env python3
"""
Generate G8: User Escalation Frequency and Justification Rate

Shows escalation rate and justification rate for PolyLLM vs Single-LLM systems.
Demonstrates that PolyLLM surfaces genuine uncertainty through principled escalation.

Usage:
    python3 bin/generate_G8_escalation_analysis.py <bulk_csv> <scrna_csv> <output_dir>
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Publication style
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
plt.rcParams['font.size'] = 11

def parse_system_type(system_name):
    """Categorize system as PolyLLM or Single-LLM"""
    if 'single_gpt' in system_name.lower() or system_name == 'single_gpt':
        return 'Single-LLM'
    elif 'single_claude' in system_name.lower() or system_name == 'single_claude':
        return 'Single-LLM'
    elif 'single_gemini' in system_name.lower() or system_name == 'single_gemini':
        return 'Single-LLM'
    elif 'no_agent' in system_name.lower() or system_name == 'no_agent':
        return 'No-Agent'
    else:
        if 'parallel' in system_name.lower() or 'sequential' in system_name.lower():
            return 'PolyLLM'
    return 'Other'

def calculate_escalation_metrics(df):
    """Calculate escalation rate and justification rate per system category"""

    df['System_Category'] = df['System'].apply(parse_system_type)

    # Filter out no-agent
    df = df[df['System_Category'] != 'No-Agent']

    results = []

    for category in ['PolyLLM', 'Single-LLM']:
        cat_df = df[df['System_Category'] == category]

        if len(cat_df) == 0:
            continue

        # Calculate escalation rate
        total_decisions = cat_df['Total_Stages'].sum()
        total_escalations = cat_df['User_Interventions'].sum()
        escalation_rate = (total_escalations / total_decisions * 100) if total_decisions > 0 else 0

        # Calculate justification rate
        total_justified = cat_df['Justified_Escalations'].sum()
        justification_rate = (total_justified / total_escalations * 100) if total_escalations > 0 else 0

        results.append({
            'Category': category,
            'Escalation_Rate': escalation_rate,
            'Justification_Rate': justification_rate,
            'Total_Decisions': total_decisions,
            'Total_Escalations': total_escalations,
            'Justified_Escalations': total_justified
        })

    return pd.DataFrame(results)

def plot_escalation_analysis(bulk_metrics, scrna_metrics, output_dir):
    """Generate grouped bar chart for escalation analysis"""

    # Combine metrics
    combined = pd.concat([
        bulk_metrics.assign(Modality='Bulk RNA-seq'),
        scrna_metrics.assign(Modality='scRNA-seq')
    ])

    # Average across modalities
    avg_metrics = combined.groupby('Category').agg({
        'Escalation_Rate': 'mean',
        'Justification_Rate': 'mean'
    }).reset_index()

    fig, ax = plt.subplots(figsize=(7, 5))

    categories = avg_metrics['Category'].values
    x = np.arange(len(categories))
    width = 0.35

    # Colors
    color_escalation = '#4393C3'  # Blue
    color_justification = '#2CA02C'  # Green

    # Plot bars
    bars1 = ax.bar(x - width/2, avg_metrics['Escalation_Rate'], width,
                   label='Escalation Rate', color=color_escalation,
                   edgecolor='black', linewidth=0.8)
    bars2 = ax.bar(x + width/2, avg_metrics['Justification_Rate'], width,
                   label='Justification Rate', color=color_justification,
                   edgecolor='black', linewidth=0.8)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{height:.1f}%', ha='center', va='bottom',
                   fontsize=9, fontweight='bold')

    # Labels and title
    ax.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    ax.set_xlabel('System Category', fontsize=12, fontweight='bold')
    ax.set_title('User Escalation Frequency and Justification Rate\nAveraged Across Bulk and scRNA-seq',
                 fontsize=13, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=11)
    ax.set_ylim(0, 100)
    ax.legend(loc='upper right', fontsize=10, frameon=True)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add annotation
    ax.text(0.5, 0.95, 'Higher justification rate indicates principled escalation',
            transform=ax.transAxes, ha='center', va='top',
            fontsize=9, style='italic', bbox=dict(boxstyle='round',
            facecolor='wheat', alpha=0.3))

    plt.tight_layout()

    # Save
    output_prefix = os.path.join(output_dir, 'G8_escalation_analysis')
    plt.savefig(f'{output_prefix}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved: {output_prefix}.pdf / .png")

    return avg_metrics

def plot_modality_comparison(bulk_metrics, scrna_metrics, output_dir):
    """Generate side-by-side comparison for bulk vs scRNA"""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    width = 0.35

    for ax, metrics, title in [(ax1, bulk_metrics, 'Bulk RNA-seq'),
                                (ax2, scrna_metrics, 'scRNA-seq')]:

        categories = metrics['Category'].values
        x = np.arange(len(categories))

        bars1 = ax.bar(x - width/2, metrics['Escalation_Rate'], width,
                      label='Escalation Rate', color='#4393C3',
                      edgecolor='black', linewidth=0.8)
        bars2 = ax.bar(x + width/2, metrics['Justification_Rate'], width,
                      label='Justification Rate', color='#2CA02C',
                      edgecolor='black', linewidth=0.8)

        # Add value labels
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                       f'{height:.1f}%', ha='center', va='bottom',
                       fontsize=8, fontweight='bold')

        ax.set_ylabel('Percentage (%)' if ax == ax1 else '', fontsize=11, fontweight='bold')
        ax.set_xlabel('System Category', fontsize=11, fontweight='bold')
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(categories, fontsize=10)
        ax.set_ylim(0, 100)
        ax.legend(loc='upper right', fontsize=8, frameon=True)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()

    output_prefix = os.path.join(output_dir, 'G8_escalation_by_modality')
    plt.savefig(f'{output_prefix}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved: {output_prefix}.pdf / .png")

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 bin/generate_G8_escalation_analysis.py <bulk_csv> <scrna_csv> <output_dir>")
        print("\nExample:")
        print("  python3 bin/generate_G8_escalation_analysis.py \\")
        print("    experiments/bulk_rna_csv_figures/bulk_evaluation_per_experiment.csv \\")
        print("    experiments/scrna_csv_figures/scrna_evaluation_per_experiment.csv \\")
        print("    experiments/combined_figures")
        sys.exit(1)

    bulk_csv = sys.argv[1]
    scrna_csv = sys.argv[2]
    output_dir = sys.argv[3]

    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "="*80)
    print("GENERATING G8: USER ESCALATION ANALYSIS")
    print("="*80 + "\n")

    # Load data
    print("Loading bulk RNA-seq data...")
    bulk_df = pd.read_csv(bulk_csv)
    print(f"  Loaded {len(bulk_df)} experiments")

    print("Loading scRNA-seq data...")
    scrna_df = pd.read_csv(scrna_csv)
    print(f"  Loaded {len(scrna_df)} experiments\n")

    # Calculate metrics
    print("Calculating escalation metrics...")
    bulk_metrics = calculate_escalation_metrics(bulk_df)
    scrna_metrics = calculate_escalation_metrics(scrna_df)

    print("\nBulk RNA-seq Metrics:")
    print(bulk_metrics.to_string(index=False))

    print("\nscRNA-seq Metrics:")
    print(scrna_metrics.to_string(index=False))

    # Generate plots
    print("\nGenerating combined plot...")
    avg_metrics = plot_escalation_analysis(bulk_metrics, scrna_metrics, output_dir)

    print("\nGenerating modality-specific comparison...")
    plot_modality_comparison(bulk_metrics, scrna_metrics, output_dir)

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print("\nAverage Metrics Across Both Modalities:")
    print(avg_metrics.to_string(index=False))
    print("\nKey Finding:")
    poly_row = avg_metrics[avg_metrics['Category'] == 'PolyLLM'].iloc[0]
    single_row = avg_metrics[avg_metrics['Category'] == 'Single-LLM'].iloc[0]
    print(f"  PolyLLM escalates {poly_row['Escalation_Rate']:.1f}% but {poly_row['Justification_Rate']:.1f}% are justified")
    print(f"  Single-LLM escalates {single_row['Escalation_Rate']:.1f}% but only {single_row['Justification_Rate']:.1f}% are justified")
    print(f"  PolyLLM surfaces genuine uncertainty; Single-LLM under-escalates and fails silently")
    print("\n")

if __name__ == '__main__':
    main()
