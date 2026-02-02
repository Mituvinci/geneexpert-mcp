#!/usr/bin/env python3
"""
Generate Stage-Wise Accuracy Plot

Shows accuracy by pipeline stage to identify where poly-foundational
validation provides the greatest benefit (e.g., complex interpretation stages
like PCA/batch detection).

Usage:
    python3 bin/plot_stage_wise_accuracy.py bulk_evaluation_per_experiment.csv bulk
    python3 bin/plot_stage_wise_accuracy.py scrna_evaluation_per_experiment.csv scrna
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Publication style
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

def parse_system_type(system_name):
    """Categorize system"""
    if 'single_gpt' in system_name.lower() or system_name == 'single_gpt':
        return 'single_gpt'
    elif 'single_claude' in system_name.lower() or system_name == 'single_claude':
        return 'single_claude'
    elif 'single_gemini' in system_name.lower() or system_name == 'single_gemini':
        return 'single_gemini'
    elif 'no_agent' in system_name.lower() or system_name == 'no_agent':
        return 'no_agent'
    else:
        if 'parallel' in system_name.lower():
            return 'poly_parallel'
        elif 'sequential' in system_name.lower():
            return 'poly_sequential'
        else:
            return 'poly_other'

def calculate_stage_accuracy(df, pipeline_type):
    """Calculate accuracy per stage per system category"""

    # Define stages
    if pipeline_type == 'bulk':
        stages = [
            ('Stage1_Match', 'Stage 1\nFASTQ QC'),
            ('Stage2_Match', 'Stage 2\nAlignment QC'),
            ('Stage3_Match', 'Stage 3\nPCA + Batch\nDetection'),
            ('Stage4_Match', 'Stage 4\nDE Analysis')
        ]
    else:  # scrna
        stages = [
            ('Stage1_Match', 'Stage 1\nLoad Data'),
            ('Stage2_Match', 'Stage 2\nQC Filtering'),
            ('Stage3A_Match', 'Stage 3A\nCell Cycle\nDetection'),
            ('Stage3B_Match', 'Stage 3B\nNormalization'),
            ('Stage4_Match', 'Stage 4\nPCA'),
            ('Stage5_Match', 'Stage 5\nClustering')
        ]

    results = []

    for stage_col, stage_label in stages:
        # Aggregate by system category

        # Single-foundational (homogeneous)
        single_llm = df[df['System'].apply(parse_system_type).isin(['single_gpt', 'single_claude', 'single_gemini'])]
        if len(single_llm) > 0 and stage_col in single_llm.columns:
            valid = single_llm[stage_col].notna()
            if valid.sum() > 0:
                accuracy = (single_llm[valid][stage_col].sum() / valid.sum()) * 100
                results.append({
                    'Stage': stage_label,
                    'Category': 'Single-Foundational\n(Homogeneous)',
                    'Accuracy': accuracy,
                    'N': valid.sum()
                })

        # Poly-foundational
        poly = df[df['System'].apply(parse_system_type).isin(['poly_parallel', 'poly_sequential'])]
        if len(poly) > 0 and stage_col in poly.columns:
            valid = poly[stage_col].notna()
            if valid.sum() > 0:
                accuracy = (poly[valid][stage_col].sum() / valid.sum()) * 100
                results.append({
                    'Stage': stage_label,
                    'Category': 'Poly-Foundational\n(Diverse Models)',
                    'Accuracy': accuracy,
                    'N': valid.sum()
                })

    return pd.DataFrame(results)

def plot_stage_accuracy(stage_df, output_prefix, pipeline_type):
    """Generate grouped bar chart of stage-wise accuracy"""

    fig, ax = plt.subplots(figsize=(12, 6))

    stages = stage_df['Stage'].unique()
    categories = stage_df['Category'].unique()

    x = np.arange(len(stages))
    width = 0.35

    # Colors: Blue for poly-foundational, Red for single-foundational
    colors = {
        'Single-Foundational\n(Homogeneous)': '#D6604D',  # Red
        'Poly-Foundational\n(Diverse Models)': '#4393C3'  # Blue
    }

    for i, category in enumerate(categories):
        cat_data = stage_df[stage_df['Category'] == category]
        accuracies = [cat_data[cat_data['Stage'] == s]['Accuracy'].values[0]
                     if len(cat_data[cat_data['Stage'] == s]) > 0 else 0
                     for s in stages]

        offset = (i - len(categories)/2 + 0.5) * width
        bars = ax.bar(x + offset, accuracies, width, label=category,
                     color=colors.get(category, '#999999'),
                     edgecolor='black', linewidth=1)

        # Add value labels on bars
        for j, (bar, acc) in enumerate(zip(bars, accuracies)):
            if acc > 0:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                       f'{acc:.1f}%', ha='center', va='bottom',
                       fontsize=9, fontweight='bold')

    # Highlight complex stages
    if pipeline_type == 'bulk':
        # Stage 3 (PCA/Batch) is most complex
        ax.axvspan(1.5, 2.5, alpha=0.1, color='yellow', zorder=-1)
        ax.text(2, 105, 'Complex\nInterpretation', ha='center', va='bottom',
               fontsize=10, style='italic', bbox=dict(boxstyle='round',
               facecolor='yellow', alpha=0.3))
    else:  # scrna
        # Stage 3A (Cell Cycle) and Stage 5 (Clustering) are most complex
        ax.axvspan(1.5, 2.5, alpha=0.1, color='yellow', zorder=-1)
        ax.axvspan(4.5, 5.5, alpha=0.1, color='yellow', zorder=-1)
        ax.text(2, 105, 'Complex', ha='center', va='bottom',
               fontsize=9, style='italic', bbox=dict(boxstyle='round',
               facecolor='yellow', alpha=0.3))
        ax.text(5, 105, 'Complex', ha='center', va='bottom',
               fontsize=9, style='italic', bbox=dict(boxstyle='round',
               facecolor='yellow', alpha=0.3))

    ax.set_ylabel('Decision Accuracy (%)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Pipeline Stage', fontsize=13, fontweight='bold')
    ax.set_title(f'Stage-Wise Validation Accuracy ({pipeline_type.upper()}-seq)\n' +
                 'Poly-Foundational Advantage at Complex Interpretation Stages',
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(stages, fontsize=10)
    ax.set_ylim(0, 115)
    ax.legend(loc='upper right', fontsize=11, frameon=True, shadow=True)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_stage_accuracy.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_stage_accuracy.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_stage_accuracy.pdf / .png")

def plot_stage_improvement(stage_df, output_prefix, pipeline_type):
    """Plot improvement (poly - single) by stage"""

    stages = stage_df['Stage'].unique()

    improvements = []
    for stage in stages:
        stage_data = stage_df[stage_df['Stage'] == stage]

        poly_acc = stage_data[stage_data['Category'].str.contains('Poly')]['Accuracy'].values
        single_acc = stage_data[stage_data['Category'].str.contains('Single')]['Accuracy'].values

        if len(poly_acc) > 0 and len(single_acc) > 0:
            improvement = poly_acc[0] - single_acc[0]
            improvements.append({'Stage': stage, 'Improvement': improvement})

    imp_df = pd.DataFrame(improvements)

    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(imp_df))
    colors = ['#D73027' if imp < 0 else '#1A9850' for imp in imp_df['Improvement']]

    bars = ax.bar(x, imp_df['Improvement'], color=colors, edgecolor='black', linewidth=1.2)

    # Add value labels
    for i, (bar, imp) in enumerate(zip(bars, imp_df['Improvement'])):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.,
               height + (0.5 if height > 0 else -0.5),
               f'{imp:.1f}', ha='center', va='bottom' if height > 0 else 'top',
               fontsize=10, fontweight='bold')

    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.set_ylabel('Accuracy Improvement\n(Poly-Foundational - Single-Foundational, %)',
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('Pipeline Stage', fontsize=13, fontweight='bold')
    ax.set_title(f'Stage-Wise Improvement with Poly-Foundational System ({pipeline_type.upper()}-seq)',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(imp_df['Stage'], fontsize=10)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_stage_improvement.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_stage_improvement.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_stage_improvement.pdf / .png")

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 bin/plot_stage_wise_accuracy.py <csv_file> <bulk|scrna> [output_prefix]")
        sys.exit(1)

    csv_file = sys.argv[1]
    pipeline_type = sys.argv[2].lower()

    # Determine output directory based on pipeline type
    if pipeline_type == 'bulk':
        default_prefix = 'experiments/bulk_rna_csv_figures/bulk_stage_analysis'
    else:
        default_prefix = f'experiments/scrna_csv_figures/{pipeline_type}_stage_analysis'

    output_prefix = sys.argv[3] if len(sys.argv) > 3 else default_prefix

    # Ensure output directory exists
    output_dir = os.path.dirname(output_prefix) if '/' in output_prefix else '.'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"📁 Created output directory: {output_dir}\n")

    print("\n" + "="*80)
    print(f"STAGE-WISE ACCURACY ANALYSIS: {pipeline_type.upper()}-seq")
    print("="*80 + "\n")

    # Load data
    print(f"📖 Loading: {csv_file}")
    df = pd.read_csv(csv_file)
    print(f"   ✓ Loaded {len(df)} experiments\n")

    # Calculate stage accuracy
    print("🔍 Calculating stage-wise accuracy...")
    stage_df = calculate_stage_accuracy(df, pipeline_type)
    print(f"   ✓ Analyzed {len(stage_df['Stage'].unique())} stages\n")

    print(stage_df.to_string(index=False))
    print()

    # Save results
    stage_csv = f'{output_prefix}_stage_accuracy.csv'
    stage_df.to_csv(stage_csv, index=False)
    print(f"✓ Saved results: {stage_csv}\n")

    # Generate plots
    print("📊 Generating Figure 1: Stage-Wise Accuracy Comparison...")
    plot_stage_accuracy(stage_df, output_prefix, pipeline_type)

    print("📊 Generating Figure 2: Stage-Wise Improvement...")
    plot_stage_improvement(stage_df, output_prefix, pipeline_type)

    print("\n" + "="*80)
    print("✨ Stage-wise analysis complete!")
    print("="*80)
    print(f"\nOutput files:")
    print(f"  • {output_prefix}_stage_accuracy.pdf / .png")
    print(f"  • {output_prefix}_stage_improvement.pdf / .png")
    print(f"  • {output_prefix}_stage_accuracy.csv")
    print()

if __name__ == '__main__':
    main()
