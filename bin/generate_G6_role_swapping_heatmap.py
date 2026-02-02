#!/usr/bin/env python3
"""
Generate G6: Role-Swapping Performance Matrix
Shows overall accuracy (%) for all 6 role permutations tested
Outputs: 3×3 heatmap (GPT role × Claude role) with Gemini role implicit
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
import os

def parse_role_configuration(system):
    """Extract role assignments from system name"""
    if 'gp_st_cl_pl_gm_bl' in system:
        return {'gpt': 'Statistics', 'claude': 'Pipeline', 'gemini': 'Biology'}
    elif 'gp_st_cl_bl_gm_pl' in system:
        return {'gpt': 'Statistics', 'claude': 'Biology', 'gemini': 'Pipeline'}
    elif 'gp_pl_cl_st_gm_bl' in system:
        return {'gpt': 'Pipeline', 'claude': 'Statistics', 'gemini': 'Biology'}
    elif 'gp_pl_cl_bl_gm_st' in system:
        return {'gpt': 'Pipeline', 'claude': 'Biology', 'gemini': 'Statistics'}
    elif 'gp_bl_cl_st_gm_pl' in system:
        return {'gpt': 'Biology', 'claude': 'Statistics', 'gemini': 'Pipeline'}
    elif 'gp_bl_cl_pl_gm_st' in system:
        return {'gpt': 'Biology', 'claude': 'Pipeline', 'gemini': 'Statistics'}
    else:
        return None

def generate_role_swapping_heatmap(csv_file, output_dir, modality):
    """Generate 3×3 role-swapping performance matrix"""

    print(f"Reading {csv_file}...")
    df = pd.read_csv(csv_file)

    # Filter to only parallel role-swapping configurations (6 permutations)
    role_systems = [
        'parallel_gp_st_cl_pl_gm_bl',
        'parallel_gp_st_cl_bl_gm_pl',
        'parallel_gp_pl_cl_st_gm_bl',
        'parallel_gp_pl_cl_bl_gm_st',
        'parallel_gp_bl_cl_st_gm_pl',
        'parallel_gp_bl_cl_pl_gm_st'
    ]

    df_roles = df[df['System'].isin(role_systems)]

    if df_roles.empty:
        print(f"Error: No role-swapping configurations found in {csv_file}")
        print("Available systems:", df['System'].unique())
        return

    # Calculate mean accuracy for each system across all datasets
    accuracy_by_system = df_roles.groupby('System')['Accuracy_%'].mean()

    # Create 3×3 matrix: rows=GPT role, columns=Claude role
    roles = ['Statistics', 'Biology', 'Pipeline']
    matrix = np.full((3, 3), np.nan)
    annotations = [['' for _ in range(3)] for _ in range(3)]

    role_to_idx = {'Statistics': 0, 'Biology': 1, 'Pipeline': 2}

    for system, accuracy in accuracy_by_system.items():
        config = parse_role_configuration(system)
        if config:
            gpt_idx = role_to_idx[config['gpt']]
            claude_idx = role_to_idx[config['claude']]
            matrix[gpt_idx, claude_idx] = accuracy
            annotations[gpt_idx][claude_idx] = f"{accuracy:.1f}%\nGemini={config['gemini']}"

    # Convert to DataFrame for better labeling
    matrix_df = pd.DataFrame(
        matrix,
        index=[f'GPT:\n{role}' for role in roles],
        columns=[f'Claude:\n{role}' for role in roles]
    )

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create heatmap with red-yellow-green colormap
    sns.heatmap(
        matrix_df,
        annot=np.array(annotations),
        fmt='',
        cmap='RdYlGn',
        vmin=70,
        vmax=95,
        cbar_kws={'label': 'Accuracy (%)', 'shrink': 0.8},
        linewidths=2,
        linecolor='white',
        ax=ax,
        cbar=True,
        square=True
    )

    # Formatting
    title = f'Role-Swapping Performance Matrix ({modality})\nAll 6 Role Permutations (Parallel Mode)'
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Claude Role Assignment', fontsize=12, fontweight='bold')
    ax.set_ylabel('GPT Role Assignment', fontsize=12, fontweight='bold')

    # Rotate labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='center', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)

    plt.tight_layout()

    # Save figure
    output_prefix = f'G6_role_swapping_matrix_{modality.lower().replace(" ", "_")}'
    png_path = os.path.join(output_dir, f'{output_prefix}.png')
    pdf_path = os.path.join(output_dir, f'{output_prefix}.pdf')

    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")

    plt.close()

    # Print summary statistics
    print(f"\n{modality} Role-Swapping Summary:")
    print(f"  Total configurations: {len(accuracy_by_system)}")
    print(f"  Mean accuracy across all permutations: {accuracy_by_system.mean():.1f}%")
    print(f"  Best configuration: {accuracy_by_system.idxmax().replace('parallel_', '').replace('_', ' ').upper()}")
    print(f"    - Accuracy: {accuracy_by_system.max():.1f}%")
    print(f"  Worst configuration: {accuracy_by_system.idxmin().replace('parallel_', '').replace('_', ' ').upper()}")
    print(f"    - Accuracy: {accuracy_by_system.min():.1f}%")

    # Additional analysis by role
    print("\n  Performance by GPT role assignment:")
    for role in roles:
        role_configs = [s for s in role_systems if parse_role_configuration(s)['gpt'] == role]
        role_acc = df_roles[df_roles['System'].isin(role_configs)]['Accuracy_%'].mean()
        print(f"    GPT={role}: {role_acc:.1f}%")

    print("\n  Performance by Claude role assignment:")
    for role in roles:
        role_configs = [s for s in role_systems if parse_role_configuration(s)['claude'] == role]
        role_acc = df_roles[df_roles['System'].isin(role_configs)]['Accuracy_%'].mean()
        print(f"    Claude={role}: {role_acc:.1f}%")

    print("\n  Performance by Gemini role assignment:")
    for role in roles:
        role_configs = [s for s in role_systems if parse_role_configuration(s)['gemini'] == role]
        role_acc = df_roles[df_roles['System'].isin(role_configs)]['Accuracy_%'].mean()
        print(f"    Gemini={role}: {role_acc:.1f}%")

def generate_alternative_bar_chart(csv_file, output_dir, modality):
    """Generate alternative bar chart showing all 6 permutations side-by-side"""

    print(f"\nGenerating alternative bar chart for {modality}...")
    df = pd.read_csv(csv_file)

    # Filter to only parallel role-swapping configurations
    role_systems = [
        'parallel_gp_st_cl_pl_gm_bl',
        'parallel_gp_st_cl_bl_gm_pl',
        'parallel_gp_pl_cl_st_gm_bl',
        'parallel_gp_pl_cl_bl_gm_st',
        'parallel_gp_bl_cl_st_gm_pl',
        'parallel_gp_bl_cl_pl_gm_st'
    ]

    df_roles = df[df['System'].isin(role_systems)]
    accuracy_by_system = df_roles.groupby('System')['Accuracy_%'].mean().sort_values(ascending=False)

    # Create readable labels
    labels = []
    for system in accuracy_by_system.index:
        config = parse_role_configuration(system)
        label = f"GPT={config['gpt'][:4]}\nClaude={config['claude'][:4]}\nGemini={config['gemini'][:4]}"
        labels.append(label)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))

    colors = ['#2ecc71' if acc >= 85 else '#f39c12' if acc >= 80 else '#e74c3c'
              for acc in accuracy_by_system.values]

    bars = ax.bar(range(len(accuracy_by_system)), accuracy_by_system.values, color=colors, alpha=0.8, edgecolor='black')

    # Add value labels on bars
    for i, (bar, acc) in enumerate(zip(bars, accuracy_by_system.values)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{acc:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Formatting
    ax.set_xticks(range(len(accuracy_by_system)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Accuracy (%)', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_title(f'Role-Swapping Performance Comparison ({modality})\nAll 6 Role Permutations',
                 fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()

    # Save figure
    output_prefix = f'G6_role_swapping_barchart_{modality.lower().replace(" ", "_")}'
    png_path = os.path.join(output_dir, f'{output_prefix}.png')
    pdf_path = os.path.join(output_dir, f'{output_prefix}.pdf')

    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")

    plt.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python generate_G6_role_swapping_heatmap.py <csv_file>")
        print("\nExample:")
        print("  python bin/generate_G6_role_swapping_heatmap.py experiments/bulk_rna_csv_figures/bulk_evaluation_per_experiment.csv")
        print("  python bin/generate_G6_role_swapping_heatmap.py experiments/scrna_csv_figures/scrna_evaluation_per_experiment.csv")
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

    # Generate both heatmap and bar chart
    generate_role_swapping_heatmap(csv_file, output_dir, modality)
    generate_alternative_bar_chart(csv_file, output_dir, modality)
    print("\nDone!")
