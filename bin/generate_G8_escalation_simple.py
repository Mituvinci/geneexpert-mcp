#!/usr/bin/env python3
"""
Generate G8: User Escalation Frequency and Justification Rate

Uses provided statistics to show escalation patterns.
"""

import matplotlib.pyplot as plt
import numpy as np
import os

# Publication style
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
plt.rcParams['font.size'] = 11

def generate_escalation_plot(output_dir):
    """Generate escalation analysis bar chart"""

    # Data from analysis
    categories = ['PolyLLM', 'Single-LLM']
    escalation_rates = [18.2, 7.9]
    justification_rates = [83.4, 45.2]

    os.makedirs(output_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5))

    x = np.arange(len(categories))
    width = 0.35

    # Colors
    color_escalation = '#4393C3'  # Blue
    color_justification = '#2CA02C'  # Green

    # Plot bars
    bars1 = ax.bar(x - width/2, escalation_rates, width,
                   label='Escalation Rate', color=color_escalation,
                   edgecolor='black', linewidth=0.8)
    bars2 = ax.bar(x + width/2, justification_rates, width,
                   label='Justification Rate', color=color_justification,
                   edgecolor='black', linewidth=0.8)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{height:.1f}%', ha='center', va='bottom',
                   fontsize=10, fontweight='bold')

    # Labels and title
    ax.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    ax.set_xlabel('System Category', fontsize=12, fontweight='bold')
    ax.set_title('User Escalation Frequency and Justification Rate\nAveraged Across Bulk and scRNA-seq Pipelines',
                 fontsize=13, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=11)
    ax.set_ylim(0, 95)
    ax.legend(loc='upper right', fontsize=10, frameon=True, shadow=False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    # Save
    output_prefix = os.path.join(output_dir, 'G8_escalation_analysis')
    plt.savefig(f'{output_prefix}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved: {output_prefix}.pdf")
    print(f"Saved: {output_prefix}.png")

    # Summary
    print("\nEscalation Analysis Summary:")
    print(f"  PolyLLM:")
    print(f"    - Escalates {escalation_rates[0]:.1f}% of decisions to users")
    print(f"    - {justification_rates[0]:.1f}% of escalations are justified (led to correct decisions)")
    print(f"  Single-LLM:")
    print(f"    - Escalates {escalation_rates[1]:.1f}% of decisions to users")
    print(f"    - {justification_rates[1]:.1f}% of escalations are justified")
    print(f"\nKey Finding:")
    print(f"  PolyLLM escalates {escalation_rates[0] - escalation_rates[1]:.1f}% more often than Single-LLM,")
    print(f"  but achieves {justification_rates[0] - justification_rates[1]:.1f}% higher justification rate.")
    print(f"  This demonstrates that PolyLLM surfaces genuine uncertainty through")
    print(f"  principled escalation, while Single-LLM under-escalates and fails silently.")

if __name__ == '__main__':
    import sys
    output_dir = sys.argv[1] if len(sys.argv) > 1 else 'experiments/combined_figures'
    print(f"Generating G8 escalation analysis in {output_dir}...")
    generate_escalation_plot(output_dir)
    print("\nDone!")
