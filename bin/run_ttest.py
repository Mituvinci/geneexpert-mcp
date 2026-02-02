#!/usr/bin/env python3
"""
T-test Analysis for Bulk RNA-seq Multi-Agent Evaluation

Compares decision accuracy between different system types:
- Multi-foundational multi-agent (default_parallel)
- Single-LLM multi-agent (single_agent_gpt, single_agent_claude, single_agent_gemini)
- No-agent baseline (no_agent)

Usage:
    python3 bin/run_ttest.py bulk_evaluation_per_experiment.csv
"""

import sys
import pandas as pd
from scipy import stats
import numpy as np

def run_ttest_analysis(csv_file):
    """Run t-tests comparing multi-agent vs single-agent systems"""

    print("\n" + "="*80)
    print("BULK RNA-SEQ MULTI-AGENT T-TEST ANALYSIS")
    print("="*80 + "\n")

    # Load data
    df = pd.read_csv(csv_file)

    print(f"Loaded {len(df)} experiments\n")

    # Define system groups
    multi_foundational_systems = ['default_parallel', 'default_sequential',
                                   'role_perm_1', 'role_perm_3', 'role_perm_4',
                                   'role_perm_5', 'role_perm_6']

    single_llm_systems = ['single_agent_gpt', 'single_agent_claude', 'single_agent_gemini']

    # Extract accuracy scores
    multi_foundational_scores = df[df['System'].isin(multi_foundational_systems)]['Accuracy_%'].values
    single_llm_scores = df[df['System'].isin(single_llm_systems)]['Accuracy_%'].values
    no_agent_scores = df[df['System'] == 'no_agent']['Accuracy_%'].values

    # Print descriptive statistics
    print("DESCRIPTIVE STATISTICS")
    print("-" * 80)

    print(f"\n1. Multi-foundational Multi-Agent (n={len(multi_foundational_scores)}):")
    print(f"   Mean: {np.mean(multi_foundational_scores):.2f}%")
    print(f"   SD:   {np.std(multi_foundational_scores, ddof=1):.2f}%")
    print(f"   Min:  {np.min(multi_foundational_scores):.2f}%")
    print(f"   Max:  {np.max(multi_foundational_scores):.2f}%")

    print(f"\n2. Single-LLM Multi-Agent (n={len(single_llm_scores)}):")
    print(f"   Mean: {np.mean(single_llm_scores):.2f}%")
    print(f"   SD:   {np.std(single_llm_scores, ddof=1):.2f}%")
    print(f"   Min:  {np.min(single_llm_scores):.2f}%")
    print(f"   Max:  {np.max(single_llm_scores):.2f}%")

    print(f"\n3. No-Agent Baseline (n={len(no_agent_scores)}):")
    print(f"   Mean: {np.mean(no_agent_scores):.2f}%")
    print(f"   SD:   {np.std(no_agent_scores, ddof=1):.2f}%")

    # Run t-tests
    print("\n\n" + "="*80)
    print("T-TEST RESULTS")
    print("="*80)

    # Test 1: Multi-foundational vs Single-LLM
    print("\n1. Multi-foundational Multi-Agent vs Single-LLM Multi-Agent")
    print("-" * 80)

    t_stat_1, p_value_1 = stats.ttest_ind(multi_foundational_scores, single_llm_scores)

    print(f"   Null Hypothesis: No difference in decision accuracy")
    print(f"   Alternative: Multi-foundational has different accuracy than single-LLM")
    print(f"   t-statistic: {t_stat_1:.4f}")
    print(f"   p-value:     {p_value_1:.6f}")

    if p_value_1 < 0.05:
        print(f"   ✓ SIGNIFICANT (p < 0.05)")
        if np.mean(multi_foundational_scores) > np.mean(single_llm_scores):
            print(f"   → Multi-foundational is BETTER than single-LLM")
        else:
            print(f"   → Single-LLM is BETTER than multi-foundational")
    else:
        print(f"   ✗ NOT SIGNIFICANT (p ≥ 0.05)")
        print(f"   → No statistically significant difference")

    # Test 2: Multi-foundational vs No-Agent
    print("\n2. Multi-foundational Multi-Agent vs No-Agent Baseline")
    print("-" * 80)

    t_stat_2, p_value_2 = stats.ttest_ind(multi_foundational_scores, no_agent_scores)

    print(f"   Null Hypothesis: No difference in decision accuracy")
    print(f"   Alternative: Multi-foundational has different accuracy than no-agent")
    print(f"   t-statistic: {t_stat_2:.4f}")
    print(f"   p-value:     {p_value_2:.6f}")

    if p_value_2 < 0.05:
        print(f"   ✓ SIGNIFICANT (p < 0.05)")
        improvement = np.mean(multi_foundational_scores) - np.mean(no_agent_scores)
        print(f"   → Multi-foundational is {improvement:.1f}% more accurate than no-agent baseline")
    else:
        print(f"   ✗ NOT SIGNIFICANT (p ≥ 0.05)")

    # Test 3: Single-LLM vs No-Agent
    print("\n3. Single-LLM Multi-Agent vs No-Agent Baseline")
    print("-" * 80)

    t_stat_3, p_value_3 = stats.ttest_ind(single_llm_scores, no_agent_scores)

    print(f"   Null Hypothesis: No difference in decision accuracy")
    print(f"   Alternative: Single-LLM has different accuracy than no-agent")
    print(f"   t-statistic: {t_stat_3:.4f}")
    print(f"   p-value:     {p_value_3:.6f}")

    if p_value_3 < 0.05:
        print(f"   ✓ SIGNIFICANT (p < 0.05)")
        improvement = np.mean(single_llm_scores) - np.mean(no_agent_scores)
        print(f"   → Single-LLM is {improvement:.1f}% more accurate than no-agent baseline")
    else:
        print(f"   ✗ NOT SIGNIFICANT (p ≥ 0.05)")

    # Per-system breakdown
    print("\n\n" + "="*80)
    print("PER-SYSTEM BREAKDOWN")
    print("="*80 + "\n")

    for system in sorted(df['System'].unique()):
        system_scores = df[df['System'] == system]['Accuracy_%'].values
        n = len(system_scores)
        mean = np.mean(system_scores)
        sd = np.std(system_scores, ddof=1) if n > 1 else 0

        print(f"{system:25} n={n:3}   Mean: {mean:5.1f}%   SD: {sd:5.2f}%")

    # Completion rate analysis
    print("\n\n" + "="*80)
    print("COMPLETION RATE ANALYSIS")
    print("="*80 + "\n")

    for group_name, systems in [('Multi-foundational', multi_foundational_systems),
                                 ('Single-LLM', single_llm_systems),
                                 ('No-agent', ['no_agent'])]:
        group_df = df[df['System'].isin(systems)]
        completion_rate = (group_df['Completed'].sum() / len(group_df)) * 100
        print(f"{group_name:25} Completion Rate: {completion_rate:5.1f}%")

    print("\n" + "="*80 + "\n")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 bin/run_ttest.py bulk_evaluation_per_experiment.csv")
        sys.exit(1)

    csv_file = sys.argv[1]
    run_ttest_analysis(csv_file)
