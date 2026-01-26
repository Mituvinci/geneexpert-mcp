#!/usr/bin/env python3
"""
Aggregate all experiment decision metrics into comprehensive CSV files
for presentation to PI and analysis.
"""

import pandas as pd
import glob
import os
from pathlib import Path

# Base directory for results
RESULTS_DIR = "/users/ha00014/Halimas_projects/multi_llm_mcp/experiments/results/completed_result"
OUTPUT_DIR = "/users/ha00014/Halimas_projects/multi_llm_mcp/experiments/results/completed_result"

def extract_experiment_info(filepath):
    """Extract dataset and system from filepath"""
    folder_name = Path(filepath).parent.name
    parts = folder_name.split('_')

    # Format: 1_GSE52778_pe_clean_parallel
    dataset_num = parts[0]
    dataset_id = parts[1]
    sequencing = parts[2]
    quality = parts[3]
    system = '_'.join(parts[4:])

    return {
        'dataset_num': dataset_num,
        'dataset_id': dataset_id,
        'sequencing_type': sequencing,
        'dataset_quality': quality,
        'system_config': system,
        'full_dataset': f"{dataset_num}_{dataset_id}_{sequencing}_{quality}"
    }

def main():
    print("=" * 70)
    print("AGGREGATING ALL EXPERIMENT DECISION METRICS")
    print("=" * 70)
    print()

    # Find all decision metrics CSV files
    csv_pattern = os.path.join(RESULTS_DIR, "*", "staged_analysis_agent_decisions_metrics.csv")
    csv_files = sorted(glob.glob(csv_pattern))

    print(f"Found {len(csv_files)} experiment files:")
    for f in csv_files:
        print(f"  - {Path(f).parent.name}")
    print()

    # Read and aggregate all CSVs
    all_data = []

    for csv_file in csv_files:
        exp_info = extract_experiment_info(csv_file)
        df = pd.read_csv(csv_file)

        # Add experiment metadata columns
        for key, value in exp_info.items():
            df.insert(0, key, value)

        all_data.append(df)

    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)

    # Save full detailed CSV
    detailed_output = os.path.join(OUTPUT_DIR, "ALL_EXPERIMENTS_DETAILED.csv")
    combined_df.to_csv(detailed_output, index=False)
    print(f"✓ Saved detailed data: {detailed_output}")
    print(f"  Total rows: {len(combined_df)}")
    print()

    # Create summary view for PI (high-level metrics)
    summary_data = []

    for csv_file in csv_files:
        exp_info = extract_experiment_info(csv_file)
        df = pd.read_csv(csv_file)

        # Handle empty dataframes (no-agent baseline)
        if df.empty:
            summary_data.append({
                'Dataset': exp_info['full_dataset'],
                'System': exp_info['system_config'],
                'Organism': 'N/A',
                'Total_Stages': 0,
                'Stages_Passed': 0,
                'Stages_with_Issues': 0,
                'Stages_Failed': 0,
                'User_Inputs_Required': 0,
                'Unanimous_Approve': 0,
                'Unanimous_Reject': 0,
                'Majority_2of3': 0,
                'Split_Decision': 0,
                'Total_Cost_USD': 0.0,
                'Duration_Seconds': 0.0,
                'Final_Stage': 'N/A',
                'Final_Status': 'EMPTY',
                'Completed_Successfully': False
            })
            continue

        # Calculate summary metrics
        total_stages = len(df)
        stages_passed = df['stage_status'].isin(['PASS', 'SUCCESS']).sum()
        stages_issues = df['stage_status'].str.contains('ISSUES', na=False).sum()
        stages_failed = df['stage_status'].isin(['FAILED', 'ABORT']).sum()

        user_inputs_required = df['user_input_required'].sum()

        # Agent agreement (stages where all 3 agents agreed)
        if 'votes_approve' in df.columns:
            unanimous_approve = (df['votes_approve'] == 3).sum()
            unanimous_reject = (df['votes_reject'] == 3).sum()
            majority_2of3 = ((df['votes_approve'] == 2) | (df['votes_reject'] == 2)).sum()
            split_decision = ((df['votes_approve'] == 1) & (df['votes_reject'] == 1) & (df['votes_uncertain'] == 1)).sum()
        else:
            unanimous_approve = 0
            unanimous_reject = 0
            majority_2of3 = 0
            split_decision = 0

        # Total cost
        total_cost = df['stage_total_cost_usd'].sum()

        # Duration
        duration = df['duration_seconds'].iloc[0] if 'duration_seconds' in df.columns and not df.empty else 0

        # Final outcome
        final_row = df.iloc[-1]
        final_stage = final_row['stage_name']
        final_status = final_row['stage_status']
        completed_successfully = final_status in ['SUCCESS', 'PASS']

        summary_data.append({
            'Dataset': exp_info['full_dataset'],
            'System': exp_info['system_config'],
            'Organism': df['organism'].iloc[0] if 'organism' in df.columns and not df.empty else 'N/A',
            'Total_Stages': total_stages,
            'Stages_Passed': stages_passed,
            'Stages_with_Issues': stages_issues,
            'Stages_Failed': stages_failed,
            'User_Inputs_Required': user_inputs_required,
            'Unanimous_Approve': unanimous_approve,
            'Unanimous_Reject': unanimous_reject,
            'Majority_2of3': majority_2of3,
            'Split_Decision': split_decision,
            'Total_Cost_USD': round(total_cost, 4),
            'Duration_Seconds': round(duration, 2),
            'Final_Stage': final_stage,
            'Final_Status': final_status,
            'Completed_Successfully': completed_successfully
        })

    summary_df = pd.DataFrame(summary_data)

    # Save summary CSV
    summary_output = os.path.join(OUTPUT_DIR, "ALL_EXPERIMENTS_SUMMARY.csv")
    summary_df.to_csv(summary_output, index=False)
    print(f"✓ Saved summary for PI: {summary_output}")
    print(f"  Total experiments: {len(summary_df)}")
    print()

    # Print quick stats for PI update
    print("=" * 70)
    print("QUICK STATS FOR PI UPDATE")
    print("=" * 70)
    print()

    print(f"Total Experiments Run: {len(summary_df)}")
    print(f"  - Datasets: {summary_df['Dataset'].nunique()}")
    print(f"  - System Configurations: {summary_df['System'].nunique()}")
    print()

    print("Success Rate by System:")
    success_by_system = summary_df.groupby('System')['Completed_Successfully'].agg(['sum', 'count'])
    success_by_system['rate'] = (success_by_system['sum'] / success_by_system['count'] * 100).round(1)
    for system, row in success_by_system.iterrows():
        print(f"  {system:20s}: {row['sum']}/{row['count']} ({row['rate']}%)")
    print()

    print("Cost Analysis:")
    cost_by_system = summary_df.groupby('System')['Total_Cost_USD'].agg(['mean', 'sum'])
    for system, row in cost_by_system.iterrows():
        print(f"  {system:20s}: ${row['mean']:.4f}/experiment (total: ${row['sum']:.4f})")
    print()

    print("Agent Agreement Rates (Multi-Agent Systems Only):")
    multi_agent = summary_df[summary_df['System'].str.contains('parallel|sequential', case=False)]
    if not multi_agent.empty:
        total_decisions = multi_agent['Total_Stages'].sum()
        unanimous = multi_agent['Unanimous_Approve'].sum() + multi_agent['Unanimous_Reject'].sum()
        majority = multi_agent['Majority_2of3'].sum()
        split = multi_agent['Split_Decision'].sum()

        print(f"  Total Decisions: {total_decisions}")
        print(f"  Unanimous (3/3): {unanimous} ({unanimous/total_decisions*100:.1f}%)")
        print(f"  Majority (2/3):  {majority} ({majority/total_decisions*100:.1f}%)")
        print(f"  Split (1/1/1):   {split} ({split/total_decisions*100:.1f}%)")
    print()

    print("User Intervention Required:")
    print(f"  Total across all experiments: {summary_df['User_Inputs_Required'].sum()}")
    print(f"  Average per experiment: {summary_df['User_Inputs_Required'].mean():.2f}")
    print()

    # Show detailed comparison for each dataset
    print("=" * 70)
    print("PER-DATASET SYSTEM COMPARISON")
    print("=" * 70)
    print()

    for dataset in sorted(summary_df['Dataset'].unique()):
        print(f"Dataset: {dataset}")
        dataset_data = summary_df[summary_df['Dataset'] == dataset]
        for _, row in dataset_data.iterrows():
            status_icon = "✓" if row['Completed_Successfully'] else "✗"
            print(f"  {status_icon} {row['System']:20s}: {row['Stages_Passed']}/{row['Total_Stages']} stages | "
                  f"Cost: ${row['Total_Cost_USD']:.4f} | "
                  f"User inputs: {row['User_Inputs_Required']}")
        print()

    print("=" * 70)
    print("FILES GENERATED:")
    print(f"  1. {detailed_output}")
    print(f"     → Full data with all agent decisions, reasoning, costs")
    print(f"  2. {summary_output}")
    print(f"     → High-level summary for PI presentation")
    print("=" * 70)

if __name__ == "__main__":
    main()
