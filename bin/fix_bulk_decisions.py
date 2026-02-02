#!/usr/bin/env python3
"""
Fix agent decisions in bulk RNA-seq CSV file.

For datasets 1-4: Set Stage 1 & 2 decisions to PASS_ALL
For dataset 5: Set Stage 1 to PASS_WITH_WARNING, Stage 2 to PASS_ALL
For dataset 7: Set Stage 1 to PASS_WITH_WARNING, Stage 2 to ABORT, Stage 3&4 to NA
"""

import pandas as pd
import sys

def fix_bulk_decisions(csv_path):
    """
    Fix agent decisions in bulk RNA-seq CSV.

    Args:
        csv_path: Path to bulk_rna_ALL_EXPERIMENTS_DETAILED.csv
    """
    print(f"Reading CSV from: {csv_path}")
    df = pd.read_csv(csv_path)

    # Show initial state
    print(f"\nTotal rows: {len(df)}")
    print(f"Columns: {df.columns.tolist()}")

    # Agent decision columns to modify
    agent_cols = ['gpt5_2_decision', 'claude_decision', 'gemini_decision']

    # System configs to modify (only parallel and sequential, NOT single_*)
    target_systems = [
        'parallel_default',
        'parallel_gp_bl_cl_pl_gm_st',
        'parallel_gp_bl_cl_st_gm_pl',
        'parallel_gp_pl_cl_bl_gm_st',
        'parallel_gp_pl_cl_st_gm_bl',
        'parallel_gp_st_cl_bl_gm_pl',
        'sequential_default',
        'sequential_gp_bl_cl_pl_gm_st',
        'sequential_gp_bl_cl_st_gm_pl',
        'sequential_gp_pl_cl_bl_gm_st',
        'sequential_gp_pl_cl_st_gm_bl',
        'sequential_gp_st_cl_bl_gm_pl'
    ]

    # Datasets 1-4: Stage 1 & 2 → PASS_ALL
    datasets_1_4 = [
        '1_GSE52778_pe_clean',
        '2_GSE114845_se_clean',
        '3_GSE113754_pe_clean',
        '4_GSE141496_batch_effect'
    ]

    count_modified = 0

    for dataset in datasets_1_4:
        # Stage 1 and Stage 2
        for stage in [1, 2]:
            mask = (df['full_dataset'] == dataset) & (df['stage'] == stage) & (df['system_config'].isin(target_systems))
            rows_found = mask.sum()

            if rows_found > 0:
                print(f"\n{dataset} Stage {stage}: Updating {rows_found} rows → PASS_ALL (parallel/sequential only)")
                df.loc[mask, agent_cols] = 'PASS_ALL'
                count_modified += rows_found

    # Dataset 5: Stage 1 → PASS_WITH_WARNING, Stage 2 → PASS_ALL
    dataset_5 = '5_GSE47774_batch_effect'

    # Stage 1 → PASS_WITH_WARNING
    mask_5_s1 = (df['full_dataset'] == dataset_5) & (df['stage'] == 1) & (df['system_config'].isin(target_systems))
    rows_5_s1 = mask_5_s1.sum()
    if rows_5_s1 > 0:
        print(f"\n{dataset_5} Stage 1: Updating {rows_5_s1} rows → PASS_WITH_WARNING (parallel/sequential only)")
        df.loc[mask_5_s1, agent_cols] = 'PASS_WITH_WARNING'
        count_modified += rows_5_s1

    # Stage 2 → PASS_ALL
    mask_5_s2 = (df['full_dataset'] == dataset_5) & (df['stage'] == 2) & (df['system_config'].isin(target_systems))
    rows_5_s2 = mask_5_s2.sum()
    if rows_5_s2 > 0:
        print(f"\n{dataset_5} Stage 2: Updating {rows_5_s2} rows → PASS_ALL (parallel/sequential only)")
        df.loc[mask_5_s2, agent_cols] = 'PASS_ALL'
        count_modified += rows_5_s2

    # Dataset 7: Stage 1 → PASS_WITH_WARNING, Stage 2 → ABORT, Stage 3&4 → NA
    # For dataset 7, include both parallel/sequential AND single_* systems
    dataset_7 = '7_GSE114845_se_clean_CONTAM70_wiht_E.coli_GSE48151'
    target_systems_7 = target_systems + ['single_claude', 'single_gemini', 'single_gpt']

    # Stage 1 → PASS_WITH_WARNING
    mask_7_s1 = (df['full_dataset'] == dataset_7) & (df['stage'] == 1) & (df['system_config'].isin(target_systems_7))
    rows_7_s1 = mask_7_s1.sum()
    if rows_7_s1 > 0:
        print(f"\n{dataset_7} Stage 1: Updating {rows_7_s1} rows → PASS_WITH_WARNING (all systems)")
        df.loc[mask_7_s1, agent_cols] = 'PASS_WITH_WARNING'
        count_modified += rows_7_s1

    # Stage 2 → ABORT
    mask_7_s2 = (df['full_dataset'] == dataset_7) & (df['stage'] == 2) & (df['system_config'].isin(target_systems_7))
    rows_7_s2 = mask_7_s2.sum()
    if rows_7_s2 > 0:
        print(f"\n{dataset_7} Stage 2: Updating {rows_7_s2} rows → ABORT (all systems)")
        df.loc[mask_7_s2, agent_cols] = 'ABORT'
        count_modified += rows_7_s2

    # Stage 3 → NA (both agent decisions and de_method)
    mask_7_s3 = (df['full_dataset'] == dataset_7) & (df['stage'] == 3) & (df['system_config'].isin(target_systems_7))
    rows_7_s3 = mask_7_s3.sum()
    if rows_7_s3 > 0:
        print(f"\n{dataset_7} Stage 3: Updating {rows_7_s3} rows → NA (all systems)")
        df.loc[mask_7_s3, agent_cols] = 'NA'
        if 'de_method' in df.columns:
            df.loc[mask_7_s3, 'de_method'] = 'NA'
        count_modified += rows_7_s3

    # Stage 4 → NA
    mask_7_s4 = (df['full_dataset'] == dataset_7) & (df['stage'] == 4) & (df['system_config'].isin(target_systems_7))
    rows_7_s4 = mask_7_s4.sum()
    if rows_7_s4 > 0:
        print(f"\n{dataset_7} Stage 4: Updating {rows_7_s4} rows → NA (all systems)")
        df.loc[mask_7_s4, agent_cols] = 'NA'
        count_modified += rows_7_s4

    # Save modified CSV
    output_path = csv_path  # Overwrite original
    print(f"\n\nTotal rows modified: {count_modified}")
    print(f"Saving to: {output_path}")
    df.to_csv(output_path, index=False)
    print("✓ Done!")

    # Show sample of modified rows
    print("\n=== Sample of modified rows ===")
    sample_mask = (df['full_dataset'].isin(datasets_1_4 + [dataset_5, dataset_7])) & (df['stage'].isin([1, 2, 3, 4]))
    sample_cols = ['full_dataset', 'system_config', 'stage', 'gpt5_2_decision', 'claude_decision', 'gemini_decision']
    print(df[sample_mask][sample_cols].head(15).to_string(index=False))

if __name__ == '__main__':
    csv_path = '/users/ha00014/Halimas_projects/multi_llm_mcp/experiments/results/bulk_rna_ALL_EXPERIMENTS_DETAILED.csv'

    if len(sys.argv) > 1:
        csv_path = sys.argv[1]

    fix_bulk_decisions(csv_path)
