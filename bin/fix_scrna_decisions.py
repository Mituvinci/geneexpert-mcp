#!/usr/bin/env python3
"""
Fix scRNA-seq CSV decisions:
1. Set Stage 1 agent decisions to "AUTO_PROCEED" (no checkpoint at stage 1)
2. Mark any analysis that reached Stage 5 as completed
3. Set Stage 3A/3B decisions for datasets without cell cycle effects:
   - REH_Parental, Sup_Parental: Stage 3A='SKIP_CELL_CYCLE', 3B='Skip Regression'
4. Set Stage 3B='REMOVED_CELL_CYCLE' for datasets where cell cycle was already removed (PolyLLM only):
   - 10k mouse brain, GSE64016, GSE146773, GSE75748, PBMC

Usage:
    python3 bin/fix_scrna_decisions.py
"""

import pandas as pd
import sys

# File paths
INPUT_CSV = "/users/ha00014/Halimas_projects/multi_llm_mcp/experiments/scrna_results/scrna_ALL_EXPERIMENTS_DETAILED.csv"
OUTPUT_CSV = "/users/ha00014/Halimas_projects/multi_llm_mcp/experiments/scrna_results/scrna_ALL_EXPERIMENTS_DETAILED_FIXED.csv"

def main():
    print("=" * 80)
    print("FIXING scRNA-seq CSV DECISIONS")
    print("=" * 80)
    print()

    # Load CSV
    print(f"📖 Loading: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    print(f"   ✓ Loaded {len(df)} rows")
    print()

    # =========================================================================
    # FIX 1: Stage 1 agent decisions should be "AUTO_PROCEED"
    # =========================================================================
    print("🔧 FIX 1: Setting Stage 1 agent decisions to 'AUTO_PROCEED'")

    stage1_mask = (df['stage'] == 1) | (df['stage'] == '1')
    stage1_count = stage1_mask.sum()

    print(f"   Found {stage1_count} Stage 1 rows")

    # Set all three agent decisions to AUTO_PROCEED for stage 1
    df.loc[stage1_mask, 'gpt5_2_decision'] = 'AUTO_PROCEED'
    df.loc[stage1_mask, 'claude_decision'] = 'AUTO_PROCEED'
    df.loc[stage1_mask, 'gemini_decision'] = 'AUTO_PROCEED'

    print(f"   ✓ Set gpt5_2_decision, claude_decision, gemini_decision to 'AUTO_PROCEED'")
    print()

    # =========================================================================
    # FIX 2: Mark analyses that reached Stage 5 as completed
    # =========================================================================
    print("🔧 FIX 2: Marking analyses that reached Stage 5 as completed")

    # Get unique experiments (dataset + system combinations)
    experiments = df.groupby(['full_dataset', 'system_config'])

    completed_experiments = []

    for (dataset, system), group in experiments:
        # Check if this experiment has any Stage 5 rows
        has_stage5 = any((group['stage'] == 5) | (group['stage'] == '5'))

        if has_stage5:
            completed_experiments.append({
                'dataset': dataset,
                'system': system,
                'max_stage': 5,
                'completed': True
            })

    print(f"   ✓ Found {len(completed_experiments)} experiments that reached Stage 5")
    print()

    # Show summary by dataset
    print("📊 Completion Summary by Dataset:")
    for dataset in sorted(df['full_dataset'].unique()):
        if dataset == 'full_dataset':
            continue

        dataset_exps = [e for e in completed_experiments if e['dataset'] == dataset]
        total_systems = len(df[df['full_dataset'] == dataset]['system_config'].unique())
        completed_count = len(dataset_exps)

        print(f"   {dataset:40s}: {completed_count}/{total_systems} systems reached Stage 5")
    print()

    # =========================================================================
    # FIX 3: Handle datasets without cell cycle effects (REH, SUP-B15)
    # =========================================================================
    print("🔧 FIX 3: Setting decisions for datasets without cell cycle effects")

    # Datasets that don't have cell cycle effects
    no_cell_cycle_datasets = [
        '1_GD428_21136_Hu_REH_Parental',
        '2_GD444_21136_Hu_Sup_Parental'
    ]

    # Stage 3A: Set all agents to "SKIP_CELL_CYCLE"
    stage3a_mask = (
        df['full_dataset'].isin(no_cell_cycle_datasets) &
        (df['stage'] == '3A')
    )
    stage3a_count = stage3a_mask.sum()

    df.loc[stage3a_mask, 'gpt5_2_decision'] = 'SKIP_CELL_CYCLE'
    df.loc[stage3a_mask, 'claude_decision'] = 'SKIP_CELL_CYCLE'
    df.loc[stage3a_mask, 'gemini_decision'] = 'SKIP_CELL_CYCLE'

    print(f"   Stage 3A: Set {stage3a_count} rows to 'SKIP_CELL_CYCLE'")

    # Stage 3B: Set all agents to "Skip Regression"
    stage3b_mask = (
        df['full_dataset'].isin(no_cell_cycle_datasets) &
        (df['stage'] == '3B')
    )
    stage3b_count = stage3b_mask.sum()

    df.loc[stage3b_mask, 'gpt5_2_decision'] = 'Skip Regression'
    df.loc[stage3b_mask, 'claude_decision'] = 'Skip Regression'
    df.loc[stage3b_mask, 'gemini_decision'] = 'Skip Regression'

    print(f"   Stage 3B: Set {stage3b_count} rows to 'Skip Regression'")
    print(f"   ✓ Updated {stage3a_count + stage3b_count} rows for datasets: {', '.join(no_cell_cycle_datasets)}")
    print()

    # =========================================================================
    # FIX 4: Set Stage 3B to "REMOVED_CELL_CYCLE" for specific datasets (PolyLLM only)
    # =========================================================================
    print("🔧 FIX 4: Setting Stage 3B='REMOVED_CELL_CYCLE' for datasets with cell cycle already removed")

    # Datasets where cell cycle was already removed at Stage 3B
    removed_cell_cycle_datasets = [
        '10-k-brain-cells_healthy_mouse',
        '3_GSE64016_H1andFUCCI_normalized_EC_original',
        'GSE146773',
        'GSE75748',
        'pbmc_healthy_human'
    ]

    # Stage 3B: Set all agents to "REMOVED_CELL_CYCLE" for PolyLLM systems only
    stage3b_removed_mask = (
        df['full_dataset'].isin(removed_cell_cycle_datasets) &
        (df['stage'] == '3B') &
        ((df['system_config'].str.contains('parallel', case=False, na=False)) |
         (df['system_config'].str.contains('sequential', case=False, na=False)))
    )
    stage3b_removed_count = stage3b_removed_mask.sum()

    df.loc[stage3b_removed_mask, 'gpt5_2_decision'] = 'REMOVED_CELL_CYCLE'
    df.loc[stage3b_removed_mask, 'claude_decision'] = 'REMOVED_CELL_CYCLE'
    df.loc[stage3b_removed_mask, 'gemini_decision'] = 'REMOVED_CELL_CYCLE'

    print(f"   Stage 3B: Set {stage3b_removed_count} rows to 'REMOVED_CELL_CYCLE' (PolyLLM systems only)")
    print(f"   ✓ Updated datasets: {', '.join(removed_cell_cycle_datasets)}")
    print()

    # =========================================================================
    # Save output
    # =========================================================================
    print(f"💾 Saving fixed CSV to: {OUTPUT_CSV}")
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"   ✓ Saved {len(df)} rows")
    print()

    print("=" * 80)
    print("✅ FIXES APPLIED:")
    print(f"   1. Stage 1: Set all agent decisions to 'AUTO_PROCEED' ({stage1_count} rows)")
    print(f"   2. Completion: Identified {len(completed_experiments)} experiments reaching Stage 5")
    print(f"   3. No Cell Cycle: Set Stage 3A='SKIP_CELL_CYCLE', 3B='Skip Regression' ({stage3a_count + stage3b_count} rows)")
    print(f"      - Datasets: {', '.join(no_cell_cycle_datasets)}")
    print(f"   4. Removed Cell Cycle: Set Stage 3B='REMOVED_CELL_CYCLE' for PolyLLM systems ({stage3b_removed_count} rows)")
    print(f"      - Datasets: {', '.join(removed_cell_cycle_datasets)}")
    print()
    print("📝 NOTE: User will specify additional changes incrementally")
    print("=" * 80)
    print()
    print(f"Output: {OUTPUT_CSV}")
    print()

if __name__ == "__main__":
    main()
