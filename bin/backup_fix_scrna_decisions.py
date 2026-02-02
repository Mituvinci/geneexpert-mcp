#!/usr/bin/env python3
"""
Fix scRNA-seq CSV decisions:
1. Set Stage 1 agent decisions to "AUTO_PROCEED" (no checkpoint at stage 1)
2. Mark any analysis that reached Stage 5 as completed
3. User will specify additional changes incrementally

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
    print()
    print("📝 NOTE: User will specify additional changes incrementally")
    print("=" * 80)
    print()
    print(f"Output: {OUTPUT_CSV}")
    print()

if __name__ == "__main__":
    main()
