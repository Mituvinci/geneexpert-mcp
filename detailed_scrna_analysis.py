#!/usr/bin/env python3
"""
Detailed analysis of scRNA-seq experimental results
Shows stage-by-stage progress, cell cycle detection, and agent decisions
"""

import os
import json
from pathlib import Path
import re

RESULTS_DIR = Path("/users/ha00014/Halimas_projects/multi_llm_mcp/experiments/scrna_results")

def parse_directory_name(dirname):
    """Extract dataset name and architecture from directory name"""
    if dirname.endswith('_no_agent'):
        arch = 'no_agent'
        dataset = dirname.replace('_no_agent', '')
    elif dirname.endswith('_single_gpt'):
        arch = 'single_gpt'
        dataset = dirname.replace('_single_gpt', '')
    elif dirname.endswith('_single_claude'):
        arch = 'single_claude'
        dataset = dirname.replace('_single_claude', '')
    elif dirname.endswith('_parallel'):
        arch = 'parallel'
        dataset = dirname.replace('_parallel', '')
    elif dirname.endswith('_sequential'):
        arch = 'sequential'
        dataset = dirname.replace('_sequential', '')
    else:
        arch = 'unknown'
        dataset = dirname
    return dataset, arch

def check_stage_completion(result_dir):
    """Check which stages were completed"""
    stage_info = {
        'stage1': False,
        'stage2': False,
        'stage3': False,
        'stage4': False,
        'stage5': False,
        'cell_cycle_detected': False,
        'cell_cycle_corrected': False,
        'stages_completed': 0
    }

    stage_names = ['load_qc', 'filter_qc', 'normalize_hvg', 'pca', 'cluster_markers']

    for i, name in enumerate(stage_names, 1):
        stage_dir = result_dir / f"stage{i}_{name}"
        if stage_dir.exists():
            stage_info[f'stage{i}'] = True
            stage_info['stages_completed'] += 1

            # Check for cell cycle analysis in stage 3
            if i == 3:
                # Check for cell cycle files
                cc_before = stage_dir / "cell_cycle_before.pdf"
                cc_after = stage_dir / "cell_cycle_after.pdf"
                cc_summary = stage_dir / "cell_cycle_summary.txt"

                if cc_summary.exists():
                    stage_info['cell_cycle_detected'] = True
                    try:
                        with open(cc_summary, 'r') as f:
                            content = f.read()
                            if 'regress' in content.lower() or 'removed' in content.lower():
                                stage_info['cell_cycle_corrected'] = True
                    except:
                        pass

                if cc_before.exists() and cc_after.exists():
                    stage_info['cell_cycle_corrected'] = True

    return stage_info

def parse_log_file(log_file):
    """Extract key decisions from log file"""
    decisions = {
        'stage2_decision': 'N/A',
        'stage2_user_input': False,
        'stage4_decision': 'N/A',
        'stage4_user_input': False,
        'stage5_decision': 'N/A',
        'stage5_user_input': False,
        'reanalysis_requested': False,
        'errors': []
    }

    if not log_file.exists():
        return decisions

    try:
        with open(log_file, 'r') as f:
            content = f.read()

            # Check for user input required
            if 'User decision required' in content or 'user input required' in content.lower():
                decisions['stage2_user_input'] = True

            # Check for reanalysis requests
            if 'REANALYZ' in content or 'ADJUST' in content or 'STOP_AND_REVIEW' in content:
                decisions['reanalysis_requested'] = True

            # Look for stage 2 decisions
            if 'PROCEED' in content:
                decisions['stage2_decision'] = 'PROCEED'
            elif 'STOP_AND_REVIEW' in content:
                decisions['stage2_decision'] = 'STOP_AND_REVIEW'

            # Look for stage 4 decisions
            if 'USE_DEFAULT' in content:
                decisions['stage4_decision'] = 'USE_DEFAULT'
            elif 'SELECT_PC_RANGE' in content:
                decisions['stage4_decision'] = 'SELECT_PC_RANGE'

            # Look for stage 5 decisions
            if 'ACCEPT_CLUSTERING' in content:
                decisions['stage5_decision'] = 'ACCEPT_CLUSTERING'
            elif 'ADJUST_RESOLUTION' in content:
                decisions['stage5_decision'] = 'ADJUST_RESOLUTION'
            elif 'FLAG_SUSPICIOUS' in content:
                decisions['stage5_decision'] = 'FLAG_SUSPICIOUS'

            # Check for errors
            if 'ERROR' in content or 'Error' in content:
                decisions['errors'].append('Error found in log')

    except Exception as e:
        decisions['errors'].append(f'Error parsing log: {str(e)}')

    return decisions

def main():
    all_results = []

    # Scan all result directories
    for entry in sorted(RESULTS_DIR.iterdir()):
        if not entry.is_dir():
            continue

        dataset, arch = parse_directory_name(entry.name)
        stage_info = check_stage_completion(entry)

        # Parse log file
        log_file = entry / "scrna_analysis_log.txt"
        decisions = parse_log_file(log_file)

        # Check final completion
        stage5_dir = entry / "stage5_cluster_markers"
        completed = False
        if stage5_dir.exists():
            seurat_file = stage5_dir / "seurat_stage5_clustered.rds"
            if seurat_file.exists():
                completed = True

        all_results.append({
            'dataset': dataset,
            'architecture': arch,
            'completed': completed,
            'stages_completed': stage_info['stages_completed'],
            'stage1': stage_info['stage1'],
            'stage2': stage_info['stage2'],
            'stage3': stage_info['stage3'],
            'stage4': stage_info['stage4'],
            'stage5': stage_info['stage5'],
            'cell_cycle_detected': stage_info['cell_cycle_detected'],
            'cell_cycle_corrected': stage_info['cell_cycle_corrected'],
            'stage2_decision': decisions['stage2_decision'],
            'stage4_decision': decisions['stage4_decision'],
            'stage5_decision': decisions['stage5_decision'],
            'reanalysis_requested': decisions['reanalysis_requested'],
            'user_input': decisions['stage2_user_input'] or decisions['stage4_user_input'] or decisions['stage5_user_input']
        })

    # Print detailed table
    print("\n" + "="*180)
    print("Detailed scRNA-seq Analysis Results - Stage-by-Stage Progress")
    print("="*180)
    print(f"{'Dataset':<45} {'Arch':<12} {'S1':<4} {'S2':<4} {'S3':<4} {'S4':<4} {'S5':<4} {'CC-Det':<7} {'CC-Cor':<7} {'S2-Dec':<18} {'S4-Dec':<18} {'S5-Dec':<20} {'Complete':<10}")
    print("-"*180)

    for r in all_results:
        s1 = '✓' if r['stage1'] else '✗'
        s2 = '✓' if r['stage2'] else '✗'
        s3 = '✓' if r['stage3'] else '✗'
        s4 = '✓' if r['stage4'] else '✗'
        s5 = '✓' if r['stage5'] else '✗'
        cc_det = '✓' if r['cell_cycle_detected'] else '✗'
        cc_cor = '✓' if r['cell_cycle_corrected'] else '✗'
        complete = '✓ YES' if r['completed'] else '✗ NO'

        print(f"{r['dataset']:<45} {r['architecture']:<12} {s1:<4} {s2:<4} {s3:<4} {s4:<4} {s5:<4} {cc_det:<7} {cc_cor:<7} {r['stage2_decision']:<18} {r['stage4_decision']:<18} {r['stage5_decision']:<20} {complete:<10}")

    print("="*180)
    print("\nLegend:")
    print("  S1-S5: Stage 1-5 completion")
    print("  CC-Det: Cell Cycle Detected (Stage 3)")
    print("  CC-Cor: Cell Cycle Corrected (Stage 3)")
    print("  S2-Dec: Stage 2 QC Filtering Decision")
    print("  S4-Dec: Stage 4 PCA Decision")
    print("  S5-Dec: Stage 5 Clustering Decision")
    print("="*180)

    # Summary by architecture
    print("\n" + "="*120)
    print("Summary by Architecture")
    print("="*120)
    print(f"{'Architecture':<20} {'Total':<10} {'Completed':<12} {'Cell Cycle OK':<15} {'Success Rate':<15}")
    print("-"*120)

    architectures = ['no_agent', 'single_gpt', 'single_claude', 'parallel', 'sequential']
    for arch in architectures:
        arch_results = [r for r in all_results if r['architecture'] == arch]
        total = len(arch_results)
        if total == 0:
            continue
        completed = sum(1 for r in arch_results if r['completed'])
        cc_ok = sum(1 for r in arch_results if r['cell_cycle_corrected'])
        success_rate = f"{(completed/total*100):.1f}%"

        print(f"{arch:<20} {total:<10} {completed:<12} {cc_ok:<15} {success_rate:<15}")

    print("="*120)

    # Export detailed CSV
    csv_file = RESULTS_DIR / "detailed_results.csv"
    with open(csv_file, 'w') as f:
        f.write("dataset,architecture,completed,stages_completed,stage1,stage2,stage3,stage4,stage5,")
        f.write("cell_cycle_detected,cell_cycle_corrected,stage2_decision,stage4_decision,stage5_decision,")
        f.write("reanalysis_requested,user_input\n")

        for r in all_results:
            f.write(f"{r['dataset']},{r['architecture']},{r['completed']},{r['stages_completed']},")
            f.write(f"{r['stage1']},{r['stage2']},{r['stage3']},{r['stage4']},{r['stage5']},")
            f.write(f"{r['cell_cycle_detected']},{r['cell_cycle_corrected']},")
            f.write(f"{r['stage2_decision']},{r['stage4_decision']},{r['stage5_decision']},")
            f.write(f"{r['reanalysis_requested']},{r['user_input']}\n")

    print(f"\nDetailed results exported to: {csv_file}")

    # Dataset-specific summary
    datasets = sorted(set(r['dataset'] for r in all_results))
    dataset_summary_file = RESULTS_DIR / "dataset_summary.csv"
    with open(dataset_summary_file, 'w') as f:
        f.write("dataset,no_agent,single_gpt,single_claude,parallel,sequential,best_architecture\n")

        for dataset in datasets:
            dataset_results = [r for r in all_results if r['dataset'] == dataset]
            row = [dataset]

            arch_status = {}
            for arch in architectures:
                arch_result = next((r for r in dataset_results if r['architecture'] == arch), None)
                if arch_result:
                    status = 'COMPLETE' if arch_result['completed'] else f"{arch_result['stages_completed']}/5"
                    arch_status[arch] = status
                    row.append(status)
                else:
                    row.append('N/A')

            # Determine best architecture
            completed_archs = [arch for arch in architectures if arch in arch_status and arch_status[arch] == 'COMPLETE']
            if completed_archs:
                if 'parallel' in completed_archs and 'sequential' in completed_archs:
                    best = 'parallel+sequential'
                elif 'parallel' in completed_archs:
                    best = 'parallel'
                elif 'sequential' in completed_archs:
                    best = 'sequential'
                else:
                    best = completed_archs[0]
            else:
                # Find highest partial completion
                max_stages = 0
                best = 'none'
                for arch in architectures:
                    if arch in arch_status and '/' in arch_status[arch]:
                        stages = int(arch_status[arch].split('/')[0])
                        if stages > max_stages:
                            max_stages = stages
                            best = arch

            row.append(best)
            f.write(','.join(row) + '\n')

    print(f"Dataset summary exported to: {dataset_summary_file}\n")

if __name__ == "__main__":
    main()
