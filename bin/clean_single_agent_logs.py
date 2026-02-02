#!/usr/bin/env python3
"""
Clean up confusing multi-agent references from single-agent experiment logs.

This script removes:
1. ROLE SWAPPING messages (misleading in single-agent mode)
2. References to unused models (e.g., Claude/Gemini in GPT-only runs)
3. Confusing role assignment displays

Usage:
    python bin/clean_single_agent_logs.py
    python bin/clean_single_agent_logs.py --folder experiments/results/GSE52778_single_gpt
"""

import os
import re
import sys
import glob
from pathlib import Path

# Base directory for experiments
BASE_DIR = Path("/users/ha00014/Halimas_projects/multi_llm_mcp/experiments/results")

def detect_agent_type(folder_name):
    """Detect which single-agent type based on folder name."""
    if "_single_gpt" in folder_name:
        return "gpt"
    elif "_single_claude" in folder_name:
        return "claude"
    elif "_single_gemini" in folder_name:
        return "gemini"
    return None

def clean_log_content(content, agent_type):
    """
    Remove confusing multi-agent references from log content.

    Args:
        content: Original log file content
        agent_type: 'gpt', 'claude', or 'gemini'

    Returns:
        Cleaned log content
    """
    # Model name mappings
    model_names = {
        'gpt': 'GPT-5.2',
        'claude': 'Claude Opus',
        'gemini': 'Gemini Pro'
    }

    active_model = model_names.get(agent_type, 'Unknown')

    # Remove entire ROLE SWAPPING sections (multi-line)
    # Pattern: [ROLE SWAPPING] ... until next non-indented line
    content = re.sub(
        r'\[ROLE SWAPPING\].*?\n(?=\[)',
        '',
        content,
        flags=re.DOTALL
    )

    # Remove role assignment blocks
    content = re.sub(
        r'\[ROLE SWAPPING\] Current agent role assignments:.*?\n\n',
        '',
        content,
        flags=re.DOTALL
    )

    content = re.sub(
        r'\[ROLE SWAPPING\] Prompts assigned to each agent:.*?\n\n',
        '',
        content,
        flags=re.DOTALL
    )

    # Remove lines mentioning unused models
    if agent_type == 'gpt':
        # Remove Claude and Gemini references
        content = re.sub(r'^.*Claude.*→.*Agent.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^.*Gemini.*→.*Agent.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^.*Claude.*received:.*prompt.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^.*Gemini.*received:.*prompt.*\n', '', content, flags=re.MULTILINE)
        # Remove cost lines for unused models
        content = re.sub(r'^\s+Claude:\s+\$[\d.]+.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^\s+Gemini:\s+\$[\d.]+.*\n', '', content, flags=re.MULTILINE)

    elif agent_type == 'claude':
        # Remove GPT and Gemini references
        content = re.sub(r'^.*GPT-5\.2.*→.*Agent.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^.*Gemini.*→.*Agent.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^.*GPT-5\.2.*received:.*prompt.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^.*Gemini.*received:.*prompt.*\n', '', content, flags=re.MULTILINE)
        # Remove cost lines for unused models
        content = re.sub(r'^\s+GPT-5\.2:\s+\$[\d.]+.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^\s+Gemini:\s+\$[\d.]+.*\n', '', content, flags=re.MULTILINE)

    elif agent_type == 'gemini':
        # Remove GPT and Claude references
        content = re.sub(r'^.*GPT-5\.2.*→.*Agent.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^.*Claude.*→.*Agent.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^.*GPT-5\.2.*received:.*prompt.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^.*Claude.*received:.*prompt.*\n', '', content, flags=re.MULTILINE)
        # Remove cost lines for unused models
        content = re.sub(r'^\s+GPT-5\.2:\s+\$[\d.]+.*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'^\s+Claude:\s+\$[\d.]+.*\n', '', content, flags=re.MULTILINE)

    # Clean up empty line runs (more than 2 consecutive newlines)
    content = re.sub(r'\n{3,}', '\n\n', content)

    return content

def clean_log_file(log_path, agent_type):
    """Clean a single log file."""
    try:
        with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        original_size = len(content)
        cleaned_content = clean_log_content(content, agent_type)
        cleaned_size = len(cleaned_content)

        # Only write if changes were made
        if original_size != cleaned_size:
            with open(log_path, 'w', encoding='utf-8') as f:
                f.write(cleaned_content)

            chars_removed = original_size - cleaned_size
            print(f"  ✓ Cleaned {log_path.name}: removed {chars_removed} chars")
            return True
        else:
            print(f"  - {log_path.name}: no changes needed")
            return False

    except Exception as e:
        print(f"  ✗ Error cleaning {log_path}: {e}")
        return False

def clean_folder(folder_path):
    """Clean all log files in a single-agent folder."""
    folder_path = Path(folder_path)

    if not folder_path.exists():
        print(f"✗ Folder not found: {folder_path}")
        return 0

    # Detect agent type
    agent_type = detect_agent_type(folder_path.name)
    if not agent_type:
        print(f"✗ Not a single-agent folder: {folder_path.name}")
        return 0

    print(f"\n{'='*70}")
    print(f"Cleaning: {folder_path.name}")
    print(f"Agent Type: {agent_type.upper()}")
    print(f"{'='*70}")

    # Find log files
    log_patterns = [
        "*_log.txt",
        "*_agent_conversations.txt",
        "staged_analysis_log.txt",
        "staged_analysis_agent_conversations.txt"
    ]

    cleaned_count = 0
    for pattern in log_patterns:
        for log_file in folder_path.glob(pattern):
            if clean_log_file(log_file, agent_type):
                cleaned_count += 1

    print(f"\n✓ Cleaned {cleaned_count} files in {folder_path.name}")
    return cleaned_count

def main():
    """Main function."""
    print("="*70)
    print("Single-Agent Log Cleaner")
    print("="*70)
    print("\nRemoving confusing multi-agent references from single-agent logs...")

    # Check if specific folder provided
    if len(sys.argv) > 1:
        folder = sys.argv[1]
        if folder.startswith("--folder="):
            folder = folder.split("=", 1)[1]

        total_cleaned = clean_folder(folder)

    else:
        # Find all single-agent folders
        patterns = [
            "*_single_gpt",
            "*_single_claude",
            "*_single_gemini"
        ]

        folders_found = []
        for pattern in patterns:
            folders_found.extend(BASE_DIR.glob(pattern))

        if not folders_found:
            print(f"\n✗ No single-agent folders found in {BASE_DIR}")
            return

        print(f"\nFound {len(folders_found)} single-agent folders")

        total_cleaned = 0
        for folder in sorted(folders_found):
            total_cleaned += clean_folder(folder)

        print("\n" + "="*70)
        print(f"✓ COMPLETE: Cleaned {total_cleaned} log files total")
        print("="*70)

if __name__ == "__main__":
    main()
