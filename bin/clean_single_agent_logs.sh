#!/bin/bash
#
# Clean up confusing multi-agent references from single-agent experiment logs
#
# Usage:
#   bash bin/clean_single_agent_logs.sh                    # Clean all single-agent folders
#   bash bin/clean_single_agent_logs.sh FOLDER_PATH       # Clean specific folder
#

python3 bin/clean_single_agent_logs.py "$@"
