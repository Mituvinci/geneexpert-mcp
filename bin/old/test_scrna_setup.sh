#!/bin/bash

#
# Test scRNA-seq Setup
# Quick verification that everything is ready
#

echo "============================================================"
echo "Testing scRNA-seq Pipeline Setup"
echo "============================================================"
echo ""

ERRORS=0

# Test 1: Check R version
echo "[1/5] Checking R version..."
R_VERSION=$(R --version 2>&1 | head -1 | grep -oP 'version \K[0-9.]+' || echo "0.0.0")
R_MAJOR=$(echo $R_VERSION | cut -d. -f1)
if [ "$R_MAJOR" -ge 4 ]; then
  echo "  ✓ R $R_VERSION (OK)"
else
  echo "  ✗ R $R_VERSION (Need R 4.0+)"
  ERRORS=$((ERRORS+1))
fi
echo ""

# Test 2: Check essential R packages
echo "[2/5] Checking R packages..."
Rscript -e "
  essential <- c('jsonlite', 'dplyr', 'Seurat')
  missing <- c()
  for (pkg in essential) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    } else {
      cat(sprintf('  ✓ %s\n', pkg))
    }
  }
  if (length(missing) > 0) {
    cat(sprintf('\n  ✗ Missing: %s\n', paste(missing, collapse=', ')))
    quit(status=1)
  }
" 2>&1
if [ $? -ne 0 ]; then
  ERRORS=$((ERRORS+1))
fi
echo ""

# Test 3: Check Node.js
echo "[3/5] Checking Node.js..."
if command -v node &> /dev/null; then
  NODE_VERSION=$(node --version)
  echo "  ✓ Node.js $NODE_VERSION"
else
  echo "  ✗ Node.js not found"
  ERRORS=$((ERRORS+1))
fi
echo ""

# Test 4: Check scRNA JS files
echo "[4/5] Checking scRNA JavaScript files..."
if [ -f "src/scrna_executor/scrna_executor.js" ]; then
  echo "  ✓ scRNA executor exists"
else
  echo "  ✗ scRNA executor missing"
  ERRORS=$((ERRORS+1))
fi

if [ -f "bin/scrna_geneexpert.js" ]; then
  echo "  ✓ scRNA CLI exists"
else
  echo "  ✗ scRNA CLI missing"
  ERRORS=$((ERRORS+1))
fi

STAGE_COUNT=$(ls src/scrna_stages/*.js 2>/dev/null | wc -l)
if [ "$STAGE_COUNT" -eq 5 ]; then
  echo "  ✓ All 5 stage modules exist"
else
  echo "  ✗ Stage modules incomplete ($STAGE_COUNT/5)"
  ERRORS=$((ERRORS+1))
fi
echo ""

# Test 5: Check test data
echo "[5/5] Checking test data..."
if [ -d "data/scRNA_data/10-k-brain-cells_healthy_mouse" ]; then
  if [ -f "data/scRNA_data/10-k-brain-cells_healthy_mouse/filtered_feature_bc_matrix.h5" ]; then
    echo "  ✓ Test dataset exists (.h5 file found)"
  elif [ -d "data/scRNA_data/10-k-brain-cells_healthy_mouse/filtered_feature_bc_matrix" ]; then
    echo "  ✓ Test dataset exists (matrix directory found)"
  else
    echo "  ⚠ Test dataset directory exists but no data files found"
  fi
else
  echo "  ⚠ Test dataset not found (optional)"
fi
echo ""

# Summary
echo "============================================================"
if [ $ERRORS -eq 0 ]; then
  echo "✅ SUCCESS: Setup complete! Ready to run scRNA-seq pipeline."
  echo "============================================================"
  echo ""
  echo "Try running:"
  echo "  node bin/scrna_geneexpert.js analyze \\"
  echo "    data/scRNA_data/10-k-brain-cells_healthy_mouse \\"
  echo "    --organism mouse --output results/scrna_test"
  echo ""
else
  echo "❌ ERRORS FOUND: $ERRORS issues detected"
  echo "============================================================"
  echo ""
  echo "To fix:"
  echo "  1. Install missing R packages manually (see commands above)"
  echo "  2. Check SCRNA_SETUP.md for full installation guide"
  echo ""
  exit 1
fi
