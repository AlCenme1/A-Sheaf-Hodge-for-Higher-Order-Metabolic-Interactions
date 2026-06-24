#!/bin/bash
# setup_repo.sh
# -------------
# Run this script from the directory where you have all the files.
# It organizes everything into the correct folder structure and
# pushes to your GitHub repository in one shot.
#
# Usage:
#   1. Place all files in one folder alongside this script.
#   2. chmod +x setup_repo.sh
#   3. ./setup_repo.sh

set -e

REPO_URL="https://github.com/AlCenme1/HIGHER-ORDER-METABOLIC-OBSTRUCTIONS-IN-MICROBIAL-COMMUNITIES-SHEAF--HODGE-FRAMEWORK.git"
REPO_DIR="sheaf-hodge-repo"

echo "=== Setting up repository structure ==="

# Initialize git repo
mkdir -p "$REPO_DIR"
cd "$REPO_DIR"
git init
git remote add origin "$REPO_URL"

# Create folder structure
mkdir -p figures notebooks code latex

# ── Figures ──────────────────────────────────────────────────────────────────
echo "Copying figures..."
for f in ima1.png fig2_power_curve.png fig3_coupling_decomposition.png \
          fig4_interface_emergence.png figS1_qqplot.png figS2_power_multiD.png \
          figS3_locoh_detection.png figS4_complex.png figS5_greedy_design.png; do
    [ -f "../$f" ] && cp "../$f" figures/ && echo "  + figures/$f"
done

# ── Notebooks ────────────────────────────────────────────────────────────────
echo "Copying notebooks..."
for nb in coupling_decomposition.ipynb Real-data_application.ipynb; do
    [ -f "../$nb" ] && cp "../$nb" notebooks/ && echo "  + notebooks/$nb"
done

# ── Code ─────────────────────────────────────────────────────────────────────
echo "Copying code..."
for py in fig2_power_curve.py fig3_from_synthetic.py supp_figures.py; do
    [ -f "../$py" ] && cp "../$py" code/ && echo "  + code/$py"
done

# ── LaTeX ────────────────────────────────────────────────────────────────────
echo "Copying LaTeX..."
for tex in NewIEEE.tex SupplementaryMaterial.tex; do
    [ -f "../$tex" ] && cp "../$tex" latex/ && echo "  + latex/$tex"
done

# ── README ───────────────────────────────────────────────────────────────────
[ -f "../README.md" ] && cp "../README.md" . && echo "  + README.md"

# ── .gitignore ───────────────────────────────────────────────────────────────
cat > .gitignore << 'EOF'
# Python
__pycache__/
*.pyc
*.pyo
.ipynb_checkpoints/

# LaTeX auxiliary files
*.aux
*.log
*.out
*.toc
*.bbl
*.blg
*.synctex.gz
*.fdb_latexmk
*.fls
*.nav
*.snm
*.vrb

# OS
.DS_Store
Thumbs.db
EOF
echo "  + .gitignore"

# ── LICENSE ──────────────────────────────────────────────────────────────────
cat > LICENSE << 'EOF'
MIT License

Copyright (c) 2026 Alex Centeno-Mejía

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOF
echo "  + LICENSE"

# ── requirements.txt ─────────────────────────────────────────────────────────
cat > requirements.txt << 'EOF'
numpy>=1.24
scipy>=1.10
networkx>=3.0
matplotlib>=3.7
jupyter>=1.0
EOF
echo "  + requirements.txt"

# ── Commit and push ──────────────────────────────────────────────────────────
echo ""
echo "=== Committing and pushing ==="
git add .
git commit -m "Initial commit: paper, supplementary, figures, notebooks, and code

- NewIEEE.tex: main paper (14 pages, IEEEtran)
- SupplementaryMaterial.tex: extended supplementary (5 sections, S1-S5)
- figures/: all main (Fig 1-4) and supplementary figures (Fig S1-S5)
- notebooks/: coupling_decomposition.ipynb, Real-data_application.ipynb
- code/: fig2_power_curve.py, fig3_from_synthetic.py, supp_figures.py
- README.md with mathematical overview and reproduction instructions"

git branch -M main
git push -u origin main

echo ""
echo "=== Done! Repository pushed to GitHub ==="
echo "    $REPO_URL"
