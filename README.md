# Code for the revised manuscript (TCBB-2026-09-0796, major revision)

This folder reproduces every number, table and figure that is new or changed in the revision of
*A Sheaf–Hodge Framework for Higher-Order Metabolic Interactions in Microbial Communities*.
The original code of the first submission is kept unchanged in `python_code.zip`.

## Quick start (Jupyter)

```bash
pip install -r requirements.txt
```

1. Copy the SBML models into `data/` (see `data/README.md`).
2. Open the notebooks in `notebooks/` and run them top to bottom.

| Notebook | Manuscript | Runtime (2 cores) |
|---|---|---|
| `01_triad_glycerol_DHA.ipynb` | Section VII-C, Table IV, Fig. 6 | ~8 min |
| `02_community_coherence.ipynb` | Sections VII-A, VII-B, VII-D, Table III, Section V-B design | ~30 min |
| `03_figures.ipynb` | Figs. 1–5 and S1–S4 (schematics and synthetic validations) | < 1 min |

Every script can also be run on its own (`python triad_analysis.py`). Paths are set in one place,
`paths.py`, or through the environment variables `TCBB_DATA`, `TCBB_RESULTS` and `TCBB_FIGS`.
`results/` ships the JSON outputs used in the manuscript, so the figures and the noise sweep can be
regenerated without re-running the flux analyses.

## What changed with respect to the first submission

- **Exact test.** `nerve_lib.gls_test` uses the pseudoinverse of δ₀Σδ₀ᵀ and d = rank δ₀; for diagonal Σ
  the statistic reduces to a per-coordinate weighted dispersion, which gives the coordinate contributions
  used for localization.
- **Community FBA.** Minimal hypersaline medium (glycerol 10, CO₂ uptake 0, trace metals ≤ 0.1); the shared
  pool is mass-balanced even though the SBML declares it as boundary species; cooperative trade-off
  (α = 0.5, MICOM-style) followed by pFBA within 0.1 %; transfers by proportional allocation in the pool.
- **Model-derived Σ.** Standard deviations from FVA at 99 % of growth with total flux ≤ 1.1 × pFBA
  (floor 10⁻³), instead of an arbitrary coefficient of variation (`noise_sweep.py` shows why).
- **Positive control** on the community (`comm_pc.py`) and a curated **glycerol–DHA triad**
  (`triad_*.py`) with D. salina as a glycerol source.
- **Design.** Regularized D-optimal criterion log det(εI + 𝓘) − log det(εI), ε = 1, compared with the
  pseudo-determinant, the minimum-overlap heuristic and random designs (`design_rev.py`).
- **Figures.** Interface-only residual for Fig. S3, f(∅) = 0 for Fig. S4, panel labels outside the axes.

## Figure files

| Manuscript | Script | Output in `figures/` |
|---|---|---|
| Fig. 1 | `fig1.py` | `Fig1.pdf` |
| Fig. 2 | `gen_FIG_nerve_schematic.py` | `FIG_nerve_schematic_v4.png` (`Fig2.png` in the LaTeX source) |
| Fig. 3 | `gen_figs_rev.py` | `FIG4_regime_decomposition.png` |
| Fig. 4 | `gen_figs_rev.py` | `FIG3_factorization_residual.png` |
| Fig. 5 | `gen_figs_rev.py` | `FIG2_power_curve.png` |
| Fig. 6 | `fig_triad.py` | `FIG_triad.pdf` |
| Figs. S1–S4 | `gen_FIGS1-4_supplementary.py` | `FIGS1_…` to `FIGS4_…` |

The file names keep those of the first submission, where the figures were numbered differently.
All of them are reproduced pixel for pixel from the scripts.

## Files

| File | Role |
|---|---|
| `paths.py` | data, results and figure directories |
| `nerve_lib.py` | context nerve, δ₀, exact GLS test, H⁰ projection, regime decomposition |
| `parse_sbml.py` | SBML parser for the community model |
| `triad_supplements.py` | essential supplements of iSRU780 |
| `triad_lib.py` | joint cobra model for {D, S, H} contexts with a shared pool |
| `triad_analysis.py` | baseline test for the triad → `triad_base.json` |
| `triad_sens.py`, `triad_slope.py`, `triad_tradeoff.py`, `triad_pfbafactor.py` | sensitivity, transfer operator, objective and flux-bound robustness |
| `comm_lib.py` | HiGHS LP engine for community contexts (medium, trade-off, pFBA, FVA) |
| `comm_coords.py`, `comm_pipeline.py` | pairwise transfers, coordinates, coherence test |
| `comm_medium.py` | organism models and growth alone on the minimal medium (Supplementary Table) |
| `comm_stage1.py`, `comm_fva_par.py` | stage 1 (165 contexts), parallel FVA and the coherence test |
| `comm_pc.py` | positive control (lactate sink in Halobacterium) |
| `noise_sweep.py` | arbitrary-CV Σ against model-derived Σ |
| `design_rev.py` | experimental design on the community |
| `plotstyle.py`, `fig1.py`, `fig_triad.py`, `gen_*.py` | figures |

## Citation

If you use this code, please cite the TCBB article (under review) and
A. Centeno-Mejía, "Global metabolic coherence as a derived property of microbial communities:
A discrete Hodge framework," *J. Theor. Biol.*, vol. 625, p. 112438, 2026, doi:10.1016/j.jtbi.2026.112438.
