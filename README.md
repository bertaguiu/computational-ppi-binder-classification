# computational-ppi-binder-classification
Statistical analysis of computational metrics for binder/non-binder classification in de novo protein-protein design

````
# Prediction of Binding Classification in Computational Protein–Protein Design

**Treball de Final de Grau — Grau en Biotecnologia**
Barcelona Supercomputing Center (BSC-CNS)

---

## Overview

Code and data from a TFG on *de novo* protein-protein binding prediction. The question driving the work: do AlphaFold3 confidence metrics and Rosetta energy terms carry enough signal to separate binders from non-binders, before running a single experiment?

Two targets were studied:

- **BBF-14**: a 122-residue *de novo* protein from the BindCraft paper (Pacesa et al., 2025). 14 complexes (6 binders, 8 non-binders). The dataset is small but the labels are solid — all experimentally validated.
- **IL-7Rα** (*Interleukin-7 Receptor alpha chain*, ~460 aa): a receptor involved in T and B cell development. 90 complexes from the RFdiffusion IL-7Rα screen, with Kd values from Adaptyv Bio's independent replication (28 binders, 62 non-binders).

Complexes were predicted with **AlphaFold3** (iPTM, pTM, ranking score, cross-chain PAE) and relaxed with **Rosetta** (REU energy terms: total score, van der Waals, electrostatics, hydrogen bonds, etc.). Interface metrics (ΔG, ΔSASA, packstat, shape complementarity, interface hydrogen bonds) came from the `proteinModels` module of `prepare_proteins`. All jobs ran on **MareNostrum 5** at BSC — AlphaFold3 on ACC (GPU), Rosetta on GPP (CPU), via SLURM.

Analysis structure:

- **Bloc 0** — Dataset overview and class balance
- **Bloc 1** — Univariate discrimination: Mann-Whitney U test, Cliff's delta, per-metric AUC, FDR correction (Benjamini-Hochberg), forest plot and volcano plot
- **Bloc 2** — Redundancy and collinearity: Pearson correlation matrix (hierarchical clustering), VIF, down-selection to 14 non-redundant metrics
- **Bloc 3** — PCA, biplot, scree plot
- **Bloc 4** — L1-regularised logistic regression (grid search over C, stratified 5×10 CV or LOOCV for BBF-14), Random Forest feature importances, permutation test (n = 250)
- **Bloc 5** — Robustness: outlier detection (IQR + robust z-score), sensitivity of Mann-Whitney results and model performance to outlier removal, bootstrap CIs for top-6 metric AUCs, sensitivity to Kd threshold (IL-7Rα only)

---

## Repository Structure

```
.
├── bbf14/
│   ├── bindcraft_bbf14.ipynb              # BindCraft dataset exploration
│   ├── bindcraft_dataset.ipynb            # Dataset preparation
│   ├── create_final_csv_bbf14.ipynb       # Merge AF3 + Rosetta outputs into final CSV
│   ├── final_analysis_bbf14.ipynb         # Main analysis
│   ├── post_af3_bbf14.ipynb               # AlphaFold3 post-processing
│   ├── post_rosetta_bbf14.ipynb           # Rosetta post-processing
│   ├── prepare_jobs_af3_bbf14.ipynb       # SLURM job setup for AlphaFold3
│   ├── prepare_jobs_rosetta_bbf14.ipynb   # SLURM job setup for Rosetta
│   ├── bindcraft_dataset/                 # Raw BindCraft dataset
│   ├── final_bbf14/                       # Final merged CSV
│   ├── bindcraft_dataset_outputs/         # Plots — dataset exploration
│   ├── final_analysis_bbf14_outputs/      # Plots and CSVs — main analysis
│   └── bindcraft_bbf14_outputs/           # Plots — additional BBF-14 outputs
│
└── il7ra/
    ├── create_final_csv_il7ra.ipynb       # Merge AF3 + Rosetta outputs into final CSV
    ├── final_analysis_il7ra.ipynb         # Main analysis
    ├── final_analysis_il7ra_defB.ipynb    # Stricter binder definition (see note below)
    ├── il7ra_binders_preparation.ipynb    # Dataset preparation
    ├── post_af3_il7ra.ipynb               # AlphaFold3 post-processing
    ├── post_rosetta_il7ra.ipynb           # Rosetta post-processing
    ├── prepare_jobs_af3_il7ra.ipynb       # SLURM job setup for AlphaFold3
    ├── prepare_jobs_rosetta_il7ra.ipynb   # SLURM job setup for Rosetta
    ├── final_il7ra.csv                    # Final merged CSV
    ├── il7ra_binders_dataset.csv          # Experimental labels and Kd values
    └── final_analysis_il7ra_outputs/      # Plots and CSVs — main analysis
```

---

## Computational Pipeline

Same pipeline for both targets:

```
FASTA input (target:binder sequences)
        │
        ▼
AlphaFold3 (ACC partition, GPU)
  → structural confidence metrics per complex
        │
        ▼
Rosetta Relax (GPP partition, CPU, nstruct=100)
  → folding energy terms per structure
        │
        ▼
Rosetta Interface Analysis (proteinModels module)
  → binding interface metrics
        │
        ▼
create_final_csv_*.ipynb
  → single merged CSV, values averaged across seeds/structures
        │
        ▼
final_analysis_*.ipynb
  → statistical analysis (Blocs 0–5)
```

Only the top-ranked AF3 model (by `ranking_score`) is used as input to Rosetta.

---

## Key Files

| File | Description |
|------|-------------|
| `bbf14/final_bbf14/` | Merged dataset for BBF-14 (14 complexes, 48 metrics) |
| `il7ra/final_il7ra.csv` | Merged dataset for IL-7Rα (90 complexes, 60 columns) |
| `il7ra/il7ra_binders_dataset.csv` | Experimental binding labels and Kd values for IL-7Rα |
| `*/final_analysis_*_outputs/` | All output figures and CSVs |

---

## Dependencies

```
numpy, pandas, matplotlib, seaborn, scipy, statsmodels, scikit-learn
```

Three internal BSC libraries handle job preparation and post-processing but are not included in this repository:

- `prepare_proteins` — sequence prep, AF3 and Rosetta job setup, interface analysis (`proteinModels`)
- `bsc_calculations` — SLURM specifications for MareNostrum 5
- `bioprospecting` — auxiliary utilities

---

*Berta Guiu Gorgas — Universitat de Vic / Barcelona Supercomputing Center, 2026*
````
