# Multi-objective supplier selection for rare-earth-oxide production — dataset and code

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21481074-blue.svg)](https://doi.org/10.5281/zenodo.21481074)
[![License: MIT](https://img.shields.io/badge/Code%20License-MIT-yellow.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA.md)

Open dataset and reproducible code for the bi-objective (cost vs. Global Warming
Potential) supplier-selection problem studied in the master's dissertation of
**Mateus Crepaldi da Silva** (Graduate Program in Computer Science — PPGCC,
State University of Maringá — UEM, Brazil), supervised by
**Prof. Dr. Rodrigo Clemente Thom de Souza**.

> **Título (PT):** *Otimização Multiobjetivo da Cadeia de Suprimentos de Terras
> Raras: Uma Abordagem Baseada em Inteligência de Enxame para o Balanço entre
> Custos e Descarbonização.*

---

## Overview

The problem selects one supplier for each of **13 intermediate inputs** used in the
production of rare-earth oxides (REO), minimising simultaneously the **total cost**
and the **total Global Warming Potential** (GWP, kg CO₂-eq). Two regimes are studied:

- **Unconstrained (separable):** each input's supplier is chosen independently. The
  exact Pareto front is obtained by decomposition (Minkowski sum of per-input local
  fronts with dominance pruning): **1918 non-dominated solutions**, cost ∈ [635.39,
  7155.25], GWP ∈ [2.5956, 4.7994].
- **Consolidation-constrained (NP-hard):** at most **K** distinct suppliers are
  allowed in the whole chain (a bi-objective *p*-median variant). Runs for
  K ∈ {2, 3, 5, 8, 10, 13} are provided.

A broad set of techniques is benchmarked against the exact front — metaheuristics
(MOPSO, NSGA-II, MOGA, MOEA/D and variants, SI-CDC), multi-criteria methods (AHP,
PROMETHEE, ELECTRE), exact/deterministic methods (MILP *p*-median, branch-and-bound,
Lagrangian relaxation, lexicographic), plus three recent techniques adapted here to
the multi-objective setting: **MO-MCDE** and **MO-MC-SHADE** (multi-child
differential evolution) and **MO-SMAC-PE** (Bayesian algorithm configuration).

## Repository structure

```
.
├── data/          Input dataset (synthetic 13×500 and original 13×4) + metadata
│   └── DATA_DICTIONARY.md
├── code/
│   ├── matlab/    Optimisers (benchmark + statistical-validation harness)
│   ├── python/    Metrics, exact front, and the MO-SMAC run
│   └── README.md  How to run each script
├── results/       Pareto fronts and quality-indicator tables produced by the code
├── figures/       Selected figures of the proposed algorithms
├── LICENSE        MIT (code)
├── LICENSE-DATA.md  CC BY 4.0 (data)
└── CITATION.cff   Citation metadata
```

## The dataset

- `data/Custos.csv`, `data/GWP.csv` — the **expanded synthetic** matrices: 13 rows
  (inputs) × 500 columns (candidate suppliers). Separator `,`; decimal `.`; no quotes.
- `data/Custos_Original.csv`, `data/GWP_Original.csv` — the **original 13×4** matrix
  (four supply regions) that anchors the expansion.
- `data/Metadados_Geracao.xlsx` — generation configuration (random seeds, correlation
  target, region layout, preserved anchors, injected outliers).

The 500 suppliers are organised in four regions of 125 each — `CN-NM`, `RER`, `RoW`
and `CN-NM ADAPTADO BR` — with the region encoded in each column name. The first
supplier of each region (columns 1, 126, 251, 376) reproduces exactly the
corresponding value of the original 13×4 matrix. The expansion is statistically
controlled (log-normal in log-space, adaptive σ, cost–GWP correlation ρ = −0.75 by
bivariate conditional sampling, intentional outliers). Reproducibility seeds:
**cost = 42, GWP = 24**. See `data/DATA_DICTIONARY.md` for full column semantics.

**Provenance note.** The original 13×4 values derive from a Life Cycle Assessment
inventory (EcoInvent database via SimaPro, IPCC GWP-100 method). Only the aggregated
cost/GWP figures are distributed here; the raw EcoInvent database is **not** included
and requires a separate license.

## Reproducing the results

1. **Benchmark (both regimes)** — MATLAB (Parallel Computing, Optimization and
   Statistics toolboxes). Place `Custos.csv` and `GWP.csv` next to the script and run
   `code/matlab/Master_Optimizer_REO_V5_MC.m`. The constraint knob `K_MAX` (line ~55)
   selects the regime; sweep K to reproduce the restricted runs. Output:
   `Fronteira_Pareto_Global*.xlsx`.
2. **Statistical validation** — run `code/matlab/Estatistica_REO_V5.m`; it executes
   the main stochastic techniques 30× per seed for K ∈ {2, 5, 13}, resumable via
   per-checkpoint `.mat` files, and writes `estatistica_out/estatistica_resultados.csv`.
3. **Metrics and exact front** — Python 3.12+ (`pip install -r
   code/python/requirements.txt`). Run `code/python/analise_v5.py` to recompute the
   exact front (by decomposition) and the HV/IGD/GD/purity indicators against it.
4. **MO-SMAC** — `code/python/mosmac_reo.py` (uses the official `smac` package,
   ParEGO variant).

See `code/README.md` for command-line details.

## Results (this release)

Base-regime quality indicators (hypervolume as a fraction of the exact HV) are in
`results/_metricas_v5.csv`; the restricted K=5 indicators in
`results/_metricas_k5_run1.csv`. The single-run benchmark tables of the dissertation
are reproduced here; the multi-seed statistical-significance study (Wilcoxon /
Friedman–Nemenyi with Holm–Bonferroni correction) is produced by the harness in
step 2.

## How to cite

If you use this dataset or code, please cite this record via its DOI:
**[10.5281/zenodo.21481074](https://doi.org/10.5281/zenodo.21481074)**. A
machine-readable citation is in [`CITATION.cff`](CITATION.cff).

## License

- **Code** — [MIT](LICENSE).
- **Data** — [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE-DATA.md).

## Acknowledgments

Third-party algorithms referenced but **not redistributed** here: SHADE (Tanabe &
Fukunaga, 2013), JADE (Zhang & Sanderson, 2009), GDE3 (Kukkonen & Lampinen, 2005),
MCDE/MC-SHADE (Storn & Price, 2025/2026), MO-SMAC (Rook, Bossek et al., 2026) and
ParEGO (Knowles, 2006). The MO-SMAC run relies on the authors' official `smac`
package.

---

*Autor: Mateus Crepaldi da Silva · Orientador: Prof. Dr. Rodrigo Clemente Thom de
Souza · PPGCC/UEM, 2026.*
