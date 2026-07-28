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
  fronts with dominance pruning — the *Pareto sum* operation): **1918 non-dominated
  solutions**, cost ∈ [635.39, 7155.25], GWP ∈ [2.5956, 4.7994].
- **Consolidation-constrained (NP-hard):** at most **K** distinct suppliers are
  allowed in the whole chain (a bi-objective *p*-median variant). Runs for
  K ∈ {2, 3, 5, 8, 10, 13} are provided.

A broad set of techniques is benchmarked against the exact front — metaheuristics
(MOPSO, NSGA-II, MOGA, MOEA/D and variants, SI-CDC), multi-criteria methods (AHP,
PROMETHEE, ELECTRE), exact/deterministic methods (MILP *p*-median, branch-and-bound,
Lagrangian relaxation, lexicographic), plus three recent techniques: multi-child
differential evolution in multi-objective form (**MO-MCDE**, **MO-MC-SHADE**) and
multi-objective Bayesian algorithm configuration (**MO-SMAC-PE**).

All techniques are the work of their original authors. MCDE and MC-SHADE (Storn) are
single-objective optimisers, run here in multi-objective form by replacing the scalar
parent–child selection with GDE3-style Pareto-dominance selection (Kukkonen &
Lampinen, 2005); MO-SMAC is multi-objective by design (Rook et al.) and is used here
as a direct solver rather than in its native algorithm-configuration role.

## Repository structure

```
.
├── data/          Input dataset (synthetic 13×500 and original 13×4) + metadata
│   └── DATA_DICTIONARY.md
├── code/
│   ├── matlab/    Optimisers (benchmark + statistical-validation harness)
│   ├── python/    Metrics, exact front, MO-SMAC run, significance tests
│   └── README.md  How to run each script
├── results/       Pareto fronts and quality-indicator tables
│   └── estatistica_out/   30-run statistical study (raw data + summaries)
├── figures/       Figures of the recent techniques and of the statistical study
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
   `code/matlab/Master_Optimizer_REO_V5_MC.m`. The constraint knob `K_MAX` selects
   the regime; sweep K to reproduce the restricted runs. Output:
   `Fronteira_Pareto_Global*.xlsx`.
2. **Statistical validation** — run `code/matlab/Estatistica_REO_V5.m`; it executes
   the five main stochastic techniques 30× (fixed seeds) for K ∈ {2, 5, 13},
   resumable via per-checkpoint `.mat` files, and writes
   `estatistica_out/estatistica_resultados.csv`.
3. **Metrics and exact front** — Python 3.12+ (`pip install -r
   code/python/requirements.txt`). Run `code/python/analise_v5.py` to recompute the
   exact front (by decomposition) and the HV/IGD/GD/purity indicators against it.
4. **Significance tests** — `code/python/testes_significancia.py` consumes the CSV
   from step 2 and produces the Friedman / paired-Wilcoxon (Holm–Bonferroni) results
   and the box-plot figure.
5. **MO-SMAC** — `code/python/mosmac_reo.py` (uses the official `smac` package,
   ParEGO variant).

See `code/README.md` for command-line details.

## Results in this release

**Single-run benchmark.** Base-regime quality indicators (hypervolume as a fraction
of the exact HV) are in `results/_metricas_v5.csv`; the restricted K = 5 indicators in
`results/_metricas_k5_run1.csv`. Pareto fronts are in `results/Fronteira_Pareto_*`.

**30-run statistical study** (`results/estatistica_out/`): five stochastic techniques
× 30 independent seeds × K ∈ {2, 5, 13}. Median hypervolume, as a fraction of the
exact front:

| Technique | K = 2 | K = 5 | K = 13 |
|---|---|---|---|
| **MO-MC-SHADE** | **0.9311** | **0.9836** | **0.9923** |
| NSGA-II | 0.9270 | 0.9756 | 0.9833 |
| MOPSO | 0.9299 | 0.9686 | 0.9818 |
| MO-MCDE | 0.9232 | 0.9645 | 0.9806 |
| MOEA/D-AV | 0.9219 | 0.9513 | 0.9639 |

Friedman rejects equality of the techniques at every K (p = 7.6e-20, 1.6e-17,
3.0e-13). Paired Wilcoxon with Holm–Bonferroni correction confirms that
**MO-MC-SHADE is significantly better than every other technique** at every K.
Notably, MOEA/D-AV — the leader of the single-run benchmark — shows the largest
variance and the worst mean rank under replication, i.e. its first place was an
artefact of one favourable run. Raw per-solution data are in
`estatistica_resultados.csv`; summaries in `_estatistica_resumo.csv`,
`_estatistica_friedman.csv` and `_estatistica_pairwise.csv`; the distribution plot is
`figures/fig11_estatistica_boxplot.png`.

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
