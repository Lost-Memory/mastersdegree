# Code

## MATLAB (`matlab/`)

Requires MATLAB with the Parallel Computing, Optimization and Statistics toolboxes.
Run from a folder that contains `Custos.csv` and `GWP.csv` (copy them from `../data/`).

### `Master_Optimizer_REO_V5_MC.m` — full benchmark
Runs the comparative benchmark of all techniques in both regimes.
- Constraint knob **`K_MAX`** (near line 55): `inf` (or ≥13) → unconstrained
  (separable) regime; an integer 2..13 → consolidation to at most K distinct
  suppliers. Sweep `K_MAX` ∈ {2,3,5,8,10,13} to reproduce the restricted runs.
- Headless mode (for `matlab -batch`): environment variables `HEADLESS=1`,
  `KMAX=<n|inf>`, `TECHS=<comma-list>` skip the GUI.
- Output: `Fronteira_Pareto_Global.xlsx` (unconstrained) or
  `Fronteira_Pareto_Global_K<K>.xlsx` (constrained), with the columns
  `Algoritmo, Custo_Total, Impacto_Ambiental, Tempo_Total_Segundos,
  Num_Fornecedores_Distintos, Selecao_Fornecedores`.

### `Estatistica_REO_V5.m` — statistical-validation harness
Runs the main stochastic techniques (MO-MC-SHADE, MO-MCDE, NSGA-II, MOEA/D-AV,
MOPSO) **30 times** (fixed seeds) for each K ∈ {2, 5, 13}, saving every run's final
Pareto front. Just run:
```matlab
>> Estatistica_REO_V5
```
- **Resumable:** each (K, technique) pair is a `.mat` checkpoint under
  `estatistica_out/`; re-running skips finished pairs.
- **Fast & stable:** one technique at a time, its 30 seeds in `parfor` on a bounded
  pool; uses vectorised dominance/non-dominated-sort (self-tested to be identical to
  the O(n²) originals).
- Output: `estatistica_out/estatistica_resultados.csv` (one row per solution:
  `K, Algoritmo, Rep, Custo_Total, Impacto_Ambiental, Tempo_Total_Segundos`).
- Test overrides (do **not** use for the final run): env `ESTAT_NSEEDS`,
  `ESTAT_KLIST`, `ESTAT_TECHS`, `ESTAT_MAXIT`.

## Python (`python/`)

```bash
pip install -r requirements.txt
```

### `analise_v5.py` — exact front + quality indicators
Recomputes the exact Pareto front by decomposition (Minkowski sum of per-input local
fronts with dominance pruning) and computes HV / IGD / GD / purity against it, for a
given benchmark `.xlsx` (and optionally a MO-SMAC CSV). Run from a folder with
`Custos.csv`, `GWP.csv` and the results file:
```bash
python analise_v5.py --xlsx Fronteira_Pareto_Global.xlsx \
       --smac Fronteira_Pareto_MOSMAC_base.csv --out _metricas_v5.csv
```

### `analise_k5_run1.py` — restricted-regime metrics (single run)
Computes the K=5 indicators for the recent techniques using one run each (fair
comparison with the published restricted table).

### `testes_significancia.py` — significance tests over the 30 runs
Consumes `estatistica_out/estatistica_resultados.csv` (produced by
`Estatistica_REO_V5.m`) and computes, per K: the Friedman omnibus test, the paired
Wilcoxon signed-rank post-hoc with Holm–Bonferroni correction, and the Nemenyi
critical-difference ranks. The design is paired because all techniques share the
same seeds (common random numbers).
```bash
python testes_significancia.py
```
- Output: `_estatistica_resumo.csv` (median and IQR per technique and K),
  `_estatistica_friedman.csv`, `_estatistica_pairwise.csv` and the box-plot
  `fig11_estatistica_boxplot.png`.

### `mosmac_reo.py` — MO-SMAC (ParEGO variant)
Runs the Bayesian algorithm-configuration method on the problem (13 categorical
parameters × 500 levels). Requires the official `smac` package (see requirements).
```bash
python mosmac_reo.py --n-trials 500 --kmax inf --seed 42 --tag base
```

## Metrics convention

All indicators are normalised by the range of the **exact** front; HV is reported as
a fraction of the exact HV; reference point (1.1, 1.1) in normalised space. This is
consistent across the MATLAB benchmark and the Python analysis (the deterministic
MILP *p*-median value cross-checks the convention).
