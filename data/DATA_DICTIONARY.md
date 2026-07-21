# Data dictionary

All matrices are plain CSV: separator `,`, decimal point `.`, UTF-8, no quoting.
Do **not** open-and-save these files in a locale that uses `;`/`,` (e.g. Excel in
pt-BR), as it will corrupt the separators.

## `Custos.csv` and `GWP.csv` — expanded dataset (13 × 500)

| Element | Description |
|---|---|
| Rows | 1 header + **13 inputs** (materials/processes), one row each. |
| Column 1 | `Material` — input name (see list below). |
| Columns 2–501 | **500 candidate suppliers**, named `Fornecedor NNN (REGION)`. |
| Values in `Custos.csv` | Unit **cost** of the input from that supplier (relative monetary units). |
| Values in `GWP.csv` | **Global Warming Potential** of the input from that supplier (kg CO₂-eq), IPCC GWP-100. |

**Suppliers / regions.** The 500 suppliers are four regions of 125 each, encoded in
the column name:

| Columns | Region tag | Meaning |
|---|---|---|
| 2–126   | `CN-NM` | China (baseline) |
| 127–251 | `RER` | Europe |
| 252–376 | `RoW` | Rest of world |
| 377–501 | `CN-NM ADAPTADO BR` | China profile adapted to a Brazilian electricity mix |

**Anchors.** The first supplier of each region (columns 2, 127, 252, 377 → suppliers
001, 126, 251, 376) reproduces exactly the corresponding column of the original 13×4
matrix. These are the real LCA anchor points; the remaining 496 suppliers are
synthetic (statistically-controlled expansion).

**The 13 inputs.** Blasting · Fatty acid · Fatty alcohol · Hard coal · Mine
infrastructure – iron · Rare earth tailings · Recultivation – iron mine · Sodium
hydroxide (50% solution) · Sodium silicate (spray powder 80%) · Steel – chromium
steel 18/8 · Diesel · Electricity · Wastewater.

Two inputs are **degenerate by design**, fixing the extremes of the Pareto front:
*Rare earth tailings* has GWP = 0 for all suppliers, and *Wastewater* has cost = 0.
These are intentional and must not be "cleaned".

## `Custos_Original.csv` and `GWP_Original.csv` — original dataset (13 × 4)

Same layout, but only the four real supply regions (columns `CN-NM`, `RER`, `RoW`,
`CN-NM ADAPTADO BR`). This is the baseline the expansion is anchored on.

## `Metadados_Geracao.xlsx` — generation metadata

Sheets document the expansion: global configuration (generation timestamp; random
seeds **cost = 42**, **GWP = 24**; target correlation **ρ = −0.75** within region;
125 suppliers/region; anchor positions 1/126/251/376), per-material summary,
preserved anchors, and the intentionally injected cost and GWP outliers.

## Generation method (summary)

The 496 non-anchor suppliers per material are drawn log-normally in log-space with an
adaptive standard deviation (σ ∈ [0.08, 0.30]); cost and GWP are correlated
(ρ = −0.75) by bivariate conditional sampling within each region; ten outliers per
material are injected at ~3σ on purpose. The full method is described in the
dissertation (Chapter 3).

## Provenance / licensing of the underlying LCA figures

The original 13×4 values derive from a Life Cycle Assessment inventory (EcoInvent
database via SimaPro, IPCC GWP-100). Only the **aggregated** cost/GWP figures are
distributed here under CC BY 4.0; the **raw EcoInvent database is not included** and
requires its own license. Users intending to reproduce or extend the LCA inventory
must obtain EcoInvent/SimaPro access independently.
