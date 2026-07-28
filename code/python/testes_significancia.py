# -*- coding: utf-8 -*-
"""Testes de significância estatística sobre estatistica_out/estatistica_resultados.csv.

Desenho: 5 técnicas estocásticas x 30 sementes (números aleatórios comuns:
a mesma semente inicia todas as técnicas -> desenho PAREADO/bloqueado) para cada
K in {2,5,13}. Métrica primária: hipervolume (HV) como fração do HV da frente
exata (mesma convenção do benchmark).

Testes (por K):
  - Friedman (omnibus, blocos = sementes) -> há diferença entre técnicas?
  - Pós-teste par a par: Wilcoxon (postos sinalizados, pareado) + Holm-Bonferroni.
  - Nemenyi: diferença crítica (CD) sobre postos médios (para o diagrama de CD).
Saídas: resumo CSV, tabela LaTeX (tab_estatistica.tex) e figura boxplot.
"""
import csv, itertools
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
OUTDIR = BASE / "estatistica_out"
FIGDIR = Path(r"C:\Users\engmc\OneDrive - Fatecie\Inteligência de Informações e Dados\Inteligência\inbox\Dissertacao_Final\figuras")

TECHS = ["MO_MCSHADE", "MO_MCDE", "NSGA2", "MOEAD_AV", "MOPSO"]
LABEL = {"MO_MCSHADE": "MO-MC-SHADE", "MO_MCDE": "MO-MCDE", "NSGA2": "NSGA-II",
         "MOEAD_AV": "MOEA/D-AV", "MOPSO": "MOPSO"}
KS = [2, 5, 13]
NSEEDS = 30

# paleta não-azul (consistente com as figuras dos novos algoritmos)
CORES = {"MO_MCSHADE": "#E8730C", "MO_MCDE": "#0B7A6B", "NSGA2": "#4D4D4D",
         "MOEAD_AV": "#8C6D1F", "MOPSO": "#7A5C99"}

# ---------- frente exata (para normalizar o HV) ----------
def ler(nome):
    L = [r for r in csv.reader(open(BASE / nome, encoding="utf-8-sig")) if r]
    return np.array([[float(x) for x in r[1:]] for r in L[1:]])

def pf(P):
    if len(P) == 0: return P
    p = P[np.lexsort((P[:, 1], P[:, 0]))]; out = []; b = np.inf
    for q in p:
        if q[1] < b - 1e-12: out.append(q); b = q[1]
    return np.array(out)

MC, MG = ler("Custos.csv"), ler("GWP.csv")
front = np.array([[0., 0.]])
for i in range(MC.shape[0]):
    loc = pf(np.column_stack([MC[i], MG[i]]))
    front = pf((front[:, None, :] + loc[None, :, :]).reshape(-1, 2))
cmin, cmax = front[:, 0].min(), front[:, 0].max()
gmin, gmax = front[:, 1].min(), front[:, 1].max()
def norm(P):
    o = P.astype(float).copy(); o[:, 0] = (o[:, 0] - cmin) / (cmax - cmin); o[:, 1] = (o[:, 1] - gmin) / (gmax - gmin); return o
frontN = norm(front); REF = np.array([1.1, 1.1])
def hv(P):
    P = P[(P[:, 0] < REF[0]) & (P[:, 1] < REF[1])]
    if len(P) == 0: return 0.0
    P = pf(P); P = P[np.argsort(P[:, 0])]; h = 0.0; ly = REF[1]
    for x, y in P: h += (REF[0] - x) * (ly - y); ly = y
    return h
HVX = hv(frontN)

# ---------- HV por (K, técnica, semente) ----------
df = pd.read_csv(OUTDIR / "estatistica_resultados.csv")
HV = {}   # (K, tech) -> array de 30 HV
for K in KS:
    for t in TECHS:
        sub = df[(df["K"] == K) & (df["Algoritmo"] == t)]
        vals = np.full(NSEEDS, np.nan)
        for rep, g in sub.groupby("Rep"):
            P = norm(g[["Custo_Total", "Impacto_Ambiental"]].values.astype(float))
            vals[int(rep) - 1] = hv(P) / HVX
        HV[(K, t)] = vals

# ---------- estatísticas descritivas + testes ----------
resumo = []       # linhas do resumo
friedman = []     # por K
pairwise = []     # por K, por par
for K in KS:
    M = np.column_stack([HV[(K, t)] for t in TECHS])   # 30 x 5
    # postos por bloco (1 = melhor HV, i.e., maior). ranked ascending of -HV
    ranks = np.apply_along_axis(lambda r: stats.rankdata(-r), 1, M)  # 30 x 5
    meanrank = ranks.mean(axis=0)
    # Friedman
    chi2, p = stats.friedmanchisquare(*[M[:, j] for j in range(len(TECHS))])
    friedman.append({"K": K, "Friedman_chi2": chi2, "Friedman_p": p})
    for j, t in enumerate(TECHS):
        v = HV[(K, t)]
        resumo.append({"K": K, "Tecnica": LABEL[t], "HV_mediana": np.median(v),
                       "HV_IQR": np.percentile(v, 75) - np.percentile(v, 25),
                       "HV_media": v.mean(), "HV_std": v.std(ddof=1),
                       "HV_min": v.min(), "HV_max": v.max(), "posto_medio": meanrank[j]})
    # pós-teste par a par: Wilcoxon signed-rank + Holm-Bonferroni
    pares = list(itertools.combinations(range(len(TECHS)), 2))
    praw = []
    for (a, b) in pares:
        try:
            _, pw = stats.wilcoxon(M[:, a], M[:, b], zero_method="wilcox")
        except ValueError:
            pw = 1.0  # todas as diferenças nulas
        praw.append(pw)
    # Holm-Bonferroni
    order = np.argsort(praw); m = len(praw); padj = np.empty(m)
    run_max = 0.0
    for rank, idx in enumerate(order):
        val = (m - rank) * praw[idx]
        run_max = max(run_max, min(val, 1.0)); padj[idx] = run_max
    # Nemenyi CD (k=5, alfa=0.05): q=2.728
    q = {2: 1.960, 3: 2.343, 4: 2.569, 5: 2.728, 6: 2.850}[len(TECHS)]
    CD = q * np.sqrt(len(TECHS) * (len(TECHS) + 1) / (6.0 * NSEEDS))
    for k, (a, b) in enumerate(pares):
        pairwise.append({"K": K, "A": LABEL[TECHS[a]], "B": LABEL[TECHS[b]],
                         "Wilcoxon_p": praw[k], "p_Holm": padj[k],
                         "signif_0.05": padj[k] < 0.05,
                         "difPostos": abs(meanrank[a] - meanrank[b]),
                         "Nemenyi_CD": CD, "signif_Nemenyi": abs(meanrank[a] - meanrank[b]) > CD})

Rresumo = pd.DataFrame(resumo); Rfried = pd.DataFrame(friedman); Rpair = pd.DataFrame(pairwise)
Rresumo.to_csv(OUTDIR / "_estatistica_resumo.csv", index=False)
Rfried.to_csv(OUTDIR / "_estatistica_friedman.csv", index=False)
Rpair.to_csv(OUTDIR / "_estatistica_pairwise.csv", index=False)

pd.set_option("display.float_format", lambda x: f"{x:.4f}")
print("=== Frente exata:", len(front), "pts | HV exato", round(HVX, 4), "===\n")
print("RESUMO (mediana HV, IQR, posto medio) por K e tecnica:\n")
print(Rresumo.to_string(index=False))
print("\nFRIEDMAN (omnibus) por K:\n"); print(Rfried.to_string(index=False))
print("\nPARES significativos (Holm < 0.05):\n")
print(Rpair[Rpair["signif_0.05"]][["K", "A", "B", "Wilcoxon_p", "p_Holm"]].to_string(index=False))

# ---------- tabela LaTeX ----------
def fmt(v, iqr):
    return f"{v:.4f}\\,({iqr:.4f})"
lines = []
for t in TECHS:
    cells = []
    for K in KS:
        r = Rresumo[(Rresumo["K"] == K) & (Rresumo["Tecnica"] == LABEL[t])].iloc[0]
        best = Rresumo[Rresumo["K"] == K]["HV_mediana"].max()
        cell = fmt(r["HV_mediana"], r["HV_IQR"])
        if abs(r["HV_mediana"] - best) < 1e-9: cell = "\\textbf{" + cell + "}"
        cells.append(cell)
    lines.append(f"{LABEL[t]} & " + " & ".join(cells) + "\\\\")
fried_row = "Friedman $p$ & " + " & ".join(
    f"{Rfried[Rfried['K']==K]['Friedman_p'].iloc[0]:.2e}" for K in KS) + "\\\\"
tex = f"""\\begin{{table}}[htbp]
\\centering
\\caption{{Validação estatística (30 execuções independentes por técnica e por $K$): mediana do hipervolume (fração do exato) com o intervalo interquartil entre parênteses. Melhor mediana por $K$ em negrito. A linha final traz o $p$-valor do teste de Friedman (omnibus). Os pós-testes par a par (Wilcoxon pareado com correção de Holm--Bonferroni) constam no texto.}}
\\label{{tab:estatistica}}
\\small
\\begin{{tabular}}{{lrrr}}
\\hline
Técnica & $K=2$ & $K=5$ & $K=13$\\\\ \\hline
{chr(10).join(lines)}
\\hline
{fried_row}
\\hline
\\end{{tabular}}
\\end{{table}}
"""
Path(BASE / "tab_estatistica.tex").write_text(tex, encoding="utf-8", newline="\n")
print("\n[gerado tab_estatistica.tex]")

# ---------- figura boxplot ----------
fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), sharey=False)
for ax, K in zip(axes, KS):
    data = [HV[(K, t)] for t in TECHS]
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=True)
    for patch, t in zip(bp["boxes"], TECHS):
        patch.set_facecolor(CORES[t]); patch.set_alpha(0.75)
    for med in bp["medians"]: med.set_color("black")
    ax.set_xticks(range(1, len(TECHS) + 1))
    ax.set_xticklabels([LABEL[t] for t in TECHS], rotation=40, ha="right", fontsize=8)
    ax.set_title(f"$K = {K}$"); ax.grid(True, axis="y", alpha=0.25)
    if ax is axes[0]: ax.set_ylabel("Hipervolume (fração do exato)")
fig.suptitle("Distribuição do hipervolume em 30 execuções independentes", fontsize=11)
fig.tight_layout()
fig.savefig(FIGDIR / "fig11_estatistica_boxplot.png", dpi=200)
plt.close(fig)
print("[gerado fig11_estatistica_boxplot.png]")
