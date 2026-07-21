# -*- coding: utf-8 -*-
"""Análise de métricas V5 — 17 técnicas MATLAB + MO-SMAC.

Recomputa a frente exata por decomposição (soma de Minkowski das frentes
locais por material, com poda de dominância) — mesma convenção do pacote
(_analise_temp/_analise_metricas): normalização pela FAIXA da frente exata,
HV 2D com ponto de referência (1,1;1,1) reportado como FRAÇÃO do HV exato,
IGD/GD/pureza idênticos. Assim os novos números são comparáveis aos da V3.

Uso:
  python analise_v5.py --xlsx Fronteira_Pareto_Global.xlsx \
      --smac Fronteira_Pareto_MOSMAC_base.csv --out _metricas_v5.csv
"""
import argparse
import csv
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parent


# ----------------------------------------------------------------- dados
def ler_csv_matriz(nome):
    linhas = []
    with open(BASE / nome, encoding="utf-8-sig") as f:
        for row in csv.reader(f):
            if row:
                linhas.append(row)
    return np.array([[float(x) for x in r[1:]] for r in linhas[1:]])


def pareto_filter(points):
    """Não-dominados (minimização 2D)."""
    if len(points) == 0:
        return points
    pts = points[np.lexsort((points[:, 1], points[:, 0]))]
    out = []
    best = np.inf
    for p in pts:
        if p[1] < best - 1e-12:
            out.append(p)
            best = p[1]
    return np.array(out)


def frente_exata(MC, MG):
    """Decomposição: soma de Minkowski incremental das frentes locais."""
    nVar, nForn = MC.shape
    front = np.array([[0.0, 0.0]])
    for i in range(nVar):
        local = pareto_filter(np.column_stack([MC[i], MG[i]]))
        # Minkowski: front (x) local, depois poda
        soma = (front[:, None, :] + local[None, :, :]).reshape(-1, 2)
        front = pareto_filter(soma)
    return front


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", default="Fronteira_Pareto_Global.xlsx")
    ap.add_argument("--smac", default=None, help="CSV do MO-SMAC (opcional)")
    ap.add_argument("--out", default="_metricas_v5.csv")
    args = ap.parse_args()

    MC = ler_csv_matriz("Custos.csv")
    MG = ler_csv_matriz("GWP.csv")
    front = frente_exata(MC, MG)
    print(f"Frente exata: {len(front)} pontos | "
          f"custo [{front[:,0].min():.2f}; {front[:,0].max():.2f}] | "
          f"GWP [{front[:,1].min():.4f}; {front[:,1].max():.4f}]")

    cmin, cmax = front[:, 0].min(), front[:, 0].max()
    gmin, gmax = front[:, 1].min(), front[:, 1].max()

    def norm(P):
        out = P.astype(float).copy()
        out[:, 0] = (out[:, 0] - cmin) / (cmax - cmin)
        out[:, 1] = (out[:, 1] - gmin) / (gmax - gmin)
        return out

    frontN = norm(front)
    REF = np.array([1.1, 1.1])

    def hv2d(P, ref=REF):
        P = P[(P[:, 0] < ref[0]) & (P[:, 1] < ref[1])]
        if len(P) == 0:
            return 0.0
        P = pareto_filter(P)
        P = P[np.argsort(P[:, 0])]
        hv = 0.0
        last_y = ref[1]
        for x, y in P:
            hv += (ref[0] - x) * (last_y - y)
            last_y = y
        return hv

    def igd(aN):
        d = np.sqrt(((frontN[:, None, :] - aN[None, :, :]) ** 2).sum(-1))
        return d.min(axis=1).mean()

    def gd(aN):
        d = np.sqrt(((aN[:, None, :] - frontN[None, :, :]) ** 2).sum(-1))
        return d.min(axis=1).mean()

    hv_exact = hv2d(frontN)
    print(f"HV exato (normalizado, ref={REF}): {hv_exact:.4f}\n")
    front_set = set((round(c, 2), round(g, 4)) for c, g in front)

    # ---- carrega técnicas MATLAB ----
    df = pd.read_excel(BASE / args.xlsx)
    df["Alg"] = df["Algoritmo"].str.replace(r"\s*\(Run \d+\)", "", regex=True)
    frames = [df[["Alg", "Custo_Total", "Impacto_Ambiental", "Tempo_Total_Segundos"]]]

    # ---- MO-SMAC (CSV próprio) ----
    if args.smac and (BASE / args.smac).exists():
        s = pd.read_csv(BASE / args.smac)
        s = s.rename(columns={"Tempo_Total_Segundos": "Tempo_Total_Segundos"})
        s["Alg"] = s["Algoritmo"]
        frames.append(s[["Alg", "Custo_Total", "Impacto_Ambiental", "Tempo_Total_Segundos"]])

    allpts = pd.concat(frames, ignore_index=True)

    rows = []
    for alg, g in allpts.groupby("Alg"):
        P = g[["Custo_Total", "Impacto_Ambiental"]].values.astype(float)
        tempo = g["Tempo_Total_Segundos"].mean()
        Pu = np.unique(P, axis=0)
        PN = norm(Pu)
        nd = pareto_filter(Pu)
        on = sum(1 for c, gg in Pu if (round(c, 2), round(gg, 4)) in front_set)
        rows.append(dict(Alg=alg, N=len(P), Nuniq=len(Pu), Nnd=len(nd),
                         HVratio=hv2d(PN) / hv_exact, IGD=igd(PN), GD=gd(PN),
                         PurityPct=100 * on / len(Pu), Tempo=tempo))
    R = pd.DataFrame(rows).sort_values("HVratio", ascending=False)
    pd.set_option("display.width", 220)
    pd.set_option("display.float_format", lambda x: f"{x:.4f}")
    print("RANKING POR HIPERVOLUME (fração do HV exato):\n")
    print(R[["Alg", "N", "Nuniq", "Nnd", "HVratio", "IGD", "GD", "PurityPct", "Tempo"]].to_string(index=False))

    R.to_csv(BASE / args.out, index=False)
    print(f"\n[salvo {args.out}]  | frente exata: {len(front)} pts")


if __name__ == "__main__":
    main()
