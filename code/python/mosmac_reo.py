# -*- coding: utf-8 -*-
"""
mosmac_reo.py — MO-SMAC (variante PE / ParEGO) aplicado ao problema bi-objetivo
de seleção de fornecedores REO da dissertação (13 materiais x 500 fornecedores).

Algoritmo: MO-SMAC — Rook, Benjamins, Bossek, Trautmann, Hoos & Lindauer,
"MO-SMAC: Multiobjective Sequential Model-Based Algorithm Configuration",
Evolutionary Computation 34(1):29-52, 2025 (doi:10.1162/evco_a_00371).

Implementação usada: pacote OFICIAL `smac` 2.4.0 (SMAC3), que incorpora os
componentes multiobjetivo do MO-SMAC publicados pelos autores:
  - agregação ParEGO (Tchebycheff aumentado, rho=0.05) => variante MO-SMAC-PE;
  - conjunto de INCUMBENTES não-dominados mantido pelo intensificador, com
    poda por crowding distance (utils/pareto_front.py) — o "MO-Intensify".
A variante PHVI exige o fork experimental dos autores (SMAC3 branch mosmac +
random_forest_run C++ compilado), que não builda em Windows; a escolha da PE
é registrada como limitação de reprodução.

Espaço de busca: 13 parâmetros CATEGÓRICOS (1 por material) com 500 opções
cada (o índice do fornecedor não tem ordem semântica -> categórico honesto).
Objetivos: custo total e GWP total (kg CO2-eq), ambos minimizados.
Restrição K opcional (consolidação): reparo Baldwiniano IDÊNTICO ao ApplyK
do Master_Optimizer (guloso: derruba o fornecedor menos usado e realoca pelos
Cn+Gn normalizados por material).

Saída: CSV com TODAS as avaliações + frente não-dominada + incumbentes,
no formato de colunas do harness MATLAB (Algoritmo, Custo_Total,
Impacto_Ambiental, Tempo_Total_Segundos, Num_Fornecedores_Distintos,
Selecao_Fornecedores) para consumo pelos scripts de métricas.

Uso:
  python mosmac_reo.py --n-trials 5000 --kmax inf --seed 42 --tag full
"""
import argparse
import csv
import time
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


# ----------------------------------------------------------------------
# Dados
# ----------------------------------------------------------------------
def carregar_dados():
    """Lê Custos.csv e GWP.csv (separador vírgula, decimal ponto, 1ª col = material)."""
    def ler(nome):
        linhas = []
        with open(BASE / nome, encoding="utf-8-sig") as f:
            for row in csv.reader(f):
                if row:
                    linhas.append(row)
        # remove cabeçalho; 1ª coluna = nome do material
        mats = [r[0] for r in linhas[1:]]
        vals = np.array([[float(x) for x in r[1:]] for r in linhas[1:]])
        return mats, vals

    mats, C = ler("Custos.csv")
    _, G = ler("GWP.csv")
    assert C.shape == G.shape, "Custos e GWP com dimensões diferentes"
    # normalização min-max POR MATERIAL (idêntica a Dados_Processo_REO da V4/V5)
    rc = C.max(axis=1) - C.min(axis=1)
    rc[rc == 0] = 1.0
    rg = G.max(axis=1) - G.min(axis=1)
    rg[rg == 0] = 1.0
    Cn = (C - C.min(axis=1, keepdims=True)) / rc[:, None]
    Gn = (G - G.min(axis=1, keepdims=True)) / rg[:, None]
    return mats, C, G, Cn, Gn


MATS, MC, MG, CN, GN = carregar_dados()
N_VAR, N_FORN = MC.shape  # 13 x 500


def apply_k(idx, K):
    """Porta fiel do ApplyK do Master_Optimizer (V4/V5, MATLAB).

    idx: array 1-based (como no MATLAB) com o fornecedor de cada material.
    Enquanto houver mais de K fornecedores distintos: elimina o MENOS usado
    (primeiro em caso de empate, como o min do MATLAB) e realoca cada material
    dele para o fornecedor remanescente de menor Cn+Gn daquele material.
    """
    idx = idx.copy()
    u = np.unique(idx)  # ordenado, como no MATLAB
    while u.size > K:
        counts = np.array([(idx == s).sum() for s in u])
        drop = u[np.argmin(counts)]          # 1º mínimo, igual ao MATLAB
        keep = u[u != drop]
        for m in np.where(idx == drop)[0]:
            sc = CN[m, keep - 1] + GN[m, keep - 1]
            idx[m] = keep[np.argmin(sc)]     # 1º mínimo, igual ao MATLAB
        u = np.unique(idx)
    return idx


def avaliar(idx_1based, K):
    """Custo e GWP totais de uma seleção (com reparo K se ativo)."""
    idx = np.asarray(idx_1based, dtype=int)
    if K < N_VAR:
        idx = apply_k(idx, K)
    cols = idx - 1
    linhas = np.arange(N_VAR)
    return float(MC[linhas, cols].sum()), float(MG[linhas, cols].sum()), idx


# ----------------------------------------------------------------------
# MO-SMAC
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-trials", type=int, default=5000)
    ap.add_argument("--kmax", default="inf")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tag", default="run")
    ap.add_argument("--max-incumbents", type=int, default=10,
                    help="tamanho máx. do conjunto de incumbentes (paper: 10)")
    args = ap.parse_args()

    K = float("inf") if str(args.kmax).lower() in ("inf", "") else int(args.kmax)
    K_efetivo = K if K < N_VAR else N_VAR

    from ConfigSpace import Categorical, ConfigurationSpace
    from smac import HyperparameterOptimizationFacade as Facade
    from smac import Scenario
    from smac.multi_objective.parego import ParEGO

    cs = ConfigurationSpace(seed=args.seed)
    for i in range(N_VAR):
        # índices de fornecedor como CATEGÓRICOS: 1..500 sem ordem semântica
        cs.add(Categorical(f"forn_{i:02d}", items=list(range(1, N_FORN + 1))))

    def target(config, seed: int = 0):
        x = np.array([config[f"forn_{i:02d}"] for i in range(N_VAR)], dtype=int)
        custo, gwp, _ = avaliar(x, K)
        return {"custo": custo, "gwp": gwp}

    outdir = BASE / f"smac3_out_{args.tag}"
    scenario = Scenario(
        cs,
        name=f"mosmac_reo_{args.tag}",
        objectives=["custo", "gwp"],
        deterministic=True,
        n_trials=args.n_trials,
        seed=args.seed,
        output_directory=outdir,
    )

    mo = ParEGO(scenario)  # rho=0.05 default, como no paper (Knowles 2006)
    intensifier = Facade.get_intensifier(scenario, max_incumbents=args.max_incumbents)

    smac = Facade(
        scenario,
        target,
        multi_objective_algorithm=mo,
        intensifier=intensifier,
        overwrite=True,
        logging_level=30,
    )

    t0 = time.perf_counter()
    incumbents = smac.optimize()
    tempo = time.perf_counter() - t0

    if not isinstance(incumbents, list):
        incumbents = [incumbents]

    # ------------------------------------------------------------------
    # Extrai TODAS as avaliações do runhistory e monta a frente global
    # ------------------------------------------------------------------
    rh = smac.runhistory
    evals = []  # (custo, gwp, x_reparado)
    for tk, tv in rh.items():
        if tv.cost is None:
            continue
        cfg = rh.get_config(tk.config_id)
        x = np.array([cfg[f"forn_{i:02d}"] for i in range(N_VAR)], dtype=int)
        custo, gwp, xrep = avaliar(x, K)   # reavalia p/ obter seleção reparada
        evals.append((custo, gwp, xrep))

    pts = np.array([[e[0], e[1]] for e in evals])
    nd = np.ones(len(evals), dtype=bool)
    order = np.lexsort((pts[:, 1], pts[:, 0]))  # varredura por custo p/ frente 2D
    best_g = np.inf
    for j in order:
        if pts[j, 1] < best_g - 1e-12:
            best_g = pts[j, 1]
        else:
            nd[j] = False
    # correção para empates exatos de custo: manter apenas o de menor gwp
    # (a varredura acima já garante isso, pois ordena por (custo, gwp))

    inc_set = set()
    for cfg in incumbents:
        x = tuple(int(cfg[f"forn_{i:02d}"]) for i in range(N_VAR))
        inc_set.add(x)

    out_csv = BASE / f"Fronteira_Pareto_MOSMAC_{args.tag}.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["Algoritmo", "Custo_Total", "Impacto_Ambiental",
                       "Tempo_Total_Segundos", "Num_Fornecedores_Distintos",
                       "Selecao_Fornecedores", "Eh_Incumbente"])
        for j in np.where(nd)[0]:
            custo, gwp, xrep = evals[j]
            sel = "[" + " ".join(str(int(v)) for v in xrep) + "]"
            ehinc = 1 if tuple(int(v) for v in xrep) in inc_set else 0
            wcsv.writerow(["MO-SMAC-PE", f"{custo:.6f}", f"{gwp:.6f}",
                           f"{tempo:.2f}", int(np.unique(xrep).size), sel, ehinc])

    ndist_max = max(int(np.unique(e[2]).size) for e in evals)
    print(f"\n===== MO-SMAC-PE concluido =====")
    print(f"K = {K_efetivo} | trials = {len(evals)} | tempo = {tempo:.1f}s")
    print(f"frente nao-dominada: {int(nd.sum())} pontos | incumbentes: {len(incumbents)}")
    print(f"max fornecedores distintos observado: {ndist_max} (deve ser <= {K_efetivo})")
    print(f"saida: {out_csv}")


if __name__ == "__main__":
    main()
