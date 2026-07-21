# -*- coding: utf-8 -*-
"""Métricas K=5 dos novos algoritmos usando APENAS a Run 1 (execução única),
para comparabilidade justa com tab_restrito (single-run) e com o regime base.
HV como fração do HV exato irrestrito (mesma convenção — MILP cross-check ok).
"""
import numpy as np, pandas as pd, csv
from pathlib import Path
BASE = Path(__file__).resolve().parent

def ler(nome):
    L = [r for r in csv.reader(open(BASE/nome, encoding="utf-8-sig")) if r]
    return np.array([[float(x) for x in r[1:]] for r in L[1:]])

def pf(P):
    if len(P)==0: return P
    p = P[np.lexsort((P[:,1],P[:,0]))]; out=[]; b=np.inf
    for q in p:
        if q[1] < b-1e-12: out.append(q); b=q[1]
    return np.array(out)

MC, MG = ler("Custos.csv"), ler("GWP.csv")
front = np.array([[0.,0.]])
for i in range(MC.shape[0]):
    loc = pf(np.column_stack([MC[i],MG[i]]))
    front = pf((front[:,None,:]+loc[None,:,:]).reshape(-1,2))
cmin,cmax = front[:,0].min(),front[:,0].max(); gmin,gmax = front[:,1].min(),front[:,1].max()
def norm(P):
    o=P.astype(float).copy(); o[:,0]=(o[:,0]-cmin)/(cmax-cmin); o[:,1]=(o[:,1]-gmin)/(gmax-gmin); return o
frontN=norm(front); REF=np.array([1.1,1.1])
def hv(P):
    P=P[(P[:,0]<REF[0])&(P[:,1]<REF[1])]
    if len(P)==0: return 0.
    P=pf(P); P=P[np.argsort(P[:,0])]; h=0.; ly=REF[1]
    for x,y in P: h+=(REF[0]-x)*(ly-y); ly=y
    return h
def igd(aN):
    d=np.sqrt(((frontN[:,None,:]-aN[None,:,:])**2).sum(-1)); return d.min(1).mean()
hvx=hv(frontN)

df = pd.read_excel(BASE/"Fronteira_Pareto_Global_K5.xlsx")
# só Run 1 (MILP é determinístico; Run 1 basta para todos)
df = df[df["Algoritmo"].str.contains(r"\(Run 1\)")].copy()
df["Alg"]=df["Algoritmo"].str.replace(r"\s*\(Run \d+\)","",regex=True)

# MO-SMAC k5 (execução única por construção)
s = pd.read_csv(BASE/"Fronteira_Pareto_MOSMAC_k5.csv"); s["Alg"]=s["Algoritmo"]
allp = pd.concat([df[["Alg","Custo_Total","Impacto_Ambiental","Tempo_Total_Segundos"]],
                  s[["Alg","Custo_Total","Impacto_Ambiental","Tempo_Total_Segundos"]]], ignore_index=True)

rows=[]
for alg,g in allp.groupby("Alg"):
    P=np.unique(g[["Custo_Total","Impacto_Ambiental"]].values.astype(float),axis=0)
    rows.append(dict(Alg=alg, Nnd=len(pf(P)), HVratio=hv(norm(P))/hvx, IGD=igd(norm(P)),
                     Tempo=g["Tempo_Total_Segundos"].mean()))
R=pd.DataFrame(rows).sort_values("HVratio",ascending=False)
pd.set_option("display.float_format",lambda x:f"{x:.4f}")
print(f"Frente exata: {len(front)} pts | HV exato {hvx:.4f}\n")
print("K=5, EXECUCAO UNICA (Run 1):\n")
print(R.to_string(index=False))
R.to_csv(BASE/"_metricas_k5_run1.csv",index=False)
