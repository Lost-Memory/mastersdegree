# Disponibilidade de Dados e Código — texto para a dissertação

Cole um dos blocos abaixo na dissertação (fim da metodologia ou seção própria),
**substituindo `<DOI>` e `<URL>`** pelos valores gerados ao arquivar no Zenodo.

## Versão em prosa (PT)

> **Disponibilidade de dados e código.** Em conformidade com os princípios de
> ciência aberta, o conjunto de dados (13 insumos × 500 fornecedores, incluindo a
> base original de 13 × 4 e os metadados de geração), o código-fonte (otimizadores
> em MATLAB e rotinas de análise em Python) e os resultados (frentes de Pareto e
> indicadores de qualidade) estão disponíveis publicamente em
> \url{https://github.com/Lost-Memory/mastersdegree}, sob licença
> Creative Commons Attribution 4.0 (dados) e MIT (código), e arquivados de forma
> permanente sob o identificador DOI \texttt{10.5281/zenodo.21481074}. Os valores de custo e GWP das quatro
> regiões-âncora derivam de inventário de Avaliação de Ciclo de Vida (base
> EcoInvent, \textit{software} SimaPro, método IPCC GWP-100); apenas os valores
> agregados são disponibilizados, não a base EcoInvent, cujo uso requer licença
> própria.

## Versão LaTeX (para `\input`)

```latex
\section*{Disponibilidade de Dados e Código}
\addcontentsline{toc}{section}{Disponibilidade de Dados e Código}
Em conformidade com os princípios de ciência aberta, o conjunto de dados
(13~insumos $\times$ 500~fornecedores, incluindo a base original de $13 \times 4$ e
os metadados de geração), o código-fonte (otimizadores em MATLAB e rotinas de
análise em Python) e os resultados (frentes de Pareto e indicadores de qualidade)
estão disponíveis publicamente em \url{https://github.com/Lost-Memory/mastersdegree},
sob licença Creative Commons Attribution~4.0 (dados) e MIT (código), e arquivados de
forma permanente sob o identificador DOI \texttt{10.5281/zenodo.21481074}. Os valores de custo e GWP das quatro regiões-âncora
derivam de inventário de Avaliação de Ciclo de Vida (base EcoInvent, \textit{software}
SimaPro, método IPCC GWP-100); apenas os valores agregados são disponibilizados, não
a base EcoInvent, cujo uso requer licença própria.
```

## Passo a passo para publicar (resumo)

1. Crie um repositório no **GitHub** e suba esta pasta (`git init`, `git add .`,
   `git commit`, `git push`).
2. Em **zenodo.org**, faça login com a conta do GitHub e ative o repositório na aba
   *GitHub* do Zenodo.
3. No GitHub, crie um **Release** (ex.: `v1.0.0`). O Zenodo arquiva automaticamente e
   **gera o DOI**.
4. Copie o DOI para: (a) `CITATION.cff` (campo `doi`), (b) o `README.md`, e (c) o
   parágrafo acima na dissertação.
5. **Antes de publicar**, confirme os termos da sua licença EcoInvent/SimaPro quanto
   à divulgação dos valores agregados das âncoras (ver nota em `LICENSE-DATA.md`).
