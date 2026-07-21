function Estatistica_REO_V5
% =====================================================================
%  Estatistica_REO_V5 — EXPERIMENTO DE VALIDAÇÃO ESTATÍSTICA
% =====================================================================
%  Roda as técnicas ESTOCÁSTICAS principais NSEEDS vezes (sementes fixas
%  e reprodutíveis) para cada valor de K, guardando a frente de Pareto
%  final de CADA execução. As métricas (HV/IGD) e os testes estatísticos
%  (Wilcoxon/Friedman) são calculados depois, em Python, sobre estas saídas.
%
%  POR QUE ESTE ARQUIVO SEPARADO
%   - O Master_Optimizer roda as 17 técnicas de uma vez; com 14 núcleos
%     isso satura o pool e trava. Aqui roda-se UMA técnica por vez, com as
%     30 sementes em parfor sobre um pool LIMITADO -> memória previsível,
%     sem travamento.
%   - Usa DOMINÂNCIA VETORIZADA (10-50x mais rápida que a O(n^2); idêntica
%     em resultado — verificada por auto-teste no início). É o que torna
%     viável rodar ~900 execuções.
%   - RESUMÍVEL: cada par (K, técnica) vira um checkpoint .mat. Pode fechar
%     o MATLAB e reabrir; ao rodar de novo, ele PULA o que já terminou.
%
%  COMO RODAR (fire-and-forget)
%   1. Abra o MATLAB nesta pasta (com Custos.csv e GWP.csv ao lado).
%   2. Comando:  >> Estatistica_REO_V5
%   3. Deixe rodar. Ao terminar, a pasta estatistica_out/ terá:
%        - ck_K*_*.mat            (checkpoints por par K×técnica)
%        - estatistica_resultados.csv  (frentes de todas as execuções)
%        - estatistica_log.txt    (diário da execução, com ETA)
%
%  AJUSTES por variável de ambiente (para TESTE, sem editar o arquivo):
%     set ESTAT_NSEEDS=3        %% menos sementes
%     set ESTAT_KLIST=13        %% só K=13 (irrestrito)
%     set ESTAT_TECHS=NSGA2     %% só uma técnica
%     set ESTAT_MAXIT=100       %% MaxIt reduzido (NUNCA no experimento final)
%  (no experimento definitivo, NÃO defina nenhuma dessas variáveis.)
% =====================================================================
clc;

%% ---------------- Configuração ----------------
NSEEDS = getenvnum('ESTAT_NSEEDS', 30);
KLIST  = getenvvec('ESTAT_KLIST',  [2 5 13]);
TECHS  = getenvcell('ESTAT_TECHS', {'MO_MCSHADE','MO_MCDE','NSGA2','MOEAD_AV','MOPSO'});
MAXIT_OVR = getenvnum('ESTAT_MAXIT', 0);

OUTDIR = fullfile(pwd, 'estatistica_out');
if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end
diary(fullfile(OUTDIR, 'estatistica_log.txt')); diary on;
fprintf('=== Estatistica_REO_V5 ===\n');
fprintf('%d sementes | K = %s | tecnicas: %s\n', NSEEDS, mat2str(KLIST), strjoin(TECHS, ','));

%% ---------------- Parâmetros (IGUAIS ao benchmark principal) ----------------
P.MOPSO      = struct('MaxIt',2000,'nPop',500,'nRep',500,'w',1.0,'wdamp',0.99,'c1',1.5,'c2',1.5,'mu',0.5);
P.NSGA2      = struct('MaxIt',2000,'nPop',500,'pCrossover',0.8,'pMutation',0.2);
P.MOEAD_AV   = struct('MaxIt',2000,'nPopInit',500,'T',10);
P.MO_MCDE    = struct('MaxIt',500,'nPop',500,'M',4,'F_L',0.3,'F_U',1.6,'Cr',0.85,'nRep',500);
P.MO_MCSHADE = struct('MaxIt',500,'nPop',500,'M',4,'H',round(500/3),'c',0.1,'nRep',500);
if MAXIT_OVR > 0
    fn = fieldnames(P);
    for i = 1:numel(fn), P.(fn{i}).MaxIt = MAXIT_OVR; end
    fprintf('*** ATENCAO: MaxIt reduzido para %d (modo TESTE) ***\n', MAXIT_OVR);
end

%% ---------------- Dados ----------------
fprintf('Lendo Custos.csv / GWP.csv...\n');
dados = Dados_Processo_REO();

%% ---------------- Auto-teste da dominância vetorizada ----------------
selftest_domination();

%% ---------------- Pool paralelo (LIMITADO — evita o travamento) ----------------
nW = max(1, min(feature('numcores') - 2, 12));
pool = gcp('nocreate');
if isempty(pool) || pool.NumWorkers ~= nW
    delete(gcp('nocreate')); pool = parpool(nW);
end
fprintf('Pool com %d workers.\n', nW);
dq = parallel.pool.DataQueue; afterEach(dq, @(~) []);   % silencioso (evita spam de 900 execucoes)

%% ---------------- Loop principal (RESUMÍVEL) ----------------
totalPairs = numel(KLIST) * numel(TECHS); donePairs = 0; tExp = tic;
for K = KLIST
    d = dados; d.K = K;
    for ti = 1:numel(TECHS)
        tech = TECHS{ti};
        ckpt = fullfile(OUTDIR, sprintf('ck_K%d_%s.mat', K, tech));
        if exist(ckpt, 'file')
            fprintf('[pular] K=%2d  %-11s (checkpoint ja existe)\n', K, tech);
            donePairs = donePairs + 1; continue;
        end
        fprintf('==> K=%2d  %-11s  rodando %d sementes ... ', K, tech, NSEEDS);
        par = P.(tech);
        fronts = cell(NSEEDS, 1); tempos = zeros(NSEEDS, 1); nds = zeros(NSEEDS, 1);
        tPair = tic;
        parfor s = 1:NSEEDS
            rng(s, 'twister');                 % semente fixa e reprodutível
            [rep, tt] = run_tech(tech, d, par, dq);
            C = costs_of(rep);
            fronts{s} = C; tempos(s) = tt; nds(s) = size(C, 1);
        end
        save(ckpt, 'fronts', 'tempos', 'nds', 'K', 'tech', 'NSEEDS', 'par', '-v7');
        donePairs = donePairs + 1;
        eta = (toc(tExp) / donePairs) * (totalPairs - donePairs);
        fprintf('%.0f s | mediana N_nd = %d | ETA restante ~ %.1f h\n', ...
                toc(tPair), round(median(nds)), eta/3600);
    end
end

%% ---------------- Consolida em CSV longo ----------------
consolidate_csv(OUTDIR);
fprintf('\n=== CONCLUIDO === resultados em: %s\n', OUTDIR);
diary off;
end

% =====================================================================
%   Infra do experimento
% =====================================================================
function [rep, tt] = run_tech(tech, d, par, dq)
% Despacho por nome (parfor exige função, não handle de campo de struct).
switch tech
    case 'MOPSO',      [rep, tt] = RunMOPSO(d, par, dq);
    case 'NSGA2',      [rep, tt] = RunNSGA2(d, par, dq);
    case 'MOEAD_AV',   [rep, tt] = RunMOEAD_AV(d, par, dq);
    case 'MO_MCDE',    [rep, tt] = RunMO_MCDE(d, par, dq);
    case 'MO_MCSHADE', [rep, tt] = RunMO_MCSHADE(d, par, dq);
    otherwise, error('Tecnica desconhecida: %s', tech);
end
end

function C = costs_of(rep)
% Extrai a matriz nSol x 2 (custo, GWP) de um repositório de soluções.
if isempty(rep) || ~isfield(rep, 'Cost'), C = zeros(0, 2); return; end
C = horzcat(rep.Cost)';
C = C(all(~isnan(C), 2), :);
end

function selftest_domination()
% Confirma que as rotinas VETORIZADAS (usadas aqui) são idênticas às O(n^2)
% originais do Master_Optimizer. Se divergir, ABORTA (não vale rodar um
% experimento longo com uma rotina incorreta). Testa dois cenários: custos
% contínuos (sem empates) e custos inteiros (COM empates, como no dado real).
rng(12345, 'twister');
for caso = 1:2
    n = 300; pop = repmat(struct('Cost', [], 'IsDominated', false, 'Rank', []), n, 1);
    for i = 1:n
        if caso == 1, pop(i).Cost = rand(2, 1) * 100;        % contínuo
        else,         pop(i).Cost = round(rand(2, 1) * 6);   % inteiro -> muitos empates
        end
    end
    % (1) DetermineDomination vetorizada x O(n^2)
    a = DetermineDomination(pop); b = DetermineDomination_ref(pop);
    if ~isequal(logical([a.IsDominated]), logical([b.IsDominated]))
        error('AUTO-TESTE FALHOU (caso %d): DetermineDomination vetorizada != O(n^2).', caso);
    end
    % (2) NonDominatedSort vetorizada x O(n^2) -> os RANKS devem ser idênticos
    c = NonDominatedSort(pop); d = NonDominatedSort_ref(pop);
    if ~isequal([c.Rank], [d.Rank])
        error('AUTO-TESTE FALHOU (caso %d): NonDominatedSort vetorizada != O(n^2).', caso);
    end
end
fprintf('Auto-teste (DetermineDomination + NonDominatedSort, com e sem empates): OK.\n');
end

function consolidate_csv(OUTDIR)
% Junta todos os checkpoints em um CSV longo: uma linha por solução.
files = dir(fullfile(OUTDIR, 'ck_K*_*.mat'));
fid = fopen(fullfile(OUTDIR, 'estatistica_resultados.csv'), 'w');
fprintf(fid, 'K,Algoritmo,Rep,Custo_Total,Impacto_Ambiental,Tempo_Total_Segundos\n');
nLin = 0;
for f = 1:numel(files)
    S = load(fullfile(OUTDIR, files(f).name));
    for s = 1:numel(S.fronts)
        C = S.fronts{s};
        for r = 1:size(C, 1)
            fprintf(fid, '%d,%s,%d,%.6f,%.6f,%.3f\n', S.K, S.tech, s, C(r,1), C(r,2), S.tempos(s));
            nLin = nLin + 1;
        end
    end
end
fclose(fid);
fprintf('CSV consolidado: %d checkpoints, %d linhas -> estatistica_resultados.csv\n', numel(files), nLin);
end

% ---- leitura de variáveis de ambiente (para o modo teste) ----
function v = getenvnum(name, def)
s = getenv(name); if isempty(s), v = def; else, v = str2double(s); if isnan(v), v = def; end; end
end
function v = getenvvec(name, def)
s = getenv(name); if isempty(s), v = def; else, v = str2num(strrep(s, ',', ' ')); if isempty(v), v = def; end; end %#ok<ST2NM>
end
function v = getenvcell(name, def)
s = getenv(name); if isempty(s), v = def; else, v = strtrim(strsplit(s, ',')); end
end

% =====================================================================
%   Dominância VETORIZADA (10-50x mais rápida; substitui a O(n^2)).
%   A O(n^2) original está preservada como DetermineDomination_ref para
%   o auto-teste de equivalência.
% =====================================================================
function pop = DetermineDomination(pop)
n = numel(pop);
for i = 1:n, pop(i).IsDominated = false; end
if n < 2, return; end
Costs = horzcat(pop.Cost)';                 % n x nObj
for i = 1:n
    le = all(Costs <= Costs(i, :), 2);      % outros <= i em todos os objetivos
    lt = any(Costs <  Costs(i, :), 2);      % outros <  i em algum objetivo
    dom = le & lt; dom(i) = false;
    if any(dom), pop(i).IsDominated = true; end
end
end

% =====================================================================
%   Ordenação por não-dominância VETORIZADA (fronts do NSGA-II).
%   Substitui a NonDominatedSort O(n^2) (o real gargalo do NSGA-II e das
%   variantes multi-child, que também a usam na seleção ambiental).
%   Produz RANKS idênticos aos da versão O(n^2) (preservada como
%   NonDominatedSort_ref para o auto-teste). Complexidade: uma matriz de
%   dominância n x n construída em O(n^2) mas VETORIZADA (sem acesso a
%   struct em laço duplo), seguida de peeling de fronts sobre vetores.
% =====================================================================
function pop = NonDominatedSort(pop)
n = numel(pop);
if n == 0, return; end
Costs = horzcat(pop.Cost)';                 % n x nObj
% Matriz de dominância D(i,j) = "i domina j" (minimização).
D = false(n, n);
for i = 1:n
    le = all(Costs(i, :) <= Costs, 2);      % i <= cada j em todos os objetivos
    lt = any(Costs(i, :) <  Costs, 2);      % i <  cada j em algum objetivo
    D(i, :) = (le & lt)';
end
domCount = sum(D, 1)';                       % nº de soluções que dominam j
rank = zeros(n, 1); remaining = true(n, 1); cur = 1;
while any(remaining)
    front = remaining & (domCount == 0);
    if ~any(front), front = remaining; end   % salvaguarda (dominância não tem ciclos)
    rank(front) = cur;
    remaining(front) = false;
    domCount = domCount - sum(D(front, :), 1)';   % remove a pressão dos que saíram
    cur = cur + 1;
end
for i = 1:n, pop(i).Rank = rank(i); end
end

% =====================================================================
%   Funções de técnica e auxiliares — EXTRAÍDAS VERBATIM do
%   Master_Optimizer_REO_V5_MC.m (garante fidelidade ao benchmark).
% =====================================================================

function [repValido, tempo_execucao] = RunNSGA2(dados, params, dq)
% NSGA-II — ordenação por não-dominância (fronts) + distância de aglomeração
% (crowding) para diversidade. Seleção por torneio, cruzamento uniforme e
% mutação; elitismo via combinação pais+filhos e truncamento por (rank, crowding).
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar];
    empty.Position = []; empty.Cost = []; empty.IsDominated = false; empty.Rank = []; empty.CrowdingDistance = []; pop = repmat(empty, params.nPop, 1);
    for i = 1:params.nPop, pop(i).Position = round(unifrnd(VarMin, VarMax, VarSize)); pop(i).Cost = CostFunction(pop(i).Position, dados); end
    pop = NonDominatedSort(pop); pop = CalcCrowdingDistance(pop);
    for it = 1:params.MaxIt
        popc = repmat(empty, params.nPop, 1);
        for i = 1:2:params.nPop
            % Torneio binário por (rank, crowding); cruzamento uniforme por máscara.
            p1 = TournamentSelection(pop); p2 = TournamentSelection(pop);
            if rand < params.pCrossover, mask = rand(VarSize) > 0.5; popc(i).Position = pop(p1).Position; popc(i).Position(mask) = pop(p2).Position(mask); popc(i+1).Position = pop(p2).Position; popc(i+1).Position(mask) = pop(p1).Position(mask); else, popc(i).Position = pop(p1).Position; popc(i+1).Position = pop(p2).Position; end
            if rand < params.pMutation, popc(i).Position = Mutate(popc(i).Position, 0.1, VarMin, VarMax); end; if rand < params.pMutation, popc(i+1).Position = Mutate(popc(i+1).Position, 0.1, VarMin, VarMax); end
            popc(i).Position = round(popc(i).Position); popc(i+1).Position = round(popc(i+1).Position); popc(i).Cost = CostFunction(popc(i).Position, dados); popc(i+1).Cost = CostFunction(popc(i+1).Position, dados); popc(i).IsDominated = false; popc(i+1).IsDominated = false;
        end
        % Elitismo: junta pais+filhos, reordena e trunca para nPop.
        popTotal = [pop; popc]; popTotal = NonDominatedSort(popTotal); popTotal = CalcCrowdingDistance(popTotal); pop = SortAndTruncateNSGA2(popTotal, params.nPop);

        if mod(it, 100) == 0 || it == params.MaxIt
            send(dq, sprintf('NSGA-II    : Iter %4d / %d concluida', it, params.MaxIt));
        end
    end
    tempo_execucao = toc(t0); repValido = pop([pop.Rank] == 1);   % retorna o 1º front
end

function [repValido, tempo_execucao] = RunMOPSO(dados, params, dq)
% MOPSO — Multi-Objective Particle Swarm Optimization com GRADE ADAPTATIVA.
% Mantém um repositório externo (rep) de soluções não-dominadas; o líder de
% cada partícula é sorteado de regiões POUCO povoadas da grade (preserva
% diversidade ao longo da frente). Mutação/turbulência evita ótimos locais.
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar];
    % Estrutura de uma partícula (posição, velocidade, custo, melhor pessoal e índices de grade).
    empty.Position = []; empty.Velocity = []; empty.Cost = []; empty.Best.Position = []; empty.Best.Cost = []; empty.IsDominated = false; empty.GridIndex = []; empty.GridSubIndex = [];
    pop = repmat(empty, params.nPop, 1);
    % Inicialização aleatória uniforme no domínio [1, nForn].
    for i = 1:params.nPop, pop(i).Position = unifrnd(VarMin, VarMax, VarSize); pop(i).Velocity = zeros(VarSize); pop(i).Cost = CostFunction(pop(i).Position, dados); pop(i).Best = pop(i); end
    % Inicializa o repositório com os não-dominados e cria a grade adaptativa.
    pop = DetermineDomination(pop); novos = pop(~[pop.IsDominated]); rep = repmat(empty, params.nRep, 1); rep(1:numel(novos)) = novos; repCount = numel(novos);
    Grid = CreateGrid(rep, repCount, 7, 0.1); for i = 1:repCount, rep(i) = FindGridIndex(rep(i), Grid); end
    w_atual = params.w;
    for it = 1:params.MaxIt
        % pm: probabilidade de mutação decrescente (mais exploração no início).
        pm = (1 - (it-1)/(params.MaxIt-1))^(1/params.mu);
        for i = 1:params.nPop
            leader = SelectLeader(rep, repCount, 2);   % líder de região esparsa
            % Atualização clássica de velocidade: inércia + componente cognitivo + social.
            pop(i).Velocity = w_atual * pop(i).Velocity + params.c1 * rand(VarSize) .* (pop(i).Best.Position - pop(i).Position) + params.c2 * rand(VarSize) .* (leader.Position - pop(i).Position);
            % MELHORIA SUGERIDA (Vmax) — evita explosão de velocidade com w=1,0:
            %   Vmax = 0.2*(VarMax-VarMin); pop(i).Velocity = max(min(pop(i).Velocity, Vmax), -Vmax);
            % (desativada para não alterar os resultados já validados; ative se notar instabilidade)
            pop(i).Position = min(max(pop(i).Position + pop(i).Velocity, VarMin), VarMax); pop(i).Cost = CostFunction(pop(i).Position, dados);
            % Operador de mutação/turbulência (aceita por dominância).
            if rand < pm
                newPos = Mutate(pop(i).Position, pm, VarMin, VarMax); newCost= CostFunction(newPos, dados);
                if Dominates(newCost, pop(i).Cost), pop(i).Position = newPos; pop(i).Cost = newCost; else, if ~Dominates(pop(i).Cost, newCost) && rand < 0.5, pop(i).Position = newPos; pop(i).Cost = newCost; end; end
            end
            % Atualiza o melhor pessoal (Pbest) por dominância (com desempate aleatório).
            if Dominates(pop(i), pop(i).Best), pop(i).Best = pop(i); else, if ~Dominates(pop(i).Best, pop(i)) && rand < 0.5, pop(i).Best = pop(i); end; end
        end
        % Insere novos não-dominados no repositório, re-filtra e poda excesso.
        pop = DetermineDomination(pop); novos = pop(~[pop.IsDominated]); k = numel(novos);
        if repCount + k > params.nRep, k = params.nRep - repCount; end
        if k > 0, rep(repCount+1:repCount+k) = novos(1:k); repCount = repCount + k; end
        tmp = rep(1:repCount); tmp = DetermineDomination(tmp); tmp = tmp(~[tmp.IsDominated]); rep(1:numel(tmp)) = tmp; repCount = numel(tmp);
        % [CORRIGIDO] agora recebe 'rep' de volta — a poda por densidade de
        % grade realmente altera o repositório (ver nota em DeleteOneRepMember).
        while repCount > params.nRep, [rep, repCount] = DeleteOneRepMember(rep, repCount, 2); end
        Grid = CreateGrid(rep, repCount, 7, 0.1); for i = 1:repCount, rep(i) = FindGridIndex(rep(i), Grid); end
        w_atual = w_atual * params.wdamp;   % decaimento da inércia

        if mod(it, 100) == 0 || it == params.MaxIt
            send(dq, sprintf('MOPSO      : Iter %4d / %d concluida', it, params.MaxIt));
        end
    end
    tempo_execucao = toc(t0); repValido = rep(1:repCount);
end

function [repValido, tempo_execucao] = RunMOEAD_AV(dados, params, dq)
% MOEA/D-AV — variante com população ADAPTATIVA (Adaptive Volume): N cresce
% se o arquivo está rico e encolhe se está pobre, recalculando pesos/vizinhança.
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar]; nObj = 2; N = params.nPopInit; W = zeros(N, nObj);
    for i = 1:N, W(i, 1) = (i-1)/(N-1); W(i, 2) = 1 - W(i, 1); end; T = min(params.T, N); sp = struct('V', [], 'Neighbors', []); sp = repmat(sp, N, 1);
    for i = 1:N, D = pdist2(W(i,:), W); [~, idx] = sort(D); sp(i).Neighbors = idx(1:T); end; empty_ind.Position = []; empty_ind.Cost = []; empty_ind.IsDominated = false; pop = repmat(empty_ind, N, 1); z = inf(nObj, 1);
    for i = 1:N, pop(i).Position = round(unifrnd(VarMin, VarMax, VarSize)); pop(i).Cost = CostFunction(pop(i).Position, dados); z = min(z, pop(i).Cost); end; EP = pop;
    for it = 1:params.MaxIt
        for i = 1:N
            neighbors = sp(i).Neighbors; p1 = neighbors(randi(T)); p2 = neighbors(randi(T)); y = empty_ind; y.Position = pop(p1).Position; mask = rand(VarSize) > 0.5; y.Position(mask) = pop(p2).Position(mask); y.Position = round(Mutate(y.Position, 0.1, VarMin, VarMax));
            y.Cost = CostFunction(y.Position, dados); y.IsDominated = false; z = min(z, y.Cost);
            for j = 1:T, k = neighbors(j); if k <= size(W, 1) && k <= length(pop), g_old = max(W(k,:)' .* abs(pop(k).Cost - z)); g_new = max(W(k,:)' .* abs(y.Cost - z)); if g_new <= g_old, pop(k) = y; end; end; end
            EP = [EP; y]; %#ok<AGROW>
        end
        EP = DetermineDomination(EP); EP = EP(~[EP.IsDominated]); if ~isempty(EP), costs = round(horzcat(EP.Cost)', 4); [~, uniqueIdx] = unique(costs, 'rows'); EP = EP(uniqueIdx); end
        % Ajuste adaptativo do tamanho da população a cada 20 iterações.
        if mod(it, 20) == 0 && N > 5 && N < 500
            if length(EP) > N * 0.5, N = N + round(N * 0.1); elseif length(EP) < N * 0.2, N = N - round(N * 0.1); end
            W = zeros(N, nObj); for i = 1:N, W(i, 1) = (i-1)/(N-1); W(i, 2) = 1 - W(i, 1); end; T = min(params.T, N); sp = repmat(struct('V', [], 'Neighbors', []), N, 1);
            for i = 1:N, D = pdist2(W(i,:), W); [~, idx] = sort(D); sp(i).Neighbors = idx(1:T); end
            if length(pop) < N, new_pop = repmat(empty_ind, N - length(pop), 1); for m = 1:length(new_pop), new_pop(m).Position = round(unifrnd(VarMin, VarMax, VarSize)); new_pop(m).Cost = CostFunction(new_pop(m).Position, dados); end; pop = [pop; new_pop]; else, pop = pop(1:N); end
        end

        if mod(it, 100) == 0 || it == params.MaxIt
            send(dq, sprintf('MOEA/D-AV  : Iter %4d / %d concluida', it, params.MaxIt));
        end
    end
    tempo_execucao = toc(t0); repValido = EP;
end

function [repValido, tempo_execucao] = RunMO_MCDE(dados, params, dq)
% MO-MCDE â€” Multi-Child DE multiobjetivo.
% NÃºcleo MCDE preservado (eqs. 1-5 do paper CEC 2026):
%   - mutaÃ§Ã£o DE/rand/1 com r0/r1/r2 SORTEADOS POR FILHO;
%   - dithering de F por competiÃ§Ã£o pai-filhos (eq. 2):
%       F_i = F_L + rand*(F_U - F_L), F_L=0,3, F_U=1,6;
%   - crossover binomial com Cr fixo (0,85) e j_rand (eq. 3);
%   - M filhos por pai; prÃ©-seleÃ§Ã£o -> campeÃ£o w_i (eq. 4);
%   - seleÃ§Ã£o final pai x campeÃ£o (eq. 5), aqui por dominÃ¢ncia.
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar];
    N = params.nPop; M = params.M;
    empty.Position = []; empty.Cost = []; empty.IsDominated = false; empty.Rank = []; empty.CrowdingDistance = [];
    pop = repmat(empty, N, 1);
    for i = 1:N, pop(i).Position = round(unifrnd(VarMin, VarMax, VarSize)); pop(i).Cost = CostFunction(pop(i).Position, dados); end
    popIni = NonDominatedSort(pop);
    rep = UpdateRepositorio(repmat(empty, 0, 1), popIni([popIni.Rank] == 1), params.nRep);
    for it = 1:params.MaxIt
        overflow = repmat(empty, 0, 1);
        newPop = pop;
        for i = 1:N
            % dithering de F: um F por competiÃ§Ã£o pai-filhos (eq. 2)
            Fi = params.F_L + rand*(params.F_U - params.F_L);
            filhos = repmat(empty, M, 1);
            for m = 1:M
                r = DistinctIdx(N, i, 3);           % r0, r1, r2 por FILHO
                v = pop(r(1)).Position + Fi*(pop(r(2)).Position - pop(r(3)).Position);
                v = min(max(v, VarMin), VarMax);
                u = pop(i).Position; jrand = randi(nVar);
                cmask = rand(VarSize) <= params.Cr; cmask(jrand) = true;
                u(cmask) = v(cmask); u = round(u);
                filhos(m).Position = u; filhos(m).Cost = CostFunction(u, dados);
            end
            w = PreSelecaoMO(filhos, pop(i));       % campeÃ£o entre os M filhos
            if Dominates(w.Cost, pop(i).Cost)
                newPop(i) = w;                      % campeÃ£o substitui o pai
            elseif ~Dominates(pop(i).Cost, w.Cost)
                overflow(end+1, 1) = w;             %#ok<AGROW> % nÃ£o-dominados entre si
            end
        end
        % seleÃ§Ã£o ambiental (GDE3): pop + excedentes -> nPop por (rank, crowding)
        pool_ = [newPop; overflow];
        pool_ = NonDominatedSort(pool_); pool_ = CalcCrowdingDistance(pool_);
        pop = SortAndTruncateNSGA2(pool_, N);
        % repositÃ³rio: basta inserir os rank-1 da geraÃ§Ã£o (os demais seriam
        % filtrados de qualquer forma dentro de UpdateRepositorio)
        rep = UpdateRepositorio(rep, pool_([pool_.Rank] == 1), params.nRep);
        if mod(it, 50) == 0 || it == params.MaxIt
            send(dq, sprintf('MO-MCDE    : Iter %4d / %d concluida (rep = %d)', it, params.MaxIt, numel(rep)));
        end
    end
    tempo_execucao = toc(t0); repValido = rep;
end

function [repValido, tempo_execucao] = RunMO_MCSHADE(dados, params, dq)
% MO-MC-SHADE â€” Multi-Child SHADE multiobjetivo.
% NÃºcleo MC-SHADE preservado (eqs. 6-14 do paper CEC 2026):
%   - mutaÃ§Ã£o DE/current-to-pbest (eq. 6): v = x_i + F*(x_pbest - x_i + x_r1 - x_r2),
%     com x_pbest sorteado entre os p melhores da geraÃ§Ã£o anterior,
%     p ~ U{2, floor(0,2N)}, mantido para os M filhos do mesmo pai;
%   - Cr_i ~ Normal(mu_Cr, 0,1) e F_i ~ Cauchy(mu_F, 0,1) POR COMPETIÃ‡ÃƒO (eq. 8);
%   - buffer circular de tamanho H = N/3 (anel), inicializado em 0,5; a cada
%     geraÃ§Ã£o um slot aleatÃ³rio fornece (mu_Cr, mu_F) e o slot corrente do
%     anel recebe a atualizaÃ§Ã£o:
%       mu_Cr novo = mu_Cr*(1-c) + c * [mÃ©dia aritmÃ©tica de Cr ponderada por delta] (eqs. 9-10)
%       mu_F  novo = mu_F *(1-c) + c * [mÃ©dia de Lehmer quadrÃ¡tica dos F de sucesso] (eqs. 12-13)
%     com c = 0,1 e delta_i pela regra de Peng (eq. 11);
%   - SEM archive externo de derrotados (fiel ao MC-SHADE original);
%   - tratamento de borda padrÃ£o SHADE: ponto mÃ©dio entre o limite e o pai.
% AdaptaÃ§Ã£o multiobjetivo: "p melhores" = ordenaÃ§Ã£o (rank, crowding);
% sucesso = campeÃ£o domina o pai; delta_i = norma da melhoria normalizada.
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar];
    N = params.nPop; M = params.M; H = max(2, params.H); c = params.c;
    empty.Position = []; empty.Cost = []; empty.IsDominated = false; empty.Rank = []; empty.CrowdingDistance = [];
    pop = repmat(empty, N, 1);
    for i = 1:N, pop(i).Position = round(unifrnd(VarMin, VarMax, VarSize)); pop(i).Cost = CostFunction(pop(i).Position, dados); end
    popIni = NonDominatedSort(pop);
    rep = UpdateRepositorio(repmat(empty, 0, 1), popIni([popIni.Rank] == 1), params.nRep);
    bufCr = 0.5*ones(H, 1); bufF = 0.5*ones(H, 1); hPos = 1;
    for it = 1:params.MaxIt
        % ordenaÃ§Ã£o (rank, crowding) da geraÃ§Ã£o anterior p/ o operador pbest
        pop = NonDominatedSort(pop); pop = CalcCrowdingDistance(pop);
        [~, ordem] = sortrows([[pop.Rank]', -[pop.CrowdingDistance]']);
        % (mu_Cr, mu_F) da geraÃ§Ã£o: slot aleatÃ³rio do buffer circular
        rIdx = randi(H); mu_Cr = bufCr(rIdx); mu_F = bufF(rIdx);
        % faixas correntes dos objetivos p/ normalizar o peso de sucesso
        C_all = horzcat(pop.Cost)'; rngObj = max(C_all, [], 1) - min(C_all, [], 1); rngObj(rngObj == 0) = 1;
        S_F = []; S_Cr = []; S_d = [];
        overflow = repmat(empty, 0, 1); newPop = pop;
        pmax = max(2, floor(0.2*N));
        for i = 1:N
            Cr_i = min(1, max(0, mu_Cr + 0.1*randn));
            F_i = 0; while F_i <= 0, F_i = mu_F + 0.1*tan(pi*(rand - 0.5)); end
            F_i = min(F_i, 1);
            psz = randi([2, pmax]); pb = ordem(randi(psz));   % pbest do pai (vale p/ os M filhos)
            filhos = repmat(empty, M, 1);
            for m = 1:M
                rr = DistinctIdx(N, i, 2);
                v = pop(i).Position + F_i*(pop(pb).Position - pop(i).Position + pop(rr(1)).Position - pop(rr(2)).Position);
                low = v < VarMin;  v(low)  = (VarMin + pop(i).Position(low))/2;
                high = v > VarMax; v(high) = (VarMax + pop(i).Position(high))/2;
                u = pop(i).Position; jrand = randi(nVar);
                cmask = rand(VarSize) <= Cr_i; cmask(jrand) = true;
                u(cmask) = v(cmask); u = round(u);
                filhos(m).Position = u; filhos(m).Cost = CostFunction(u, dados);
            end
            w = PreSelecaoMO(filhos, pop(i));
            if Dominates(w.Cost, pop(i).Cost)
                dvec = (pop(i).Cost(:).' - w.Cost(:).') ./ rngObj;   % melhoria normalizada
                S_F(end+1, 1) = F_i; S_Cr(end+1, 1) = Cr_i; S_d(end+1, 1) = norm(dvec); %#ok<AGROW>
                newPop(i) = w;
            elseif ~Dominates(pop(i).Cost, w.Cost)
                overflow(end+1, 1) = w; %#ok<AGROW>
            end
        end
        % atualizaÃ§Ã£o do buffer circular (sÃ³ se houve sucesso â€” eq. 10/13)
        if ~isempty(S_F) && sum(S_d) > 0
            CrW = sum(S_d .* S_Cr) / sum(S_d);                    % eq. 10 (ponderada por delta)
            FLeh = sum(S_F.^2) / max(1e-12, sum(S_F));            % eq. 13 (Lehmer quadrÃ¡tica)
            bufCr(hPos) = mu_Cr*(1 - c) + c*CrW;                  % eq. 9
            bufF(hPos)  = mu_F *(1 - c) + c*FLeh;                 % eq. 12
            hPos = hPos + 1; if hPos > H, hPos = 1; end           % anel
        end
        pool_ = [newPop; overflow];
        pool_ = NonDominatedSort(pool_); pool_ = CalcCrowdingDistance(pool_);
        pop = SortAndTruncateNSGA2(pool_, N);
        rep = UpdateRepositorio(rep, pool_([pool_.Rank] == 1), params.nRep);
        if mod(it, 50) == 0 || it == params.MaxIt
            send(dq, sprintf('MO-MC-SHADE: Iter %4d / %d concluida (rep = %d)', it, params.MaxIt, numel(rep)));
        end
    end
    tempo_execucao = toc(t0); repValido = rep;
end

function z = CostFunction(x, dados)
% Avalia os dois objetivos de uma solução x (vetor de índices de fornecedor
% por material). Passos:
%   1) Discretiza x para inteiros válidos em [1, nForn] (PSO opera contínuo).
%   2) [NOVO] Se a restrição K está ativa, REPARA para <= K fornecedores
%      distintos (reparo Baldwiniano: vale para todas as técnicas).
%   3) Soma o custo e o GWP do fornecedor escolhido de cada material.
    idx_fornecedores = max(1, min(round(x), dados.nFornecedores));
    if isfield(dados,'K') && ~isempty(dados.K) && dados.K < dados.nVar
        idx_fornecedores = ApplyK(idx_fornecedores, dados);
    end
    linhas = 1:dados.nVar;
    indices_lineares = sub2ind(size(dados.MatrizCustos), linhas, idx_fornecedores);
    z = [sum(dados.MatrizCustos(indices_lineares)); sum(dados.MatrizGWP(indices_lineares))];
end

function idx = ApplyK(idx, dados)
% [NOVO] Operador de REPARO da restrição de consolidação: reduz o número de
% fornecedores DISTINTOS em 'idx' para no máximo dados.K.
% Heurística (gulosa e determinística):
%   enquanto há mais de K fornecedores distintos:
%     - identifica o fornecedor MENOS usado (candidato a ser eliminado);
%     - reatribui cada material que o usava ao fornecedor REMANESCENTE de
%       menor custo+GWP NORMALIZADOS para aquele material (Cn+Gn).
% Motivo da normalização: custo (~10^2-10^3) e GWP (~10^0) têm escalas muito
% diferentes; normalizar por material evita que o custo domine a reatribuição.
    K = dados.K;
    u = unique(idx);
    while numel(u) > K
        counts = arrayfun(@(s) sum(idx == s), u);   % uso de cada fornecedor
        [~, mi] = min(counts); drop = u(mi);         % elimina o menos usado
        keep = u(u ~= drop);
        mats = find(idx == drop);
        for m = mats
            sc = dados.Cn(m, keep) + dados.Gn(m, keep);  % custo normalizado agregado
            [~, bi] = min(sc); idx(m) = keep(bi);
        end
        u = unique(idx);
    end
    % NOTA: este é o reparo "Baldwiniano" (usado na AVALIAÇÃO). Para a variante
    % "Lamarckiana" (escrever o reparo de volta no genótipo, melhor convergência),
    % chamar ApplyK também após cada atualização de posição dentro das técnicas.
end

function dados = Dados_Processo_REO()
% Lê Custos.csv e GWP.csv (formato: 1ª coluna = nome do material; demais =
% fornecedores). Converte vírgula decimal para ponto, descarta linhas com
% valores inválidos (NaN) e monta a struct 'dados'.
% [NOVO] Calcula também Cn/Gn (matrizes normalizadas por material em [0,1]),
% usadas pelo operador de reparo ApplyK.
    opts = detectImportOptions('Custos.csv'); opts.VariableNamingRule = 'preserve'; opts = setvartype(opts, opts.VariableNames(2:end), 'string'); tabCusto = readtable('Custos.csv', opts);
    if ~isfile('GWP.csv'), error('GWP.csv nao encontrado.'); end
    optsGWP = detectImportOptions('GWP.csv'); optsGWP.VariableNamingRule = 'preserve'; optsGWP = setvartype(optsGWP, optsGWP.VariableNames(2:end), 'string'); tabGWP = readtable('GWP.csv', optsGWP);
    minHeight = min(height(tabCusto), height(tabGWP)); tabCusto = tabCusto(1:minHeight, :); tabGWP = tabGWP(1:minHeight, :);
    matCusto = str2double(strrep(tabCusto{:, 2:end}, ',', '.')); matGWP = str2double(strrep(tabGWP{:, 2:end}, ',', '.'));
    linhasInvalidas = any(isnan(matCusto), 2) | any(isnan(matGWP), 2); matCusto = matCusto(~linhasInvalidas, :); matGWP = matGWP(~linhasInvalidas, :);
    dados.nVar = size(matCusto, 1); dados.nFornecedores = width(tabCusto) - 1; dados.MatrizCustos = matCusto; dados.MatrizGWP = matGWP;
    % Normalização min-max POR MATERIAL (linha); evita divisão por zero.
    rangeC = max(matCusto, [], 2) - min(matCusto, [], 2); rangeC(rangeC == 0) = 1;
    rangeG = max(matGWP,   [], 2) - min(matGWP,   [], 2); rangeG(rangeG == 0) = 1;
    dados.Cn = (matCusto - min(matCusto, [], 2)) ./ rangeC;
    dados.Gn = (matGWP   - min(matGWP,   [], 2)) ./ rangeG;
end

function b = Dominates(a, b_)
% Dominância de Pareto (minimização): 'a' domina 'b_' se é <= em todos os
% objetivos e < em pelo menos um. Aceita structs (usa .Cost) ou vetores.
    if isstruct(a), a = a.Cost; end; if isstruct(b_), b_ = b_.Cost; end
    b = all(all(a <= b_)) && any(any(a < b_));
end

function pop = DetermineDomination_ref(pop)
% Marca pop(i).IsDominated = true se existe alguma outra solução que o domina.
% Complexidade O(n^2) com acesso a structs — é o principal gargalo de tempo.
    n = numel(pop); for i = 1:n, pop(i).IsDominated = false; end
    for i = 1:n-1
        for j = i+1:n
            if Dominates(pop(i), pop(j)), pop(j).IsDominated = true; end
            if Dominates(pop(j), pop(i)), pop(i).IsDominated = true; end
        end
    end
end

function pop = NonDominatedSort_ref(pop)
% Ordenação rápida por não-dominância (fast non-dominated sort do NSGA-II):
% atribui pop(i).Rank = nº do front (1 = melhor). Base do elitismo do NSGA-II.
    n = numel(pop); domCount = zeros(n,1); dominatedSet = cell(n,1);
    for i = 1:n
        for j = i+1:n, if Dominates(pop(i), pop(j)), dominatedSet{i} = [dominatedSet{i}, j]; domCount(j) = domCount(j) + 1; elseif Dominates(pop(j), pop(i)), dominatedSet{j} = [dominatedSet{j}, i]; domCount(i) = domCount(i) + 1; end; end
    end
    F = find(domCount == 0)'; f_idx = 1;
    while ~isempty(F)
        nextF = []; for i = F, pop(i).Rank = f_idx; for j = dominatedSet{i}, domCount(j) = domCount(j) - 1; if domCount(j) == 0, nextF = [nextF, j]; end; end; end
        F = nextF; f_idx = f_idx + 1;
    end
end

function pop = CalcCrowdingDistance(pop)
% Distância de aglomeração (crowding) por front: mede o "espaço" ao redor de
% cada solução; extremos recebem inf. Preserva diversidade no NSGA-II.
    n = numel(pop); nObj = numel(pop(1).Cost); for i = 1:n, pop(i).CrowdingDistance = 0; end; maxRank = max([pop.Rank]);
    for r = 1:maxRank
        F = find([pop.Rank] == r); l = numel(F); if l == 0, continue; end
        if l <= 2, for i = F, pop(i).CrowdingDistance = inf; end; continue; end
        costs = horzcat(pop(F).Cost)';
        for m = 1:nObj
            [~, idx] = sort(costs(:, m)); F_sorted = F(idx); pop(F_sorted(1)).CrowdingDistance = inf; pop(F_sorted(end)).CrowdingDistance = inf;
            minC = costs(idx(1), m); maxC = costs(idx(end), m); if maxC == minC, continue; end
            for i = 2:l-1, pop(F_sorted(i)).CrowdingDistance = pop(F_sorted(i)).CrowdingDistance + (costs(idx(i+1),m) - costs(idx(i-1),m)) / (maxC - minC); end
        end
    end
end

function popNew = SortAndTruncateNSGA2(pop, nPop)
% Seleção ambiental do NSGA-II: ordena por (Rank crescente, Crowding
% decrescente) e mantém os nPop melhores.
    ranks = [pop.Rank]; cds = [pop.CrowdingDistance]; [~, order] = sortrows([ranks', -cds']); popNew = pop(order(1:nPop));
end

function idx = TournamentSelection(pop)
% Torneio binário: sorteia 2 indivíduos e vence o de menor rank (desempate
% pela maior distância de aglomeração).
    n = numel(pop); i1 = randi(n); i2 = randi(n);
    if pop(i1).Rank < pop(i2).Rank, idx = i1; elseif pop(i2).Rank < pop(i1).Rank, idx = i2; else, if pop(i1).CrowdingDistance > pop(i2).CrowdingDistance, idx = i1; else, idx = i2; end; end
end

function xnew = Mutate(x, pm, VarMin, VarMax)
% Mutação: perturba UMA variável aleatória dentro de uma janela proporcional
% a pm (amplitude da perturbação). Mantém a solução no domínio.
    j = randi(numel(x)); dx = pm*(VarMax-VarMin); lb = max(VarMin, x(j)-dx); ub = min(VarMax, x(j)+dx); xnew = x; xnew(j) = unifrnd(lb, ub);
end

function idx = RouletteWheelSelection(P)
% Seleção por roleta: sorteia um índice com probabilidade proporcional a P.
    r = rand; C = cumsum(P); idx = find(r <= C, 1, 'first');
end

function part = FindGridIndex(part, Grid)
% Mapeia uma solução do repositório para uma célula da grade adaptativa
% (índice de hipercubo no espaço de objetivos). Usado pelo MOPSO.
    nObj = numel(part.Cost); nGrid = numel(Grid(1).LB); part.GridSubIndex = zeros(1, nObj);
    for j = 1:nObj, val = part.Cost(j); if isnan(val), idx = nGrid; else, idx = find(val < Grid(j).UB, 1, 'first'); if isempty(idx), idx = numel(Grid(j).UB); end; end; part.GridSubIndex(j) = idx; end
    part.GridIndex = part.GridSubIndex(1); for j = 2:nObj, part.GridIndex = nGrid*(part.GridIndex-1) + part.GridSubIndex(j); end
end

function leader = SelectLeader(rep, repCount, beta)
% Seleciona o líder global do MOPSO: células MENOS povoadas têm MAIOR
% probabilidade (P ~ exp(-beta*N)), empurrando o enxame para regiões esparsas
% da frente (promove diversidade/cobertura).
    GI = [rep(1:repCount).GridIndex]; OC = unique(GI); N = zeros(size(OC)); for k = 1:numel(OC), N(k) = sum(GI == OC(k)); end
    P = exp(-beta*N); P = P./sum(P); sci = RouletteWheelSelection(P); if isempty(sci), sci = randi(numel(OC)); end
    sc = OC(sci); SCM = find(GI == sc); leader = rep(SCM(randi(numel(SCM))));
end

function [rep, repCount] = DeleteOneRepMember(rep, repCount, gamma)
% Poda do repositório quando excede a capacidade: células MAIS povoadas têm
% MAIOR probabilidade de perder um membro (P ~ exp(gamma*(N-max))), mantendo
% a frente uniforme.
%
% [BUG CORRIGIDO NA V4] Na V3 a assinatura era "function repCount = ...": como
% o MATLAB passa arrays POR VALOR, a remoção "rep(remIdx) = rep(repCount)" era
% feita numa CÓPIA local e PERDIDA ao retornar (só repCount voltava). O efeito
% real era apenas truncar o último elemento — a poda por densidade de grade
% descrita na metodologia (Seção 3.1.4.3) NÃO acontecia. Agora 'rep' também é
% retornado, fazendo a poda funcionar como documentado.
% Para reproduzir EXATAMENTE o comportamento da V3, reverta a assinatura para
% "function repCount = ..." e o chamador no MOPSO para "repCount = DeleteOne...".
    GI = [rep(1:repCount).GridIndex]; OC = unique(GI); N = zeros(size(OC)); for k = 1:numel(OC), N(k) = sum(GI == OC(k)); end
    P = exp(gamma*(N-max(N))); P = P./sum(P); sci = RouletteWheelSelection(P); if isempty(sci), sci = randi(numel(OC)); end
    sc = OC(sci); SCM = find(GI == sc); remIdx = SCM(randi(numel(SCM)));
    rep(remIdx) = rep(repCount);   % move o último para a posição removida (swap-remove)
    repCount = repCount - 1;
end

function Grid = CreateGrid(pop, repCount, nGrid, alpha)
% Cria a grade adaptativa no espaço de objetivos: divide cada eixo em nGrid
% faixas, com folga 'alpha' nas bordas. Recalculada a cada iteração conforme
% a frente evolui (por isso "adaptativa").
    if repCount == 0, Grid = repmat(struct('LB', [-inf 0], 'UB', [0 inf]), 2, 1); return; end
    nObj = numel(pop(1).Cost); c = reshape([pop(1:repCount).Cost], nObj, []); c(isnan(c)) = max(c(~isnan(c)));
    cmin = min(c, [], 2); cmax = max(c, [], 2); dc = cmax-cmin; cmin = cmin-alpha*dc; cmax = cmax+alpha*dc;
    empty.LB = []; empty.UB = []; Grid = repmat(empty, nObj, 1);
    for j = 1:nObj, cj = linspace(cmin(j), cmax(j), nGrid+1); Grid(j).LB = [-inf cj]; Grid(j).UB = [cj +inf]; end
end

function w = PreSelecaoMO(filhos, pai)
% PrÃ©-seleÃ§Ã£o multiobjetivo do multi-child: devolve UM campeÃ£o entre os M
% filhos (equivalente Pareto do argmin da eq. 4):
%   1) se algum filho DOMINA o pai, sorteia entre os que dominam o pai e
%      sÃ£o nÃ£o-dominados entre si (maior pressÃ£o seletiva rumo Ã  frente);
%   2) senÃ£o, sorteia entre os filhos nÃ£o-dominados entre si.
    Mloc = numel(filhos);
    domPai = false(Mloc, 1);
    for m = 1:Mloc, domPai(m) = Dominates(filhos(m).Cost, pai.Cost); end
    if any(domPai)
        idx = find(domPai);
        sub = DetermineDomination(filhos(idx));
        nd = find(~[sub.IsDominated]);
        pick = idx(nd(randi(numel(nd))));
    else
        sub = DetermineDomination(filhos);
        nd = find(~[sub.IsDominated]);
        pick = nd(randi(numel(nd)));
    end
    w = filhos(pick);
    w.IsDominated = false;
end

function r = DistinctIdx(N, excl, k)
% k Ã­ndices distintos em 1..N, todos diferentes entre si e de 'excl'.
    r = zeros(1, k);
    for a = 1:k
        cand = randi(N);
        while cand == excl || any(cand == r(1:a-1)), cand = randi(N); end
        r(a) = cand;
    end
end

function rep = UpdateRepositorio(rep, novos, nRep)
% RepositÃ³rio externo de nÃ£o-dominados (a frente devolvida no final):
% junta, refiltra por dominÃ¢ncia, deduplica em (custo, GWP) e, se exceder a
% capacidade, poda os de MENOR crowding (preserva cobertura da frente).
% Poda determinÃ­stica por crowding â€” anÃ¡loga em intenÃ§Ã£o Ã  poda por grade do
% MOPSO, reutilizando CalcCrowdingDistance.
    rep = [rep; novos];
    if isempty(rep), return; end
    rep = DetermineDomination(rep); rep = rep(~[rep.IsDominated]);
    costs = round(horzcat(rep.Cost)', 6); [~, uIdx] = unique(costs, 'rows'); rep = rep(uIdx);
    if numel(rep) > nRep
        for k = 1:numel(rep), rep(k).Rank = 1; end
        rep = CalcCrowdingDistance(rep);
        [~, ord] = sort([rep.CrowdingDistance], 'descend');
        rep = rep(ord(1:nRep));
    end
end
