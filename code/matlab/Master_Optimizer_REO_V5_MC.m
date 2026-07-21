function Master_Optimizer_REO_V5_MC
% =====================================================================
%  Master_Optimizer_REO  —  VERSÃO 5 (V4 + MO-MCDE + MO-MC-SHADE)
% =====================================================================
%  O QUE MUDOU DA V4 PARA A V5 -----------------------------------------
%  [NOVO] Duas técnicas multi-child de Evolução Diferencial, criadas como
%         versões MULTIOBJETIVO dos algoritmos de Storn:
%           - MO-MCDE     (base: Storn & Price, "Multi-Child DE – A Massively
%                          Parallel Differential Evolution Algorithm",
%                          IEEE CIES/SSCI 2025)
%           - MO-MC-SHADE (base: Storn, "Comparison of MCDE and MC-SHADE for
%                          Massively Parallel Optimization", IEEE CEC 2026)
%         A "multiobjetivação" usa seleção por dominância de Pareto estilo
%         GDE3 (Kukkonen & Lampinen, 2005) + repositório externo podado por
%         crowding. Detalhes nos cabeçalhos de RunMO_MCDE e RunMO_MCSHADE.
%  [NOVO] Modo HEADLESS para execução via "matlab -batch":
%           - env HEADLESS=1  -> pula a GUI (usa defaults / env TECHS);
%           - env KMAX=<n|inf> -> define K_MAX sem editar o arquivo;
%           - env TECHS=MOPSO,NSGA2,... -> roda só as técnicas listadas;
%           - plotagem desligada quando HEADLESS=1.
%  [FIX]  Pool paralelo: parpool agora é limitado ao nº de núcleos e o
%         parfeval ENFILEIRA o excedente (na V4, numSelected > núcleos
%         derrubava o parpool — provável causa das ausências históricas de
%         MOGA no caso base e SI-CDC nas rodadas K: 15 técnicas > 14 núcleos
%         obrigavam a desmarcar uma técnica).
% =====================================================================
%  OBJETIVO
%  Resolver, de forma comparativa, o problema BI-OBJETIVO de seleção de
%  fornecedores para a produção de óxidos de terras raras (REO):
%     minimizar  f1 = custo total   e   f2 = GWP total (kg CO2-eq).
%  Cada material (linha das matrizes) recebe UM fornecedor (coluna).
%  15 técnicas em 3 famílias rodam em paralelo e exportam suas frentes
%  de Pareto para Excel.
%
%  O QUE MUDOU DA V3 PARA A V4 (resumo da revisão) -------------------
%  [NOVO]   Restrição de CONSOLIDAÇÃO (Opção A): no máximo K fornecedores
%           DISTINTOS em toda a cadeia (eq. |{x_1..x_M}| <= K). Quebra a
%           separabilidade do problema e o torna NP-difícil (variante de
%           p-medianas) — é o regime que JUSTIFICA as metaheurísticas.
%           Implementada por:
%             (a) operador de reparo ApplyK() chamado dentro de CostFunction
%                 -> vale automaticamente para TODAS as 15 técnicas
%                    ("reparo Baldwiniano": repara na avaliação);
%             (b) RunMILP reformulada como p-medianas EXATO (intlinprog)
%                 quando K está ativo -> baseline exato correto do problema
%                 restrito (a versão original de atribuição vira o ramo
%                 "sem restrição").
%  [NOVO]   Exportação enriquecida: nº de fornecedores distintos por solução
%           e a própria seleção (para a análise de "padrões de fornecedor"
%           e para AUDITAR que a restrição foi respeitada).
%  [NOVO]   dados.Cn / dados.Gn (matrizes normalizadas por material) usadas
%           pelo reparo para escolher o melhor fornecedor remanescente.
%  [MELHORIA SUGERIDA, comentada no código, NÃO ativada p/ preservar os
%           resultados já validados]:
%             - Limitação de velocidade (Vmax) em MOPSO/SI-CDC (evita
%               explosão de velocidade com w inicial = 1,0);
%             - DetermineDomination vetorizada (10-50x mais rápida; a
%               versão O(n^2) original é o gargalo dos ~30-80 min de execução).
%
%  COMO DESLIGAR A RESTRIÇÃO: definir K_MAX = inf (ou K_MAX >= nVar). Nesse
%  caso o código se comporta exatamente como a V3 (problema separável).
%  COMO VARRER K (experimento do gráfico "custo da consolidação x K"):
%  rodar o script alterando K_MAX para 2, 3, 5, 8, 13 e guardar cada Excel.
%
%  Requisitos: Parallel Computing Toolbox (parpool/parfeval),
%              Optimization Toolbox (intlinprog), Statistics Toolbox (pdist2).
%  Arquivos de entrada (na pasta de trabalho): Custos.csv, GWP.csv.
% =====================================================================

clc; clear; close all;

%% 0. PARÂMETRO DA RESTRIÇÃO DE CONSOLIDAÇÃO  <<< AJUSTE AQUI >>>
% K_MAX = nº máximo de fornecedores DISTINTOS permitido em toda a cadeia.
% Motivo: cada fornecedor extra implica custo de qualificação, contrato,
% auditoria e gestão; consolidar a base é exigência operacional realista e,
% computacionalmente, acopla as decisões (torna o problema NP-difícil).
K_MAX = inf;   % use inf (ou >=13) para DESLIGAR a restrição (= comportamento V3)
% [V5] Sobrescrita por variável de ambiente (para varreduras via -batch):
%      set KMAX=5 && matlab -batch Master_Optimizer_REO_V5_MC
kEnv = getenv('KMAX');
if ~isempty(kEnv)
    kNum = str2double(kEnv);
    if ~isnan(kNum), K_MAX = kNum; end
    fprintf('K_MAX definido via env KMAX = %g\n', K_MAX);
end

%% 1. Carregamento dos Dados
% Lê as matrizes de custo e GWP (materiais x fornecedores) dos CSVs e monta
% a struct 'dados'. Em caso de erro de leitura, aborta com mensagem clara.
fprintf('Lendo arquivos CSV...\n');
try
    dados = Dados_Processo_REO();
    if dados.nVar == 0
        error('Nao existem materiais validos para otimizacao.');
    end
catch ME
    errordlg(sprintf('Erro critico na leitura:\n%s', ME.message), 'Erro');
    error('Execucao interrompida: %s', ME.message);
end

% Injeta o parâmetro de consolidação na struct para que CHEGUE A TODOS os
% workers paralelos (dados é copiado para cada parfeval). Se K_MAX >= nVar,
% a restrição é inócua (não há como usar mais fornecedores que materiais).
dados.K = K_MAX;
if dados.K < dados.nVar
    fprintf('Restricao de consolidacao ATIVA: no maximo K = %d fornecedores distintos (de %d materiais).\n', dados.K, dados.nVar);
else
    fprintf('Restricao de consolidacao DESLIGADA (K = %g >= %d materiais): problema separavel.\n', dados.K, dados.nVar);
end

%% 2. Interface Grafica e Parametros
% GUI com checkboxes para escolher as técnicas e os parâmetros de cada uma.
% Tem um timer de 15 s que inicia o "teste completo" (todas as técnicas) se
% o usuário não interagir.
[flags, params] = InterfaceSelecaoHeuristicas();

% Conta quantas técnicas foram selecionadas (define o dimensionamento do pool).
numSelected = flags.MOPSO + flags.NSGA2 + flags.MOGA + ...
              flags.MOEAD + flags.LEXICOGRAFICO + flags.MOEAD_AV + ...
              flags.MOEAD_PS + flags.MOEAD_STN + flags.SI_CDC + ...
              flags.AHP + flags.PROMETHEE + flags.ELECTRE + ...
              flags.MILP + flags.BNB + flags.LAGRANGIANA + ...
              flags.MO_MCDE + flags.MO_MCSHADE;

if numSelected == 0
    fprintf('Nenhuma tecnica selecionada. Encerrando...\n');
    return;
end

%% 3. Configuracao do Pool Paralelo Dinamico
% Distribui os núcleos disponíveis entre as técnicas: cada técnica recebe
% 'runsPerTech' execuções independentes (réplicas) para captar a
% estocasticidade. totalWorkers = técnicas x réplicas.
hwCores = feature('numcores');
runsPerTech = max(1, floor(hwCores / numSelected));
totalWorkers = numSelected * runsPerTech;
totalTasks = totalWorkers;

fprintf('\nIniciando pool com %d workers (%d execucoes por tecnica baseando-se em %d nucleos detectados)...\n', totalWorkers, runsPerTech, hwCores);

% [V5-FIX] O pool é limitado ao nº de núcleos físicos; o parfeval enfileira
% as tarefas excedentes automaticamente (na V4, parpool(totalWorkers) com
% totalWorkers > núcleos ERRAVA, forçando a desmarcar técnicas na mão).
poolSize = min(totalWorkers, hwCores);
pool = gcp('nocreate');
if isempty(pool) || pool.NumWorkers ~= poolSize
    delete(gcp('nocreate'));
    pool = parpool(poolSize);
end

% Fila de mensagens dos workers -> imprime progresso com timestamp [HH:MM:SS].
% Motivo: parfeval não imprime direto no console; a DataQueue centraliza os logs.
dq = parallel.pool.DataQueue;
afterEach(dq, @(msg) fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), msg));

%% 4. Despacho das Tarefas Paralelas (parfeval)
% Cada técnica selecionada é submetida como tarefa assíncrona (parfeval).
% 'taskInfo' guarda nome e índice de cor para a plotagem/legenda.
f = parallel.FevalFuture.empty;
taskInfo = struct('Name', {}, 'ColorIdx', {});

% --- Metaheuristicas (busca populacional baseada em dominância de Pareto) ---
if flags.MOPSO
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunMOPSO, 2, dados, params.MOPSO, dq); taskInfo(end+1).Name = sprintf('MOPSO (Run %d)', r); taskInfo(end).ColorIdx = 1; end
end
if flags.NSGA2
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunNSGA2, 2, dados, params.NSGA2, dq); taskInfo(end+1).Name = sprintf('NSGA-II (Run %d)', r); taskInfo(end).ColorIdx = 2; end
end
if flags.MOGA
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunMOGA, 2, dados, params.MOGA, dq); taskInfo(end+1).Name = sprintf('MOGA (Run %d)', r); taskInfo(end).ColorIdx = 3; end
end
if flags.MOEAD
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunMOEAD, 2, dados, params.MOEAD, dq); taskInfo(end+1).Name = sprintf('MOEA/D (Run %d)', r); taskInfo(end).ColorIdx = 4; end
end
if flags.MOEAD_AV
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunMOEAD_AV, 2, dados, params.MOEAD_AV, dq); taskInfo(end+1).Name = sprintf('MOEA/D-AV (Run %d)', r); taskInfo(end).ColorIdx = 5; end
end
if flags.MOEAD_PS
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunMOEAD_PS, 2, dados, params.MOEAD_PS, dq); taskInfo(end+1).Name = sprintf('MOEA/D-PS (Run %d)', r); taskInfo(end).ColorIdx = 6; end
end
if flags.MOEAD_STN
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunMOEAD_STN, 2, dados, params.MOEAD_STN, dq); taskInfo(end+1).Name = sprintf('MOEA/D-STN (Run %d)', r); taskInfo(end).ColorIdx = 7; end
end
if flags.SI_CDC
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunSI_CDC, 2, dados, params.SI_CDC, dq); taskInfo(end+1).Name = sprintf('SI-CDC (Run %d)', r); taskInfo(end).ColorIdx = 8; end
end

% --- [V5] Multi-child DE multiobjetivo (novas técnicas 16 e 17) ---
if flags.MO_MCDE
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunMO_MCDE, 2, dados, params.MO_MCDE, dq); taskInfo(end+1).Name = sprintf('MO-MCDE (Run %d)', r); taskInfo(end).ColorIdx = 16; end
end
if flags.MO_MCSHADE
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunMO_MCSHADE, 2, dados, params.MO_MCSHADE, dq); taskInfo(end+1).Name = sprintf('MO-MC-SHADE (Run %d)', r); taskInfo(end).ColorIdx = 17; end
end

% --- MCDM (multicritério: constroem solução por material, varrendo pesos) ---
if flags.AHP
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunAHP, 2, dados, params.AHP, dq); taskInfo(end+1).Name = sprintf('AHP (Run %d)', r); taskInfo(end).ColorIdx = 9; end
end
if flags.PROMETHEE
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunPROMETHEE, 2, dados, params.PROMETHEE, dq); taskInfo(end+1).Name = sprintf('PROMETHEE (Run %d)', r); taskInfo(end).ColorIdx = 10; end
end
if flags.ELECTRE
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunELECTRE, 2, dados, params.ELECTRE, dq); taskInfo(end+1).Name = sprintf('ELECTRE (Run %d)', r); taskInfo(end).ColorIdx = 11; end
end

% --- Sequenciais / Deterministicos / Exatos ---
if flags.LEXICOGRAFICO
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunLexicografico, 2, dados, params.LEXICOGRAFICO, dq); taskInfo(end+1).Name = sprintf('Lexicografico (Run %d)', r); taskInfo(end).ColorIdx = 12; end
end
if flags.MILP
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunMILP, 2, dados, params.MILP, dq); taskInfo(end+1).Name = sprintf('MILP (Run %d)', r); taskInfo(end).ColorIdx = 13; end
end
if flags.BNB
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunBNB, 2, dados, params.BNB, dq); taskInfo(end+1).Name = sprintf('Branch and Bound (Run %d)', r); taskInfo(end).ColorIdx = 14; end
end
if flags.LAGRANGIANA
    for r = 1:runsPerTech, f(end+1) = parfeval(pool, @RunLagrangiana, 2, dados, params.LAGRANGIANA, dq); taskInfo(end+1).Name = sprintf('Relax. Lagrangiana (Run %d)', r); taskInfo(end).ColorIdx = 15; end
end

%% 5. Coleta de Resultados, Plotagem e Exportacao Excel
% Conforme cada tarefa termina (fetchNext), seus pontos (custo, GWP) são
% plotados e acumulados para o Excel. NOVO: para cada solução também se
% registra o nº de fornecedores distintos e a seleção, permitindo auditar
% a restrição e analisar os padrões de escolha.
cores = lines(17);
% [V5] Plotagem desligada em modo headless (matlab -batch não tem display útil).
doPlot = isempty(getenv('HEADLESS'));
if doPlot
    figure('Name', 'Comparacao das Fronteiras de Pareto', 'Position', [100 100 900 650]);
    hold on; grid on;
    xlabel('Custo Total ($)'); ylabel('Impacto Ambiental (CO2-eq)');
    title('Comparativo de Fronteiras de Pareto - Multi-Metodologias');
end

fprintf('\nAguardando algoritmos finalizarem...\n\n');

techNames = {'MOPSO', 'NSGA-II', 'MOGA', 'MOEA/D', ...
             'MOEA/D-AV', 'MOEA/D-PS', 'MOEA/D-STN', 'SI-CDC', ...
             'AHP', 'PROMETHEE', 'ELECTRE', ...
             'Lexicografico', 'MILP (intlinprog)', 'Branch and Bound', 'Relax. Lagrangiana', ...
             'MO-MCDE', 'MO-MC-SHADE'};
legendAdded = false(1, 17);

% Vetores-coluna que vão virar a tabela do Excel.
xls_Algoritmo = {};
xls_Custo = [];
xls_GWP = [];
xls_Tempo = [];
xls_NumForn = [];      % NOVO: nº de fornecedores distintos por solução
xls_Selecao = {};      % NOVO: vetor de fornecedores escolhidos (string)
maxDistintosGlobal = 0; % NOVO: para o relatório de viabilidade da restrição

for i = 1:totalTasks
    % Bloqueia até a PRÓXIMA tarefa concluir (ordem de término, não de envio).
    [completedIdx, paretoResult, tempo] = fetchNext(f);
    nomeAlg = taskInfo(completedIdx).Name;
    corIdx = taskInfo(completedIdx).ColorIdx;

    fprintf('>>> [%s] FINALIZADO em %.2f segundos! <<<\n', nomeAlg, tempo);

    if ~isempty(paretoResult)
        if isstruct(paretoResult) && isfield(paretoResult, 'Cost') && ~isempty(paretoResult(1).Cost)
            custos = horzcat(paretoResult.Cost)';   % nSol x 2 (col 1 = custo, col 2 = GWP)

            % Alimenta as colunas do Excel, solução a solução.
            for sol = 1:size(custos, 1)
                xls_Algoritmo{end+1, 1} = nomeAlg;
                xls_Custo(end+1, 1) = custos(sol, 1);
                xls_GWP(end+1, 1) = custos(sol, 2);
                xls_Tempo(end+1, 1) = tempo;

                % NOVO: reconstrói a seleção de fornecedores da solução,
                % aplicando a MESMA discretização+reparo da avaliação, para
                % que a seleção reportada seja coerente com o custo/GWP
                % (que já foram calculados sobre a versão reparada).
                posSol = [];
                if isfield(paretoResult(sol), 'Position') && ~isempty(paretoResult(sol).Position)
                    posSol = max(1, min(round(paretoResult(sol).Position), dados.nFornecedores));
                    if isfield(dados,'K') && ~isempty(dados.K) && dados.K < dados.nVar
                        posSol = ApplyK(posSol, dados);
                    end
                end
                if isempty(posSol)
                    xls_NumForn(end+1, 1) = NaN; xls_Selecao{end+1, 1} = '';
                else
                    nd = numel(unique(posSol));
                    xls_NumForn(end+1, 1) = nd;
                    xls_Selecao{end+1, 1} = mat2str(posSol(:)');
                    maxDistintosGlobal = max(maxDistintosGlobal, nd);
                end
            end

            % Plotagem: adiciona à legenda apenas na 1ª aparição da técnica.
            if doPlot
                if ~legendAdded(corIdx)
                    plot(custos(:,1), custos(:,2), 'o', 'MarkerFaceColor', cores(corIdx,:), ...
                        'MarkerEdgeColor', 'k', 'DisplayName', techNames{corIdx}, ...
                        'MarkerSize', 8);
                    legendAdded(corIdx) = true;
                else
                    plot(custos(:,1), custos(:,2), 'o', 'MarkerFaceColor', cores(corIdx,:), ...
                        'MarkerEdgeColor', 'k', 'HandleVisibility', 'off', ...
                        'MarkerSize', 8);
                end
            end
        end
    end
end

if doPlot
    legend('Location', 'best');
    hold off;
end
fprintf('\nTodas as otimizacoes foram concluidas.\n');

% NOVO: relatório de viabilidade da restrição (sanity check). Se a restrição
% está ativa, o nº máximo de fornecedores distintos observado DEVE ser <= K.
if dados.K < dados.nVar
    if maxDistintosGlobal <= dados.K
        fprintf('Verificacao da restricao OK: max de fornecedores distintos = %d (<= K = %d).\n', maxDistintosGlobal, dados.K);
    else
        warning('Restricao VIOLADA em alguma solucao: max distintos = %d > K = %d. Verifique ApplyK.', maxDistintosGlobal, dados.K);
    end
end

% --- Gravação do Excel (nome do arquivo embute o K para não sobrescrever
%     execuções de Ks diferentes — útil no experimento de varredura de K). ---
if ~isempty(xls_Custo)
    T_Export = table(xls_Algoritmo, xls_Custo, xls_GWP, xls_Tempo, xls_NumForn, xls_Selecao, ...
        'VariableNames', {'Algoritmo', 'Custo_Total', 'Impacto_Ambiental', 'Tempo_Total_Segundos', 'Num_Fornecedores_Distintos', 'Selecao_Fornecedores'});

    if dados.K < dados.nVar
        nomeArquivo = sprintf('Fronteira_Pareto_Global_K%d.xlsx', dados.K);
    else
        nomeArquivo = 'Fronteira_Pareto_Global.xlsx';
    end
    writetable(T_Export, nomeArquivo);
    fprintf('\n>>> DADOS EXPORTADOS COM SUCESSO PARA "%s" <<<\n', nomeArquivo);
end

end

% =====================================================================
%   INTERFACE GRAFICA & PARAMETROS (AJUSTADOS PARA TESTE COMPLETO)
% =====================================================================
% Define os hiperparâmetros de cada técnica e exibe checkboxes para
% seleção. OBS (revisão): estes valores são os EFETIVAMENTE executados —
% ao escrever a dissertação, alinhe o texto a eles (MaxIt=2000, nPop=500,
% c1=c2=1,5, mu=0,5 para o MOPSO; etc.).
function [flags, params] = InterfaceSelecaoHeuristicas()
    % --- Parâmetros das metaheurísticas ---
    % MaxIt: nº de iterações; nPop: tamanho da população; nRep: tamanho do
    % repositório externo (arquivo de não-dominados); w/wdamp: inércia e seu
    % decaimento; c1/c2: coef. cognitivo/social; mu: controle da mutação.
    params.MOPSO = struct('MaxIt', 2000, 'nPop', 500, 'nRep', 500, 'w', 1.0, 'wdamp', 0.99, 'c1', 1.5, 'c2', 1.5, 'mu', 0.5);
    params.NSGA2 = struct('MaxIt', 2000, 'nPop', 500, 'pCrossover', 0.8, 'pMutation', 0.2);
    params.MOGA  = struct('MaxIt', 2000, 'nPop', 500, 'pCrossover', 0.8, 'pMutation', 0.2);
    params.MOEAD = struct('MaxIt', 2000, 'nPop', 500, 'T', 10);                 % T: tamanho da vizinhança
    params.LEXICOGRAFICO = struct('MaxIt', 2000, 'nPop', 500, 'Tolerancia', 1.05);
    params.MOEAD_AV = struct('MaxIt', 2000, 'nPopInit', 500, 'T', 10);          % AV: população adaptativa
    params.MOEAD_PS = struct('MaxIt', 2000, 'nPop', 500, 'T', 10, 'nRatio', 0.1);% PS: atualização parcial
    params.MOEAD_STN = struct('MaxIt', 2000, 'nPop', 500, 'T', 10);             % STN: rastreia trajetória
    params.SI_CDC = struct('MaxIt', 2000, 'nPop', 500, 'nRep', 50, 'w', 0.8, 'c1', 1.5, 'c2', 1.5);

    % --- MCDM e exatos: nPop = nº de pesos (cenários) varridos em [0,1] ---
    params.AHP = struct('nPop', 500);
    params.PROMETHEE = struct('nPop', 500);
    params.ELECTRE = struct('nPop', 500);
    params.MILP = struct('nPop', 500);   % p/ K ativo, RunMILP limita internamente o nº de pesos (frente restrita tem poucos pontos suportados)
    params.BNB = struct('nPop', 500);
    params.LAGRANGIANA = struct('nPop', 500, 'MaxIt', 2000, 'step', 0.05);

    % --- [V5] Novas técnicas multi-child DE (multiobjetivo) ---
    % Orçamento de avaliações IGUAL ao das demais metaheurísticas:
    %   demais: MaxIt=2000 x nPop=500        = 1.000.000 avaliações
    %   V5:     MaxIt= 500 x nPop=500 x M=4  = 1.000.000 avaliações
    % (o "multi-child" troca gerações sequenciais por filhos paralelizáveis —
    %  é exatamente a proposta do MCDE; ver cabeçalho de RunMO_MCDE.)
    % F_L/F_U (dithering) e Cr do paper do CEC 2026; H = N/3 do MC-SHADE.
    params.MO_MCDE    = struct('MaxIt', 500, 'nPop', 500, 'M', 4, ...
                               'F_L', 0.3, 'F_U', 1.6, 'Cr', 0.85, 'nRep', 500);
    params.MO_MCSHADE = struct('MaxIt', 500, 'nPop', 500, 'M', 4, ...
                               'H', round(500/3), 'c', 0.1, 'nRep', 500);

    % Por padrão, TODAS as técnicas ligadas (modo "teste completo").
    flags.MOPSO = 1; flags.NSGA2 = 1; flags.MOGA = 1; flags.MOEAD = 1;
    flags.LEXICOGRAFICO = 1; flags.MOEAD_AV = 1; flags.MOEAD_PS = 1; flags.MOEAD_STN = 1;
    flags.SI_CDC = 1; flags.AHP = 1; flags.PROMETHEE = 1; flags.ELECTRE = 1;
    flags.MILP = 1; flags.BNB = 1; flags.LAGRANGIANA = 1;
    flags.MO_MCDE = 1; flags.MO_MCSHADE = 1;

    % --- [V5] Modo HEADLESS: pula a GUI por completo (p/ matlab -batch). ---
    % env TECHS (opcional): lista separada por vírgula das técnicas a rodar,
    % com os nomes dos campos de 'flags' (ex.: TECHS=MO_MCDE,MO_MCSHADE).
    if ~isempty(getenv('HEADLESS'))
        tsel = getenv('TECHS');
        if ~isempty(tsel)
            fn = fieldnames(flags);
            for kf = 1:numel(fn), flags.(fn{kf}) = 0; end
            partes = strsplit(tsel, ',');
            for kp = 1:numel(partes)
                key = strtrim(partes{kp});
                if isfield(flags, key)
                    flags.(key) = 1;
                else
                    warning('TECHS: tecnica desconhecida "%s" ignorada.', key);
                end
            end
        end
        ligadas = ''; fn = fieldnames(flags);
        for kf = 1:numel(fn), if flags.(fn{kf}), ligadas = [ligadas ' ' fn{kf}]; end; end %#ok<AGROW>
        fprintf('Modo HEADLESS: GUI ignorada. Tecnicas ativas:%s\n', ligadas);
        return;
    end

    fig = figure('Name', 'Configuracao de Otimizacao (Modo Teste Completo)', 'Position', [400 100 450 755], ...
                 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off', 'Resize', 'off');

    uicontrol('Style', 'text', 'String', 'Metaheuristicas:', 'Position', [30 715 200 20], ...
              'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 9, 'ForegroundColor', [0 0.3 0.6]);

    chkMOPSO     = uicontrol('Style', 'checkbox', 'String', 'MOPSO', 'Position', [40 690 300 20], 'Value', 1);
    chkNSGA2     = uicontrol('Style', 'checkbox', 'String', 'NSGA-II', 'Position', [40 665 300 20], 'Value', 1);
    chkMOGA      = uicontrol('Style', 'checkbox', 'String', 'MOGA', 'Position', [40 640 300 20], 'Value', 1);
    chkMOEAD     = uicontrol('Style', 'checkbox', 'String', 'MOEA/D', 'Position', [40 615 300 20], 'Value', 1);
    chkMOEAD_AV  = uicontrol('Style', 'checkbox', 'String', 'MOEA/D-AV', 'Position', [40 590 300 20], 'Value', 1);
    chkMOEAD_PS  = uicontrol('Style', 'checkbox', 'String', 'MOEA/D-PS', 'Position', [40 565 300 20], 'Value', 1);
    chkMOEAD_STN = uicontrol('Style', 'checkbox', 'String', 'MOEA/D-STN', 'Position', [40 540 300 20], 'Value', 1);
    chkSI_CDC    = uicontrol('Style', 'checkbox', 'String', 'SI-CDC', 'Position', [40 515 300 20], 'Value', 1);

    uicontrol('Style', 'text', 'String', 'Metodos MCDM:', 'Position', [30 480 300 20], ...
              'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 9, 'ForegroundColor', [0 0.3 0.6]);

    chkAHP       = uicontrol('Style', 'checkbox', 'String', 'AHP/ANP', 'Position', [40 455 300 20], 'Value', 1);
    chkPROMETHEE = uicontrol('Style', 'checkbox', 'String', 'PROMETHEE', 'Position', [40 430 300 20], 'Value', 1);
    chkELECTRE   = uicontrol('Style', 'checkbox', 'String', 'ELECTRE', 'Position', [40 405 300 20], 'Value', 1);

    uicontrol('Style', 'text', 'String', 'Metodos Exatos / Deterministicos:', 'Position', [30 370 300 20], ...
              'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 9, 'ForegroundColor', [0 0.3 0.6]);

    chkLEXI      = uicontrol('Style', 'checkbox', 'String', 'Lexicografico (Artigo REO)', 'Position', [40 345 300 20], 'Value', 1);
    chkMILP      = uicontrol('Style', 'checkbox', 'String', 'MILP (Programacao Linear Inteira Mista)', 'Position', [40 320 300 20], 'Value', 1);
    chkBNB       = uicontrol('Style', 'checkbox', 'String', 'Branch and Bound', 'Position', [40 295 300 20], 'Value', 1);
    chkLAGR      = uicontrol('Style', 'checkbox', 'String', 'Relaxacao Lagrangiana', 'Position', [40 270 300 20], 'Value', 1);

    uicontrol('Style', 'text', 'String', 'Multi-Child DE (novas V5):', 'Position', [30 235 300 20], ...
              'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 9, 'ForegroundColor', [0 0.3 0.6]);

    chkMO_MCDE    = uicontrol('Style', 'checkbox', 'String', 'MO-MCDE', 'Position', [40 210 140 20], 'Value', 1);
    chkMO_MCSHADE = uicontrol('Style', 'checkbox', 'String', 'MO-MC-SHADE', 'Position', [190 210 160 20], 'Value', 1);

    txtTimer = uicontrol('Style', 'text', 'String', 'Iniciando padrao em 15s...', ...
                         'Position', [30 140 340 20], 'ForegroundColor', 'r', 'FontSize', 10);

    btnConfig = uicontrol('Style', 'pushbutton', 'String', 'Configurar Parametros', ...
                           'Position', [40 80 150 40], 'Callback', 'set(gcbf, ''UserData'', 2); uiresume(gcbf);'); %#ok<NASGU>

    btnStart  = uicontrol('Style', 'pushbutton', 'String', 'Executar Teste Agora', ...
                           'Position', [210 80 150 40], 'FontWeight', 'bold', ...
                           'Callback', 'set(gcbf, ''UserData'', 1); uiresume(gcbf);'); %#ok<NASGU>

    % Timer de 15 s: se ninguém clicar, dispara o teste completo automaticamente.
    set(fig, 'UserData', 0);
    t = 15;
    while t > 0 && ishandle(fig) && get(fig, 'UserData') == 0
        set(txtTimer, 'String', sprintf('Iniciando Teste com parametros basicos em %ds...', t));
        pause(1); t = t - 1;
    end
    if ~ishandle(fig), return; end

    acao = get(fig, 'UserData'); if acao == 0, acao = 1; end %#ok<NASGU>

    % Lê o estado final dos checkboxes para montar 'flags'.
    flags.MOPSO = get(chkMOPSO, 'Value'); flags.NSGA2 = get(chkNSGA2, 'Value'); flags.MOGA  = get(chkMOGA, 'Value');
    flags.MOEAD = get(chkMOEAD, 'Value'); flags.LEXICOGRAFICO = get(chkLEXI, 'Value');
    flags.MOEAD_AV = get(chkMOEAD_AV, 'Value'); flags.MOEAD_PS = get(chkMOEAD_PS, 'Value'); flags.MOEAD_STN = get(chkMOEAD_STN, 'Value');
    flags.SI_CDC = get(chkSI_CDC, 'Value'); flags.AHP = get(chkAHP, 'Value'); flags.PROMETHEE = get(chkPROMETHEE, 'Value');
    flags.ELECTRE = get(chkELECTRE, 'Value'); flags.MILP = get(chkMILP, 'Value'); flags.BNB = get(chkBNB, 'Value');
    flags.LAGRANGIANA = get(chkLAGR, 'Value');
    flags.MO_MCDE = get(chkMO_MCDE, 'Value'); flags.MO_MCSHADE = get(chkMO_MCSHADE, 'Value');
    close(fig);
end

% =====================================================================
%   MÓDULOS 1-8: METAHEURÍSTICAS
%   (todas avaliam via CostFunction, logo TODAS já respeitam a restrição K
%    pelo reparo Baldwiniano embutido em CostFunction.)
% =====================================================================

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

function [repValido, tempo_execucao] = RunMOGA(dados, params, dq)
% MOGA — algoritmo genético multiobjetivo com aptidão baseada no rank de
% dominância (1/rank) e seleção por roleta; elitismo dos 10% melhores.
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar];
    empty.Position = []; empty.Cost = []; empty.Rank = []; empty.Fitness = []; empty.IsDominated = false; pop = repmat(empty, params.nPop, 1);
    for i = 1:params.nPop, pop(i).Position = round(unifrnd(VarMin, VarMax, VarSize)); pop(i).Cost = CostFunction(pop(i).Position, dados); end
    for it = 1:params.MaxIt
        % Rank por contagem de dominadores; aptidão = 1/rank; probabilidade ~ aptidão.
        pop = DetermineDominationCount(pop); fits = 1 ./ [pop.Rank]; probs = fits / sum(fits); popc = repmat(empty, params.nPop, 1);
        for i = 1:2:params.nPop
            p1 = RouletteWheelSelection(probs); p2 = RouletteWheelSelection(probs);
            if rand < params.pCrossover, mask = rand(VarSize) > 0.5; popc(i).Position = pop(p1).Position; popc(i).Position(mask) = pop(p2).Position(mask); popc(i+1).Position = pop(p2).Position; popc(i+1).Position(mask) = pop(p1).Position(mask); else, popc(i).Position = pop(p1).Position; popc(i+1).Position = pop(p2).Position; end
            if rand < params.pMutation, popc(i).Position = Mutate(popc(i).Position, 0.1, VarMin, VarMax); end; if rand < params.pMutation, popc(i+1).Position = Mutate(popc(i+1).Position, 0.1, VarMin, VarMax); end
            popc(i).Position = round(popc(i).Position); popc(i+1).Position = round(popc(i+1).Position); popc(i).Cost = CostFunction(popc(i).Position, dados); popc(i+1).Cost = CostFunction(popc(i+1).Position, dados); popc(i).IsDominated = false; popc(i+1).IsDominated = false;
        end
        % Elitismo: preserva os melhores não-dominados na nova geração.
        nd_indices = find([pop.Rank] == 1); qtdElitismo = min(length(nd_indices), round(params.nPop * 0.1)); popc(1:qtdElitismo) = pop(nd_indices(1:qtdElitismo)); pop = popc;

        if mod(it, 100) == 0 || it == params.MaxIt
            send(dq, sprintf('MOGA       : Iter %4d / %d concluida', it, params.MaxIt));
        end
    end
    pop = DetermineDominationCount(pop); tempo_execucao = toc(t0); repValido = pop([pop.Rank] == 1);
end

function [repValido, tempo_execucao] = RunMOEAD(dados, params, dq)
% MOEA/D — decomposição: N subproblemas escalares (pesos W) resolvidos
% cooperativamente usando a VIZINHANÇA de pesos (T). z = ponto de referência
% ideal; abordagem de Tchebycheff (max W.*|f - z|). EP = arquivo externo.
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar]; nObj = 2; N = params.nPop; W = zeros(N, nObj);
    for i = 1:N, W(i, 1) = (i-1)/(N-1); W(i, 2) = 1 - W(i, 1); end; T = min(params.T, N); sp = struct('V', [], 'Neighbors', []); sp = repmat(sp, N, 1);
    % Vizinhança: os T pesos mais próximos (em distância euclidiana no espaço de pesos).
    for i = 1:N, D = pdist2(W(i,:), W); [~, idx] = sort(D); sp(i).Neighbors = idx(1:T); end; empty_ind.Position = []; empty_ind.Cost = []; empty_ind.IsDominated = false; pop = repmat(empty_ind, N, 1); z = inf(nObj, 1);
    for i = 1:N, pop(i).Position = round(unifrnd(VarMin, VarMax, VarSize)); pop(i).Cost = CostFunction(pop(i).Position, dados); z = min(z, pop(i).Cost); end; EP = pop;
    for it = 1:params.MaxIt
        for i = 1:N
            % Recombina dois vizinhos, muta, atualiza z e os vizinhos que melhorarem.
            neighbors = sp(i).Neighbors; p1 = neighbors(randi(T)); p2 = neighbors(randi(T)); y = empty_ind; y.Position = pop(p1).Position; mask = rand(VarSize) > 0.5; y.Position(mask) = pop(p2).Position(mask); y.Position = round(Mutate(y.Position, 0.1, VarMin, VarMax));
            y.Cost = CostFunction(y.Position, dados); y.IsDominated = false; z = min(z, y.Cost);
            for j = 1:T, k = neighbors(j); g_old = max(W(k,:)' .* abs(pop(k).Cost - z)); g_new = max(W(k,:)' .* abs(y.Cost - z)); if g_new <= g_old, pop(k) = y; end; end
            EP = [EP; y]; %#ok<AGROW>
        end
        % Atualiza o arquivo externo (não-dominados únicos).
        EP = DetermineDomination(EP); EP = EP(~[EP.IsDominated]); if ~isempty(EP), costs = round(horzcat(EP.Cost)', 4); [~, uniqueIdx] = unique(costs, 'rows'); EP = EP(uniqueIdx); end

        if mod(it, 100) == 0 || it == params.MaxIt
            send(dq, sprintf('MOEA/D     : Iter %4d / %d concluida', it, params.MaxIt));
        end
    end
    tempo_execucao = toc(t0); repValido = EP;
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

function [repValido, tempo_execucao] = RunMOEAD_PS(dados, params, dq)
% MOEA/D-PS — variante com SELEÇÃO PARCIAL (Partial Update): a cada iteração
% só uma fração dos subproblemas é atualizada (+ os extremos 1 e N),
% reduzindo custo computacional por iteração.
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar]; nObj = 2; N = params.nPop; n_update = max(1, round(N * params.nRatio)); W = zeros(N, nObj);
    for i = 1:N, W(i, 1) = (i-1)/(N-1); W(i, 2) = 1 - W(i, 1); end; T = min(params.T, N); sp = struct('V', [], 'Neighbors', []); sp = repmat(sp, N, 1);
    for i = 1:N, D = pdist2(W(i,:), W); [~, idx] = sort(D); sp(i).Neighbors = idx(1:T); end; empty_ind.Position = []; empty_ind.Cost = []; empty_ind.IsDominated = false; pop = repmat(empty_ind, N, 1); z = inf(nObj, 1);
    for i = 1:N, pop(i).Position = round(unifrnd(VarMin, VarMax, VarSize)); pop(i).Cost = CostFunction(pop(i).Position, dados); z = min(z, pop(i).Cost); end; EP = pop;
    for it = 1:params.MaxIt
        % Seleciona quais subproblemas atualizar nesta iteração.
        if N > 2, cand = 2:(N-1); if length(cand) > n_update, selecionados = [1, N, cand(randperm(length(cand), n_update))]; else, selecionados = 1:N; end; else, selecionados = 1:N; end
        for idx = 1:length(selecionados)
            i = selecionados(idx); neighbors = sp(i).Neighbors; p1 = neighbors(randi(T)); p2 = neighbors(randi(T)); y = empty_ind; y.Position = pop(p1).Position; mask = rand(VarSize) > 0.5; y.Position(mask) = pop(p2).Position(mask); y.Position = round(Mutate(y.Position, 0.1, VarMin, VarMax));
            y.Cost = CostFunction(y.Position, dados); y.IsDominated = false; z = min(z, y.Cost);
            for j = 1:T, k = neighbors(j); g_old = max(W(k,:)' .* abs(pop(k).Cost - z)); g_new = max(W(k,:)' .* abs(y.Cost - z)); if g_new <= g_old, pop(k) = y; end; end
            EP = [EP; y]; %#ok<AGROW>
        end
        EP = DetermineDomination(EP); EP = EP(~[EP.IsDominated]); if ~isempty(EP), costs = round(horzcat(EP.Cost)', 4); [~, uniqueIdx] = unique(costs, 'rows'); EP = EP(uniqueIdx); end

        if mod(it, 100) == 0 || it == params.MaxIt
            send(dq, sprintf('MOEA/D-PS  : Iter %4d / %d concluida', it, params.MaxIt));
        end
    end
    tempo_execucao = toc(t0); repValido = EP;
end

function [repValido, tempo_execucao] = RunMOEAD_STN(dados, params, dq)
% MOEA/D-STN — variante que rastreia a TRAJETÓRIA (Search Trajectory Network)
% do subproblema central, registrando como a solução central evolui (útil
% para diagnóstico de convergência). A otimização em si segue o MOEA/D base.
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar]; nObj = 2; N = params.nPop; W = zeros(N, nObj);
    for i = 1:N, W(i, 1) = (i-1)/(N-1); W(i, 2) = 1 - W(i, 1); end; T = min(params.T, N); sp = struct('V', [], 'Neighbors', []); sp = repmat(sp, N, 1);
    for i = 1:N, D = pdist2(W(i,:), W); [~, idx] = sort(D); sp(i).Neighbors = idx(1:T); end; empty_ind.Position = []; empty_ind.Cost = []; empty_ind.IsDominated = false; pop = repmat(empty_ind, N, 1); z = inf(nObj, 1);
    for i = 1:N, pop(i).Position = round(unifrnd(VarMin, VarMax, VarSize)); pop(i).Cost = CostFunction(pop(i).Position, dados); z = min(z, pop(i).Cost); end; EP = pop; center_idx = max(1, round(N/2)); stn_trajectory_nodes = {pop(center_idx).Position};
    for it = 1:params.MaxIt
        for i = 1:N
            neighbors = sp(i).Neighbors; p1 = neighbors(randi(T)); p2 = neighbors(randi(T)); y = empty_ind; y.Position = pop(p1).Position; mask = rand(VarSize) > 0.5; y.Position(mask) = pop(p2).Position(mask); y.Position = round(Mutate(y.Position, 0.1, VarMin, VarMax));
            y.Cost = CostFunction(y.Position, dados); y.IsDominated = false; z = min(z, y.Cost);
            for j = 1:T, k = neighbors(j); g_old = max(W(k,:)' .* abs(pop(k).Cost - z)); g_new = max(W(k,:)' .* abs(y.Cost - z)); if g_new <= g_old, pop(k) = y; if k == center_idx && ~isequal(stn_trajectory_nodes{end}, y.Position), stn_trajectory_nodes{end+1} = y.Position; end; end; end %#ok<AGROW>
            EP = [EP; y]; %#ok<AGROW>
        end
        EP = DetermineDomination(EP); EP = EP(~[EP.IsDominated]); if ~isempty(EP), costs = round(horzcat(EP.Cost)', 4); [~, uniqueIdx] = unique(costs, 'rows'); EP = EP(uniqueIdx); end

        if mod(it, 100) == 0 || it == params.MaxIt
            send(dq, sprintf('MOEA/D-STN : Iter %4d / %d concluida', it, params.MaxIt));
        end
    end
    tempo_execucao = toc(t0); repValido = EP;
end

function [repValido, tempo_execucao] = RunSI_CDC(dados, params, dq)
% SI-CDC — enxame simplificado (Swarm Intelligence) com líder sorteado do
% arquivo externo EP e reflexão de velocidade nas bordas (rebote). Arquivo
% limitado a nRep (poda aleatória se exceder).
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar];
    empty.Position = []; empty.Velocity = []; empty.Cost = []; empty.Best.Position = []; empty.Best.Cost = []; empty.IsDominated = false; pop = repmat(empty, params.nPop, 1);
    for i = 1:params.nPop, pop(i).Position = unifrnd(VarMin, VarMax, VarSize); pop(i).Velocity = zeros(VarSize); pop(i).Cost = CostFunction(pop(i).Position, dados); pop(i).Best = pop(i); end
    pop = DetermineDomination(pop); EP = pop(~[pop.IsDominated]);
    for it = 1:params.MaxIt
        for i = 1:params.nPop
            if isempty(EP), leader = pop(randi(params.nPop)); else, leader = EP(randi(numel(EP))); end
            pop(i).Velocity = params.w * pop(i).Velocity + params.c1 * rand(VarSize) .* (pop(i).Best.Position - pop(i).Position) + params.c2 * rand(VarSize) .* (leader.Position - pop(i).Position); pop(i).Position = pop(i).Position + pop(i).Velocity;
            % MELHORIA SUGERIDA (Vmax): clampe a velocidade aqui, como no MOPSO.
            % Tratamento de borda: clampa a posição e reflete (amortecida) a velocidade.
            for d = 1:nVar, if pop(i).Position(d) < VarMin || pop(i).Position(d) > VarMax, pop(i).Position(d) = max(VarMin, min(VarMax, pop(i).Position(d))); pop(i).Velocity(d) = -pop(i).Velocity(d) * 0.5; end; end
            pop(i).Cost = CostFunction(pop(i).Position, dados); if Dominates(pop(i).Cost, pop(i).Best.Cost), pop(i).Best = pop(i); else, if ~Dominates(pop(i).Best.Cost, pop(i).Cost) && rand < 0.5, pop(i).Best = pop(i); end; end
        end
        EP = [EP; pop]; EP = DetermineDomination(EP); EP = EP(~[EP.IsDominated]);
        if numel(EP) > params.nRep, costs = horzcat(EP.Cost)'; [~, idx] = unique(round(costs, 4), 'rows'); EP = EP(idx); if numel(EP) > params.nRep, EP = EP(randperm(numel(EP), params.nRep)); end; end

        if mod(it, 100) == 0 || it == params.MaxIt
            send(dq, sprintf('SI-CDC     : Iter %4d / %d concluida', it, params.MaxIt));
        end
    end
    tempo_execucao = toc(t0); repValido = EP;
end

% =====================================================================
%   MÓDULOS 9-11: MCDM (Tomada de Decisão Multicritério)
%   Constroem UMA solução por peso w1 (custo) / w2=1-w1 (GWP), escolhendo o
%   fornecedor de cada material independentemente. OBS (revisão): por
%   decidirem material a material, naturalmente usam muitos fornecedores;
%   sob a restrição K, o reparo em CostFunction os torna viáveis, porém a
%   escolha NÃO é ótima para o problema restrito (só o p-median MILP é).
% =====================================================================

function [repValido, tempo_execucao] = RunAHP(dados, params, dq)
% AHP — prioridades por normalização do inverso (quanto menor custo/GWP,
% maior a prioridade). Para cada peso, escolhe o fornecedor de maior
% prioridade global por material.
    t0 = tic; nVar = dados.nVar; nFornecedores = dados.nFornecedores; empty.Position = zeros(1, nVar); empty.Cost = []; empty.IsDominated = false; pop = repmat(empty, params.nPop, 1);
    for it = 1:params.nPop
        w1 = (it - 1) / max(1, params.nPop - 1); w2 = 1 - w1; pos = zeros(1, nVar);
        for i = 1:nVar
            C = dados.MatrizCustos(i, :); G = dados.MatrizGWP(i, :); C = max(C, 1e-6); G = max(G, 1e-6);
            P_C = (1./C) / sum(1./C); P_G = (1./G) / sum(1./G); P_Global = w1 * P_C + w2 * P_G; [~, bestIdx] = max(P_Global); pos(i) = bestIdx;
        end
        pop(it).Position = pos; pop(it).Cost = CostFunction(pos, dados);
    end
    pop = DetermineDomination(pop); repValido = pop(~[pop.IsDominated]); if ~isempty(repValido), costs = round(horzcat(repValido.Cost)', 4); [~, uniqueIdx] = unique(costs, 'rows'); repValido = repValido(uniqueIdx); end; tempo_execucao = toc(t0);
    send(dq, sprintf('AHP        : Processamento de %d cenarios concluido', params.nPop));
end

function [repValido, tempo_execucao] = RunPROMETHEE(dados, params, dq)
% PROMETHEE — fluxo líquido de preferência: para cada material, compara cada
% fornecedor com todos os outros (preferência proporcional à diferença
% normalizada) e escolhe o de maior fluxo Phi.
    t0 = tic; nVar = dados.nVar; nFornecedores = dados.nFornecedores; empty.Position = zeros(1, nVar); empty.Cost = []; empty.IsDominated = false; pop = repmat(empty, params.nPop, 1);
    for it = 1:params.nPop
        w1 = (it - 1) / max(1, params.nPop - 1); w2 = 1 - w1; pos = zeros(1, nVar);
        for i = 1:nVar
            C = dados.MatrizCustos(i, :); G = dados.MatrizGWP(i, :); maxC = max(C); minC = min(C); maxG = max(G); minG = min(G); Phi = zeros(1, nFornecedores);
            for j = 1:nFornecedores
                for k = 1:nFornecedores
                    if j == k, continue; end
                    P_C_jk = max(0, C(k) - C(j)) / max(1e-6, maxC - minC); P_G_jk = max(0, G(k) - G(j)) / max(1e-6, maxG - minG); Pi_jk = w1 * P_C_jk + w2 * P_G_jk;
                    P_C_kj = max(0, C(j) - C(k)) / max(1e-6, maxC - minC); P_G_kj = max(0, G(j) - G(k)) / max(1e-6, maxG - minG); Pi_kj = w1 * P_C_kj + w2 * P_G_kj; Phi(j) = Phi(j) + Pi_jk - Pi_kj;
                end
            end
            [~, bestIdx] = max(Phi); pos(i) = bestIdx;
        end
        pop(it).Position = pos; pop(it).Cost = CostFunction(pos, dados);
    end
    pop = DetermineDomination(pop); repValido = pop(~[pop.IsDominated]); if ~isempty(repValido), costs = round(horzcat(repValido.Cost)', 4); [~, uniqueIdx] = unique(costs, 'rows'); repValido = repValido(uniqueIdx); end; tempo_execucao = toc(t0);
    send(dq, sprintf('PROMETHEE  : Processamento de %d cenarios concluido', params.nPop));
end

function [repValido, tempo_execucao] = RunELECTRE(dados, params, dq)
% ELECTRE — concordância: para cada material, soma a concordância (pesos dos
% critérios em que j é melhor) menos a do reverso, e escolhe o de maior score.
    t0 = tic; nVar = dados.nVar; nFornecedores = dados.nFornecedores; empty.Position = zeros(1, nVar); empty.Cost = []; empty.IsDominated = false; pop = repmat(empty, params.nPop, 1);
    for it = 1:params.nPop
        w1 = (it - 1) / max(1, params.nPop - 1); w2 = 1 - w1; pos = zeros(1, nVar);
        for i = 1:nVar
            C = dados.MatrizCustos(i, :); G = dados.MatrizGWP(i, :); Score = zeros(1, nFornecedores);
            for j = 1:nFornecedores
                for k = 1:nFornecedores
                    if j == k, continue; end
                    c_jk = 0; if C(j) <= C(k), c_jk = c_jk + w1; end; if G(j) <= G(k), c_jk = c_jk + w2; end
                    c_kj = 0; if C(k) <= C(j), c_kj = c_kj + w1; end; if G(k) <= G(j), c_kj = c_kj + w2; end; Score(j) = Score(j) + c_jk - c_kj;
                end
            end
            [~, bestIdx] = max(Score); pos(i) = bestIdx;
        end
        pop(it).Position = pos; pop(it).Cost = CostFunction(pos, dados);
    end
    pop = DetermineDomination(pop); repValido = pop(~[pop.IsDominated]); if ~isempty(repValido), costs = round(horzcat(repValido.Cost)', 4); [~, uniqueIdx] = unique(costs, 'rows'); repValido = repValido(uniqueIdx); end; tempo_execucao = toc(t0);
    send(dq, sprintf('ELECTRE    : Processamento de %d cenarios concluido', params.nPop));
end

% =====================================================================
%   MÓDULOS 12-15: EXATOS E DETERMINÍSTICOS
% =====================================================================

function [repValido, tempo_execucao] = RunLexicografico(dados, params, dq)
% LEXICOGRÁFICO — otimização em duas fases: (1) minimiza o custo; (2) dentro
% de uma tolerância de custo (Tolerancia), minimiza o GWP. Retorna a melhor
% solução de compromisso encontrada.
    t0 = tic; nVar = dados.nVar; VarMin = 1; VarMax = dados.nFornecedores; VarSize = [1 nVar]; empty.Position = []; empty.Cost = []; empty.IsDominated = false; pop = repmat(empty, params.nPop, 1);
    for i = 1:params.nPop, pop(i).Position = round(unifrnd(VarMin, VarMax, VarSize)); pop(i).Cost = CostFunction(pop(i).Position, dados); end
    bestObj1 = inf; maxItPhase1 = floor(params.MaxIt / 2);
    % Fase 1: busca local para reduzir o custo (objetivo 1).
    for it = 1:maxItPhase1
        for i = 1:params.nPop
            newPos = round(Mutate(pop(i).Position, 0.2, VarMin, VarMax)); newCost = CostFunction(newPos, dados);
            if newCost(1) < pop(i).Cost(1), pop(i).Position = newPos; pop(i).Cost = newCost; end; if newCost(1) < bestObj1, bestObj1 = newCost(1); end
        end
        if mod(it, 100) == 0, send(dq, sprintf('LEXICOGRA. : (Fase 1) Iter %4d / %d concluida', it, maxItPhase1)); end
    end
    % Fase 2: respeitando custo <= bestObj1 * Tolerancia, minimiza o GWP (objetivo 2).
    bestObj2 = inf; bestSol = empty;
    for it = maxItPhase1 + 1 : params.MaxIt
        for i = 1:params.nPop
            newPos = round(Mutate(pop(i).Position, 0.2, VarMin, VarMax)); newCost = CostFunction(newPos, dados); limiteCusto = bestObj1 * params.Tolerancia;
            if newCost(1) <= limiteCusto, if newCost(2) < pop(i).Cost(2), pop(i).Position = newPos; pop(i).Cost = newCost; else, if pop(i).Cost(1) > limiteCusto, pop(i).Position = newPos; pop(i).Cost = newCost; end; end
                if newCost(2) < bestObj2, bestObj2 = newCost(2); bestSol = pop(i); end
            end
        end
        if mod(it, 100) == 0 || it == params.MaxIt, send(dq, sprintf('LEXICOGRA. : (Fase 2) Iter %4d / %d concluida', it, params.MaxIt)); end
    end
    tempo_execucao = toc(t0); if isempty(bestSol.Cost), repValido = pop(1); else, repValido = bestSol; end
end

function [repValido, tempo_execucao] = RunMILP(dados, params, dq)
% MILP — Programação Linear Inteira Mista (intlinprog), varrendo pesos.
%
%  >>> REVISÃO/MELHORIA PRINCIPAL <<<
%  Há DOIS ramos:
%   (A) SE a restrição K está ativa -> resolve o p-MEDIANAS EXATO:
%       escolhe <= K fornecedores (y_j) e atribui cada material ao melhor
%       fornecedor ABERTO. Este é o baseline EXATO correto do problema
%       RESTRITO (recupera apenas os pontos SUPORTADOS, pois usa soma
%       ponderada; a frente restrita é não-convexa).
%   (B) SENÃO -> MILP de ATRIBUIÇÃO original (problema separável/irrestrito):
%       cada material escolhe seu melhor fornecedor (trivial).
%  Variáveis no ramo (A): z_ij (material i -> fornecedor j) e y_j (forn. ativo).
    t0 = tic; nVar = dados.nVar; nFornecedores = dados.nFornecedores;
    temK = isfield(dados,'K') && ~isempty(dados.K) && dados.K < nVar;
    empty.Position = zeros(1, nVar); empty.Cost = []; empty.IsDominated = false;
    options = optimoptions('intlinprog', 'Display', 'off');

    if temK
        % ----- Ramo (A): p-medianas exato com no máximo K fornecedores -----
        nZ = nVar * nFornecedores; numVars = nZ + nFornecedores;   % [z(:); y(:)]
        % Aeq: cada material atendido por exatamente 1 fornecedor (sum_j z_ij = 1).
        ri = repelem((1:nVar)', nFornecedores); ci = (1:nZ)';
        Aeq = sparse(ri, ci, 1, nVar, numVars); beq = ones(nVar, 1);
        % A1: z_ij - y_j <= 0  (só pode usar fornecedor aberto).
        jY = nZ + (mod((0:nZ-1)', nFornecedores) + 1);            % var y_j de cada z_ij
        A1 = sparse([ci; ci], [ci; jY], [ones(nZ,1); -ones(nZ,1)], nZ, numVars);
        % A2: sum_j y_j <= K  (no máximo K fornecedores ativos).
        A2 = sparse(ones(nFornecedores,1), (nZ+1:numVars)', 1, 1, numVars);
        A = [A1; A2]; b = [zeros(nZ,1); dados.K];
        lb = zeros(numVars, 1); ub = ones(numVars, 1); intcon = (1:numVars)';
        % A frente restrita é não-convexa e tem POUCOS pontos suportados; não
        % compensa resolver centenas de pesos. Limita-se a um número modesto.
        nW = min(params.nPop, 50); pop = repmat(empty, nW, 1);
        for it = 1:nW
            w1 = (it - 1) / max(1, nW - 1); w2 = 1 - w1;
            Dt = (w1 * dados.MatrizCustos + w2 * dados.MatrizGWP)';  % nForn x nVar
            f = [Dt(:); zeros(nFornecedores,1)];                    % custo de z; y sem custo
            [x_opt, ~, flag] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, options);
            pos = ones(1, nVar);
            if ~isempty(x_opt) && flag > 0
                Zmat = reshape(round(x_opt(1:nZ)), nFornecedores, nVar)';  % nVar x nForn
                [~, pos] = max(Zmat, [], 2); pos = pos';
            end
            pop(it).Position = pos; pop(it).Cost = CostFunction(pos, dados);
        end
        send(dq, sprintf('MILP       : p-medianas exato (K=%d) — %d pesos resolvidos', dados.K, nW));
    else
        % ----- Ramo (B): MILP de atribuição (irrestrito/separável) -----
        numVars = nVar * nFornecedores;
        ri = repelem((1:nVar)', nFornecedores); ci = (1:numVars)';
        Aeq = sparse(ri, ci, 1, nVar, numVars); beq = ones(nVar, 1);   % sum_j z_ij = 1
        lb = zeros(numVars, 1); ub = ones(numVars, 1); intcon = (1:numVars)';
        pop = repmat(empty, params.nPop, 1);
        for it = 1:params.nPop
            w1 = (it - 1) / max(1, params.nPop - 1); w2 = 1 - w1;
            Dt = (w1 * dados.MatrizCustos + w2 * dados.MatrizGWP)';
            f = Dt(:);
            [x_opt, ~] = intlinprog(f, intcon, [], [], Aeq, beq, lb, ub, options);
            pos = zeros(1, nVar);
            if ~isempty(x_opt), x_mat = reshape(round(x_opt), nFornecedores, nVar)'; [~, pos] = max(x_mat, [], 2); pos = pos'; end
            pop(it).Position = pos; pop(it).Cost = CostFunction(pos, dados);
        end
        send(dq, sprintf('MILP       : Processamento de %d variacoes de peso concluido', params.nPop));
    end

    % Filtra os não-dominados e remove duplicatas em (custo, GWP).
    pop = DetermineDomination(pop); repValido = pop(~[pop.IsDominated]);
    if ~isempty(repValido), costs = round(horzcat(repValido.Cost)', 4); [~, uniqueIdx] = unique(costs, 'rows'); repValido = repValido(uniqueIdx); end
    tempo_execucao = toc(t0);
end

function [repValido, tempo_execucao] = RunBNB(dados, params, dq)
% BRANCH AND BOUND (versão construtiva) — para cada peso, escolhe por material
% o "ramo" (fornecedor) de menor custo ponderado. OBS (revisão): no problema
% IRRESTRITO isto coincide com o ótimo (separável); sob a restrição K NÃO é
% exato (o reparo em CostFunction garante viabilidade, mas não otimalidade).
% Para um B&B exato do problema restrito, usar o p-medianas do RunMILP.
    t0 = tic; nVar = dados.nVar; nFornecedores = dados.nFornecedores; empty.Position = zeros(1, nVar); empty.Cost = []; empty.IsDominated = false; pop = repmat(empty, params.nPop, 1);
    for it = 1:params.nPop
        w1 = (it - 1) / max(1, params.nPop - 1); w2 = 1 - w1; pos = zeros(1, nVar);
        for i = 1:nVar
            custo_ramos = w1 * dados.MatrizCustos(i,:) + w2 * dados.MatrizGWP(i,:); [~, bestBranch] = min(custo_ramos); pos(i) = bestBranch;
        end
        pop(it).Position = pos; pop(it).Cost = CostFunction(pos, dados);
    end
    pop = DetermineDomination(pop); repValido = pop(~[pop.IsDominated]);
    if ~isempty(repValido), costs = round(horzcat(repValido.Cost)', 4); [~, uniqueIdx] = unique(costs, 'rows'); repValido = repValido(uniqueIdx); end; tempo_execucao = toc(t0);
    send(dq, sprintf('BNB        : Processamento de %d buscas concluido', params.nPop));
end

function [repValido, tempo_execucao] = RunLagrangiana(dados, params, dq)
% RELAXAÇÃO LAGRANGIANA — para cada peso, resolve o subproblema relaxado por
% material. OBS (revisão): como não há restrição acoplante no modelo base, o
% gradiente subgradiente é nulo (grad=0) e o método reduz-se à escolha por
% material (igual ao BNB). Sob a restrição K, seria o lugar natural para
% DUALIZAR a restrição de cardinalidade (multiplicador penalizando abrir
% fornecedores novos) — fica como extensão; aqui o reparo garante viabilidade.
    t0 = tic; nVar = dados.nVar; nFornecedores = dados.nFornecedores; empty.Position = zeros(1, nVar); empty.Cost = []; empty.IsDominated = false; pop = repmat(empty, params.nPop, 1);
    for it = 1:params.nPop
        w1 = (it - 1) / max(1, params.nPop - 1); w2 = 1 - w1; lambda = zeros(nVar, 1); step_size = params.step; best_pos = zeros(1, nVar);
        for k = 1:params.MaxIt
            pos = zeros(1, nVar);
            for i = 1:nVar, custo_relaxado = (w1 * dados.MatrizCustos(i,:) + w2 * dados.MatrizGWP(i,:)) - lambda(i); [~, pos(i)] = min(custo_relaxado); end
            grad = zeros(nVar, 1); lambda = lambda + step_size * grad; best_pos = pos;
            if all(grad == 0), break; end   % sem restrição acoplante -> converge em 1 passo
        end
        pop(it).Position = best_pos; pop(it).Cost = CostFunction(best_pos, dados);
        if mod(it, 20) == 0 || it == params.nPop, send(dq, sprintf('LAGRANGIANA: Peso %d / %d resolvido', it, params.nPop)); end
    end
    pop = DetermineDomination(pop); repValido = pop(~[pop.IsDominated]);
    if ~isempty(repValido), costs = round(horzcat(repValido.Cost)', 4); [~, uniqueIdx] = unique(costs, 'rows'); repValido = repValido(uniqueIdx); end; tempo_execucao = toc(t0);
end

% =====================================================================
%   FUNCOES AUXILIARES DE CUSTO, RESTRIÇÃO E DOMINANCIA
% =====================================================================

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

function pop = DetermineDomination(pop)
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
% -------------------------------------------------------------------------
% MELHORIA SUGERIDA (vetorizada, 10-50x mais rápida) — valide e, se desejar,
% substitua a função acima por esta. Reduz drasticamente o tempo das
% metaheurísticas (que chamam DetermineDomination a cada iteração).
% function pop = DetermineDomination(pop)
%     n = numel(pop); for i = 1:n, pop(i).IsDominated = false; end
%     if n < 2, return; end
%     Costs = horzcat(pop.Cost)';                 % n x nObj
%     for i = 1:n
%         le = all(Costs <= Costs(i,:), 2);       % outros <= i em todos os obj
%         lt = any(Costs <  Costs(i,:), 2);       % outros <  i em algum obj
%         dom = le & lt; dom(i) = false;
%         if any(dom), pop(i).IsDominated = true; end
%     end
% end
% -------------------------------------------------------------------------

function pop = DetermineDominationCount(pop)
% Atribui a cada solução um 'Rank' = 1 + nº de soluções que a dominam
% (rank 1 = não-dominada). Usado pelo MOGA.
    n = numel(pop); for i = 1:n, pop(i).Rank = 1; end
    for i = 1:n, for j = 1:n, if i ~= j && Dominates(pop(j), pop(i)), pop(i).Rank = pop(i).Rank + 1; end; end; end
end

function pop = NonDominatedSort(pop)
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

% =====================================================================
%   [V5] MÃ“DULOS 16-17: MULTI-CHILD DIFFERENTIAL EVOLUTION MULTIOBJETIVO
%   Fontes dos nÃºcleos mono-objetivo:
%     [MCDE]     Storn & Price, "Multi-Child DE â€” A Massively Parallel
%                Differential Evolution Algorithm", IEEE CIES/SSCI 2025.
%     [MC-SHADE] Storn, "Comparison of MCDE and MC-SHADE for Massively
%                Parallel Optimization", IEEE CEC 2026 (eqs. 1-14).
%   AdaptaÃ§Ã£o multiobjetivo (contribuiÃ§Ã£o desta dissertaÃ§Ã£o):
%     - a prÃ©-seleÃ§Ã£o (argmin f) e a seleÃ§Ã£o final (f(w)<=f(x)) do original
%       sÃ£o substituÃ­das por regras de DOMINÃ‚NCIA DE PARETO no estilo GDE3
%       (Kukkonen & Lampinen, CEC 2005): campeÃ£o que domina o pai substitui;
%       mutuamente nÃ£o-dominados vÃ£o para um pool e a populaÃ§Ã£o volta a nPop
%       por seleÃ§Ã£o ambiental NSGA-II (rank + crowding), reutilizando os
%       helpers NonDominatedSort/CalcCrowdingDistance/SortAndTruncateNSGA2;
%     - repositÃ³rio externo de nÃ£o-dominados com capacidade nRep, podado
%       por crowding (Ã© a frente devolvida â€” mesma convenÃ§Ã£o do harness);
%     - no MO-MC-SHADE, o "sucesso" da adaptaÃ§Ã£o de F/Cr passa a ser
%       "campeÃ£o DOMINA o pai" e o peso delta_i da regra de Peng (eq. 11)
%       vira a norma da melhoria NORMALIZADA pela faixa corrente dos
%       objetivos na populaÃ§Ã£o (substituto vetorial de |f(w)-f(x)|).
%   As avaliaÃ§Ãµes continuam via CostFunction => a restriÃ§Ã£o K (reparo
%   Baldwiniano) vale automaticamente, como nas demais tÃ©cnicas.
% =====================================================================

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

% ---------------------------------------------------------------------
%  Helpers das tÃ©cnicas multi-child (V5)
% ---------------------------------------------------------------------

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

