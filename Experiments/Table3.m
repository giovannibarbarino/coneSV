% Numerical comparison for Gurobi, BFAS, E-AO and
%  SRPL for the problem of finding the maximum edge biclique
%   in four different bipartite graphs
%
% Reports the maximum edge biclique found in the timelimit
%  (10 seconds) for Gurobi and BFAS
% The reported numbers for E-AO and SRPL are instead the 
%  values found at 10 seconds in 100 runs and their averages
% Gurobi cannot be executed on the last graph due to its 
%  excessive size

% The largest edge bicliques found in the i-th graph in the 
%  100 iterations are reported in EAO_all(i,:) and 
%   SRPL_all(i,:), while their averages are reported in 
%    EAO_mean(i) and SRPL_mean(i)

% The largest edge bicliques found by Gur and BFAS in the
%  i-th graph are reported in BFAS(i) and Gur(i)
% Notice that Gur(4) does not exists since Gurobi crashes 
%  on the last graph

% The four bipartite graphs are taken from the dataset in
%  https://github.com/shahamer/maximum-biclique-benchmark 
%   Shaham, E.: maximum biclique benchmark. (Dic 2019)

% include algorithms
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename)); cd('..'); 
addpath(genpath('./'));

timelimit = 10;
Ninitpoint = 100;

EAO_all = zeros(4,Ninitpoint); EAO_mean = zeros(1,4);
SRPL_all = zeros(4,Ninitpoint); SRPL_mean = zeros(1,4);
BFAS = zeros(1,4); Gur = zeros(1,3);

for i = 1:4
    Adj = readmatrix(strcat("Biclique_matrix_",num2str(i),".txt"));
    n = max(size(Adj)); A = -Adj + (1-Adj)*n;
    
    % E-AO
    options.cone.P = 'nonnegort'; options.cone.Q = 'nonnegort';
    options.accuracy = 1e-10; options.maxiter = 500;
    options.beta = .5; tol = 1e-8;
    for j = 1:Ninitpoint % for all random initial points
        emin = 0; T = 0;
        while T < timelimit
            tic;
            options.v0 = rand(n,1);
            options.v0 = options.v0/norm(options.v0);
            [~,~,e] = AlternatingOptimization(A,options); %E-AO
            T = T + toc;
            if e(end) < emin - tol, emin = e(end); end
        end
        EAO_all(i,j) = emin^2; % the optimal value is in EAO_all
    end
    %average of all runs
    EAO_mean(i) = mean(EAO_all(i,:));

    % SRPL
    G = eye(size(A,1)); H = eye(size(A,2));
    tol = 1e-8; mu1 = 0.25; mu2 = 0.01;
    for j = 1:Ninitpoint % for all random initial points
        emin = 0; T = 0;
        while T < timelimit
            tic;
            x0 = rand(size(A,1),1); x0 = x0/sum(x0);
            y0 = rand(size(A,2),1); y0 = y0/sum(y0);
            [~,~,la,~,~] = srpl_poly(A,G,H,x0,y0,mu1,mu2); %SRPL
            T = T + toc;
            if la < emin - tol, emin = la; end
        end
        SRPL_all(i,j) = emin^2; % the optimal value is in SRPL_all
    end
    %average of all runs
    SRPL_mean(i) = mean(SRPL_all(i,:));


    % BFAS
    [lamvec,~,~,~] = bfas_timestamps_test(G,H,A,timelimit);
    BFAS(i) = lamvec(end)^2;


    % Gurobi
    if i < 4
        [val,~] = generators_gurobi_uv_test(G,H,A,timelimit);
        Gur(i) = val^2;
    end

end
