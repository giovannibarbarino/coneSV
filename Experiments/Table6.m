% Numerical comparison for Gurobi, BFAS, E-AO and
%  SRPL for the problem of finding the maximum angle between
%   the PSD cone and the nonnegative symmetric cone, both
%    restricted to the subalgebra of circulant matrices
%
% Reports the maximum angle found in the timelimit
%  (10 seconds) for Gurobi and BFAS
% The reported numbers for E-AO and SRPL are instead the 
%  values found at 10 seconds in 100 runs and their averages


% The largest angle expressed as a fraction of pi found in
%  the 100 iterations are reported in EAO_all(i,:) and 
%   SRPL_all(i,:), while their averages are reported in 
%    EAO_mean(i) and SRPL_mean(i)

% The largest angle expressed as a fraction of pi found by
%  Gur and BFAS are reported in BFAS(i) and Gur(i)




% include algorithms
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename)); cd('..'); 
addpath(genpath('./'));

timelimit = 10;
Ninitpoint = 100;

ndim = 17:2:27; nexp = length(ndim);
EAO_all = zeros(nexp,Ninitpoint); EAO_mean = zeros(1,nexp);
SRPL_all = zeros(nexp,Ninitpoint); SRPL_mean = zeros(1,nexp);
BFAS = zeros(1,nexp); Gur = zeros(1,nexp);

for i = 1:nexp

    n = ndim(i); N = n; 
    m = (N-1)/2; Aux = fft(eye(N))/sqrt(N); 
    Aux = real(Aux); A = 2*Aux(2:m+1,2:m+1);
    G = eye(m); H = G;

    
    % E-AO
    options.cone.P = 'nonnegort'; options.cone.Q = 'nonnegort';
    options.accuracy = 1e-10; options.maxiter = 500;
    options.beta = .5; tol = 1e-8;
    for j = 1:Ninitpoint % for all random initial points
        emin = 0; T = 0;
        while T < timelimit
            tic;
            options.v0 = rand(m,1);
            options.v0 = options.v0/norm(options.v0);
            [~,~,e] = AlternatingOptimization(A,options); %E-AO
            T = T + toc;
            if e(end) < emin - tol, emin = e(end); end
        end
        EAO_all(i,j) = acos(emin)/pi; % the optimal value is in EAO_all
    end
    %average of all runs
    EAO_mean(i) = mean(EAO_all(i,:));

    % SRPL
    tol = 1e-8; mu1 = 0.25; mu2 = 0.01;
    for j = 1:Ninitpoint % for all random initial points
        emin = 0; T = 0;
        while T < timelimit
            tic;
            x0 = rand(m,1); x0 = x0/sum(x0);
            y0 = rand(m,1); y0 = y0/sum(y0);
            [~,~,la,~,~] = srpl_poly(A,G,H,x0,y0,mu1,mu2); %SRPL
            T = T + toc;
            if la < emin - tol, emin = la; end
        end
        SRPL_all(i,j) = acos(emin)/pi; % the optimal value is in SRPL_all
    end
    %average of all runs
    SRPL_mean(i) = mean(SRPL_all(i,:));


    % BFAS
    [lamvec,~,~,~] = bfas_timestamps_test(G,H,A,timelimit);
    BFAS(i) = acos(lamvec(end))/pi;


    % Gurobi

    [val,~] = generators_gurobi_uv_test(G,H,A,timelimit);
    Gur(i) = acos(val)/pi;


end
