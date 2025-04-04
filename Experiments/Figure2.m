% E-AO, SRPL and Gurobi compared on the problem of finding the 
%  maximum angle between the Schur cone and the nonnegative orthant,
%   and between the Schur cone and itself in dimension n = 200, with a 
%    timelimit of 10 seconds

% For E-AO and SRPL, 100 iterations from random generated points 
%  are performed, and their performance over time is computed together
%   with their average
% The data for the 100 iterations are reported in EAO_all and SRPL_all
%  where EAO_all(i,t), SRPL_all(i,t) are the best values found on iteration
%   i at time t/100 seconds, expressed as fractions of the angle pi
% The averages at time t/100 seconds are reported in EAO_mean(t) and
%  SRPL_mean(t), expressed as fractions of the angle pi

% The data for Gurobi are logged on screen


%% include algorithms
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename)); cd('..'); 
addpath(genpath('./'));

n = 200; timelimit = 10; Ninitpoint = 100;
% timemult is how many timestamps we log each second for E-AO and SRPL 
timemult = 100; Ndisctime = timelimit*timemult;

%% Schur - Positive Orthant

A = eye(n); H = eye(n);
% G generates the Schur cone
G = (diag(ones(1,n-1),1) - diag(ones(1,n)))/sqrt(2);
G = G(:,2:end);

% E-AO
options.cone.P = 'generator'; options.cone.Q = 'nonnegort';
options.G = G;
options.accuracy = 1e-6; options.maxiter = 500;
options.beta = .5; tol = 1e-8;

EAO_all = zeros(Ninitpoint,Ndisctime+1);
for j = 1:Ninitpoint % for all random initial points
    T = 0; EAO = []; EAO_times = []; emin = 1;
    while T < timelimit 
        tic;
        options.v0 = rand(n,1);
        options.v0 = options.v0/norm(options.v0);
        [~,~,e] = AlternatingOptimization(A,options); %E-AO
        T = T + toc;
        if e(end) < emin - tol % if the objective has improved
            EAO = [EAO e(end)]; emin = e(end);
            EAO_times = [EAO_times T];
        end
    end
    % for any timestamp, report the best value found
    pip = 0;
    for i = 1:length(EAO_times)
        if EAO_times(i) > timelimit, break; end
        pipend = floor(EAO_times(i)*timemult);
        EAO_all(j,pip+1:pipend+1) = EAO(i);
        pip = pipend+1;
    end
    if isempty(EAO), EAO_all(j,:) = ones(1,Ndisctime+1);
    else, EAO_all(j,pip+1:end) = EAO(i); end
end
EAO_all = acos(EAO_all)/pi;
%average of all runs
EAO_mean = mean(EAO_all,1);


% SRPL
tol = 1e-8; mu1 = 0.25; mu2 = 0.01;

SRPL_all = zeros(Ninitpoint,Ndisctime+1);
for j = 1:Ninitpoint % for all random initial points
    T = 0; SRPL = []; SRPL_times = []; emin = 1;
    while T < timelimit 
        tic;
        x0 = rand(n-1,1); x0 = x0/sum(x0);
        y0 = rand(n,1); y0 = y0/sum(y0);
        [~,~,la,~,~] = srpl_poly(A,G,H,x0,y0,mu1,mu2); %SRPL
        T = T + toc;
        if la < emin - tol % if the objective has improved
            SRPL = [SRPL la]; emin = la;
            SRPL_times = [SRPL_times T];
        end
    end
    % for any timestamp, report the best value found
    pip = 0;
    for i = 1:length(SRPL_times)
        if SRPL_times(i) > timelimit, break; end
        pipend = floor(SRPL_times(i)*timemult);
        SRPL_all(j,pip+1:pipend+1) = SRPL(i);
        pip = pipend+1;
    end
    if isempty(SRPL), SRPL_all(j,:) = ones(1,Ndisctime+1);
    else, SRPL_all(j,pip+1:end) = SRPL(i); end
end
SRPL_all = acos(SRPL_all)/pi;
%average of all runs
SRPL_mean = mean(SRPL_all,1);



% Gurobi
generators_gurobi_uv_test(G,H,A,timelimit,1);




%% Schur - Schur

A = eye(n);
% G generates the Schur cone
G = (diag(ones(1,n-1),1) - diag(ones(1,n)))/sqrt(2);
G = G(:,2:end); H = G;

% E-AO
options.cone.P = 'generator'; options.cone.Q = 'generator';
options.G = G; options.H = H;
options.accuracy = 1e-10; options.maxiter = 500;
options.beta = .5; tol = 1e-8;

EAO_all = zeros(Ninitpoint,Ndisctime+1);
for j = 1:Ninitpoint % for all random initial points
    T = 0; EAO = []; EAO_times = []; emin = 1;
    while T < timelimit 
        tic;
        options.v0 = rand(n,1);
        options.v0 = options.v0/norm(options.v0);
        [~,~,e] = AlternatingOptimization(A,options); %E-AO
        T = T + toc;
        if e(end) < emin - tol % if the objective has improved
            EAO = [EAO e(end)]; emin = e(end);
            EAO_times = [EAO_times T];
        end
    end
    % for any timestamp, report the best value found
    pip = 0;
    for i = 1:length(EAO_times)
        if EAO_times(i) > timelimit, break; end
        pipend = floor(EAO_times(i)*timemult);
        EAO_all(j,pip+1:pipend+1) = EAO(i);
        pip = pipend+1;
    end
    if isempty(EAO), EAO_all(j,:) = ones(1,Ndisctime+1);
    else, EAO_all(j,pip+1:end) = EAO(i); end
end
EAO_all = acos(EAO_all)/pi;
%average of all runs
EAO_mean = mean(EAO_all,1);


% SRPL
tol = 1e-8; mu1 = 1; mu2 = 1;

SRPL_all = zeros(Ninitpoint,Ndisctime+1);
for j = 1:Ninitpoint % for all random initial points
    T = 0; SRPL = []; SRPL_times = []; emin = 1;
    while T < timelimit 
        tic;
        x0 = rand(n-1,1); x0 = x0/sum(x0);
        y0 = rand(n-1,1); y0 = y0/sum(y0);
        [~,~,la,~,~] = srpl_poly(A,G,H,x0,y0,mu1,mu2);
        T = T + toc;
        if la < emin - tol % if the objective has improved
            SRPL = [SRPL la]; emin = la;
            SRPL_times = [SRPL_times T];
        end
    end
    % for any timestamp, report the best value found
    pip = 0;
    for i = 1:length(SRPL_times)
        if SRPL_times(i) > timelimit, break; end
        pipend = floor(SRPL_times(i)*timemult);
        SRPL_all(j,pip+1:pipend+1) = SRPL(i);
        pip = pipend+1;
    end
    if isempty(SRPL), SRPL_all(j,:) = ones(1,Ndisctime+1);
    else, SRPL_all(j,pip+1:end) = SRPL(i); end
end
SRPL_all = acos(SRPL_all)/pi;
%average of all runs
SRPL_mean = mean(SRPL_all,1);



% Gurobi
A = 9*eye(n); 
generators_gurobi_uv_test(G,H,A,timelimit,1);
