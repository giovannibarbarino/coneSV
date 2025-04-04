% Numerical comparison for Gurobi and BFAS for 
%  the problem of finding the maximum angle between the
%   Schur cone and the positive orthant cone with dimensions
%    in the vector ndim

% Returns the optimal angle found by Gurobi and BFAS
%  in the timelimit (60 seconds) and the actual elapsed time

% The values BFAS(1,:) and Gur(1,:) are the optimal angles
%  found as a fraction of pi

% The elapsed times are in BFAS(2,:) and Gur(2,:)


% include algorithms
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename)); cd('..'); 
addpath(genpath('./'));

ndim = [5 10 20 50 100 200 500]; nexp = length(ndim);
BFAS = zeros(2,nexp); Gur = zeros(2,nexp);
timelimit = 60;

for i = 1:nexp

    n = ndim(i);
    A = eye(n); H = eye(n);
    % G generates the Schur cone
    G = (diag(ones(1,n-1),1) - diag(ones(1,n)))/sqrt(2);
    G = G(:,2:end);
    % BFAS
    [lamvec,~,~,timest] = bfas_timestamps_test(G,H,A,timelimit);
    BFAS(1,i) = acos(lamvec(end))/pi; BFAS(2,i) = timest(end);
    % Gurobi
    [val, time] = generators_gurobi_uv_test(G,H,A,timelimit);
    Gur(1,i) = acos(val)/pi; Gur(2,i) = time;

end

