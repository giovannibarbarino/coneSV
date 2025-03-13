% Numerical comparison for E-AO and SRPL for 
%  different dimensions for the problem of finding
%   the maximum angle between the PSD cone and the 
%    nonnegative symmetric cone
% 
% In EAO_best, EAO_ave, SRPL_best, SRPL_ave are 
%  repported the best and average angles found over 
%   10000 random initializations for EAO and SRPL, 
%    expressed as fractions of pi
% 
% In EAO_meantime and SRPL_meantime are reported the 
%  average elapsed time for EAO and SRPL

% In EAO_all(i,j,1) and SRPL_all(i,j,1) are reported 
%  the angles found on randum input number j for the
%   problem with dimension nlist(i)

% In EAO_all(i,j,2) and SRPL_all(i,j,2) are reported 
%  the running times on randum input number j for the
%   problem with dimension nlist(i)

% include algorithms
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename)); cd('..'); 
addpath(genpath('./'));

nlist = 20:10:60;
Ninitpoint = 1;  % must be at least 2

% E-AO
EAO_all = zeros(length(nlist),Ninitpoint,2);
indnlist = 0;
options.cone.P = 'semidefin';
options.cone.Q = 'nonnegort';
options.beta = .5;
options.accuracy = 1e-6;
options.maxiter = 500;
for n = nlist
    indnlist = indnlist + 1;
    A = eye(n^2);
    for j = 1:Ninitpoint
        tic;
        options.v0 = rand(n^2,1);
        options.v0 = options.v0/norm(options.v0);
        [u,v,e] = AlternatingOptimization(A,options);
        EAO_all(indnlist,j,2) = toc;
        EAO_all(indnlist,j,1) = acos(e(end))/pi;
    end
end
EAO_best = max(EAO_all(:,:,1)');
EAO_ave = mean(EAO_all(:,:,1)');
EAO_meantime = mean(EAO_all(:,:,2)');

% SRPL

mu1 = 0.1; 
mu2 = 5;
SRPL_all = zeros(length(nlist),Ninitpoint,2);
indnlist = 0;
for nn = nlist
    indnlist = indnlist + 1;
    for j = 1:Ninitpoint
        tic;
        M = ToRun_SDP_N_stat(nn,mu1,mu2);
        e = min(M(:,1));
        SRPL_all(indnlist,j,2) = toc;
        SRPL_all(indnlist,j,1) = acos(e(end))/pi;
    end
end
SRPL_best = max(SRPL_all(:,:,1)');
SRPL_ave = mean(SRPL_all(:,:,1)');
SRPL_meantime = mean(SRPL_all(:,:,2)');

