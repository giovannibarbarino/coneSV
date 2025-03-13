function [val,time] = generators_gurobi_uv_test(G,H,A,timelimit,OF)

% [val,time] = generators_gurobi_uv_test(G,H,A,timelimit)
%
% Gurobi for solving
%
%           min_{u in P, v in Q, ||u||=||v||=1} <u, Av>     (1)
%
% where
%       P = {u | u = G*x, x >= 0},  Q = {u | u = H*x, x >= 0}

% *** Input ***
%   A              : m-by-n matrix
%
%   G, H           : The columns of G (resp. H) are the rays of P (resp. Q)
%                     hence G (m-by-p) u = G*x in P for some x >= 0.
%                      and H (n-by-q) v = H*y in Q for some y >= 0.
%   timelimit      : time limit for the algorithm
%   OF             : Output Flag. If it is not zero, it prints the partial
%                     data found during the run by Gurobi, including 
%                      timestamps and upper/lower bounds. The default is 0
%
% *** Output ***
%   val            : approximate optimal value for <u, Av> in Problem (1)
%   time           : runtime of the algorithm

if nargin<=4 % default values for OF
    OF = 0;
end

G = normc(G); H = normc(H);
[m1,n1] = size(G); 
[m2,n2] = size(H); 
n = m1+m2+n1+n2; 
GtAH = G'*A*H; % size n1 x n2

[minGAH,indminGAH] = min(GtAH(:)); 
if minGAH >= 0 
    disp('Since G''*A*H >= 0, the optimal solution are extreme rays.'); 
    j = floor((indminGAH-1)/n2)+1;
    i = indminGAH-(j-1)*n2; 
    x = zeros(n1,1); 
    x(i) = 1; 
    y = zeros(n2,1); 
    y(j) = 1; 
    u = G*x; 
    y = H*y; 
    results = 'trivial case'; 
    return; 
end
%% Construct model and solve the Gurobi model for (Problem 1) above
model.modelsense = 'min';
model.vtype = 'C';
model.Q = sparse(n,n); 
model.Q(1:m1,m1+1:m1+m2) = 0.5*A; % u,v x,y 
model.Q(m1+1:m1+m2,1:m1) = 0.5*A'; 
% Equalities: u = Gx; v = Hy 
model.A = [speye(m1)    sparse(m1,m2) -G           sparse(m1,n2); ... 
           zeros(m2,m1) speye(m2)     zeros(m2,n1) -H ]; 
model.rhs = zeros(m1+m2,1);
model.sense = '='; 
model.lb(1:m1+m2) = -1; % Bounds for u,v: -1 <= u,v <= 1 
model.ub(1:m1+m2) = 1;
model.lb(m1+m2+1:n) = zeros(n1+n2,1); % Bounds for x,y >= 0 
model.ub(m1+m2+1:n) = Inf; 
% This is a non-convex QCQP
params.NonConvex = 2;
% ||u|| <= 1 
model.quadcon(1).Qc = sparse(n,n); 
model.quadcon(1).Qc(1:m1,1:m1) = eye(m1); 
model.quadcon(1).q  = zeros(n,1);
model.quadcon(1).rhs = 1.0;
model.quadcon(1).sense = '<'; 
% ||v|| <= 1 
model.quadcon(2).Qc = sparse(n,n); 
model.quadcon(2).Qc(m1+1:m1+m2,m1+1:m1+m2) = eye(m2); 
model.quadcon(2).q  = zeros(n,1);
model.quadcon(2).rhs = 1.0;
model.quadcon(2).sense = '<'; 
if nargin >= 4
    params.TimeLimit = timelimit;
end
params.Outputflag = OF; % display on/off
% gurobi_write(model, 'qcp.lp');
results = gurobi(model, params);
time = results.runtime;
val = results.objval;