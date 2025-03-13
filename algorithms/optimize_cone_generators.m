function [u,results] = optimize_cone_generators(G,Av,timelimit);

% Given G and Av, solve 
% 
% min_{u=Gx, x>=0} <u,Av> such that ||u|| <= 1. 

[m,n] = size(G); 
%% Construct model and solve the Gurobi model for (Problem 1) above
model.modelsense = 'min';
model.vtype = 'C';
model.obj = G'*Av; 
model.A = sparse(1,n);
model.rhs = 0;
model.sense = '=';
model.lb = zeros(n,1);
% ||Gx|| <= 1 
model.quadcon(1).Qc = sparse([G'*G]);
model.quadcon(1).q  = zeros(n,1);
model.quadcon(1).rhs = 1.0;
model.quadcon(1).sense = '<'; 
if nargin >= 4
    params.TimeLimit = timelimit;
end
params.Outputflag = 0; % display on/off
%gurobi_write(model, 'qcp.lp');
results = gurobi(model, params);
u = G*results.x(1:n); 
if norm(u) < 0.1 % meaning ||u|| is zero because if it is positive, it is one 
    % pick column of G that minimizes G^T Av
    [~,b] = min(model.obj); 
    u = G(:,b); 
    u = u/norm(u); 
end
