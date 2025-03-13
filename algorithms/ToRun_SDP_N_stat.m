function M = ToRun_SDP_N_stat(n,rho_x,rho_y)

% [x,y,la,it,status] = srpl_poly(A,G,H,xk,yk,rho_x,rho_y)
%
% Sequential Regularized Partial Linearization for solving
%
%           min_{u in P, v in Q, u,v \ne 0} <u, Av> / ||u||||v||    (1)
%
% where P is the cone of SPD matrices and Q is the cone of nonnegative
% symmetric matrices 

% The algorithm starts from a randomly generated feasible point

% *** Input ***
%   n              : size of the matrices
%
%   rho_x, rho_y   : Penalization parameters. 
%
% *** Output ***
%   M              : vector containing 
%                    [fval_rspl, it_rspl, cputime_rspl]
%                    fval_rspl is the optimal value found
%                    it_rspl is the number of iterations
%                    cputime_rspl is the runtime
                
%% Initial point
x0 = rand(n,1);
x0 = x0/sum(x0);
X0 = sparse(diag(x0));
N = 0.5*n*(n+1);
y0 = generator_simplex(N,2);

%% Solve with SRPL
tt = tic;
[~,~,la,it,~] = srpl_SDP_N(X0,y0,rho_x,rho_y);
cputime_rspl  = toc(tt);
fval_rspl     = la;
it_rspl       = it;

M = [fval_rspl, it_rspl, cputime_rspl];
end


function x0 = generator_simplex(n,method)

switch method
    case 1
        
        x0    = zeros(n,1);
        somma = 0;
        for j = 1:n-1
            x0(j,1) = rand(1)*(1-somma);
            somma = somma+x0(j,1);
        end
        x0(n,1) = 1-somma;
        
    case 2
        
        x0 = rand(n,1);
        x0 = x0/sum(x0);
        
end

end
