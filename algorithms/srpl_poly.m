function [x,y,la,it,status] = srpl_poly(A,G,H,xk,yk,rho_x,rho_y)

% [x,y,la,it,status] = srpl_poly(A,G,H,xk,yk,rho_x,rho_y)
%
% Sequential Regularized Partial Linearization for solving
%
%           min_{u in P, v in Q, u,v \ne 0} <u, Av> / ||u||||v||    (1)
%
% where
%       P = {u | u = G*x, x >= 0},  Q = {u | u = H*x, x >= 0}

% *** Input ***
%   A              : m-by-n matrix
%
%   G, H           : The columns of G (resp. H) are the rays of P (resp. Q)
%                    hence G (m-by-k) u = G*x in P for some x >= 0.
%                      and H (n-by-p) v = H*y in Q for some y >= 0.
%                    G, H must have positively independent columns.
%   xk, yk         : initializations for the algorithm.
%   rho_x, rho_y   : Penalization parameters. Default values = 1.
%
% *** Output ***
%   (x,y)          : approximate solution to Problem (1)
%   la             : final value of u'*A*v / ||u||||v||
%   it             : number of iterations.
%   status         : 0 if it surpasses itmax. 1 otherwise.
%                    default itmax = 5000

p = length(xk);
q = length(yk);

tol_1 = 1e-6;
tol_2 = 1e-6;

it = 0;
itmax = 5000;
status = 1;

% n = size(G,1);
% p = size(G,2);
% q = size(H,2);

if nargin<=4 % default values for prox-parameters
    rho_x = 1;
    rho_y = 1;
end


Phi = @(x,y)( x'*G'*A*H*y/(norm(G*x)*norm(H*y)) ); 
% fla = @(x,y,la)( x'*G'*A*H*y-la*(norm(G*x)*norm(H*y)));
proj_simplex = @(y) max(y-max((cumsum(sort(y,1,'descend'),1)-1)./(1:size(y,1))'),0);


lak_o = 1e6;
count = 0;


while (it < itmax)
    gx = G*xk; 
    hy = H*yk;
    Atgx = A'*gx;
    Ahy = A*hy;
    ngx = norm(gx);
    nhy = norm(hy);
    lak = (gx'*Ahy)/(ngx*nhy);
    c1 = G'*(Ahy-lak*(nhy/ngx)*gx);
    c2 = H'*(Atgx-lak*(ngx/nhy)*hy);
    if rho_x == 0
        Aineqx = -eye(p); bineqx = zeros(p,1);
        Aeqx = ones(1,p); beqx=1;
        xk1 = linprog(c1,Aineqx,bineqx,Aeqx,beqx);
    else
    % Find xk1 by solving the regularized problem given (xk,yk)
        b_x = -(c1/rho_x-xk);
        xk1 = proj_simplex(b_x);
    end
    if rho_y == 0
        Aineqy = -eye(q); bineqy = zeros(q,1);
        Aeqy = ones(1,q); beqy=1;
        yk1 = linprog(c2,Aineqy,bineqy,Aeqy,beqy);
    else
    % Find yk1 by solving the regularized problem given (xk,yk)
        b_y = -(c2/rho_y-yk);
        yk1 = proj_simplex(b_y);
    end

%-----------------------------------------------------------------    
    d1k = xk1-xk;
    d2k = yk1-yk;
    
    L1_d1 = d1k'*c1;
    L2_d2 = d2k'*c2;
    
    gradPhi_d = (L1_d1+L2_d2)/(ngx*nhy);
    
    
    if (abs(lak_o-lak)<1e-5)
        count = count+1;
    else 
        count = 0; % count algorithm stanilize over the last five iterations
        lak_o = lak;
    end
        
    if((abs(L1_d1)< tol_1)&&(abs(L2_d2)< tol_2))  &&count>=5   % Crit 1 
        x = xk;
        y = yk;
        la = lak;      
        return;
    else
         
        [xk1,yk1] = LSearch(Phi,xk,d1k,yk,d2k,gradPhi_d); 
        % flak =  fla(xk1,yk1,lak);
        it = it + 1;
        xk = xk1;
        yk = yk1;
    end
    
end

if(it >= itmax)
    status = 0;
    x      = [];
    y      = [];
    la     = lak;
end


end

function [xkk,ykk,i] = LSearch(f,xk,d1,yk,d2,grad_d)
        rho   = 0.2; %rho \in (0,1);      
        alpha = 1e-3; % alpha \in (0,1);
        t     = 1; % initial t
        fk    = f(xk,yk);
        i     = 0;
        nd2   = grad_d;
        while (i<20) && (f(xk+t*d1,yk+t*d2)>fk+alpha*t*nd2)
            t = t*rho; i=i+1;
%               fprintf('i: %d ,f(x+td): %2.2e, fk+alpha*t*grad: %2.2e, gk: %2.2e \n',i, f(xk+t*d1,yk+t*d2),fk+alpha*t*nd2,grad_d);      
        end
        if i ==20
            disp('LS failed')
        end
        xkk = xk + t*d1;
        ykk = yk + t*d2;
end

