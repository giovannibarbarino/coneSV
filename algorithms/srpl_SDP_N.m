function [X,y,la,it,status] = srpl_SDP_N(Xk,yk,rho_x,rho_y)
% 

tol_1 = 1e-6;
tol_2 = 1e-6;

it     = 0;
itmax  = 5e3; 
status = 1;

n = size(Xk,1);
N = numel(yk);

if nargin<4 % default values for prox-parameters
    rho_x = 1;
    rho_y = .1;
end


Phi = @(X,Hy)( trace(X'*Hy)/(norm(X,'fro')*norm(Hy,'fro')) ); 
fla = @(X,Hy,la)( trace(X'*Hy)-la*(norm(X,'fro')*norm(Hy,'fro')));
proj_simplex = @(y) max(y-max((cumsum(sort(y,1,'descend'),1)-1)./(1:size(y,1))'),0);

lak_o = 1e6;
count = 0;
while (it < itmax)
    
    Hyk  = matrixize(yk,n);     
    nXk  = norm(Xk,'fro');
    nHyk = norm(Hyk,'fro');
    
    lak = trace(Xk*Hyk)/(nXk*nHyk); % Xk is symmetric

%----------------------------------------------------------------   
% Find Xk1 by solving the regularized problem given (Xk,yk)
    
    aux1    = (Hyk-lak*(1/nXk)*nHyk*Xk);
    Bk      = Xk-(1/rho_x)*aux1;   
    [Q,laB] = eig(0.5*(Bk+Bk')); %instead of eig(Bk) to avoid the problem the it would compute compl eig although the matrix is symm  
    la_bar  = proj_simplex(diag(laB));
    Xk1     = Q*diag(la_bar)*Q';

  
% Find yk1 by solving the regularized problem given (xk,yk)
    aux2 = Xk-lak*nXk*(1/nHyk)*Hyk;
    HTX  = vectorize(aux2,n,N);
    bk   = yk-1/rho_y*(HTX);
    yk1  = proj_simplex(bk);

%-----------------------------------------------------------------    
    d1k = Xk1-Xk;
    d2k = yk1-yk;
    
  
    L1_d1 = trace(d1k'*aux1);
    L2_d2 =  d2k'*HTX; %same thing than trace(matrixize_v2(d2k,n)'*aux2);
    
    gradPhi_d = (L1_d1+L2_d2)/(nXk*nHyk);

    if (abs(lak_o-lak)<1e-7) 
        count = count+1;
    else 
        count = 0;
        lak_o = lak;
    end
        
    if((abs(L1_d1)< tol_1)&&(abs(L2_d2)< tol_2))&&count>=5  % Crit 1 
        X  = Xk;
        y  = yk;
        la = lak;      
%         fprintf('it: %d, L1: %2.2e, L2: %2.2e, lak: %2.6e, flak: %2.2e, ||d1||_inf: %2.2e, ||d2||_inf: %2.2e   \n',it,L1_d1,L2_d2,lak,flak,norm(d1k,inf),norm(d2k,inf));
        return;
    else
         
        [Xk1,yk1] = LSearch(Phi,Xk,d1k,yk,d2k,gradPhi_d,n); 
        flak      =  fla(Xk1,matrixize(yk1,n),lak);

%         if(mod(it,1)==0)
%             fprintf('it: %d, L1: %2.2e, L2: %2.2e, lak: %2.6e, flak: %2.2e, ||d1||_inf: %2.2e, ||d2||_inf: %2.2e \n',it,L1_d1,L2_d2,lak,flak,norm(vectorize(d1k,n,N),inf),norm(d2k,inf));
%         end
        
        it = it + 1;
        Xk = Xk1;
        yk = yk1;
    end
    
end

if(it >= itmax)
    status = 0;
    X = [];
    y = [];
    la = lak;
end


end

function [xkk,ykk,i] = LSearch(f,Xk,d1,yk,d2,grad_d,n)
        rho   = 0.2;  %rho \in (0,1);      
        alpha = 1e-3; % alpha \in (0,1);
        t     = 1;    % initial t
        Hyk   = matrixize(yk,n);
        fk    = f(Xk,Hyk);
        i     = 0;
        nd2   = grad_d;
        
        Hykd2  = matrixize(yk+t*d2,n); 
        while (i<20) && (f(Xk+t*d1,Hykd2)>fk+alpha*t*nd2)
            t = t*rho; i=i+1;
            Hykd2  = matrixize(yk+t*d2,n); 
%                fprintf('i: %d ,f(x+td): %2.4e, fk+alpha*t*grad: %2.4e, gk: %2.4e \n',i, f(Xk+t*d1,Hykd2),fk+alpha*t*nd2,grad_d);      
        end
        if i ==20
            disp('LS failed')
        end
        xkk = Xk + t*d1;
        ykk = yk + t*d2;
end

function HTX = vectorize(X,n,N)
% Generate a column vector HTX from a symmetric matrix X of size n x n.

    X        = X+ triu(X, 1);
    HTX      = zeros(N,1);
    indices  = triu(true(n));
    HTX(1:N) = X(indices);

end

function Hy = matrixize(y, n)
% Generate a symmetric matrix Hy from a column vector y of length n(n+1)/2.
Hy          = zeros(n);
indices     = triu(true(n));
Hy(indices) = y;
Hy          = Hy + triu(Hy, 1)';

end
