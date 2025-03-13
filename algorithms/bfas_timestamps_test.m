function [lamvec,u,v,timest] = bfas_timestamps_test(G,H,A,timelimit)
%% BFAS_TIMESTAMPS_TEST finds the minimum singular values of A with conic constraints
% [lamvec,u,v,timest] = bfas_timestamps_test(G,H,A,timelimit)
%
% Brute Force Active Set method for solving
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
%                    G, H must have positively independent columns.
%                    the columns of G, H must have unit norm. 
%                    It is supposed that the solution to (1) is in 
%                     [ -||A||, 0 ]
%   timelimit      : time limit for the algorithm
%
% *** Output ***
%   (u,v)          : approximate solution to Problem (1)
%   lamvec         : vector of values of u'*A*v for each iteration where 
%                     it improves.
%   timest         : vector of time in seconds corrisponding to 
%                     the entries in lamvec


% The algorithm tries every subset of indices I<[p], J<[q] such that 
% |I| + |J| <= n + m - r such that the respective reduced G, H are 
% full column rank. r is the multiplicity of ||A|| as a singular value 
% of A.


timest = []; lamvec = [];

tic;
flag_overtime = 0;
     
tol = 1e-10;
m = size(A,1); n = size(A,2); p = size(G,2); q = size(H,2);  

gA = G'*A*H; [iu,iv] = find(gA==min(gA(:))); lam = gA(iu(1),iv(1));
lamvec = [lamvec lam]; time = toc; timest = [timest time];
u = G(:,iu(1)); v = H(:,iv(1));
% start from the minimum element i,j in G^TAH and the
% associated couple u=g_i, v=h_j
s = svd(A); nA = s(1); r = length( s(s>nA-tol) ); 
% compute r as the the multiplicity of ||A|| as a singular value of A


for p1 = 1:min(p,m) % |I|, number of face generator in P
    if flag_overtime, break; end
    time = toc; if time > timelimit, flag_overtime = 1; break; end
    for q1 = 1:min(q,n) % |J|, number of face generator in Q
        if flag_overtime, break; end
        time = toc; if time > timelimit, flag_overtime = 1; break; end
        if q1+p1 == 2, continue; end % the best value for 1x1 is min(G^TAH)
        if q1+p1 > n+m-r, break; end % |I| + |J| <= n + m - r
        idxp = nchoosek(1:p,p1); idxq = nchoosek(1:q,q1); 
        % generate all possible I,J with cardinality p1, q1
        for II = idxp', sG = G(:,II); G2 = sG'*sG; 
            if flag_overtime, break; end
            time = toc; if time > timelimit, flag_overtime = 1; break; end
            if det(G2)<tol, continue; end % the restricted G must have full column rank
            for JJ = idxq', sH = H(:,JJ); H2 = sH'*sH;
                time = toc; if time > timelimit, flag_overtime = 1; break; end
                if det(H2)<tol, continue; end % the restricted H must have full column rank
                sA = gA(II,JJ); Ax = H2\sA'; Ay = G2\sA;
                if p1 < q1, Alam = Ay*Ax; Aoth = Ax; s=1; else, Alam = Ax*Ay; Aoth = Ay; s=2; end
                % choose the matrix with least dimension
                if norm(Alam) < (lam-tol)^2, continue; end 
                [sl,sx,sy] = greatest_nonnegative_eig(Alam,Aoth,lam); 
                % check if x,y nonnegative nonzero and lam < 0 exist 
                % such that Ax x = lam y, Ay y = lam x
                if sl > lam-tol, continue; end 
                % if best sv is worse than current best, then skip
                if s==1, x = sx; y = sy; else, y = sx; x = sy; end
                su = sG*x;
                sv = sH*y;
                rho = norm(su);
                assert(rho>tol); assert(abs(rho-norm(sv))<tol);
                lam = sl; u = su/rho; v = sv/rho; 
                lamvec = [lamvec lam]; time = toc; timest = [timest time];
                % otherwise compute u=Gx/||Gx||, v=Hy/||Hy|| 
                % and update the best with the newfound lam
            end
        end
    end
end

end