function [maxLambda,x,y] = greatest_nonnegative_eig(AyAx,Ax,maxlam)
%% GREATEST_NONNEGATIVE_EIG biggest eig of AyAx constrained
% AyAx is pxp and Ax is qxp with p<=q

% looks for the greatest Lambda such that
%     -Lambda < maxlam,  Lambda < 0
%     Lambda^2 is an eigenvalue of AyAx with eigenvector x>=0 and nonzero
%     y = Ax x / Lambda <= 0 
% Takes the eigendecomposition of AyAx and checks the positive eigenvalues in descending order
% If the right eigenspace of Lambda has dimension 1, checks that the eigenvector w is either 
%  nonnegative or nonpositive, and in case sets x = |w| and then checks if y = Ax x / Lambda >= 0
% If the right eigenspace of Lambda has dimension >1, takes a basis of the eigenspace U, 
%  computes a slim QR decomposition of [Ax U/Lambda; U] = VR, solves 
%       min ||(VV^T - I)z|| = 0 : z >= 0, e^T z=1
%  and then sets z = [y; x].

% If x,y have been set, then returns (Lambda,x,y). Otherwise go to the next positive eigenvalue of
% AyAx until they finish. 

% If x,y have never been set, returns Lambda = maxlam and random x,y

tol = 1e-10; maxLambda = maxlam;
q = size(Ax,1); p = size(Ax,2); x = zeros(p,1); y = zeros(q,1);
if norm(AyAx) < maxlam^2+tol, return; end 
% if the norm is less then maxlam^2, then it cannot improve the bound
[U,D] = eig(AyAx); % AyAx U = U D
lam_list = diag(D); assert(norm(imag(lam_list)) < tol); lam_list = real(lam_list); 
[lam_list, Ind_sort] = sort(lam_list,"descend"); 
lam_list(lam_list <= maxlam^2 + tol) = []; % remove the eig that cannot improve the bound
lam_list = - sqrt(lam_list); % take the square root and change sign
s = length(lam_list); U = U(:,Ind_sort); U = U(:,1:s); % remove the correspective eigenvectors

eig_ind = 1;
while eig_ind <= s
    Lambda = lam_list(eig_ind);
    trun_ds = lam_list(lam_list < Lambda+tol); eind_end = length(trun_ds);
    % Lambda has multiplicity  1 + eind_end - eig_ind
    U_lam = U(:,eig_ind:eind_end); V_lam = Ax*U_lam/Lambda;
    assert(norm(U_lam'*U_lam)>tol); % U_lam should be full rank, but just in case
    W_lam = [V_lam;U_lam];
    if eind_end == eig_ind % if the dimension of the eigenspace is one
         if ~any(diff(sign(W_lam(W_lam~=0)))) % it is enough to check if W_lam has sign changes
             maxLambda = Lambda; x = abs(U_lam); y = abs(V_lam); return;
             % if there are no sign changes, we have found Lambda, x and y
         end
    else % if the dimension is greater than one
        V = orth(W_lam,tol);
        % try
        %     [V,~] = qr(W_lam,"econ"); % slim QR
        % catch err
        %     [V,~] = qr(W_lam,0); % "old" versions of MATLAB 
        % end
        C = V*V' - eye(p+q); d = zeros(p+q,1); lb = zeros(p+q,1); ub = ones(p+q,1); 
        Aeq = ones(1,p+q); beq = 1; Att = []; b = []; 
        options = optimoptions(@lsqlin,'Display','none');
        assert(norm(imag(C)) < tol); C = real(C); 
        [z,resnorm] = lsqlin(C,d,Att,b,Aeq,beq,lb,ub,zeros(p+q,1),options);
        % solves resnorm = min ||Cz|| : 1 >= z >= 0, e^T z = 1
        if resnorm < tol % if Cz = 0, we have found Lambda, x and y
            maxLambda = Lambda; x = z(q+1:q+p); y = z(1:q); return;
        end
    end
    eig_ind = eind_end+1; % go to the next eigenvalue
end
end