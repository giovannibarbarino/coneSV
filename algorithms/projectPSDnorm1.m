% Project the matrix Q onto the PSD cone 
% 
% This requires an eigendecomposition and then setting the negative
% eigenvalues to zero, 
% or all eigenvalues in the interval [epsilon,delta] if specified. 

function Qp = projectPSDnorm1(Q) 
n = size(Q,1); 
if isempty(Q)
    Qp = Q;
    return;
end
Q = (Q+Q')/2; 
if max(max(isnan(Q))) == 1 || max(max(isinf(Q))) == 1
    error('Input matrix has infinite of NaN entries');
end
[V,e] = eig(Q); 
e = projectNonnegOrthnorm1(diag(e)); 
Qp = V * diag(e) * V'; 