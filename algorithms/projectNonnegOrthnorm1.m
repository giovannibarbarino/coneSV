% Project on norm1 nonnegative orthant
% max_{u>=0,||u||=1} u^T Av

function u = projectNonnegOrthnorm1(Av);

u = max(0, -Av);
if norm(u) <= 1e-9 
    [~,b] = max(-Av);
    u(b) = 1;
else
    u = u/norm(u);
end