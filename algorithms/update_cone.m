% Solve min_{u in P, ||u||=1} <u, Av>
% 
% coneP is the type of cone: 
% 1) 'generator': P = { u | u = options.G*x, x >= 0}
%               In that case, options.G provides the generators
% 2) 'nonnegort': P = { x | x >= 0}
% 3) 'semidefin': P = { x | mat(x) is PSD}
% 4) 'facetsrep': P = { x | G*x >= g}

function u = update_cone(Av,coneP,G,g);

if coneP == 'generator'
    % Solve max u^T*Av such that ||u||=1, u=G*x, x>= 0
    u = optimize_cone_generators(G,Av); 
elseif coneP == 'facetsrep'
    % Solve max u^T*Av such that ||u||=1, G*u >= 0 
    u = optimize_cone_facets(G,g,Av);     
elseif coneP == 'nonnegort'
    u = projectNonnegOrthnorm1(Av); 
elseif coneP == 'semidefin'
    n = length(Av); 
    u = vec(projectPSDnorm1(reshape(Av,sqrt(n),sqrt(n))));  
end
end 



