% [u,v,e] = AlternatingOptimization(A,options)
%
% Alternating optimization for solving
%
%           min_{u in P, v in Q, ||u||=||v||=1} <u, Av>    (1)
%
% where
%       P = {u | u = G*x, x >= 0}, or
%         = {u | G*u >= g}, or
%         = {x | x >= 0}, (defaut) or
%         = {x | mat(x) is positive semidefinite (PSD)}.
%
% *** Input ***
%   A              : m-by-n matrix
%
% options:
%   maxiter        : maximum number of iterations (default = 500)
%   accuracy       : continue when  [ ||uk-1 - uk|| >= accuracy or
%                                     ||vk-1 - vk|| <= accuracy ]
%                               and objetive(k-1)-objective(k) >= accuracy
%                    (default = 1e-4)
%   cone.P         : type of cones, we have implemented 3 cases:
%                    1) 'generator': P = { u | u = options.G*x, x >= 0}
%                       In that case, options.G provides the generators
%                    2) 'nonnegort': P = { u | u >= 0}
%                    3) 'semidefin': P = { u | mat(u) is PSD}
%                    4) 'facetsrep': P = { u | options.G*u >= options.g}
%   cone.Q         : same structure as .P
%   G, H           : When options.cone.P = 'generator' or 'facetsrep'
%                     G is required, and similarly for Q and H.
%                    In the case 'generator':
%                    The columns of G (resp. H) are the rays of P (resp. Q)
%                    hence G (m-by-k) u = G*x in P for some x >= 0.
%                      and H (n-by-p) v = H*y in Q for some y >= 0.
%   g,h              In the case 'facetsrep':
%                    G provides the inequalities for P = {u | G*u >= g}
%                    H provides the inequalities for Q = {v | H*v >= h}
%                    default: g=0, h=0
%   v0             : initializations for the algorithm.
%                     default: v <-- argmin_{v in Q} (u^TA) v
%                               where u is randn(m,1).
%
% *** Output ***
%   (u,v) in PxQ   : approximate solution to Problem (1)
%    e             : evolution of u'*A*v

function [u,v,e] = AlternatingOptimization(A,options)

curtime = tic;
[m,n] = size(A);
if nargin <= 1
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 500;
end
if ~isfield(options,'accuracy')
    options.accuracy = 1e-6;
end
% Representation for P and Q
if ~isfield(options,'cone')
    warning('No cones specified: nonnegative orthant is used');
    if ~isfield(options.cone,'P')
        options.cone.P = 'nonnegort';
    end
    if ~isfield(options.cone,'Q')
        options.cone.P = 'nonnegort';
    end
end
if ~isfield(options,'G')
    if options.cone.P == 'generator'
        error('No generators specified for P.');
    elseif options.cone.P == 'facetsrep'
        error('No facets specified for P.');
    else
        options.G = [];
    end
end
if ~isfield(options,'H')
    if options.cone.Q == 'generator'
        error('No generators specified for Q.');
    elseif options.cone.Q == 'facetsrep'
        error('No facets specified for Q.');
    else
        options.H = [];
    end
end
if options.cone.P == 'facetsrep'
    if ~isfield(options,'g')
        warning('No right hand side g specified for P: zero is used');
        options.g = zeros(size(options.G,1),1);
    end
else
    options.g = [];
end
if options.cone.Q == 'facetsrep'
    if ~isfield(options,'h')
        warning('No right hand side h specified for Q: zero is used');
        options.h = zeros(size(options.H,1),1);
    end
else
    options.h = [];
end
if ~isfield(options,'display')
    options.display = 0;
end
% Extrapoaltion parameters
if ~isfield(options,'beta')
    options.beta = 0.5;
end
if ~isfield(options,'eta')
    options.eta = 2;
end
if ~isfield(options,'gamma')
    options.gamma = 1.05;
end
% Initialization: random by default
if ~isfield(options,'v0')
    u0 = randn(m,1);
    options.v0 = update_cone(A'*u0,options.cone.Q,options.H,options.h);
end
v0 = options.v0; % initialization
v = 0; ve = v0; vp = 0; % v sequences
u = 0; ue = 0; up = 0; % u sequences
i = 1; % iteration count
beta_p = options.beta; rs = 0; % restart
if options.display == 1
    disp('Evolution of iteration number and the objective:');
    tdelay = 0.05; numprint = 0;
end
while i <= options.maxiter && ...
        (rs == 1 || norm(v-vp) >= options.accuracy || ...
        norm(u-up) >= options.accuracy || ...
        (i <= 3 || abs(e(i-2)-e(i-1)) >= options.accuracy*abs(e(i-2))))
    rs = 0; % set restart to zero
    up = u; vp = v; % previous iterates
    % Update u
    u = update_cone(A*ve,options.cone.P,options.G,options.g);
    % Extrapolate u
    if i >= 2
        ue = u+options.beta*(u-up);
    else
        ue = u;
    end
    % Update v
    Atu = A'*ue;
    v = update_cone(Atu,options.cone.Q,options.H,options.h);
    e(i) = u'*A*v; % error at iteration i
    % Extrapolate v
    if i >= 2
        ve = v+options.beta*(v-vp);
    else
        ve = v;
    end
    % Update extrapolation parameters depending on the error
    if i>=2 && e(i-1) < e(i) && options.beta > 0
        beta_p = options.beta/options.eta;
        options.beta = 0;
        u = up; v = vp; ve = vp; % ue = up; 
        e(i) = e(i-1); rs = 1; % restart is one
    else
        options.beta = min(1,beta_p*options.gamma);
        beta_p = options.beta;
    end
    % Display evolution of the error
    if options.display == 1 && toc >= tdelay*(2^numprint)
        fprintf('%3.0f: %1.4f...',i,e(i));
        numprint = numprint+1;
        if mod(numprint,5) == 0
            fprintf('\n');
        end
    end
    i = i+1;
end
if options.display == 1
    fprintf('\n');
end