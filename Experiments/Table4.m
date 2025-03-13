% Largest angle between circulant positive definite matrices
%  and circulant nonnegative symmetric matrices computed by
%   Gurobi for dinemsions from 5 to 23

% The angle, expressed as a fraction of pi, is stored in the
%  array Gur, where Gur(n) is referred to dimension n+4

% WARNING: for high dimensions n Gurobi may take more than
% 24h, so a timelimit is set at 24h

% include algorithms
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename)); cd('..'); 
addpath(genpath('./'));

timelimit = 24*3600;
Gur = zeros(1,19);

for n = 5:23
    N = n
    if mod(n,2) == 0
        m = N/2; Aux = fft(eye(N))/sqrt(N); Aux = real(Aux); 
        A = 2*Aux(2:m+1,2:m+1);
        A(m,:) = A(m,:)/sqrt(2); 
        A(:,m) = A(:,m)/sqrt(2);
    else
        m = (N-1)/2; Aux = fft(eye(N))/sqrt(N); 
        Aux = real(Aux); A = 2*Aux(2:m+1,2:m+1);
    end
    G = eye(m); H = G;

    [val,~] = generators_gurobi_uv_test(G,H,A,timelimit);
    Gur(n-4) = acos(val)/pi; 

end