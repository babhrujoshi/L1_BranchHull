function [h, m, xi] = L1BH_ADMM(B, C, y, t, params)
% function [h, m, xi] = L1BH_ADMM(B, C, y, t, params)
% A function that performs the proposed ADMM scheme with L_1 penalties on
% both h and m. In this program:
%
% inputs:
% B: is an L-by-K matrix, representative of the w subspace
% C: is an L-by-N matrix, representative of the x subspace
% y: is an L-by-1 vector of the observations
% t: is an L-by-1 vector of the known signs
% params.lambda: the regularization parameter for noisy case
% params.maxIter: maximum number of ADMM iterations
% params.rho: the ADMM penalty parameter (propoer choice can significantly improve
% the speed)
%
% outputs:
% h: a vector of length K
% m: a vector of length N
% xi: a vector of length L, representative of the noise

global L K N lamb Bb Cc;

[L,K] = size(B);
[Ltemp, N] = size(C);
lamb = params.lambda;
Bb = B;
Cc = C;
rho = params.rho;

if (Ltemp ~= L)
    error('Matrices B and C should have similar number of rows!');
end

% To address the linear solve (E'E + I)\d we use the Cholesky factorization
BB = B'*B + speye(K,K);
CC = C'*C + speye(N,N);

global Bchol Cchol;

Bchol = sparse(chol(BB));
Cchol = sparse(chol(CC));

% initializing all the parameters with some random vectors
u = zeros(3*L,1);
v = zeros(K+N+L,1);
z = v;
beta = v;
alpha = u;

h = zeros(K,1);
m = zeros(N,1);
i = 1;
succ_error = 1;
tol = 1e-7;
while i <= params.maxIter && succ_error>tol
    
    arg = applyE(z) - alpha;
    [term1, term2, term3] = projC(arg(1:L), arg((L+1):(2*L)), arg((2*L+1):(3*L)), y, t);
    u = [term1;term2;term3];
    v = softTh(z - beta, 1/rho);
    z = E_LinSolve( applyEt(alpha+u) + beta + v );
    alpha = alpha + u -applyE(z);
    beta = beta + v - z;
    m_old = m;
    h_old = h;
    
    m = z(1:N);
    h = z((N+1):(N+K));
    succ_error = sqrt(norm(h-h_old,2)^2+norm(m-m_old)^2)/sqrt(K+N);
%     if round(i/1000) == i/1000
%         fprintf('ADMM iterations performed = %d, successive error = %e\n',i,succ_error);
%     end
    i = i+1;
end
m = z(1:N);
h = z((N+1):(N+K));
xi = z((N+K+1):end)/lamb;

%% Bunch of functions listed here
    function xs = softTh(x,rho)
        % function xs = softTh(x,rho)
        % this function performs the soft-thresholding operator
        xs = sign(x).*max(abs(x) - 1/rho, 0);
    end

    function sol = E_LinSolve(d)
        % function sol = E_LinSolve(d)
        % since E is a large matrix and has a nice diagonal structure,
        % instead of explicitly forming E and calculating inv(E'E + I)*d
        % at each iteration, we apply it in a blockwise manner
        solC = Cchol\(Cchol'\d(1:N,1));
        solB = Bchol\(Bchol'\d((N+1):(K+N),1));
        solLamb = d((K+N+1):end,1)/(1+lamb^(-2));
        sol = [solC;solB;solLamb];     
    end

    function u = applyE(v)
        % function u = applyE(v)
        % since E is a large matrix and has a nice diagonal structure,
        % instead of forming E and calculating E*v, we use this function to
        % apply E block by block
        uC = Cc*v(1:N,1);
        uB = Bb*v((N+1):(K+N),1);
        uLamb = v((K+N+1):end,1)/lamb;
        u = [uC;uB;uLamb];
    end

    function v = applyEt(u)
        % function v = applyEt(u)
        % since E is a large matrix and has a nice diagonal structure,
        % instead of forming E and calculating E'*u, we use this function 
        % to apply E' block by block
        vC = Cc'*u(1:L,1);
        vB = Bb'*u((L+1):(2*L),1);
        vLamb = u((2*L+1):end)/lamb;
        v = [vC;vB;vLamb];
    end

end




