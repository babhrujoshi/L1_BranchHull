function [h, m, xi] = TV2DBH_ADMM(B, C, y, t, params, signal_size)
% function [h, m, xi] = TV2DBH_ADMM(B, C, y, t, params,)
% A function that performs the proposed ADMM scheme with TV penalty on
% h and L_1 penalty on m. In this program:
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

global L K N lamb Bb Cc DhB p gamma;
p = signal_size.row;
q = signal_size.col;

[L,K] = size(B);
[Ltemp, N] = size(C);
lamb = params.lambda;
gamma = params.gamma;
Bb = B;
Cc = C;
rho = params.rho;

% difference in vertical direction
Dv = spdiags([-ones(L,1), ones(L,1)],[0 1],L,L);
Dv(L,L)=1;
Dv(L,L-1) = -1;
for i = 1:q
    Dv(i*p,:) = Dv(i*p-1,:);
end

%difference in horizontal direction
Dh = spdiags([-ones(L,1), ones(L,1)],[0,p],L-p,L);
D = [Dh;Dv];



DhB = D*B; %size of DhB is (2L-p) by L

if (Ltemp ~= L)
    error('Matrices B and C should have similar number of rows!');
end

% To address the linear solve (E'E + D'D)\d we use the Cholesky factorization
BB = B'*B + DhB'*DhB + 1E-13*speye(K,K);
CC = gamma^(-2)*(C'*C) + speye(N,N);

global Bchol Cchol;

Bchol = sparse(chol(BB));
Cchol = sparse(chol(CC));

% initializing all the parameters with some random vectors
u = zeros(3*L,1);
v = zeros(N+L+2*L-p,1);
z = zeros(L+K+N,1);
beta = v;
alpha = u;
h = zeros(K,1);

for i = 1 : params.maxIter 
    arg = applyE(z) - alpha;
    [term1, term2, term3] = projC(arg(1:L), arg((L+1):(2*L)), arg((2*L+1):(3*L)), y, t);
    u = [term1;term2;term3];
    v = softTh(applyD(z) - beta, 1/rho);
    z = E_LinSolve( applyEt(alpha+u) + applyDt(beta + v) );
        alpha = alpha + u - applyE(z);
    beta = beta + v - applyD(z);
    
    h_new = z((N+1):(N+K));  
    suc_error = norm(h_new-h)/sqrt(K);
    h = h_new;
    if round(i/1000) == i/1000    
        fprintf('ADMM iterations performed %d with error %e\n',i,suc_error);        
    end
    
end

m = z(1:N)/gamma;
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
        uC = Cc*v(1:N,1)/gamma;
        uB = Bb*v((N+1):(K+N),1);
        uLamb = v((K+N+1):end,1)/lamb;
        u = [uC;uB;uLamb];
    end

    function v = applyEt(u)
        % function v = applyEt(u)
        % since E is a large matrix and has a nice diagonal structure,
        % instead of forming E and calculating E'*u, we use this function
        % to apply E' block by block
        vC = Cc'*u(1:L,1)/gamma;
        vB = Bb'*u((L+1):(2*L),1);
        vLamb = u((2*L+1):end)/lamb;
        v = [vC;vB;vLamb];
    end

    function u = applyD(v)
        % function u = applyE(v)
        % since E is a large matrix and has a nice diagonal structure,
        % instead of forming E and calculating E*v, we use this function to
        % apply E block by block
        uC = v(1:N,1);
        uB = DhB*v((N+1):(K+N),1);
        uLamb = v((K+N+1):end,1);
        u = [uC;uB;uLamb];
    end

    function v = applyDt(u)
        % function v = applyEt(u)
        % since E is a large matrix and has a nice diagonal structure,
        % instead of forming E and calculating E'*u, we use this function
        % to apply E' block by block
        vC = u(1:N,1);
        vB = DhB'*u((N+1):(2*L-p+N),1);
        vLamb = u((2*L-p+N+1):end,1);
        v = [vC;vB;vLamb];
    end

end
