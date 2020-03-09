% Test script for solving sparse rank-1 matrix recovery from its linear
% measurements using algorithms (1), (2), and (3) in the paper "Near Optimal Compressed Sensing
% of a Class of Sparse Low-Rank Matrices via Sparse Power Factorization" by
% Lee, Wu, and Bressler. 

L = 1000; %number of measurements
K = 140; %size of ho
N = 100; % size of mo

B = randn(L,K); % dictionary for ho
C = randn(L,N); % dictionary for mo

ho = zeros(K,1);
mo = zeros(N,1);
Sh = 0.05*K;
Sm = 0.05*N;

Sm_ind = randperm(N);
Sh_ind = randperm(K);

ho(Sh_ind(1:Sh)) = ones(Sh,1);
mo(Sm_ind(1:Sm)) = ones(Sm,1);

% normalizing the signal produces better result

ho = ho/norm(ho);
mo = mo/norm(mo);

% noise level
alpha = 0.01;

xi = randn(L,1);
xi = alpha*xi/norm(xi,2);

% noisy observation
y = (B*ho).*(C*mo).*(ones(L,1) + xi );

% gather the rank-1 measurement matrix
M  = {};
for i = 1:L
    M{i} = B(i,:)'*C(i,:);
end

% compute the initialzation for SPF-HTP
mo_init = thres_init(M, y, K, N , Sh, Sm);

[hk,mk] = SPF_HTP(M, y, K, N, Sh, Sm, mo_init);

fprintf("relative error = %f\n",norm(hk*mk'-ho*mo','fro')/norm(ho*mo','fro'))
