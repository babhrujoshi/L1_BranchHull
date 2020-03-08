close all
clear all

L = 100;
K = 140;
N = 140;

B = randn(L,K);
C = randn(L,N);

ho = zeros(K,1);
mo = zeros(N,1);
Sh = 0.05*K;
Sm = 0.05*N;

Sm_ind = randperm(N);
Sh_ind = randperm(K);

ho(Sh_ind(1:Sh)) = ones(Sh,1);
mo(Sm_ind(1:Sm)) = ones(Sm,1);

alpha = 0;

xi = randn(L,1);
xi = alpha*xi/norm(xi,2);

y = (B*ho).*(C*mo).*(ones(L,1) + xi );

M  = {};
for i = 1:L
    M{i} = B(i,:)'*C(i,:);
end

[hk,mk] = SPF_HTP(M, y, K, N, Sh, Sm, mo+randn(N,1));



function Fy = rightF(M, y)
    [N,~] = size(M{1});
    Fy = zeros(length(M),N);
    for i = 1:length(M)
        Fy(i,:) = y'*M{i}';
    end
end

function Gx = leftG(M,x)
    [~,K] = size(M{1});
    Gx = zeros(length(M),K);
    for i = 1:length(M)
        Gx(i,:) = x'*M{i};
    end
end

function xk = HTP(M, b, s)
    [~,P] = size(M);  
    tol = 1e-7;
    max_iter = 1000;
    k = 1;
    step = 0.01;
    succ_err = 1;
    maxsInd = 1:s;
    xk = zeros(P,1);
    
    while k<max_iter && succ_err >tol
        xk_maxs = xk(maxsInd);
        M_maxs = M(:,maxsInd);
        xk = xk+step*M'*(b-M_maxs*xk_maxs);
        [~,sortInd] = sort(xk, 'descend');
        maxsInd = sortInd(1:s);
        x_next = zeros(P,1);
        x_next(maxsInd) = M(:,maxsInd)\b;
        succ_err = norm(xk - x_next,2)/norm(x_next);
        xk = x_next;
%         fprintf("iteration = %d, error = %f\n",k, succ_err)
        k = k+1;
    end
end


function [hk,mk] = SPF_HTP(M, y, K, N, S1, S2, mo_init)
    k = 1;
    tol = 1e-7;
    max_iter = 1000;
    mk = mo_init;
    hk = zeros(K,1);
    succ_error = 1;
    while k < max_iter && succ_error > tol
        mk_1 = mk;
        hk_1 = hk;
        
        mk = mk/norm(mk);
        Fm = rightF(M,mk);
        
        if S1 < K
            hk = HTP(Fm, y, S1);
        else
            hk = Fm \ y;
        end
        
        hk = hk/norm(hk);
        Gh = leftG(M,hk);
        
        if S2 < N
            mk = HTP(Gh, y, S2);
        else
            mk = Gh \ y;
        end
        
        succ_error = norm([hk_1;mk_1] - [hk;mk])/norm([hk;mk]);
        fprintf("iteration = %d, succ error = %f\n", k, succ_error)                  
        k = k+1;
    end
        
end





