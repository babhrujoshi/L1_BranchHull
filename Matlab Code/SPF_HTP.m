function [hk,mk] = SPF_HTP(M, y, K, N, S1, S2, mo_init)
% A function that performs Sparse power factorization with Hard
% Thresholding Pursuit from the paper "Near Optimal Compressed Sensing
% of a Class of Sparse Low-Rank Matrices via Sparse Power Factorization" by
% Lee, Wu, and Bressler. This is algorithm (1) in the paper.
%
% input
% M = is an array where each entry is the rank-1 measurement matrix
% y = is the observation where y_i = <M{i}, ho*mo'> 
% K = is the dimension of ho
% N = is the dimension o mo
% S1 = is the sparsity level of ho
% S2 = is the sparsity level of mo
% mo_init = is the initial estimate of mo used to start SPF-HTP
%
% output
% hk = is the estimate of ho
% mk = is the estimate of mo

    % parameter 
    k = 1;
    tol = 1e-14;
    max_iter = 1000;
    succ_error = 1;

    % initialization
    mk = mo_init;
    hk = zeros(K,1);
    
    while k < max_iter && succ_error > tol
        mk_1 = mk;
        hk_1 = hk;
        
        mk = mk/norm(mk);
        
        Fm = rightF(M,mk); % equation (6) in the paper
        
        if S1 < K
            hk = HTP(Fm, y, S1);
        else
            hk = Fm \ y;
        end
        
        hk = hk/norm(hk);
        
        Gh = leftG(M,hk); % equation (6) in the paper
        
        if S2 < N
            mk = HTP(Gh, y, S2);
        else
            mk = Gh \ y;
        end
        
        % compute successive error
        succ_error = norm([hk_1;mk_1] - [hk;mk])/norm([hk;mk]);
        fprintf("iteration = %d, succ error = %f\n", k, succ_error)                  
        k = k+1;
    end        
end

function Fy = rightF(M, y)
% A function that performs the linear operator F in equation (6) of the
% paper
    [N,~] = size(M{1});
    Fy = zeros(length(M),N);
    for i = 1:length(M)
        Fy(i,:) = y'*M{i}';
    end
end

function Gx = leftG(M,x)
% A function that performs the linear operator G in equation (6) of the
% paper
    [~,K] = size(M{1});
    Gx = zeros(length(M),K);
    for i = 1:length(M)
        Gx(i,:) = x'*M{i};
    end
end








