function mo_init = thres_init(M, y, K, N , S1, S2)
% A function that computes the initialization for Sparse power factorization with Hard
% Thresholding Pursuit from the paper "Near Optimal Compressed Sensing
% of a Class of Sparse Low-Rank Matrices via Sparse Power Factorization" by
% Lee, Wu, and Bressler. This is algorithm (3) in the paper.
%
% input
% M = is an array where each entry is the rank-1 measurement matrix
% y = is the observation where y_i = <M{i}, ho*mo'> 
% K = is the dimension of ho
% N = is the dimension o mo
% S1 = is the sparsity level of ho
% S2 = is the sparsity level of mo
% 
% output
% mo_init = is the initial estimate of mo used to start SPF-HTP


    M_adj = zeros(size(M{1}));
    
    for i  = 1:length(M)
        M_adj = y(i)*M{i}+M_adj;
    end
    
    [~, sortIndM_adj] = sort(abs(M_adj)',"descend");
    xi = zeros(K,1);
    for j = 1:K
        xi(j) = norm(M_adj(j,sortIndM_adj(1:S2)));
    end
    
    [~, indxi] = sort(xi, "descend");
    J1 = indxi(1:S1);
    
    ProJ1 = zeros(K,K);
    for i = 1:S1
        ProJ1(J1(i),J1(i)) = 1;
    end
    M_adj = ProJ1*M_adj;
    
    [~,indJ2] = sort(vecnorm(M_adj,2,1), "descend");
    J2 = indJ2(1:S2);
    
    ProJ2 = zeros(N,N);
    for i = 1:S2
        ProJ2(J2(i),J2(i)) = 1;
    end
    
    M_adj = M_adj*ProJ2;
    
    [~,~,mo_init] = svd(M_adj);
    mo_init = mo_init(:,1);
    
end