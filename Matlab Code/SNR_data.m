
L = 150;
K = 140;
N = 140;

B = randn(L,K);
C = randn(L,N);

M  = {};
for i = 1:L
    M{i} = B(i,:)'*C(i,:);
end

ho = zeros(K,1);
mo = zeros(N,1);
Sh = 0.05*K;
Sm = 0.05*N;

Sm_ind = randperm(N);
Sh_ind = randperm(K);

ho(Sh_ind(1:Sh)) = ones(Sh,1);
mo(Sm_ind(1:Sm)) = ones(Sm,1);



alpha_list = 0:.1:.8;

rel_error_mean = alpha_list;
rel_error_std = alpha_list;

max_sup_noise = alpha_list;

trials = 2;


params.lambda = 0;
params.maxIter = 10000;
params.rho = 1;



for i = 1:length(alpha_list)
    alpha= alpha_list(i);
    rel_error_list = [];
    max_sup_noise_list = [];
    
    for j = 1:trials
        xi = randn(L,1);
        xi = alpha*xi/norm(xi,2);
        
        max_sup_noise_list = horzcat(max_sup_noise_list, max(abs(xi)));

        y = (B*ho).*(C*mo).*(ones(L,1) + xi );
        

        t= sign(B*ho);

        c = sqrt(norm(mo,1)/norm(ho,1));
        ho_hat = ho*c;
        mo_hat = mo/c;
        [h_rec,m_rec] = SPF_HTP(M, y, K, N, Sh, Sm);

%         [h_rec,m_rec] = L1BH_ADMM_noslack(B, C, y, t, ho_hat, mo_hat, params);

%         rel_error_list = horzcat(rel_error_list, norm([h_rec;m_rec]-[ho_hat;mo_hat])/norm([ho_hat;mo_hat]));  
        rel_error_list = horzcat(rel_error_list, norm(h_rec*m_rec'-ho_hat*mo_hat','fro')/norm(ho_hat*mo_hat','fro'));  
    end
    rel_error_mean(i) = mean(rel_error_list);
    rel_error_std(i) = std(rel_error_list);
    max_sup_noise(i) = max(max_sup_noise_list);
    fprintf("alpha = %f, relative error (mean, std) = (%f,%f)\n", alpha,rel_error_mean(i),rel_error_std(i))
    
end

result = [alpha_list; rel_error_mean; rel_error_std; max_sup_noise];
% save('L1_BH_SNR.mat', 'result')
