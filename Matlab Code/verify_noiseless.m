% this script generates tensor M of size (trials, L, N) where 
% trials is the number of trials for each (L,N). This checks the recovery
% performace in the noiseless case.

for L = 4:4:140
    for N = 20:20:300
        K = N;
        for i = 1:10    
            % generate data
            B = randn(L,N)/sqrt(L);
            C = randn(L,K)/sqrt(L);
            
            Sm = .05*N;
            mo = zeros(N,1);
            permm = randperm(N);
            mo(permm(1:(Sm))) = ones(Sm,1);
            t = randn(N,1);
            t(t<0) = -1; t(t>0) = 1;
            mo = mo.*t;

            Sh = .05*K;
            ho = zeros(K,1);
            permh = randperm(K);
            ho(permh(1:(Sh))) = ones(Sh,1);
            t = randn(K,1);
            t(t<0) = -1; t(t>0) = 1;
            ho = ho.*t;
            
            xo = C*mo;
            wo = B*ho;
            t = sign(wo);
            yo = (xo).*(wo); 
            y= yo;
            
            scale = sqrt(norm(ho,1)/norm(mo,1));
            h_hat = ho/scale;
            m_hat = mo*scale;
            
            
            params.lambda = 1;
            params.rho = 1;
            params.maxIter = 100000;
            
            [h,m] = L1BH_ADMM_noslack(B,C,y,t,h_hat,m_hat,params);
            
            error = (norm(h-h_hat)^2 + norm(m-m_hat)^2)^(1/2);
            M(i,cast(35-(L-4)/4,'int8'),cast((N-20)/20+1,'int8')) = error;
        end
        mean_error=sum(M(:,cast(35-(L-4)/4,'int8'),cast((N-20)/20+1,'int8')))/10;
        fprintf('total success for L = %d, N =K = %d and Sm = Sh = %d is %d\n',L,K,Sm,mean_error); 
    end
end