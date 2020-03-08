close all
new = 1;
p = 256/2;
q = 256/2;
L = p*q;
N = 50;
B = speye(L,L);

g = linspace(-9,4,L)';
C = ones(L,N);
rng(1);
for i = 2:N
   C(:,i) = besselj(g/(6+.1*abs(randn))+5*abs(randn),.1+10*abs(randn)); 
end
C(:,1) = max(C(:,2))*C(:,1);
B = B*norm(C,'fro')/norm(B,'fro');

I_unt = imread('rice.png');
I = im2double(I_unt);
I = I';
delta = 256/2;

params.lambda = 1000;
params.rho = 1e-7;
params.maxIter = 400000;
params.gamma = 1;
signal_size.row = p;
signal_size.col = q;
h = zeros(256,256);

for k = 1:1
    for l = 1:1
    
        y_til = I(delta*(k-1)+1:delta*k,delta*(l-1)+1:delta*l);
        y = y_til(:);
        t = sign(y);
        [h_hat,m_hat,xi]=TV2DBH_ADMM(B, C, y, t, params, signal_size);
        what = B*h_hat;
        xhat = C*m_hat;
        yhat = xhat.*what;
        h(delta*(k-1)+1:delta*k,delta*(l-1)+1:delta*l)=reshape(what,p,q);
        fprintf('error in y at (%d,%d) is %f\n',k,l,norm(y-B*h_hat.*C*m_hat)/sqrt(size(B,1)))
    end
end
