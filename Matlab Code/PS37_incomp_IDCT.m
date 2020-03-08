
% read the distorted image and resize it for processing
I_unt = imread('PS37_illum.jpg');
I_unt = rgb2gray(I_unt);
I_unt = imresize(I_unt,.04);
I = im2double(I_unt);
[row,col] = size(I);

% 
p = row;
q = col;
L = p*q;
B = speye(L,L);

N1 = 200;
N2 = 200;
N = 400;
C = ones(L,N);
temp = dctmtx(p);

ind = randperm(L);
ind = ind(1:N-1);

for i = 1:N-1
    [a, b] = quorem(sym(ind(i)-1),sym(p));
    C(:,i+1) = kron(temp(:,a+1),ones(p,1)).*repmat(temp(:,b+1),p,1);
end

% g = linspace(-9,4,L)';
% C1 = ones(L,N1);
% rng(1);
% for i = 2:N1
%    C1(:,i) = besselj(g/(6+.1*abs(randn))+5*abs(randn),.1+10*abs(randn)); 
% end
% C2 = ones(L,N1);
% rng(1);
% for i = 2:N1
%    C2(:,i) = besselj(g/(6+.1*abs(randn))+5*abs(randn),.1+10*abs(randn)); 
% end
% C2_temp = zeros(p*q,N1); 
% for i = 0:p*q-1
%         C2_temp(floor(i/p)+1+mod(i,p)*p,:) = C2(i+1,:);
% end
% keyboard
% C2 = C2_temp;
% 
% C = [C1,C2];
% C(:,1) = max(C(:,2))*C(:,1);

% C(:,1) = max(C(:,2))*C(:,1);


C = C*norm(B,'fro')/norm(C,'fro');

params.lambda = 1000;
params.rho = 7e-4;
params.maxIter = 10000;%64000;
params.gamma = 1;
signal_size.row = p;
signal_size.col = q;
    
y_til = I;
y = y_til(:);
t = sign(y);

[h_hat,m_hat,xi]=TV2DBH_ADMM(B, C, y, t, params, signal_size);

what = B*h_hat;
xhat = C*m_hat;
yhat = xhat.*what;
h =reshape(what,p,q);
imshow(h)