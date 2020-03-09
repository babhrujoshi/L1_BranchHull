% this script outputs the undistorted image of the rice grain in the paper.

% load the distorted rice grain image
I_unt = imread('rice.png');
I = im2double(I_unt);
I = I';

p = 256/2;
q = 256/2;
L = p*q;
N = 50;

% use the Identity matrix as the dictionary for the clean rice grain image
B = speye(L,L);

% use the bessel function as the dictionary for the distortion image
g = linspace(-9,4,L)';
C = ones(L,N);
rng(1);
for i = 2:N
   C(:,i) = besselj(g/(6+.1*abs(randn))+5*abs(randn),.1+10*abs(randn)); 
end
C(:,1) = max(C(:,2))*C(:,1);
B = B*norm(C,'fro')/norm(B,'fro');

% parameters
params.lambda = 1000;
params.rho = 1e-5;
params.maxIter = 10000;%400000;
params.gamma = 1;
signal_size.row = p;
signal_size.col = q;

% sign of the unknown signal
y_til = I(1:p,1:q);
y = y_til(:);
t = sign(y);


[h_hat,m_hat,xi]=TV2DBH_ADMM(B, C, y, t, params, signal_size);
what = B*h_hat;
xhat = C*m_hat;
yhat = xhat.*what;
h=reshape(what,p,q);

% recovery error
fprintf('error in recovery is %f\n',norm(y-B*h_hat.*C*m_hat)/sqrt(size(B,1)))
imshow(h)


