% this script outputs the undistorted image of the mousepad in the paper.

% load the figure and resize it
I_unt = imread('PS37_illum.jpg');
I_unt = rgb2gray(I_unt);
I_unt = imresize(I_unt,.04);
I = im2double(I_unt);
[row,col] = size(I);


p = row;
q = col;
L = p*q;

% use the Identity matrix as the dictionary for the mouspad 
B = speye(L,L);

% use 400 columns of DCT matrix as the dictionary for the distortion
N = 400;
C = ones(L,N);
temp = dctmtx(p);

ind = randperm(L);
ind = ind(1:N-1);

for i = 1:N-1
    [a, b] = quorem(sym(ind(i)-1),sym(p));
    C(:,i+1) = kron(temp(:,a+1),ones(p,1)).*repmat(temp(:,b+1),p,1);
end

% normalize the dictionaries so that they are of equal weight
C = C*norm(B,'fro')/norm(C,'fro');

% parameters for the TVL1_BH program
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

% output the clean image
imshow(h)



