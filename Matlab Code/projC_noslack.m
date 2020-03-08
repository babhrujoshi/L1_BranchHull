function [x, w] = projC_noslack(xp, wp, y, t)
% function [x, w, xi] = projC(xp, wp, y, t)
% projects a vector in the form of (xp,wp) in R^{2L} onto the set C
% defined in the paper. In this function all inputs should have identical
% sizes of L-by-1

L = size(xp,1);
if L == 1
    error('Make sure all the inputs are vertical vectors of the same size!');
end

% making sure all inputs are of the same size
AllSizes = [size(xp,1), size(wp,1), size(y,1), size(t,1)];
if (norm(diff(AllSizes))>0)
    error('Make sure all the inputs are vertical vectors of the same size!');
end

% Our assumptions are based on having nonzero entries for y, so if an entry
% of y is zero we replace it with a very small number
y(y==0) = 1e-12;

mac = 1e-13;

% setting s to be the sign of y
s = sign(y);

% checking if a point meets any of the inequality constraints
C1 = abs(y) - s.*(xp.*wp);
C2 = t.*wp;

inSet = (C1 <= mac)&(C2>=-mac);

% setting the default projection to the point itself
x  = xp;
w  = wp;

% forming a matrix of L-by-5 as the coefficient matrix for all L
% polynomials

% pCoef = [4*ones(L,1), -4*wp, xp.^2, zeros(L,1), -(y.^2)];
pCoef = [ones(L,1), -wp, zeros(L,1), y.*xp, -(y.^2)];


% lopping to find the projections one after the other

for i = 1 : L
    % skipping the points with trivial projection (projection of the points
    % that are already in the set is trivial)
    if (inSet(i)>0)
        continue;
    end
    pRoots = roots(pCoef(i,:));
    pRoots = pRoots(abs(imag(pRoots))<1e-15);
    
    isFeasible = ((t(i)*pRoots)>=-mac)&( (abs(y(i)) - s(i)*(x(i))*pRoots) >=-mac );
    % if none of the roots meet the required conditions we let the
    % projection be the point itself
    if sum(isFeasible) == 0
        continue;
    end
    % if for some numerical error multiple roots meet the condition we pick
    % one of them
    if sum(isFeasible) > 1
        I = find(isFeasible);
        isFeasible(I(2:end)) = 0;
    end
    pRoots = pRoots(isFeasible);
    w(i) = pRoots;
%     mu = (abs(y(i)) - s(i)*(x(i))*pRoots)/(2*pRoots^2);
    mu = (abs(y(i)) - s(i)*(x(i))*pRoots)/(pRoots^2);

    x(i) = x(i) + mu*s(i)*pRoots;
end




