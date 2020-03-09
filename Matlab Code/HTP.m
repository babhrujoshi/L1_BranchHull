function xk = HTP(M, b, s)
% A function that performs hard thresholding pursuit
%
% input:
% M = is the meaurement matrix
% b = is the observation
% s = sparsity level of the unknown signal
%
% output: xk = is the estimate of the sparse signal that solves Mx = b

    [~,P] = size(M);  

    % parameter of the function   
    tol = 1e-7;
    max_iter = 1000;
    k = 1;
    step = 1;
    succ_err = 1;
    
    % initialize the algorithm with zero vector with maxsInd the indicies of
    % non-zero elements
    xk = zeros(P,1);
    maxsInd = 1:s;    
    
    
    while k<max_iter && succ_err >tol
        
        % xk and M restricted to the support 
        xk_maxs = xk(maxsInd);
        M_maxs = M(:,maxsInd);
        
        xk = xk+step*M'*(b-M_maxs*xk_maxs);
        [~,sortInd] = sort(abs(xk), 'descend');
        maxsInd = sortInd(1:s);
        x_next = zeros(P,1);
        
        % perform least squares on the support
        x_next(maxsInd) = M(:,maxsInd)\b;
        
        % compute successive error
        succ_err = norm(xk - x_next,2)/norm(x_next);
        xk = x_next;
%         fprintf("iteration = %d, error = %f\n",k, succ_err)
        k = k+1;
    end
end