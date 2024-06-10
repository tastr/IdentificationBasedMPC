function [A,B,e,condNumb] = lsq_ext_reg(X0, X1, U,n,m,p, reg_param)
    % p = output dimension, i.e. for state p = n 
    
    [~,N] = size(X0);
    
    Z = [X0; U; ones(1,N)]; 
    

    M = X1*transpose(Z)*pinv(Z*transpose(Z) + reg_param*eye(size(Z*transpose(Z)))); 

    A = M(1:p,1:n); 
    B = M(1:p, n+1:n+m); 
    e = M(1:p, end);
    condNumb = cond(Z*transpose(Z)); 
    
end