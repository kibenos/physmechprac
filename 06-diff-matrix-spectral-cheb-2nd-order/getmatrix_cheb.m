function [D, x] = getmatrix_cheb(N)
    if N==0, D=0; x=1; return, end
    
    x = cos(pi*(0:N)/N)'; 
    %t  = (2 * x - (x0 + x1)) / (x1 - x0);
    c  = [2; ones(N-1,1); 2].*(-1).^(0:N)';
    X  = repmat(x,1,N+1);
    dX = X-X';                  
    D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
    D  = D - diag(sum(D'));                 % diagonal entries
end