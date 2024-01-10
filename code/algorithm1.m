function [Ar,Br,Cr,sigma] = algorithm1(A,B,C,r,phi,init,maxiter)
% [Ar,Br,Cr] = IRKAcom(A,B,C,r,phi,init)
% Computes the local optimal reduced order model in the H2(\barA^c) norm
% 
% INPUTS
% A,B,C = full order model system matrices
% r = reduced order 
% phi = function for the interpolation points
% init = initial guess for the interpolation points
% maxiter = maximum number of iterations
% 
% OUTPUTS
% Ar,Br,Cr = reduced order model system matrices
% sigma = final interpolation points


n = size(A,1);
clear V W
V = zeros(n,r);
W = zeros(n,r);
iter = 0;
conv_crit = inf;
conv_tol = 1e-6;
s = init;
while(conv_crit > conv_tol && iter < maxiter)
    iter = iter+1;   
    
    for i = 1:r
        V(:,i) = (s(i)*speye(n)-A)\B;
        W(:,i) = ((s(i)*speye(n)-A)')\C';
    end
    [V,~] = qr(V,0);
    [W,~] = qr(W,0);
    Ar = (W'*V)\(W'*A*V);
    Br = (W'*V)\(W'*B);
    Cr = C*V;   

    s_old = s;

    s = eig(Ar);
    for i=1:r
        snew(i) = phi(s(i));
    end
    s = snew.';
    conv_crit = norm(sort(s)-sort(s_old))/norm(s_old);
    fprintf('Iteration %d - Convergence %f \n', iter, conv_crit);
end
sigma = s;

end