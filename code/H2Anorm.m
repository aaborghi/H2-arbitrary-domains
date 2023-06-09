function [abserror,relerror] = H2Anorm(A,B,C,Ar,Br,Cr,psicom,dpsi)
% [abserror,relerror] = H2Anorm(A,B,C,Ar,Br,Cr,psicom,dpsi)
% computes the H2Abar error norm 
% 
% INPUTS
% A,B,C = full order model system matrices
% Ar,Br,Cr = reduced order model system matrices
% psicom = conformal map
% dpsi = derivative of the conformal map
% 
% OUTPUTS
% abserror = absolute H2(\bar A^c) error norm
% relerror = relative H2(\bar A^c) error norm

n = length(A);
r = length(Ar);
fom = @(z) (C*((psicom(1i.*z)*eye(n)-A)\B)) .* sqrt(dpsi(1i.*z));
rom = @(z) (Cr*((psicom(1i.*z)*eye(r)-Ar)\Br)) .* sqrt(dpsi(1i.*z));
funerror = @(z) (fom(z)-rom(z))*conj(fom(z)-rom(z));
funfom = @(z) (fom(z)*conj(fom(z)));
funrom = @(z) (rom(z)*conj(rom(z)));
abserror = sqrt((1/(2*pi))*integral(funerror,-Inf,Inf,'ArrayValued',true));
relerror = abserror/sqrt((1/(2*pi))*integral(funfom,-Inf,Inf,'ArrayValued',true));
end