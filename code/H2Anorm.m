function [abserror,relerror] = H2Anorm(A,B,C,Ar,Br,Cr,psicom,dpsi)
% out = H2Abarquad(A,B,C,Ar,Br,Cr,psi,dpsi)
% computes the H2Abar error norm 
n = length(A);
r = length(Ar);
fom = @(z) (C*((psicom(1i.*z)*eye(n)-A)\B)) .* sqrt(dpsi(1i.*z));
rom = @(z) (Cr*((psicom(1i.*z)*eye(r)-Ar)\Br)) .* sqrt(dpsi(1i.*z));
funerror = @(z) (fom(z)-rom(z))*conj(fom(z)-rom(z));
funfom = @(z) (fom(z)*conj(fom(z)));
funrom = @(z) (rom(z)*conj(rom(z)));
abserror = sqrt((1/(2*pi))*integral(funerror,-Inf,Inf,'ArrayValued',true));%,'RelTol',1e-6,'AbsTol',1e-10)));
relerror = abserror/sqrt((1/(2*pi))*integral(funfom,-Inf,Inf,'ArrayValued',true));%,'RelTol',1e-6,'AbsTol',1e-10)));
% fomnorm = sqrt((1/(2*pi))*integral(funfom,-Inf,Inf,'ArrayValued',true));
% romnorm = sqrt((1/(2*pi))*integral(funrom,-Inf,Inf,'ArrayValued',true));


end