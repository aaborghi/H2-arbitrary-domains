clear; clc;
rng(4)

%% Constructing the discretized wave equation
nx = 50;
n = 2*nx;
xa = 0;
xb = 1;
nu = 1;
damping = 0;
hx = (xb-xa)/(nx+1);
xd = xa+hx:hx:xb-hx;
ex = ones(nx,1);
I = speye(nx);
Laplace_x = 1/hx^2*spdiags([ex -2*ex ex], -1:1, nx, nx);
e1x = I(:,1);
enx = I(:,nx);
O = sparse(nx,nx);
A = nu*Laplace_x;
AA = [O,I;A,-damping*I];
B = zeros(nx,1);
C = zeros(1,nx);
for i = 1:nx
    if i*hx >= 0.1 && i*hx <= 0.4
        C(i) = hx;
    end
    if i*hx >= 0.6 && i*hx <= 0.7
        B(i) = 1;
    end
end
BB = [zeros(nx,1);B];
CC = [C,zeros(1,nx)];
II = speye(nx);

%% Running Algorithm conformalBT with ellipse
r = 20;
R = 1+5e-3;
M = 3e2;
c = -1e-2;
r1 = 0.5*(R+inv(R));
r2= 0.5*(R-inv(R));
t = linspace(0,2*pi,500);

psi = @(x) c + 0.5*1i*M.*(R.*(x+1)./(x-1)+(x-1)./(R.*(x+1)));
dpsi = @(x) 1i*M.*(-R./((x-1).^2)+1/(R.*(x+1).^2));
psiinv = @(x) ((-1i.*(x-c)./M+sqrt((-1i.*(x-c)./M).^2-1))./R+1)/((-1i.*(x-c)./M+sqrt((-1i.*(x-c)./M).^2-1))./R-1);
% psi = @(x) c + x;
% dpsi = @(x) 1;
con = @(z) (psi(1i.*z)*eye(n)-AA)\BB * BB'/(psi(1i.*z)*eye(n)-AA)' * sqrt(dpsi(1i.*z)) * sqrt(dpsi(1i.*z))' ;
obs = @(z) (psi(1i.*z)*eye(n)-AA)'\CC' * CC/(psi(1i.*z)*eye(n)-AA) * sqrt(dpsi(1i.*z)) * sqrt(dpsi(1i.*z))' ;


P = integral(con,-Inf,Inf,'ArrayValued',true)/(2*pi);
P = 0.5*(P + P')+1e-12*eye(n);
Q = integral(obs,-Inf,Inf,'ArrayValued',true)/(2*pi);
Q = 0.5*(Q + Q')+1e-12*eye(n);
U = chol(P); U = U'; 
L = chol(Q); L = L';
[Z,S,Y] = svd(L'*U, 'econ');



Z1 = Z(:,1:r);
Y1 = Y(:,1:r);
S1 = S(1:r,1:r);  S1half = sqrt(S1);

Wr = L*Z1/S1half;
Vr = U*Y1/S1half;

Ar = (Wr'*Vr)\(Wr'*AA*Vr);
Br = (Wr'*Vr)\(Wr'*BB);
Cr = CC*Vr;

eigAr = eig(Ar);
eigA = eig(full(AA));
for i = 1:500
    plot(real(c) + M*r2*sin(t(i)),imag(c) + M*r1*cos(t(i)),'.k')
    hold on
end
plot(real(eigA),imag(eigA),'ob'); hold on
plot(real(eig(Ar)),imag(eig(Ar)),'rx');
plot(real(psiinv(eigA)),imag(psiinv(eigA)),'kx');

conr = @(z) (psi(1i.*z)*eye(r)-Ar)\Br * Br'/(psi(1i.*z)*eye(r)-Ar)' * sqrt(dpsi(1i.*z)) * sqrt(dpsi(1i.*z))';
obsr = @(z) (psi(1i.*z)*eye(r)-Ar)'\Cr' * Cr/(psi(1i.*z)*eye(r)-Ar) * sqrt(dpsi(1i.*z)) * sqrt(dpsi(1i.*z))';
% Pr = integral(conr,-Inf,Inf,'ArrayValued',true)/(2*pi);
% Qr = integral(obsr,-Inf,Inf,'ArrayValued',true)/(2*pi);

% % test
% psiinv = @(x) ((-1i*(x-eye(n)*c)/M+sqrtm((-1i*(x-eye(n)*c)/M)^2-eye(n)))/R+eye(n))/((-1i*(x-eye(n)*c)/M+sqrtm((-1i*(x-eye(n)*c)/M)^2-eye(n)))/R-eye(n));
% dpsiinv = @(x) -2*(inv(sqrtm(x*(x^2-eye(n))))*x+eye(n))/((x+sqrtm(x^2-eye(n)))*inv(R)-eye(n))^2;
% Alyap = psiinv(AA);
% Blyap = sqrt(dpsiinv(AA))' * BB*BB' * sqrt(dpsiinv(AA));
% Xc = lyap(Alyap,Blyap);
% Clyap = sqrt(dpsiinv(AA)) * CC'*CC * sqrt(dpsiinv(AA))';
% Xo = lyap(Alyap,Clyap);

% 
% Xc = 0.5*(Xc + Xc')+1e-16*eye(n);
% Xo = 0.5*(Xo + Xo')+1e-16*eye(n);
% U = chol(Xc); U = U'; 
% L = chol(Xo); L = L';
% [Z,S,Y] = svd(L'*U, 'econ');
% 
% 
% 
% Z1 = Z(:,1:r);
% Y1 = Y(:,1:r);
% S1 = S(1:r,1:r);  S1half = sqrt(S1);
% 
% Wr = L*Z1/S1half;
% Vr = U*Y1/S1half;
% 
% Ar = (Wr'*Vr)\(Wr'*AA*Vr);
% Br = (Wr'*Vr)\(Wr'*BB);
% Cr = CC*Vr;

%% Running Algorithm conformalBT with cubed square root
% r = 13; %13
% psi = @(x) (-x).^(5/4);
% dpsi = @(x) -(5/4)*(-x).^(5/4);
% psiinv = @(x) -(x).^(5/4);
% con = @(z) (psi(1i.*z)*eye(n)-AA)\BB * BB'/(psi(1i.*z)*eye(n)-AA)' * sqrt(dpsi(1i.*z)) * sqrt(dpsi(1i.*z))' ;
% obs = @(z) (psi(1i.*z)*eye(n)-AA)'\CC' * CC/(psi(1i.*z)*eye(n)-AA) * sqrt(dpsi(1i.*z)) * sqrt(dpsi(1i.*z))' ;
% 
% 
% P = integral(con,-Inf,Inf,'ArrayValued',true)/(2*pi);
% P = 0.5*(P + P')+1e-12*eye(n);
% Q = integral(obs,-Inf,Inf,'ArrayValued',true)/(2*pi);
% Q = 0.5*(Q + Q')+1e-12*eye(n);
% U = chol(P); U = U'; 
% L = chol(Q); L = L';
% [Z,S,Y] = svd(L'*U, 'econ');
% 
% 
% 
% Z1 = Z(:,1:r);
% Y1 = Y(:,1:r);
% S1 = S(1:r,1:r);  S1half = sqrt(S1);
% 
% Wr = L*Z1/S1half;
% Vr = U*Y1/S1half;
% 
% Ar = (Wr'*Vr)\(Wr'*AA*Vr);
% Br = (Wr'*Vr)\(Wr'*BB);
% Cr = CC*Vr;
% 
% eigAr = eig(Ar);
% eigA = eig(full(AA));
% figure()
% plot(real(eigA),imag(eigA),'ob'); hold on
% plot(real(eig(Ar)),imag(eig(Ar)),'rx')
% plot(real(psiinv(eigA)),imag(psiinv(eigA)),'k^'); hold off
% legend('$\Lambda(A)$','$\Lambda(A_r)$','$\psi^{-1}(\Lambda(A))$','interpreter','latex','fontsize',20)
% %test
% 
% % Alyap = psiinv(AA);
% % Blyap = -sqrt(dpsiinv(AA))' * BB*BB' * sqrt(dpsiinv(AA));
% % Xc = lyap(Alyap,Blyap);
% % Xc = 0.5*(Xc + Xc')+1e-12*eye(n);
% % Clyap = -sqrt(dpsiinv(AA)) * CC'*CC * sqrt(dpsiinv(AA))';
% % Xo = lyap(Alyap,Clyap);
% % Xo = 0.5*(Xo + Xo')+1e-12*eye(n);
% conr = @(z) (psi(1i.*z)*eye(r)-Ar)\Br * Br'/(psi(1i.*z)*eye(r)-Ar)' * sqrt(dpsi(1i.*z)) * sqrt(dpsi(1i.*z))';
% obsr = @(z) (psi(1i.*z)*eye(r)-Ar)'\Cr' * Cr/(psi(1i.*z)*eye(r)-Ar) * sqrt(dpsi(1i.*z)) * sqrt(dpsi(1i.*z))';
% Pr = integral(conr,-Inf,Inf,'ArrayValued',true)/(2*pi);
% Qr = integral(obsr,-Inf,Inf,'ArrayValued',true)/(2*pi);
% 

%% Computing systems output trajectories for impulse response
dynamics = @(t,x,A,B) A*x;
[t1,x] = ode23(dynamics,linspace(0,10,5000),BB,[],AA,BB);
y1 = CC*x';
[t2,xr] = ode23(dynamics,linspace(0,10,5000),Br,[],Ar,Br);
y2 = Cr*xr.';
y3 = y1-y2;

%% Plots
figure()
set(gcf,'position',[100,100,1100,500])
subplot(2,1,1)
plot(t1,real(y1),'r-', 'Linewidth', 3); hold on
plot(t2,real(y2),'b--', 'Linewidth', 3); hold off
ax = gca;
ax.FontSize = 14;
legend('$y(t)$','$\widehat{y}_r(t)$','fontsize',20, 'interpreter','latex', 'Location', 'southeast')
subplot(2,1,2)
plot(t2,abs(real(y3)),'k', 'Linewidth', 1.5)
ax = gca;
ax.FontSize = 14;
xlabel('time [s]','fontsize',20,'interpreter','latex')
legend('$|y(t)-\widehat{y}_r(t)|$','fontsize',20, 'interpreter','latex', 'Location', 'northwest')

% cleanfigure;
% matlab2tikz('betterconfBTshift.tex')
