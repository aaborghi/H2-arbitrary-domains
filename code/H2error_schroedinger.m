clear; clc; 
rng(4)

%% Constructing the discretized Schroedinger equation
nx = 1000;
xa = 0;
xb = 1;
nu = 1;
hx = (xb-xa)/(nx+1);
xd = xa+hx:hx:xb-hx;
ex = ones(nx,1);
I = speye(nx);
Laplace_x = 1/hx^2*spdiags([ex -2*ex ex], -1:1, nx, nx);
e1x = I(:,1);
enx = I(:,nx);
O = sparse(nx,nx);
A = nu*Laplace_x;
B = zeros(nx,1);
B(end) = -1i*nu*1/hx^2;
C = ones(1,nx)*hx;
n = size(A,1);
A=-1i*A;

%% Running Algorithm1 and computing the H2(\barA^c) error norm
m = 25;
maxiter = 100;
H2Arel = zeros(1,m);
H2Arel_irka = zeros(1,m);
n = size(A,1);
  
rinit = 5;
for r = rinit:5:m
    init = 100*randn(r,1)-1e2i;
    phi = @(z) (conj(z));
    % Algorithm 1
    [Ar,Br,Cr,s] = algorithm1(A,B,C,r,phi,init,maxiter); 
    init = 10*1i*[1:r]+1i*1e4;
    phi2 = @(z) (-conj(z));
    % IRKA
    [Ar_irka,Br_irka,Cr_irka,~] = algorithm1(A,B,C,r,phi2,init,maxiter); 
    psi = @(z) (-1i*z);
    dpsi = @(z) (-1i);
    [~,H2Arel(1,r)] = H2Anorm(A,B,C,Ar,Br,Cr,psi,dpsi);
    [~,H2Arel_irka(1,r)] = H2Anorm(A,B,C,Ar_irka,Br_irka,Cr_irka,psi,dpsi);
    
end

%% Plots
figure()
set(gcf,'position',[100,100,550,500])
semilogy(rinit:5:m,H2Arel(1,rinit:5:m),'k-o','Linewidth', 2); hold on
semilogy(rinit:5:m,H2Arel_irka(1,rinit:5:m),'b--x','Linewidth', 2); 
ax = gca;
ax.FontSize = 14; 
xlabel('$r$','fontsize',20,'interpreter','latex')
ylabel('Relative error norm','fontsize',20,'interpreter','latex')
legend('Algorithm 1', 'IRKA', 'fontsize',20, 'interpreter','latex','Location', 'southwest')
xlim([5,25])
saveas(gcf,'schroedingerH2error.eps', 'epsc')
